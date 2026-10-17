//! Index loading and data structures.

pub mod bifrost;
pub mod builder;
mod graph_build;
mod workers;
mod writer;

pub use builder::{IndexBuildOptions, IndexBuildReport, build_index, build_index_with_report};

use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
use std::path::Path;

use croaring_sys::{
    roaring_bitmap_deserialize, roaring_bitmap_free, roaring_bitmap_get_cardinality,
    roaring_bitmap_to_uint32_array,
};

use crate::ec::EcList;
use crate::{Error, Result};

const SIZE_T_BYTES: usize = 8;
const INDEX_VERSION_CURRENT: u64 = 13;
const KMER_BYTES_CANDIDATES: [usize; 4] = [8, 16, 24, 32];
const BFG_METABIN_FORMAT_HEADER: u64 = 0x267c3d5d;

/// Basic transcript metadata.
#[derive(Debug, Clone)]
pub struct TranscriptInfo {
    pub name: String,
    pub length: u32,
}

/// Loaded index metadata and lookup tables.
#[derive(Debug, Clone)]
pub struct Index {
    pub k: u32,
    pub index_version: u64,
    pub minimizer_len: Option<u32>,
    pub unitigs: Option<u64>,
    pub kmers: Option<u64>,
    pub transcripts: Vec<TranscriptInfo>,
}

impl Index {
    /// Load a kallisto index with partial parsing (metadata only).
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
        let mut reader = BufReader::new(file);

        let index_version = read_u64_le(&mut reader)?;
        if index_version != INDEX_VERSION_CURRENT {
            return Err(Error::UnsupportedFeature(format!(
                "unsupported index version {index_version}"
            )));
        }
        let dbg_size = read_u64_le(&mut reader)?;

        let mut k: Option<u32> = None;
        let mut minimizer_len: Option<u32> = None;
        let mut unitigs: Option<u64> = None;
        let mut kmers: Option<u64> = None;
        let mut kmer_bytes: Option<usize> = None;

        if dbg_size > 0 {
            let dbg_start = reader.stream_position()?;
            let meta =
                parse_graph_section(&mut reader, dbg_start, dbg_size, &KMER_BYTES_CANDIDATES)?;
            k = Some(meta.k);
            minimizer_len = Some(meta.g);
            unitigs = Some(meta.unitigs);
            kmers = Some(meta.kmers);
            kmer_bytes = Some(meta.kmer_bytes);

            let mphf_size = read_u64_le(&mut reader)?;
            reader.seek(SeekFrom::Current(mphf_size as i64))?;
        }

        let dlist_size = read_u64_le(&mut reader)?;
        let _dlist_overhang = read_u64_le(&mut reader)?;

        let file_len = reader.get_ref().metadata()?.len();
        let after_dlist_pos = reader.stream_position()?;

        let parse = parse_after_dlist(
            &mut reader,
            after_dlist_pos,
            file_len,
            dlist_size,
            k,
            kmer_bytes,
        )?;

        let ParseResult {
            kmer_len,
            transcripts,
        } = parse;

        let k = k.ok_or_else(|| Error::UnsupportedFeature("k-mer length missing".into()))?;
        let _ = kmer_len;
        Ok(Index {
            k,
            index_version,
            minimizer_len,
            unitigs,
            kmers,
            transcripts,
        })
    }
}

/// Extract ECs from a kallisto index by scanning node EC blocks.
pub fn extract_ec_list(path: &Path) -> Result<EcList> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let index_version = read_u64_le(&mut reader)?;
    if index_version != INDEX_VERSION_CURRENT {
        return Err(Error::UnsupportedFeature(format!(
            "unsupported index version {index_version}"
        )));
    }
    let dbg_size = read_u64_le(&mut reader)?;

    let mut k: Option<u32> = None;
    let mut kmer_bytes: Option<usize> = None;

    if dbg_size > 0 {
        let dbg_start = reader.stream_position()?;
        let meta = parse_graph_section(&mut reader, dbg_start, dbg_size, &KMER_BYTES_CANDIDATES)?;
        k = Some(meta.k);
        kmer_bytes = Some(meta.kmer_bytes);

        let mphf_size = read_u64_le(&mut reader)?;
        reader.seek(SeekFrom::Current(mphf_size as i64))?;
    }

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let kmer_bytes =
        kmer_bytes.ok_or_else(|| Error::InvalidFormat("missing k-mer width".into()))?;
    let dlist_skip = dlist_size
        .checked_mul(kmer_bytes as u64)
        .ok_or_else(|| Error::InvalidFormat("d-list size overflow".into()))?;
    reader.seek(SeekFrom::Current(dlist_skip as i64))?;

    let node_count = read_u64_le(&mut reader)?;
    let k = k.ok_or_else(|| Error::InvalidFormat("missing k-mer length".into()))?;

    let mut classes: Vec<Vec<u32>> = Vec::new();
    let mut map: std::collections::HashMap<Vec<u32>, usize> = std::collections::HashMap::new();

    for _ in 0..node_count {
        reader.seek(SeekFrom::Current(k as i64))?;
        let node_size = read_u32_le(&mut reader)? as usize;
        let mut buf = vec![0u8; node_size];
        reader.read_exact(&mut buf)?;

        let mut cur = Cursor::new(buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let block_ecs = read_block_array_ecs(&mut cur)?;
        for ec in block_ecs {
            if !map.contains_key(&ec) {
                let id = classes.len();
                map.insert(ec.clone(), id);
                classes.push(ec);
            }
        }
    }

    Ok(EcList { classes })
}

struct ParseResult {
    kmer_len: u32,
    transcripts: Vec<TranscriptInfo>,
}

fn parse_after_dlist<R: Read + Seek>(
    reader: &mut R,
    after_dlist_pos: u64,
    file_len: u64,
    dlist_size: u64,
    k: Option<u32>,
    kmer_bytes_hint: Option<usize>,
) -> Result<ParseResult> {
    let kmer_len = k.ok_or_else(|| Error::UnsupportedFeature("k-mer length missing".into()))?;
    let mut last_err: Option<Error> = None;

    let candidates: Vec<usize> = match kmer_bytes_hint {
        Some(v) => vec![v],
        None => KMER_BYTES_CANDIDATES.to_vec(),
    };

    for &kmer_bytes in &candidates {
        if let Err(err) = reader.seek(SeekFrom::Start(after_dlist_pos)) {
            return Err(err.into());
        }
        let skip = match dlist_size.checked_mul(kmer_bytes as u64) {
            Some(v) => v,
            None => {
                last_err = Some(Error::InvalidFormat("d-list size overflow".into()));
                continue;
            }
        };
        if after_dlist_pos + skip > file_len {
            last_err = Some(Error::InvalidFormat("d-list exceeds file size".into()));
            continue;
        }
        reader.seek(SeekFrom::Current(skip as i64))?;

        let node_count = match read_u64_le(reader) {
            Ok(v) => v,
            Err(e) => {
                last_err = Some(e);
                continue;
            }
        };

        if let Err(err) = skip_nodes(reader, kmer_len, node_count, file_len) {
            last_err = Some(err);
            continue;
        }

        let num_trans = match read_i32_le(reader) {
            Ok(v) => v,
            Err(e) => {
                last_err = Some(e);
                continue;
            }
        };
        if num_trans < 0 {
            last_err = Some(Error::InvalidFormat("negative transcript count".into()));
            continue;
        }
        let num_trans = num_trans as usize;

        let mut lengths = Vec::with_capacity(num_trans);
        let mut ok = true;
        for _ in 0..num_trans {
            let len = match read_i32_le(reader) {
                Ok(v) => v,
                Err(e) => {
                    last_err = Some(e);
                    ok = false;
                    break;
                }
            };
            if len < 0 {
                last_err = Some(Error::InvalidFormat("negative transcript length".into()));
                ok = false;
                break;
            }
            lengths.push(len as u32);
        }
        if !ok {
            continue;
        }

        let mut transcripts = Vec::with_capacity(num_trans);
        for (idx, len) in lengths.iter().copied().enumerate() {
            if idx >= num_trans {
                break;
            }
            let name_len = match read_u64_le(reader) {
                Ok(v) => v as usize,
                Err(e) => {
                    last_err = Some(e);
                    ok = false;
                    break;
                }
            };
            let mut buf = vec![0u8; name_len];
            if let Err(e) = reader.read_exact(&mut buf) {
                last_err = Some(e.into());
                ok = false;
                break;
            }
            let name = String::from_utf8_lossy(&buf).to_string();
            transcripts.push(TranscriptInfo { name, length: len });
        }
        if !ok {
            continue;
        }

        let onlist_size = match read_u64_le(reader) {
            Ok(v) => v,
            Err(e) => {
                last_err = Some(e);
                continue;
            }
        };
        if onlist_size > 0
            && let Err(e) = reader.seek(SeekFrom::Current(onlist_size as i64))
        {
            last_err = Some(e.into());
            continue;
        }

        return Ok(ParseResult {
            kmer_len,
            transcripts,
        });
    }

    Err(last_err.unwrap_or_else(|| Error::InvalidFormat("failed to parse index".into())))
}

fn skip_nodes<R: Read + Seek>(
    reader: &mut R,
    kmer_len: u32,
    node_count: u64,
    file_len: u64,
) -> Result<()> {
    for _ in 0..node_count {
        let pos = reader.stream_position()?;
        let header_end = pos + kmer_len as u64 + 4;
        if header_end > file_len {
            return Err(Error::InvalidFormat("node header exceeds file size".into()));
        }
        reader.seek(SeekFrom::Current(kmer_len as i64))?;
        let node_size = read_u32_le(reader)?;
        let end = header_end + node_size as u64;
        if end > file_len {
            return Err(Error::InvalidFormat(
                "node payload exceeds file size".into(),
            ));
        }
        reader.seek(SeekFrom::Current(node_size as i64))?;
    }
    Ok(())
}

pub(crate) struct GraphMeta {
    pub k: u32,
    pub g: u32,
    pub unitigs: u64,
    pub kmers: u64,
    pub kmer_bytes: usize,
}

pub(crate) fn parse_graph_section<R: Read + Seek>(
    reader: &mut R,
    dbg_start: u64,
    dbg_size: u64,
    kmer_candidates: &[usize],
) -> Result<GraphMeta> {
    let file_format_version = read_u64_le(reader)?;
    let kmer_len = read_i32_le(reader)?;
    let minimizer_len = read_i32_le(reader)?;
    let _ = file_format_version;

    if kmer_len <= 0 || minimizer_len <= 0 {
        return Err(Error::InvalidFormat("invalid graph header".into()));
    }

    let header_bytes = (SIZE_T_BYTES + 4 + 4) as u64;
    if dbg_size < header_bytes {
        return Err(Error::InvalidFormat("corrupt dBG header".into()));
    }

    let dbg_end = dbg_start + dbg_size;
    let mut last_err: Option<Error> = None;

    for &kmer_bytes in kmer_candidates {
        reader.seek(SeekFrom::Start(dbg_start + header_bytes))?;

        let v_unitigs_sz = match read_u64_le(reader) {
            Ok(v) => v,
            Err(e) => {
                last_err = Some(e);
                continue;
            }
        };

        let mut kmers_count: u64 = 0;
        for _ in 0..v_unitigs_sz {
            let len = match read_u64_le(reader) {
                Ok(v) => v,
                Err(e) => {
                    last_err = Some(e);
                    break;
                }
            };
            let data_sz = len.div_ceil(4);
            let end = reader.stream_position()? + data_sz;
            if end > dbg_end {
                last_err = Some(Error::InvalidFormat(
                    "unitig sequence exceeds graph size".into(),
                ));
                break;
            }
            reader.seek(SeekFrom::Current(data_sz as i64))?;
            if len >= kmer_len as u64 {
                kmers_count += len - kmer_len as u64 + 1;
            }
        }
        if let Some(err) = last_err.as_ref() {
            let _ = err;
            continue;
        }

        let km_unitigs_sz = match read_u64_le(reader) {
            Ok(v) => v,
            Err(e) => {
                last_err = Some(e);
                continue;
            }
        };
        let km_bytes = km_unitigs_sz
            .checked_mul(kmer_bytes as u64)
            .ok_or_else(|| Error::InvalidFormat("km_unitigs size overflow".into()))?;
        if reader.stream_position()? + km_bytes > dbg_end {
            last_err = Some(Error::InvalidFormat("km_unitigs exceeds graph size".into()));
            continue;
        }
        reader.seek(SeekFrom::Current(km_bytes as i64))?;
        kmers_count += km_unitigs_sz;

        let h_kmers_ccov_sz = match read_u64_le(reader) {
            Ok(v) => v,
            Err(e) => {
                last_err = Some(e);
                continue;
            }
        };
        let h_bytes = h_kmers_ccov_sz
            .checked_mul(kmer_bytes as u64)
            .ok_or_else(|| Error::InvalidFormat("h_kmers size overflow".into()))?;
        if reader.stream_position()? + h_bytes > dbg_end {
            last_err = Some(Error::InvalidFormat("h_kmers exceeds graph size".into()));
            continue;
        }
        reader.seek(SeekFrom::Current(h_bytes as i64))?;
        kmers_count += h_kmers_ccov_sz;

        let pos_after_graph = reader.stream_position()?;
        if pos_after_graph > dbg_end {
            last_err = Some(Error::InvalidFormat("graph section exceeds size".into()));
            continue;
        }

        if pos_after_graph + SIZE_T_BYTES as u64 > dbg_end {
            last_err = Some(Error::InvalidFormat("missing graph index header".into()));
            continue;
        }

        let meta_header = read_u64_le(reader)?;
        if (meta_header >> 32) != BFG_METABIN_FORMAT_HEADER {
            last_err = Some(Error::InvalidFormat("invalid graph index header".into()));
            continue;
        }

        reader.seek(SeekFrom::Start(dbg_end))?;

        let unitigs = v_unitigs_sz + km_unitigs_sz + h_kmers_ccov_sz;
        return Ok(GraphMeta {
            k: kmer_len as u32,
            g: minimizer_len as u32,
            unitigs,
            kmers: kmers_count,
            kmer_bytes,
        });
    }

    Err(last_err.unwrap_or_else(|| Error::InvalidFormat("failed to parse graph section".into())))
}

fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u32_le<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_i32_le<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

pub(crate) struct EcBlock {
    pub lb: u32,
    pub ub: u32,
    pub ec: Vec<u32>,
    pub positions: Option<Vec<PositionSet>>,
    pub strands: Option<Vec<u8>>,
}

pub(crate) enum PositionSet {
    Single(u32),
    Many(Vec<u32>),
}

impl PositionSet {
    fn from_values(mut values: Vec<u32>) -> Self {
        match values.len() {
            0 => Self::Many(values),
            1 => Self::Single(values[0]),
            _ => {
                values.sort_unstable();
                values.dedup();
                Self::Many(values)
            }
        }
    }

    pub(crate) fn contains(&self, value: u32) -> bool {
        match self {
            Self::Single(raw) => *raw == value,
            Self::Many(values) => values.binary_search(&value).is_ok(),
        }
    }

    pub(crate) fn min_max(&self) -> Option<(u32, u32)> {
        match self {
            Self::Single(raw) => Some((*raw, *raw)),
            Self::Many(values) => Some((*values.first()?, *values.last()?)),
        }
    }

    pub(crate) fn for_each(&self, mut f: impl FnMut(u32)) {
        match self {
            Self::Single(raw) => f(*raw),
            Self::Many(values) => {
                for &raw in values {
                    f(raw);
                }
            }
        }
    }
}

fn read_block_array_ecs<R: Read>(reader: &mut R) -> Result<Vec<Vec<u32>>> {
    let blocks = read_block_array_blocks(reader)?;
    Ok(blocks.into_iter().map(|b| b.ec).collect())
}

pub(crate) fn read_block_array_blocks<R: Read>(reader: &mut R) -> Result<Vec<EcBlock>> {
    let flag = read_u8(reader)?;
    match flag {
        0 => Ok(Vec::new()),
        1 => {
            let lb = read_u32_le(reader)?;
            let ub = read_u32_le(reader)?;
            let ec = read_sparse_vector_ec(reader)?;
            Ok(vec![EcBlock {
                lb,
                ub,
                ec,
                positions: None,
                strands: None,
            }])
        }
        2 => {
            let count = read_u64_le(reader)?;
            let mut ecs = Vec::with_capacity(count as usize);
            for _ in 0..count {
                let lb = read_u32_le(reader)?;
                let ub = read_u32_le(reader)?;
                let ec = read_sparse_vector_ec(reader)?;
                ecs.push(EcBlock {
                    lb,
                    ub,
                    ec,
                    positions: None,
                    strands: None,
                });
            }
            Ok(ecs)
        }
        _ => Err(Error::InvalidFormat("invalid block array flag".into())),
    }
}

pub(crate) fn read_block_array_blocks_with_positions<R: Read>(
    reader: &mut R,
) -> Result<Vec<EcBlock>> {
    let flag = read_u8(reader)?;
    match flag {
        0 => Ok(Vec::new()),
        1 => {
            let lb = read_u32_le(reader)?;
            let ub = read_u32_le(reader)?;
            let (ec, positions, strands) = read_sparse_vector_with_positions(reader)?;
            Ok(vec![EcBlock {
                lb,
                ub,
                ec,
                positions: Some(positions),
                strands: Some(strands),
            }])
        }
        2 => {
            let count = read_u64_le(reader)?;
            let mut ecs = Vec::with_capacity(count as usize);
            for _ in 0..count {
                let lb = read_u32_le(reader)?;
                let ub = read_u32_le(reader)?;
                let (ec, positions, strands) = read_sparse_vector_with_positions(reader)?;
                ecs.push(EcBlock {
                    lb,
                    ub,
                    ec,
                    positions: Some(positions),
                    strands: Some(strands),
                });
            }
            Ok(ecs)
        }
        _ => Err(Error::InvalidFormat("invalid block array flag".into())),
    }
}

pub(crate) fn read_sparse_vector_ec<R: Read>(reader: &mut R) -> Result<Vec<u32>> {
    let roaring_size = read_u64_le(reader)? as usize;
    let mut buf = vec![0u8; roaring_size];
    reader.read_exact(&mut buf)?;
    let ids = unsafe { deserialize_roaring_to_vec(&buf) }
        .ok_or_else(|| Error::InvalidFormat("failed to deserialize roaring bitmap".into()))?;

    let v_size = read_u64_le(reader)? as usize;
    let mut skip = [0u8; 1024];
    for _ in 0..v_size {
        let mut remaining = read_u64_le(reader)? as usize;
        while remaining != 0 {
            let chunk = remaining.min(skip.len());
            reader.read_exact(&mut skip[..chunk])?;
            remaining -= chunk;
        }
    }

    if ids.len() != v_size {
        let mut ids = ids;
        ids.sort_unstable();
        return Ok(ids);
    }
    Ok(ids)
}

pub(crate) fn read_sparse_vector_with_positions<R: Read>(
    reader: &mut R,
) -> Result<(Vec<u32>, Vec<PositionSet>, Vec<u8>)> {
    let roaring_size = read_u64_le(reader)? as usize;
    let mut buf = vec![0u8; roaring_size];
    reader.read_exact(&mut buf)?;
    let ids = unsafe { deserialize_roaring_to_vec(&buf) }
        .ok_or_else(|| Error::InvalidFormat("failed to deserialize roaring bitmap".into()))?;

    let v_size = read_u64_le(reader)? as usize;
    let mut positions = Vec::with_capacity(v_size);
    let mut strands = Vec::with_capacity(v_size);
    let mut pbuf = Vec::new();
    for _ in 0..v_size {
        let p_size = read_u64_le(reader)? as usize;
        pbuf.resize(p_size, 0);
        reader.read_exact(&mut pbuf)?;
        let raw_positions = unsafe { deserialize_roaring_to_vec(&pbuf) }
            .ok_or_else(|| Error::InvalidFormat("failed to deserialize roaring bitmap".into()))?;
        let raw_positions = PositionSet::from_values(raw_positions);
        let (min, max) = raw_positions
            .min_max()
            .ok_or_else(|| Error::InvalidFormat("empty position bitmap".into()))?;
        let min_sense = (min & 0x7fff_ffff) == min;
        let max_sense = (max & 0x7fff_ffff) == max;
        let strand = if min_sense == max_sense {
            if min_sense { 1 } else { 0 }
        } else {
            2
        };
        positions.push(raw_positions);
        strands.push(strand);
    }

    let mut ids = ids;
    if ids.len() != v_size {
        ids.sort_unstable();
    }
    if ids.len() != positions.len() {
        return Err(Error::InvalidFormat("sparse vector size mismatch".into()));
    }
    if ids.len() != strands.len() {
        return Err(Error::InvalidFormat("sparse vector size mismatch".into()));
    }
    Ok((ids, positions, strands))
}

unsafe fn deserialize_roaring_to_vec(buf: &[u8]) -> Option<Vec<u32>> {
    let mut out = Vec::new();
    unsafe { append_roaring_to_vec(buf, &mut out) }?;
    Some(out)
}

pub(crate) unsafe fn append_roaring_to_vec(buf: &[u8], out: &mut Vec<u32>) -> Option<()> {
    if buf.first() == Some(&1) {
        let count = u32::from_le_bytes(buf.get(1..5)?.try_into().ok()?) as usize;
        let end = count.checked_mul(4)?.checked_add(5)?;
        let payload = buf.get(5..end)?;
        let start = out.len();
        out.extend(
            payload
                .as_chunks::<4>()
                .0
                .iter()
                .map(|word| u32::from_le_bytes(*word)),
        );
        if out[start..].windows(2).any(|pair| pair[0] >= pair[1]) {
            out[start..].sort_unstable();
            let mut write = start;
            for read in start..out.len() {
                let value = out[read];
                if write == start || out[write - 1] != value {
                    out[write] = value;
                    write += 1;
                }
            }
            out.truncate(write);
        }
        return Some(());
    }
    if buf.is_empty() {
        return None;
    }
    let ptr = unsafe { roaring_bitmap_deserialize(buf.as_ptr().cast()) };
    if ptr.is_null() {
        return None;
    }
    let card = unsafe { roaring_bitmap_get_cardinality(ptr) } as usize;
    let start = out.len();
    out.resize(start + card, 0);
    unsafe {
        roaring_bitmap_to_uint32_array(ptr, out.as_mut_ptr().add(start));
        roaring_bitmap_free(ptr);
    }
    Some(())
}

fn node_eof() -> Error {
    std::io::Error::new(
        std::io::ErrorKind::UnexpectedEof,
        "failed to fill whole buffer",
    )
    .into()
}

pub(crate) fn read_block_array_slice(
    bytes: &[u8],
    mut block: impl FnMut(u32, u32, &[u8]) -> Result<()>,
) -> Result<()> {
    let mut reader = Cursor::new(bytes);
    let count = match read_u8(&mut reader)? {
        0 => 0,
        1 => 1,
        2 => read_u64_le(&mut reader)?,
        _ => return Err(Error::InvalidFormat("invalid block array flag".into())),
    };
    for _ in 0..count {
        let lb = read_u32_le(&mut reader)?;
        let ub = read_u32_le(&mut reader)?;
        let size = read_u64_le(&mut reader)? as usize;
        let start = reader.position() as usize;
        let end = start.checked_add(size).ok_or_else(node_eof)?;
        let bitmap = bytes.get(start..end).ok_or_else(node_eof)?;
        block(lb, ub, bitmap)?;
        reader.set_position(end as u64);
        let positions = read_u64_le(&mut reader)?;
        for _ in 0..positions {
            let size = read_u64_le(&mut reader)? as usize;
            let start = reader.position() as usize;
            let end = start.checked_add(size).ok_or_else(node_eof)?;
            bytes.get(start..end).ok_or_else(node_eof)?;
            reader.set_position(end as u64);
        }
    }
    Ok(())
}

#[cfg(test)]
mod skipped_position_tests {
    use super::*;

    #[test]
    fn direct_array_decode_equals_c_bitmap_set_semantics() {
        for ids in [
            vec![],
            vec![7u32],
            vec![7, 7],
            vec![9, 3, 9, 1],
            vec![0, u32::MAX, 100_000, 65_535, 65_536],
        ] {
            let mut bytes = vec![1u8];
            bytes.extend_from_slice(&(ids.len() as u32).to_le_bytes());
            for id in ids {
                bytes.extend_from_slice(&id.to_le_bytes());
            }
            let expected = unsafe {
                let bitmap = roaring_bitmap_deserialize(bytes.as_ptr().cast());
                assert!(!bitmap.is_null());
                let mut expected = vec![0; roaring_bitmap_get_cardinality(bitmap) as usize];
                roaring_bitmap_to_uint32_array(bitmap, expected.as_mut_ptr());
                roaring_bitmap_free(bitmap);
                expected
            };
            assert_eq!(
                unsafe { deserialize_roaring_to_vec(&bytes) }.unwrap(),
                expected
            );
            bytes.extend_from_slice(&[19, 7, 3]);
            assert_eq!(
                unsafe { deserialize_roaring_to_vec(&bytes) }.unwrap(),
                expected
            );
        }
        for bytes in [
            vec![],
            vec![1],
            vec![1, 1, 0, 0, 0],
            vec![1, 2, 0, 0, 0, 5, 0, 0, 0],
        ] {
            assert!(unsafe { deserialize_roaring_to_vec(&bytes) }.is_none());
        }
    }

    fn sparse_bytes(sizes: &[usize]) -> Vec<u8> {
        let mut bitmap = vec![1u8];
        bitmap.extend_from_slice(&3u32.to_le_bytes());
        for id in [1u32, 7, 33] {
            bitmap.extend_from_slice(&id.to_le_bytes());
        }
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&(bitmap.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&bitmap);
        bytes.extend_from_slice(&(sizes.len() as u64).to_le_bytes());
        for &size in sizes {
            bytes.extend_from_slice(&(size as u64).to_le_bytes());
            bytes.extend(std::iter::repeat_n(91u8, size));
        }
        bytes
    }

    #[test]
    fn slice_blocks_match_generic_blocks_and_preserve_existing_prefix() {
        for count in 0..=3u64 {
            let mut bytes = if count == 0 {
                vec![0]
            } else if count == 1 {
                vec![1]
            } else {
                let mut bytes = vec![2];
                bytes.extend_from_slice(&count.to_le_bytes());
                bytes
            };
            for i in 0..count {
                bytes.extend_from_slice(&(i as u32).to_le_bytes());
                bytes.extend_from_slice(&((i + 1) as u32).to_le_bytes());
                bytes.extend_from_slice(&sparse_bytes(&[1, 1025, 7]));
            }
            let expected: Vec<_> = read_block_array_blocks(&mut Cursor::new(&bytes))
                .unwrap()
                .into_iter()
                .map(|b| (b.lb, b.ub, b.ec))
                .collect();
            let mut actual = Vec::new();
            read_block_array_slice(&bytes, |lb, ub, bitmap| {
                let mut ids = vec![9, 9];
                unsafe { append_roaring_to_vec(bitmap, &mut ids) }.unwrap();
                assert_eq!(&ids[..2], &[9, 9]);
                actual.push((lb, ub, ids[2..].to_vec()));
                Ok(())
            })
            .unwrap();
            assert_eq!(actual, expected);
            for cut in [0, bytes.len() - 1] {
                let old = read_block_array_blocks(&mut Cursor::new(&bytes[..cut]));
                let new = read_block_array_slice(&bytes[..cut], |_, _, bitmap| {
                    unsafe { deserialize_roaring_to_vec(bitmap) }.ok_or_else(|| {
                        Error::InvalidFormat("failed to deserialize roaring bitmap".into())
                    })?;
                    Ok(())
                });
                assert_eq!(old.is_err(), new.is_err());
            }
        }
    }

    #[test]
    fn skipped_positions_consume_exact_bytes_and_reject_truncated_payloads() {
        for sizes in [[0usize, 0, 0], [1, 1024, 1025], [4099, 3, 5]] {
            let bytes = sparse_bytes(&sizes);
            let mut input = bytes.clone();
            input.push(173);
            let mut reader = std::io::Cursor::new(input);
            assert_eq!(read_sparse_vector_ec(&mut reader).unwrap(), [1, 7, 33]);
            assert_eq!(reader.position() as usize, bytes.len());
            assert_eq!(read_u8(&mut reader).unwrap(), 173);
            for cut in [0, 7, bytes.len() - 1] {
                assert!(read_sparse_vector_ec(&mut std::io::Cursor::new(&bytes[..cut])).is_err());
            }
        }
    }
}

mod packed_transcripts;

mod scratch_pages;

mod transcript_runs;

mod block_positions;
