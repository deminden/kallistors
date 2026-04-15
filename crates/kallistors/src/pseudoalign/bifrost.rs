use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::ops::Index;
use std::path::Path;

use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::index::bifrost::{BooPhf, encode_minimizer_rep, read_bitcontainer_values};
use crate::index::{parse_graph_section, read_block_array_blocks};
use crate::{Error, Result};

use super::ec::encode_kmer_pair;
use super::io::{read_i32_le, read_u32_le, read_u64_le};
use super::{KMER_BYTES_CANDIDATES, KmerEcIndex, read_onlist_with_lengths};

pub type KmerPosIndex = HashMap<u64, Vec<(usize, usize, bool)>>;

pub struct MinzPositionIndex {
    offsets: Vec<usize>,
    values: Vec<u64>,
}

impl MinzPositionIndex {
    fn from_counts(counts: Vec<usize>) -> Self {
        let mut offsets = Vec::with_capacity(counts.len() + 1);
        offsets.push(0);
        let mut total = 0usize;
        for count in counts {
            total += count;
            offsets.push(total);
        }
        Self {
            offsets,
            values: vec![0; total],
        }
    }

    fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    fn cursors(&self) -> Vec<usize> {
        self.offsets[..self.len()].to_vec()
    }

    pub fn get(&self, idx: usize) -> &[u64] {
        &self.values[self.offsets[idx]..self.offsets[idx + 1]]
    }

    fn push_at(&mut self, cursors: &mut [usize], idx: usize, pos_id: u64) {
        let cursor = &mut cursors[idx];
        self.values[*cursor] = pos_id;
        *cursor += 1;
    }
}

impl Index<usize> for MinzPositionIndex {
    type Output = [u64];

    fn index(&self, index: usize) -> &Self::Output {
        self.get(index)
    }
}

pub struct BifrostIndex {
    pub k: usize,
    pub g: usize,
    pub unitigs: Vec<Vec<u8>>,
    pub km_unitigs: Vec<Vec<u8>>,
    pub h_kmers: Vec<Vec<u8>>,
    pub h_kmer_map: Option<HashMap<u64, usize>>,
    pub dlist: Option<HashSet<u64>>,
    pub kmer_pos_index: Option<KmerPosIndex>,
    pub kmer_index: Option<KmerEcIndex>,
    pub(crate) ec_blocks: Vec<Vec<crate::index::EcBlock>>,
    pub minz_positions: MinzPositionIndex,
    pub mphf: BooPhf,
    pub onlist: Option<Vec<bool>>,
    pub transcript_lengths: Vec<u32>,
    pub transcript_names: Vec<String>,
    pub shade_to_color: Vec<Option<u32>>,
    pub shade_sequences: Vec<bool>,
    pub use_shade: bool,
}

pub fn build_bifrost_index(path: &Path) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, false, 1)
}

pub fn build_bifrost_index_with_kmer_pos(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index_with_positions_threaded(path, load_positional_info, 1)?;
    index.kmer_pos_index = Some(build_kmer_pos_index(&index));
    Ok(index)
}

pub fn build_bifrost_index_with_positions(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, load_positional_info, 1)
}

pub fn build_bifrost_index_threaded(path: &Path, threads: usize) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, false, threads)
}

pub fn build_bifrost_index_with_kmer_pos_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    let mut index =
        build_bifrost_index_with_positions_threaded(path, load_positional_info, threads)?;
    index.kmer_pos_index = Some(build_kmer_pos_index(&index));
    Ok(index)
}

pub fn build_bifrost_index_with_positions_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    let kmer_bytes = detect_kmer_bytes(path)?;
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;

    let file_format_version = read_u64_le(&mut reader)?;
    let k = read_i32_le(&mut reader)?;
    let g = read_i32_le(&mut reader)?;
    let _ = file_format_version;

    if k <= 0 || g <= 0 {
        return Err(Error::InvalidFormat("invalid k/g".into()));
    }
    let k = k as usize;
    let g = g as usize;
    if k > 32 || g > 32 {
        return Err(Error::UnsupportedFeature(
            "k/g > 32 not supported in Bifrost path".into(),
        ));
    }

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut unitigs = Vec::with_capacity(v_unitigs_sz);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        unitigs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut km_unitigs = Vec::with_capacity(km_unitigs_sz);
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        km_unitigs.push(decode_kmer(&buf, k));
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    let mut h_kmers = Vec::with_capacity(h_kmers_ccov_sz);
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        h_kmers.push(decode_kmer(&buf, k));
    }
    let h_kmer_map = if h_kmers.is_empty() {
        None
    } else {
        let mut map = HashMap::with_capacity(h_kmers.len());
        let base = unitigs.len() + km_unitigs.len();
        for (idx, seq) in h_kmers.iter().enumerate() {
            if let Some((fwd, rev)) = encode_kmer_pair(seq) {
                let code = fwd.min(rev);
                map.insert(code, base + idx);
            }
        }
        Some(map)
    };

    let _bfg_header = read_u64_le(&mut reader)?;
    let _checksum = read_u64_le(&mut reader)?;
    let _v_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _km_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _h_kmers_sz2 = read_u64_le(&mut reader)?;
    let hmap_min_unitigs_sz = read_u64_le(&mut reader)? as usize;

    let nb_bmp_unitigs = read_u64_le(&mut reader)? as usize;
    let mut bmp_unitigs = Vec::with_capacity(nb_bmp_unitigs);
    for _ in 0..nb_bmp_unitigs {
        bmp_unitigs.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_bmp_km = read_u64_le(&mut reader)? as usize;
    let mut bmp_km = Vec::with_capacity(nb_bmp_km);
    for _ in 0..nb_bmp_km {
        bmp_km.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_special_minz = read_u64_le(&mut reader)? as usize;
    let nb_bmp_special = (nb_special_minz >> 32) + 1;
    let mut abundant = HashSet::new();
    let mut overcrowded = HashSet::new();
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            abundant.insert(id as usize);
        }
    }
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            overcrowded.insert(id as usize);
        }
    }

    let mut special_minz = Vec::with_capacity(nb_special_minz);
    for i in 0..nb_special_minz {
        let minz = read_u64_le(&mut reader)?;
        let is_abundant = abundant.contains(&i);
        let is_overcrowded = overcrowded.contains(&i);
        let pos_id = if is_abundant || is_overcrowded {
            let mut id = 0xffff_ffff_0000_0000 | if is_overcrowded { 0x8000_0000 } else { 0 };
            if is_abundant {
                let count = read_u32_le(&mut reader)? as u64;
                id |= count;
            }
            id
        } else {
            read_u64_le(&mut reader)?
        };
        special_minz.push((minz, pos_id));
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;

    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mphf = BooPhf::load(&mut mphf_cur)?;

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let mut dlist = if dlist_size > 0 {
        HashSet::with_capacity(dlist_size as usize)
    } else {
        HashSet::new()
    };
    for _ in 0..dlist_size {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        let seq = decode_kmer(&buf, k);
        if let Some((fwd, rev)) = encode_kmer_pair(&seq) {
            dlist.insert(fwd.min(rev));
        }
    }
    let dlist = if dlist.is_empty() { None } else { Some(dlist) };

    let node_count = read_u64_le(&mut reader)? as usize;
    let mut ec_blocks = Vec::with_capacity(node_count);
    for _ in 0..node_count {
        reader.seek(SeekFrom::Current(k as i64))?;
        let node_size = read_u32_le(&mut reader)? as usize;
        let mut buf = vec![0u8; node_size];
        reader.read_exact(&mut buf)?;
        let mut cur = std::io::Cursor::new(buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let blocks = if load_positional_info {
            crate::index::read_block_array_blocks_with_positions(&mut cur)?
        } else {
            read_block_array_blocks(&mut cur)?
        };
        ec_blocks.push(blocks);
    }

    let load_threads = threads.max(1);
    let rayon_pool = if load_threads > 1 {
        Some(
            ThreadPoolBuilder::new()
                .num_threads(load_threads)
                .build()
                .map_err(|e| Error::InvalidFormat(format!("failed to build rayon pool: {e}")))?,
        )
    } else {
        None
    };
    let mut minz_counts = vec![0usize; hmap_min_unitigs_sz];
    populate_unitig_minimizer_counts(
        &mut minz_counts,
        &bmp_unitigs,
        &unitigs,
        g,
        &mphf,
        rayon_pool.as_ref(),
    );

    let km_glen = k.saturating_sub(g) + 1;
    let kmer_load_ctx = KmerLoadContext {
        km_unitigs: &km_unitigs,
        g,
        km_glen,
        mphf: &mphf,
        pool: rayon_pool.as_ref(),
    };
    populate_kmer_minimizer_counts(&mut minz_counts, &bmp_km, &kmer_load_ctx);

    for (minz, _pos_id) in special_minz.iter().copied() {
        let bytes = minz.to_le_bytes();
        if let Some(idx) = mphf.lookup(&bytes) {
            minz_counts[idx as usize] += 1;
        }
    }

    let mut minz_positions = MinzPositionIndex::from_counts(minz_counts);
    let mut minz_cursors = minz_positions.cursors();
    populate_unitig_minimizer_storage(
        &mut minz_positions,
        &mut minz_cursors,
        &bmp_unitigs,
        &unitigs,
        g,
        &mphf,
        rayon_pool.as_ref(),
    );
    populate_kmer_minimizer_storage(
        &mut minz_positions,
        &mut minz_cursors,
        &bmp_km,
        &kmer_load_ctx,
    );

    for (minz, pos_id) in special_minz {
        let bytes = minz.to_le_bytes();
        if let Some(idx) = mphf.lookup(&bytes) {
            minz_positions.push_at(&mut minz_cursors, idx as usize, pos_id);
        }
    }

    let onlist = read_onlist_with_lengths(&mut reader)?;

    Ok(BifrostIndex {
        k,
        g,
        unitigs,
        km_unitigs,
        h_kmers,
        h_kmer_map,
        dlist,
        kmer_pos_index: None,
        kmer_index: None,
        ec_blocks,
        minz_positions,
        mphf,
        onlist: onlist.onlist,
        transcript_lengths: onlist.lengths,
        transcript_names: onlist.names,
        shade_to_color: onlist.shade_to_color,
        shade_sequences: onlist.shade_sequences,
        use_shade: onlist.use_shade,
    })
}

const LOAD_BATCH_SIZE: usize = 1 << 20;

struct KmerLoadContext<'a> {
    km_unitigs: &'a [Vec<u8>],
    g: usize,
    km_glen: usize,
    mphf: &'a BooPhf,
    pool: Option<&'a ThreadPool>,
}

fn batch_capacity(len: usize) -> usize {
    len.clamp(1, LOAD_BATCH_SIZE)
}

fn populate_unitig_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_unitigs: &[Vec<u32>],
    unitigs: &[Vec<u8>],
    g: usize,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
) {
    if bmp_unitigs.is_empty() || unitigs.is_empty() {
        return;
    }

    let mut unitig_id = 0usize;
    let mut total_unitig_len = 0u64;
    let mut current_unitig_len = unitigs[0].len().saturating_sub(g) as u64 + 1;

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        let mut batch = Vec::with_capacity(batch_capacity(bitmap.len()));
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            while pos >= total_unitig_len + current_unitig_len {
                unitig_id += 1;
                total_unitig_len += current_unitig_len;
                current_unitig_len = unitigs[unitig_id].len().saturating_sub(g) as u64 + 1;
            }
            let rel = pos - total_unitig_len;
            batch.push(((unitig_id as u64) << 32) | rel);
            if batch.len() >= LOAD_BATCH_SIZE {
                flush_unitig_count_batch(minz_counts, &batch, unitigs, g, mphf, pool);
                batch.clear();
            }
        }
        if !batch.is_empty() {
            flush_unitig_count_batch(minz_counts, &batch, unitigs, g, mphf, pool);
        }
    }
}

fn flush_unitig_count_batch(
    minz_counts: &mut [usize],
    batch: &[u64],
    unitigs: &[Vec<u8>],
    g: usize,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
) {
    let pairs = parallel_pairs(batch, pool, |pos_id| {
        let pos_id = *pos_id;
        let unitig_id = (pos_id >> 32) as usize;
        let rel = (pos_id & 0xffff_ffff) as usize;
        let seq = &unitigs[unitig_id];
        encode_minimizer_rep(&seq[rel..rel + g])
            .and_then(|bytes| mphf.lookup(&bytes))
            .map(|idx| (idx as usize, pos_id))
    });
    for (idx, _) in pairs {
        minz_counts[idx] += 1;
    }
}

fn populate_unitig_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_unitigs: &[Vec<u32>],
    unitigs: &[Vec<u8>],
    g: usize,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
) {
    if bmp_unitigs.is_empty() || unitigs.is_empty() {
        return;
    }

    let mut unitig_id = 0usize;
    let mut total_unitig_len = 0u64;
    let mut current_unitig_len = unitigs[0].len().saturating_sub(g) as u64 + 1;

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        let mut batch = Vec::with_capacity(batch_capacity(bitmap.len()));
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            while pos >= total_unitig_len + current_unitig_len {
                unitig_id += 1;
                total_unitig_len += current_unitig_len;
                current_unitig_len = unitigs[unitig_id].len().saturating_sub(g) as u64 + 1;
            }
            let rel = pos - total_unitig_len;
            batch.push(((unitig_id as u64) << 32) | rel);
            if batch.len() >= LOAD_BATCH_SIZE {
                flush_unitig_storage_batch(
                    minz_positions,
                    minz_cursors,
                    &batch,
                    unitigs,
                    g,
                    mphf,
                    pool,
                );
                batch.clear();
            }
        }
        if !batch.is_empty() {
            flush_unitig_storage_batch(
                minz_positions,
                minz_cursors,
                &batch,
                unitigs,
                g,
                mphf,
                pool,
            );
        }
    }
}

fn flush_unitig_storage_batch(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    batch: &[u64],
    unitigs: &[Vec<u8>],
    g: usize,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
) {
    let pairs = parallel_pairs(batch, pool, |pos_id| {
        let pos_id = *pos_id;
        let unitig_id = (pos_id >> 32) as usize;
        let rel = (pos_id & 0xffff_ffff) as usize;
        let seq = &unitigs[unitig_id];
        encode_minimizer_rep(&seq[rel..rel + g])
            .and_then(|bytes| mphf.lookup(&bytes))
            .map(|idx| (idx as usize, pos_id))
    });
    for (idx, pos_id) in pairs {
        minz_positions.push_at(minz_cursors, idx, pos_id);
    }
}

fn populate_kmer_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_km: &[Vec<u32>],
    ctx: &KmerLoadContext<'_>,
) {
    for (i, bitmap) in bmp_km.iter().enumerate() {
        let base = (i as u64) << 32;
        let mut batch = Vec::with_capacity(batch_capacity(bitmap.len()));
        for &pos_bmp in bitmap {
            batch.push(base + pos_bmp as u64);
            if batch.len() >= LOAD_BATCH_SIZE {
                flush_kmer_count_batch(minz_counts, &batch, ctx);
                batch.clear();
            }
        }
        if !batch.is_empty() {
            flush_kmer_count_batch(minz_counts, &batch, ctx);
        }
    }
}

fn flush_kmer_count_batch(minz_counts: &mut [usize], batch: &[u64], ctx: &KmerLoadContext<'_>) {
    let pairs = parallel_pairs(batch, ctx.pool, |pos| {
        let pos = *pos;
        let km_id = (pos / ctx.km_glen as u64) as usize;
        let km_pos = (pos % ctx.km_glen as u64) as usize;
        if km_id >= ctx.km_unitigs.len() || km_pos + ctx.g > ctx.km_unitigs[km_id].len() {
            return None;
        }
        let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | (km_pos as u64);
        encode_minimizer_rep(&ctx.km_unitigs[km_id][km_pos..km_pos + ctx.g])
            .and_then(|bytes| ctx.mphf.lookup(&bytes))
            .map(|idx| (idx as usize, pos_id))
    });
    for (idx, _) in pairs {
        minz_counts[idx] += 1;
    }
}

fn populate_kmer_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_km: &[Vec<u32>],
    ctx: &KmerLoadContext<'_>,
) {
    for (i, bitmap) in bmp_km.iter().enumerate() {
        let base = (i as u64) << 32;
        let mut batch = Vec::with_capacity(batch_capacity(bitmap.len()));
        for &pos_bmp in bitmap {
            batch.push(base + pos_bmp as u64);
            if batch.len() >= LOAD_BATCH_SIZE {
                flush_kmer_storage_batch(minz_positions, minz_cursors, &batch, ctx);
                batch.clear();
            }
        }
        if !batch.is_empty() {
            flush_kmer_storage_batch(minz_positions, minz_cursors, &batch, ctx);
        }
    }
}

fn flush_kmer_storage_batch(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    batch: &[u64],
    ctx: &KmerLoadContext<'_>,
) {
    let pairs = parallel_pairs(batch, ctx.pool, |pos| {
        let pos = *pos;
        let km_id = (pos / ctx.km_glen as u64) as usize;
        let km_pos = (pos % ctx.km_glen as u64) as usize;
        if km_id >= ctx.km_unitigs.len() || km_pos + ctx.g > ctx.km_unitigs[km_id].len() {
            return None;
        }
        let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | (km_pos as u64);
        encode_minimizer_rep(&ctx.km_unitigs[km_id][km_pos..km_pos + ctx.g])
            .and_then(|bytes| ctx.mphf.lookup(&bytes))
            .map(|idx| (idx as usize, pos_id))
    });
    for (idx, pos_id) in pairs {
        minz_positions.push_at(minz_cursors, idx, pos_id);
    }
}

fn parallel_pairs<T, F>(items: &[T], pool: Option<&ThreadPool>, map: F) -> Vec<(usize, u64)>
where
    T: Sync + Send,
    F: Fn(&T) -> Option<(usize, u64)> + Sync + Send,
{
    if pool.is_none() || items.len() < 4096 {
        return items.iter().filter_map(map).collect();
    }

    pool.expect("checked above").install(|| {
        items
            .par_iter()
            .with_min_len(4096)
            .filter_map(map)
            .collect()
    })
}

fn build_kmer_pos_index(index: &BifrostIndex) -> KmerPosIndex {
    let mut map: KmerPosIndex = HashMap::new();
    let k = index.k;
    for (unitig_id, seq) in index.unitigs.iter().enumerate() {
        if seq.len() < k {
            continue;
        }
        for pos in 0..=seq.len() - k {
            if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + k]) {
                map.entry(fwd).or_default().push((unitig_id, pos, false));
                if rev != fwd {
                    map.entry(rev).or_default().push((unitig_id, pos, true));
                }
            }
        }
    }
    let base_km = index.unitigs.len();
    for (km_id, seq) in index.km_unitigs.iter().enumerate() {
        if seq.len() != k {
            continue;
        }
        if let Some((fwd, rev)) = encode_kmer_pair(seq) {
            let unitig_id = base_km + km_id;
            map.entry(fwd).or_default().push((unitig_id, 0, false));
            if rev != fwd {
                map.entry(rev).or_default().push((unitig_id, 0, true));
            }
        }
    }
    let base_h = index.unitigs.len() + index.km_unitigs.len();
    for (hk_id, seq) in index.h_kmers.iter().enumerate() {
        if seq.len() != k {
            continue;
        }
        if let Some((fwd, rev)) = encode_kmer_pair(seq) {
            let unitig_id = base_h + hk_id;
            map.entry(fwd).or_default().push((unitig_id, 0, false));
            if rev != fwd {
                map.entry(rev).or_default().push((unitig_id, 0, true));
            }
        }
    }
    map
}

fn detect_kmer_bytes(path: &Path) -> Result<usize> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    if dbg_size == 0 {
        return Err(Error::InvalidFormat("missing graph section".into()));
    }

    let dbg_start = reader.stream_position()?;
    let meta = parse_graph_section(&mut reader, dbg_start, dbg_size, &KMER_BYTES_CANDIDATES)?;
    Ok(meta.kmer_bytes)
}

pub(super) fn read_unitig_sequences<R: Read + Seek>(
    reader: &mut R,
    dbg_start: u64,
    dbg_size: u64,
    k: u32,
    kmer_bytes: usize,
) -> Result<Vec<Vec<u8>>> {
    reader.seek(SeekFrom::Start(dbg_start))?;
    let _file_format_version = read_u64_le(reader)?;
    let _k = read_i32_le(reader)?;
    let _g = read_i32_le(reader)?;

    let v_unitigs_sz = read_u64_le(reader)?;
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(v_unitigs_sz as usize);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        seqs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(reader)?;
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    let h_kmers_ccov_sz = read_u64_le(reader)?;
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    Ok(seqs)
}

fn decode_compressed_sequence(data: &[u8], len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let byte = data[i >> 2];
        let shift = (i & 0x3) << 1;
        let v = (byte >> shift) & 0x03;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}

fn decode_kmer(bytes: &[u8], k: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(k);
    let mut longs: Vec<u64> = Vec::with_capacity(bytes.len() / 8);
    for chunk in bytes.chunks_exact(8) {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(chunk);
        longs.push(u64::from_le_bytes(arr));
    }
    for i in 0..k {
        let idx = i >> 5;
        let shift = 62 - ((i & 0x1F) << 1);
        let v = (longs[idx] >> shift) & 0x3;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}
