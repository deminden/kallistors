use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::ops::Index;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::index::bifrost::{BooPhf, read_bitcontainer_values};
use crate::index::{parse_graph_section, read_block_array_blocks};
use crate::timing::{self, Stage};
use crate::{Error, Result};

use super::ec::encode_kmer_pair;
use super::io::{read_i32_le, read_u32_le, read_u64_le};
use super::{KMER_BYTES_CANDIDATES, KmerEcIndex, read_onlist_with_lengths};

pub type KmerPosIndex = HashMap<u64, Vec<(usize, usize, bool)>>;

const LOAD_BATCH_POSITIONS: usize = 1 << 20;
const LOAD_MIN_PAR_CHUNK: usize = 1 << 14;

pub(crate) struct EncodedMinimizerStore {
    offsets: Vec<usize>,
    values: Vec<[u8; 8]>,
}

impl EncodedMinimizerStore {
    fn from_sequences(seqs: &[Vec<u8>], g: usize, pool: Option<&ThreadPool>) -> Self {
        let mut offsets = Vec::with_capacity(seqs.len() + 1);
        offsets.push(0);
        let mut total = 0usize;
        for seq in seqs {
            total += minimizer_count(seq.len(), g);
            offsets.push(total);
        }

        let build_sequential = || {
            let mut values = vec![[0u8; 8]; total];
            for (seq_idx, seq) in seqs.iter().enumerate() {
                let start = offsets[seq_idx];
                let end = offsets[seq_idx + 1];
                if start == end {
                    continue;
                }
                compute_minimizer_reps_into(seq, g, &mut values[start..end]);
            }
            values
        };

        let values = if let Some(pool) = pool {
            let mut values = vec![[0u8; 8]; total];
            let base_ptr = values.as_mut_ptr() as usize;
            pool.install(|| {
                (0..seqs.len())
                    .into_par_iter()
                    .with_min_len(2048)
                    .for_each(|seq_idx| {
                        let start = offsets[seq_idx];
                        let end = offsets[seq_idx + 1];
                        if start == end {
                            return;
                        }
                        let out = unsafe {
                            std::slice::from_raw_parts_mut(
                                (base_ptr as *mut [u8; 8]).add(start),
                                end - start,
                            )
                        };
                        compute_minimizer_reps_into(&seqs[seq_idx], g, out);
                    });
            });
            values
        } else {
            build_sequential()
        };

        Self { offsets, values }
    }

    #[inline]
    fn len(&self, seq_idx: usize) -> usize {
        self.offsets[seq_idx + 1] - self.offsets[seq_idx]
    }

    #[inline]
    fn get(&self, seq_idx: usize, pos: usize) -> Option<&[u8; 8]> {
        let start = *self.offsets.get(seq_idx)?;
        let end = *self.offsets.get(seq_idx + 1)?;
        self.values.get(start..end)?.get(pos)
    }
}

pub(crate) struct EncodedBaseStore {
    offsets: Vec<usize>,
    values: Vec<u8>,
}

impl EncodedBaseStore {
    fn from_sequences(seqs: &[Vec<u8>], pool: Option<&ThreadPool>) -> Self {
        let mut offsets = Vec::with_capacity(seqs.len() + 1);
        offsets.push(0);
        let mut total = 0usize;
        for seq in seqs {
            total += seq.len();
            offsets.push(total);
        }

        let build_sequential = || {
            let mut values = vec![0u8; total];
            for (seq_idx, seq) in seqs.iter().enumerate() {
                let start = offsets[seq_idx];
                let end = offsets[seq_idx + 1];
                if start == end {
                    continue;
                }
                encode_bases_into(seq, &mut values[start..end]);
            }
            values
        };

        let values = if let Some(pool) = pool {
            let mut values = vec![0u8; total];
            let base_ptr = values.as_mut_ptr() as usize;
            pool.install(|| {
                (0..seqs.len())
                    .into_par_iter()
                    .with_min_len(2048)
                    .for_each(|seq_idx| {
                        let start = offsets[seq_idx];
                        let end = offsets[seq_idx + 1];
                        if start == end {
                            return;
                        }
                        // SAFETY: each sequence owns a disjoint prefix-summed slice.
                        let out = unsafe {
                            std::slice::from_raw_parts_mut(
                                (base_ptr as *mut u8).add(start),
                                end - start,
                            )
                        };
                        encode_bases_into(&seqs[seq_idx], out);
                    });
            });
            values
        } else {
            build_sequential()
        };

        Self { offsets, values }
    }

    #[inline]
    pub(crate) fn encode_kmer_at(&self, seq_idx: usize, start: usize, k: usize) -> Option<u64> {
        if k > 32 {
            return None;
        }
        let seq_start = *self.offsets.get(seq_idx)?;
        let seq_end = *self.offsets.get(seq_idx + 1)?;
        let seq_len = seq_end.checked_sub(seq_start)?;
        if start.checked_add(k)? > seq_len {
            return None;
        }
        let mut out = 0u64;
        for &base in &self.values[seq_start + start..seq_start + start + k] {
            out = (out << 2) | u64::from(base);
        }
        Some(out)
    }
}

pub(crate) struct EncodedUnitigs {
    pub(crate) unitig_minimizers: EncodedMinimizerStore,
    pub(crate) km_minimizers: EncodedMinimizerStore,
    pub(crate) unitig_bases: EncodedBaseStore,
    pub(crate) km_bases: EncodedBaseStore,
}

pub(crate) struct FlatEcIndex {
    unitig_block_offsets: Vec<usize>,
    block_lbs: Vec<u32>,
    block_ubs: Vec<u32>,
    block_ec_offsets: Vec<usize>,
    ec_values: Vec<u32>,
}

impl FlatEcIndex {
    fn new(unitig_capacity: usize) -> Self {
        let mut unitig_block_offsets = Vec::with_capacity(unitig_capacity + 1);
        unitig_block_offsets.push(0);
        let block_ec_offsets = vec![0];
        Self {
            unitig_block_offsets,
            block_lbs: Vec::new(),
            block_ubs: Vec::new(),
            block_ec_offsets,
            ec_values: Vec::new(),
        }
    }

    fn push_unitig_blocks(&mut self, blocks: &[crate::index::EcBlock]) {
        for block in blocks {
            self.block_lbs.push(block.lb);
            self.block_ubs.push(block.ub);
            self.ec_values.extend_from_slice(&block.ec);
            self.block_ec_offsets.push(self.ec_values.len());
        }
        self.unitig_block_offsets.push(self.block_lbs.len());
    }

    pub(crate) fn block_index_for_position(&self, unitig_id: usize, pos: usize) -> Option<usize> {
        let start = *self.unitig_block_offsets.get(unitig_id)?;
        let end = *self.unitig_block_offsets.get(unitig_id + 1)?;
        if start == end {
            return None;
        }

        let pos = pos as u32;
        let mut lo = start;
        let mut hi = end;
        while lo < hi {
            let mid = (lo + hi) / 2;
            let lb = self.block_lbs[mid];
            let ub = self.block_ubs[mid];
            if pos < lb {
                hi = mid;
            } else if pos >= ub {
                lo = mid + 1;
            } else {
                return Some(mid - start);
            }
        }
        None
    }

    pub(crate) fn ec(&self, unitig_id: usize, block_idx: usize) -> &[u32] {
        let global_idx = self.unitig_block_offsets[unitig_id] + block_idx;
        &self.ec_values[self.block_ec_offsets[global_idx]..self.block_ec_offsets[global_idx + 1]]
    }

    pub(crate) fn block_bounds(&self, unitig_id: usize, block_idx: usize) -> Option<(u32, u32)> {
        let global_idx = self.unitig_block_offsets.get(unitig_id)? + block_idx;
        Some((
            *self.block_lbs.get(global_idx)?,
            *self.block_ubs.get(global_idx)?,
        ))
    }
}

pub(crate) struct SelectiveFallbackKmerIndex {
    entries: HashMap<u64, Vec<(usize, u64)>>,
}

impl SelectiveFallbackKmerIndex {
    fn build(unitigs_len: usize, km_unitigs: &[Vec<u8>], h_kmers: &[Vec<u8>]) -> Self {
        let mut entries = HashMap::<u64, Vec<(usize, u64)>>::new();
        for (km_id, seq) in km_unitigs.iter().enumerate() {
            if let Some((fwd, rev)) = encode_kmer_pair(seq) {
                entries
                    .entry(fwd.min(rev))
                    .or_default()
                    .push((unitigs_len + km_id, fwd));
            }
        }
        let base_h = unitigs_len + km_unitigs.len();
        for (hk_id, seq) in h_kmers.iter().enumerate() {
            if let Some((fwd, rev)) = encode_kmer_pair(seq) {
                entries
                    .entry(fwd.min(rev))
                    .or_default()
                    .push((base_h + hk_id, fwd));
            }
        }
        Self { entries }
    }

    pub(crate) fn lookup(&self, code: u64) -> Option<&[(usize, u64)]> {
        self.entries.get(&code).map(Vec::as_slice)
    }
}

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
    #[allow(dead_code)]
    pub(crate) encoded_unitigs: EncodedUnitigs,
    pub kmer_pos_index: Option<KmerPosIndex>,
    pub kmer_index: Option<KmerEcIndex>,
    pub(crate) ec_blocks: Vec<Vec<crate::index::EcBlock>>,
    pub(crate) flat_ec: FlatEcIndex,
    pub minz_positions: MinzPositionIndex,
    pub mphf: BooPhf,
    pub(crate) selective_fallback_kmers: SelectiveFallbackKmerIndex,
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

pub fn build_bifrost_index_quant_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(path, load_positional_info, threads, false)
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
    build_bifrost_index_with_positions_threaded_impl(path, load_positional_info, threads, true)
}

fn build_bifrost_index_with_positions_threaded_impl(
    path: &Path,
    load_positional_info: bool,
    _threads: usize,
    store_ec_blocks: bool,
) -> Result<BifrostIndex> {
    let _header_timing = timing::scoped(Stage::IndexHeaderParse);
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;
    let graph_meta = parse_graph_section(&mut reader, dbg_start, dbg_size, &KMER_BYTES_CANDIDATES)?;
    let kmer_bytes = graph_meta.kmer_bytes;
    reader.seek(SeekFrom::Start(dbg_start))?;
    drop(_header_timing);

    let _graph_timing = timing::scoped(Stage::GraphDecode);
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
    let mut ec_blocks = if store_ec_blocks {
        Vec::with_capacity(node_count)
    } else {
        Vec::new()
    };
    let mut flat_ec = FlatEcIndex::new(node_count);
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
        flat_ec.push_unitig_blocks(&blocks);
        if store_ec_blocks {
            ec_blocks.push(blocks);
        }
    }
    let selective_fallback_kmers =
        SelectiveFallbackKmerIndex::build(unitigs.len(), &km_unitigs, &h_kmers);
    let load_threads = _threads.max(1);
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
    let encoded_unitigs = precompute_encoded_unitigs(&unitigs, &km_unitigs, g, rayon_pool.as_ref());

    let mut minz_counts = vec![0usize; hmap_min_unitigs_sz];
    drop(_graph_timing);
    let _count_timing = timing::scoped(Stage::MinimizerCountPass);
    populate_unitig_minimizer_counts(
        &mut minz_counts,
        &bmp_unitigs,
        &encoded_unitigs,
        &mphf,
        rayon_pool.as_ref(),
    );

    let km_glen = k.saturating_sub(g) + 1;
    let kmer_load_ctx = KmerLoadContext {
        km_minimizers: &encoded_unitigs.km_minimizers,
        km_glen,
        mphf: &mphf,
    };
    populate_kmer_minimizer_counts(
        &mut minz_counts,
        &bmp_km,
        &kmer_load_ctx,
        rayon_pool.as_ref(),
    );

    for (minz, _pos_id) in special_minz.iter().copied() {
        let bytes = minz.to_le_bytes();
        if let Some(idx) = mphf.lookup(&bytes) {
            minz_counts[idx as usize] += 1;
        }
    }

    let mut minz_positions = MinzPositionIndex::from_counts(minz_counts);
    let mut minz_cursors = minz_positions.cursors();
    drop(_count_timing);
    let _fill_timing = timing::scoped(Stage::MinimizerFillPass);
    populate_unitig_minimizer_storage(
        &mut minz_positions,
        &mut minz_cursors,
        &bmp_unitigs,
        &encoded_unitigs,
        &mphf,
        rayon_pool.as_ref(),
    );
    populate_kmer_minimizer_storage(
        &mut minz_positions,
        &mut minz_cursors,
        &bmp_km,
        &kmer_load_ctx,
        rayon_pool.as_ref(),
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
        encoded_unitigs,
        kmer_pos_index: None,
        kmer_index: None,
        ec_blocks,
        flat_ec,
        minz_positions,
        mphf,
        selective_fallback_kmers,
        onlist: onlist.onlist,
        transcript_lengths: onlist.lengths,
        transcript_names: onlist.names,
        shade_to_color: onlist.shade_to_color,
        shade_sequences: onlist.shade_sequences,
        use_shade: onlist.use_shade,
    })
}

struct KmerLoadContext<'a> {
    km_minimizers: &'a EncodedMinimizerStore,
    km_glen: usize,
    mphf: &'a BooPhf,
}

fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn encode_word_pair(seq: &[u8]) -> Option<(u64, u64)> {
    let mut fwd = 0u64;
    for (i, &base) in seq.iter().enumerate() {
        let code = base_code(base)?;
        fwd |= code << (62 - ((i & 0x1f) << 1));
    }

    let mut rev = 0u64;
    for (i, &base) in seq.iter().rev().enumerate() {
        let code = base_code(base)?;
        rev |= (3 - code) << (62 - ((i & 0x1f) << 1));
    }
    Some((fwd, rev))
}

#[inline]
fn minimizer_count(seq_len: usize, g: usize) -> usize {
    if g == 0 || seq_len < g {
        0
    } else {
        seq_len - g + 1
    }
}

fn compute_minimizer_reps_into(seq: &[u8], g: usize, out: &mut [[u8; 8]]) {
    if out.is_empty() {
        return;
    }

    debug_assert_eq!(out.len(), minimizer_count(seq.len(), g));

    let Some((mut fwd, mut rev)) = encode_word_pair(&seq[..g]) else {
        return;
    };
    let tail_shift = 64usize.saturating_sub(g * 2);
    let high_mask = if g == 32 {
        u64::MAX
    } else {
        (!0u64) << tail_shift
    };

    out[0] = fwd.min(rev).to_le_bytes();

    for (dst, &base) in out.iter_mut().skip(1).zip(&seq[g..]) {
        let Some(code) = base_code(base) else {
            return;
        };
        let comp = 3 - code;
        fwd = ((fwd << 2) & high_mask) | (code << tail_shift);
        rev = ((rev >> 2) & high_mask) | (comp << 62);
        *dst = fwd.min(rev).to_le_bytes();
    }
}

fn encode_bases_into(seq: &[u8], out: &mut [u8]) {
    debug_assert_eq!(seq.len(), out.len());
    for (src, dst) in seq.iter().zip(out.iter_mut()) {
        *dst = match *src {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        };
    }
}

fn precompute_encoded_unitigs(
    unitigs: &[Vec<u8>],
    km_unitigs: &[Vec<u8>],
    g: usize,
    pool: Option<&ThreadPool>,
) -> EncodedUnitigs {
    EncodedUnitigs {
        unitig_minimizers: EncodedMinimizerStore::from_sequences(unitigs, g, pool),
        km_minimizers: EncodedMinimizerStore::from_sequences(km_unitigs, g, pool),
        unitig_bases: EncodedBaseStore::from_sequences(unitigs, pool),
        km_bases: EncodedBaseStore::from_sequences(km_unitigs, pool),
    }
}

#[derive(Clone, Copy)]
struct UnitigBlockState {
    unitig_id: usize,
    total_unitig_len: u64,
    current_unitig_len: u64,
}

fn advance_unitig_state(state: &mut UnitigBlockState, encoded_unitigs: &EncodedUnitigs, pos: u64) {
    while pos >= state.total_unitig_len + state.current_unitig_len {
        state.unitig_id += 1;
        state.total_unitig_len += state.current_unitig_len;
        state.current_unitig_len = encoded_unitigs.unitig_minimizers.len(state.unitig_id) as u64;
    }
}

fn precompute_unitig_block_states(
    bmp_unitigs: &[Vec<u32>],
    encoded_unitigs: &EncodedUnitigs,
) -> Vec<UnitigBlockState> {
    let mut out = Vec::with_capacity(bmp_unitigs.len());
    if encoded_unitigs.unitig_minimizers.offsets.len() <= 1 {
        return out;
    }

    let mut state = UnitigBlockState {
        unitig_id: 0,
        total_unitig_len: 0,
        current_unitig_len: encoded_unitigs.unitig_minimizers.len(0) as u64,
    };
    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        advance_unitig_state(&mut state, encoded_unitigs, base);
        out.push(state);
        if let Some(&last_pos) = bitmap.last() {
            advance_unitig_state(&mut state, encoded_unitigs, base + last_pos as u64);
        }
    }
    out
}

fn populate_unitig_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_unitigs: &[Vec<u32>],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
) {
    if bmp_unitigs.is_empty() || encoded_unitigs.unitig_minimizers.offsets.len() <= 1 {
        return;
    }

    let Some(pool) = pool else {
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizers.len(0) as u64,
        };
        for (i, bitmap) in bmp_unitigs.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizers.get(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(bytes)
                {
                    minz_counts[idx as usize] += 1;
                }
            }
        }
        return;
    };

    let states = precompute_unitig_block_states(bmp_unitigs, encoded_unitigs);
    let mut batch_pos_ids = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut flush_batch = |batch_pos_ids: &mut Vec<u64>| {
        if batch_pos_ids.is_empty() {
            return;
        }
        let grouped = collect_unitig_batch_counts(batch_pos_ids, encoded_unitigs, mphf, pool);
        merge_count_batches(minz_counts, grouped);
        batch_pos_ids.clear();
    };

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let mut state = states[i];
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            let rel = (pos - state.total_unitig_len) as usize;
            batch_pos_ids.push(((state.unitig_id as u64) << 32) | rel as u64);
            if batch_pos_ids.len() >= LOAD_BATCH_POSITIONS {
                flush_batch(&mut batch_pos_ids);
            }
        }
    }
    flush_batch(&mut batch_pos_ids);
}

fn populate_unitig_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_unitigs: &[Vec<u32>],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
) {
    if bmp_unitigs.is_empty() || encoded_unitigs.unitig_minimizers.offsets.len() <= 1 {
        return;
    }

    let Some(pool) = pool else {
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizers.len(0) as u64,
        };
        for (i, bitmap) in bmp_unitigs.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizers.get(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(bytes)
                {
                    let pos_id = ((state.unitig_id as u64) << 32) | rel as u64;
                    minz_positions.push_at(minz_cursors, idx as usize, pos_id);
                }
            }
        }
        return;
    };

    let states = precompute_unitig_block_states(bmp_unitigs, encoded_unitigs);
    let cursors = minz_cursors
        .iter()
        .copied()
        .map(AtomicUsize::new)
        .collect::<Vec<_>>();
    let values = minz_positions
        .values
        .iter()
        .copied()
        .map(AtomicU64::new)
        .collect::<Vec<_>>();
    pool.install(|| {
        bmp_unitigs.par_iter().enumerate().for_each(|(i, bitmap)| {
            let mut state = states[i];
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizers.get(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(bytes)
                {
                    let pos_id = ((state.unitig_id as u64) << 32) | rel as u64;
                    let slot = cursors[idx as usize].fetch_add(1, Ordering::Relaxed);
                    values[slot].store(pos_id, Ordering::Relaxed);
                }
            }
        });
    });
    for (dst, src) in minz_positions.values.iter_mut().zip(values) {
        *dst = src.into_inner();
    }
    for (dst, src) in minz_cursors.iter_mut().zip(cursors) {
        *dst = src.into_inner();
    }
}

fn populate_kmer_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_km: &[Vec<u32>],
    ctx: &KmerLoadContext<'_>,
    pool: Option<&ThreadPool>,
) {
    let Some(pool) = pool else {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_minimizers.get(km_id, km_pos)
                    && let Some(idx) = ctx.mphf.lookup(bytes)
                {
                    minz_counts[idx as usize] += 1;
                }
            }
        }
        return;
    };

    let mut batch_pos_ids = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut flush_batch = |batch_pos_ids: &mut Vec<u64>| {
        if batch_pos_ids.is_empty() {
            return;
        }
        let grouped = collect_kmer_batch_counts(batch_pos_ids, ctx, pool);
        merge_count_batches(minz_counts, grouped);
        batch_pos_ids.clear();
    };

    for (i, bitmap) in bmp_km.iter().enumerate() {
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            let km_id = (pos / ctx.km_glen as u64) as usize;
            let km_pos = (pos % ctx.km_glen as u64) as usize;
            batch_pos_ids.push(((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64);
            if batch_pos_ids.len() >= LOAD_BATCH_POSITIONS {
                flush_batch(&mut batch_pos_ids);
            }
        }
    }
    flush_batch(&mut batch_pos_ids);
}

fn populate_kmer_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_km: &[Vec<u32>],
    ctx: &KmerLoadContext<'_>,
    pool: Option<&ThreadPool>,
) {
    let Some(pool) = pool else {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_minimizers.get(km_id, km_pos)
                    && let Some(idx) = ctx.mphf.lookup(bytes)
                {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64;
                    minz_positions.push_at(minz_cursors, idx as usize, pos_id);
                }
            }
        }
        return;
    };

    let cursors = minz_cursors
        .iter()
        .copied()
        .map(AtomicUsize::new)
        .collect::<Vec<_>>();
    let values = minz_positions
        .values
        .iter()
        .copied()
        .map(AtomicU64::new)
        .collect::<Vec<_>>();
    pool.install(|| {
        bmp_km.par_iter().enumerate().for_each(|(i, bitmap)| {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_minimizers.get(km_id, km_pos)
                    && let Some(idx) = ctx.mphf.lookup(bytes)
                {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64;
                    let slot = cursors[idx as usize].fetch_add(1, Ordering::Relaxed);
                    values[slot].store(pos_id, Ordering::Relaxed);
                }
            }
        });
    });
    for (dst, src) in minz_positions.values.iter_mut().zip(values) {
        *dst = src.into_inner();
    }
    for (dst, src) in minz_cursors.iter_mut().zip(cursors) {
        *dst = src.into_inner();
    }
}

fn collect_unitig_batch_counts(
    pos_ids: &[u64],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: &ThreadPool,
) -> Vec<Vec<usize>> {
    let chunk_len =
        (pos_ids.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    pool.install(|| {
        pos_ids
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut indices = Vec::with_capacity(chunk.len());
                for &pos_id in chunk {
                    let unitig_id = (pos_id >> 32) as usize;
                    let rel = (pos_id & 0xffff_ffff) as usize;
                    if let Some(bytes) = encoded_unitigs.unitig_minimizers.get(unitig_id, rel)
                        && let Some(idx) = mphf.lookup(bytes)
                    {
                        indices.push(idx as usize);
                    }
                }
                indices.sort_unstable();
                indices
            })
            .collect()
    })
}

fn collect_kmer_batch_counts(
    pos_ids: &[u64],
    ctx: &KmerLoadContext<'_>,
    pool: &ThreadPool,
) -> Vec<Vec<usize>> {
    let chunk_len =
        (pos_ids.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    pool.install(|| {
        pos_ids
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut indices = Vec::with_capacity(chunk.len());
                for &pos_id in chunk {
                    let km_id = (pos_id >> 32) as usize;
                    let km_pos = (pos_id & 0x7fff_ffff) as usize;
                    if let Some(bytes) = ctx.km_minimizers.get(km_id, km_pos)
                        && let Some(idx) = ctx.mphf.lookup(bytes)
                    {
                        indices.push(idx as usize);
                    }
                }
                indices.sort_unstable();
                indices
            })
            .collect()
    })
}

fn merge_count_batches(minz_counts: &mut [usize], grouped: Vec<Vec<usize>>) {
    for indices in grouped {
        let mut iter = indices.into_iter();
        let Some(mut idx) = iter.next() else {
            continue;
        };
        let mut count = 1usize;
        for next_idx in iter {
            if next_idx == idx {
                count += 1;
            } else {
                minz_counts[idx] += count;
                idx = next_idx;
                count = 1;
            }
        }
        minz_counts[idx] += count;
    }
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
