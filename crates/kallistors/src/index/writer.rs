use super::workers::IndexWorkers;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use croaring_sys::{
    roaring_bitmap_add_many, roaring_bitmap_add_range, roaring_bitmap_create, roaring_bitmap_free,
    roaring_bitmap_portable_serialize, roaring_bitmap_portable_size_in_bytes,
    roaring_bitmap_run_optimize, roaring_bitmap_serialize, roaring_bitmap_size_in_bytes,
};

use super::bifrost::{ReducedDomain, minimizer_hash, wyhash};
use super::graph_build::{BuildGraph, NodeBlock, NodePayload, encode_kmer};
use crate::{Error, Result};

const INDEX_VERSION: u64 = 13;
const BFG_GRAPHBIN_FORMAT_HEADER: u64 = 0x7e21_5f3f;
const BFG_METABIN_FORMAT_HEADER: u64 = 0x267c_3d5d;
const BFG_FORMAT_VERSION: u64 = 1;
const FLAG_PTR_BITMAP: u64 = 0x3;
const SHIFT_MASK_BITS: u64 = 3;
const BOO_GAMMA: f64 = 2.0;
const BOO_LEVELS: usize = 25;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct IndexWriteReport {
    pub(super) total: Duration,
    pub(super) mphf: Duration,
}

pub(super) fn write_index_with_report(
    path: &Path,
    graph: &BuildGraph,
    threads: usize,
) -> Result<IndexWriteReport> {
    let total_start = Instant::now();
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, file);
    write_u64(&mut writer, INDEX_VERSION)?;

    let mut graph_section = Vec::new();
    write_bifrost_graph_and_meta(&mut graph_section, graph)?;
    write_u64(&mut writer, graph_section.len() as u64)?;
    writer.write_all(&graph_section)?;
    super::scratch_pages::retire(graph_section);

    let mphf_start = Instant::now();
    let mut mphf = Vec::new();
    write_boo_mphf(&mut mphf, &graph.minimizer_keys, threads)?;
    write_u64(&mut writer, mphf.len() as u64)?;
    writer.write_all(&mphf)?;
    let mphf = mphf_start.elapsed();

    write_u64(&mut writer, 0)?;
    write_u64(&mut writer, 1)?;

    write_nodes(&mut writer, graph, threads)?;
    write_transcript_metadata(&mut writer, graph)?;
    writer.flush()?;
    Ok(IndexWriteReport {
        total: total_start.elapsed(),
        mphf,
    })
}

fn write_bifrost_graph_and_meta(out: &mut Vec<u8>, graph: &BuildGraph) -> Result<()> {
    write_u64(out, (BFG_GRAPHBIN_FORMAT_HEADER << 32) | BFG_FORMAT_VERSION)?;
    write_i32(out, graph.k as i32)?;
    write_i32(out, graph.g as i32)?;

    write_u64(out, graph.unitigs.len() as u64)?;
    let mut checksum = wyhash(&(graph.k as u64).to_le_bytes(), 0);
    checksum = wyhash(&(graph.g as u64).to_le_bytes(), checksum);
    for seq in &graph.unitigs {
        write_u64(out, seq.len() as u64)?;
        let range = append_compressed_sequence(out, seq);
        checksum = wyhash(&out[range], checksum);
    }

    write_u64(out, graph.km_unitigs.len() as u64)?;
    for seq in &graph.km_unitigs {
        let word = encode_kmer(seq)?;
        write_kmer_word(out, word)?;
        checksum = wyhash(&word.to_le_bytes(), checksum);
    }
    write_u64(out, 0)?;

    write_u64(out, (BFG_METABIN_FORMAT_HEADER << 32) | BFG_FORMAT_VERSION)?;
    write_u64(out, checksum)?;
    write_u64(out, graph.unitigs.len() as u64)?;
    write_u64(out, graph.km_unitigs.len() as u64)?;
    write_u64(out, 0)?;
    write_u64(out, graph.minimizer_keys.len() as u64)?;

    write_u64(out, graph.unitig_bitmap_blocks.len() as u64)?;
    for vals in &graph.unitig_bitmap_blocks {
        write_bitcontainer(out, vals)?;
    }

    write_u64(out, graph.km_bitmap_blocks.len() as u64)?;
    for vals in &graph.km_bitmap_blocks {
        write_bitcontainer(out, vals)?;
    }

    write_u64(out, 0)?;
    write_bitcontainer(out, &[])?;
    write_bitcontainer(out, &[])?;
    Ok(())
}

struct BooBitVector {
    size: u64,
    words: Vec<u64>,
    ranks: Vec<u64>,
}

impl BooBitVector {
    fn new(size: u64) -> Result<Self> {
        let nchar = size
            .checked_div(64)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| Error::InvalidFormat("MPHF bitvector too large".into()))?;
        let nchar = usize::try_from(nchar)
            .map_err(|_| Error::InvalidFormat("MPHF bitvector too large".into()))?;
        Ok(Self {
            size,
            words: vec![0; nchar],
            ranks: Vec::new(),
        })
    }

    fn set(&mut self, pos: u64) {
        self.words[(pos >> 6) as usize] |= 1u64 << (pos & 63);
    }

    fn build_ranks(&mut self, offset: u64) -> u64 {
        self.ranks.clear();
        let mut rank = offset;
        for (idx, &word) in self.words.iter().enumerate() {
            if (idx as u64 * 64).is_multiple_of(512) {
                self.ranks.push(rank);
            }
            rank += u64::from(word.count_ones());
        }
        rank
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        write_u64(writer, self.size)?;
        write_u64(writer, self.words.len() as u64)?;
        for &word in &self.words {
            write_u64(writer, word)?;
        }
        write_u64(writer, self.ranks.len() as u64)?;
        for &rank in &self.ranks {
            write_u64(writer, rank)?;
        }
        Ok(())
    }
}

trait MphfIndex: Copy + Send + Sync {
    fn from_index(index: usize) -> Self;
    fn index(self) -> usize;
}
impl MphfIndex for u32 {
    fn from_index(index: usize) -> Self {
        index as u32
    }
    fn index(self) -> usize {
        self as usize
    }
}
impl MphfIndex for usize {
    fn from_index(index: usize) -> Self {
        index
    }
    fn index(self) -> usize {
        self
    }
}

trait MphfBucket: Copy + Send + Sync {
    fn from_bucket(bucket: u64) -> Self;
    fn bucket(self) -> u64;
}
impl MphfBucket for u32 {
    fn from_bucket(bucket: u64) -> Self {
        bucket as u32
    }
    fn bucket(self) -> u64 {
        self as u64
    }
}
impl MphfBucket for u64 {
    fn from_bucket(bucket: u64) -> Self {
        bucket
    }
    fn bucket(self) -> u64 {
        self
    }
}

fn write_boo_mphf(out: &mut Vec<u8>, keys: &[u64], threads: usize) -> Result<()> {
    if boo_level_domains(keys.len())?
        .iter()
        .all(|&domain| domain <= u32::MAX as u64)
    {
        write_boo_mphf_ids::<u32, u32>(out, keys, threads)
    } else if keys.len() <= u32::MAX as usize {
        write_boo_mphf_ids::<u32, u64>(out, keys, threads)
    } else {
        write_boo_mphf_ids::<usize, u64>(out, keys, threads)
    }
}

#[derive(Clone, Copy, Default)]
struct MphfMarks {
    seen: u64,
    collisions: u64,
}

fn write_boo_mphf_ids<Id: MphfIndex, Bucket: MphfBucket>(
    out: &mut Vec<u8>,
    keys: &[u64],
    threads: usize,
) -> Result<()> {
    let pool = if threads > 1 && keys.len() >= 1024 {
        Some(IndexWorkers::new(threads)?)
    } else {
        None
    };
    let domains = boo_level_domains(keys.len())?;
    let mut survivors: Vec<Id> = (0..keys.len()).map(Id::from_index).collect();
    let mut bitvectors = Vec::with_capacity(BOO_LEVELS);
    let mut offset = 0u64;

    for (level, &domain) in domains.iter().enumerate() {
        let mut bitvector = BooBitVector::new(domain)?;
        if level + 1 < BOO_LEVELS && !survivors.is_empty() {
            if let Some(pool) = pool.as_ref().filter(|_| survivors.len() >= 1024) {
                let chunk = survivors
                    .len()
                    .div_ceil(threads.min(survivors.len().div_ceil(256)).max(1));
                let reduced = ReducedDomain::new(domain);
                let shards = pool.install(|| {
                    survivors
                        .par_chunks(chunk)
                        .map(|ids| {
                            // Both fields are touched for each hashed bit. Keep
                            // them on the same cache line without atomic updates.
                            let mut marks = vec![MphfMarks::default(); bitvector.words.len()];
                            let buckets = ids
                                .iter()
                                .map(|&key_idx| {
                                    let bucket = reduced.remainder(boo_level_hash(
                                        &keys[key_idx.index()].to_be_bytes(),
                                        level,
                                    ));
                                    let word = (bucket >> 6) as usize;
                                    let mask = 1 << (bucket & 63);
                                    let mark = &mut marks[word];
                                    mark.collisions |= mark.seen & mask;
                                    mark.seen |= mask;
                                    Bucket::from_bucket(bucket)
                                })
                                .collect::<Vec<_>>();
                            (marks, buckets)
                        })
                        .collect::<Vec<_>>()
                });
                pool.install(|| {
                    bitvector
                        .words
                        .par_iter_mut()
                        .enumerate()
                        .for_each(|(word, dest)| {
                            let mut seen = 0;
                            let mut collisions = 0;
                            for (marks, _) in &shards {
                                let mark = marks[word];
                                collisions |= mark.collisions | (seen & mark.seen);
                                seen |= mark.seen;
                            }
                            *dest = seen & !collisions;
                        })
                });
                let next = pool.install(|| {
                    survivors
                        .par_chunks(chunk)
                        .zip(shards.par_iter())
                        .map(|(ids, (_, buckets))| {
                            ids.iter()
                                .zip(buckets)
                                .filter_map(|(&key_idx, &bucket)| {
                                    let bucket = bucket.bucket();
                                    (bitvector.words[(bucket >> 6) as usize] & (1 << (bucket & 63))
                                        == 0)
                                        .then_some(key_idx)
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>()
                });
                survivors = next.into_iter().flatten().collect();
            } else {
                let domain_usize = usize::try_from(domain)
                    .map_err(|_| Error::InvalidFormat("MPHF domain too large".into()))?;
                let mut counts = vec![0u8; domain_usize];
                let mut buckets = Vec::with_capacity(survivors.len());
                let reduced_domain = ReducedDomain::new(domain);
                for &key_idx in &survivors {
                    let bucket = reduced_domain
                        .remainder(boo_level_hash(&keys[key_idx.index()].to_be_bytes(), level));
                    let count = &mut counts[bucket as usize];
                    *count = count.saturating_add(1).min(2);
                    buckets.push(Bucket::from_bucket(bucket));
                }

                let mut next_survivors = Vec::new();
                for (&key_idx, bucket) in survivors.iter().zip(buckets) {
                    let bucket = bucket.bucket();
                    if counts[bucket as usize] == 1 {
                        bitvector.set(bucket);
                    } else {
                        next_survivors.push(key_idx);
                    }
                }
                survivors = next_survivors;
            }
        }
        offset = bitvector.build_ranks(offset);
        bitvectors.push(bitvector);
    }

    write_f64(out, BOO_GAMMA)?;
    write_i32(out, BOO_LEVELS as i32)?;
    write_u64(out, offset)?;
    write_u64(out, keys.len() as u64)?;
    for bitvector in &bitvectors {
        bitvector.write(out)?;
    }

    write_u64(out, survivors.len() as u64)?;
    for (idx, &key_idx) in survivors.iter().enumerate() {
        out.write_all(&keys[key_idx.index()].to_be_bytes())?;
        write_u64(out, idx as u64)?;
    }
    Ok(())
}

fn boo_level_domains(n: usize) -> Result<Vec<u64>> {
    let nelem = n as f64;
    let hash_domain = (nelem * BOO_GAMMA).ceil();
    let proba_collision = if n == 0 {
        1.0
    } else {
        1.0 - (((BOO_GAMMA * nelem - 1.0) / (BOO_GAMMA * nelem)).powf(nelem - 1.0))
    };
    let mut out = Vec::with_capacity(BOO_LEVELS);
    for level in 0..BOO_LEVELS {
        let raw = (hash_domain * proba_collision.powi(level as i32)) as u64;
        let domain = raw
            .checked_add(63)
            .map(|v| (v / 64) * 64)
            .ok_or_else(|| Error::InvalidFormat("MPHF domain too large".into()))?;
        out.push(domain.max(64));
    }
    Ok(out)
}

fn boo_level_hash(key: &[u8; 8], level: usize) -> u64 {
    if level == 0 {
        return minimizer_hash(key, 0xAAAA_AAAA_5555_5555);
    }
    let mut s0 = minimizer_hash(key, 0xAAAA_AAAA_5555_5555);
    let mut s1 = minimizer_hash(key, 0x3333_3333_CCCC_CCCC);
    if level == 1 {
        return s1;
    }
    let mut out = 0;
    for _ in 2..=level {
        out = xorshift_next(&mut s0, &mut s1);
    }
    out
}

fn xorshift_next(s0: &mut u64, s1: &mut u64) -> u64 {
    let mut x = *s0;
    let y = *s1;
    *s0 = y;
    x ^= x << 23;
    *s1 = x ^ y ^ (x >> 17) ^ (y >> 26);
    s1.wrapping_add(y)
}

fn write_nodes<W: Write>(writer: &mut W, graph: &BuildGraph, threads: usize) -> Result<()> {
    write_u64(writer, graph.nodes.len() as u64)?;
    if threads > 1 && graph.nodes.len() >= 1024 {
        let pool = IndexWorkers::new(threads)?;
        // Keep only a bounded ordered batch of serialized nodes in memory.
        for first in (0..graph.nodes.len()).step_by(2048 * threads) {
            let end = (first + 2048 * threads).min(graph.nodes.len());
            let records = pool.install(|| {
                graph.nodes[first..end]
                    .par_chunks(128)
                    .enumerate()
                    .map(|(chunk, nodes)| {
                        let mut record = Vec::with_capacity(nodes.len() * (graph.k + 128));
                        for (relative, node) in nodes.iter().enumerate() {
                            let idx = first + chunk * 128 + relative;
                            let seq = if idx < graph.unitigs.len() {
                                &graph.unitigs[idx]
                            } else {
                                &graph.km_unitigs[idx - graph.unitigs.len()]
                            };
                            record.extend_from_slice(&seq[..graph.k]);
                            let length_offset = record.len();
                            record.extend_from_slice(&[0; 4]);
                            write_u32(&mut record, node.id)?;
                            write_block_array(&mut record, node)?;
                            let length = (record.len() - length_offset - 4) as u32;
                            record[length_offset..length_offset + 4]
                                .copy_from_slice(&length.to_le_bytes());
                        }
                        Ok::<_, Error>(record)
                    })
                    .collect::<Result<Vec<_>>>()
            })?;
            for record in records {
                writer.write_all(&record)?;
            }
        }
        return Ok(());
    }
    for (seq, node) in graph
        .unitigs
        .iter()
        .chain(&graph.km_unitigs)
        .zip(&graph.nodes)
    {
        writer.write_all(&seq[..graph.k])?;
        let mut payload = Vec::new();
        write_u32(&mut payload, node.id)?;
        write_block_array(&mut payload, node)?;
        write_u32(writer, payload.len() as u32)?;
        writer.write_all(&payload)?;
    }
    Ok(())
}

fn write_block_array(out: &mut Vec<u8>, node: &NodePayload) -> Result<()> {
    match node.blocks.as_slice() {
        [] => write_u8(out, 0),
        [block] => {
            write_u8(out, 1)?;
            write_u32(out, block.lb)?;
            write_u32(out, block.ub)?;
            write_block_sparse_vector(out, block)
        }
        blocks => {
            write_u8(out, 2)?;
            write_u64(out, blocks.len() as u64)?;
            for block in blocks {
                write_u32(out, block.lb)?;
                write_u32(out, block.ub)?;
                write_block_sparse_vector(out, block)?;
            }
            Ok(())
        }
    }
}

fn write_block_sparse_vector(out: &mut Vec<u8>, block: &NodeBlock) -> Result<()> {
    write_prefixed_bitmap::<false>(out, &block.ec)?;
    write_u64(out, block.positions.len() as u64)?;
    for vals in block.positions.iter() {
        write_prefixed_bitmap::<true>(out, vals)?;
    }
    Ok(())
}

// The locked native library chooses runs only on a strict size win. For a
// single run container its portable header omits offsets; the native wrapper
// chooses this payload only when smaller than the plain-u32 representation.
fn write_single_run_bitmap(out: &mut Vec<u8>, values: &[u32]) -> bool {
    if !cfg!(target_endian = "little")
        || values.len() < 3
        || values[0] >> 16 != values[values.len() - 1] >> 16
    {
        return false;
    }
    let mut runs = 1usize;
    for pair in values.windows(2) {
        if pair[0] >= pair[1] {
            return false;
        }
        runs += usize::from(pair[1] - pair[0] != 1);
    }
    let container = 2 + 4 * runs;
    let original = if values.len() <= 4096 {
        2 * values.len()
    } else {
        8192
    };
    let portable = 9 + container;
    if container >= original || portable >= 4 + 4 * values.len() {
        return false;
    }
    out.reserve(8 + 1 + portable);
    out.extend_from_slice(&((1 + portable) as u64).to_le_bytes());
    out.push(2);
    out.extend_from_slice(&12347u32.to_le_bytes());
    out.push(1);
    out.extend_from_slice(&((values[0] >> 16) as u16).to_le_bytes());
    out.extend_from_slice(&((values.len() - 1) as u16).to_le_bytes());
    out.extend_from_slice(&(runs as u16).to_le_bytes());
    let mut start = values[0];
    let mut previous = start;
    for &value in &values[1..] {
        if value != previous + 1 {
            out.extend_from_slice(&(start as u16).to_le_bytes());
            out.extend_from_slice(&((previous - start) as u16).to_le_bytes());
            start = value;
        }
        previous = value;
    }
    out.extend_from_slice(&(start as u16).to_le_bytes());
    out.extend_from_slice(&((previous - start) as u16).to_le_bytes());
    true
}

#[inline]
fn write_prefixed_bitmap<const RUN: bool>(out: &mut Vec<u8>, values: &[u32]) -> Result<()> {
    if RUN && write_single_run_bitmap(out, values) {
        return Ok(());
    }
    let run_optimize = RUN;
    // A run costs at least6bytes, so at most3values cannot strictly beat
    // an array container even when adjacent. Preserve native plain-u32 bytes.
    let ordered = |require_gaps: bool| {
        values
            .windows(2)
            .all(|pair| pair[0] < pair[1] && (!require_gaps || pair[1] - pair[0] > 1))
    };
    if cfg!(target_endian = "little")
        && values.len() <= 6
        && ordered(run_optimize && values.len() > 3)
    {
        write_u64(out, (5 + 4 * values.len()) as u64)?;
        write_u8(out, 1)?;
        write_u32(out, values.len() as u32)?;
        for &value in values {
            write_u32(out, value)?;
        }
    } else if cfg!(target_endian = "little")
        && (7..=4096).contains(&values.len())
        && values[0] >> 16 == values[values.len() - 1] >> 16
        && ordered(run_optimize)
    {
        // One sorted array container (<=4096 values). With run optimization,
        // nonadjacent values cannot benefit from a run container. CRoaring
        // selects this exact portable payload over its larger plain-u32 array:
        // native kind 2, no-run cookie, one key/cardinality, offset16, u16 data.
        write_u64(out, (17 + 2 * values.len()) as u64)?;
        write_u8(out, 2)?;
        write_u32(out, 12346)?;
        write_u32(out, 1)?;
        out.extend_from_slice(&((values[0] >> 16) as u16).to_le_bytes());
        out.extend_from_slice(&((values.len() - 1) as u16).to_le_bytes());
        write_u32(out, 16)?;
        for &value in values {
            out.extend_from_slice(&(value as u16).to_le_bytes());
        }
    } else {
        let bitmap = serialize_roaring(values, run_optimize, false)?;
        write_u64(out, bitmap.len() as u64)?;
        out.write_all(&bitmap)?;
    }
    Ok(())
}

fn write_transcript_metadata<W: Write>(writer: &mut W, graph: &BuildGraph) -> Result<()> {
    write_i32(writer, graph.transcript_names.len() as i32)?;
    for &len in &graph.transcript_lengths {
        write_i32(writer, len as i32)?;
    }
    for name in &graph.transcript_names {
        write_u64(writer, name.len() as u64)?;
        writer.write_all(name.as_bytes())?;
    }
    let onlist = serialize_roaring_range(graph.transcript_names.len() as u32)?;
    write_u64(writer, onlist.len() as u64)?;
    writer.write_all(&onlist)?;
    Ok(())
}

fn write_bitcontainer(out: &mut Vec<u8>, values: &[u32]) -> Result<()> {
    let bitmap = serialize_roaring(values, true, true)?;
    write_u64(
        out,
        ((bitmap.len() as u64) << SHIFT_MASK_BITS) | FLAG_PTR_BITMAP,
    )?;
    out.write_all(&bitmap)?;
    Ok(())
}

fn append_compressed_sequence(out: &mut Vec<u8>, seq: &[u8]) -> std::ops::Range<usize> {
    let start = out.len();
    out.reserve(seq.len().div_ceil(4));
    let code = |base| match base {
        b'A' => 0u8,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    };
    #[cfg(target_arch = "aarch64")]
    let seq = {
        let (chunks, remainder) = seq.as_chunks::<16>();
        for chunk in chunks {
            if let Some(packed) = pack_sixteen(chunk) {
                out.extend_from_slice(&packed);
            } else {
                out.extend(chunk.as_chunks::<4>().0.iter().map(|b| {
                    code(b[0]) | (code(b[1]) << 2) | (code(b[2]) << 4) | (code(b[3]) << 6)
                }));
            }
        }
        remainder
    };
    let (chunks, remainder) = seq.as_chunks::<4>();
    out.extend(chunks.iter().map(|bases| {
        code(bases[0]) | (code(bases[1]) << 2) | (code(bases[2]) << 4) | (code(bases[3]) << 6)
    }));
    if !remainder.is_empty() {
        let mut tail = 0u8;
        for (i, &base) in remainder.iter().enumerate() {
            tail |= code(base) << (2 * i);
        }
        out.push(tail);
    }
    start..out.len()
}

#[cfg(target_arch = "aarch64")]
fn pack_sixteen(bases: &[u8]) -> Option<[u8; 4]> {
    use std::arch::aarch64::*;
    debug_assert_eq!(bases.len(), 16);
    // The slice contains16bytes. These NEON operations are available on every
    // aarch64 target; reject non-ACGT bytes before applying the compact mapping.
    unsafe {
        let x = vld1q_u8(bases.as_ptr());
        let valid = vorrq_u8(
            vorrq_u8(vceqq_u8(x, vdupq_n_u8(b'A')), vceqq_u8(x, vdupq_n_u8(b'C'))),
            vorrq_u8(vceqq_u8(x, vdupq_n_u8(b'G')), vceqq_u8(x, vdupq_n_u8(b'T'))),
        );
        if vminvq_u8(valid) == 0 {
            return None;
        }
        let codes = vandq_u8(
            veorq_u8(vshrq_n_u8::<1>(x), vshrq_n_u8::<2>(x)),
            vdupq_n_u8(3),
        );
        let weights = [1, 4, 16, 64, 1, 4, 16, 64, 1, 4, 16, 64, 1, 4, 16, 64];
        let groups = vpaddlq_u16(vpaddlq_u8(vmulq_u8(codes, vld1q_u8(weights.as_ptr()))));
        Some([
            vgetq_lane_u32::<0>(groups) as u8,
            vgetq_lane_u32::<1>(groups) as u8,
            vgetq_lane_u32::<2>(groups) as u8,
            vgetq_lane_u32::<3>(groups) as u8,
        ])
    }
}

fn serialize_roaring(values: &[u32], run_optimize: bool, portable: bool) -> Result<Vec<u8>> {
    let ptr = unsafe { roaring_bitmap_create() };
    if ptr.is_null() {
        return Err(Error::InvalidFormat(
            "failed to create roaring bitmap".into(),
        ));
    }
    if !values.is_empty() {
        unsafe {
            roaring_bitmap_add_many(ptr, values.len(), values.as_ptr());
        }
    }
    if run_optimize {
        unsafe {
            roaring_bitmap_run_optimize(ptr);
        }
    }
    let size = if portable {
        unsafe { roaring_bitmap_portable_size_in_bytes(ptr) }
    } else {
        unsafe { roaring_bitmap_size_in_bytes(ptr) }
    };
    let mut out = vec![0u8; size];
    unsafe {
        if portable {
            roaring_bitmap_portable_serialize(ptr, out.as_mut_ptr().cast());
        } else {
            roaring_bitmap_serialize(ptr, out.as_mut_ptr().cast());
        }
        roaring_bitmap_free(ptr);
    }
    Ok(out)
}

fn serialize_roaring_range(end: u32) -> Result<Vec<u8>> {
    let ptr = unsafe { roaring_bitmap_create() };
    if ptr.is_null() {
        return Err(Error::InvalidFormat(
            "failed to create roaring bitmap".into(),
        ));
    }
    unsafe {
        roaring_bitmap_add_range(ptr, 0, u64::from(end));
        roaring_bitmap_run_optimize(ptr);
    }
    let size = unsafe { roaring_bitmap_portable_size_in_bytes(ptr) };
    let mut out = vec![0u8; size];
    unsafe {
        roaring_bitmap_portable_serialize(ptr, out.as_mut_ptr().cast());
        roaring_bitmap_free(ptr);
    }
    Ok(out)
}

fn write_kmer_word<W: Write>(writer: &mut W, word: u64) -> Result<()> {
    writer.write_all(&word.to_le_bytes())?;
    Ok(())
}

fn write_u8<W: Write>(writer: &mut W, v: u8) -> Result<()> {
    writer.write_all(&[v])?;
    Ok(())
}

fn write_u32<W: Write>(writer: &mut W, v: u32) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_i32<W: Write>(writer: &mut W, v: i32) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_u64<W: Write>(writer: &mut W, v: u64) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_f64<W: Write>(writer: &mut W, v: f64) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::io::Cursor;

    use super::*;
    use crate::index::bifrost::BooPhf;

    #[test]
    fn compressed_sequence_append_preserves_native_bits_and_tail() {
        for len in 0..128 {
            let sequence = b"ACGTNacgu-"
                .iter()
                .copied()
                .cycle()
                .take(len)
                .collect::<Vec<_>>();
            let mut expected = vec![0u8; len.div_ceil(4)];
            for (i, &base) in sequence.iter().enumerate() {
                let code = match base {
                    b'A' => 0u8,
                    b'C' => 1,
                    b'G' => 2,
                    _ => 3,
                };
                expected[i / 4] |= code << (2 * (i % 4));
            }
            let mut actual = vec![0x12, 0x34];
            let range = append_compressed_sequence(&mut actual, &sequence);
            assert_eq!(&actual[..2], &[0x12, 0x34]);
            assert_eq!(&actual[range], expected);
        }
    }

    #[test]
    fn small_bitmap_records_match_croaring() {
        let mut cases = vec![
            vec![],
            vec![0],
            vec![u32::MAX],
            vec![0, u32::MAX],
            vec![1, 1],
            vec![2, 1],
            vec![1, 2, 3],
            (0..6).collect(),
            (0..7).collect(),
            vec![0, 65536, 131072, 196608, 262144, 327680],
        ];
        for len in [3, 4, 5, 6, 7, 8, 16, 31, 255, 1024, 4096, 4097] {
            cases.push((0..len).map(|i| 65536 + 2 * i as u32).collect());
            cases.push((0..len).map(|i| i as u32).collect());
            cases.push((0..len).map(|i| (i as u32) << 16).collect());
        }
        cases.push(vec![65534, 65535, 65536]);
        cases.push(vec![65533, 65535, 65537, u32::MAX]);
        cases.push(vec![u32::MAX - 4, u32::MAX - 2, u32::MAX]);
        for size in [3, 4, 5, 6, 7, 8, 31, 4096, 4097, 65536] {
            cases.push((0..size).collect());
            cases.push((0..size).map(|i| i + i / 3).collect());
            cases.push((0..size.min(32000)).map(|i| 65536 + i + i / 7).collect());
        }
        for values in cases {
            for optimize in [false, true] {
                let mut actual = Vec::new();
                if optimize {
                    write_prefixed_bitmap::<true>(&mut actual, &values).unwrap();
                } else {
                    write_prefixed_bitmap::<false>(&mut actual, &values).unwrap();
                }
                let native = serialize_roaring(&values, optimize, false).unwrap();
                let mut expected = (native.len() as u64).to_le_bytes().to_vec();
                expected.extend(native);
                assert_eq!(actual, expected, "values={values:?}, optimize={optimize}");
            }
        }
    }

    #[test]
    fn boo_mphf_roundtrips_all_keys_without_all_final_hash() {
        let keys = (0..2000u64)
            .map(|i| u64::from_be_bytes(i.wrapping_mul(0x9e37_79b9_7f4a_7c15).to_le_bytes()))
            .collect::<Vec<_>>();

        let mut buf = Vec::new();
        write_boo_mphf(&mut buf, &keys, 1).expect("write mphf");
        let mut parallel = Vec::new();
        write_boo_mphf(&mut parallel, &keys, 8).expect("parallel mphf");
        assert_eq!(parallel, buf);
        for threads in [1, 8] {
            let mut wide = Vec::new();
            write_boo_mphf_ids::<usize, u64>(&mut wide, &keys, threads).expect("wide survivor IDs");
            assert_eq!(wide, buf);
            let mut wide_buckets = Vec::new();
            write_boo_mphf_ids::<u32, u64>(&mut wide_buckets, &keys, threads).unwrap();
            assert_eq!(wide_buckets, buf);
        }
        assert_eq!(
            <u32 as MphfBucket>::from_bucket(u32::MAX as u64).bucket(),
            u32::MAX as u64
        );
        assert_eq!(
            <u64 as MphfBucket>::from_bucket(u64::MAX).bucket(),
            u64::MAX
        );
        assert_eq!(
            <u32 as MphfIndex>::from_index(u32::MAX as usize).index(),
            u32::MAX as usize
        );
        assert_eq!(
            <usize as MphfIndex>::from_index(usize::MAX).index(),
            usize::MAX
        );

        let mphf = BooPhf::load(&mut Cursor::new(buf)).expect("load mphf");
        assert_eq!(mphf.size(), keys.len() as u64);
        assert!(mphf.final_hash_entries().len() < keys.len());

        let mut ranks = HashSet::new();
        for key in &keys {
            let rank = mphf.lookup(&key.to_be_bytes()).expect("lookup key");
            assert!(rank < keys.len() as u64);
            assert!(ranks.insert(rank), "duplicate rank {rank}");
        }
        assert_eq!(ranks.len(), keys.len());
    }
}
