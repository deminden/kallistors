//! Serialize the kallisto v13 layout, including upstream bitmap tags and widths.
//! In-memory packing choices must not leak into these on-disk representations.

mod mphf;
use mphf::write_boo_mphf;
#[cfg(test)]
use mphf::{MphfBucket, MphfIndex, write_boo_mphf_ids, write_boo_mphf_ids_first};

use super::binary::{write_f64, write_i32, write_u8, write_u32, write_u64};
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

use super::bifrost::wyhash;
use super::graph_build::{BuildGraph, NodeBlock, NodePayload, encode_kmer};
use crate::{Error, Result};

const INDEX_VERSION: u64 = 13;
const BFG_GRAPHBIN_FORMAT_HEADER: u64 = 0x7e21_5f3f;
const BFG_METABIN_FORMAT_HEADER: u64 = 0x267c_3d5d;
const BFG_FORMAT_VERSION: u64 = 1;
const FLAG_PTR_BITMAP: u64 = 0x3;
const SHIFT_MASK_BITS: u64 = 3;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct IndexWriteReport {
    pub(super) total: Duration,
    pub(super) mphf: Duration,
}

pub(super) fn write_index_with_report(
    path: &Path,
    graph: &mut BuildGraph,
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
    // These positions are now completely serialized. Node rows use sequence
    // prefixes, not the minimizer bitmaps, so they need not overlap the MPHF.
    for block in std::mem::take(&mut graph.unitig_bitmap_blocks) {
        super::scratch_pages::retire(block);
    }
    for block in std::mem::take(&mut graph.km_bitmap_blocks) {
        super::scratch_pages::retire(block);
    }

    let mphf_start = Instant::now();
    let mut mphf = Vec::new();
    write_boo_mphf(&mut mphf, &graph.minimizer_keys, threads)?;
    write_u64(&mut writer, mphf.len() as u64)?;
    writer.write_all(&mphf)?;
    let mphf = mphf_start.elapsed();

    // The supported builder emits no D-list entries but still writes both
    // required header fields: size and overhang.
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

#[cfg(test)]
mod tests;
