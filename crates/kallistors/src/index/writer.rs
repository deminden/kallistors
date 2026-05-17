use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

use croaring_sys::{
    roaring_bitmap_add_many, roaring_bitmap_add_range, roaring_bitmap_create, roaring_bitmap_free,
    roaring_bitmap_portable_serialize, roaring_bitmap_portable_size_in_bytes,
    roaring_bitmap_run_optimize, roaring_bitmap_serialize, roaring_bitmap_size_in_bytes,
};

use super::bifrost::{minimizer_hash, wyhash};
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

pub(super) fn write_index_with_report(path: &Path, graph: &BuildGraph) -> Result<IndexWriteReport> {
    let total_start = Instant::now();
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_u64(&mut writer, INDEX_VERSION)?;

    let mut graph_section = Vec::new();
    write_bifrost_graph_and_meta(&mut graph_section, graph)?;
    write_u64(&mut writer, graph_section.len() as u64)?;
    writer.write_all(&graph_section)?;

    let mphf_start = Instant::now();
    let mut mphf = Vec::new();
    write_boo_mphf(&mut mphf, &graph.minimizer_keys)?;
    write_u64(&mut writer, mphf.len() as u64)?;
    writer.write_all(&mphf)?;
    let mphf = mphf_start.elapsed();

    write_u64(&mut writer, 0)?;
    write_u64(&mut writer, 1)?;

    write_nodes(&mut writer, graph)?;
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
    for seq in &graph.unitigs {
        write_compressed_sequence(out, seq)?;
    }

    write_u64(out, graph.km_unitigs.len() as u64)?;
    for seq in &graph.km_unitigs {
        write_kmer_word(out, encode_kmer(seq)?)?;
    }
    write_u64(out, 0)?;

    write_u64(out, (BFG_METABIN_FORMAT_HEADER << 32) | BFG_FORMAT_VERSION)?;
    write_u64(out, graph_checksum(graph)?)?;
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

fn graph_checksum(graph: &BuildGraph) -> Result<u64> {
    let mut checksum = wyhash(&(graph.k as u64).to_le_bytes(), 0);
    checksum = wyhash(&(graph.g as u64).to_le_bytes(), checksum);
    for seq in &graph.unitigs {
        checksum = wyhash(&compressed_sequence_bytes(seq), checksum);
    }
    for seq in &graph.km_unitigs {
        let word = encode_kmer(seq)?;
        checksum = wyhash(&word.to_le_bytes(), checksum);
    }
    Ok(checksum)
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

fn write_boo_mphf(out: &mut Vec<u8>, keys: &[[u8; 8]]) -> Result<()> {
    let domains = boo_level_domains(keys.len())?;
    let mut survivors: Vec<usize> = (0..keys.len()).collect();
    let mut bitvectors = Vec::with_capacity(BOO_LEVELS);
    let mut offset = 0u64;

    for (level, &domain) in domains.iter().enumerate() {
        let mut bitvector = BooBitVector::new(domain)?;
        if level + 1 < BOO_LEVELS && !survivors.is_empty() {
            let domain_usize = usize::try_from(domain)
                .map_err(|_| Error::InvalidFormat("MPHF domain too large".into()))?;
            let mut counts = vec![0u8; domain_usize];
            let mut buckets = Vec::with_capacity(survivors.len());
            for &key_idx in &survivors {
                let bucket = boo_level_hash(&keys[key_idx], level) % domain;
                let count = &mut counts[bucket as usize];
                *count = count.saturating_add(1).min(2);
                buckets.push(bucket);
            }

            let mut next_survivors = Vec::new();
            for (&key_idx, bucket) in survivors.iter().zip(buckets) {
                if counts[bucket as usize] == 1 {
                    bitvector.set(bucket);
                } else {
                    next_survivors.push(key_idx);
                }
            }
            survivors = next_survivors;
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
        out.write_all(&keys[key_idx])?;
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

fn write_nodes<W: Write>(writer: &mut W, graph: &BuildGraph) -> Result<()> {
    write_u64(writer, graph.nodes.len() as u64)?;
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
    let ec_bitmap = serialize_roaring(&block.ec, false, false)?;
    write_u64(out, ec_bitmap.len() as u64)?;
    out.write_all(&ec_bitmap)?;

    write_u64(out, block.positions.len() as u64)?;
    for vals in &block.positions {
        let bitmap = serialize_roaring(vals, true, false)?;
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

fn write_compressed_sequence<W: Write>(writer: &mut W, seq: &[u8]) -> Result<()> {
    write_u64(writer, seq.len() as u64)?;
    writer.write_all(&compressed_sequence_bytes(seq))?;
    Ok(())
}

fn compressed_sequence_bytes(seq: &[u8]) -> Vec<u8> {
    let mut data = vec![0u8; seq.len().div_ceil(4)];
    for (idx, &base) in seq.iter().enumerate() {
        let code = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        };
        data[idx >> 2] |= code << ((idx & 0x3) << 1);
    }
    data
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
    fn boo_mphf_roundtrips_all_keys_without_all_final_hash() {
        let keys = (0..2000u64)
            .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15).to_le_bytes())
            .collect::<Vec<_>>();

        let mut buf = Vec::new();
        write_boo_mphf(&mut buf, &keys).expect("write mphf");
        let mphf = BooPhf::load(&mut Cursor::new(buf)).expect("load mphf");
        assert_eq!(mphf.size(), keys.len() as u64);
        assert!(mphf.final_hash_entries().len() < keys.len());

        let mut ranks = HashSet::new();
        for key in &keys {
            let rank = mphf.lookup(key).expect("lookup key");
            assert!(rank < keys.len() as u64);
            assert!(ranks.insert(rank), "duplicate rank {rank}");
        }
        assert_eq!(ranks.len(), keys.len());
    }
}
