//! Pseudoalignment algorithms.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use boomphf::Mphf;

use io::{read_u32_le, read_u64_le};

use crate::bias::BiasCounts;
use crate::index::{parse_graph_section, read_block_array_blocks};
use crate::{Error, Result};

mod bias;
mod bifrost;
mod core;
mod debug;
mod ec;
mod io;
mod minimizers;
mod onlist;
mod paired;
mod single;
mod threaded;
mod utils;
use bias::bias_hexamer_for_match;
pub use bifrost::{BifrostIndex, build_bifrost_index, build_bifrost_index_with_positions};
use core::{
    apply_strand_filter, block_index_for_position, ec_block_at, ec_for_read_bifrost,
    filter_ec_by_fragment,
};
pub use debug::{DebugFailReason, DebugReport, ReadTrace};
use debug::{ReadDebugState, debug_reason};
pub(super) use ec::{
    ec_has_offlist, encode_kmer_pair, filter_shades, filter_shades_in_place, intersect_sorted,
    lookup_ec, merge_sorted_unique, merge_sorted_unique_vec,
};
use onlist::{read_onlist, read_onlist_with_lengths};
pub use paired::{
    pseudoalign_paired_bifrost_debug_with_options, pseudoalign_paired_bifrost_debug_with_strand,
    pseudoalign_paired_bifrost_with_options, pseudoalign_paired_bifrost_with_strand,
    pseudoalign_paired_naive,
};
pub use single::{
    pseudoalign_single_end, pseudoalign_single_end_bifrost, pseudoalign_single_end_bifrost_debug,
    pseudoalign_single_end_bifrost_debug_with_options,
    pseudoalign_single_end_bifrost_debug_with_strand,
    pseudoalign_single_end_bifrost_debug_with_strand_and_filter,
    pseudoalign_single_end_bifrost_with_options, pseudoalign_single_end_bifrost_with_strand,
    pseudoalign_single_end_bifrost_with_strand_and_filter,
};
pub use threaded::{
    pseudoalign_paired_bifrost_with_options_threaded,
    pseudoalign_single_end_bifrost_with_options_threaded,
};
pub(super) use utils::merge_ec_counts;

const KMER_BYTES_CANDIDATES: [usize; 4] = [8, 16, 24, 32];
const BATCH_SIZE: usize = 10_000;
const MAX_FRAG_LEN: i64 = 1000;

/// A naive k-mer -> EC mapping built from a kallisto index.
pub struct KmerEcIndex {
    pub k: usize,
    pub mphf: Mphf<u64>,
    pub keys_by_index: Vec<u64>,
    pub ecs: Vec<Vec<u32>>,
    pub onlist: Option<Vec<bool>>,
}

/// Pseudoalignment EC counts.
pub struct EcCounts {
    pub ec_list: Vec<Vec<u32>>,
    pub counts: Vec<u32>,
    pub reads_processed: u64,
    pub reads_aligned: u64,
    pub bias: Option<BiasCounts>,
    pub fragment_length_stats: Option<FragmentLengthStats>,
    pub fragment_length_hist: Option<Vec<u32>>,
}

impl EcCounts {
    fn new(with_bias: bool) -> Self {
        Self {
            ec_list: Vec::new(),
            counts: Vec::new(),
            reads_processed: 0,
            reads_aligned: 0,
            bias: with_bias.then(BiasCounts::new),
            fragment_length_stats: None,
            fragment_length_hist: None,
        }
    }
}

/// Running fragment length statistics for paired-end reads.
#[derive(Debug, Clone, Copy)]
pub struct FragmentLengthStats {
    count: u64,
    mean: f64,
    m2: f64,
}

impl FragmentLengthStats {
    pub fn new() -> Self {
        Self {
            count: 0,
            mean: 0.0,
            m2: 0.0,
        }
    }

    pub fn add(&mut self, value: f64) {
        self.count += 1;
        let delta = value - self.mean;
        self.mean += delta / self.count as f64;
        let delta2 = value - self.mean;
        self.m2 += delta * delta2;
    }

    pub fn merge(&mut self, other: &Self) {
        if other.count == 0 {
            return;
        }
        if self.count == 0 {
            *self = *other;
            return;
        }
        let total = self.count + other.count;
        let delta = other.mean - self.mean;
        self.mean += delta * (other.count as f64 / total as f64);
        self.m2 +=
            other.m2 + delta * delta * (self.count as f64 * other.count as f64) / total as f64;
        self.count = total;
    }

    pub fn count(&self) -> u64 {
        self.count
    }

    pub fn mean(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.mean)
        }
    }

    pub fn sd(&self) -> Option<f64> {
        if self.count < 2 {
            None
        } else {
            Some((self.m2 / (self.count as f64 - 1.0)).sqrt())
        }
    }
}

impl Default for FragmentLengthStats {
    fn default() -> Self {
        Self::new()
    }
}

/// Count reads assigned to a unique transcript (EC size 1).
pub fn unique_pseudoaligned_reads(counts: &EcCounts) -> u64 {
    let mut total = 0u64;
    for (idx, ec) in counts.ec_list.iter().enumerate() {
        if ec.len() == 1 {
            total = total.saturating_add(counts.counts.get(idx).copied().unwrap_or(0) as u64);
        }
    }
    total
}

#[derive(Debug, Clone, Copy)]
pub struct FragmentFilter {
    pub fragment_length: u32,
    pub single_overhang: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct PseudoalignOptions {
    pub min_range: usize,
    pub do_union: bool,
    pub dfk_onlist: bool,
    pub strand_specific: Option<StrandSpecific>,
    pub no_jump: bool,
    pub bias: bool,
    pub max_bias: usize,
}

impl Default for PseudoalignOptions {
    fn default() -> Self {
        Self {
            min_range: 1,
            do_union: false,
            dfk_onlist: false,
            strand_specific: None,
            no_jump: false,
            bias: false,
            max_bias: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandSpecific {
    FR,
    RF,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Unstranded,
    Forward,
    Reverse,
}

#[derive(Debug, Clone, Copy)]
struct MatchInfo {
    unitig_id: usize,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
}

#[derive(Debug, Clone)]
struct Hit {
    unitig_id: usize,
    read_pos: usize,
    block_idx: usize,
    used_revcomp: bool,
}

#[derive(Debug)]
struct ReadEc {
    ec: Vec<u32>,
    best_match: Option<MatchInfo>,
    first_hit: Option<Hit>,
    hits: Vec<Hit>,
    had_offlist: bool,
    shade_union: Vec<u32>,
}

/// Build a naive k-mer -> EC map by scanning unitigs and EC blocks in the index.
pub fn build_kmer_ec_index(path: &Path) -> Result<KmerEcIndex> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;

    let mut k: Option<u32> = None;
    let mut kmer_bytes: Option<usize> = None;

    let sequences = if dbg_size > 0 {
        let dbg_start = reader.stream_position()?;
        let meta = parse_graph_section(&mut reader, dbg_start, dbg_size, &KMER_BYTES_CANDIDATES)?;
        k = Some(meta.k);
        kmer_bytes = Some(meta.kmer_bytes);

        reader.seek(SeekFrom::Start(dbg_start))?;
        let seqs = bifrost::read_unitig_sequences(
            &mut reader,
            dbg_start,
            dbg_size,
            meta.k,
            meta.kmer_bytes,
        )?;
        reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;

        let mphf_size = read_u64_le(&mut reader)?;
        reader.seek(SeekFrom::Current(mphf_size as i64))?;
        seqs
    } else {
        Vec::new()
    };

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let kmer_bytes =
        kmer_bytes.ok_or_else(|| Error::InvalidFormat("missing k-mer width".into()))?;
    let dlist_skip = dlist_size
        .checked_mul(kmer_bytes as u64)
        .ok_or_else(|| Error::InvalidFormat("d-list size overflow".into()))?;
    reader.seek(SeekFrom::Current(dlist_skip as i64))?;

    let node_count = read_u64_le(&mut reader)?;
    let k = k.ok_or_else(|| Error::InvalidFormat("missing k-mer length".into()))? as usize;
    if k > 32 {
        return Err(Error::UnsupportedFeature(
            "k-mer length > 32 not supported in MPHF path".into(),
        ));
    }

    if sequences.len() != node_count as usize {
        return Err(Error::InvalidFormat(
            "unitig and node counts mismatch".into(),
        ));
    }

    let mut map: HashMap<u64, Vec<u32>> = HashMap::new();

    for (idx, seq) in sequences.iter().enumerate() {
        let _ = idx;
        reader.seek(SeekFrom::Current(k as i64))?;
        let node_size = read_u32_le(&mut reader)? as usize;
        let mut buf = vec![0u8; node_size];
        reader.read_exact(&mut buf)?;

        let mut cur = std::io::Cursor::new(buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let blocks = read_block_array_blocks(&mut cur)?;
        for block in blocks {
            let start = block.lb as usize;
            let end = block.ub as usize;
            for pos in start..end {
                if pos + k > seq.len() {
                    break;
                }
                if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + k]) {
                    map.entry(fwd)
                        .and_modify(|existing| merge_sorted_unique(existing, &block.ec))
                        .or_insert_with(|| block.ec.clone());
                    if rev != fwd {
                        map.entry(rev)
                            .and_modify(|existing| merge_sorted_unique(existing, &block.ec))
                            .or_insert_with(|| block.ec.clone());
                    }
                }
            }
        }
    }

    let mut keys: Vec<u64> = map.keys().cloned().collect();
    keys.sort_unstable();
    let mphf = Mphf::new(1.7, &keys);
    let mut keys_by_index = vec![0u64; keys.len()];
    let mut ecs = vec![Vec::new(); keys.len()];
    for (key, ec) in map {
        let idx = mphf.hash(&key) as usize;
        keys_by_index[idx] = key;
        ecs[idx] = ec;
    }

    let onlist = read_onlist(&mut reader)?;

    Ok(KmerEcIndex {
        k,
        mphf,
        keys_by_index,
        ecs,
        onlist,
    })
}
