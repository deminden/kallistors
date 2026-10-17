//! Pseudoalignment algorithms.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::OnceLock;

use boomphf::Mphf;

use io::{read_u32_le, read_u64_le};

use crate::bias::BiasCounts;
use crate::index::{parse_graph_section, read_block_array_blocks};
use crate::{Error, Result};

mod bias;
mod bifrost;
mod core;
mod debug;
mod diagnostics;
mod ec;
mod io;
mod minimizers;
mod onlist;
mod paired;
mod single;
mod threaded;
mod utils;
use bias::bias_hexamer_for_match;
pub use bifrost::{
    BifrostIndex, build_bifrost_index, build_bifrost_index_quant_threaded,
    build_bifrost_index_threaded, build_bifrost_index_with_kmer_pos,
    build_bifrost_index_with_kmer_pos_threaded, build_bifrost_index_with_positions,
    build_bifrost_index_with_positions_threaded,
};
pub use core::local_kmer_hits;
use core::{
    apply_strand_filter, ec_for_read_bifrost, filter_ec_by_fragment, reset_thread_local_caches,
};
pub use debug::{
    DebugFailReason, DebugReport, DroppedHitTrace, JumpDecisionTrace, MinimizerCandidateTrace,
    ReadTrace,
};
use debug::{ReadDebugState, debug_reason};
use diagnostics::placement_for_match;
pub(crate) use diagnostics::trace_result_from_debug_state;
pub use diagnostics::{
    debug_minimizer_lookup, trace_read_bifrost, trace_read_bifrost_with_hits, write_index_ec_dump,
    write_index_ec_dump_filtered,
};
pub(super) use ec::{
    ec_has_offlist, encode_kmer_pair, filter_shades, filter_shades_in_place, intersect_sorted,
    lookup_ec, merge_sorted_unique, merge_sorted_unique_vec,
};
use onlist::{read_onlist, read_onlist_with_lengths};
pub use paired::{
    pseudoalign_paired_bifrost_debug_with_options, pseudoalign_paired_bifrost_debug_with_strand,
    pseudoalign_paired_bifrost_with_options, pseudoalign_paired_bifrost_with_strand,
    pseudoalign_paired_naive, trace_read_pair_bifrost,
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
pub(super) use utils::{add_ec_count, merge_ec_counts};

const KMER_BYTES_CANDIDATES: [usize; 4] = [8, 16, 24, 32];
const BATCH_SIZE: usize = 12_000;
pub const MAX_FRAG_LEN: i64 = 1000;
static RESET_ALL_CACHES_PER_READ: OnceLock<bool> = OnceLock::new();

#[inline]
fn reset_all_caches_per_read() -> bool {
    *RESET_ALL_CACHES_PER_READ
        .get_or_init(|| std::env::var_os("KALLISTORS_RESET_ALL_CACHES_PER_READ").is_some())
}

pub fn build_bifrost_index_with_kmer(path: &Path) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index(path)?;
    index.kmer_index = Some(build_kmer_ec_index(path)?);
    Ok(index)
}

pub fn build_bifrost_index_with_kmer_threaded(path: &Path, threads: usize) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index_threaded(path, threads)?;
    index.kmer_index = Some(build_kmer_ec_index(path)?);
    Ok(index)
}

pub fn build_bifrost_index_with_positions_and_kmer(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index_with_positions(path, load_positional_info)?;
    index.kmer_index = Some(build_kmer_ec_index(path)?);
    Ok(index)
}

pub fn build_bifrost_index_with_positions_and_kmer_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    let mut index =
        build_bifrost_index_with_positions_threaded(path, load_positional_info, threads)?;
    index.kmer_index = Some(build_kmer_ec_index(path)?);
    Ok(index)
}

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

#[derive(Debug, Clone)]
pub struct ReadTraceResult {
    pub ec_before_filters: Option<Vec<u32>>,
    pub ec_after_fragment_filter: Option<Vec<u32>>,
    pub ec_after_strand_filter: Option<Vec<u32>>,
    pub reason: Option<DebugFailReason>,
    pub kmer_pos: Option<usize>,
    pub min_pos: Option<usize>,
    pub kmer_seq: Option<String>,
    pub min_seq: Option<String>,
    pub positions: Option<usize>,
    pub sample_positions: Option<String>,
    pub used_revcomp: bool,
    pub positions_visited: Option<Vec<usize>>,
    pub dropped_hits: Option<Vec<debug::DroppedHitTrace>>,
    pub minimizer_candidates: Option<Vec<debug::MinimizerCandidateTrace>>,
    pub jump_decisions: Option<Vec<debug::JumpDecisionTrace>>,
    pub placement: Option<PseudoalignPlacement>,
}

#[derive(Debug, Clone)]
pub struct PairedTraceResult {
    pub left: ReadTraceResult,
    pub right: ReadTraceResult,
    pub merged_ec: Option<Vec<u32>>,
    pub merged_reason: Option<DebugFailReason>,
    pub hard_reject_pair: bool,
    pub had_offlist: bool,
    pub estimated_fragment_length: Option<i64>,
    pub placement: Option<PseudoalignPlacement>,
}

#[derive(Debug, Clone, Copy)]
pub struct PseudoalignPlacement {
    pub transcript_id: u32,
    pub start: usize,
    pub reverse: bool,
}

#[derive(Debug, Clone)]
pub struct HitTrace {
    pub unitig_id: usize,
    pub block_idx: usize,
    pub read_pos: usize,
    pub used_revcomp: bool,
    pub ec: Vec<u32>,
    pub offlist: bool,
}

#[derive(Debug, Clone)]
pub struct MinimizerDebugRow {
    pub read_pos: usize,
    pub kmer: String,
    pub kmer_canon: String,
    pub used_revcomp: bool,
    pub min_pos: Option<usize>,
    pub min_seq: Option<String>,
    pub mphf_key: Option<String>,
    pub mphf_level: Option<u32>,
    pub mphf_hash_raw: Option<u64>,
    pub mphf_hash_domain: Option<u64>,
    pub mphf_bucket: Option<u64>,
    pub mphf_rank: Option<u64>,
    pub mphf_final_hash_hit: Option<bool>,
    pub mphf_hit: bool,
    pub positions: Option<usize>,
    pub sample_positions: Option<String>,
    pub matched: bool,
    pub match_unitig_id: Option<usize>,
    pub match_unitig_pos: Option<usize>,
    pub match_block_idx: Option<usize>,
    pub is_km: Option<bool>,
    pub special_code: Option<u64>,
    pub special_in_dlist: Option<bool>,
    pub special_in_hmap: Option<bool>,
    pub special_uid: Option<usize>,
    pub note: String,
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

/// Pseudoalign one sequence and return its equivalence class.
pub fn ec_for_sequence_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<Vec<u32>> {
    let mut read_ec = ec_for_read_bifrost(index, seq, strand, None, options)?;
    if let Some(mode) = options.strand_specific {
        let comprehensive = options.do_union || options.no_jump || index.use_shade;
        if let Some(filtered) =
            apply_strand_filter(index, &read_ec, mode, true, comprehensive, &mut None, b"")
        {
            read_ec.ec = filtered;
        }
    }
    (!read_ec.ec.is_empty()).then_some(read_ec.ec)
}

pub fn unmapped_kmer_ratio_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    options: PseudoalignOptions,
) -> f64 {
    let total = seq.len().saturating_sub(index.k).saturating_add(1);
    if total == 0 {
        return 1.0;
    }
    let hits = local_kmer_hits(
        index,
        seq,
        options.kallisto_enum,
        options.kallisto_strict,
        options.skip_overcrowded_minimizer,
        options.kallisto_bifrost_find,
        usize::MAX,
    )
    .into_iter()
    .filter(|(_, _, hit, _)| *hit)
    .count();
    (total.saturating_sub(hits)) as f64 / total as f64
}

#[derive(Debug, Clone, Copy, Default)]
pub struct InvestigationOptions {
    pub trace_fast_path: bool,
    pub bounded_incremental_scan: bool,
    pub faster_probe_state: bool,
    pub altered_hit_bookkeeping: bool,
    pub candidate_reuse: bool,
}

impl InvestigationOptions {
    pub fn any_enabled(self) -> bool {
        self.trace_fast_path
            || self.bounded_incremental_scan
            || self.faster_probe_state
            || self.altered_hit_bookkeeping
            || self.candidate_reuse
    }
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
    pub kallisto_enum: bool,
    pub kallisto_strict: bool,
    pub kallisto_local_fallback: bool,
    pub kallisto_fallback: bool,
    pub discard_special_only: bool,
    pub skip_overcrowded_minimizer: bool,
    pub kallisto_direct_kmer: bool,
    pub kallisto_bifrost_find: bool,
    pub kallisto_sparse_hits: bool,
    pub bias: bool,
    pub max_bias: usize,
    /// Paired fragment observations used for estimation; zero collects all.
    pub fragment_length_samples: usize,
    pub investigation: InvestigationOptions,
}

impl Default for PseudoalignOptions {
    fn default() -> Self {
        Self {
            min_range: 1,
            do_union: false,
            dfk_onlist: false,
            strand_specific: None,
            no_jump: false,
            kallisto_enum: false,
            kallisto_strict: false,
            kallisto_local_fallback: false,
            kallisto_fallback: false,
            discard_special_only: false,
            skip_overcrowded_minimizer: false,
            kallisto_direct_kmer: false,
            kallisto_bifrost_find: false,
            kallisto_sparse_hits: false,
            bias: false,
            max_bias: 1_000_000,
            fragment_length_samples: 10_000,
            investigation: InvestigationOptions::default(),
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
    hard_reject_pair: bool,
}

/// Build a naive k-mer -> EC map by scanning unitigs and EC blocks in the index.
pub fn build_kmer_ec_index(path: &Path) -> Result<KmerEcIndex> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::new(file);

    let _index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;

    let mut k: Option<u32> = None;
    let mut kmer_bytes: Option<usize> = None;

    let (sequences, abundant_base) = if dbg_size > 0 {
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
        (Vec::new(), 0)
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

    let mut node_kmer = vec![0; k];
    for (idx, graph_seq) in sequences.iter().enumerate() {
        let seq = if idx >= abundant_base {
            reader.read_exact(&mut node_kmer)?;
            node_kmer.as_slice()
        } else {
            reader.seek(SeekFrom::Current(k as i64))?;
            graph_seq.as_slice()
        };
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
