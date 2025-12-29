//! Pseudoalignment algorithms.

use std::collections::{HashMap, HashSet};
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
pub use bifrost::{
    BifrostIndex, build_bifrost_index, build_bifrost_index_with_kmer_pos,
    build_bifrost_index_with_positions,
};
pub use core::local_kmer_hits;
use core::{
    apply_strand_filter, block_index_for_position, ec_block_at, ec_for_read_bifrost,
    filter_ec_by_fragment, special_unitig_for_kmer, special_unitig_for_kmer_raw,
};
use debug::format_positions_sample;
pub use debug::{DebugFailReason, DebugReport, ReadTrace};
use debug::{ReadDebugState, debug_reason};
pub(super) use ec::{
    ec_has_offlist, encode_kmer_pair, filter_shades, filter_shades_in_place, intersect_sorted,
    lookup_ec, merge_sorted_unique, merge_sorted_unique_vec,
};
use minimizers::{
    fill_revcomp, minhash_candidates_for_kmer, minimizer_for_kmer_strict, minimizers_for_kmer,
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
use std::io::Write;
pub use threaded::{
    pseudoalign_paired_bifrost_with_options_threaded,
    pseudoalign_single_end_bifrost_with_options_threaded,
};
pub(super) use utils::merge_ec_counts;

const KMER_BYTES_CANDIDATES: [usize; 4] = [8, 16, 24, 32];
const BATCH_SIZE: usize = 10_000;
const MAX_FRAG_LEN: i64 = 1000;

pub fn build_bifrost_index_with_kmer(path: &Path) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index(path)?;
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

pub fn write_index_ec_dump(index_path: &Path, out: &Path, limit: usize) -> Result<()> {
    write_index_ec_dump_filtered(index_path, out, limit, None)
}

pub fn write_index_ec_dump_filtered(
    index_path: &Path,
    out: &Path,
    limit: usize,
    unitig_filter: Option<&HashSet<usize>>,
) -> Result<()> {
    let index = build_bifrost_index(index_path)?;
    let mut writer = std::io::BufWriter::new(std::fs::File::create(out)?);
    writeln!(
        writer,
        "unitig_id\tunitig_pos\tkmer\tblock_id\tblock_lb\tblock_ub\tec"
    )?;
    let mut rows = 0usize;
    let mut remaining = unitig_filter.map(|v| v.len());
    for (unitig_id, blocks) in index.ec_blocks.iter().enumerate() {
        if let Some(filter) = unitig_filter
            && !filter.contains(&unitig_id)
        {
            continue;
        }
        let seq = if unitig_id < index.unitigs.len() {
            &index.unitigs[unitig_id]
        } else {
            let km_idx = unitig_id - index.unitigs.len();
            if km_idx < index.km_unitigs.len() {
                &index.km_unitigs[km_idx]
            } else {
                continue;
            }
        };
        for (block_id, block) in blocks.iter().enumerate() {
            let pos = block.lb as usize;
            let kmer = if pos + index.k <= seq.len() {
                String::from_utf8_lossy(&seq[pos..pos + index.k]).into_owned()
            } else if pos < seq.len() {
                String::from_utf8_lossy(&seq[pos..]).into_owned()
            } else {
                "-".to_string()
            };
            let ec_text = if block.ec.is_empty() {
                "-".to_string()
            } else {
                block
                    .ec
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                unitig_id, pos, kmer, block_id, block.lb, block.ub, ec_text
            )?;
            rows += 1;
            if limit > 0 && rows >= limit {
                return Ok(());
            }
        }
        if let Some(ref mut left) = remaining
            && *left > 0
        {
            *left -= 1;
            if *left == 0 {
                return Ok(());
            }
        }
    }
    Ok(())
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

/// Debug minimizer lookups and unitig matches for specific read positions.
pub fn debug_minimizer_lookup(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    positions: &[usize],
    kallisto_bifrost_find: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
) -> Vec<MinimizerDebugRow> {
    let allow_forward_base = strand != Strand::Reverse;
    let allow_rev_base = strand != Strand::Forward;
    let diff = index.k.saturating_sub(index.g);
    let mut rev_buf: Vec<u8> = Vec::new();
    let mut out = Vec::new();
    let format_mphf_key = |key: &[u8; 8]| -> String {
        let mut out = String::with_capacity(16);
        for byte in key {
            use std::fmt::Write;
            let _ = write!(out, "{byte:02x}");
        }
        out
    };

    for &pos in positions {
        if pos + index.k > seq.len() {
            out.push(MinimizerDebugRow {
                read_pos: pos,
                kmer: "-".to_string(),
                kmer_canon: "-".to_string(),
                used_revcomp: false,
                min_pos: None,
                min_seq: None,
                mphf_key: None,
                mphf_level: None,
                mphf_hash_raw: None,
                mphf_hash_domain: None,
                mphf_bucket: None,
                mphf_rank: None,
                mphf_final_hash_hit: None,
                mphf_hit: false,
                positions: None,
                sample_positions: None,
                matched: false,
                match_unitig_id: None,
                match_unitig_pos: None,
                match_block_idx: None,
                is_km: None,
                special_code: None,
                special_in_dlist: None,
                special_in_hmap: None,
                special_uid: None,
                note: "out_of_range".to_string(),
            });
            continue;
        }

        let kmer = &seq[pos..pos + index.k];
        let kmer_text = String::from_utf8_lossy(kmer).into_owned();
        if !fill_revcomp(kmer, &mut rev_buf) {
            out.push(MinimizerDebugRow {
                read_pos: pos,
                kmer: kmer_text.clone(),
                kmer_canon: kmer_text,
                used_revcomp: false,
                min_pos: None,
                min_seq: None,
                mphf_key: None,
                mphf_level: None,
                mphf_hash_raw: None,
                mphf_hash_domain: None,
                mphf_bucket: None,
                mphf_rank: None,
                mphf_final_hash_hit: None,
                mphf_hit: false,
                positions: None,
                sample_positions: None,
                matched: false,
                match_unitig_id: None,
                match_unitig_pos: None,
                match_block_idx: None,
                is_km: None,
                special_code: None,
                special_in_dlist: None,
                special_in_hmap: None,
                special_uid: None,
                note: "invalid_kmer".to_string(),
            });
            continue;
        }

        let use_revcomp = rev_buf.as_slice() < kmer;

        let kmer_canon = if use_revcomp {
            rev_buf.as_slice()
        } else {
            kmer
        };
        let kmer_canon_text = String::from_utf8_lossy(kmer_canon).into_owned();
        let min_input = if kallisto_bifrost_find || kallisto_strict {
            kmer_canon
        } else {
            kmer
        };
        let mut special_code: Option<u64> = None;
        let mut special_in_dlist: Option<bool> = None;
        let mut special_in_hmap: Option<bool> = None;
        let mut special_uid: Option<usize> = None;
        if let Some((fwd, rev)) = encode_kmer_pair(kmer) {
            let code = fwd.min(rev);
            special_code = Some(code);
            if let Some(dlist) = index.dlist.as_ref() {
                special_in_dlist = Some(dlist.contains(&code));
            }
            if let Some(map) = index.h_kmer_map.as_ref() {
                if let Some(&uid) = map.get(&code) {
                    special_in_hmap = Some(true);
                    special_uid = Some(uid);
                } else {
                    special_in_hmap = Some(false);
                }
            }
        }
        let (mut min_candidates, mut fallback_min) = if kallisto_bifrost_find {
            match minhash_candidates_for_kmer(min_input, index.g) {
                Some(v) => v,
                None => {
                    out.push(MinimizerDebugRow {
                        read_pos: pos,
                        kmer: kmer_text.clone(),
                        kmer_canon: kmer_canon_text,
                        used_revcomp: use_revcomp,
                        min_pos: None,
                        min_seq: None,
                        mphf_key: None,
                        mphf_level: None,
                        mphf_hash_raw: None,
                        mphf_hash_domain: None,
                        mphf_bucket: None,
                        mphf_rank: None,
                        mphf_final_hash_hit: None,
                        mphf_hit: false,
                        positions: None,
                        sample_positions: None,
                        matched: false,
                        match_unitig_id: None,
                        match_unitig_pos: None,
                        match_block_idx: None,
                        is_km: None,
                        special_code,
                        special_in_dlist,
                        special_in_hmap,
                        special_uid,
                        note: "no_minimizer".to_string(),
                    });
                    continue;
                }
            }
        } else {
            let candidates = if kallisto_strict {
                match minimizer_for_kmer_strict(min_input, index.g) {
                    Some(v) => vec![v],
                    None => {
                        out.push(MinimizerDebugRow {
                            read_pos: pos,
                            kmer: kmer_text.clone(),
                            kmer_canon: kmer_canon_text,
                            used_revcomp: use_revcomp,
                            min_pos: None,
                            min_seq: None,
                            mphf_key: None,
                            mphf_level: None,
                            mphf_hash_raw: None,
                            mphf_hash_domain: None,
                            mphf_bucket: None,
                            mphf_rank: None,
                            mphf_final_hash_hit: None,
                            mphf_hit: false,
                            positions: None,
                            sample_positions: None,
                            matched: false,
                            match_unitig_id: None,
                            match_unitig_pos: None,
                            match_block_idx: None,
                            is_km: None,
                            special_code,
                            special_in_dlist,
                            special_in_hmap,
                            special_uid,
                            note: "no_minimizer".to_string(),
                        });
                        continue;
                    }
                }
            } else {
                match minimizers_for_kmer(min_input, index.g) {
                    Some(v) => v,
                    None => {
                        out.push(MinimizerDebugRow {
                            read_pos: pos,
                            kmer: kmer_text.clone(),
                            kmer_canon: kmer_canon_text,
                            used_revcomp: use_revcomp,
                            min_pos: None,
                            min_seq: None,
                            mphf_key: None,
                            mphf_level: None,
                            mphf_hash_raw: None,
                            mphf_hash_domain: None,
                            mphf_bucket: None,
                            mphf_rank: None,
                            mphf_final_hash_hit: None,
                            mphf_hit: false,
                            positions: None,
                            sample_positions: None,
                            matched: false,
                            match_unitig_id: None,
                            match_unitig_pos: None,
                            match_block_idx: None,
                            is_km: None,
                            special_code,
                            special_in_dlist,
                            special_in_hmap,
                            special_uid,
                            note: "no_minimizer".to_string(),
                        });
                        continue;
                    }
                }
            };
            (candidates, None)
        };

        let mut fallback_active = false;
        let mut special_seen = false;
        'min_loop: loop {
            let mut skip_minimizer = false;
            for (min_bytes, min_pos) in min_candidates.iter().copied() {
                let min_pos_fwd = if use_revcomp {
                    diff.saturating_sub(min_pos)
                } else {
                    min_pos
                };
                let min_pos_rev = if use_revcomp {
                    min_pos
                } else {
                    diff.saturating_sub(min_pos)
                };
                let min_seq = kmer_canon
                    .get(min_pos..min_pos + index.g)
                    .map(|v| String::from_utf8_lossy(v).into_owned());
                let (min_idx, mphf_debug) = index.mphf.debug_lookup(&min_bytes);
                let mphf_key = Some(format_mphf_key(&min_bytes));
                let Some(min_idx) = min_idx else {
                    let mut note = "mphf_miss".to_string();
                    if fallback_active {
                        note.push_str("|minhash_fallback");
                    }
                    out.push(MinimizerDebugRow {
                        read_pos: pos,
                        kmer: kmer_text.clone(),
                        kmer_canon: kmer_canon_text.clone(),
                        used_revcomp: use_revcomp,
                        min_pos: Some(min_pos),
                        min_seq,
                        mphf_key,
                        mphf_level: Some(mphf_debug.level),
                        mphf_hash_raw: Some(mphf_debug.hash_raw),
                        mphf_hash_domain: Some(mphf_debug.hash_domain),
                        mphf_bucket: mphf_debug.bucket,
                        mphf_rank: mphf_debug.rank,
                        mphf_final_hash_hit: Some(mphf_debug.final_hash_hit),
                        mphf_hit: false,
                        positions: None,
                        sample_positions: None,
                        matched: false,
                        match_unitig_id: None,
                        match_unitig_pos: None,
                        match_block_idx: None,
                        is_km: None,
                        special_code,
                        special_in_dlist,
                        special_in_hmap,
                        special_uid,
                        note,
                    });
                    continue;
                };

                let positions = &index.minz_positions[min_idx as usize];
                if positions.is_empty() {
                    let mut note = "no_positions".to_string();
                    if fallback_active {
                        note.push_str("|minhash_fallback");
                    }
                    out.push(MinimizerDebugRow {
                        read_pos: pos,
                        kmer: kmer_text.clone(),
                        kmer_canon: kmer_canon_text.clone(),
                        used_revcomp: use_revcomp,
                        min_pos: Some(min_pos),
                        min_seq,
                        mphf_key,
                        mphf_level: Some(mphf_debug.level),
                        mphf_hash_raw: Some(mphf_debug.hash_raw),
                        mphf_hash_domain: Some(mphf_debug.hash_domain),
                        mphf_bucket: mphf_debug.bucket,
                        mphf_rank: mphf_debug.rank,
                        mphf_final_hash_hit: Some(mphf_debug.final_hash_hit),
                        mphf_hit: true,
                        positions: Some(0),
                        sample_positions: Some("-".to_string()),
                        matched: false,
                        match_unitig_id: None,
                        match_unitig_pos: None,
                        match_block_idx: None,
                        is_km: None,
                        special_code,
                        special_in_dlist,
                        special_in_hmap,
                        special_uid,
                        note,
                    });
                    continue;
                }

                let mut matched = None;
                let mut special_note: Option<&str> = None;
                for &pos_id in positions {
                    let unitig_id_raw = (pos_id >> 32) as u32;
                    if unitig_id_raw == u32::MAX {
                        special_seen = true;
                        if kallisto_strict {
                            special_note = Some("special_strict_skip");
                            continue;
                        }
                        let overcrowded = (pos_id & 0x8000_0000) != 0;
                        if overcrowded {
                            special_note = Some("special_overcrowded");
                            if skip_overcrowded_minimizer {
                                if kallisto_bifrost_find {
                                    skip_minimizer = true;
                                    break;
                                }
                                continue;
                            }
                        }
                        let special_uid = if kallisto_bifrost_find {
                            special_unitig_for_kmer_raw(index, kmer)
                        } else {
                            special_unitig_for_kmer(index, kmer)
                        };
                        if let Some(uid) = special_uid {
                            let used_revcomp = kmer_canon != kmer;
                            if (used_revcomp && !allow_rev_base)
                                || (!used_revcomp && !allow_forward_base)
                            {
                                special_note = Some("special_abundant_strand_skip");
                                continue;
                            }
                            let block_idx =
                                block_index_for_position(index.ec_blocks[uid].as_slice(), 0);
                            matched = Some((uid, 0usize, block_idx, None, used_revcomp));
                            special_note = Some("special_abundant_match");
                            break;
                        }
                        special_note = Some("special_abundant_miss");
                        continue;
                    }
                    let (unitig_id, rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
                    if is_km {
                        let uid = unitig_id as usize;
                        if uid >= index.km_unitigs.len() {
                            continue;
                        }
                        let km_pos = rel_pos as usize;
                        let km_seq = index.km_unitigs[uid].as_slice();
                        let (used_revcomp, min_pos_match) = if km_seq == kmer {
                            (false, min_pos_fwd)
                        } else if km_seq == rev_buf.as_slice() {
                            (true, min_pos_rev)
                        } else {
                            continue;
                        };
                        if !kallisto_bifrost_find
                            && min_pos_match != km_pos
                            && min_pos_match + km_pos != diff
                        {
                            continue;
                        }
                        if (used_revcomp && !allow_rev_base)
                            || (!used_revcomp && !allow_forward_base)
                        {
                            continue;
                        }
                        let unitig_id = index.unitigs.len() + uid;
                        let block_idx =
                            block_index_for_position(index.ec_blocks[unitig_id].as_slice(), 0);
                        matched = Some((unitig_id, 0usize, block_idx, Some(true), used_revcomp));
                        break;
                    } else {
                        let uid = unitig_id as usize;
                        if uid >= index.unitigs.len() {
                            continue;
                        }
                        let rel_pos = rel_pos as isize;
                        let start_fwd = rel_pos - min_pos_fwd as isize;
                        if allow_forward_base && start_fwd >= 0 {
                            let start = start_fwd as usize;
                            if start + index.k <= index.unitigs[uid].len()
                                && &index.unitigs[uid][start..start + index.k] == kmer
                            {
                                let block_idx = block_index_for_position(
                                    index.ec_blocks[uid].as_slice(),
                                    start,
                                );
                                matched = Some((uid, start, block_idx, Some(false), false));
                                break;
                            }
                        }
                        let start_rev = rel_pos - diff as isize + min_pos_rev as isize;
                        if allow_rev_base && start_rev >= 0 {
                            let start = start_rev as usize;
                            if start + index.k <= index.unitigs[uid].len()
                                && &index.unitigs[uid][start..start + index.k] == rev_buf.as_slice()
                            {
                                let block_idx = block_index_for_position(
                                    index.ec_blocks[uid].as_slice(),
                                    start,
                                );
                                matched = Some((uid, start, block_idx, Some(false), true));
                                break;
                            }
                        }
                        if kallisto_bifrost_find && special_seen {
                            let unitig = &index.unitigs[uid];
                            if unitig.len() >= index.k {
                                let rel_pos = rel_pos as usize;
                                let start_low = rel_pos.saturating_sub(diff);
                                let start_high = rel_pos.min(unitig.len().saturating_sub(index.k));
                                for start in start_low..=start_high {
                                    if allow_forward_base && &unitig[start..start + index.k] == kmer
                                    {
                                        let block_idx = block_index_for_position(
                                            index.ec_blocks[uid].as_slice(),
                                            start,
                                        );
                                        matched = Some((uid, start, block_idx, Some(false), false));
                                        special_note = Some("match_relaxed");
                                        break;
                                    }
                                    if allow_rev_base
                                        && &unitig[start..start + index.k] == rev_buf.as_slice()
                                    {
                                        let block_idx = block_index_for_position(
                                            index.ec_blocks[uid].as_slice(),
                                            start,
                                        );
                                        matched = Some((uid, start, block_idx, Some(false), true));
                                        special_note = Some("match_relaxed");
                                        break;
                                    }
                                }
                            }
                            if matched.is_some() {
                                break;
                            }
                        }
                    }
                }

                let (match_unitig_id, match_unitig_pos, match_block_idx, is_km, note, match_rev) =
                    if let Some((uid, start, block_idx, is_km, match_rev)) = matched {
                        (
                            Some(uid),
                            Some(start),
                            block_idx,
                            is_km,
                            "match",
                            Some(match_rev),
                        )
                    } else if skip_minimizer {
                        (None, None, None, None, "special_overcrowded_skip", None)
                    } else if let Some(note) = special_note {
                        (None, None, None, None, note, None)
                    } else {
                        (None, None, None, None, "no_match", None)
                    };
                let used_revcomp = match_rev.unwrap_or(use_revcomp);
                let mut note_text = note.to_string();
                if fallback_active {
                    note_text.push_str("|minhash_fallback");
                }

                out.push(MinimizerDebugRow {
                    read_pos: pos,
                    kmer: kmer_text.clone(),
                    kmer_canon: kmer_canon_text.clone(),
                    used_revcomp,
                    min_pos: Some(min_pos),
                    min_seq,
                    mphf_key,
                    mphf_level: Some(mphf_debug.level),
                    mphf_hash_raw: Some(mphf_debug.hash_raw),
                    mphf_hash_domain: Some(mphf_debug.hash_domain),
                    mphf_bucket: mphf_debug.bucket,
                    mphf_rank: mphf_debug.rank,
                    mphf_final_hash_hit: Some(mphf_debug.final_hash_hit),
                    mphf_hit: true,
                    positions: Some(positions.len()),
                    sample_positions: Some(format_positions_sample(positions)),
                    matched: match_unitig_id.is_some(),
                    match_unitig_id,
                    match_unitig_pos,
                    match_block_idx,
                    is_km,
                    special_code,
                    special_in_dlist,
                    special_in_hmap,
                    special_uid,
                    note: note_text,
                });
                if skip_minimizer {
                    break;
                }
            }
            if !kallisto_bifrost_find {
                break;
            }
            if skip_minimizer && let Some((bytes, pos)) = fallback_min.take() {
                min_candidates = vec![(bytes, pos)];
                fallback_active = true;
                continue 'min_loop;
            }
            break;
        }
    }
    out
}

pub fn trace_read_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> ReadTraceResult {
    let (_read_ec, trace) = trace_read_bifrost_inner(index, seq, strand, filter, options);
    trace
}

pub fn trace_read_bifrost_with_hits(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> (ReadTraceResult, Vec<HitTrace>) {
    let (read_ec, trace) = trace_read_bifrost_inner(index, seq, strand, filter, options);
    let mut read_ec = read_ec;
    if read_ec.is_none() && trace.reason == Some(DebugFailReason::IntersectionEmpty) {
        let mut union_options = options;
        union_options.do_union = true;
        union_options.min_range = 1;
        read_ec = ec_for_read_bifrost(index, seq, strand, None, union_options);
    }
    let Some(read_ec) = read_ec else {
        return (trace, Vec::new());
    };
    let use_shade = index.use_shade && !options.do_union;
    let mut shade_scratch = Vec::new();
    let mut hits = Vec::with_capacity(read_ec.hits.len());
    for hit in &read_ec.hits {
        let Some(blocks) = index.ec_blocks.get(hit.unitig_id) else {
            continue;
        };
        let Some(block) = blocks.get(hit.block_idx) else {
            continue;
        };
        let raw_ec = &block.ec;
        let ec_vec = if use_shade {
            filter_shades(raw_ec, &index.shade_sequences, &mut shade_scratch);
            shade_scratch.clone()
        } else {
            raw_ec.to_vec()
        };
        let offlist = if let Some(onlist) = index.onlist.as_deref() {
            ec_has_offlist(&ec_vec, onlist)
        } else {
            false
        };
        hits.push(HitTrace {
            unitig_id: hit.unitig_id,
            block_idx: hit.block_idx,
            read_pos: hit.read_pos,
            used_revcomp: hit.used_revcomp,
            ec: ec_vec,
            offlist,
        });
    }
    (trace, hits)
}

fn trace_read_bifrost_inner(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> (Option<ReadEc>, ReadTraceResult) {
    let mut dbg = ReadDebugState::default();
    let read_ec = ec_for_read_bifrost(index, seq, strand, Some(&mut dbg), options);

    if let Some(read_ec) = read_ec {
        let mut read_ec = read_ec;
        let ec_before = read_ec.ec.clone();
        let mut ec_fragment = ec_before.clone();
        if let Some(filter) = filter
            && !filter.single_overhang
            && filter.fragment_length > 0
            && let Some(best_match) = read_ec.best_match
        {
            ec_fragment = filter_ec_by_fragment(
                index,
                &ec_fragment,
                best_match,
                filter.fragment_length as i64,
            );
        }
        read_ec.ec = ec_fragment.clone();
        let mut ec_strand = ec_fragment.clone();
        if let Some(mode) = options.strand_specific {
            let comprehensive = options.do_union || options.no_jump || index.use_shade;
            if let Some(filtered) =
                apply_strand_filter(index, &read_ec, mode, true, comprehensive, &mut None, b"")
            {
                ec_strand = filtered;
            } else {
                ec_strand.clear();
            }
        }
        let reason = if ec_strand.is_empty() {
            Some(DebugFailReason::Unknown)
        } else {
            None
        };
        return (
            Some(read_ec),
            ReadTraceResult {
                ec_before_filters: Some(ec_before),
                ec_after_fragment_filter: Some(ec_fragment),
                ec_after_strand_filter: Some(ec_strand),
                reason,
                kmer_pos: None,
                min_pos: None,
                kmer_seq: None,
                min_seq: None,
                positions: None,
                sample_positions: None,
                used_revcomp: dbg.used_revcomp,
                positions_visited: Some(dbg.visited_positions.clone()),
            },
        );
    }

    let (reason, kmer_pos, min_pos) = debug_reason(&dbg);
    let kmer_seq = kmer_pos.and_then(|pos| {
        if pos + index.k <= seq.len() {
            Some(String::from_utf8_lossy(&seq[pos..pos + index.k]).into_owned())
        } else {
            None
        }
    });
    let min_seq = match (kmer_pos, min_pos) {
        (Some(kp), Some(mp)) => {
            let start = kp + mp;
            let end = start + index.g;
            if end <= seq.len() {
                Some(String::from_utf8_lossy(&seq[start..end]).into_owned())
            } else {
                None
            }
        }
        _ => None,
    };
    let (positions, sample_positions) = dbg
        .first_no_match_positions
        .as_ref()
        .map(|v| (Some(v.2), Some(v.3.clone())))
        .unwrap_or((None, None));
    (
        None,
        ReadTraceResult {
            ec_before_filters: None,
            ec_after_fragment_filter: None,
            ec_after_strand_filter: None,
            reason: Some(reason),
            kmer_pos,
            min_pos,
            kmer_seq,
            min_seq,
            positions,
            sample_positions,
            used_revcomp: dbg.used_revcomp,
            positions_visited: Some(dbg.visited_positions.clone()),
        },
    )
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
    pub kallisto_enum: bool,
    pub kallisto_strict: bool,
    pub kallisto_local_fallback: bool,
    pub kallisto_fallback: bool,
    pub discard_special_only: bool,
    pub skip_overcrowded_minimizer: bool,
    pub kallisto_direct_kmer: bool,
    pub kallisto_bifrost_find: bool,
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
            kallisto_enum: false,
            kallisto_strict: false,
            kallisto_local_fallback: false,
            kallisto_fallback: false,
            discard_special_only: false,
            skip_overcrowded_minimizer: false,
            kallisto_direct_kmer: false,
            kallisto_bifrost_find: false,
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
