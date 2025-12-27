use std::collections::HashMap;

use crate::bias::BiasCounts;
use crate::io::ReadSource;
use crate::{Error, Result};

use super::{
    BifrostIndex, DebugFailReason, DebugReport, EcCounts, FragmentLengthStats, KmerEcIndex,
    PseudoalignOptions, ReadDebugState, Strand,
};

pub fn pseudoalign_paired_naive<R1: ReadSource, R2: ReadSource>(
    index: &KmerEcIndex,
    reader1: &mut R1,
    reader2: &mut R2,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;

    loop {
        let r1 = reader1.next_record();
        let r2 = reader2.next_record();
        match (r1, r2) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()));
            }
            (Some(a), Some(b)) => {
                let a = a?;
                let b = b?;
                reads_processed += 1;
                let ec1 = ec_for_read_naive(index, &a.seq);
                let ec2 = ec_for_read_naive(index, &b.seq);
                let (Some(ec1), Some(ec2)) = (ec1, ec2) else {
                    continue;
                };
                let mut merged = Vec::new();
                super::intersect_sorted(&ec1, &ec2, &mut merged);
                if merged.is_empty() {
                    continue;
                }
                reads_aligned += 1;
                let ec_id = match ec_map.get(&merged) {
                    Some(id) => *id,
                    None => {
                        let id = ec_list.len();
                        ec_map.insert(merged.clone(), id);
                        ec_list.push(merged.clone());
                        counts.push(0);
                        id
                    }
                };
                counts[ec_id] += 1;
            }
        }
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias: None,
        fragment_length_stats: None,
        fragment_length_hist: None,
    })
}

pub fn pseudoalign_paired_bifrost_with_strand<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
) -> Result<EcCounts> {
    pseudoalign_paired_bifrost_inner(
        index,
        reader1,
        reader2,
        strand,
        None,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_paired_bifrost_debug_with_strand<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    max_traces: usize,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_paired_bifrost_inner(
        index,
        reader1,
        reader2,
        strand,
        Some(&mut report),
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_paired_bifrost_with_options<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    pseudoalign_paired_bifrost_inner(index, reader1, reader2, strand, None, options)
}

pub fn pseudoalign_paired_bifrost_debug_with_options<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    max_traces: usize,
    options: PseudoalignOptions,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_paired_bifrost_inner(
        index,
        reader1,
        reader2,
        strand,
        Some(&mut report),
        options,
    )?;
    Ok((counts, report))
}

fn pseudoalign_paired_bifrost_inner<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    mut report: Option<&mut DebugReport>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;
    let mut frag_stats = FragmentLengthStats::new();
    let mut frag_hist = vec![0u32; super::MAX_FRAG_LEN as usize];
    let mut bias = if options.bias {
        Some(BiasCounts::new())
    } else {
        None
    };

    loop {
        let r1 = reader1.next_record();
        let r2 = reader2.next_record();
        match (r1, r2) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()));
            }
            (Some(a), Some(b)) => {
                let a = a?;
                let b = b?;
                reads_processed += 1;

                let mut dbg1 = report.as_ref().map(|_| ReadDebugState::default());
                let mut dbg2 = report.as_ref().map(|_| ReadDebugState::default());
                let ec1 = super::ec_for_read_bifrost(index, &a.seq, strand, dbg1.as_mut(), options);
                let ec2 = super::ec_for_read_bifrost(index, &b.seq, strand, dbg2.as_mut(), options);

                if ec1.is_none() || ec2.is_none() {
                    if let Some(r) = report.as_deref_mut() {
                        let (state, header) = if ec1.is_none() {
                            (dbg1, &a.header)
                        } else {
                            (dbg2, &b.header)
                        };
                        if let Some(state) = state {
                            let (reason, kmer_pos, min_pos) = super::debug_reason(&state);
                            let (positions, sample_positions) = state
                                .first_no_match_positions
                                .as_ref()
                                .map(|v| (Some(v.2), Some(v.3.clone())))
                                .unwrap_or((None, None));
                            r.record(
                                header,
                                reason,
                                kmer_pos,
                                min_pos,
                                None,
                                None,
                                positions,
                                sample_positions,
                                state.used_revcomp,
                            );
                        }
                    }
                    continue;
                }
                let ec1 = ec1.unwrap();
                let ec2 = ec2.unwrap();
                let mut ec1 = ec1;
                let mut ec2 = ec2;
                if let Some(mode) = options.strand_specific {
                    let comprehensive = options.do_union || options.no_jump || index.use_shade;
                    if let Some(filtered) = super::apply_strand_filter(
                        index,
                        &ec1,
                        mode,
                        true,
                        comprehensive,
                        &mut report,
                        &a.header,
                    ) {
                        ec1.ec = filtered;
                    }
                    if let Some(filtered) = super::apply_strand_filter(
                        index,
                        &ec2,
                        mode,
                        false,
                        comprehensive,
                        &mut report,
                        &b.header,
                    ) {
                        ec2.ec = filtered;
                    }
                }
                if ec1.ec.is_empty() || ec2.ec.is_empty() {
                    continue;
                }

                let use_shade = index.use_shade;
                if use_shade && options.do_union {
                    super::filter_shades_in_place(&mut ec1.ec, &index.shade_sequences);
                    super::filter_shades_in_place(&mut ec2.ec, &index.shade_sequences);
                }

                let mut merged = Vec::new();
                super::intersect_sorted(&ec1.ec, &ec2.ec, &mut merged);
                if options.dfk_onlist
                    && (ec1.had_offlist || ec2.had_offlist)
                    && !merged.is_empty()
                    && let Some(onlist) = index.onlist.as_deref()
                {
                    let dummy = onlist.len() as u32;
                    if merged.last().copied() != Some(dummy) {
                        merged.push(dummy);
                    }
                }
                if merged.is_empty() {
                    if let Some(r) = report.as_deref_mut() {
                        let mut header = Vec::new();
                        header.extend_from_slice(&a.header);
                        header.extend_from_slice(b"|");
                        header.extend_from_slice(&b.header);
                        r.record(
                            &header,
                            DebugFailReason::IntersectionEmpty,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            dbg1.as_ref().map(|d| d.used_revcomp).unwrap_or(false)
                                || dbg2.as_ref().map(|d| d.used_revcomp).unwrap_or(false),
                        );
                    }
                    continue;
                }

                if use_shade && !merged.is_empty() {
                    let mut shade_candidates = Vec::new();
                    super::merge_sorted_unique_vec(&mut shade_candidates, &ec1.shade_union);
                    super::merge_sorted_unique_vec(&mut shade_candidates, &ec2.shade_union);
                    if !shade_candidates.is_empty() {
                        for shade in shade_candidates {
                            let Some(color) =
                                index.shade_to_color.get(shade as usize).copied().flatten()
                            else {
                                continue;
                            };
                            if merged.binary_search(&color).is_ok() {
                                merged.push(shade);
                            }
                        }
                        merged.sort_unstable();
                        merged.dedup();
                    }
                }

                if merged.len() == 1 {
                    let tr = merged[0];
                    let is_shade = index
                        .shade_sequences
                        .get(tr as usize)
                        .copied()
                        .unwrap_or(false);
                    if !is_shade
                        && let Some(frag_len) =
                            estimate_fragment_length_for_pair(index, tr, &ec1, &ec2)
                    {
                        let idx = frag_len as usize;
                        if idx < frag_hist.len() {
                            frag_hist[idx] = frag_hist[idx].saturating_add(1);
                            frag_stats.add(frag_len as f64);
                        }
                    }
                }

                if let Some(bias_counts) = bias.as_mut()
                    && bias_counts.total < options.max_bias as u64
                    && let Some(best_match) = ec1.best_match
                    && let Some(hex) = super::bias_hexamer_for_match(index, best_match)
                {
                    bias_counts.record(hex);
                }

                reads_aligned += 1;
                let ec_id = match ec_map.get(&merged) {
                    Some(id) => *id,
                    None => {
                        let id = ec_list.len();
                        ec_map.insert(merged.clone(), id);
                        ec_list.push(merged.clone());
                        counts.push(0);
                        id
                    }
                };
                counts[ec_id] += 1;
            }
        }
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias,
        fragment_length_stats: (frag_stats.count() > 0).then_some(frag_stats),
        fragment_length_hist: (frag_stats.count() > 0).then_some(frag_hist),
    })
}

fn estimate_fragment_length_for_pair(
    index: &BifrostIndex,
    tr: u32,
    left: &super::ReadEc,
    right: &super::ReadEc,
) -> Option<i64> {
    let left_match = left.best_match?;
    let right_match = right.best_match?;
    if left_match.unitig_id != right_match.unitig_id {
        return None;
    }
    let blocks = index.ec_blocks.get(left_match.unitig_id)?;
    let block_left = super::block_index_for_position(blocks, left_match.unitig_pos)?;
    let block_right = super::block_index_for_position(blocks, right_match.unitig_pos)?;
    if block_left != block_right {
        return None;
    }
    let (_, ub_left) = super::ec_block_at(blocks, left_match.unitig_pos as u32)?;
    let (_, ub_right) = super::ec_block_at(blocks, right_match.unitig_pos as u32)?;
    if ub_left != ub_right {
        return None;
    }
    if left_match.used_revcomp == right_match.used_revcomp {
        return None;
    }
    if !blocks
        .get(block_left)
        .map(|block| block.ec.binary_search(&tr).is_ok())
        .unwrap_or(false)
    {
        return None;
    }
    let k = index.k as i64;
    let p1 = if left_match.used_revcomp {
        left_match.unitig_pos as i64 + k + left_match.read_pos as i64
    } else {
        left_match.unitig_pos as i64 - left_match.read_pos as i64
    };
    let p2 = if right_match.used_revcomp {
        right_match.unitig_pos as i64 + k + right_match.read_pos as i64
    } else {
        right_match.unitig_pos as i64 - right_match.read_pos as i64
    };
    let frag_len = (p1 - p2).abs();
    if frag_len > 0 && frag_len < super::MAX_FRAG_LEN {
        return Some(frag_len);
    }
    None
}

fn ec_for_read_naive(index: &KmerEcIndex, seq: &[u8]) -> Option<Vec<u32>> {
    if seq.len() < index.k {
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut has_hit = false;
    for pos in 0..=seq.len() - index.k {
        if let Some((fwd, rev)) = super::encode_kmer_pair(&seq[pos..pos + index.k])
            && let Some(ec) = super::lookup_ec(index, fwd).or_else(|| super::lookup_ec(index, rev))
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                super::intersect_sorted(&current, ec, &mut next);
                std::mem::swap(&mut current, &mut next);
                if current.is_empty() {
                    break;
                }
            }
        }
    }
    if !has_hit || current.is_empty() {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
    }
    if current.is_empty() {
        None
    } else {
        Some(current)
    }
}
