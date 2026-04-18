use std::collections::HashMap;

use crate::bias::BiasCounts;
use crate::io::{PackedFastqBatch, ReadSource};
use crate::{Error, Result};

use super::{
    BifrostIndex, DebugFailReason, DebugReport, EcCounts, FragmentLengthStats, KmerEcIndex,
    PairedTraceResult, PseudoalignOptions, ReadDebugState, Strand,
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
                if std::env::var_os("KALLISTORS_RESET_ALL_CACHES_PER_READ").is_some() {
                    super::reset_thread_local_caches();
                }
                reads_processed += 1;
                let ec1 = ec_for_read_naive(index, &a.seq);
                let ec2 = ec_for_read_naive(index, &b.seq);
                let mut merged = Vec::new();
                match (ec1.as_ref(), ec2.as_ref()) {
                    (None, None) => continue,
                    (Some(ec), None) | (None, Some(ec)) => merged.extend_from_slice(ec),
                    (Some(left), Some(right)) => {
                        super::intersect_sorted(left, right, &mut merged);
                    }
                }
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

pub fn trace_read_pair_bifrost(
    index: &BifrostIndex,
    left_seq: &[u8],
    right_seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> PairedTraceResult {
    let mut dbg1 = ReadDebugState::default();
    let mut dbg2 = ReadDebugState::default();
    let mut ec1 = super::ec_for_read_bifrost(index, left_seq, strand, Some(&mut dbg1), options);
    let mut ec2 = super::ec_for_read_bifrost(index, right_seq, strand, Some(&mut dbg2), options);

    if let Some(mode) = options.strand_specific {
        let comprehensive = options.do_union || options.no_jump || index.use_shade;
        if let Some(read_ec) = ec1.as_mut()
            && let Some(filtered) = super::apply_strand_filter(
                index,
                read_ec,
                mode,
                true,
                comprehensive,
                &mut None,
                b"",
            )
        {
            read_ec.ec = filtered;
        }
        if let Some(read_ec) = ec2.as_mut()
            && let Some(filtered) = super::apply_strand_filter(
                index,
                read_ec,
                mode,
                false,
                comprehensive,
                &mut None,
                b"",
            )
        {
            read_ec.ec = filtered;
        }
    }

    let mut left_reason_override = None;
    let mut right_reason_override = None;
    if ec1.as_ref().is_some_and(|v| v.ec.is_empty())
        && ec1.as_ref().is_some_and(|v| !v.hard_reject_pair)
    {
        left_reason_override = Some(DebugFailReason::EmptyEc);
        ec1 = None;
    }
    if ec2.as_ref().is_some_and(|v| v.ec.is_empty())
        && ec2.as_ref().is_some_and(|v| !v.hard_reject_pair)
    {
        right_reason_override = Some(DebugFailReason::EmptyEc);
        ec2 = None;
    }

    let hard_reject_pair = ec1.as_ref().is_some_and(|v| v.hard_reject_pair)
        || ec2.as_ref().is_some_and(|v| v.hard_reject_pair);

    let mut merged = Vec::new();
    let mut had_offlist = false;
    let merged_reason = if hard_reject_pair {
        Some(DebugFailReason::IntersectionEmpty)
    } else {
        match (ec1.as_ref(), ec2.as_ref()) {
            (None, None) => {
                if left_reason_override.is_none() && right_reason_override.is_none() {
                    left_reason_override = Some(DebugFailReason::NoKmerMatch);
                    right_reason_override = Some(DebugFailReason::NoKmerMatch);
                }
                None
            }
            (Some(ec), None) | (None, Some(ec)) => {
                merged.extend_from_slice(&ec.ec);
                had_offlist = ec.had_offlist;
                None
            }
            (Some(left), Some(right)) => {
                super::intersect_sorted(&left.ec, &right.ec, &mut merged);
                had_offlist = left.had_offlist || right.had_offlist;
                if merged.is_empty() {
                    Some(DebugFailReason::IntersectionEmpty)
                } else {
                    None
                }
            }
        }
    };

    let left_before = ec1.as_ref().map(|v| v.ec.clone());
    let right_before = ec2.as_ref().map(|v| v.ec.clone());
    let left_trace = super::trace_result_from_debug_state(
        index,
        left_seq,
        &dbg1,
        ec1.is_some(),
        left_before.clone(),
        left_before.clone(),
        left_before,
        left_reason_override,
    );
    let right_trace = super::trace_result_from_debug_state(
        index,
        right_seq,
        &dbg2,
        ec2.is_some(),
        right_before.clone(),
        right_before.clone(),
        right_before,
        right_reason_override,
    );

    PairedTraceResult {
        left: left_trace,
        right: right_trace,
        merged_ec: (!merged.is_empty()).then_some(merged),
        merged_reason,
        hard_reject_pair,
        had_offlist,
    }
}

pub(crate) fn pseudoalign_paired_bifrost_batch_into(
    index: &BifrostIndex,
    left_batch: PackedFastqBatch,
    right_batch: PackedFastqBatch,
    strand: Strand,
    options: PseudoalignOptions,
    counts: &mut EcCounts,
    ec_map: &mut HashMap<Vec<u32>, usize>,
) -> Result<()> {
    if left_batch.len() != right_batch.len() {
        return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()));
    }
    let mut frag_stats = counts.fragment_length_stats.take().unwrap_or_default();
    let mut frag_hist = counts
        .fragment_length_hist
        .take()
        .unwrap_or_else(|| vec![0u32; super::MAX_FRAG_LEN as usize]);

    for (a, b) in left_batch.records().zip(right_batch.records()) {
        if std::env::var_os("KALLISTORS_RESET_ALL_CACHES_PER_READ").is_some() {
            super::reset_thread_local_caches();
        }
        counts.reads_processed = counts.reads_processed.saturating_add(1);

        let mut ec1 = super::ec_for_read_bifrost(index, a.seq, strand, None, options);
        let mut ec2 = super::ec_for_read_bifrost(index, b.seq, strand, None, options);
        if ec1.is_none() && ec2.is_none() {
            continue;
        }
        if let Some(mode) = options.strand_specific {
            let comprehensive = options.do_union || options.no_jump || index.use_shade;
            if let Some(read_ec) = ec1.as_mut()
                && let Some(filtered) = super::apply_strand_filter(
                    index,
                    read_ec,
                    mode,
                    true,
                    comprehensive,
                    &mut None,
                    b"",
                )
            {
                read_ec.ec = filtered;
            }
            if let Some(read_ec) = ec2.as_mut()
                && let Some(filtered) = super::apply_strand_filter(
                    index,
                    read_ec,
                    mode,
                    false,
                    comprehensive,
                    &mut None,
                    b"",
                )
            {
                read_ec.ec = filtered;
            }
        }
        if ec1.as_ref().is_some_and(|v| v.ec.is_empty())
            && ec1.as_ref().is_some_and(|v| !v.hard_reject_pair)
        {
            ec1 = None;
        }
        if ec2.as_ref().is_some_and(|v| v.ec.is_empty())
            && ec2.as_ref().is_some_and(|v| !v.hard_reject_pair)
        {
            ec2 = None;
        }
        if ec1.as_ref().is_some_and(|v| v.hard_reject_pair)
            || ec2.as_ref().is_some_and(|v| v.hard_reject_pair)
        {
            continue;
        }

        let use_shade = index.use_shade;
        if use_shade && options.do_union {
            if let Some(read_ec) = ec1.as_mut() {
                super::filter_shades_in_place(&mut read_ec.ec, &index.shade_sequences);
            }
            if let Some(read_ec) = ec2.as_mut() {
                super::filter_shades_in_place(&mut read_ec.ec, &index.shade_sequences);
            }
        }

        let mut merged = Vec::new();
        match (ec1.as_ref(), ec2.as_ref()) {
            (None, None) => continue,
            (Some(ec), None) | (None, Some(ec)) => merged.extend_from_slice(&ec.ec),
            (Some(left), Some(right)) => super::intersect_sorted(&left.ec, &right.ec, &mut merged),
        }
        if options.dfk_onlist
            && (ec1.as_ref().is_some_and(|v| v.had_offlist)
                || ec2.as_ref().is_some_and(|v| v.had_offlist))
            && !merged.is_empty()
            && let Some(onlist) = index.onlist.as_deref()
        {
            let dummy = onlist.len() as u32;
            if merged.last().copied() != Some(dummy) {
                merged.push(dummy);
            }
        }
        if merged.is_empty() {
            continue;
        }

        if use_shade && !merged.is_empty() {
            let mut shade_candidates = Vec::new();
            if let Some(read_ec) = ec1.as_ref() {
                super::merge_sorted_unique_vec(&mut shade_candidates, &read_ec.shade_union);
            }
            if let Some(read_ec) = ec2.as_ref() {
                super::merge_sorted_unique_vec(&mut shade_candidates, &read_ec.shade_union);
            }
            if !shade_candidates.is_empty() {
                for shade in shade_candidates {
                    let Some(color) = index.shade_to_color.get(shade as usize).copied().flatten()
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
                && ec1.is_some()
                && ec2.is_some()
                && let Some(frag_len) = estimate_fragment_length_for_pair(
                    index,
                    tr,
                    ec1.as_ref().expect("checked"),
                    ec2.as_ref().expect("checked"),
                )
            {
                let idx = frag_len as usize;
                if idx < frag_hist.len() {
                    frag_hist[idx] = frag_hist[idx].saturating_add(1);
                    frag_stats.add(frag_len as f64);
                }
            }
        }

        if let Some(bias_counts) = counts.bias.as_mut()
            && bias_counts.total < options.max_bias as u64
            && let Some(best_match) = ec1
                .as_ref()
                .and_then(|v| v.best_match)
                .or_else(|| ec2.as_ref().and_then(|v| v.best_match))
            && let Some(hex) = super::bias_hexamer_for_match(index, best_match)
        {
            bias_counts.record(hex);
        }

        counts.reads_aligned = counts.reads_aligned.saturating_add(1);
        super::add_ec_count(counts, ec_map, merged);
    }
    counts.fragment_length_stats = (frag_stats.count() > 0).then_some(frag_stats);
    counts.fragment_length_hist = counts.fragment_length_stats.as_ref().map(|_| frag_hist);
    Ok(())
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

                if ec1.is_none() && ec2.is_none() {
                    if let Some(r) = report.as_deref_mut() {
                        let (state, header) = if dbg1
                            .as_ref()
                            .and_then(|s| s.first_no_match_positions.as_ref())
                            .is_some()
                        {
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
                let mut ec1 = ec1;
                let mut ec2 = ec2;
                if let Some(mode) = options.strand_specific {
                    let comprehensive = options.do_union || options.no_jump || index.use_shade;
                    if let Some(read_ec) = ec1.as_mut()
                        && let Some(filtered) = super::apply_strand_filter(
                            index,
                            read_ec,
                            mode,
                            true,
                            comprehensive,
                            &mut report,
                            &a.header,
                        )
                    {
                        read_ec.ec = filtered;
                    }
                    if let Some(read_ec) = ec2.as_mut()
                        && let Some(filtered) = super::apply_strand_filter(
                            index,
                            read_ec,
                            mode,
                            false,
                            comprehensive,
                            &mut report,
                            &b.header,
                        )
                    {
                        read_ec.ec = filtered;
                    }
                }
                if ec1.as_ref().is_some_and(|v| v.ec.is_empty())
                    && ec1.as_ref().is_some_and(|v| !v.hard_reject_pair)
                {
                    ec1 = None;
                }
                if ec2.as_ref().is_some_and(|v| v.ec.is_empty())
                    && ec2.as_ref().is_some_and(|v| !v.hard_reject_pair)
                {
                    ec2 = None;
                }
                if ec1.as_ref().is_some_and(|v| v.hard_reject_pair)
                    || ec2.as_ref().is_some_and(|v| v.hard_reject_pair)
                {
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

                let use_shade = index.use_shade;
                if use_shade && options.do_union {
                    if let Some(read_ec) = ec1.as_mut() {
                        super::filter_shades_in_place(&mut read_ec.ec, &index.shade_sequences);
                    }
                    if let Some(read_ec) = ec2.as_mut() {
                        super::filter_shades_in_place(&mut read_ec.ec, &index.shade_sequences);
                    }
                }

                let mut merged = Vec::new();
                match (ec1.as_ref(), ec2.as_ref()) {
                    (None, None) => continue,
                    (Some(ec), None) | (None, Some(ec)) => merged.extend_from_slice(&ec.ec),
                    (Some(left), Some(right)) => {
                        super::intersect_sorted(&left.ec, &right.ec, &mut merged);
                    }
                }
                if options.dfk_onlist
                    && (ec1.as_ref().is_some_and(|v| v.had_offlist)
                        || ec2.as_ref().is_some_and(|v| v.had_offlist))
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
                    if let Some(read_ec) = ec1.as_ref() {
                        super::merge_sorted_unique_vec(&mut shade_candidates, &read_ec.shade_union);
                    }
                    if let Some(read_ec) = ec2.as_ref() {
                        super::merge_sorted_unique_vec(&mut shade_candidates, &read_ec.shade_union);
                    }
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
                        && ec1.is_some()
                        && ec2.is_some()
                        && let Some(frag_len) = estimate_fragment_length_for_pair(
                            index,
                            tr,
                            ec1.as_ref().expect("checked"),
                            ec2.as_ref().expect("checked"),
                        )
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
                    && let Some(best_match) = ec1
                        .as_ref()
                        .and_then(|v| v.best_match)
                        .or_else(|| ec2.as_ref().and_then(|v| v.best_match))
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
    let block_left = index
        .flat_ec
        .block_index_for_position(left_match.unitig_id, left_match.unitig_pos)?;
    let block_right = index
        .flat_ec
        .block_index_for_position(right_match.unitig_id, right_match.unitig_pos)?;
    if block_left != block_right {
        return None;
    }
    let (_, ub_left) = index
        .flat_ec
        .block_bounds(left_match.unitig_id, block_left)?;
    let (_, ub_right) = index
        .flat_ec
        .block_bounds(right_match.unitig_id, block_right)?;
    if ub_left != ub_right {
        return None;
    }
    if left_match.used_revcomp == right_match.used_revcomp {
        return None;
    }
    if index
        .flat_ec
        .ec(left_match.unitig_id, block_left)
        .binary_search(&tr)
        .is_err()
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
