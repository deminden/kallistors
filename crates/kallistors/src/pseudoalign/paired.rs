pub use super::diagnostics::trace_read_pair_bifrost;
use super::diagnostics::{PairedFinalTraceRow, PairedFinalTraceWriter};
use super::utils::EcCountMap;

use crate::bias::BiasCounts;
use crate::io::{PackedSeqBatch, ReadSource};
use crate::{Error, Result};

use super::{
    BifrostIndex, DebugFailReason, DebugReport, EcCounts, FragmentLengthStats, KmerEcIndex,
    PseudoalignOptions, ReadDebugState, Strand,
};

/// Upstream-style observation goal shared by paired workers. Each batch snapshots
/// the remaining quota, then publishes its observations when that batch finishes.
/// Concurrent batches can overshoot the goal, as in original kallisto.
pub(crate) struct FragmentLengthSampler {
    goal: usize,
    observations: std::sync::atomic::AtomicUsize,
}

impl FragmentLengthSampler {
    pub(crate) fn new(goal: usize) -> Self {
        Self {
            goal,
            observations: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    pub(crate) fn remaining(&self) -> usize {
        self.goal
            .saturating_sub(self.observations.load(std::sync::atomic::Ordering::Relaxed))
    }

    pub(crate) fn submit(&self, count: usize) {
        self.observations
            .fetch_add(count, std::sync::atomic::Ordering::Relaxed);
    }
}

pub fn pseudoalign_paired_naive<R1: ReadSource, R2: ReadSource>(
    index: &KmerEcIndex,
    reader1: &mut R1,
    reader2: &mut R2,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: EcCountMap = EcCountMap::default();
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
                if super::reset_all_caches_per_read() {
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
                let ec_id = match ec_map.get(&merged, &ec_list) {
                    Some(id) => id,
                    None => {
                        let id = ec_list.len();
                        ec_map.insert(&merged, id);
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn pseudoalign_paired_bifrost_batch_into<'index>(
    index: &'index BifrostIndex,
    left_batch: &PackedSeqBatch,
    right_batch: &PackedSeqBatch,
    read_range: std::ops::Range<usize>,
    strand: Strand,
    options: PseudoalignOptions,
    counts: &mut EcCounts,
    ec_map: &mut EcCountMap,
    sampling: Option<&FragmentLengthSampler>,
    pairs: &mut super::core::ScopedEcPairs<'index>,
) -> Result<()> {
    if left_batch.len() != right_batch.len() {
        return Err(Error::InvalidFormat("paired FASTQ length mismatch".into()));
    }
    if read_range.start > read_range.end || read_range.end > left_batch.len() {
        return Err(Error::InvalidFormat(
            "paired FASTQ batch range out of bounds".into(),
        ));
    }
    let mut frag_stats = counts.fragment_length_stats.take().unwrap_or_default();
    let mut frag_hist = counts
        .fragment_length_hist
        .take()
        .unwrap_or_else(|| vec![0u32; super::MAX_FRAG_LEN as usize]);

    let fragment_budget = sampling.map_or(usize::MAX, FragmentLengthSampler::remaining);
    if fragment_budget == 0
        && !super::bias_sampling_pending(counts.bias.as_ref(), options.max_bias)
        && super::counting_filters_supported(index, options)
        && !index.use_shade
        && !options.dfk_onlist
    {
        count_paired_borrowed_batch(
            index,
            left_batch,
            right_batch,
            read_range,
            strand,
            options,
            counts,
            ec_map,
            pairs,
        );
        if let Some(sampler) = sampling {
            sampler.submit(0);
        }
        counts.fragment_length_stats = (frag_stats.count() > 0).then_some(frag_stats);
        counts.fragment_length_hist = counts.fragment_length_stats.as_ref().map(|_| frag_hist);
        return Ok(());
    }

    let mut sampled = 0usize;
    let mut merged = Vec::new();
    for idx in read_range {
        let a_seq = unsafe { left_batch.seq_unchecked(idx) };
        let b_seq = unsafe { right_batch.seq_unchecked(idx) };
        if super::reset_all_caches_per_read() {
            super::reset_thread_local_caches();
        }
        counts.reads_processed = counts.reads_processed.saturating_add(1);

        let count_only = fragment_budget == 0
            && !super::bias_sampling_pending(counts.bias.as_ref(), options.max_bias)
            && super::counting_filters_supported(index, options);
        let scan = |seq: &[u8]| {
            if count_only {
                super::core::ec_for_read_bifrost_counting(index, seq, strand, options)
            } else {
                super::ec_for_read_bifrost(index, seq, strand, None, options)
            }
        };
        let mut ec1 = scan(a_seq);
        // A hard rejection of either mate always rejects the pair. Both FASTQ
        // records have already been parsed and this pair is already counted.
        if ec1.as_ref().is_some_and(|ec| ec.hard_reject_pair) {
            continue;
        }
        let mut ec2 = scan(b_seq);
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
                read_ec.hard_reject_pair |= filtered.is_empty();
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
                read_ec.hard_reject_pair |= filtered.is_empty();
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

        merged.clear();
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

        if sampled < fragment_budget && merged.len() == 1 {
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
                    sampled += 1;
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
        super::utils::add_ec_count_slice(counts, ec_map, &merged);
    }
    if let Some(sampler) = sampling {
        sampler.submit(sampled);
    }
    counts.fragment_length_stats = (frag_stats.count() > 0).then_some(frag_stats);
    counts.fragment_length_hist = counts.fragment_length_stats.as_ref().map(|_| frag_hist);
    Ok(())
}

// This caller has validated the paired range and exhausted the placement quota.
// The index stays immutably borrowed while both mate EC views are consumed.
#[allow(clippy::too_many_arguments)]
fn count_paired_borrowed_batch<'index>(
    index: &'index BifrostIndex,
    left: &PackedSeqBatch,
    right: &PackedSeqBatch,
    range: std::ops::Range<usize>,
    strand: Strand,
    options: PseudoalignOptions,
    counts: &mut EcCounts,
    ec_map: &mut EcCountMap,
    pairs: &mut super::core::ScopedEcPairs<'index>,
) {
    let mut merged = Vec::new();
    let mut count = |a: Option<super::BorrowedReadEc<'_>>, b: Option<super::BorrowedReadEc<'_>>| {
        counts.reads_processed = counts.reads_processed.saturating_add(1);
        if a.as_ref().is_some_and(|read| read.read.hard_reject_pair)
            || b.as_ref().is_some_and(|read| read.read.hard_reject_pair)
        {
            return;
        }
        let a = a.as_ref().filter(|read| !read.ec().is_empty());
        let b = b.as_ref().filter(|read| !read.ec().is_empty());
        let ec = match (a, b) {
            (None, None) => return,
            (Some(read), None) | (None, Some(read)) => read.ec(),
            (Some(a), Some(b)) if std::ptr::eq(a.ec(), b.ec()) => a.ec(),
            (Some(a), Some(b)) => {
                super::intersect_sorted(a.ec(), b.ec(), &mut merged);
                &merged
            }
        };
        if ec.is_empty() {
            return;
        }
        counts.reads_aligned = counts.reads_aligned.saturating_add(1);
        super::utils::add_ec_count_slice(counts, ec_map, ec);
    };
    if super::core::try_paired_counting_batch(
        index,
        left,
        right,
        range.clone(),
        strand,
        options,
        pairs,
        &mut count,
    ) {
        return;
    }
    for idx in range {
        if super::reset_all_caches_per_read() {
            super::reset_thread_local_caches();
        }
        // SAFETY: caller proved both batch lengths equal and this range in bounds.
        let a_seq = unsafe { left.seq_unchecked(idx) };
        let b_seq = unsafe { right.seq_unchecked(idx) };
        let a = super::core::ec_for_read_bifrost_counting_borrowed(index, a_seq, strand, options)
            .map(|read| {
                if let Some(mode) = options.strand_specific {
                    read.filter_strand(index, mode, true)
                } else {
                    read
                }
            });
        if a.as_ref().is_some_and(|read| read.read.hard_reject_pair) {
            count(a, None);
            continue;
        }
        let b = super::core::ec_for_read_bifrost_counting_borrowed(index, b_seq, strand, options)
            .map(|read| {
                if let Some(mode) = options.strand_specific {
                    read.filter_strand(index, mode, false)
                } else {
                    read
                }
            });
        count(a, b);
    }
}

fn pseudoalign_paired_bifrost_inner<R1: ReadSource, R2: ReadSource>(
    index: &BifrostIndex,
    reader1: &mut R1,
    reader2: &mut R2,
    strand: Strand,
    mut report: Option<&mut DebugReport>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    let _cache_budget = super::core::MatchCacheBudgetScope::enter(
        options.match_cache_size.allocation(1)?.entries_per_worker(),
    )?;
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: EcCountMap = EcCountMap::default();
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
    let mut final_trace = PairedFinalTraceWriter::from_env()?;

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
                // Once the placement quota is complete, ordinary counting needs
                // the exact EC fold but no hit collection or best-match metadata.
                // Diagnostics, bias and stranded filters retain the complete scan.
                let count_only = options.fragment_length_samples != 0
                    && frag_stats.count() >= options.fragment_length_samples as u64
                    && !super::bias_sampling_pending(bias.as_ref(), options.max_bias)
                    && super::counting_filters_supported(index, options)
                    && report.is_none()
                    && final_trace.is_none();
                let (ec1, ec2) = if count_only {
                    (
                        super::core::ec_for_read_bifrost_counting(index, &a.seq, strand, options),
                        super::core::ec_for_read_bifrost_counting(index, &b.seq, strand, options),
                    )
                } else {
                    (
                        super::ec_for_read_bifrost(index, &a.seq, strand, dbg1.as_mut(), options),
                        super::ec_for_read_bifrost(index, &b.seq, strand, dbg2.as_mut(), options),
                    )
                };

                if ec1.is_none() && ec2.is_none() {
                    if let Some(trace) = final_trace.as_mut() {
                        trace.write(PairedFinalTraceRow {
                            read_id: reads_processed - 1,
                            left: &a,
                            right: &b,
                            merged_ec: None,
                            reason: "no_hits",
                            left_ec: None,
                            right_ec: None,
                        })?;
                    }
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
                        read_ec.hard_reject_pair |= filtered.is_empty();
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
                        read_ec.hard_reject_pair |= filtered.is_empty();
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
                    if let Some(trace) = final_trace.as_mut() {
                        trace.write(PairedFinalTraceRow {
                            read_id: reads_processed - 1,
                            left: &a,
                            right: &b,
                            merged_ec: None,
                            reason: "hard_reject_pair",
                            left_ec: ec1.as_ref().map(|v| v.ec.as_slice()),
                            right_ec: ec2.as_ref().map(|v| v.ec.as_slice()),
                        })?;
                    }
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
                    if let Some(trace) = final_trace.as_mut() {
                        trace.write(PairedFinalTraceRow {
                            read_id: reads_processed - 1,
                            left: &a,
                            right: &b,
                            merged_ec: None,
                            reason: "intersection_empty",
                            left_ec: ec1.as_ref().map(|v| v.ec.as_slice()),
                            right_ec: ec2.as_ref().map(|v| v.ec.as_slice()),
                        })?;
                    }
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

                if (options.fragment_length_samples == 0
                    || frag_stats.count() < options.fragment_length_samples as u64)
                    && merged.len() == 1
                {
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
                if let Some(trace) = final_trace.as_mut() {
                    trace.write(PairedFinalTraceRow {
                        read_id: reads_processed - 1,
                        left: &a,
                        right: &b,
                        merged_ec: Some(merged.as_slice()),
                        reason: "ok",
                        left_ec: ec1.as_ref().map(|v| v.ec.as_slice()),
                        right_ec: ec2.as_ref().map(|v| v.ec.as_slice()),
                    })?;
                }
                let ec_id = match ec_map.get(&merged, &ec_list) {
                    Some(id) => id,
                    None => {
                        let id = ec_list.len();
                        ec_map.insert(&merged, id);
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

pub(super) fn estimate_fragment_length_for_pair(
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

#[cfg(test)]
mod fragment_sampling_tests {
    use super::*;

    #[test]
    fn concurrent_batches_use_remaining_goal_and_may_overshoot() {
        let sampler = FragmentLengthSampler::new(10_000);
        sampler.submit(9_700);
        let first = sampler.remaining();
        let concurrent = sampler.remaining();
        assert_eq!(first, 300);
        assert_eq!(concurrent, 300);
        sampler.submit(first);
        assert_eq!(sampler.remaining(), 0);
        sampler.submit(concurrent);
        assert_eq!(
            sampler
                .observations
                .load(std::sync::atomic::Ordering::Relaxed),
            10_300
        );
        assert_eq!(sampler.remaining(), 0);
    }

    #[test]
    fn short_batches_preserve_unused_quota() {
        let sampler = FragmentLengthSampler::new(10_000);
        sampler.submit(12);
        assert_eq!(sampler.remaining(), 9_988);
        sampler.submit(0);
        assert_eq!(sampler.remaining(), 9_988);
    }
}
