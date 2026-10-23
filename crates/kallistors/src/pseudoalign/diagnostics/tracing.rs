//! Single-read trace construction and placement reporting.

use super::super::debug::{ReadDebugState, debug_reason};
use super::super::{
    BifrostIndex, DebugFailReason, FragmentFilter, HitTrace, MatchInfo, PseudoalignOptions,
    PseudoalignPlacement, ReadEc, ReadTraceResult, Strand, apply_strand_filter,
    ec_for_read_bifrost, ec_has_offlist, filter_ec_by_fragment, filter_shades,
};

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
        let Some(blocks) = index.flat_ec.owned_ec_blocks.get(hit.unitig_id) else {
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn trace_result_from_debug_state(
    index: &BifrostIndex,
    seq: &[u8],
    dbg: &ReadDebugState,
    read_ec_present: bool,
    ec_before_filters: Option<Vec<u32>>,
    ec_after_fragment_filter: Option<Vec<u32>>,
    ec_after_strand_filter: Option<Vec<u32>>,
    reason_override: Option<DebugFailReason>,
) -> ReadTraceResult {
    if let Some(reason) = reason_override {
        return ReadTraceResult {
            ec_before_filters,
            ec_after_fragment_filter,
            ec_after_strand_filter,
            reason: Some(reason),
            kmer_pos: None,
            min_pos: None,
            kmer_seq: None,
            min_seq: None,
            positions: None,
            sample_positions: None,
            used_revcomp: dbg.used_revcomp,
            positions_visited: Some(dbg.visited_positions.clone()),
            dropped_hits: Some(dbg.dropped_hits.clone()),
            minimizer_candidates: Some(dbg.minimizer_candidates.clone()),
            jump_decisions: Some(dbg.jump_decisions.clone()),
            placement: None,
        };
    }
    if read_ec_present {
        return ReadTraceResult {
            ec_before_filters,
            ec_after_fragment_filter,
            ec_after_strand_filter,
            reason: None,
            kmer_pos: None,
            min_pos: None,
            kmer_seq: None,
            min_seq: None,
            positions: None,
            sample_positions: None,
            used_revcomp: dbg.used_revcomp,
            positions_visited: Some(dbg.visited_positions.clone()),
            dropped_hits: Some(dbg.dropped_hits.clone()),
            minimizer_candidates: Some(dbg.minimizer_candidates.clone()),
            jump_decisions: Some(dbg.jump_decisions.clone()),
            placement: None,
        };
    }

    let (reason, kmer_pos, min_pos) = debug_reason(dbg);
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
    ReadTraceResult {
        ec_before_filters,
        ec_after_fragment_filter,
        ec_after_strand_filter,
        reason: Some(reason),
        kmer_pos,
        min_pos,
        kmer_seq,
        min_seq,
        positions,
        sample_positions,
        used_revcomp: dbg.used_revcomp,
        positions_visited: Some(dbg.visited_positions.clone()),
        dropped_hits: Some(dbg.dropped_hits.clone()),
        minimizer_candidates: Some(dbg.minimizer_candidates.clone()),
        jump_decisions: Some(dbg.jump_decisions.clone()),
        placement: None,
    }
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
        let mut trace = trace_result_from_debug_state(
            index,
            seq,
            &dbg,
            true,
            Some(ec_before),
            Some(ec_fragment),
            Some(ec_strand),
            reason,
        );
        trace.placement = read_ec
            .best_match
            .and_then(|best_match| placement_for_match(&read_ec.ec, best_match, seq.len()));
        return (Some(read_ec), trace);
    }

    (
        None,
        trace_result_from_debug_state(index, seq, &dbg, false, None, None, None, None),
    )
}

pub(in crate::pseudoalign) fn placement_for_match(
    ec: &[u32],
    best_match: MatchInfo,
    read_len: usize,
) -> Option<PseudoalignPlacement> {
    let transcript_id = *ec.first()?;
    let start = if best_match.used_revcomp {
        best_match
            .unitig_pos
            .saturating_add(read_len)
            .saturating_sub(best_match.read_pos)
            .saturating_sub(1)
    } else {
        best_match.unitig_pos.saturating_sub(best_match.read_pos)
    };
    Some(PseudoalignPlacement {
        transcript_id,
        start,
        reverse: best_match.used_revcomp,
    })
}
