//! General read scan, accepted-hit stream and jump/backoff decisions.

use super::{finish, probe};
use crate::pseudoalign::core::matcher::{
    jump_distance_for_match, match_kmer_at_pos, match_kmer_direct,
};
use crate::pseudoalign::core::runtime::{
    ACCEPT_RELAXED_STREAM, BACKOFF_DIRECT_RESCUE, block_index_for_position_fast, ec_slice, env_flag,
};
use crate::pseudoalign::debug::{DroppedHitTrace, JumpDecisionTrace, ReadDebugState};
use crate::pseudoalign::ec::{
    encode_kmer_pair, filter_shades, intersect_sorted, valid_kmer_start_at_or_after,
};
use crate::pseudoalign::minimizers::fill_revcomp;
use crate::pseudoalign::{BifrostIndex, Hit, MatchInfo, PseudoalignOptions, ReadEc, Strand};

pub(super) fn scan_read(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    let mut rev_buf: Vec<u8> = Vec::new();
    let allow_forward_base = strand != Strand::Reverse;
    let allow_rev_base = strand != Strand::Forward;
    let mut has_hit = false;
    let diff = index.k.saturating_sub(index.g);
    let mut best_match: Option<MatchInfo> = None;
    let mut first_hit: Option<Hit> = None;
    let mut hits: Vec<Hit> = Vec::new();
    let use_shade = index.use_shade && !options.do_union;
    let mut shade_scratch: Vec<u32> = Vec::new();
    let mut next_shade_scratch: Vec<u32> = Vec::new();
    let mut online_intersection: Vec<u32> = Vec::new();
    let mut has_online_intersection = false;
    let mut saw_special_hit = false;
    let mut saw_non_special_hit = false;
    // Delay the reference-wide count until a synthetic D-list hit needs its ID.
    let dlist_onlist = index.onlist.as_deref();
    let mut any_dlist_kmer = false;
    let mut dlist_early_return = false;
    let mut append_dummy_hit = false;
    // Kallisto falls back to incremental scanning after a failed jump; track that window.
    let mut backoff_until: Option<usize> = None;
    let mut fallback_probe_pending = false;

    let mut pos = 0usize;
    let last_pos = seq.len() - index.k;
    while pos <= last_pos {
        if let Some(state) = dbg.as_deref_mut()
            && state.visited_positions.len() < 4096
        {
            state.visited_positions.push(pos);
        }
        let mut in_backoff = false;
        if let Some(until) = backoff_until {
            if pos > until {
                backoff_until = None;
            } else {
                in_backoff = true;
            }
        }
        let kmer = &seq[pos..pos + index.k];
        let kmer_in_dlist = index.dlist.as_ref().is_some_and(|dlist| {
            encode_kmer_pair(kmer).is_some_and(|(fwd, rev)| dlist.contains(&fwd.min(rev)))
        });
        if kmer_in_dlist {
            any_dlist_kmer = true;
        }
        if !allow_forward_base && !allow_rev_base {
            pos += 1;
            continue;
        }
        if !fill_revcomp(kmer, &mut rev_buf) {
            pos += 1;
            continue;
        }
        let fallback_probe = std::mem::take(&mut fallback_probe_pending);
        let probe::ProbeOutcome::Candidates(probe::ProbeCandidates {
            mut matched,
            min_candidates,
            any_mphf,
            any_positions,
        }) = probe::probe_kmer(
            index,
            kmer,
            &rev_buf,
            pos,
            strand,
            options,
            in_backoff,
            dbg.as_deref_mut(),
        )?
        else {
            pos += 1;
            continue;
        };
        if matched.is_none()
            && !options.kallisto_bifrost_find
            && !options.kallisto_strict
            && !options.kallisto_enum
            && any_mphf
            && any_positions
            && in_backoff
            && min_candidates.len() == 1
            && let Some((uid, start, used_revcomp, _block_idx, matched_relaxed)) = match_kmer_at_pos(
                index,
                kmer,
                allow_forward_base,
                allow_rev_base,
                diff,
                &mut rev_buf,
                options.kallisto_direct_kmer,
                options.kallisto_enum,
                options.kallisto_strict,
                options.skip_overcrowded_minimizer,
                true,
                true,
            )
        {
            matched = Some((
                uid,
                start,
                pos,
                0usize,
                used_revcomp,
                false,
                false,
                matched_relaxed,
                !used_revcomp,
            ));
        }

        if matched.is_none()
            && in_backoff
            && has_hit
            && env_flag(&BACKOFF_DIRECT_RESCUE, "KALLISTORS_BACKOFF_DIRECT_RESCUE")
            && let Some((uid, start, used_revcomp)) =
                match_kmer_direct(index, kmer, allow_forward_base, allow_rev_base)
            && let Some(block_idx) = block_index_for_position_fast(index, uid, start)
        {
            matched = Some((
                uid,
                start,
                pos,
                0usize,
                used_revcomp,
                false,
                false,
                false,
                !used_revcomp,
            ));
            if let Some(state) = dbg.as_deref_mut()
                && state.dropped_hits.len() < 256
            {
                state.dropped_hits.push(DroppedHitTrace {
                    read_pos: pos,
                    min_pos: 0,
                    unitig_id: uid,
                    unitig_pos: start,
                    block_idx: Some(block_idx),
                    reason: "backoff_direct_rescue",
                    used_revcomp,
                    is_special: false,
                });
            }
        }

        let Some((
            uid,
            start,
            kmer_pos,
            min_pos,
            used_revcomp,
            is_special,
            _minimizer_overcrowded,
            matched_relaxed,
            forward_strand,
        )) = matched
        else {
            pos += 1;
            continue;
        };
        // Keep EC intersection on kallisto-like accepted hits only.
        // Relaxed/fallback matches are treated as probes and should not
        // enter the accepted hit stream used for intersection.
        let mut accept_for_stream = true;
        let mut stream_drop_reason = "";
        if matched_relaxed
            && has_hit
            && !env_flag(&ACCEPT_RELAXED_STREAM, "KALLISTORS_ACCEPT_RELAXED_STREAM")
        {
            accept_for_stream = false;
            stream_drop_reason = if in_backoff {
                "backoff_relaxed_probe"
            } else {
                "relaxed_probe"
            };
        }
        if !accept_for_stream
            && let Some(state) = dbg.as_deref_mut()
            && state.dropped_hits.len() < 256
        {
            state.dropped_hits.push(DroppedHitTrace {
                read_pos: kmer_pos,
                min_pos,
                unitig_id: uid,
                unitig_pos: start,
                block_idx: block_index_for_position_fast(index, uid, start),
                reason: stream_drop_reason,
                used_revcomp,
                is_special,
            });
        }
        let Some(block_idx) = block_index_for_position_fast(index, uid, start) else {
            if let Some(state) = dbg.as_deref_mut()
                && state.dropped_hits.len() < 256
            {
                state.dropped_hits.push(DroppedHitTrace {
                    read_pos: kmer_pos,
                    min_pos,
                    unitig_id: uid,
                    unitig_pos: start,
                    block_idx: None,
                    reason: "block_idx_missing",
                    used_revcomp,
                    is_special,
                });
            }
            pos += 1;
            continue;
        };
        let ec = ec_slice(index, uid, block_idx);
        let ec = if use_shade {
            filter_shades(ec, &index.shade_sequences, &mut shade_scratch);
            shade_scratch.as_slice()
        } else {
            ec
        };
        // Match kallisto's partial intersection behavior: maintain a live running
        // intersection and stop the read as soon as it collapses.
        if accept_for_stream && !options.do_union && !ec.is_empty() {
            if has_online_intersection {
                let mut next = Vec::new();
                intersect_sorted(&online_intersection, ec, &mut next);
                if next.is_empty() {
                    if let Some(state) = dbg.as_deref_mut() {
                        state.intersection_empty = true;
                        if state.dropped_hits.len() < 256 {
                            state.dropped_hits.push(DroppedHitTrace {
                                read_pos: kmer_pos,
                                min_pos,
                                unitig_id: uid,
                                unitig_pos: start,
                                block_idx: Some(block_idx),
                                reason: "intersection_empty",
                                used_revcomp,
                                is_special,
                            });
                        }
                    }
                    return Some(ReadEc {
                        ec: Vec::new(),
                        best_match,
                        first_hit,
                        hits,
                        had_offlist: false,
                        shade_union: Vec::new(),
                        hard_reject_pair: true,
                    });
                }
                online_intersection = next;
            } else {
                online_intersection.clear();
                online_intersection.extend_from_slice(ec);
                has_online_intersection = true;
            }
        }
        if is_special {
            saw_special_hit = true;
        } else {
            saw_non_special_hit = true;
        }
        if let Some(state) = dbg.as_deref_mut() {
            state.saw_match = true;
        }
        if best_match.map(|m| kmer_pos < m.read_pos).unwrap_or(true) {
            best_match = Some(MatchInfo {
                unitig_id: uid,
                unitig_pos: start,
                read_pos: kmer_pos,
                used_revcomp,
            });
        }
        if ec.is_empty() {
            if let Some(state) = dbg.as_deref_mut()
                && state.first_empty_ec.is_none()
            {
                state.first_empty_ec = Some((kmer_pos, min_pos));
            }
        } else if let Some(state) = dbg.as_deref_mut() {
            state.saw_ec = true;
        }
        hits.push(Hit {
            unitig_id: uid,
            read_pos: kmer_pos,
            block_idx,
            used_revcomp,
        });
        if first_hit
            .as_ref()
            .map(|h| kmer_pos < h.read_pos)
            .unwrap_or(true)
        {
            first_hit = Some(Hit {
                unitig_id: uid,
                read_pos: kmer_pos,
                block_idx,
                used_revcomp,
            });
        }
        has_hit = true;
        // Complete the next valid k-mer's fallback probe without starting another
        // jump, matching kallisto's inner probe followed by its outer jump loop.
        if fallback_probe {
            pos += 1;
            continue;
        }
        // Excluding a probe from EC intersection does not invalidate its exact
        // sequence match for confirming a jump, including on reverse unitigs.
        let mut jump_to: Option<usize> = None;
        let mut force_break = false;
        let mut synthetic_hit: Option<Hit> = None;
        let mut next_pos;
        if !options.no_jump
            && let Some(dist) =
                jump_distance_for_match(index, uid, block_idx, start, forward_strand)
        {
            next_pos = pos + dist;
            if next_pos > last_pos {
                next_pos = last_pos;
            }
            if next_pos > pos {
                let Some(next_query_pos) = valid_kmer_start_at_or_after(seq, next_pos, index.k)
                else {
                    break;
                };
                let kmer_next = &seq[next_query_pos..next_query_pos + index.k];
                let next_hit = match_kmer_at_pos(
                    index,
                    kmer_next,
                    allow_forward_base,
                    allow_rev_base,
                    diff,
                    &mut rev_buf,
                    options.kallisto_direct_kmer,
                    options.kallisto_enum,
                    options.kallisto_strict,
                    options.skip_overcrowded_minimizer,
                    options.kallisto_bifrost_find,
                    true,
                );
                let current_ec = ec_slice(index, uid, block_idx);
                let mut jump_reason = "no_jump";
                let mut next_hit_found = false;
                let mut next_hit_relaxed = false;
                let mut next_hit_same_unitig = false;
                let mut next_hit_same_ec = false;
                let mut mid_hit_found = false;
                let mut mid_hit_relaxed = false;
                let mut mid_hit_matches_either = false;
                let mut use_backoff = false;
                if let Some((uid2, _start2, _rev2, block_idx2, next_relaxed)) = next_hit {
                    next_hit_found = true;
                    next_hit_relaxed = next_relaxed;
                    let next_ec_raw = ec_slice(index, uid2, block_idx2);
                    let next_ec = if use_shade {
                        filter_shades(next_ec_raw, &index.shade_sequences, &mut next_shade_scratch);
                        next_shade_scratch.as_slice()
                    } else {
                        next_ec_raw
                    };
                    next_hit_same_unitig = uid2 == uid;
                    next_hit_same_ec = next_ec_raw == current_ec;
                    if next_hit_same_unitig && next_hit_same_ec {
                        if next_pos >= last_pos {
                            synthetic_hit = Some(Hit {
                                unitig_id: uid,
                                read_pos: last_pos,
                                block_idx,
                                used_revcomp,
                            });
                            force_break = true;
                            jump_reason = "next_same_ec_direct_end";
                        } else {
                            synthetic_hit = Some(Hit {
                                unitig_id: uid,
                                read_pos: next_pos,
                                block_idx,
                                used_revcomp,
                            });
                            jump_to = Some(next_query_pos);
                            jump_reason = "next_same_ec_direct";
                        }
                    } else if dist > 4 {
                        let middle_pos = (pos + next_pos) / 2;
                        if let Some(middle_query_pos) =
                            valid_kmer_start_at_or_after(seq, middle_pos, index.k)
                        {
                            let kmer_mid = &seq[middle_query_pos..middle_query_pos + index.k];
                            if let Some((uid3, _start3, rev3, block_idx3, mid_relaxed)) =
                                match_kmer_at_pos(
                                    index,
                                    kmer_mid,
                                    allow_forward_base,
                                    allow_rev_base,
                                    diff,
                                    &mut rev_buf,
                                    options.kallisto_direct_kmer,
                                    options.kallisto_enum,
                                    options.kallisto_strict,
                                    options.skip_overcrowded_minimizer,
                                    options.kallisto_bifrost_find,
                                    true,
                                )
                            {
                                mid_hit_found = true;
                                mid_hit_relaxed = mid_relaxed;
                                let mid_ec = ec_slice(index, uid3, block_idx3);
                                let mid_matches_current = uid3 == uid && mid_ec == current_ec;
                                let mid_matches_next = uid3 == uid2 && mid_ec == next_ec;
                                mid_hit_matches_either = mid_matches_current || mid_matches_next;
                                if mid_matches_current || mid_matches_next {
                                    // Kallisto accepts this middle hit before final jump decision.
                                    if accept_for_stream && !options.do_union && !mid_ec.is_empty()
                                    {
                                        if has_online_intersection {
                                            let mut next = Vec::new();
                                            intersect_sorted(
                                                &online_intersection,
                                                mid_ec,
                                                &mut next,
                                            );
                                            if next.is_empty() {
                                                if let Some(state) = dbg.as_deref_mut() {
                                                    state.intersection_empty = true;
                                                    if state.dropped_hits.len() < 256 {
                                                        state.dropped_hits.push(DroppedHitTrace {
                                                            read_pos: middle_pos,
                                                            min_pos,
                                                            unitig_id: uid3,
                                                            unitig_pos: _start3,
                                                            block_idx: Some(block_idx3),
                                                            reason: "intersection_empty",
                                                            used_revcomp: rev3,
                                                            is_special: false,
                                                        });
                                                    }
                                                }
                                                return Some(ReadEc {
                                                    ec: Vec::new(),
                                                    best_match,
                                                    first_hit,
                                                    hits,
                                                    had_offlist: false,
                                                    shade_union: Vec::new(),
                                                    hard_reject_pair: true,
                                                });
                                            }
                                            online_intersection = next;
                                        } else {
                                            online_intersection.clear();
                                            online_intersection.extend_from_slice(mid_ec);
                                            has_online_intersection = true;
                                        }
                                    }
                                    let found3pos = if mid_matches_current {
                                        middle_pos
                                    } else {
                                        pos + dist
                                    };
                                    synthetic_hit = Some(Hit {
                                        unitig_id: uid3,
                                        read_pos: found3pos,
                                        block_idx: block_idx3,
                                        used_revcomp: rev3,
                                    });
                                    if next_pos >= last_pos {
                                        force_break = true;
                                        jump_reason = "mid_match_direct_end";
                                    } else {
                                        jump_to = Some(next_query_pos);
                                        jump_reason = "mid_match_direct";
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // Kallisto jumps even when the next k-mer is missing.
                    jump_reason = "next_missing";
                    let next_in_dlist = !options.do_union
                        && index.dlist.as_ref().is_some_and(|dlist| {
                            encode_kmer_pair(kmer_next)
                                .is_some_and(|(fwd, rev)| dlist.contains(&fwd.min(rev)))
                        });
                    if next_in_dlist {
                        dlist_early_return = true;
                        jump_reason = "next_missing_dlist_early_return";
                        append_dummy_hit = dlist_onlist.is_some();
                    }
                    if next_pos >= last_pos {
                        synthetic_hit = Some(Hit {
                            unitig_id: uid,
                            read_pos: last_pos,
                            block_idx,
                            used_revcomp,
                        });
                        force_break = true;
                    } else {
                        // In kallisto this corresponds to found2pos = pos.
                        synthetic_hit = Some(Hit {
                            unitig_id: uid,
                            read_pos: kmer_pos,
                            block_idx,
                            used_revcomp,
                        });
                        jump_to = Some(next_query_pos);
                    }
                    if dlist_early_return {
                        force_break = true;
                        jump_to = None;
                    }
                }
                if jump_to.is_none() && !force_break {
                    use_backoff = true;
                }
                if jump_to.is_none() && !force_break && next_hit_found {
                    jump_reason = "next_or_mid_not_direct";
                }
                if let Some(state) = dbg.as_deref_mut()
                    && state.jump_decisions.len() < 512
                {
                    state.jump_decisions.push(JumpDecisionTrace {
                        read_pos: kmer_pos,
                        matched_unitig_id: uid,
                        matched_block_idx: block_idx,
                        jump_distance: dist,
                        next_pos,
                        in_backoff,
                        next_hit_found,
                        next_hit_relaxed,
                        next_hit_same_unitig,
                        next_hit_same_ec,
                        mid_hit_found,
                        mid_hit_relaxed,
                        mid_hit_matches_either,
                        jumped: jump_to.is_some() || force_break,
                        reason: jump_reason,
                    });
                }
                if let Some(hit) = synthetic_hit.take() {
                    hits.push(hit);
                }
                if use_backoff {
                    backoff_until = Some(next_pos);
                    fallback_probe_pending = true;
                }
            }
        }
        if force_break {
            break;
        }
        if let Some(next_pos) = jump_to {
            pos = next_pos.saturating_add(1);
            continue;
        }
        pos += 1;
    }

    // Mirror kallisto's post-loop D-list dummy append in partial mode.
    if !dlist_early_return && any_dlist_kmer && (has_hit || options.do_union) {
        append_dummy_hit = dlist_onlist.is_some();
    }

    finish::finish_read(
        index,
        seq,
        dbg,
        options,
        finish::ScannedRead {
            shade_scratch,
            best_match,
            first_hit,
            hits,
            has_hit,
            saw_special_hit,
            saw_non_special_hit,
            append_dummy_hit,
        },
    )
}
