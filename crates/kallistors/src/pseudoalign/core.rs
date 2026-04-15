use super::debug::{
    DroppedHitTrace, JumpDecisionTrace, MinimizerCandidateTrace, ReadDebugState,
    format_positions_sample,
};
use super::ec::encode_kmer_pair;
use super::minimizers::{
    fill_revcomp, minhash_next_after_hash, minhash_primary_for_kmer, minimizer_for_kmer_strict,
    minimizer_tail_for_kmer, minimizers_for_kmer, minimizers_ranked_for_kmer,
};
use super::{BifrostIndex, KmerEcIndex, PseudoalignOptions};
use super::{DebugFailReason, DebugReport, Hit, MatchInfo, ReadEc, Strand, StrandSpecific};
use super::{
    ec_has_offlist, filter_shades, intersect_sorted, merge_sorted_unique, merge_sorted_unique_vec,
};
use std::collections::{HashMap, HashSet};

fn minimizer_hex(bytes: [u8; 8]) -> String {
    let mut out = String::with_capacity(16);
    for b in bytes {
        out.push_str(&format!("{:02x}", b));
    }
    out
}

#[inline]
fn onlist_cardinality(onlist: &[bool]) -> u32 {
    onlist.iter().filter(|&&v| v).count() as u32
}

#[allow(clippy::too_many_arguments)]
fn match_unitig_candidate(
    unitig: &[u8],
    k: usize,
    rel_pos: usize,
    diff: usize,
    min_pos_fwd: usize,
    min_pos_rev: usize,
    kmer: &[u8],
    rev_kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
    allow_relaxed: bool,
) -> Option<(usize, bool, bool, bool)> {
    if unitig.len() < k {
        return None;
    }

    let rel = rel_pos as isize;
    let start_fwd = rel - min_pos_fwd as isize;
    if allow_forward && start_fwd >= 0 {
        let start = start_fwd as usize;
        if start + k <= unitig.len() && &unitig[start..start + k] == kmer {
            return Some((start, false, false, true));
        }
    }

    let start_rev = rel - diff as isize + min_pos_rev as isize;
    if allow_rev && start_rev >= 0 {
        let start = start_rev as usize;
        if start + k <= unitig.len() && &unitig[start..start + k] == rev_kmer {
            return Some((start, true, false, false));
        }
    }

    if !allow_relaxed {
        return None;
    }

    // Relaxed fallback: allow any k-mer start whose window can contain
    // the selected minimizer occurrence (rel_pos - (k-g) .. rel_pos).
    let max_start = unitig.len() - k;
    let start_lo = rel_pos.saturating_sub(diff).min(max_start);
    let start_hi = rel_pos.min(max_start);
    if start_lo > start_hi {
        return None;
    }
    for start in start_lo..=start_hi {
        if allow_forward && &unitig[start..start + k] == kmer {
            return Some((start, false, true, true));
        }
        if allow_rev && &unitig[start..start + k] == rev_kmer {
            return Some((start, true, true, false));
        }
    }

    None
}

pub(super) fn ec_for_read_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    if seq.len() < index.k {
        if let Some(state) = dbg {
            state.saw_valid_kmer = false;
        }
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
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
    let mut online_intersection: Option<Vec<u32>> = None;
    let mut saw_special_hit = false;
    let mut saw_non_special_hit = false;
    let dlist_dummy = index.onlist.as_deref().map(onlist_cardinality);
    let mut any_dlist_kmer = false;
    let mut dlist_early_return = false;
    let mut extra_ec_hits: Vec<Vec<u32>> = Vec::new();
    // Kallisto falls back to incremental scanning after a failed jump; track that window.
    let mut backoff_until: Option<usize> = None;
    // Legacy parity switch kept for CLI compatibility.
    let _kallisto_sparse_hits = options.kallisto_sparse_hits;

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
        let mut matched = None;
        let mut any_mphf = false;
        let mut any_positions = false;
        if !allow_forward_base && !allow_rev_base {
            pos += 1;
            continue;
        }
        if !fill_revcomp(kmer, &mut rev_buf) {
            pos += 1;
            continue;
        }
        let use_revcomp = rev_buf.as_slice() < kmer;
        let kmer_canon = if use_revcomp {
            rev_buf.as_slice()
        } else {
            kmer
        };
        let min_input = if options.kallisto_bifrost_find || options.kallisto_strict {
            kmer_canon
        } else {
            kmer
        };

        if options.kallisto_direct_kmer {
            if let Some((uid, start, used_revcomp)) =
                match_kmer_direct(index, kmer, allow_forward_base, allow_rev_base)
            {
                matched = Some((
                    uid,
                    start,
                    pos,
                    0,
                    used_revcomp,
                    false,
                    false,
                    false,
                    !used_revcomp,
                ));
            }
            if matched.is_none() {
                pos += 1;
                continue;
            }
        }

        let (mut min_candidates, mut min_hash_current) = if options.kallisto_bifrost_find {
            match minhash_primary_for_kmer(min_input, index.g) {
                Some((min_hash, candidates)) => (candidates, Some(min_hash)),
                None => {
                    pos += 1;
                    continue;
                }
            }
        } else {
            let candidates = if options.kallisto_strict {
                minimizer_for_kmer_strict(min_input, index.g).map(|v| vec![v])
            } else if options.kallisto_enum {
                minimizers_ranked_for_kmer(min_input, index.g, 2)
            } else {
                minimizers_for_kmer(min_input, index.g)
            };
            let Some(candidates) = candidates else {
                pos += 1;
                continue;
            };
            (candidates, None)
        };
        if let Some(state) = dbg.as_deref_mut() {
            state.saw_valid_kmer = true;
        }
        let mut tried_tail_minimizer = false;
        let mut saw_special_overcrowded_unresolved = false;
        'min_outer: loop {
            let mut request_next_min = false;
            for (min_bytes, min_pos) in min_candidates.iter().copied() {
                if matched.is_some() {
                    break;
                }
                let mut candidate = MinimizerCandidateTrace {
                    read_pos: pos,
                    min_pos,
                    minimizer: minimizer_hex(min_bytes),
                    mphf_hit: false,
                    positions_len: 0,
                    sample_positions: "-".to_string(),
                    has_special: false,
                    overcrowded: false,
                    matched: false,
                };
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
                let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
                    if let Some(state) = dbg.as_deref_mut()
                        && state.first_mphf_miss.is_none()
                    {
                        state.first_mphf_miss = Some((pos, min_pos));
                    }
                    if let Some(state) = dbg.as_deref_mut()
                        && state.minimizer_candidates.len() < 1024
                    {
                        state.minimizer_candidates.push(candidate);
                    }
                    continue;
                };
                candidate.mphf_hit = true;
                any_mphf = true;
                let positions = &index.minz_positions[min_idx as usize];
                if positions.is_empty() {
                    if let Some(state) = dbg.as_deref_mut()
                        && state.first_no_positions.is_none()
                    {
                        state.first_no_positions = Some((pos, min_pos));
                    }
                    if let Some(state) = dbg.as_deref_mut()
                        && state.minimizer_candidates.len() < 1024
                    {
                        state.minimizer_candidates.push(candidate);
                    }
                    continue;
                }
                let mut minimizer_overcrowded = false;
                if options.skip_overcrowded_minimizer {
                    for &pos_id in positions {
                        let unitig_id_raw = (pos_id >> 32) as u32;
                        if unitig_id_raw == u32::MAX && (pos_id & 0x8000_0000) != 0 {
                            minimizer_overcrowded = true;
                            break;
                        }
                    }
                }
                candidate.positions_len = positions.len();
                candidate.sample_positions = format_positions_sample(positions);
                candidate.overcrowded = minimizer_overcrowded;
                any_positions = true;
                if let Some(state) = dbg.as_deref_mut()
                    && state.first_no_match_positions.is_none()
                {
                    state.first_no_match_positions = Some((
                        pos,
                        min_pos,
                        positions.len(),
                        format_positions_sample(positions),
                    ));
                }
                let mut local_match = None;
                for &pos_id in positions {
                    let unitig_id_raw = (pos_id >> 32) as u32;
                    if unitig_id_raw == u32::MAX {
                        candidate.has_special = true;
                        if options.kallisto_strict {
                            continue;
                        }
                        let overcrowded = (pos_id & 0x8000_0000) != 0;
                        candidate.overcrowded = overcrowded;

                        let special_uid = if options.kallisto_bifrost_find {
                            special_unitig_for_kmer_raw(index, kmer)
                        } else {
                            special_unitig_for_kmer(index, kmer)
                        };
                        if let Some(uid) = special_uid {
                            let used_revcomp = kmer_canon != kmer;
                            if (used_revcomp && !allow_rev_base)
                                || (!used_revcomp && !allow_forward_base)
                            {
                                continue;
                            }
                            candidate.matched = true;
                            local_match = Some((
                                uid,
                                0usize,
                                pos,
                                min_pos,
                                used_revcomp,
                                true,
                                minimizer_overcrowded,
                                false,
                                !used_revcomp,
                            ));
                            break;
                        }
                        if overcrowded {
                            saw_special_overcrowded_unresolved = true;
                        }
                        if options.kallisto_bifrost_find && overcrowded {
                            request_next_min = true;
                        } else if options.skip_overcrowded_minimizer && overcrowded {
                            if let Some(state) = dbg.as_deref_mut() {
                                state.overcrowded_jump_skip = true;
                            }
                            continue;
                        }
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
                        if km_seq != kmer_canon {
                            continue;
                        }
                        if min_pos != km_pos && min_pos + km_pos != diff {
                            continue;
                        }
                        let used_revcomp = kmer_canon != kmer;
                        if (used_revcomp && !allow_rev_base)
                            || (!used_revcomp && !allow_forward_base)
                        {
                            continue;
                        }
                        candidate.matched = true;
                        local_match = Some((
                            index.unitigs.len() + uid,
                            0usize,
                            pos,
                            min_pos,
                            used_revcomp,
                            false,
                            minimizer_overcrowded,
                            false,
                            !used_revcomp,
                        ));
                        if used_revcomp && let Some(state) = dbg.as_deref_mut() {
                            state.used_revcomp = true;
                        }
                        break;
                    } else {
                        let uid = unitig_id as usize;
                        if uid >= index.unitigs.len() {
                            continue;
                        }
                        if let Some((start, used_revcomp, matched_relaxed, forward_strand)) =
                            match_unitig_candidate(
                                index.unitigs[uid].as_slice(),
                                index.k,
                                rel_pos as usize,
                                diff,
                                min_pos_fwd,
                                min_pos_rev,
                                kmer,
                                rev_buf.as_slice(),
                                allow_forward_base,
                                allow_rev_base,
                                true,
                            )
                        {
                            candidate.matched = true;
                            local_match = Some((
                                uid,
                                start,
                                pos,
                                min_pos,
                                used_revcomp,
                                false,
                                minimizer_overcrowded,
                                matched_relaxed,
                                forward_strand,
                            ));
                            if used_revcomp && let Some(state) = dbg.as_deref_mut() {
                                state.used_revcomp = true;
                            }
                            break;
                        }
                    }
                }
                if let Some(hit) = local_match {
                    matched = Some(hit);
                    if let Some(state) = dbg.as_deref_mut()
                        && state.minimizer_candidates.len() < 1024
                    {
                        state.minimizer_candidates.push(candidate);
                    }
                    break;
                }
                if let Some(state) = dbg.as_deref_mut()
                    && state.first_no_match.is_none()
                {
                    state.first_no_match = Some((pos, min_pos));
                }
                if let Some(state) = dbg.as_deref_mut()
                    && state.minimizer_candidates.len() < 1024
                {
                    state.minimizer_candidates.push(candidate);
                }
            }
            if matched.is_some() {
                break;
            }
            if !options.kallisto_bifrost_find
                && !options.kallisto_strict
                && saw_special_overcrowded_unresolved
            {
                let mut probe_rev_buf: Vec<u8> = Vec::new();
                if let Some((uid, start, used_revcomp, _block_idx, matched_relaxed)) =
                    match_kmer_at_pos(
                        index,
                        kmer,
                        allow_forward_base,
                        allow_rev_base,
                        diff,
                        &mut probe_rev_buf,
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
                    break;
                }
            }
            if !options.kallisto_bifrost_find
                && !options.kallisto_strict
                && !options.kallisto_enum
                && !in_backoff
                && !tried_tail_minimizer
                && let Some((tail_bytes, tail_pos)) = minimizer_tail_for_kmer(min_input, index.g)
                && !min_candidates
                    .iter()
                    .any(|&(min_bytes, min_pos)| min_bytes == tail_bytes && min_pos == tail_pos)
            {
                tried_tail_minimizer = true;
                min_candidates = vec![(tail_bytes, tail_pos)];
                continue 'min_outer;
            }
            if !options.kallisto_bifrost_find {
                break;
            }
            if request_next_min
                && let Some(curr_hash) = min_hash_current
                && let Some((next_hash, next_min)) =
                    minhash_next_after_hash(min_input, index.g, curr_hash)
            {
                min_hash_current = Some(next_hash);
                min_candidates = vec![next_min];
                continue 'min_outer;
            }
            break;
        }
        if let Some(state) = dbg.as_deref_mut() {
            if any_mphf {
                state.saw_mphf_hit = true;
            }
            if any_positions {
                state.saw_positions = true;
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
        if matched_relaxed && has_hit {
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
                block_idx: block_index_for_position(&index.ec_blocks[uid], start),
                reason: stream_drop_reason,
                used_revcomp,
                is_special,
            });
        }
        let Some(block_idx) = block_index_for_position(&index.ec_blocks[uid], start) else {
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
        let ec = &index.ec_blocks[uid][block_idx].ec;
        let ec = if use_shade {
            filter_shades(ec, &index.shade_sequences, &mut shade_scratch);
            shade_scratch.as_slice()
        } else {
            ec
        };
        // Match kallisto's partial intersection behavior: maintain a live running
        // intersection and stop the read as soon as it collapses.
        if accept_for_stream && !options.do_union && !ec.is_empty() {
            if let Some(current) = online_intersection.as_ref() {
                let mut next = Vec::new();
                intersect_sorted(current, ec, &mut next);
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
                online_intersection = Some(next);
            } else {
                online_intersection = Some(ec.to_vec());
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
                let kmer_next = &seq[next_pos..next_pos + index.k];
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
                let current_ec = &index.ec_blocks[uid][block_idx].ec;
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
                    let next_ec = &index.ec_blocks[uid2][block_idx2].ec;
                    next_hit_same_unitig = uid2 == uid;
                    next_hit_same_ec = next_ec == current_ec;
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
                            jump_to = Some(next_pos);
                            jump_reason = "next_same_ec_direct";
                        }
                    } else if dist > 4 {
                        let middle_pos = (pos + next_pos) / 2;
                        if middle_pos <= last_pos {
                            let kmer_mid = &seq[middle_pos..middle_pos + index.k];
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
                                let mid_ec = &index.ec_blocks[uid3][block_idx3].ec;
                                let mid_matches_current = uid3 == uid && mid_ec == current_ec;
                                let mid_matches_next = uid3 == uid2 && mid_ec == next_ec;
                                mid_hit_matches_either = mid_matches_current || mid_matches_next;
                                if mid_matches_current || mid_matches_next {
                                    // Kallisto accepts this middle hit before final jump decision.
                                    if accept_for_stream && !options.do_union && !mid_ec.is_empty()
                                    {
                                        if let Some(current) = online_intersection.as_ref() {
                                            let mut next = Vec::new();
                                            intersect_sorted(current, mid_ec, &mut next);
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
                                            online_intersection = Some(next);
                                        } else {
                                            online_intersection = Some(mid_ec.to_vec());
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
                                        jump_to = Some(next_pos);
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
                        if let Some(dummy) = dlist_dummy {
                            extra_ec_hits.push(vec![dummy]);
                        }
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
                        jump_to = Some(next_pos);
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
    if !dlist_early_return
        && any_dlist_kmer
        && (has_hit || options.do_union)
        && let Some(dummy) = dlist_dummy
    {
        extra_ec_hits.push(vec![dummy]);
    }

    if options.discard_special_only && saw_special_hit && !saw_non_special_hit {
        if let Some(state) = dbg.as_deref_mut() {
            state.special_only = true;
        }
        return None;
    }

    if !has_hit || hits.is_empty() {
        if (options.kallisto_local_fallback || options.kallisto_bifrost_find)
            && !options.kallisto_strict
        {
            let unitigs = collect_unitigs_for_read(
                index,
                seq,
                options.kallisto_enum,
                options.kallisto_strict,
                options.skip_overcrowded_minimizer,
                options.kallisto_bifrost_find,
                256,
            );
            if !unitigs.is_empty()
                && let Some(ec) =
                    ec_for_read_local_kmer(index, seq, &unitigs, true, options.discard_special_only)
            {
                if let Some(state) = dbg.as_deref_mut() {
                    state.saw_valid_kmer = true;
                    state.saw_ec = true;
                }
                return Some(ReadEc {
                    ec,
                    best_match: None,
                    first_hit: None,
                    hits: Vec::new(),
                    had_offlist: false,
                    shade_union: Vec::new(),
                    hard_reject_pair: false,
                });
            }
        }
        if options.kallisto_fallback
            && let Some(kmer_index) = index.kmer_index.as_ref()
            && let Some(ec) = ec_for_read_kmer_index(kmer_index, seq)
        {
            if let Some(state) = dbg.as_deref_mut() {
                state.saw_valid_kmer = true;
                state.saw_ec = true;
            }
            return Some(ReadEc {
                ec,
                best_match: None,
                first_hit: None,
                hits: Vec::new(),
                had_offlist: false,
                shade_union: Vec::new(),
                hard_reject_pair: false,
            });
        }
        return None;
    }

    let mut minpos = usize::MAX;
    let mut maxpos = 0usize;
    let mut last_ec: Vec<u32> = Vec::new();
    let mut last_unitig = None;
    let mut found_nonempty = false;
    let mut current_offlist = false;
    let mut ec_scratch: Vec<u32> = Vec::new();
    for hit in &hits {
        minpos = minpos.min(hit.read_pos);
        maxpos = maxpos.max(hit.read_pos);
        let ec = &index.ec_blocks[hit.unitig_id][hit.block_idx].ec;
        let ec = if use_shade {
            filter_shades(ec, &index.shade_sequences, &mut shade_scratch);
            shade_scratch.as_slice()
        } else {
            ec
        };
        let ec_offlist = if let Some(onlist) = index.onlist.as_deref() {
            ec_has_offlist(ec, onlist)
        } else {
            false
        };
        if !found_nonempty {
            if !ec.is_empty() {
                if options.do_union {
                    current.extend_from_slice(ec);
                    current.sort_unstable();
                    current.dedup();
                } else {
                    current.extend_from_slice(ec);
                    last_ec = ec.to_vec();
                    last_unitig = Some(hit.unitig_id);
                }
                current_offlist = ec_offlist;
                found_nonempty = true;
            }
            continue;
        }
        if last_unitig == Some(hit.unitig_id) && ec == last_ec.as_slice() {
            continue;
        }
        if ec.is_empty() {
            continue;
        }
        if options.do_union {
            if options.dfk_onlist
                && (current_offlist || ec_offlist)
                && let Some(dummy) = index.onlist.as_ref().map(|v| onlist_cardinality(v))
                && current.last().copied() != Some(dummy)
            {
                current.push(dummy);
            }
            merge_sorted_unique_vec(&mut current, ec);
            current_offlist = current_offlist || ec_offlist;
        } else {
            if current.is_empty() {
                current.extend_from_slice(ec);
                last_ec = ec.to_vec();
                last_unitig = Some(hit.unitig_id);
                current_offlist = ec_offlist;
                continue;
            }
            let ec_slice = if options.dfk_onlist && (current_offlist || ec_offlist) {
                if let Some(onlist) = index.onlist.as_deref() {
                    let dummy = onlist_cardinality(onlist);
                    ec_scratch.clear();
                    ec_scratch.extend_from_slice(ec);
                    if ec_scratch.last().copied() != Some(dummy) {
                        ec_scratch.push(dummy);
                    }
                    if current.last().copied() != Some(dummy) {
                        current.push(dummy);
                    }
                    ec_scratch.sort_unstable();
                    ec_scratch.dedup();
                    ec_scratch.as_slice()
                } else {
                    ec
                }
            } else {
                ec
            };
            next.clear();
            intersect_sorted(&current, ec_slice, &mut next);
            std::mem::swap(&mut current, &mut next);
            if current.is_empty() {
                if let Some(state) = dbg.as_deref_mut() {
                    state.intersection_empty = true;
                }
                return Some(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: current_offlist,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                });
            }
            last_ec = ec.to_vec();
            last_unitig = Some(hit.unitig_id);
            current_offlist = current_offlist || ec_offlist;
        }
    }

    // Apply synthetic D-list dummy hits (kallisto um_dummy behavior).
    for ec in &extra_ec_hits {
        if ec.is_empty() {
            continue;
        }
        if options.do_union {
            merge_sorted_unique_vec(&mut current, ec);
            current_offlist = true;
            continue;
        }
        if current.is_empty() {
            current.extend_from_slice(ec);
            current_offlist = true;
            continue;
        }
        next.clear();
        intersect_sorted(&current, ec, &mut next);
        std::mem::swap(&mut current, &mut next);
        if current.is_empty() {
            if let Some(state) = dbg.as_deref_mut() {
                state.intersection_empty = true;
            }
            return Some(ReadEc {
                ec: Vec::new(),
                best_match,
                first_hit,
                hits,
                had_offlist: current_offlist,
                shade_union: Vec::new(),
                hard_reject_pair: true,
            });
        }
        current_offlist = true;
    }

    if current.is_empty() {
        return None;
    }
    let min_range = options.min_range.max(1);
    if maxpos >= minpos && (maxpos - minpos + index.k) < min_range {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
        if options.dfk_onlist && current_offlist && !current.is_empty() {
            let dummy = onlist_cardinality(onlist);
            if current.last().copied() != Some(dummy) {
                current.push(dummy);
            }
        }
    }
    if current.is_empty() {
        None
    } else {
        let mut shade_union = Vec::new();
        if index.use_shade {
            for hit in &hits {
                let ec = &index.ec_blocks[hit.unitig_id][hit.block_idx].ec;
                for &tr in ec {
                    let idx = tr as usize;
                    if idx < index.shade_sequences.len() && index.shade_sequences[idx] {
                        shade_union.push(tr);
                    }
                }
            }
            shade_union.sort_unstable();
            shade_union.dedup();
        }
        Some(ReadEc {
            ec: current,
            best_match,
            first_hit,
            hits,
            had_offlist: current_offlist,
            shade_union,
            hard_reject_pair: false,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn match_kmer_at_pos(
    index: &BifrostIndex,
    kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
    _diff: usize,
    rev_buf: &mut Vec<u8>,
    kallisto_direct_kmer: bool,
    kallisto_enum: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_bifrost_find: bool,
    allow_relaxed: bool,
) -> Option<(usize, usize, bool, usize, bool)> {
    if !fill_revcomp(kmer, rev_buf) {
        return None;
    }
    let use_revcomp = rev_buf.as_slice() < kmer;
    let kmer_canon = if use_revcomp {
        rev_buf.as_slice()
    } else {
        kmer
    };
    let min_input = if kallisto_bifrost_find || kallisto_strict {
        kmer_canon
    } else {
        kmer
    };
    if kallisto_direct_kmer
        && let Some((uid, start, used_revcomp)) =
            match_kmer_direct(index, kmer, allow_forward, allow_rev)
    {
        let block_idx = block_index_for_position(&index.ec_blocks[uid], start)?;
        return Some((uid, start, used_revcomp, block_idx, false));
    }
    let (mut min_candidates, mut min_hash_current) = if kallisto_bifrost_find {
        let (min_hash, candidates) = minhash_primary_for_kmer(min_input, index.g)?;
        (candidates, Some(min_hash))
    } else {
        let candidates = if kallisto_strict {
            vec![minimizer_for_kmer_strict(min_input, index.g)?]
        } else if kallisto_enum {
            minimizers_ranked_for_kmer(min_input, index.g, 2)?
        } else {
            minimizers_for_kmer(min_input, index.g)?
        };
        (candidates, None)
    };
    if !kallisto_bifrost_find
        && !kallisto_strict
        && min_candidates.iter().any(|(min_bytes, _)| {
            index
                .mphf
                .lookup(min_bytes)
                .map(|min_idx| {
                    index.minz_positions[min_idx as usize]
                        .iter()
                        .any(|&pos_id| {
                            (pos_id >> 32) as u32 == u32::MAX && (pos_id & 0x8000_0000) != 0
                        })
                })
                .unwrap_or(false)
        })
    {
        return match_kmer_at_pos(
            index,
            kmer,
            allow_forward,
            allow_rev,
            _diff,
            rev_buf,
            kallisto_direct_kmer,
            kallisto_enum,
            kallisto_strict,
            skip_overcrowded_minimizer,
            true,
            allow_relaxed,
        );
    }
    let mut tried_tail_minimizer = false;
    'min_outer: loop {
        let mut request_next_min = false;
        for (min_bytes, min_pos) in min_candidates.iter().copied() {
            let min_pos_fwd = if use_revcomp {
                _diff.saturating_sub(min_pos)
            } else {
                min_pos
            };
            let min_pos_rev = if use_revcomp {
                min_pos
            } else {
                _diff.saturating_sub(min_pos)
            };
            let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
                continue;
            };
            let positions = &index.minz_positions[min_idx as usize];
            if positions.is_empty() {
                continue;
            }
            for &pos_id in positions {
                let unitig_id_raw = (pos_id >> 32) as u32;
                if unitig_id_raw == u32::MAX {
                    if kallisto_strict {
                        continue;
                    }
                    let overcrowded = (pos_id & 0x8000_0000) != 0;
                    if !kallisto_bifrost_find && overcrowded && skip_overcrowded_minimizer {
                        continue;
                    }
                    let special_uid = if kallisto_bifrost_find {
                        special_unitig_for_kmer_raw(index, kmer)
                    } else {
                        special_unitig_for_kmer(index, kmer)
                    };
                    if let Some(uid) = special_uid {
                        let used_revcomp = kmer_canon != kmer;
                        if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                            continue;
                        }
                        let block_idx = block_index_for_position(&index.ec_blocks[uid], 0)?;
                        return Some((uid, 0usize, used_revcomp, block_idx, false));
                    }
                    if kallisto_bifrost_find && overcrowded {
                        request_next_min = true;
                    }
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
                    if km_seq != kmer_canon {
                        continue;
                    }
                    if min_pos != km_pos && min_pos + km_pos != _diff {
                        continue;
                    }
                    let used_revcomp = kmer_canon != kmer;
                    if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                        continue;
                    }
                    let unitig_id = index.unitigs.len() + uid;
                    let block_idx = block_index_for_position(&index.ec_blocks[unitig_id], 0)?;
                    return Some((unitig_id, 0usize, used_revcomp, block_idx, false));
                } else {
                    let uid = unitig_id as usize;
                    if uid >= index.unitigs.len() {
                        continue;
                    }
                    if let Some((start, used_revcomp, matched_relaxed, _forward_strand)) =
                        match_unitig_candidate(
                            index.unitigs[uid].as_slice(),
                            index.k,
                            rel_pos as usize,
                            _diff,
                            min_pos_fwd,
                            min_pos_rev,
                            kmer,
                            rev_buf.as_slice(),
                            allow_forward,
                            allow_rev,
                            allow_relaxed,
                        )
                    {
                        let block_idx = block_index_for_position(&index.ec_blocks[uid], start)?;
                        return Some((uid, start, used_revcomp, block_idx, matched_relaxed));
                    }
                }
            }
        }
        if !kallisto_bifrost_find
            && !kallisto_strict
            && !kallisto_enum
            && !tried_tail_minimizer
            && let Some((tail_bytes, tail_pos)) = minimizer_tail_for_kmer(min_input, index.g)
            && !min_candidates
                .iter()
                .any(|&(min_bytes, min_pos)| min_bytes == tail_bytes && min_pos == tail_pos)
        {
            tried_tail_minimizer = true;
            min_candidates = vec![(tail_bytes, tail_pos)];
            continue 'min_outer;
        }
        if !kallisto_bifrost_find {
            break;
        }
        if request_next_min
            && let Some(curr_hash) = min_hash_current
            && let Some((next_hash, next_min)) =
                minhash_next_after_hash(min_input, index.g, curr_hash)
        {
            min_hash_current = Some(next_hash);
            min_candidates = vec![next_min];
            continue 'min_outer;
        }
        break;
    }
    None
}

fn match_kmer_direct(
    index: &BifrostIndex,
    kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
) -> Option<(usize, usize, bool)> {
    let map = index.kmer_pos_index.as_ref()?;
    let (fwd, rev) = encode_kmer_pair(kmer)?;
    if let Some(entries) = map.get(&fwd) {
        for &(unitig_id, pos, used_revcomp) in entries {
            if used_revcomp {
                continue;
            }
            if allow_forward {
                return Some((unitig_id, pos, false));
            }
        }
    }
    if let Some(entries) = map.get(&rev) {
        for &(unitig_id, pos, used_revcomp) in entries {
            if !used_revcomp {
                continue;
            }
            if allow_rev {
                return Some((unitig_id, pos, true));
            }
        }
    }
    None
}

fn collect_unitigs_for_read(
    index: &BifrostIndex,
    seq: &[u8],
    kallisto_enum: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_bifrost_find: bool,
    max_unitigs: usize,
) -> Vec<usize> {
    if seq.len() < index.k {
        return Vec::new();
    }
    let mut out: HashSet<usize> = HashSet::new();
    let mut rev_buf: Vec<u8> = Vec::new();
    for pos in 0..=seq.len() - index.k {
        let kmer = &seq[pos..pos + index.k];
        if !fill_revcomp(kmer, &mut rev_buf) {
            continue;
        }
        let use_revcomp = rev_buf.as_slice() < kmer;
        let kmer_canon = if use_revcomp {
            rev_buf.as_slice()
        } else {
            kmer
        };
        let min_input = if kallisto_bifrost_find || kallisto_strict {
            kmer_canon
        } else {
            kmer
        };
        let (mut min_candidates, mut min_hash_current) = if kallisto_bifrost_find {
            let (min_hash, candidates) = match minhash_primary_for_kmer(min_input, index.g) {
                Some(v) => v,
                None => continue,
            };
            (candidates, Some(min_hash))
        } else {
            let candidates = if kallisto_strict {
                match minimizer_for_kmer_strict(min_input, index.g) {
                    Some(v) => vec![v],
                    None => continue,
                }
            } else if kallisto_enum {
                match minimizers_ranked_for_kmer(min_input, index.g, 2) {
                    Some(v) => v,
                    None => continue,
                }
            } else {
                match minimizers_for_kmer(min_input, index.g) {
                    Some(v) => v,
                    None => continue,
                }
            };
            (candidates, None)
        };
        'min_loop: loop {
            let mut request_next_min = false;
            for (min_bytes, _min_pos) in min_candidates.iter() {
                let Some(min_idx) = index.mphf.lookup(min_bytes) else {
                    continue;
                };
                let positions = &index.minz_positions[min_idx as usize];
                if positions.is_empty() {
                    continue;
                }
                for &pos_id in positions {
                    let unitig_id_raw = (pos_id >> 32) as u32;
                    if unitig_id_raw == u32::MAX {
                        if !kallisto_strict {
                            let overcrowded = (pos_id & 0x8000_0000) != 0;
                            if !kallisto_bifrost_find && overcrowded && skip_overcrowded_minimizer {
                                continue;
                            }
                            let special_uid = if kallisto_bifrost_find {
                                special_unitig_for_kmer_raw(index, kmer)
                            } else {
                                special_unitig_for_kmer(index, kmer)
                            };
                            if let Some(uid) = special_uid {
                                out.insert(uid);
                                if out.len() >= max_unitigs {
                                    return out.into_iter().collect();
                                }
                            }
                            if kallisto_bifrost_find && overcrowded {
                                request_next_min = true;
                            }
                        }
                        continue;
                    }
                    let (unitig_id, _rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
                    let uid = unitig_id as usize;
                    let real_id = if is_km {
                        index.unitigs.len() + uid
                    } else {
                        uid
                    };
                    out.insert(real_id);
                    if out.len() >= max_unitigs {
                        return out.into_iter().collect();
                    }
                }
            }
            if !kallisto_bifrost_find {
                break;
            }
            if request_next_min
                && let Some(curr_hash) = min_hash_current
                && let Some((next_hash, next_min)) =
                    minhash_next_after_hash(min_input, index.g, curr_hash)
            {
                min_hash_current = Some(next_hash);
                min_candidates = vec![next_min];
                continue 'min_loop;
            }
            break;
        }
    }
    out.into_iter().collect()
}

fn build_local_kmer_map(index: &BifrostIndex, unitig_ids: &[usize]) -> HashMap<u64, Vec<u32>> {
    let mut map: HashMap<u64, Vec<u32>> = HashMap::new();
    for &unitig_id in unitig_ids {
        let seq_ref: Option<&[u8]> = if unitig_id < index.unitigs.len() {
            Some(index.unitigs[unitig_id].as_slice())
        } else {
            let km_idx = unitig_id - index.unitigs.len();
            index.km_unitigs.get(km_idx).map(|v| v.as_slice())
        };
        let Some(seq_ref) = seq_ref else { continue };
        let blocks = match index.ec_blocks.get(unitig_id) {
            Some(v) => v,
            None => continue,
        };
        for block in blocks {
            let start = block.lb as usize;
            let end = block.ub as usize;
            for pos in start..end {
                if pos + index.k > seq_ref.len() {
                    break;
                }
                if let Some((fwd, rev)) = encode_kmer_pair(&seq_ref[pos..pos + index.k]) {
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
    map
}

pub fn local_kmer_hits(
    index: &BifrostIndex,
    seq: &[u8],
    kallisto_enum: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_bifrost_find: bool,
    max_unitigs: usize,
) -> Vec<(usize, String, bool, usize)> {
    if seq.len() < index.k {
        return Vec::new();
    }
    let unitigs = collect_unitigs_for_read(
        index,
        seq,
        kallisto_enum,
        kallisto_strict,
        skip_overcrowded_minimizer,
        kallisto_bifrost_find,
        max_unitigs,
    );
    let map = build_local_kmer_map(index, &unitigs);
    let mut out = Vec::new();
    for pos in 0..=seq.len() - index.k {
        let kmer = &seq[pos..pos + index.k];
        let mut hit = false;
        let mut ec_size = 0usize;
        if let Some((fwd, rev)) = encode_kmer_pair(kmer)
            && let Some(ec) = map.get(&fwd).or_else(|| map.get(&rev))
        {
            hit = true;
            ec_size = ec.len();
        } else if !kallisto_strict && let Some(ec) = special_ec_for_kmer(index, kmer) {
            hit = true;
            ec_size = ec.len();
        }
        out.push((
            pos,
            String::from_utf8_lossy(kmer).into_owned(),
            hit,
            ec_size,
        ));
    }
    out
}

pub(super) fn special_unitig_for_kmer(index: &BifrostIndex, kmer: &[u8]) -> Option<usize> {
    let map = index.h_kmer_map.as_ref()?;
    let (fwd, rev) = encode_kmer_pair(kmer)?;
    let code = fwd.min(rev);
    map.get(&code).copied()
}

pub(super) fn special_unitig_for_kmer_raw(index: &BifrostIndex, kmer: &[u8]) -> Option<usize> {
    let map = index.h_kmer_map.as_ref()?;
    let (fwd, rev) = encode_kmer_pair(kmer)?;
    let code = fwd.min(rev);
    map.get(&code).copied()
}

fn special_ec_for_kmer<'a>(index: &'a BifrostIndex, kmer: &'a [u8]) -> Option<&'a [u32]> {
    let uid = special_unitig_for_kmer(index, kmer)?;
    let blocks = index.ec_blocks.get(uid)?;
    let block_idx = block_index_for_position(blocks, 0)?;
    Some(&blocks[block_idx].ec)
}

fn ec_for_read_local_kmer(
    index: &BifrostIndex,
    seq: &[u8],
    unitig_ids: &[usize],
    allow_special: bool,
    discard_special_only: bool,
) -> Option<Vec<u32>> {
    if seq.len() < index.k {
        return None;
    }
    let map = build_local_kmer_map(index, unitig_ids);
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut has_hit = false;
    let mut saw_non_special = false;
    let mut saw_special = false;
    for pos in 0..=seq.len() - index.k {
        if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + index.k])
            && let Some(ec) = map.get(&fwd).or_else(|| map.get(&rev))
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                intersect_sorted(&current, ec, &mut next);
                std::mem::swap(&mut current, &mut next);
                if current.is_empty() {
                    break;
                }
            }
            saw_non_special = true;
        } else if allow_special
            && let Some(ec) = special_ec_for_kmer(index, &seq[pos..pos + index.k])
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                intersect_sorted(&current, ec, &mut next);
                std::mem::swap(&mut current, &mut next);
                if current.is_empty() {
                    break;
                }
            }
            saw_special = true;
        }
    }
    if discard_special_only && saw_special && !saw_non_special {
        return None;
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
fn ec_for_read_kmer_index(index: &KmerEcIndex, seq: &[u8]) -> Option<Vec<u32>> {
    if seq.len() < index.k {
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut has_hit = false;
    for pos in 0..=seq.len() - index.k {
        if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + index.k])
            && let Some(ec) = super::lookup_ec(index, fwd).or_else(|| super::lookup_ec(index, rev))
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                intersect_sorted(&current, ec, &mut next);
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

fn jump_distance_for_match(
    index: &BifrostIndex,
    unitig_id: usize,
    block_idx: usize,
    start: usize,
    forward_strand: bool,
) -> Option<usize> {
    let blocks = index.ec_blocks.get(unitig_id)?;
    let block = blocks.get(block_idx)?;
    if block.ub <= block.lb {
        return None;
    }
    let contig_start = block.lb as isize;
    let contig_len = (block.ub - block.lb) as isize;
    let start = start as isize;
    let dist = if forward_strand {
        contig_len - 1 - (start - contig_start)
    } else {
        start - contig_start
    };
    if dist >= 2 { Some(dist as usize) } else { None }
}

pub(super) fn filter_ec_by_fragment(
    index: &BifrostIndex,
    ec: &[u32],
    best_match: MatchInfo,
    fragment_length: i64,
) -> Vec<u32> {
    if fragment_length <= 0 {
        return ec.to_vec();
    }
    if index.transcript_lengths.is_empty() {
        return ec.to_vec();
    }
    if best_match.unitig_id >= index.ec_blocks.len() {
        return ec.to_vec();
    }
    let blocks = &index.ec_blocks[best_match.unitig_id];
    if blocks.is_empty() || blocks[0].positions.is_none() {
        return ec.to_vec();
    }
    let unitig_len = if best_match.unitig_id < index.unitigs.len() {
        index.unitigs[best_match.unitig_id].len()
    } else {
        index.k
    };
    let mut filtered = Vec::new();
    for &tr in ec {
        let (pos, forward) = match find_position_in_transcript(
            blocks,
            tr,
            best_match.unitig_pos,
            best_match.read_pos,
            best_match.used_revcomp,
            unitig_len,
            index.k,
        ) {
            Some((pos, forward)) => (pos, forward),
            None => {
                filtered.push(tr);
                continue;
            }
        };
        let len = index
            .transcript_lengths
            .get(tr as usize)
            .copied()
            .unwrap_or(0) as i64;
        if forward {
            if pos + fragment_length - 1 <= len {
                filtered.push(tr);
            }
        } else if pos - fragment_length >= 0 {
            filtered.push(tr);
        }
    }
    filtered
}

pub(super) fn apply_strand_filter(
    index: &BifrostIndex,
    read_ec: &ReadEc,
    mode: StrandSpecific,
    is_first_read: bool,
    comprehensive: bool,
    report: &mut Option<&mut DebugReport>,
    header: &[u8],
) -> Option<Vec<u32>> {
    let ec = &read_ec.ec;
    let target = match mode {
        StrandSpecific::FR => is_first_read,
        StrandSpecific::RF => !is_first_read,
    };

    if comprehensive {
        let mut union = Vec::new();
        for hit in &read_ec.hits {
            if let Some(filtered) = filter_ec_for_hit(index, ec, hit, target)
                && !filtered.is_empty()
            {
                super::merge_sorted_unique_vec(&mut union, &filtered);
            }
        }
        return Some(union);
    }

    let hit = read_ec.first_hit.as_ref()?;
    let filtered = filter_ec_for_hit(index, ec, hit, target)?;
    if filtered.len() < ec.len()
        && let Some(r) = report.as_deref_mut()
    {
        r.record(
            header,
            DebugFailReason::Unknown,
            None,
            None,
            None,
            None,
            None,
            None,
            read_ec
                .best_match
                .as_ref()
                .map(|m| m.used_revcomp)
                .unwrap_or(false),
        );
    }
    Some(filtered)
}

fn filter_ec_for_hit(
    index: &BifrostIndex,
    ec: &[u32],
    hit: &Hit,
    target: bool,
) -> Option<Vec<u32>> {
    let mut filtered = Vec::new();
    let block = index.ec_blocks.get(hit.unitig_id)?;
    let block = block.get(hit.block_idx)?;
    let strands = block.strands.as_ref()?;
    let um_strand = !hit.used_revcomp;
    for &tr in ec {
        let idx = match block.ec.binary_search(&tr) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let sense = strands.get(idx).copied().unwrap_or(2);
        if sense == 2 || ((um_strand == (sense == 1)) == target) {
            filtered.push(tr);
        }
    }
    if filtered.is_empty() {
        None
    } else {
        Some(filtered)
    }
}

pub(super) fn find_position_in_transcript(
    blocks: &[crate::index::EcBlock],
    tr: u32,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
    unitig_len: usize,
    k: usize,
) -> Option<(i64, bool)> {
    let idx = unitig_pos as u32;
    let mc = ec_block_at(blocks, idx)?;
    let ecs = ec_blocks_leading_vals(blocks, idx);
    if ecs.is_empty() {
        return None;
    }
    let v_ec = ecs.last().copied()?;
    let rawpos = block_min_pos(v_ec, tr)?;
    let trpos = (rawpos & 0x7fff_ffff) as i64;
    let trsense = rawpos == (rawpos & 0x7fff_ffff);

    let csense = !used_revcomp;
    let um_dist = unitig_pos as i64;
    let um_size = unitig_len as i64;
    let p = read_pos as i64;
    let k = k as i64;
    let mc_first = mc.0 as i64;
    let mc_second = mc.1 as i64;
    let _um_dist_block = um_dist - mc_first;

    if trsense {
        if csense {
            let mut padding = 0i64;
            if trpos == 0 {
                let mut mc_cur = mc;
                for block in ecs.iter().rev().skip(1) {
                    if !block_contains(block, tr) {
                        padding = mc_cur.0 as i64;
                        break;
                    }
                    mc_cur = prev_block_at(blocks, mc_cur.0)?;
                }
            }
            let pos = trpos - p + um_dist + 1 - padding;
            Some((pos, csense))
        } else {
            let mut mc_cur = mc;
            let mut right_one = 0i64;
            let mut left_one = 0i64;
            let initial = mc_second;
            for (i, block) in ecs.iter().enumerate().rev() {
                if i == ecs.len() - 1 {
                    right_one = mc_cur.1 as i64;
                }
                if !block_contains(block, tr) {
                    left_one = mc_cur.1 as i64;
                    break;
                } else if i == 0 {
                    left_one = 0;
                }
                mc_cur = prev_block_at(blocks, mc_cur.0)?;
            }
            let padding = -(left_one + right_one - um_size + k - 1);
            let pos = trpos + p + k - (um_size - k - um_dist) + initial - 1 + padding;
            Some((pos, csense))
        }
    } else if csense {
        let mut left_one = 0i64;
        let mut right_one = 0i64;
        let mut unmapped_len = 0i64;
        let mut found_first_mapped = false;
        let ecs_all = ec_blocks_leading_vals(blocks, u32::MAX);
        let mut curr_mc = 0u32;
        for block in ecs_all {
            let mc_ = ec_block_at(blocks, curr_mc)?;
            if !block_contains(block, tr) && found_first_mapped {
                if unmapped_len == 0 {
                    left_one = mc_.0 as i64;
                }
                right_one = mc_.1 as i64;
                unmapped_len += mc_.1 as i64 - mc_.0 as i64;
            }
            if block_contains(block, tr) {
                found_first_mapped = true;
            }
            curr_mc = mc_.1;
        }
        let mut start = 0i64;
        start -= right_one - left_one;
        start += um_size - k;
        let pos = trpos + (-(um_dist - start)) + k + p;
        Some((pos, !csense))
    } else {
        let mut left_one = 0i64;
        let mut right_one = 0i64;
        let mut unmapped_len = 0i64;
        let mut found_first_mapped = false;
        let ecs_all = ec_blocks_leading_vals(blocks, u32::MAX);
        let mut curr_mc = 0u32;
        for block in ecs_all {
            let mc_ = ec_block_at(blocks, curr_mc)?;
            if !block_contains(block, tr) && found_first_mapped {
                if unmapped_len == 0 {
                    left_one = mc_.0 as i64;
                }
                right_one = mc_.1 as i64;
                unmapped_len += mc_.1 as i64 - mc_.0 as i64;
            }
            if block_contains(block, tr) {
                found_first_mapped = true;
            }
            curr_mc = mc_.1;
        }
        let unmapped_len = right_one - left_one;
        let padding = um_size - um_dist - unmapped_len - k + 1;
        let pos = trpos + padding - p;
        Some((pos, !csense))
    }
}

pub(super) fn ec_blocks_leading_vals(
    blocks: &[crate::index::EcBlock],
    idx: u32,
) -> Vec<&crate::index::EcBlock> {
    if blocks.is_empty() {
        return Vec::new();
    }
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= idx {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    blocks[..lo].iter().collect()
}

pub(super) fn block_index_for_position(
    blocks: &[crate::index::EcBlock],
    pos: usize,
) -> Option<usize> {
    if blocks.is_empty() {
        return None;
    }
    if blocks.len() == 1 {
        return Some(0);
    }
    let pos = pos as u32;
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= pos {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        return None;
    }
    Some(lo - 1)
}

pub(super) fn ec_block_at(blocks: &[crate::index::EcBlock], idx: u32) -> Option<(u32, u32)> {
    if blocks.is_empty() {
        return None;
    }
    if blocks.len() == 1 {
        return Some((blocks[0].lb, blocks[0].ub));
    }
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= idx {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        return None;
    }
    let block = &blocks[lo - 1];
    Some((block.lb, block.ub))
}

pub(super) fn prev_block_at(blocks: &[crate::index::EcBlock], lb: u32) -> Option<(u32, u32)> {
    if lb == 0 {
        ec_block_at(blocks, u32::MAX)
    } else {
        ec_block_at(blocks, lb - 1)
    }
}

pub(super) fn block_contains(block: &crate::index::EcBlock, tr: u32) -> bool {
    block.ec.binary_search(&tr).is_ok()
}

pub(super) fn block_min_pos(block: &crate::index::EcBlock, tr: u32) -> Option<u32> {
    let positions = block.positions.as_ref()?;
    let idx = block.ec.binary_search(&tr).ok()?;
    positions.get(idx).copied()
}
