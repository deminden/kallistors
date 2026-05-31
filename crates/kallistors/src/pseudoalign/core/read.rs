use super::super::debug::{
    DroppedHitTrace, JumpDecisionTrace, MinimizerCandidateTrace, ReadDebugState,
    format_positions_sample,
};
use super::super::ec::{
    ec_has_offlist, encode_kmer_pair, filter_shades, intersect_sorted, merge_sorted_unique_vec,
};
use super::super::minimizers::{
    fill_revcomp, minhash_next_after_hash, minhash_primary_for_kmer, minimizer_for_kmer_strict,
    minimizer_tail_for_kmer, minimizers_for_kmer, minimizers_ranked_for_kmer,
};
use super::super::{BifrostIndex, Hit, MatchInfo, PseudoalignOptions, ReadEc, Strand};
use super::fast::ec_for_read_bifrost_fast;
use super::local::{collect_unitigs_for_read, ec_for_read_kmer_index, ec_for_read_local_kmer};
use super::matcher::{
    jump_distance_for_match, match_kmer_at_pos, match_kmer_direct, match_unitig_candidate,
    special_unitig_for_kmer, special_unitig_for_kmer_raw,
};
use super::runtime::{
    ACCEPT_RELAXED_STREAM, ALLOW_TAIL_MINIMIZER_IN_BACKOFF, BACKOFF_DIRECT_RESCUE,
    block_index_for_position_fast, ec_slice, env_flag, fast_path_enabled, minimizer_hex,
    onlist_cardinality,
};

// Main parity scanner: stream hits, keep an online EC intersection, and jump when safe.
pub(crate) fn ec_for_read_bifrost(
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
    if fast_path_enabled(index, dbg.is_some(), options) {
        return ec_for_read_bifrost_fast(index, seq, strand, dbg, options);
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
    let mut next_shade_scratch: Vec<u32> = Vec::new();
    let mut online_intersection: Vec<u32> = Vec::new();
    let mut has_online_intersection = false;
    let mut saw_special_hit = false;
    let mut saw_non_special_hit = false;
    let dlist_dummy = index.onlist.as_deref().map(onlist_cardinality);
    let mut any_dlist_kmer = false;
    let mut dlist_early_return = false;
    let mut append_dummy_hit = false;
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
                        false,
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
                && !tried_tail_minimizer
                && (!in_backoff
                    || env_flag(
                        &ALLOW_TAIL_MINIMIZER_IN_BACKOFF,
                        "KALLISTORS_ALLOW_TAIL_MINIMIZER_IN_BACKOFF",
                    ))
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
        if !accept_for_stream && in_backoff {
            let can_terminal_jump =
                jump_distance_for_match(index, uid, block_idx, start, forward_strand)
                    .is_some_and(|dist| pos.saturating_add(dist) >= last_pos);
            if !can_terminal_jump {
                pos += 1;
                continue;
            }
        }

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
                        append_dummy_hit = dlist_dummy.is_some();
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
    if !dlist_early_return && any_dlist_kmer && (has_hit || options.do_union) {
        append_dummy_hit = dlist_dummy.is_some();
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
                if let Some(state) = dbg {
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
    let mut last_unitig = None;
    let mut last_block_idx = None;
    let mut found_nonempty = false;
    let mut current_offlist = false;
    let mut ec_scratch: Vec<u32> = Vec::new();
    for hit in &hits {
        minpos = minpos.min(hit.read_pos);
        maxpos = maxpos.max(hit.read_pos);
        let ec = ec_slice(index, hit.unitig_id, hit.block_idx);
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
                    last_unitig = Some(hit.unitig_id);
                    last_block_idx = Some(hit.block_idx);
                }
                current_offlist = ec_offlist;
                found_nonempty = true;
            }
            continue;
        }
        if last_unitig == Some(hit.unitig_id) && last_block_idx == Some(hit.block_idx) {
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
                last_unitig = Some(hit.unitig_id);
                last_block_idx = Some(hit.block_idx);
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
                if let Some(state) = dbg {
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
            last_unitig = Some(hit.unitig_id);
            last_block_idx = Some(hit.block_idx);
            current_offlist = current_offlist || ec_offlist;
        }
    }

    // Apply synthetic D-list dummy hits (kallisto um_dummy behavior).
    if append_dummy_hit && let Some(dummy) = dlist_dummy {
        let ec = [dummy];
        if options.do_union {
            merge_sorted_unique_vec(&mut current, &ec);
            current_offlist = true;
        } else if current.is_empty() {
            current.extend_from_slice(&ec);
            current_offlist = true;
        } else {
            next.clear();
            intersect_sorted(&current, &ec, &mut next);
            std::mem::swap(&mut current, &mut next);
            if current.is_empty() {
                if let Some(state) = dbg {
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
                let ec = ec_slice(index, hit.unitig_id, hit.block_idx);
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
