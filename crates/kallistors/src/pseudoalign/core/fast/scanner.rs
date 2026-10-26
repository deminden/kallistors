//! Ordered read scanning, jump decisions and online EC folding.

use super::*;

pub(super) fn ec_for_read_bifrost_fast_with_cache<
    'index,
    const TRACE: bool,
    const COUNT_ONLY: bool,
    const PLACEMENT: bool,
>(
    index: &'index BifrostIndex,
    seq: &[u8],
    strand: Strand,
    dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
    cache: &mut MphfLookupCache,
    scratch: &mut FastScannerScratch,
) -> Option<BorrowedReadEc<'index>> {
    if index.has_large_candidate_cache() {
        ec_for_read_bifrost_fast_with_cache_cached::<TRACE, COUNT_ONLY, PLACEMENT, true>(
            index, seq, strand, dbg, options, cache, scratch,
        )
    } else {
        ec_for_read_bifrost_fast_with_cache_cached::<TRACE, COUNT_ONLY, PLACEMENT, false>(
            index, seq, strand, dbg, options, cache, scratch,
        )
    }
}

fn ec_for_read_bifrost_fast_with_cache_cached<
    'index,
    const TRACE: bool,
    const COUNT_ONLY: bool,
    const PLACEMENT: bool,
    const CACHED: bool,
>(
    index: &'index BifrostIndex,
    seq: &[u8],
    strand: Strand,
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
    cache: &mut MphfLookupCache,
    scratch: &mut FastScannerScratch,
) -> Option<BorrowedReadEc<'index>> {
    let mut current: Vec<u32> = Vec::new();
    let min_candidates = &mut scratch.min_candidates;
    let jump_candidates = &mut scratch.jump_candidates;
    let allow_forward_base = strand != Strand::Reverse;
    let allow_rev_base = strand != Strand::Forward;
    let mut has_hit = false;
    let diff = index.k.saturating_sub(index.g);
    let mut best_match: Option<MatchInfo> = None;
    let mut first_hit: Option<Hit> = None;
    let mut hits: Vec<Hit> = Vec::new();
    let mut counted = super::super::count::FinalHitFold::default();
    // The folds share an exact prefix until a recorded constraint was omitted
    // from the online stream. Materialize that prefix only on divergence.
    let mut counted_diverged = false;
    // Ordinary reads only need to know whether an on-list exists. Counting every
    // transcript here would rescan the whole reference for every read.
    let dlist_onlist = index.onlist.as_deref();
    let mut any_dlist_kmer = false;
    let mut dlist_early_return = false;
    let mut append_dummy_hit = false;
    let online_intersection = &mut OnlineIntersection {
        initial: None,
        owned: &mut scratch.online_intersection,
        pairs: scratch.pairs.as_deref_mut(),
    };
    let mut has_online_intersection = false;
    let mut last_online_block = None;
    let mut backoff_until: Option<usize> = None;
    let mut fallback_probe_pending = false;
    let minimizer_cache_disabled = env_flag(
        &DISABLE_MINIMIZER_CACHE,
        "KALLISTORS_DISABLE_MINIMIZER_CACHE",
    );
    let match_cache_disabled = env_flag(&DISABLE_MATCH_CACHE, "KALLISTORS_DISABLE_MATCH_CACHE");
    let negative_cache_enabled =
        dbg.is_none() && !options.investigation.any_enabled() && all_ascii_uppercase(seq);
    let mut kmer_cursor = KmerPairCursor::new(index.k);
    let mut pos = 0usize;
    let last_pos = seq.len() - index.k;
    while pos <= last_pos {
        if TRACE && let Some(state) = dbg.as_deref_mut() {
            state.saw_valid_kmer = true;
            if state.visited_positions.len() < 4096 {
                state.visited_positions.push(pos);
            }
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
        let kmer_codes = kmer_cursor.encode(kmer, pos);
        let kmer_in_dlist = index
            .dlist
            .as_ref()
            .zip(kmer_codes)
            .is_some_and(|(dlist, (fwd, rev))| dlist.contains(&fwd.min(rev)));
        if kmer_in_dlist {
            any_dlist_kmer = true;
        }
        if !allow_forward_base && !allow_rev_base {
            pos += 1;
            continue;
        }
        let Some((read_fwd, read_rev)) = kmer_codes else {
            pos += 1;
            continue;
        };
        let fallback_probe = std::mem::take(&mut fallback_probe_pending);
        if dlist_onlist.is_some() && kmer_in_dlist {
            append_dummy_hit = true;
        }
        if TRACE && let Some(state) = dbg.as_deref_mut() {
            if !minimizer_candidates_cached_into_with_code(
                cache,
                kmer,
                index.g,
                read_fwd,
                read_rev,
                min_candidates,
                minimizer_cache_disabled,
            ) {
                pos += 1;
                continue;
            }
            for (min_bytes, min_pos) in min_candidates.iter() {
                let bucket = minimizer_bucket_idx(index, cache, min_bytes);
                let positions_len = bucket
                    .map(|idx| index.minz_positions.get(idx).len())
                    .unwrap_or(0);
                let sample_positions = bucket
                    .map(|idx| format_positions_sample(index.minz_positions.get(idx)))
                    .unwrap_or_else(|| "-".to_string());
                if bucket.is_some() {
                    state.saw_mphf_hit = true;
                    if positions_len > 0 {
                        state.saw_positions = true;
                    }
                }
                state.minimizer_candidates.push(MinimizerCandidateTrace {
                    read_pos: pos,
                    min_pos: *min_pos,
                    minimizer: minimizer_hex(*min_bytes),
                    mphf_hit: bucket.is_some(),
                    positions_len,
                    sample_positions,
                    has_special: false,
                    overcrowded: positions_len > 64,
                    matched: false,
                });
            }
        }
        let matched = if let Some(matched) = match_kmer_at_pos_fast_with_codes::<CACHED>(
            index,
            kmer,
            read_fwd,
            read_rev,
            allow_forward_base,
            allow_rev_base,
            diff,
            min_candidates,
            cache,
            true,
            !in_backoff,
            minimizer_cache_disabled,
            match_cache_disabled || TRACE,
            false,
            negative_cache_enabled,
        ) {
            matched
        } else if !options.kallisto_strict
            && !options.kallisto_enum
            && in_backoff
            && min_candidates.len() == 1
            && let Some(matched) = super::rescue::match_bifrost_rescue::<CACHED>(
                index,
                kmer,
                read_fwd,
                read_rev,
                allow_forward_base,
                allow_rev_base,
                diff,
                cache,
                true,
            )
        {
            matched
        } else {
            if TRACE
                && let Some(state) = dbg.as_deref_mut()
                && state.first_no_match.is_none()
            {
                state.first_no_match = Some((pos, 0));
            }
            pos += 1;
            continue;
        };
        let uid = matched.unitig_id;
        let start = matched.start;
        let kmer_pos = pos;
        let block_idx = matched.block_idx;
        let used_revcomp = matched.used_revcomp;
        let forward_strand = matched.forward_strand;
        let matched_relaxed = matched.matched_relaxed;
        let ec = ec_slice(index, uid, block_idx);
        if TRACE && let Some(state) = dbg.as_deref_mut() {
            state.saw_match = true;
            if ec.is_empty() {
                if state.first_empty_ec.is_none() {
                    state.first_empty_ec = Some((pos, 0));
                }
            } else {
                state.saw_ec = true;
            }
            if let Some(entry) = state
                .minimizer_candidates
                .iter_mut()
                .rev()
                .find(|entry| entry.read_pos == pos)
            {
                entry.matched = true;
            }
            state.used_revcomp = used_revcomp;
        }
        let mut accept_for_stream = true;
        if matched_relaxed && has_hit {
            accept_for_stream = false;
            if TRACE
                && let Some(state) = dbg.as_deref_mut()
                && state.dropped_hits.len() < 256
            {
                state.dropped_hits.push(DroppedHitTrace {
                    read_pos: kmer_pos,
                    min_pos: 0,
                    unitig_id: uid,
                    unitig_pos: start,
                    block_idx: Some(block_idx),
                    reason: if in_backoff {
                        "backoff_relaxed_probe"
                    } else {
                        "relaxed_probe"
                    },
                    used_revcomp,
                    is_special: false,
                });
            }
        }
        if accept_for_stream && !ec.is_empty() && last_online_block != Some((uid, block_idx)) {
            if has_online_intersection {
                online_intersection.intersect(ec);
                if online_intersection.is_empty() {
                    return Some(BorrowedReadEc::owned(ReadEc {
                        ec: Vec::new(),
                        best_match,
                        first_hit,
                        hits,
                        had_offlist: false,
                        shade_union: Vec::new(),
                        hard_reject_pair: true,
                    }));
                }
            } else {
                online_intersection.initialize(ec);
                has_online_intersection = true;
            }
            last_online_block = Some((uid, block_idx));
        }
        if PLACEMENT && best_match.map(|m| kmer_pos < m.read_pos).unwrap_or(true) {
            best_match = Some(MatchInfo {
                unitig_id: uid,
                unitig_pos: start,
                read_pos: kmer_pos,
                used_revcomp,
            });
        }
        let hit = Hit {
            unitig_id: uid,
            read_pos: kmer_pos,
            block_idx,
            used_revcomp,
        };
        if COUNT_ONLY {
            if counted_diverged {
                counted.push(index, hit);
            } else if !accept_for_stream
                && !ec.is_empty()
                && last_online_block != Some((uid, block_idx))
            {
                counted = super::super::count::FinalHitFold::from_online(
                    online_intersection,
                    last_online_block,
                );
                counted.push(index, hit);
                counted_diverged = true;
            }
            if counted.rejected {
                return Some(BorrowedReadEc::owned(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: false,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                }));
            }
        } else {
            hits.push(hit);
        }
        if PLACEMENT
            && first_hit
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
        // Kallisto probes the next valid k-mer after a failed jump before returning
        // to its outer jump loop. That probe cannot itself initiate a jump.
        if fallback_probe {
            pos += 1;
            continue;
        }
        // Excluding a probe from EC intersection does not invalidate its exact
        // sequence match for confirming a jump, including on reverse unitigs.
        let mut jump_to: Option<usize> = None;
        let mut force_break = false;
        let mut synthetic_hit: Option<Hit> = None;
        if let Some(dist) = jump_distance_for_match(index, uid, block_idx, start, forward_strand) {
            let mut next_pos = pos + dist;
            if next_pos > last_pos {
                next_pos = last_pos;
            }
            if next_pos > pos {
                let Some((next_query_pos, next_codes)) =
                    encoded_kmer_at_or_after(seq, next_pos, index.k)
                else {
                    break;
                };
                let kmer_next = &seq[next_query_pos..next_query_pos + index.k];
                let kmer_next_codes = Some(next_codes);
                let mut next_hit = kmer_next_codes.and_then(|(next_fwd, next_rev)| {
                    match_kmer_at_pos_fast_with_codes::<CACHED>(
                        index,
                        kmer_next,
                        next_fwd,
                        next_rev,
                        allow_forward_base,
                        allow_rev_base,
                        diff,
                        jump_candidates,
                        cache,
                        true,
                        !in_backoff,
                        minimizer_cache_disabled,
                        match_cache_disabled,
                        true,
                        negative_cache_enabled,
                    )
                });
                if next_hit.is_none() && in_backoff {
                    next_hit = super::rescue::match_ordinary_tail_after_primary_miss::<CACHED>(
                        index,
                        kmer_next,
                        next_codes,
                        allow_forward_base,
                        allow_rev_base,
                        diff,
                        jump_candidates,
                        cache,
                        true,
                        options.skip_overcrowded_minimizer,
                        minimizer_cache_disabled,
                    );
                }
                let current_ec = ec_slice(index, uid, block_idx);
                let mut use_backoff = false;
                if let Some(next_hit) = next_hit {
                    let uid2 = next_hit.unitig_id;
                    let block_idx2 = next_hit.block_idx;
                    let next_relaxed = next_hit.matched_relaxed;
                    let next_ec = ec_slice(index, uid2, block_idx2);
                    if uid2 == uid && next_ec == current_ec {
                        if TRACE && let Some(state) = dbg.as_deref_mut() {
                            state.jump_decisions.push(JumpDecisionTrace {
                                read_pos: pos,
                                matched_unitig_id: uid,
                                matched_block_idx: block_idx,
                                jump_distance: dist,
                                next_pos,
                                in_backoff,
                                next_hit_found: true,
                                next_hit_relaxed: next_relaxed,
                                next_hit_same_unitig: true,
                                next_hit_same_ec: true,
                                mid_hit_found: false,
                                mid_hit_relaxed: false,
                                mid_hit_matches_either: false,
                                jumped: true,
                                reason: "next_same_unitig_same_ec",
                            });
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
                            synthetic_hit = Some(Hit {
                                unitig_id: uid,
                                read_pos: next_pos,
                                block_idx,
                                used_revcomp,
                            });
                            jump_to = Some(next_query_pos);
                        }
                    } else if dist > 4 {
                        let middle_pos = (pos + next_pos) / 2;
                        if let Some((middle_query_pos, mid_codes)) =
                            encoded_kmer_at_or_after(seq, middle_pos, index.k)
                        {
                            let kmer_mid = &seq[middle_query_pos..middle_query_pos + index.k];
                            let mut mid_hit = Some(mid_codes).and_then(|(mid_fwd, mid_rev)| {
                                match_kmer_at_pos_fast_with_codes::<CACHED>(
                                    index,
                                    kmer_mid,
                                    mid_fwd,
                                    mid_rev,
                                    allow_forward_base,
                                    allow_rev_base,
                                    diff,
                                    jump_candidates,
                                    cache,
                                    true,
                                    !in_backoff,
                                    minimizer_cache_disabled,
                                    match_cache_disabled,
                                    true,
                                    negative_cache_enabled,
                                )
                            });
                            if mid_hit.is_none() && in_backoff {
                                mid_hit = super::rescue::match_ordinary_tail_after_primary_miss::<
                                    CACHED,
                                >(
                                    index,
                                    kmer_mid,
                                    mid_codes,
                                    allow_forward_base,
                                    allow_rev_base,
                                    diff,
                                    jump_candidates,
                                    cache,
                                    true,
                                    options.skip_overcrowded_minimizer,
                                    minimizer_cache_disabled,
                                );
                            }
                            if let Some(mid_hit) = mid_hit {
                                let uid3 = mid_hit.unitig_id;
                                let rev3 = mid_hit.used_revcomp;
                                let block_idx3 = mid_hit.block_idx;
                                let mid_relaxed = mid_hit.matched_relaxed;
                                let mid_ec = ec_slice(index, uid3, block_idx3);
                                let mid_matches_current = uid3 == uid && mid_ec == current_ec;
                                let mid_matches_next = uid3 == uid2 && mid_ec == next_ec;
                                if mid_matches_current || mid_matches_next {
                                    if TRACE && let Some(state) = dbg.as_deref_mut() {
                                        state.jump_decisions.push(JumpDecisionTrace {
                                            read_pos: pos,
                                            matched_unitig_id: uid,
                                            matched_block_idx: block_idx,
                                            jump_distance: dist,
                                            next_pos,
                                            in_backoff,
                                            next_hit_found: true,
                                            next_hit_relaxed: next_relaxed,
                                            next_hit_same_unitig: uid2 == uid,
                                            next_hit_same_ec: next_ec == current_ec,
                                            mid_hit_found: true,
                                            mid_hit_relaxed: mid_relaxed,
                                            mid_hit_matches_either: true,
                                            jumped: true,
                                            reason: "middle_confirms_jump",
                                        });
                                    }
                                    if !mid_ec.is_empty()
                                        && last_online_block != Some((uid3, block_idx3))
                                    {
                                        if has_online_intersection {
                                            online_intersection.intersect(mid_ec);
                                            if online_intersection.is_empty() {
                                                return Some(BorrowedReadEc::owned(ReadEc {
                                                    ec: Vec::new(),
                                                    best_match,
                                                    first_hit,
                                                    hits,
                                                    had_offlist: false,
                                                    shade_union: Vec::new(),
                                                    hard_reject_pair: true,
                                                }));
                                            }
                                        } else {
                                            online_intersection.initialize(mid_ec);
                                            has_online_intersection = true;
                                        }
                                        last_online_block = Some((uid3, block_idx3));
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
                                    } else {
                                        jump_to = Some(next_query_pos);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    let next_in_dlist = index
                        .dlist
                        .as_ref()
                        .zip(kmer_next_codes)
                        .is_some_and(|(dlist, (fwd, rev))| dlist.contains(&fwd.min(rev)));
                    if next_in_dlist {
                        dlist_early_return = true;
                        append_dummy_hit = dlist_onlist.is_some();
                    }
                    if TRACE && let Some(state) = dbg.as_deref_mut() {
                        state.jump_decisions.push(JumpDecisionTrace {
                            read_pos: pos,
                            matched_unitig_id: uid,
                            matched_block_idx: block_idx,
                            jump_distance: dist,
                            next_pos,
                            in_backoff,
                            next_hit_found: false,
                            next_hit_relaxed: false,
                            next_hit_same_unitig: false,
                            next_hit_same_ec: false,
                            mid_hit_found: false,
                            mid_hit_relaxed: false,
                            mid_hit_matches_either: false,
                            jumped: !dlist_early_return,
                            reason: if dlist_early_return {
                                "next_miss_dlist_early_return"
                            } else {
                                "next_miss_jump_or_break"
                            },
                        });
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
                    if options.investigation.bounded_incremental_scan {
                        let stop = next_pos.min(last_pos);
                        let mut recovered = false;
                        let mut scan_pos = pos + 1;
                        while scan_pos <= stop {
                            let scan_kmer = &seq[scan_pos..scan_pos + index.k];
                            let scan_hit =
                                encode_kmer_pair(scan_kmer).and_then(|(scan_fwd, scan_rev)| {
                                    match_kmer_at_pos_fast_with_codes::<CACHED>(
                                        index,
                                        scan_kmer,
                                        scan_fwd,
                                        scan_rev,
                                        allow_forward_base,
                                        allow_rev_base,
                                        diff,
                                        jump_candidates,
                                        cache,
                                        options.investigation.candidate_reuse,
                                        !in_backoff,
                                        minimizer_cache_disabled,
                                        match_cache_disabled,
                                        true,
                                        negative_cache_enabled,
                                    )
                                });
                            if let Some(scan_hit) = scan_hit {
                                synthetic_hit = Some(Hit {
                                    unitig_id: scan_hit.unitig_id,
                                    read_pos: scan_pos,
                                    block_idx: scan_hit.block_idx,
                                    used_revcomp: scan_hit.used_revcomp,
                                });
                                jump_to = Some(scan_pos);
                                recovered = true;
                                if TRACE && let Some(state) = dbg.as_deref_mut() {
                                    state.jump_decisions.push(JumpDecisionTrace {
                                        read_pos: pos,
                                        matched_unitig_id: uid,
                                        matched_block_idx: block_idx,
                                        jump_distance: dist,
                                        next_pos: scan_pos,
                                        in_backoff,
                                        next_hit_found: false,
                                        next_hit_relaxed: false,
                                        next_hit_same_unitig: scan_hit.unitig_id == uid,
                                        next_hit_same_ec: false,
                                        mid_hit_found: false,
                                        mid_hit_relaxed: false,
                                        mid_hit_matches_either: false,
                                        jumped: true,
                                        reason: "bounded_incremental_recovery",
                                    });
                                }
                                break;
                            }
                            scan_pos += 1;
                        }
                        if !recovered {
                            use_backoff = true;
                        }
                    } else {
                        use_backoff = true;
                    }
                }
                if let Some(hit) = synthetic_hit.take() {
                    if COUNT_ONLY {
                        if counted_diverged {
                            counted.push(index, hit);
                        } else if last_online_block != Some((hit.unitig_id, hit.block_idx))
                            && !ec_slice(index, hit.unitig_id, hit.block_idx).is_empty()
                        {
                            counted = super::super::count::FinalHitFold::from_online(
                                online_intersection,
                                last_online_block,
                            );
                            counted.push(index, hit);
                            counted_diverged = true;
                        }
                        if counted.rejected {
                            return Some(BorrowedReadEc::owned(ReadEc {
                                ec: Vec::new(),
                                best_match,
                                first_hit,
                                hits,
                                had_offlist: false,
                                shade_union: Vec::new(),
                                hard_reject_pair: true,
                            }));
                        }
                    } else {
                        hits.push(hit);
                    }
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

    if !dlist_early_return && any_dlist_kmer && has_hit {
        append_dummy_hit = dlist_onlist.is_some();
    }
    if !has_hit || (!COUNT_ONLY && hits.is_empty()) {
        return None;
    }

    let mut minpos = usize::MAX;
    let mut maxpos = 0usize;
    let mut last_unitig = None;
    let mut last_block_idx = None;
    let mut found_nonempty = false;
    if COUNT_ONLY {
        if counted.rejected {
            return Some(BorrowedReadEc::owned(ReadEc {
                ec: Vec::new(),
                best_match,
                first_hit,
                hits,
                had_offlist: false,
                shade_union: Vec::new(),
                hard_reject_pair: true,
            }));
        }
        if !counted_diverged
            && !append_dummy_hit
            && let Some(ec) = online_intersection.initial
            && !ec.is_empty()
            && index
                .onlist
                .as_deref()
                .is_none_or(|mask| ec.iter().all(|&id| mask.get(id as usize) == Some(&true)))
        {
            return Some(BorrowedReadEc {
                borrowed: Some(ec),
                read: ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: false,
                    shade_union: Vec::new(),
                    hard_reject_pair: false,
                },
            });
        }
        if counted_diverged {
            current = counted.ec;
        } else {
            current.extend_from_slice(online_intersection);
        }
    } else {
        for hit in &hits {
            minpos = minpos.min(hit.read_pos);
            maxpos = maxpos.max(hit.read_pos);
            let ec = ec_slice(index, hit.unitig_id, hit.block_idx);
            if !found_nonempty {
                if !ec.is_empty() {
                    current.extend_from_slice(ec);
                    last_unitig = Some(hit.unitig_id);
                    last_block_idx = Some(hit.block_idx);
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
            if current.is_empty() {
                current.extend_from_slice(ec);
                last_unitig = Some(hit.unitig_id);
                last_block_idx = Some(hit.block_idx);
                continue;
            }
            intersect_sorted_in_place(&mut current, ec);
            if current.is_empty() {
                return Some(BorrowedReadEc::owned(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: false,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                }));
            }
            last_unitig = Some(hit.unitig_id);
            last_block_idx = Some(hit.block_idx);
        }
    }

    if append_dummy_hit && let Some(onlist) = dlist_onlist {
        let dummy = onlist_cardinality(onlist);
        let ec = [dummy];
        if current.is_empty() {
            current.extend_from_slice(&ec);
        } else {
            intersect_sorted_in_place(&mut current, &ec);
            if current.is_empty() {
                return Some(BorrowedReadEc::owned(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: true,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                }));
            }
        }
    }

    if current.is_empty() {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
    }
    if current.is_empty() {
        return None;
    }
    let min_range = 1usize;
    if maxpos >= minpos && (maxpos - minpos + index.k) < min_range {
        return None;
    }
    Some(BorrowedReadEc::owned(ReadEc {
        ec: current,
        best_match,
        first_hit,
        hits,
        had_offlist: false,
        shade_union: Vec::new(),
        hard_reject_pair: false,
    }))
}
