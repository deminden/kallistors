use super::super::debug::{
    DroppedHitTrace, JumpDecisionTrace, MinimizerCandidateTrace, ReadDebugState,
    format_positions_sample,
};
use super::super::ec::{encode_kmer_pair, intersect_sorted};
use super::super::minimizers::minimizer_tail_for_kmer;
use super::super::{BifrostIndex, Hit, MatchInfo, PseudoalignOptions, ReadEc, Strand};
use super::matcher::{
    jump_distance_for_match, match_kmer_at_pos, match_kmer_direct, match_unitig_candidate_encoded,
    special_unitig_for_kmer,
};
use super::runtime::{
    DISABLE_MATCH_CACHE, DISABLE_MINIMIZER_CACHE, FastKmerMatch, MphfLookupCache,
    block_index_for_position_fast, ec_slice, env_flag, minimizer_bucket_has_overcrowded_marker,
    minimizer_bucket_idx, minimizer_candidates_cached_into_with_code, minimizer_hex,
    onlist_cardinality, read_fast_match_cache, write_fast_match_cache,
};

// Encoded fast path for the common unshaded, non-union pseudoalignment mode.
#[allow(clippy::too_many_arguments)]
fn match_kmer_at_pos_fast_uncached_precomputed(
    index: &BifrostIndex,
    kmer: &[u8],
    read_fwd: u64,
    read_rev: u64,
    allow_forward: bool,
    allow_rev: bool,
    diff: usize,
    min_candidates: &[([u8; 8], usize)],
    next_min_candidate: Option<([u8; 8], usize)>,
    cache: &mut MphfLookupCache,
    allow_relaxed: bool,
) -> Option<FastKmerMatch> {
    let use_revcomp = read_rev < read_fwd;
    let mut saw_overcrowded = false;
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
        let Some(bucket_idx) = minimizer_bucket_idx(index, cache, &min_bytes) else {
            continue;
        };
        for &pos_id in index.minz_positions.get(bucket_idx) {
            let unitig_id_raw = (pos_id >> 32) as u32;
            if unitig_id_raw == u32::MAX {
                if (pos_id & 0x8000_0000) != 0 {
                    saw_overcrowded = true;
                    continue;
                }
                if (pos_id & 0xffff_ffff) != 0
                    && let Some(uid) = special_unitig_for_kmer(index, kmer)
                {
                    let used_revcomp = use_revcomp;
                    if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                        continue;
                    }
                    let block_idx = block_index_for_position_fast(index, uid, 0)?;
                    return Some(FastKmerMatch {
                        unitig_id: uid,
                        start: 0,
                        block_idx,
                        used_revcomp,
                        forward_strand: !used_revcomp,
                        matched_relaxed: false,
                    });
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
                let kmer_canon = if use_revcomp { read_rev } else { read_fwd };
                if index
                    .encoded_unitigs
                    .km_bases
                    .encode_kmer_at(uid, 0, index.k)
                    != Some(kmer_canon)
                {
                    continue;
                }
                if min_pos != km_pos && min_pos + km_pos != diff {
                    continue;
                }
                let used_revcomp = use_revcomp;
                if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                    continue;
                }
                let unitig_id = index.unitigs.len() + uid;
                let block_idx = block_index_for_position_fast(index, unitig_id, 0)?;
                return Some(FastKmerMatch {
                    unitig_id,
                    start: 0,
                    block_idx,
                    used_revcomp,
                    forward_strand: !used_revcomp,
                    matched_relaxed: false,
                });
            }
            let uid = unitig_id as usize;
            if uid >= index.unitigs.len() {
                continue;
            }
            if let Some((start, used_revcomp, matched_relaxed, forward_strand)) =
                match_unitig_candidate_encoded(
                    index,
                    uid,
                    rel_pos as usize,
                    diff,
                    min_pos_fwd,
                    min_pos_rev,
                    read_fwd,
                    read_rev,
                    allow_forward,
                    allow_rev,
                    allow_relaxed,
                )
            {
                let block_idx = block_index_for_position_fast(index, uid, start)?;
                return Some(FastKmerMatch {
                    unitig_id: uid,
                    start,
                    block_idx,
                    used_revcomp,
                    forward_strand,
                    matched_relaxed,
                });
            }
        }
    }
    if saw_overcrowded
        && let Some((next_bytes, next_pos)) = next_min_candidate
        && !min_candidates
            .iter()
            .any(|&(bytes, pos)| bytes == next_bytes && pos == next_pos)
        && let Some(bucket_idx) = minimizer_bucket_idx(index, cache, &next_bytes)
    {
        let min_pos_fwd = if use_revcomp {
            diff.saturating_sub(next_pos)
        } else {
            next_pos
        };
        let min_pos_rev = if use_revcomp {
            next_pos
        } else {
            diff.saturating_sub(next_pos)
        };
        for &pos_id in index.minz_positions.get(bucket_idx) {
            let unitig_id_raw = (pos_id >> 32) as u32;
            if unitig_id_raw == u32::MAX {
                continue;
            }
            let (unitig_id, rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
            if is_km {
                let uid = unitig_id as usize;
                if uid >= index.km_unitigs.len() {
                    continue;
                }
                let km_pos = rel_pos as usize;
                let kmer_canon = if use_revcomp { read_rev } else { read_fwd };
                if index
                    .encoded_unitigs
                    .km_bases
                    .encode_kmer_at(uid, 0, index.k)
                    != Some(kmer_canon)
                {
                    continue;
                }
                if next_pos != km_pos && next_pos + km_pos != diff {
                    continue;
                }
                let used_revcomp = use_revcomp;
                if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                    continue;
                }
                let unitig_id = index.unitigs.len() + uid;
                let block_idx = block_index_for_position_fast(index, unitig_id, 0)?;
                return Some(FastKmerMatch {
                    unitig_id,
                    start: 0,
                    block_idx,
                    used_revcomp,
                    forward_strand: !used_revcomp,
                    matched_relaxed: false,
                });
            }
            let uid = unitig_id as usize;
            if uid >= index.unitigs.len() {
                continue;
            }
            if let Some((start, used_revcomp, matched_relaxed, forward_strand)) =
                match_unitig_candidate_encoded(
                    index,
                    uid,
                    rel_pos as usize,
                    diff,
                    min_pos_fwd,
                    min_pos_rev,
                    read_fwd,
                    read_rev,
                    allow_forward,
                    allow_rev,
                    allow_relaxed,
                )
            {
                let block_idx = block_index_for_position_fast(index, uid, start)?;
                return Some(FastKmerMatch {
                    unitig_id: uid,
                    start,
                    block_idx,
                    used_revcomp,
                    forward_strand,
                    matched_relaxed,
                });
            }
        }
    }
    if saw_overcrowded
        && let Some((unitig_id, start, used_revcomp)) =
            match_kmer_direct(index, kmer, allow_forward, allow_rev)
    {
        let block_idx = block_index_for_position_fast(index, unitig_id, start)?;
        return Some(FastKmerMatch {
            unitig_id,
            start,
            block_idx,
            used_revcomp,
            forward_strand: !used_revcomp,
            matched_relaxed: false,
        });
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn match_kmer_at_pos_fast_uncached_with_codes(
    index: &BifrostIndex,
    kmer: &[u8],
    read_fwd: u64,
    read_rev: u64,
    allow_forward: bool,
    allow_rev: bool,
    diff: usize,
    min_candidates: &mut Vec<([u8; 8], usize)>,
    cache: &mut MphfLookupCache,
    allow_relaxed: bool,
    allow_tail_minimizer: bool,
    minimizer_cache_disabled: bool,
) -> Option<FastKmerMatch> {
    if !minimizer_candidates_cached_into_with_code(
        kmer,
        index.g,
        read_fwd,
        min_candidates,
        minimizer_cache_disabled,
    ) {
        return None;
    }
    let tail_candidate = minimizer_tail_for_kmer(kmer, index.g);
    match_kmer_at_pos_fast_uncached_with_candidates(
        index,
        kmer,
        read_fwd,
        read_rev,
        allow_forward,
        allow_rev,
        diff,
        min_candidates,
        tail_candidate,
        cache,
        allow_relaxed,
        allow_tail_minimizer,
    )
}

#[allow(clippy::too_many_arguments)]
fn match_kmer_at_pos_fast_uncached_with_candidates(
    index: &BifrostIndex,
    kmer: &[u8],
    read_fwd: u64,
    read_rev: u64,
    allow_forward: bool,
    allow_rev: bool,
    diff: usize,
    min_candidates: &[([u8; 8], usize)],
    tail_candidate: Option<([u8; 8], usize)>,
    cache: &mut MphfLookupCache,
    allow_relaxed: bool,
    allow_tail_minimizer: bool,
) -> Option<FastKmerMatch> {
    if min_candidates.is_empty() {
        return None;
    }
    let needs_slow_path = min_candidates.iter().any(|(min_bytes, _)| {
        *min_bytes == [0u8; 8] || minimizer_bucket_has_overcrowded_marker(index, cache, min_bytes)
    });
    if needs_slow_path {
        let mut rev_buf = Vec::new();
        return match_kmer_at_pos(
            index,
            kmer,
            allow_forward,
            allow_rev,
            diff,
            &mut rev_buf,
            false,
            false,
            false,
            false,
            true,
            allow_relaxed,
        )
        .map(
            |(unitig_id, start, used_revcomp, block_idx, matched_relaxed)| FastKmerMatch {
                unitig_id,
                start,
                block_idx,
                used_revcomp,
                forward_strand: !used_revcomp,
                matched_relaxed,
            },
        );
    }
    let matched = match_kmer_at_pos_fast_uncached_precomputed(
        index,
        kmer,
        read_fwd,
        read_rev,
        allow_forward,
        allow_rev,
        diff,
        min_candidates,
        None,
        cache,
        allow_relaxed,
    );
    if matched.is_some() || !allow_tail_minimizer {
        return matched;
    }
    if let Some((tail_bytes, tail_pos)) = tail_candidate
        && !min_candidates
            .iter()
            .any(|&(min_bytes, min_pos)| min_bytes == tail_bytes && min_pos == tail_pos)
    {
        if tail_bytes == [0u8; 8]
            || minimizer_bucket_has_overcrowded_marker(index, cache, &tail_bytes)
        {
            let mut rev_buf = Vec::new();
            return match_kmer_at_pos(
                index,
                kmer,
                allow_forward,
                allow_rev,
                diff,
                &mut rev_buf,
                false,
                false,
                false,
                false,
                true,
                allow_relaxed,
            )
            .map(
                |(unitig_id, start, used_revcomp, block_idx, matched_relaxed)| FastKmerMatch {
                    unitig_id,
                    start,
                    block_idx,
                    used_revcomp,
                    forward_strand: !used_revcomp,
                    matched_relaxed,
                },
            );
        }
        let tail_candidate = [(tail_bytes, tail_pos)];
        return match_kmer_at_pos_fast_uncached_precomputed(
            index,
            kmer,
            read_fwd,
            read_rev,
            allow_forward,
            allow_rev,
            diff,
            &tail_candidate,
            None,
            cache,
            allow_relaxed,
        );
    }
    matched
}

#[inline]
fn fast_match_cache_key_from_codes(fwd: u64, rev: u64) -> u64 {
    fwd ^ rev.rotate_left(17)
}

#[allow(clippy::too_many_arguments)]
fn match_kmer_at_pos_fast_with_codes(
    index: &BifrostIndex,
    kmer: &[u8],
    read_fwd: u64,
    read_rev: u64,
    allow_forward: bool,
    allow_rev: bool,
    diff: usize,
    min_candidates: &mut Vec<([u8; 8], usize)>,
    cache: &mut MphfLookupCache,
    allow_relaxed: bool,
    allow_tail_minimizer: bool,
    minimizer_cache_disabled: bool,
    match_cache_disabled: bool,
) -> Option<FastKmerMatch> {
    let cache_key = fast_match_cache_key_from_codes(read_fwd, read_rev);
    let flags = ((allow_forward as u8) << 3)
        | ((allow_rev as u8) << 2)
        | ((allow_relaxed as u8) << 1)
        | (allow_tail_minimizer as u8);
    if !match_cache_disabled && let Some(cached) = read_fast_match_cache(cache_key, flags) {
        return cached;
    }
    let matched = match_kmer_at_pos_fast_uncached_with_codes(
        index,
        kmer,
        read_fwd,
        read_rev,
        allow_forward,
        allow_rev,
        diff,
        min_candidates,
        cache,
        allow_relaxed,
        allow_tail_minimizer,
        minimizer_cache_disabled,
    );
    if !match_cache_disabled {
        write_fast_match_cache(cache_key, flags, matched);
    }
    matched
}

pub(crate) fn ec_for_read_bifrost_fast(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut min_candidates: Vec<([u8; 8], usize)> = Vec::with_capacity(8);
    let mut jump_candidates: Vec<([u8; 8], usize)> = Vec::with_capacity(8);
    let mut cache = MphfLookupCache::default();
    let allow_forward_base = strand != Strand::Reverse;
    let allow_rev_base = strand != Strand::Forward;
    let mut has_hit = false;
    let diff = index.k.saturating_sub(index.g);
    let mut best_match: Option<MatchInfo> = None;
    let mut first_hit: Option<Hit> = None;
    let mut hits: Vec<Hit> = Vec::new();
    let dlist_dummy = index.onlist.as_deref().map(onlist_cardinality);
    let mut any_dlist_kmer = false;
    let mut dlist_early_return = false;
    let mut append_dummy_hit = false;
    let mut online_intersection: Vec<u32> = Vec::new();
    let mut has_online_intersection = false;
    let mut backoff_until: Option<usize> = None;
    let mut fallback_rev_buf: Vec<u8> = Vec::new();
    let minimizer_cache_disabled = env_flag(
        &DISABLE_MINIMIZER_CACHE,
        "KALLISTORS_DISABLE_MINIMIZER_CACHE",
    );
    let match_cache_disabled = env_flag(&DISABLE_MATCH_CACHE, "KALLISTORS_DISABLE_MATCH_CACHE");
    let mut pos = 0usize;
    let last_pos = seq.len() - index.k;
    while pos <= last_pos {
        if let Some(state) = dbg.as_deref_mut() {
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
        let kmer_codes = encode_kmer_pair(kmer);
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
        if dlist_dummy.is_some() && kmer_in_dlist {
            append_dummy_hit = true;
        }
        if let Some(state) = dbg.as_deref_mut() {
            if !minimizer_candidates_cached_into_with_code(
                kmer,
                index.g,
                read_fwd,
                &mut min_candidates,
                minimizer_cache_disabled,
            ) {
                pos += 1;
                continue;
            }
            for (min_bytes, min_pos) in &min_candidates {
                let bucket = minimizer_bucket_idx(index, &mut cache, min_bytes);
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
        let matched = if let Some(matched) = match_kmer_at_pos_fast_uncached_with_codes(
            index,
            kmer,
            read_fwd,
            read_rev,
            allow_forward_base,
            allow_rev_base,
            diff,
            &mut min_candidates,
            &mut cache,
            true,
            !in_backoff,
            minimizer_cache_disabled,
        ) {
            matched
        } else if !options.kallisto_strict
            && !options.kallisto_enum
            && in_backoff
            && min_candidates.len() == 1
            && let Some((uid, start, used_revcomp, block_idx, matched_relaxed)) = match_kmer_at_pos(
                index,
                kmer,
                allow_forward_base,
                allow_rev_base,
                diff,
                &mut fallback_rev_buf,
                options.kallisto_direct_kmer,
                options.kallisto_enum,
                options.kallisto_strict,
                options.skip_overcrowded_minimizer,
                true,
                true,
            )
        {
            FastKmerMatch {
                unitig_id: uid,
                start,
                block_idx,
                used_revcomp,
                forward_strand: !used_revcomp,
                matched_relaxed,
            }
        } else {
            if let Some(state) = dbg.as_deref_mut()
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
        if let Some(state) = dbg.as_deref_mut() {
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
            if let Some(state) = dbg.as_deref_mut()
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
        if accept_for_stream && !ec.is_empty() {
            if has_online_intersection {
                next.clear();
                intersect_sorted(&online_intersection, ec, &mut next);
                if next.is_empty() {
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
                std::mem::swap(&mut online_intersection, &mut next);
            } else {
                online_intersection.extend_from_slice(ec);
                has_online_intersection = true;
            }
        }
        if best_match.map(|m| kmer_pos < m.read_pos).unwrap_or(true) {
            best_match = Some(MatchInfo {
                unitig_id: uid,
                unitig_pos: start,
                read_pos: kmer_pos,
                used_revcomp,
            });
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
        if let Some(dist) = jump_distance_for_match(index, uid, block_idx, start, forward_strand) {
            let mut next_pos = pos + dist;
            if next_pos > last_pos {
                next_pos = last_pos;
            }
            if next_pos > pos {
                let kmer_next = &seq[next_pos..next_pos + index.k];
                let kmer_next_codes = encode_kmer_pair(kmer_next);
                let mut next_hit = kmer_next_codes.and_then(|(next_fwd, next_rev)| {
                    match_kmer_at_pos_fast_with_codes(
                        index,
                        kmer_next,
                        next_fwd,
                        next_rev,
                        allow_forward_base,
                        allow_rev_base,
                        diff,
                        &mut jump_candidates,
                        &mut cache,
                        true,
                        !in_backoff,
                        minimizer_cache_disabled,
                        match_cache_disabled,
                    )
                });
                if next_hit.is_none()
                    && in_backoff
                    && let Some((unitig_id, start, used_revcomp, block_idx, matched_relaxed)) =
                        match_kmer_at_pos(
                            index,
                            kmer_next,
                            allow_forward_base,
                            allow_rev_base,
                            diff,
                            &mut fallback_rev_buf,
                            options.kallisto_direct_kmer,
                            options.kallisto_enum,
                            options.kallisto_strict,
                            options.skip_overcrowded_minimizer,
                            options.kallisto_bifrost_find,
                            true,
                        )
                {
                    next_hit = Some(FastKmerMatch {
                        unitig_id,
                        start,
                        block_idx,
                        used_revcomp,
                        forward_strand: !used_revcomp,
                        matched_relaxed,
                    });
                }
                let current_ec = ec_slice(index, uid, block_idx);
                let mut use_backoff = false;
                if let Some(next_hit) = next_hit {
                    let uid2 = next_hit.unitig_id;
                    let block_idx2 = next_hit.block_idx;
                    let next_relaxed = next_hit.matched_relaxed;
                    let next_ec = ec_slice(index, uid2, block_idx2);
                    if uid2 == uid && next_ec == current_ec {
                        if let Some(state) = dbg.as_deref_mut() {
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
                            jump_to = Some(next_pos);
                        }
                    } else if dist > 4 {
                        let middle_pos = (pos + next_pos) / 2;
                        if middle_pos <= last_pos {
                            let kmer_mid = &seq[middle_pos..middle_pos + index.k];
                            let mut mid_hit =
                                encode_kmer_pair(kmer_mid).and_then(|(mid_fwd, mid_rev)| {
                                    match_kmer_at_pos_fast_with_codes(
                                        index,
                                        kmer_mid,
                                        mid_fwd,
                                        mid_rev,
                                        allow_forward_base,
                                        allow_rev_base,
                                        diff,
                                        &mut jump_candidates,
                                        &mut cache,
                                        true,
                                        !in_backoff,
                                        minimizer_cache_disabled,
                                        match_cache_disabled,
                                    )
                                });
                            if mid_hit.is_none()
                                && in_backoff
                                && let Some((
                                    unitig_id,
                                    start,
                                    used_revcomp,
                                    block_idx,
                                    matched_relaxed,
                                )) = match_kmer_at_pos(
                                    index,
                                    kmer_mid,
                                    allow_forward_base,
                                    allow_rev_base,
                                    diff,
                                    &mut fallback_rev_buf,
                                    options.kallisto_direct_kmer,
                                    options.kallisto_enum,
                                    options.kallisto_strict,
                                    options.skip_overcrowded_minimizer,
                                    options.kallisto_bifrost_find,
                                    true,
                                )
                            {
                                mid_hit = Some(FastKmerMatch {
                                    unitig_id,
                                    start,
                                    block_idx,
                                    used_revcomp,
                                    forward_strand: !used_revcomp,
                                    matched_relaxed,
                                });
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
                                    if let Some(state) = dbg.as_deref_mut() {
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
                                    if !mid_ec.is_empty() {
                                        if has_online_intersection {
                                            next.clear();
                                            intersect_sorted(
                                                &online_intersection,
                                                mid_ec,
                                                &mut next,
                                            );
                                            if next.is_empty() {
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
                                            std::mem::swap(&mut online_intersection, &mut next);
                                        } else {
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
                                    } else {
                                        jump_to = Some(next_pos);
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
                        append_dummy_hit = dlist_dummy.is_some();
                    }
                    if let Some(state) = dbg.as_deref_mut() {
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
                        jump_to = Some(next_pos);
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
                                    match_kmer_at_pos_fast_with_codes(
                                        index,
                                        scan_kmer,
                                        scan_fwd,
                                        scan_rev,
                                        allow_forward_base,
                                        allow_rev_base,
                                        diff,
                                        &mut jump_candidates,
                                        &mut cache,
                                        options.investigation.candidate_reuse,
                                        !in_backoff,
                                        minimizer_cache_disabled,
                                        match_cache_disabled,
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
                                if let Some(state) = dbg.as_deref_mut() {
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
        let _ = in_backoff;
        pos += 1;
    }

    if !dlist_early_return && any_dlist_kmer && has_hit {
        append_dummy_hit = dlist_dummy.is_some();
    }
    if !has_hit || hits.is_empty() {
        return None;
    }

    let mut minpos = usize::MAX;
    let mut maxpos = 0usize;
    let mut last_unitig = None;
    let mut last_block_idx = None;
    let mut found_nonempty = false;
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
        next.clear();
        intersect_sorted(&current, ec, &mut next);
        std::mem::swap(&mut current, &mut next);
        if current.is_empty() {
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
        last_unitig = Some(hit.unitig_id);
        last_block_idx = Some(hit.block_idx);
    }

    if append_dummy_hit && let Some(dummy) = dlist_dummy {
        let ec = [dummy];
        if current.is_empty() {
            current.extend_from_slice(&ec);
        } else {
            next.clear();
            intersect_sorted(&current, &ec, &mut next);
            std::mem::swap(&mut current, &mut next);
            if current.is_empty() {
                return Some(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: true,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                });
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
    Some(ReadEc {
        ec: current,
        best_match,
        first_hit,
        hits,
        had_offlist: false,
        shade_union: Vec::new(),
        hard_reject_pair: false,
    })
}
