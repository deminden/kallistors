use super::super::BorrowedReadEc;
use super::super::debug::{
    DroppedHitTrace, JumpDecisionTrace, MinimizerCandidateTrace, ReadDebugState,
    format_positions_sample,
};
use super::super::ec::{
    KmerPairCursor, encode_kmer_pair, encoded_kmer_at_or_after, intersect_sorted_in_place,
};
use super::super::minimizers::minimizer_tail_for_encoded;
use super::super::{BifrostIndex, Hit, MatchInfo, PseudoalignOptions, ReadEc, Strand};
use super::matcher::{
    jump_distance_for_match, match_kmer_at_pos, match_kmer_direct, match_unitig_candidate_encoded,
    special_unitig_for_kmer,
};
use super::runtime::{
    DISABLE_MATCH_CACHE, DISABLE_MINIMIZER_CACHE, FastKmerMatch, MphfLookupCache,
    block_index_for_position_fast, ec_slice, env_flag, minimizer_bucket_has_overcrowded_marker,
    minimizer_bucket_idx, minimizer_bucket_with_fingerprint,
    minimizer_candidates_cached_into_with_code, minimizer_hex, onlist_cardinality,
};

#[inline]
fn all_ascii_uppercase(sequence: &[u8]) -> bool {
    #[cfg(target_arch = "aarch64")]
    {
        use std::arch::aarch64::*;
        let (wide, rest) = sequence.as_chunks::<64>();
        for block in wide {
            // SAFETY: each of the four loads is contained in this 64-byte
            // array. Unsigned subtraction implements the exact A..Z range.
            let invalid = unsafe {
                let bad =
                    |bytes| vcgtq_u8(vsubq_u8(bytes, vdupq_n_u8(b'A')), vdupq_n_u8(b'Z' - b'A'));
                let a = bad(vld1q_u8(block.as_ptr()));
                let b = bad(vld1q_u8(block.as_ptr().add(16)));
                let c = bad(vld1q_u8(block.as_ptr().add(32)));
                let d = bad(vld1q_u8(block.as_ptr().add(48)));
                vmaxvq_u8(vorrq_u8(vorrq_u8(a, b), vorrq_u8(c, d)))
            };
            if invalid != 0 {
                return false;
            }
        }
        let (blocks, tail) = rest.as_chunks::<16>();
        for block in blocks {
            // SAFETY: one complete 16-byte array, baseline AArch64 NEON.
            let invalid = unsafe {
                vmaxvq_u8(vcgtq_u8(
                    vsubq_u8(vld1q_u8(block.as_ptr()), vdupq_n_u8(b'A')),
                    vdupq_n_u8(b'Z' - b'A'),
                ))
            };
            if invalid != 0 {
                return false;
            }
        }
        tail.iter().all(u8::is_ascii_uppercase)
    }
    #[cfg(not(target_arch = "aarch64"))]
    sequence.iter().all(u8::is_ascii_uppercase)
}

#[inline]
fn bucket_may_match(
    index: &BifrostIndex,
    cache: &mut MphfLookupCache,
    fingerprint: u16,
    key: &[u8; 8],
    forward: u64,
    reverse: u64,
) -> bool {
    if index.k != 31 || index.g != 23 {
        return true;
    }
    // Uninitialized, conflicting, KM and special buckets retain every probe.
    fingerprint == 0
        || fingerprint == super::super::fingerprints::DISABLED
        || fingerprint == super::super::fingerprints::word_fingerprint(u64::from_le_bytes(*key))
        || cache.query_fingerprint_matches(forward, reverse, fingerprint)
}

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
    first_positions: Option<&[u64]>,
    next_min_candidate: Option<([u8; 8], usize)>,
    cache: &mut MphfLookupCache,
    allow_relaxed: bool,
) -> Option<FastKmerMatch> {
    let use_revcomp = read_rev < read_fwd;
    let mut saw_overcrowded = false;
    for (candidate_idx, (min_bytes, min_pos)) in min_candidates.iter().copied().enumerate() {
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
        let positions = if candidate_idx == 0
            && let Some(positions) = first_positions
        {
            positions
        } else {
            let Some((bucket_idx, fingerprint)) =
                minimizer_bucket_with_fingerprint(index, cache, &min_bytes)
            else {
                continue;
            };
            if !bucket_may_match(index, cache, fingerprint, &min_bytes, read_fwd, read_rev) {
                continue;
            }
            index.minz_positions.get(bucket_idx)
        };
        for &pos_id in positions {
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
                if uid >= index.km_unitig_count() {
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
                let unitig_id = index.unitig_count() + uid;
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
            // The encoded matcher validates the node through its metadata lookup.
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
                if uid >= index.km_unitig_count() {
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
                let unitig_id = index.unitig_count() + uid;
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
            // The encoded matcher validates the node through its metadata lookup.
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
        cache,
        kmer,
        index.g,
        read_fwd,
        read_rev,
        min_candidates,
        minimizer_cache_disabled,
    ) {
        return None;
    }
    match_kmer_at_pos_fast_uncached_with_candidates(
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
    cache: &mut MphfLookupCache,
    allow_relaxed: bool,
    allow_tail_minimizer: bool,
) -> Option<FastKmerMatch> {
    if min_candidates.is_empty() {
        return None;
    }
    // Keep the first probe alive through matching instead of consulting the same
    // hash cache and decoding the same bucket a second time.
    let first_positions = if let [(min_bytes, _)] = min_candidates
        && *min_bytes != [0u8; 8]
    {
        Some(
            minimizer_bucket_with_fingerprint(index, cache, min_bytes).map_or(
                (&[][..], false),
                |(bucket, fingerprint)| {
                    if bucket_may_match(index, cache, fingerprint, min_bytes, read_fwd, read_rev) {
                        index.minz_positions.get_with_overcrowded(bucket)
                    } else {
                        (&[][..], false)
                    }
                },
            ),
        )
    } else {
        None
    };
    let needs_slow_path = if let Some((_, overcrowded)) = first_positions {
        overcrowded
    } else {
        min_candidates.iter().any(|(min_bytes, _)| {
            *min_bytes == [0u8; 8]
                || minimizer_bucket_has_overcrowded_marker(index, cache, min_bytes)
        })
    };
    if needs_slow_path {
        return match_kmer_at_pos(
            index,
            kmer,
            allow_forward,
            allow_rev,
            diff,
            &mut cache.slow_rev_buf,
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
        first_positions.map(|(positions, _)| positions),
        None,
        cache,
        allow_relaxed,
    );
    if matched.is_some() || !allow_tail_minimizer {
        return matched;
    }
    // The tail is consulted only after the primary candidates fail.
    if let Some((tail_bytes, tail_pos)) =
        minimizer_tail_for_encoded(kmer.len(), index.g, read_fwd, read_rev)
        && !min_candidates
            .iter()
            .any(|&(min_bytes, min_pos)| min_bytes == tail_bytes && min_pos == tail_pos)
    {
        let tail_positions = if tail_bytes == [0u8; 8] {
            &[][..]
        } else {
            minimizer_bucket_with_fingerprint(index, cache, &tail_bytes).map_or(
                &[][..],
                |(bucket, fingerprint)| {
                    if bucket_may_match(index, cache, fingerprint, &tail_bytes, read_fwd, read_rev)
                    {
                        index.minz_positions.get(bucket)
                    } else {
                        &[][..]
                    }
                },
            )
        };
        if tail_bytes == [0u8; 8]
            || tail_positions
                .iter()
                .any(|&pos_id| (pos_id >> 32) as u32 == u32::MAX && (pos_id & 0x8000_0000) != 0)
        {
            return match_kmer_at_pos(
                index,
                kmer,
                allow_forward,
                allow_rev,
                diff,
                &mut cache.slow_rev_buf,
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
            Some(tail_positions),
            None,
            cache,
            allow_relaxed,
        );
    }
    matched
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
    cache_misses: bool,
    negative_cache_enabled: bool,
) -> Option<FastKmerMatch> {
    // The full forward word is unique within the owning index/k context.
    let cache_key = read_fwd;
    let flags = ((allow_forward as u8) << 3)
        | ((allow_rev as u8) << 2)
        | ((allow_relaxed as u8) << 1)
        | (allow_tail_minimizer as u8);
    if !match_cache_disabled && let Some(cached) = cache.match_get(cache_key, flags) {
        if cache_misses || cached.is_some() {
            return cached;
        }
        if negative_cache_enabled {
            // Only ordinary backoff callers inspect the primary candidate list.
            if !allow_tail_minimizer {
                minimizer_candidates_cached_into_with_code(
                    cache,
                    kmer,
                    index.g,
                    read_fwd,
                    read_rev,
                    min_candidates,
                    minimizer_cache_disabled,
                );
            }
            return None;
        }
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
    if !match_cache_disabled && (cache_misses || matched.is_some() || negative_cache_enabled) {
        cache.match_put(cache_key, flags, matched);
    }
    matched
}

// The lazy final-hit fold shares this prefix. Keep an immutable graph slice
// until an actually different constraint requires writable intersection storage.
struct OnlineIntersection<'index, 'scratch> {
    initial: Option<&'index [u32]>,
    owned: &'scratch mut Vec<u32>,
}
impl std::ops::Deref for OnlineIntersection<'_, '_> {
    type Target = [u32];
    fn deref(&self) -> &[u32] {
        self.initial.unwrap_or(self.owned)
    }
}
impl<'index> OnlineIntersection<'index, '_> {
    fn initialize(&mut self, ec: &'index [u32]) {
        self.initial = Some(ec);
    }
    #[inline(always)]
    fn intersect(&mut self, ec: &'index [u32]) {
        if let Some(initial) = self.initial {
            if initial == ec {
                return;
            }
            // Nonoverlapping transcript ranges prove rejection before a copy.
            if initial.is_empty()
                || ec.is_empty()
                || initial[0] > ec[ec.len() - 1]
                || ec[0] > initial[initial.len() - 1]
            {
                self.initial = Some(&initial[..0]);
                return;
            }
            let (shorter, longer) = if initial.len() <= ec.len() {
                (initial, ec)
            } else {
                (ec, initial)
            };
            // Many adjacent transcript annotations add/remove alternatives. A
            // subset is already the exact intersection and can remain borrowed.
            if shorter.len() != longer.len() && sorted_subset(shorter, longer) {
                self.initial = Some(shorter);
                return;
            }
            self.owned.extend_from_slice(shorter);
            self.initial = None;
            intersect_sorted_in_place(self.owned, longer);
        } else {
            intersect_sorted_in_place(self.owned, ec);
        }
    }
}

fn sorted_subset(shorter: &[u32], longer: &[u32]) -> bool {
    if shorter.is_empty() {
        return true;
    }
    if longer.is_empty()
        || shorter[0] < longer[0]
        || shorter[shorter.len() - 1] > longer[longer.len() - 1]
    {
        return false;
    }
    let mut cursor = 0;
    for &value in shorter {
        if longer.len() / shorter.len() >= 16 {
            cursor += longer[cursor..].partition_point(|&candidate| candidate < value);
        } else {
            while cursor < longer.len() && longer[cursor] < value {
                cursor += 1;
            }
        }
        if longer.get(cursor) != Some(&value) {
            return false;
        }
        cursor += 1;
    }
    true
}

#[derive(Default)]
pub(super) struct FastScannerScratch {
    min_candidates: Vec<([u8; 8], usize)>,
    jump_candidates: Vec<([u8; 8], usize)>,
    online_intersection: Vec<u32>,
    fallback_rev_buf: Vec<u8>,
}

impl FastScannerScratch {
    fn clear(&mut self) {
        self.min_candidates.clear();
        self.jump_candidates.clear();
        self.online_intersection.clear();
        self.fallback_rev_buf.clear();
    }
}

pub(crate) fn ec_for_read_bifrost_fast(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    super::runtime::with_mphf_lookup_cache(|cache| {
        let mut scratch = std::mem::take(&mut cache.scan_scratch);
        scratch.clear();
        let result = if dbg.is_some() {
            ec_for_read_bifrost_fast_with_cache::<true, false>(
                index,
                seq,
                strand,
                dbg,
                options,
                cache,
                &mut scratch,
            )
        } else {
            ec_for_read_bifrost_fast_with_cache::<false, false>(
                index,
                seq,
                strand,
                dbg,
                options,
                cache,
                &mut scratch,
            )
        };
        cache.scan_scratch = scratch;
        result.map(BorrowedReadEc::into_owned)
    })
}

#[inline(always)]
pub(super) fn ec_for_read_bifrost_fast_counting_borrowed<'index>(
    index: &'index BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<BorrowedReadEc<'index>> {
    super::runtime::with_mphf_lookup_cache(|cache| {
        let mut scratch = std::mem::take(&mut cache.scan_scratch);
        scratch.clear();
        let result = ec_for_read_bifrost_fast_with_cache::<false, true>(
            index,
            seq,
            strand,
            None,
            options,
            cache,
            &mut scratch,
        );
        cache.scan_scratch = scratch;
        result
    })
}

fn ec_for_read_bifrost_fast_with_cache<'index, const TRACE: bool, const COUNT_ONLY: bool>(
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
    let mut counted = super::count::FinalHitFold::default();
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
    };
    let mut has_online_intersection = false;
    let mut last_online_block = None;
    let mut backoff_until: Option<usize> = None;
    let mut fallback_probe_pending = false;
    let fallback_rev_buf = &mut scratch.fallback_rev_buf;
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
        let matched = if let Some(matched) = match_kmer_at_pos_fast_with_codes(
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
            && let Some((uid, start, used_revcomp, block_idx, matched_relaxed)) = match_kmer_at_pos(
                index,
                kmer,
                allow_forward_base,
                allow_rev_base,
                diff,
                fallback_rev_buf,
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
        if !COUNT_ONLY && best_match.map(|m| kmer_pos < m.read_pos).unwrap_or(true) {
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
                counted =
                    super::count::FinalHitFold::from_online(online_intersection, last_online_block);
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
        if !COUNT_ONLY
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
                    match_kmer_at_pos_fast_with_codes(
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
                if next_hit.is_none()
                    && in_backoff
                    && let Some((unitig_id, start, used_revcomp, block_idx, matched_relaxed)) =
                        match_kmer_at_pos(
                            index,
                            kmer_next,
                            allow_forward_base,
                            allow_rev_base,
                            diff,
                            fallback_rev_buf,
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
                                match_kmer_at_pos_fast_with_codes(
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
                                    fallback_rev_buf,
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
                                    match_kmer_at_pos_fast_with_codes(
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
                            counted = super::count::FinalHitFold::from_online(
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

#[cfg(test)]
mod uppercase_tests {
    use super::all_ascii_uppercase;
    #[test]
    fn uppercase_vectors_match_scalar_for_all_bytes_and_record_tails() {
        for len in [
            0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 95, 100, 127, 128, 129, 257,
        ] {
            let mut input = vec![b'N'; len];
            assert!(all_ascii_uppercase(&input));
            for position in 0..len {
                for value in 0..=u8::MAX {
                    input[position] = value;
                    assert_eq!(
                        all_ascii_uppercase(&input),
                        input.iter().all(u8::is_ascii_uppercase),
                        "{len} {position} {value}"
                    );
                }
                input[position] = b'N';
            }
        }
    }
}

#[cfg(test)]
mod borrowed_prefix_tests {
    use super::*;
    #[test]
    fn borrowed_prefix_materializes_only_when_a_constraint_changes() {
        let original = [1, 1, 3, 7];
        let same = [1, 1, 3, 7];
        let mut scratch = Vec::new();
        let mut online = OnlineIntersection {
            initial: None,
            owned: &mut scratch,
        };
        online.initialize(&original);
        online.intersect(&same);
        assert!(online.initial.is_some());
        assert_eq!(&*online, &original);
        online.intersect(&[1, 3, 3, 8]);
        assert!(online.initial.is_none());
        assert_eq!(&*online, &[1, 3]);
        online.intersect(&[3]);
        assert_eq!(&*online, &[3]);
        online.intersect(&[]);
        assert!(online.is_empty());
        let mut multisets: Vec<Vec<u32>> = (0..81)
            .map(|mut code| {
                let mut values = Vec::new();
                for value in 0..4 {
                    for _ in 0..code % 3 {
                        values.push(value);
                    }
                    code /= 3;
                }
                values
            })
            .collect();
        multisets.extend([(0..256).collect(), vec![15, 15], vec![255], vec![254, 256]]);
        for left in &multisets {
            for right in &multisets {
                let mut expected = Vec::new();
                super::super::super::ec::intersect_sorted(left, right, &mut expected);
                let mut scratch = Vec::new();
                let mut intersection = OnlineIntersection {
                    initial: Some(left),
                    owned: &mut scratch,
                };
                intersection.intersect(right);
                assert_eq!(&*intersection, &expected);
                if expected == *left || expected == *right {
                    assert!(intersection.initial.is_some());
                }
            }
        }
    }
}
