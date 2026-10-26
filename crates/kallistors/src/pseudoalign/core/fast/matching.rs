//! Fast k-mer confirmation, matching scratch and hit insertion.

use super::*;

#[inline]
pub(super) fn all_ascii_uppercase(sequence: &[u8]) -> bool {
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
pub(super) fn bucket_may_match(
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
        || fingerprint == super::super::super::fingerprints::DISABLED
        || fingerprint
            == super::super::super::fingerprints::word_fingerprint(u64::from_le_bytes(*key))
        || cache.query_fingerprint_matches(forward, reverse, fingerprint)
}

// Encoded fast path for the common unshaded, non-union pseudoalignment mode.
#[allow(clippy::too_many_arguments)]
pub(super) fn match_kmer_at_pos_fast_uncached_precomputed<const CACHED: bool>(
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
        if CACHED && let Some(entries) = index.large_candidates_for(positions) {
            if let Some(found) = match_cached_candidates(
                index,
                entries,
                positions,
                diff,
                min_pos_fwd,
                min_pos_rev,
                read_fwd,
                read_rev,
                allow_forward,
                allow_rev,
                allow_relaxed,
            )
            .ok()?
            {
                return Some(found);
            }
            continue;
        }
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
                match_unitig_candidate_encoded_typed::<CACHED>(
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
                match_unitig_candidate_encoded_typed::<CACHED>(
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
pub(super) fn match_kmer_at_pos_fast_uncached_with_codes<const CACHED: bool>(
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
    let valid = if CACHED {
        minimizer_candidates_single_cell(cache, kmer, index.g, read_fwd, read_rev, min_candidates)
    } else {
        minimizer_candidates_cached_into_with_code(
            cache,
            kmer,
            index.g,
            read_fwd,
            read_rev,
            min_candidates,
            minimizer_cache_disabled,
        )
    };
    if !valid {
        return None;
    }
    match_kmer_at_pos_fast_uncached_with_candidates::<CACHED>(
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
pub(super) fn match_kmer_at_pos_fast_uncached_with_candidates<const CACHED: bool>(
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
        return super::rescue::match_bifrost_rescue::<CACHED>(
            index,
            kmer,
            read_fwd,
            read_rev,
            allow_forward,
            allow_rev,
            diff,
            cache,
            allow_relaxed,
        );
    }
    // first_positions exists only for one already resolved primary. Empty
    // positions prove the full verifier cannot match; retain the ordered tail
    // and special/overcrowded handling below and above this narrow shortcut.
    let matched = if CACHED && first_positions.is_some_and(|(positions, _)| positions.is_empty()) {
        None
    } else {
        match_kmer_at_pos_fast_uncached_precomputed::<CACHED>(
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
        )
    };
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
        let (tail_positions, tail_overcrowded) = if tail_bytes == [0u8; 8] {
            (&[][..], false)
        } else {
            minimizer_bucket_with_fingerprint(index, cache, &tail_bytes).map_or(
                (&[][..], false),
                |(bucket, fingerprint)| {
                    if bucket_may_match(index, cache, fingerprint, &tail_bytes, read_fwd, read_rev)
                    {
                        index.minz_positions.get_with_overcrowded(bucket)
                    } else {
                        (&[][..], false)
                    }
                },
            )
        };
        if tail_bytes == [0u8; 8] || tail_overcrowded {
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
        if CACHED && tail_positions.is_empty() {
            return None;
        }
        let tail_candidate = [(tail_bytes, tail_pos)];
        return match_kmer_at_pos_fast_uncached_precomputed::<CACHED>(
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
pub(super) fn match_kmer_at_pos_fast_with_codes<const CACHED: bool>(
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
                if CACHED {
                    minimizer_candidates_single_cell(
                        cache,
                        kmer,
                        index.g,
                        read_fwd,
                        read_rev,
                        min_candidates,
                    );
                } else {
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
            }
            return None;
        }
    }
    let matched = match_kmer_at_pos_fast_uncached_with_codes::<CACHED>(
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

// Cache only exact negative proofs. Potential hits retain the established
// verifier, annotation lookup and posting order, including relaxed metadata.
#[allow(clippy::too_many_arguments)]
#[inline]
pub(super) fn match_cached_candidates(
    index: &BifrostIndex,
    bucket: &crate::pseudoalign::bifrost::CachedBucket,
    positions: &[u64],
    diff: usize,
    min_forward: usize,
    min_reverse: usize,
    forward: u64,
    reverse: u64,
    allow_forward: bool,
    allow_reverse: bool,
    relaxed: bool,
) -> std::result::Result<Option<FastKmerMatch>, ()> {
    if diff == 8 && min_forward <= diff && min_reverse <= diff {
        if let Some(candidates) = bucket.matching_mask_candidates(
            positions,
            forward,
            reverse,
            allow_forward,
            allow_reverse,
        ) {
            return match_candidate_iterator(
                index,
                candidates,
                diff,
                min_forward,
                min_reverse,
                forward,
                reverse,
                allow_forward,
                allow_reverse,
                relaxed,
            );
        }
        match_candidate_iterator(
            index,
            bucket.matching_candidates_lists(
                positions,
                forward,
                reverse,
                allow_forward,
                allow_reverse,
            ),
            diff,
            min_forward,
            min_reverse,
            forward,
            reverse,
            allow_forward,
            allow_reverse,
            relaxed,
        )
    } else {
        match_candidate_iterator(
            index,
            positions.iter().copied(),
            diff,
            min_forward,
            min_reverse,
            forward,
            reverse,
            allow_forward,
            allow_reverse,
            relaxed,
        )
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn match_candidate_iterator(
    index: &BifrostIndex,
    candidates: impl Iterator<Item = u64>,
    diff: usize,
    min_forward: usize,
    min_reverse: usize,
    forward: u64,
    reverse: u64,
    allow_forward: bool,
    allow_reverse: bool,
    relaxed: bool,
) -> std::result::Result<Option<FastKmerMatch>, ()> {
    for position in candidates {
        // Directory tags select an ordered conservative superset. Exact
        // verification remains authoritative for both tagged selection and
        // the general-offset fallback, including annotation-failure stopping.
        let (uid, relative, _) = crate::index::bifrost::decode_pos_id(position);
        let uid = uid as usize;
        if let Some((start, used_revcomp, matched_relaxed, forward_strand)) =
            match_unitig_candidate_encoded_typed::<true>(
                index,
                uid,
                relative as usize,
                diff,
                min_forward,
                min_reverse,
                forward,
                reverse,
                allow_forward,
                allow_reverse,
                relaxed,
            )
        {
            let block_idx = block_index_for_position_fast(index, uid, start).ok_or(())?;
            return Ok(Some(FastKmerMatch {
                unitig_id: uid,
                start,
                block_idx,
                used_revcomp,
                matched_relaxed,
                forward_strand,
            }));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod cached_candidate_tests {
    use super::*;

    #[test]
    fn missing_ec_block_stops_query_instead_of_reporting_an_empty_bucket() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("reference.fa");
        let path = directory.path().join("reference.idx");
        std::fs::write(&fasta, ">short\nCGGGTGTGGTGGCGGGCACCTGTAGTCCCAG\n").unwrap();
        crate::index::build_index(&path, &[fasta], Default::default()).unwrap();
        let mut index = crate::pseudoalign::build_bifrost_index(&path).unwrap();
        let sequence = index.unitig_sequence(0).unwrap().into_owned();
        let (forward, reverse) = crate::dna::encode_kmer_pair(&sequence[..31]).unwrap();
        // The established traversal stops the entire query on this annotation
        // failure. A cached traversal must distinguish it from an empty bucket.
        index.flat_ec.owned_ec_blocks = (0..index.unitig_count()).map(|_| Box::default()).collect();
        super::super::super::runtime::ensure_cache_context(&index);
        let entries = [0u64];
        assert!(
            match_candidate_iterator(
                &index,
                entries.iter().copied(),
                8,
                0,
                8,
                forward,
                reverse,
                true,
                true,
                true
            )
            .is_err()
        );
    }
}
