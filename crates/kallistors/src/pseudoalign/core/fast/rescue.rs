//! Exact encoded Bifrost rescue for ordinary uppercase 31/23 probes.
use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn match_bifrost_rescue<const CACHED: bool>(
    index: &BifrostIndex,
    kmer: &[u8],
    forward: u64,
    reverse: u64,
    allow_forward: bool,
    allow_reverse: bool,
    diff: usize,
    cache: &mut MphfLookupCache,
    relaxed: bool,
) -> Option<FastKmerMatch> {
    let legacy = |cache: &mut MphfLookupCache| {
        match_kmer_at_pos(
            index,
            kmer,
            allow_forward,
            allow_reverse,
            diff,
            &mut cache.slow_rev_buf,
            false,
            false,
            false,
            false,
            true,
            relaxed,
        )
        .map(
            |(unitig_id, start, used_revcomp, block_idx, matched_relaxed)| FastKmerMatch {
                unitig_id,
                start,
                used_revcomp,
                block_idx,
                matched_relaxed,
                forward_strand: !used_revcomp,
            },
        )
    };
    if index.k != 31 || index.g != 23 || !all_ascii_uppercase(kmer) {
        return legacy(cache);
    }
    let reversed = reverse < forward;
    let (_, candidates) = super::super::super::minimizers::minhash_primary_for_encoded(
        kmer.len(),
        index.g,
        forward.min(reverse),
        forward.max(reverse),
    )?;
    for (key, position) in candidates {
        let Some((bucket, fingerprint)) = minimizer_bucket_with_fingerprint(index, cache, &key)
        else {
            continue;
        };
        if !matching::bucket_may_match(index, cache, fingerprint, &key, forward, reverse) {
            continue;
        }
        let (positions, overcrowded) = index.minz_positions.get_with_overcrowded(bucket);
        // Do not traverse a difficult bucket twice before delegating its ordered
        // next-minimum search. Packed descriptors expose this flag directly.
        if overcrowded {
            return legacy(cache);
        }
        if CACHED && positions.is_empty() {
            continue;
        }
        if CACHED && let Some(entries) = index.large_candidates_for(positions) {
            let (min_forward, min_reverse) = if reversed {
                (diff.saturating_sub(position), position)
            } else {
                (position, diff.saturating_sub(position))
            };
            if let Some(found) = matching::match_cached_candidates(
                index,
                entries,
                positions,
                diff,
                min_forward,
                min_reverse,
                forward,
                reverse,
                allow_forward,
                allow_reverse,
                relaxed,
            )
            .ok()?
            {
                return Some(found);
            }
            continue;
        }
        for &value in positions {
            if (value >> 32) as u32 == u32::MAX {
                // Marker handling includes ordered next-minimum rescues and raw
                // special-node lookup. Keep the established oracle for it.
                return legacy(cache);
            }
            let (uid, relative, km) = crate::index::bifrost::decode_pos_id(value);
            let uid = uid as usize;
            if km {
                if index
                    .encoded_unitigs
                    .km_bases
                    .encode_kmer_at(uid, 0, index.k)
                    != Some(forward.min(reverse))
                    || (position != relative as usize && position + relative as usize != diff)
                    || (reversed && !allow_reverse)
                    || (!reversed && !allow_forward)
                {
                    continue;
                }
                let unitig_id = index.unitig_count() + uid;
                let block_idx = block_index_for_position_fast(index, unitig_id, 0)?;
                return Some(FastKmerMatch {
                    unitig_id,
                    start: 0,
                    block_idx,
                    used_revcomp: reversed,
                    forward_strand: !reversed,
                    matched_relaxed: false,
                });
            }
            let (min_forward, min_reverse) = if reversed {
                (diff.saturating_sub(position), position)
            } else {
                (position, diff.saturating_sub(position))
            };
            if let Some((start, used_revcomp, matched_relaxed, forward_strand)) =
                match_unitig_candidate_encoded_typed::<CACHED>(
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
    None
}

// Jump and middle probes call this only after an exact primary-only miss.
// Preserve ordinary six-anchor-plus-tail semantics; Bifrost's seven-anchor
// rescues have different marker rules and remain a separate operation.
#[allow(clippy::too_many_arguments)]
pub(super) fn match_ordinary_tail_after_primary_miss<const CACHED: bool>(
    index: &BifrostIndex,
    kmer: &[u8],
    codes: (u64, u64),
    allow_forward: bool,
    allow_reverse: bool,
    diff: usize,
    candidates: &mut Vec<([u8; 8], usize)>,
    cache: &mut MphfLookupCache,
    relaxed: bool,
    skip_overcrowded: bool,
    minimizer_cache_disabled: bool,
) -> Option<FastKmerMatch> {
    let legacy = |cache: &mut MphfLookupCache| {
        match_kmer_at_pos(
            index,
            kmer,
            allow_forward,
            allow_reverse,
            diff,
            &mut cache.slow_rev_buf,
            false,
            false,
            false,
            skip_overcrowded,
            false,
            relaxed,
        )
        .map(
            |(unitig_id, start, used_revcomp, block_idx, matched_relaxed)| FastKmerMatch {
                unitig_id,
                start,
                used_revcomp,
                block_idx,
                matched_relaxed,
                forward_strand: !used_revcomp,
            },
        )
    };
    if index.k != 31 || index.g != 23 || !all_ascii_uppercase(kmer) {
        return legacy(cache);
    }
    let (forward, reverse) = codes;
    // A negative match-cache hit may not have filled scanner scratch. Recover
    // the exact primary candidates, normally from the bounded candidate cache.
    let valid = if CACHED {
        minimizer_candidates_single_cell(cache, kmer, index.g, forward, reverse, candidates)
    } else {
        minimizer_candidates_cached_into_with_code(
            cache,
            kmer,
            index.g,
            forward,
            reverse,
            candidates,
            minimizer_cache_disabled,
        )
    };
    if !valid {
        return legacy(cache);
    }
    if candidates.iter().any(|(key, _)| {
        *key == [0; 8]
            || minimizer_bucket_idx(index, cache, key)
                .is_some_and(|bucket| index.minz_positions.get_with_special(bucket).1)
    }) {
        // Even uncommon zero-count special markers retain the complete byte
        // oracle; a primary fast miss alone does not prove their legacy result.
        return legacy(cache);
    }
    let (key, position) = minimizer_tail_for_encoded(kmer.len(), index.g, forward, reverse)?;
    if candidates.contains(&(key, position)) {
        return None;
    }
    let (bucket, fingerprint) = minimizer_bucket_with_fingerprint(index, cache, &key)?;
    if !matching::bucket_may_match(index, cache, fingerprint, &key, forward, reverse) {
        return None;
    }
    let (positions, special) = index.minz_positions.get_with_special(bucket);
    if special {
        return legacy(cache);
    }
    if CACHED && positions.is_empty() {
        return None;
    }
    matching::match_kmer_at_pos_fast_uncached_precomputed::<CACHED>(
        index,
        kmer,
        forward,
        reverse,
        allow_forward,
        allow_reverse,
        diff,
        &[(key, position)],
        Some(positions),
        None,
        cache,
        relaxed,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoded_bifrost_rescue_matches_independent_legacy_oracle() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("reference.fa");
        let path = dir.path().join("reference.idx");
        let mut state = 79u64;
        let sequence: Vec<u8> = (0..700)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect();
        let repeat = [b'A'; 100];
        std::fs::write(
            &fasta,
            format!(
                ">a\n{}\n>b\n{}{}\n>repeat\n{}\n",
                String::from_utf8_lossy(&sequence),
                String::from_utf8_lossy(&sequence[..250]),
                String::from_utf8_lossy(&sequence[400..]),
                String::from_utf8_lossy(&repeat)
            ),
        )
        .unwrap();
        crate::index::build_index(&path, &[fasta], Default::default()).unwrap();
        let mut queries = Vec::new();
        for source in [&sequence[..], &repeat[..]] {
            for window in source.windows(31).step_by(29) {
                queries.push(window.to_vec());
                queries.push(crate::util::reverse_complement(window));
                let mut altered = window.to_vec();
                altered[13] = if altered[13] == b'A' { b'C' } else { b'A' };
                queries.push(altered);
                queries.push(window.to_ascii_lowercase());
            }
        }
        for index in [
            crate::pseudoalign::build_bifrost_index_threaded(&path, 2).unwrap(),
            crate::pseudoalign::build_bifrost_index_positional_packed_threaded(&path, 2, false)
                .unwrap(),
        ] {
            super::super::super::runtime::ensure_cache_context(&index);
            let mut cache = MphfLookupCache::default();
            let project = |value: Option<FastKmerMatch>| {
                value.map(|v| {
                    (
                        v.unitig_id,
                        v.start,
                        v.used_revcomp,
                        v.block_idx,
                        v.matched_relaxed,
                    )
                })
            };
            for query in &queries {
                let (forward, reverse) = crate::dna::encode_kmer_pair(query).unwrap();
                for skip in [false, true] {
                    for bits in 0..8 {
                        let allow_forward = bits & 1 != 0;
                        let allow_reverse = bits & 2 != 0;
                        let relaxed = bits & 4 != 0;
                        let expected = match_kmer_at_pos(
                            &index,
                            query,
                            allow_forward,
                            allow_reverse,
                            8,
                            &mut Vec::new(),
                            false,
                            false,
                            false,
                            skip,
                            true,
                            relaxed,
                        );
                        for _ in 0..2 {
                            let expected_ordinary = match_kmer_at_pos(
                                &index,
                                query,
                                allow_forward,
                                allow_reverse,
                                8,
                                &mut Vec::new(),
                                false,
                                false,
                                false,
                                skip,
                                false,
                                relaxed,
                            );
                            let mut candidates = Vec::new();
                            let primary = matching::match_kmer_at_pos_fast_with_codes::<false>(
                                &index,
                                query,
                                forward,
                                reverse,
                                allow_forward,
                                allow_reverse,
                                8,
                                &mut candidates,
                                &mut cache,
                                relaxed,
                                false,
                                false,
                                false,
                                true,
                                true,
                            );
                            if primary.is_none() || !all_ascii_uppercase(query) {
                                let actual = match_ordinary_tail_after_primary_miss::<false>(
                                    &index,
                                    query,
                                    (forward, reverse),
                                    allow_forward,
                                    allow_reverse,
                                    8,
                                    &mut candidates,
                                    &mut cache,
                                    relaxed,
                                    skip,
                                    false,
                                );
                                assert_eq!(
                                    project(actual),
                                    expected_ordinary,
                                    "ordinary skip={skip} flags={bits}"
                                );
                            }
                            let actual = match_bifrost_rescue::<false>(
                                &index,
                                query,
                                forward,
                                reverse,
                                allow_forward,
                                allow_reverse,
                                8,
                                &mut cache,
                                relaxed,
                            );
                            assert_eq!(project(actual), expected, "skip={skip} flags={bits}");
                        }
                    }
                }
            }
        }
    }
}
