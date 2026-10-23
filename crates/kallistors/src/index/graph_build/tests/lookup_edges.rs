//! Lookup edges regressions.

use super::*;
#[test]
fn bidirectional_edge_states_match_separate_orientation_oracle() {
    for parallel in [false, true] {
        for bias in [0, 1 << 28] {
            let offsets = (0u32..1537)
                .map(|i| {
                    bias + (i / 256) * 100_000 + (i % 256) * [1, 20, 300][(i / 256) as usize % 3]
                })
                .collect::<Vec<_>>();
            let ranks = RankDirectory::new(offsets.clone(), parallel);
            assert_eq!(ranks.delta_bits, if bias == 0 { 12 } else { 16 });
            for (i, &expected) in offsets.iter().enumerate() {
                assert_eq!(ranks.get(i), expected as usize);
                if i + 1 < offsets.len() {
                    assert_eq!(ranks.range(i), (expected as usize, offsets[i + 1] as usize));
                }
            }
        }
        let wide = [0, (1 << 31) - 1, 1 << 31, u32::MAX];
        let ranks = RankDirectory::new(wide.to_vec(), parallel);
        for (i, expected) in wide.into_iter().enumerate() {
            assert_eq!(ranks.get(i), expected as usize);
        }
    }

    for k in [3, 5, 7] {
        let mask = kmer_mask(k);
        for value in 0..1u64 << (2 * k) {
            let word = value << (64 - 2 * k);
            let reverse = reverse_complement_word(word, k, mask);
            let self_edges = (0..4)
                .filter(|&base| {
                    let next = shift_append_word(word, base, k, mask);
                    next == word || next == reverse
                })
                .count();
            assert!(
                self_edges <= 1,
                "multiple self successors k={k} word={word}"
            );
        }
    }
    for length in 5..=12 {
        let mut entries = (0..length - 1)
            .map(|value| (value as u64) << 2)
            .collect::<Vec<_>>();
        entries.push(1u64 << 36);
        let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
        for query in [0, 1, 4, 8, 44, 48, 1u64 << 36, 1u64 << 37] {
            assert_eq!(
                lookup.find(&entries, query),
                entries.binary_search(&query).ok()
            );
            assert_eq!(
                lookup.find_many(&entries, [query; 4]),
                [entries.binary_search(&query).ok(); 4]
            );
        }
    }

    let entries = [0, 4, 8, 1u64 << 16];
    let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
    for query in [0, 1, 4, 8, 12, 1u64 << 16, 1u64 << 17] {
        assert_eq!(
            lookup.find(&entries, query),
            entries.binary_search(&query).ok()
        );
    }

    let entries = [0, 4, 1u64 << 21];
    let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
    for query in [0, 1, 4, 8, 1u64 << 21, 1u64 << 22] {
        assert_eq!(
            lookup.find(&entries, query),
            entries.binary_search(&query).ok()
        );
    }

    // Equal guard bits still require complete-key equality, including padding.
    for entries in [vec![0, 1u64 << 32], vec![0, 4]] {
        let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
        for query in [0, 1, 4, 8, 1u64 << 32, 1u64 << 33] {
            assert_eq!(
                lookup.find(&entries, query),
                entries.binary_search(&query).ok()
            );
            assert_eq!(
                lookup.find_many(&entries, [query; 4]),
                [entries.binary_search(&query).ok(); 4]
            );
        }
    }

    for next in [0, (1 << 28) - 1] {
        let cached = CachedEdges {
            packed: (1 << EDGE_OUT_DEGREE_SHIFT)
                | (1 << (EDGE_STATE_BITS + EDGE_OUT_DEGREE_SHIFT))
                | (3 << EDGE_OUT_CODE_SHIFT)
                | (2 << (EDGE_STATE_BITS + EDGE_OUT_CODE_SHIFT))
                | (1 << (2 * EDGE_STATE_BITS)),
            next: [next, next],
        };
        let compact = compact_cached_edges(cached);
        let cache = EdgeStateCache {
            packed: EdgeStorage::Compact(vec![compact_pair_word(compact)]),
        };
        assert_eq!(cache.neighbor(0, true, 0), (next as usize, true));
        assert_eq!(cache.neighbor(0, false, 0), (next as usize, false));
    }
    let wide = EdgeStateCache {
        packed: EdgeStorage::Wide(vec![CachedEdges {
            packed: 1 << (2 * EDGE_STATE_BITS + 1),
            next: [u32::MAX, 1 << 28],
        }]),
    };
    assert_eq!(wide.neighbor(0, true, 0), (u32::MAX as usize, false));
    assert_eq!(wide.neighbor(0, false, 0), (1 << 28, true));

    assert_eq!(std::mem::size_of::<RankedEdges>(), 9);
    for count in [(1usize << 28) + 1, (1usize << 31) - 5] {
        for degrees in [[0u32, 4], [1, 1], [2, 3]] {
            for flags in 0..4u32 {
                let cached = CachedEdges {
                    packed: (degrees[0] << EDGE_OUT_DEGREE_SHIFT)
                        | (degrees[1] << (EDGE_STATE_BITS + EDGE_OUT_DEGREE_SHIFT))
                        | (3 << EDGE_OUT_CODE_SHIFT)
                        | (2 << (EDGE_STATE_BITS + EDGE_OUT_CODE_SHIFT))
                        | (1 << EDGE_HAS_NON_SELF_SHIFT)
                        | (flags << (2 * EDGE_STATE_BITS)),
                    next: [0, count as u32 - 1],
                };
                let ranked = RankedEdges::new(cached, count);
                for (orientation, &degree) in degrees.iter().enumerate() {
                    let edge = ranked.edge(orientation);
                    assert_eq!(RankedEdges::degree(edge, count), degree);
                    assert_eq!(
                        (ranked.codes.load(Ordering::Relaxed) >> (2 * orientation)) & 3,
                        3 - orientation as u8
                    );
                    if degree == 1 {
                        assert_eq!(
                            edge & 0x7fff_ffff,
                            if orientation == 0 {
                                0
                            } else {
                                count as u32 - 1
                            }
                        );
                        assert_eq!(edge >> 31, (flags >> orientation) & 1);
                    } else {
                        assert_eq!(edge & 0x7fff_ffff, count as u32 + degree);
                        assert_eq!(edge >> 31, u32::from(orientation == 0));
                    }
                }
            }
        }
    }

    for k in (3..=31).step_by(2) {
        let mask = kmer_mask(k);
        let mut entries = vec![0];
        let mut state = 42u64;
        for _ in 0..512 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let word = state & mask;
            entries.push(word.min(reverse_complement_word_slow(word, k)));
        }
        let center = entries[1];
        for code in 0..4 {
            for word in [
                shift_append_word(center, code, k, mask),
                shift_prepend_word(center, code, mask),
            ] {
                entries.push(word.min(reverse_complement_word_slow(word, k)));
            }
        }
        entries.sort_unstable();
        entries.dedup();
        let membership = entries.iter().copied().collect::<HashSet<_>>();
        let lookup = KmerLookup::with_threads(&entries, 1, k).unwrap();
        let symmetric_cache = EdgeStateCache::build(&entries, &lookup, k, 1).unwrap();
        let ranked_cache = EdgeStateCache {
            packed: EdgeStorage::Ranked(RankedEdges::build(&entries, &lookup, k, k % 4 == 1)),
        };
        let relative_cache = EdgeStateCache {
            packed: EdgeStorage::Relative(RelativeEdges::build(
                &entries,
                &lookup,
                k,
                RelativeEdges::prefixes(&entries, &lookup, k).unwrap(),
                k % 4 == 1,
            )),
        };
        let mut fallback = KmerLookup::new(&entries).unwrap();
        fallback.directory = KmerDirectory::Ranks(RankDirectory::new(
            (0..=(1usize << (64 - fallback.shift)))
                .map(|bucket| {
                    entries.partition_point(|&word| (word >> fallback.shift) < bucket as u64) as u32
                })
                .collect(),
            false,
        ));
        let packed_cache = EdgeStateCache {
            packed: EdgeStorage::Relative7(PackedRelativeEdges::build(
                &entries,
                &fallback,
                k,
                PackedRelativeEdges::prefixes(&entries, &lookup, k).unwrap(),
                true,
            )),
        };
        let dense_cache = EdgeStateCache {
            packed: EdgeStorage::Relative6(DenseRelativeEdges::build(
                &entries,
                &fallback,
                k,
                DenseRelativeEdges::prefixes(&entries, &lookup, k).unwrap(),
                true,
            )),
        };
        let adaptive_cache = EdgeStateCache {
            packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_plan(
                &entries,
                &fallback,
                k,
                adaptive_prefixes::PrefixPlan::new(&entries, &fallback, k, 1, 4),
                true,
            )),
        };
        let mut superset = entries.clone();
        for &word in &entries {
            superset.extend(edge_neighbor_words(word, k, mask));
        }
        superset.sort_unstable();
        superset.dedup();
        let filtered_caches = [&entries, &superset].map(|filter_keys| {
            let filter = negative_filter::NegativeFilter::build(filter_keys, k, 1).unwrap();
            let cache = EdgeStateCache {
                packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_filtered(
                    &entries,
                    &fallback,
                    k,
                    adaptive_prefixes::PrefixPlan::new(&entries, &fallback, k, 1, 4),
                    &filter,
                    true,
                )),
            };
            filter.retire();
            cache
        });
        for &word in &entries {
            let queries = [
                word,
                word | 1,
                word ^ (1u64 << (64 - 2 * k)),
                word ^ (1u64 << 63),
            ];
            let expected = queries.map(|query| entries.binary_search(&query).ok());
            for count in 0..=queries.len() {
                let mut partial = expected;
                partial[count..].fill(None);
                assert_eq!(lookup.find_many_count(&entries, queries, count), partial);
                assert_eq!(fallback.find_many_count(&entries, queries, count), partial);
            }

            assert_eq!(
                lookup.find_many(&entries, queries),
                expected,
                "inline k={k}"
            );
            assert_eq!(
                fallback.find_many(&entries, queries),
                expected,
                "rank k={k}"
            );
            for (query, expected) in queries.into_iter().zip(expected) {
                assert_eq!(lookup.find(&entries, query), expected);
                assert_eq!(fallback.find(&entries, query), expected);
            }
        }
        for &word in &entries {
            let oracle = |oriented| {
                let mut incoming = 0;
                let mut outgoing = 0;
                let mut out_code = 0;
                let mut non_self = false;
                let mut self_edges = 0;
                for code in 0..4 {
                    let next = shift_append_word(oriented, code, k, mask);
                    let next = next.min(reverse_complement_word_slow(next, k));
                    self_edges += usize::from(next == word);
                    if membership.contains(&next) {
                        outgoing += 1;
                        out_code = code as u32;
                        non_self |= next != word;
                    }
                    let prev = shift_prepend_word(oriented, code, mask);
                    let prev = prev.min(reverse_complement_word_slow(prev, k));
                    incoming += u32::from(membership.contains(&prev));
                }
                assert!(self_edges <= 1, "canonical self-edge bound k={k}");
                incoming
                    | (outgoing << EDGE_OUT_DEGREE_SHIFT)
                    | (out_code << EDGE_OUT_CODE_SHIFT)
                    | (u32::from(non_self) << EDGE_HAS_NON_SELF_SHIFT)
            };
            let cached = cache_edge_states(word, &entries, &lookup, k, mask);
            let idx = lookup.find(&entries, word).unwrap();
            let compact = compact_cached_edges(cache_edge_states(word, &entries, &lookup, k, mask));
            let mut compact_entries = vec![0u64; entries.len()];
            compact_entries[idx] = compact_pair_word(compact);
            let compact_cache = EdgeStateCache {
                packed: EdgeStorage::Compact(compact_entries),
            };
            let mut wide_entries = (0..entries.len())
                .map(|_| CachedEdges {
                    packed: 0,
                    next: [0; 2],
                })
                .collect::<Vec<_>>();
            wide_entries[idx] = cache_edge_states(word, &entries, &lookup, k, mask);
            let wide_cache = EdgeStateCache {
                packed: EdgeStorage::Wide(wide_entries),
            };
            for forward in [true, false] {
                let oriented = if forward {
                    word
                } else {
                    reverse_complement_word_slow(word, k)
                };
                assert_eq!(
                    compact_cache.state_bits(idx, forward, word),
                    wide_cache.state_bits(idx, forward, word)
                );
                let bits = wide_cache.state_bits(idx, forward, word);
                assert_eq!(ranked_cache.state_bits(idx, forward, word), bits);
                assert_eq!(relative_cache.state_bits(idx, forward, word), bits);
                assert_eq!(packed_cache.state_bits(idx, forward, word), bits);
                assert_eq!(dense_cache.state_bits(idx, forward, word), bits);
                assert_eq!(adaptive_cache.state_bits(idx, forward, word), bits);
                for filtered in &filtered_caches {
                    assert_eq!(filtered.state_bits(idx, forward, word), bits);
                    assert_eq!(
                        filtered.outgoing(idx, forward, oriented),
                        wide_cache.outgoing(idx, forward, oriented)
                    );
                }

                assert_eq!(
                    adaptive_cache.outgoing(idx, forward, oriented),
                    wide_cache.outgoing(idx, forward, oriented)
                );

                assert_eq!(
                    dense_cache.outgoing(idx, forward, oriented),
                    relative_cache.outgoing(idx, forward, oriented)
                );
                assert_eq!(
                    packed_cache.outgoing(idx, forward, oriented),
                    relative_cache.outgoing(idx, forward, oriented)
                );
                assert_eq!(relative_cache.is_internal(idx), bits & 0x3f == 9);
                assert_eq!(
                    relative_cache.incoming_degree(idx, forward),
                    (bits & 7) as u8
                );
                assert_eq!(
                    relative_cache.outgoing_continuation(idx, forward, oriented),
                    symmetric_cache.outgoing_continuation(idx, forward, oriented)
                );
                assert_eq!(ranked_cache.incoming_degree(idx, forward), (bits & 7) as u8);
                assert_eq!(ranked_cache.is_internal(idx), bits & 0x3f == 9);
                assert_eq!(
                    ranked_cache.outgoing_continuation(idx, forward, oriented),
                    symmetric_cache.outgoing_continuation(idx, forward, oriented)
                );
                assert_eq!(
                    symmetric_cache.state_bits(idx, forward, word),
                    bits,
                    "symmetric k={k}"
                );
                if (bits >> EDGE_OUT_DEGREE_SHIFT) & 7 == 1 {
                    assert_eq!(
                        symmetric_cache.neighbor(idx, forward, oriented),
                        wide_cache.neighbor(idx, forward, oriented)
                    );
                }

                assert_eq!(wide_cache.is_internal(idx), bits & 0x3f == 9);
                assert_eq!(compact_cache.is_internal(idx), bits & 0x3f == 9);

                assert_eq!(
                    compact_cache.incoming_degree(idx, forward),
                    (bits & 7) as u8
                );
                assert_eq!(wide_cache.incoming_degree(idx, forward), (bits & 7) as u8);
                let expected = if (bits >> EDGE_OUT_DEGREE_SHIFT) & 7 == 1 {
                    let (index, orientation) = wide_cache.neighbor(idx, forward, oriented);
                    Some((
                        index,
                        orientation,
                        u64::from((bits >> EDGE_OUT_CODE_SHIFT) & 3),
                    ))
                } else {
                    None
                };
                assert_eq!(wide_cache.outgoing(idx, forward, oriented), expected);
                assert_eq!(ranked_cache.outgoing(idx, forward, oriented), expected);
                assert_eq!(relative_cache.outgoing(idx, forward, oriented), expected);
                assert_eq!(compact_cache.outgoing(idx, forward, oriented), expected);
                assert_eq!(symmetric_cache.outgoing(idx, forward, oriented), expected);
                assert_eq!(
                    symmetric_cache.outgoing_continuation(idx, forward, oriented),
                    expected.filter(|&(target, orientation, _)| symmetric_cache
                        .incoming_degree(target, orientation)
                        == 1)
                );

                if ((wide_cache.state_bits(idx, forward, word) >> EDGE_OUT_DEGREE_SHIFT)
                    & EDGE_IN_DEGREE_MASK)
                    == 1
                {
                    assert_eq!(
                        compact_cache.neighbor(idx, forward, oriented),
                        wide_cache.neighbor(idx, forward, oriented)
                    );
                }
            }

            for (orientation, oriented) in [word, reverse_complement_word_slow(word, k)]
                .into_iter()
                .enumerate()
            {
                let bits = (cached.packed >> (orientation * EDGE_STATE_BITS))
                    & ((1 << EDGE_STATE_BITS) - 1);
                if (bits >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK == 1 {
                    let code = (bits >> EDGE_OUT_CODE_SHIFT) & 3;
                    let next = shift_append_word(oriented, u64::from(code), k, mask);
                    let canonical = next.min(reverse_complement_word_slow(next, k));
                    assert_eq!(entries[cached.next[orientation] as usize], canonical);
                    assert_eq!(
                        cached.packed & (1 << (2 * EDGE_STATE_BITS + orientation)) != 0,
                        next == canonical
                    );
                }
            }
            let expected =
                oracle(word) | (oracle(reverse_complement_word_slow(word, k)) << EDGE_STATE_BITS);
            assert_eq!(
                pack_edge_states(word, &entries, &lookup, k, mask),
                expected,
                "k={k}, word={word}"
            );
        }
    }
}
