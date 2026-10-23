//! Compaction regressions.

use super::*;
#[test]
fn oriented_start_prefix_ranges_preserve_unique_duplicate_and_missing_words() {
    let mut records = Vec::new();
    for i in 0..2048u64 {
        for copy in 0..=i % 3 {
            records.push(OrientedStart {
                word: (i << 52) | (copy & 1),
                unitig: (i * 3 + copy) as u32,
            });
        }
    }
    let identity = |v: Option<OrientedStart>| v.map(|v| (v.word(), v.unitig_id(), v.forward()));
    let mut expected = HashMap::<u64, Option<OrientedStart>>::new();
    for &record in &records {
        expected
            .entry(record.word())
            .and_modify(|entry| *entry = None)
            .or_insert(Some(record));
    }
    for threads in [1, 8] {
        let lookup = OrientedStartLookup::new(&records, threads);
        for i in 0..4096u64 {
            for word in [i << 52, (i << 52) | (1 << 40)] {
                assert_eq!(
                    identity(lookup.find(word)),
                    identity(expected.get(&word).copied().flatten())
                );
            }
        }
    }
    for records in [&records[..0], &records[..1], &records[..17]] {
        let lookup = OrientedStartLookup::new(records, 1);
        for word in [0, 1 << 52, u64::MAX] {
            assert_eq!(
                identity(lookup.find(word)),
                identity(find_unique_oriented_start(records, word))
            );
        }
    }
}

#[test]
fn shared_transcript_ends_preserve_every_seam_window() {
    // Exact differing-byte positions are the proof boundary, including
    // shorter slices and every 8/32-byte seam in either direction.
    for length in 0usize..97 {
        let a = (0..length).map(|i| (i * 37) as u8).collect::<Vec<_>>();
        for other_length in [length, length.saturating_sub(3)] {
            for change in 0..=other_length {
                let mut b = a[..other_length].to_vec();
                if change < b.len() {
                    b[change] ^= 0x80;
                }
                assert_eq!(
                    shared_prefix_bytes(&a, &b),
                    a.iter().zip(&b).take_while(|(a, b)| a == b).count()
                );
                assert_eq!(
                    shared_suffix_bytes(&a, &b),
                    a.iter()
                        .rev()
                        .zip(b.iter().rev())
                        .take_while(|(a, b)| a == b)
                        .count()
                );
            }
        }
    }
    let mut seed = 17u64;
    let mut sequence = |length| {
        (0..length)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                code_to_base(seed >> 62)
            })
            .collect::<Vec<_>>()
    };
    let interior = sequence(4096);
    let mut first = sequence(89);
    first.extend_from_slice(&interior);
    first.extend(sequence(47));
    let mut second = sequence(137);
    second.extend_from_slice(&interior);
    second.extend(sequence(23));
    let mut third = second.clone();
    third[700] = code_to_base((u64::from(BASE_CODES[third[700] as usize]) + 1) & 3);
    for k in (3..=31).step_by(2) {
        let transcripts = [first.as_slice(), second.as_slice(), third.as_slice()];
        let expected = transcripts
            .iter()
            .flat_map(|seq| seq.windows(k))
            .map(|window| {
                let (a, b) = encode_kmer_pair(window).unwrap();
                a.min(b)
            })
            .collect::<HashSet<_>>();
        let mut anchors = AnchorMap::default();
        let mut sparse = Vec::new();
        let mut every_anchor = AnchorMap::default();
        let mut dense = Vec::new();
        for (id, seq) in transcripts.into_iter().enumerate() {
            insert_interior_unique_kmers::<4>(&mut sparse, seq, k, id, &mut anchors).unwrap();
            insert_interior_unique_kmers::<0>(&mut dense, seq, k, id, &mut every_anchor).unwrap();
        }
        assert_eq!(
            sparse.iter().copied().collect::<HashSet<_>>(),
            expected,
            "sparse interior k={k}"
        );
        assert_eq!(
            dense.iter().copied().collect::<HashSet<_>>(),
            expected,
            "dense interior k={k}"
        );
        assert!(sparse.len() < first.len() + second.len() + third.len() - 3 * (k - 1));
    }
    let previous = b"ACGTTGCATGTCAGTACCTGATCGATGCTAGCATCGTACTGATGCGTACGTAGCTAGCATGCTAGT";
    for k in (5..=31).step_by(2) {
        for prefix in [0, k - 1, k, k + 1, 40] {
            for suffix in [0, k - 1, k, k + 1, 40] {
                let mut current = previous[..prefix].to_vec();
                current.extend_from_slice(b"GTACCGATTCGA");
                current.extend_from_slice(&previous[previous.len() - suffix..]);
                let words = |seq: &[u8]| {
                    seq.windows(k)
                        .map(|word| {
                            let (forward, reverse) = encode_kmer_pair(word).unwrap();
                            forward.min(reverse)
                        })
                        .collect::<HashSet<_>>()
                };
                let mut expected = words(previous);
                expected.extend(words(&current));
                let mut actual = words(previous);
                actual.extend(words(novel_transcript_region(&current, previous, k)));
                assert_eq!(actual, expected, "k={k}, prefix={prefix}, suffix={suffix}");
                let unrelated = vec![b'C'; previous.len()];
                let history = [
                    previous.to_vec(),
                    unrelated,
                    current.clone(),
                    previous.to_vec(),
                ]
                .into_iter()
                .enumerate()
                .map(|(id, seq)| Transcript {
                    name: id.to_string(),
                    original_len: seq.len() as u32,
                    seq,
                })
                .collect::<Vec<_>>();
                let mut anchored = HashSet::new();
                let mut prefixes = AnchorMap::default();
                let mut suffixes = AnchorMap::default();
                for i in 0..history.len() {
                    anchored.extend(words(anchored_transcript_region(
                        &history,
                        i,
                        k,
                        &mut prefixes,
                        &mut suffixes,
                    )));
                }
                let all = history
                    .iter()
                    .flat_map(|transcript| words(&transcript.seq))
                    .collect::<HashSet<_>>();
                assert_eq!(anchored, all, "anchor seam k={k}");
                let mut observed = HashSet::new();
                let mut complete = HashSet::new();
                for (i, transcript) in history.iter().enumerate() {
                    observed.extend(words(novel_transcript_history(
                        &transcript.seq,
                        &history[..i],
                        k,
                    )));
                    complete.extend(words(&transcript.seq));
                }
                assert_eq!(
                    observed, complete,
                    "history seam k={k},prefix={prefix},suffix={suffix}"
                );
            }
        }
        assert!(novel_transcript_region(previous, previous, k).is_empty());
    }
    let mut transcripts = (0..3073)
        .map(|id| {
            let mut seq = b"ACGT".repeat(96);
            seq[id % 384] = b"TGCA"[id & 3];
            Transcript {
                name: id.to_string(),
                original_len: seq.len() as u32,
                seq,
            }
        })
        .collect::<Vec<_>>();
    for seq in [Vec::new(), b"ACG".to_vec()] {
        transcripts.push(Transcript {
            name: transcripts.len().to_string(),
            original_len: seq.len() as u32,
            seq,
        });
    }
    let expected = collect_kmers(&transcripts, 31, 1).unwrap().sorted;
    for budget in [64, 65536, usize::MAX] {
        let mut owned = transcripts.clone();
        let mut actual = Vec::new();
        let mut identifiers = Vec::new();
        for group in bounded_transcript_groups(&mut owned, 31, budget) {
            assert!(!group.is_empty() && group.len() <= 2048);
            let before_last = group[..group.len() - 1]
                .iter()
                .map(|tr| kmer_position_count(tr.seq.len(), 31))
                .sum::<usize>();
            assert!(before_last < budget);
            let mut prefixes = AnchorMap::default();
            let mut suffixes = AnchorMap::default();
            let mut anchors = AnchorMap::default();
            for (id, transcript) in group.iter().enumerate() {
                identifiers.push(transcript.name.clone());
                let sequence =
                    anchored_transcript_region(group, id, 31, &mut prefixes, &mut suffixes);
                insert_interior_unique_kmers::<4>(&mut actual, sequence, 31, id, &mut anchors)
                    .unwrap();
            }
        }
        actual.sort_unstable();
        actual.dedup();
        assert_eq!(actual, expected, "bounded groups budget={budget}");
        assert_eq!(
            identifiers,
            transcripts
                .iter()
                .map(|tr| tr.name.clone())
                .collect::<Vec<_>>()
        );
    }

    let mut owned = transcripts.clone();
    let (actual, packed) = collect_kmers_and_pack(&mut owned, 31, 8).unwrap();
    assert_eq!(actual.sorted, expected);
    assert_eq!(packed.len(), transcripts.len());
    for ((packed, source), retired) in packed.iter().zip(&transcripts).zip(&owned) {
        assert_eq!(packed.len(), source.seq.len());
        assert_eq!(retired.name, source.name);
        assert_eq!(retired.original_len, source.original_len);
        assert!(retired.seq.is_empty());
        for (index, &base) in source.seq.iter().enumerate() {
            assert_eq!(packed.base(index), base);
        }
    }
}

#[test]
fn parallel_path_preparation_preserves_serial_claim_order() {
    assert_eq!(std::mem::size_of::<JoinQuery>(), 13);
    for source in [0u64, (1 << 33) - 1] {
        if let Ok(source) = usize::try_from(source) {
            let query = JoinQuery::new(source, u64::MAX).unwrap();
            assert_eq!(query.source(), source);
            assert_eq!(query.word(), u64::MAX);
        }
    }
    if let Ok(overflow) = usize::try_from(1u64 << 40) {
        assert!(JoinQuery::new(overflow, 0).is_err());
    }

    assert!(sparse_locations::worthwhile(1 << 30, 1 << 20, 1 << 20));
    assert!(!sparse_locations::worthwhile(1 << 30, 1 << 29, 1));

    assert_eq!(join_state(&JoinLinks::Wide(vec![usize::MAX]), 0), None);
    assert_eq!(join_state(&JoinLinks::Wide(vec![0]), 0), Some(0));
    let paths = (0..1200)
        .map(|id| UnitigPath {
            ends: [u32::MAX; 2],
            seq: if id & 1 == 0 {
                Box::from(&b"ACGTTGCATGCA"[..])
            } else {
                Box::from(&b"TGCATGCAACGT"[..])
            },
        })
        .collect::<Vec<_>>();
    for mode in 0..3 {
        let active = (0..paths.len())
            .map(|id| match mode {
                0 => true,
                1 => false,
                _ => id % 3 != 0,
            })
            .collect::<Vec<_>>();
        let mut serial = build_oriented_starts(&paths, &active, 5, 1)
            .unwrap()
            .into_iter()
            .map(|entry| (entry.word(), entry.unitig_id(), entry.forward()))
            .collect::<Vec<_>>();
        let mut parallel = build_oriented_starts(&paths, &active, 5, 8)
            .unwrap()
            .into_iter()
            .map(|entry| (entry.word(), entry.unitig_id(), entry.forward()))
            .collect::<Vec<_>>();
        serial.sort_unstable();
        parallel.sort_unstable();
        assert_eq!(parallel, serial);
    }
    let mut invalid_paths = paths
        .iter()
        .map(|path| UnitigPath {
            ends: [u32::MAX; 2],
            seq: path.seq.clone(),
        })
        .collect::<Vec<_>>();
    invalid_paths[713].seq[0] = b'N';
    assert!(build_oriented_starts(&invalid_paths, &vec![true; invalid_paths.len()], 5, 8).is_err());
    let wide_start = OrientedStart {
        word: 0x9abc_def0_1234_5671,
        unitig: u32::MAX,
    };
    assert_eq!(wide_start.unitig_id(), u32::MAX as usize);
    assert!(wide_start.forward());
    assert_eq!(wide_start.word(), 0x9abc_def0_1234_5670);

    assert_eq!(join_state(&JoinLinks::Narrow(vec![u32::MAX]), 0), None);
    assert_eq!(join_state(&JoinLinks::Narrow(vec![0]), 0), Some(0));
    assert_eq!(
        join_state(&JoinLinks::Narrow(vec![u32::MAX - 1]), 0),
        Some(u32::MAX as usize - 1)
    );
    assert_eq!(
        join_state(&JoinLinks::Wide(vec![usize::MAX - 1]), 0),
        Some(usize::MAX - 1)
    );

    for index in [0, 1, (1 << 30) - 1] {
        for base in 0..4 {
            let narrow = <u32 as PathStep>::pack(index, base);
            let wide = <u64 as PathStep>::pack(index, base);
            assert_eq!(narrow.index(), index);
            assert_eq!(wide.index(), index);
            assert_eq!(narrow.base(), base);
            assert_eq!(wide.base(), base);
        }
    }
    let wide = <u64 as PathStep>::pack(u32::MAX as usize, 3);
    assert_eq!(wide.index(), u32::MAX as usize);

    let mut cycle_rng = 987654321u64;
    let mut ring = (0..10000)
        .map(|_| {
            cycle_rng ^= cycle_rng << 13;
            cycle_rng ^= cycle_rng >> 7;
            cycle_rng ^= cycle_rng << 17;
            b"ACGT"[(cycle_rng & 3) as usize]
        })
        .collect::<Vec<_>>();
    ring.extend_from_within(..30);
    let cycle = vec![Transcript {
        name: "bounded-cycle".into(),
        original_len: ring.len() as u32,
        seq: ring,
    }];
    let words = collect_kmers(&cycle, 31, 1).unwrap();
    let lookup = KmerLookup::new(&words.sorted).unwrap();
    let edges = EdgeStateCache::build(&words.sorted, &lookup, 31, 1).unwrap();
    let topology = GraphTopology::new(&words.sorted, &lookup, &edges, 31);
    let start = choose_any_start_orientation(words.sorted[0], 0, &topology);
    assert!(
        prepare_path_group::<u64, true>(
            &[(0, start)],
            &mut [Vec::new()],
            &topology,
            &VisitedKmers::new(words.sorted.len())
        )[0]
        .is_none()
    );
    let serial = compact_unitigs(&words.sorted, &topology, 1)
        .unwrap()
        .into_iter()
        .map(|path| path.seq)
        .collect::<Vec<_>>();
    let parallel = compact_unitigs_parallel(&words.sorted, &topology, 8)
        .unwrap()
        .into_iter()
        .map(|path| path.seq)
        .collect::<Vec<_>>();
    assert_eq!(parallel, serial, "bounded cycle fallback");
    let wide = compact_unitigs_parallel_steps::<u64>(&words.sorted, &topology, 8)
        .unwrap()
        .into_iter()
        .map(|path| path.seq)
        .collect::<Vec<_>>();
    assert_eq!(wide, serial, "wide bounded cycle fallback");
    for batch in [512, 4096] {
        let mut visited = VisitedKmers::new(words.sorted.len());
        let mut paths = PathChunks::default();
        let workers = IndexWorkers::new(8).unwrap();
        finish_unitig_cycles_batched::<u64>(
            &words.sorted,
            &topology,
            &mut visited,
            &mut paths,
            &workers,
            batch,
        )
        .unwrap();
        let actual = paths
            .finish()
            .into_iter()
            .map(|path| path.seq)
            .collect::<Vec<_>>();
        assert_eq!(actual, serial, "ordered cycle batch {batch}");
    }

    let relative_edges = EdgeStateCache {
        packed: EdgeStorage::Relative(RelativeEdges::build(
            &words.sorted,
            &lookup,
            31,
            RelativeEdges::prefixes(&words.sorted, &lookup, 31).unwrap(),
            true,
        )),
    };
    let relative_topology = GraphTopology::new(&words.sorted, &lookup, &relative_edges, 31);
    let relative = compact_unitigs_parallel_steps::<u64>(&words.sorted, &relative_topology, 8)
        .unwrap()
        .into_iter()
        .map(|path| path.seq)
        .collect::<Vec<_>>();
    assert_eq!(relative, serial, "relative bounded cycle fallback");

    let mut rng = 42u64;
    let mut transcripts = Vec::new();
    for id in 0..1200 {
        let seq = (0..96)
            .map(|_| {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                b"ACGT"[(rng & 3) as usize]
            })
            .collect::<Vec<_>>();
        transcripts.push(Transcript {
            name: id.to_string(),
            original_len: seq.len() as u32,
            seq,
        });
    }
    let first = transcripts[0].seq.clone();
    for (id, seq) in [
        first[7..83].to_vec(),
        reverse_complement_seq(&first[11..89]),
    ]
    .into_iter()
    .enumerate()
    {
        transcripts.push(Transcript {
            name: format!("contained-{id}"),
            original_len: seq.len() as u32,
            seq,
        });
    }

    transcripts.push(Transcript {
        name: "reverse".into(),
        original_len: first.len() as u32,
        seq: reverse_complement_seq(&first),
    });
    transcripts.push(Transcript {
        name: "cycle".into(),
        original_len: 256,
        seq: b"ACGT".repeat(64),
    });
    for k in [5, 15, 31] {
        let words = collect_kmers(&transcripts, k, 1).unwrap();
        let lookup = KmerLookup::new(&words.sorted).unwrap();
        let edges = EdgeStateCache::build(&words.sorted, &lookup, k, 1).unwrap();
        let topology = GraphTopology::new(&words.sorted, &lookup, &edges, k);
        let expected = compact_unitigs(&words.sorted, &topology, 1)
            .unwrap()
            .into_iter()
            .map(|path| path.seq)
            .collect::<Vec<_>>();
        let relative_edges = EdgeStateCache {
            packed: EdgeStorage::Relative(RelativeEdges::build(
                &words.sorted,
                &lookup,
                k,
                RelativeEdges::prefixes(&words.sorted, &lookup, k).unwrap(),
                true,
            )),
        };
        let packed_edges = EdgeStateCache {
            packed: EdgeStorage::Relative7(PackedRelativeEdges::build(
                &words.sorted,
                &lookup,
                k,
                PackedRelativeEdges::prefixes(&words.sorted, &lookup, k).unwrap(),
                true,
            )),
        };
        let dense_edges = EdgeStateCache {
            packed: EdgeStorage::Relative6(DenseRelativeEdges::build(
                &words.sorted,
                &lookup,
                k,
                DenseRelativeEdges::prefixes(&words.sorted, &lookup, k).unwrap(),
                true,
            )),
        };
        let relative_topology = GraphTopology::new(&words.sorted, &lookup, &relative_edges, k);
        let packed_topology = GraphTopology::new(&words.sorted, &lookup, &packed_edges, k);
        let dense_topology = GraphTopology::new(&words.sorted, &lookup, &dense_edges, k);
        let adaptive_edges = EdgeStateCache {
            packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_plan(
                &words.sorted,
                &lookup,
                k,
                adaptive_prefixes::PrefixPlan::new(&words.sorted, &lookup, k, 1, 15),
                true,
            )),
        };
        let adaptive_topology = GraphTopology::new(&words.sorted, &lookup, &adaptive_edges, k);

        if k >= 15 {
            // Force the large exact key geometry on this serial-order oracle.
            // Sparse and empty prefix ranges must reconstruct full words too.
            let mut rank_lookup = KmerLookup::with_threads(&words.sorted, 1, k).unwrap();
            let bits = (2 * k).saturating_sub(40).max(12);
            rank_lookup.shift = 64 - bits;
            let mut cursor = 0;
            let offsets = (0..=(1usize << bits))
                .map(|bucket| {
                    while cursor < words.sorted.len()
                        && words.sorted[cursor] >> rank_lookup.shift < bucket as u64
                    {
                        cursor += 1;
                    }
                    cursor as u32
                })
                .collect();
            rank_lookup.directory = KmerDirectory::Ranks(RankDirectory::new(offsets, false));
            let mut keys = PrefixKeys::for_topology(words.sorted.clone(), &rank_lookup, 1).unwrap();
            let filter =
                negative_filter::NegativeFilter::build_mapped(&keys, &rank_lookup, k, 8).unwrap();
            let packed_edges = EdgeStateCache {
                packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_filtered(
                    &keys,
                    &rank_lookup,
                    k,
                    adaptive_prefixes::PrefixPlan::new(&words.sorted, &rank_lookup, k, 1, 128),
                    &filter,
                    true,
                )),
            };
            filter.retire();
            let packed_topology = GraphTopology::new(&keys, &rank_lookup, &packed_edges, k);
            for threads in [1, 8] {
                let paths = if threads == 1 {
                    compact_unitigs(&keys, &packed_topology, threads)
                } else {
                    compact_unitigs_parallel(&keys, &packed_topology, threads)
                }
                .unwrap();
                assert_eq!(
                    paths.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                    expected.iter().collect::<Vec<_>>(),
                    "packed-key paths k={k}, threads={threads}"
                );
                let joined = join_unitig_paths(paths, &packed_topology, threads).unwrap();
                let raw_paths = compact_unitigs(&words.sorted, &topology, 1).unwrap();
                let raw_joined = join_unitig_paths(raw_paths, &topology, 1).unwrap();
                assert_eq!(
                    joined.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                    raw_joined.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                    "packed-key joins k={k}, threads={threads}"
                );
            }
            let raw_paths = compact_unitigs(&words.sorted, &topology, 1).unwrap();
            let normalize = |groups: Vec<Vec<JoinQuery>>| {
                groups
                    .into_iter()
                    .flatten()
                    .map(|query| (query.source(), query.word()))
                    .collect::<Vec<_>>()
            };
            let expected_probes =
                normalize(collect_join_queries(&raw_paths, &packed_topology, 1).unwrap());
            keys.retire_word_hints();
            let retired_topology = GraphTopology::new(&keys, &rank_lookup, &packed_edges, k);
            let actual_probes =
                normalize(collect_join_queries(&raw_paths, &retired_topology, 1).unwrap());
            assert_eq!(
                actual_probes, expected_probes,
                "joins after hint retirement"
            );
        }

        for threads in [2, 8] {
            let actual = compact_unitigs_parallel(&words.sorted, &topology, threads)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(actual, expected, "k={k}, threads={threads}");
            let relative_paths =
                compact_unitigs_parallel_steps::<u64>(&words.sorted, &relative_topology, threads)
                    .unwrap();
            assert_eq!(
                relative_paths
                    .iter()
                    .map(|path| &path.seq)
                    .collect::<Vec<_>>(),
                expected.iter().collect::<Vec<_>>(),
                "relative paths k={k}, threads={threads}"
            );
            let packed_paths =
                compact_unitigs_parallel_steps::<u64>(&words.sorted, &packed_topology, threads)
                    .unwrap();
            assert_eq!(
                packed_paths
                    .iter()
                    .map(|path| &path.seq)
                    .collect::<Vec<_>>(),
                expected.iter().collect::<Vec<_>>(),
                "packed paths k={k}, threads={threads}"
            );
            let dense_paths =
                compact_unitigs_parallel_steps::<u64>(&words.sorted, &dense_topology, threads)
                    .unwrap();
            assert_eq!(
                dense_paths.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                expected.iter().collect::<Vec<_>>(),
                "dense paths k={k}, threads={threads}"
            );
            let adaptive_paths =
                compact_unitigs_parallel_steps::<u64>(&words.sorted, &adaptive_topology, threads)
                    .unwrap();
            assert_eq!(
                adaptive_paths
                    .iter()
                    .map(|path| &path.seq)
                    .collect::<Vec<_>>(),
                expected.iter().collect::<Vec<_>>(),
                "adaptive paths k={k}, threads={threads}"
            );
            let paths = expected
                .iter()
                .map(|seq| UnitigPath {
                    seq: seq.clone(),
                    ends: [u32::MAX; 2],
                })
                .collect::<Vec<_>>();
            let serial_join = join_unitig_paths(paths.clone(), &topology, 1)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            let adaptive_join = join_unitig_paths(adaptive_paths, &adaptive_topology, threads)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(
                adaptive_join, serial_join,
                "adaptive join k={k}, threads={threads}"
            );
            let relative_join = join_unitig_paths(relative_paths, &relative_topology, threads)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(
                relative_join, serial_join,
                "relative join k={k}, threads={threads}"
            );
            let packed_join = join_unitig_paths(packed_paths, &packed_topology, threads)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(
                packed_join, serial_join,
                "packed join k={k}, threads={threads}"
            );
            let dense_join = join_unitig_paths(dense_paths, &dense_topology, threads)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(
                dense_join, serial_join,
                "dense join k={k}, threads={threads}"
            );
            let parallel_join = join_unitig_paths(paths, &topology, threads)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(parallel_join, serial_join, "join k={k}, threads={threads}");
            let joined = serial_join
                .into_iter()
                .map(|seq| UnitigPath {
                    seq,
                    ends: [u32::MAX; 2],
                })
                .collect::<Vec<_>>();
            // Exercise cross-unitig batching directly even below the public
            // parallel-size threshold, including group and tail seams.
            if threads == 2 {
                let mut expected_locations = Vec::new();
                for (unitig, path) in joined.iter().enumerate() {
                    store_unitig_locations(
                        path,
                        unitig,
                        k,
                        &words.sorted,
                        &lookup,
                        |_| {},
                        |rank, loc| expected_locations.push((rank, loc.0)),
                    )
                    .unwrap();
                }
                for group_size in [1, 63, 64, 65] {
                    let mut actual_locations = Vec::new();
                    for (group, paths) in joined.chunks(group_size).enumerate() {
                        store_location_group::<64, _, _, _>(
                            paths,
                            group * group_size,
                            k,
                            &words.sorted,
                            &lookup,
                            |_| {},
                            |rank, loc| actual_locations.push((rank, loc.0)),
                        )
                        .unwrap();
                    }
                    assert_eq!(
                        actual_locations, expected_locations,
                        "location groups k={k}, size={group_size}"
                    );
                }
            }
            let locs = build_kmer_locs(&joined, &words.sorted, &lookup, k, threads).unwrap();
            let mut inputs = transcripts.clone();
            let mut packed = pack_transcript_sequences(&mut inputs, threads).unwrap();
            let sparse =
                sparse_locations::SparseLocations::build(&joined, &packed, k, threads).unwrap();
            assert!(
                sparse.find(u64::MAX).is_none(),
                "empty slot is not membership"
            );
            let dense = DenseReferenceLocations {
                keys: &words.sorted,
                lookup: &lookup,
                locations: &locs,
            };
            let dense_shards = build_trinfos(&mut packed, k, &dense, &joined, threads).unwrap();
            let mut inputs = transcripts.clone();
            let mut packed = pack_transcript_sequences(&mut inputs, threads).unwrap();
            let sparse_shards = build_trinfos(&mut packed, k, &sparse, &joined, threads).unwrap();
            let counts = joined
                .iter()
                .map(|path| (path.seq.len() - k + 1) as u32)
                .collect::<Vec<_>>();
            let dense_rows = dense_shards
                .iter()
                .flat_map(|shard| shard.iter(&counts))
                .collect::<Vec<_>>();
            let sparse_rows = sparse_shards
                .iter()
                .flat_map(|shard| shard.iter(&counts))
                .collect::<Vec<_>>();
            assert_eq!(
                sparse_rows, dense_rows,
                "sparse exact runs k={k}, threads={threads}"
            );
            for path in &joined {
                let word = encode_kmer(&path.seq[..k]).unwrap();
                assert!(
                    sparse.find(word | 1).is_none(),
                    "complete-key padding check"
                );
            }
        }
    }
}

#[test]
fn joins_cyclic_repeat_unitigs() {
    let seq = b"ACGT".repeat(20);
    let tr = Transcript {
        name: "tx".into(),
        original_len: seq.len() as u32,
        seq,
    };
    let (graph, _) = build_kmer_unitig_graph_with_report(vec![tr], 31, 23, -1, 1).expect("graph");
    assert_eq!(graph.unitigs.len(), 1);
    assert_eq!(graph.km_unitigs.len(), 0);
    assert_eq!(graph.nodes.len(), graph.unitigs.len());
}
