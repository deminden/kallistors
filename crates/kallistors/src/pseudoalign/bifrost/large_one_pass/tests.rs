use super::*;

#[test]
fn large_one_pass_codes_sorting_and_runtime_flags_preserve_wide_sources() {
    let mut events = Vec::new();
    let mut records = Vec::new();
    for bucket in [0, 1, 255, 256, 1 << 28, u32::MAX] {
        for uid in [0, (1 << 21) + 17, UID_LIMIT - 1] {
            for rel in [0, 1, 65535] {
                let source = source_code(uid, rel);
                assert_eq!(position(source), ((uid as u64) << 32) | rel as u64);
                events.push(Event::new(bucket, source));
                for fp in [0, 1, 255, 256, 4094, 4095] {
                    let item = Record::new(bucket, source, fp);
                    assert_eq!(item.source(), source);
                    assert_eq!(item.fingerprint(), fp);
                    records.push(item);
                }
            }
        }
    }
    events.reverse();
    records.reverse();
    let mut expected = events
        .iter()
        .map(|event| (event.bucket, event.source()))
        .collect::<Vec<_>>();
    expected.sort_unstable();
    events
        .radix_sort_builder()
        .with_parallel(false)
        .with_single_threaded_tuner()
        .sort();
    assert_eq!(
        events
            .iter()
            .map(|event| (event.bucket, event.source()))
            .collect::<Vec<_>>(),
        expected
    );
    records
        .radix_sort_builder()
        .with_parallel(false)
        .with_single_threaded_tuner()
        .sort();
    assert!(
        records
            .windows(2)
            .all(|items| items[0].bucket <= items[1].bucket)
    );
    let pool = ThreadPoolBuilder::new().num_threads(3).build().unwrap();
    let mut builder = Builder::new(7, 8, 8);
    let cases = [
        (2, source_code(3, 7)),
        (1, source_code(UID_LIMIT - 1, 65535)),
        (2, source_code(0, 0)),
        (2, source_code(1, 9)),
        (4, source_code(0, 0)),
        (2, source_code(0, 0)),
    ];
    for &(bucket, source) in cases.iter().rev() {
        record(
            &mut builder.states[bucket],
            &mut builder.events,
            bucket as u32,
            source,
        );
    }
    let (result, fingerprints) = builder.finish(&pool);
    assert!(
        fingerprints
            .iter()
            .all(|slot| slot.load(Ordering::Relaxed) == 0)
    );
    assert_eq!(result.get(0), &[]);
    assert_eq!(
        result.get(1),
        &[position(source_code(UID_LIMIT - 1, 65535))]
    );
    assert!(
        !result.get_with_overcrowded(1).1,
        "ordinary bit61 is data, not a multi-row flag"
    );
    assert_eq!(result.get(2), &[0, 0, (1 << 32) | 9, (3 << 32) | 7]);
    assert_eq!(
        result.get(4),
        &[0],
        "real position0 cannot be an empty-state sentinel"
    );
}

#[test]
fn large_one_pass_resource_policy_preserves_bounded_fallback() {
    assert_eq!(size_eligible(370_014_451, 406_910_794), usize::BITS == 64);
    assert!(!size_eligible(49_948_538, 53_728_714));
    assert!(
        !size_eligible(370_014_451, 500_000_000),
        "repetitive references retain the previous loader"
    );
    assert!(!size_eligible(370_014_451, INPUT_LIMIT));
    assert!(!size_eligible(370_014_451, 370_014_450));
    assert!(!size_eligible(u32::MAX as usize, INPUT_LIMIT - 1));
    let mut builder = Builder::new(3, 10, 6);
    record(
        &mut builder.states[0],
        &mut builder.events,
        0,
        source_code(0, 0),
    );
    record(
        &mut builder.states[0],
        &mut builder.events,
        0,
        source_code(1, 0),
    );
    assert_eq!(builder.events.len(), 2);
    assert_eq!(builder.events.capacity(), 6);
    builder.retire();
}

#[test]
fn large_one_pass_matches_independent_ascii_stream_across_batches_and_fallback() {
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("transcripts.fa");
    let mut random = 173u64;
    let mut sequence = |length| {
        (0..length)
            .map(|_| {
                random ^= random << 13;
                random ^= random >> 7;
                random ^= random << 17;
                b"ACGT"[(random & 3) as usize]
            })
            .collect::<Vec<_>>()
    };
    let source = sequence(5000);
    std::fs::write(
        &fasta,
        format!(">t\n{}\n", std::str::from_utf8(&source).unwrap()),
    )
    .unwrap();
    let path = directory.path().join("index.idx");
    crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default()).unwrap();
    let bytes = std::fs::read(path).unwrap();
    let graph = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let mphf = BooPhf::load(&mut std::io::Cursor::new(&bytes[16 + graph + 8..])).unwrap();
    let sequences = vec![
        Vec::new(),
        sequence(22),
        source[..31].to_vec(),
        source[..65].to_vec(),
        source[..900].to_vec(),
        source[..900].to_vec(),
        source[..1100].to_vec(),
    ];
    for (layout, workers) in [(0, 1), (0, 2), (1, 3), (2, 8)] {
        let mut encoded = precompute_encoded_unitigs(&sequences, &[], 23, None);
        if layout == 1 {
            let SequenceMetadata::Compact(items) = &encoded.unitig_bases.metadata else {
                panic!("compact fixture")
            };
            encoded.unitig_bases.metadata = SequenceMetadata::Wide(
                items
                    .iter()
                    .map(|&value| (value as u32 as usize, (value >> 32) as usize))
                    .collect(),
            );
        } else if layout == 2 {
            encoded.unitig_bases.inline_short_words();
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        let ctx = Context {
            encoded: &encoded,
            mphf: &mphf,
            pool: &pool,
            buckets: mphf.size() as usize,
            k: 31,
        };
        let mut selected = Vec::new();
        let mut ordinal = 0u32;
        let mut expected = vec![Vec::new(); ctx.buckets];
        let fp = (0..ctx.buckets)
            .map(|_| AtomicU16::new(0))
            .collect::<Vec<_>>();
        // No packed metadata/source-code helper is used for the expected
        // payloads or hashes. Repeated reference sequence at different UIDs
        // forces true multi rows spanning batches and ownership ranges.
        for (uid, bases) in sequences.iter().enumerate() {
            for (relative, query) in bases.windows(23).enumerate() {
                if ordinal % 5 != 2 {
                    selected.push(ordinal);
                    let key = crate::index::bifrost::encode_minimizer_rep(query).unwrap();
                    if let Some(bucket) = mphf.debug_lookup(&key).0 {
                        expected[bucket as usize].push(((uid as u64) << 32) | relative as u64);
                        fingerprints::observe(&fp[bucket as usize], u64::from_le_bytes(key));
                    }
                }
                ordinal += 1;
            }
        }
        assert!(selected.len() > 3 * 512);
        let bitmap = BitmapPositions::new(selected.clone(), true);
        // Internal forced-wide path avoids allocating millions of empty
        // buckets merely to exercise production batch and source-order logic.
        let (actual, actual_fp) =
            try_build_sized::<64, 512, 4>(&ctx, &bitmap, bitmap.len()).unwrap();
        for bucket in 0..ctx.buckets {
            assert_eq!(
                actual.get(bucket),
                expected[bucket],
                "layout={layout} workers={workers} bucket={bucket}"
            );
            assert!(!actual.get_with_overcrowded(bucket).1);
            assert_eq!(
                actual_fp[bucket].load(Ordering::Relaxed),
                fp[bucket].load(Ordering::Relaxed)
            );
        }
        assert_eq!(bitmap.iter().collect::<Vec<_>>(), selected);
        // The first512-position batch fits; later worst-case growth refuses
        // before exceeding the budget. Caller-owned bitmap remains usable by
        // the established two-pass fallback, with no consumed prefix or FP state.
        assert!(try_build_sized::<8, 512, 4>(&ctx, &bitmap, 700).is_none());
        assert_eq!(bitmap.iter().collect::<Vec<_>>(), selected);
    }
}

#[test]
fn large_one_pass_first_publication_matches_atomic_fingerprint_observer() {
    for first in 0..=fingerprints::DISABLED {
        for second in [0, 1, first, first ^ 1, fingerprints::DISABLED] {
            for source in [0, source_code(UID_LIMIT - 1, 65535)] {
                let expected = AtomicU16::new(0);
                let mut state = 0;
                let mut events = Vec::new();
                for (step, (position, fingerprint)) in
                    [(source, first), (source, second), (source ^ 65536, first)]
                        .into_iter()
                        .enumerate()
                {
                    // The established observer is independent of the packed
                    // descriptor and consensus representation under test.
                    fingerprints::observe_serial(&expected, fingerprint);
                    record_fingerprinted(&mut state, &mut events, 7, position, fingerprint);
                    assert_eq!(
                        ((state & FINGERPRINT_MASK) >> FINGERPRINT_SHIFT) as u16,
                        expected.load(Ordering::Relaxed)
                    );
                    if step == 0 {
                        assert_eq!(state & !FINGERPRINT_MASK, source + 1);
                    }
                }
                assert_eq!(state & !FINGERPRINT_MASK, MULTI | 3);
                assert_eq!(
                    events
                        .iter()
                        .map(|event| event.source())
                        .collect::<Vec<_>>(),
                    [source, source, source ^ 65536]
                );
            }
        }
    }
}

#[test]
fn large_one_pass_packed_fingerprints_preserve_empty_and_maximum_singletons() {
    for workers in [1, 3] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        for fingerprint in [0, 1, 4094, fingerprints::DISABLED] {
            let mut builder = Builder::new(7, 8, 8);
            let source = source_code(UID_LIMIT - 1, 65535);
            record_fingerprinted(
                &mut builder.states[1],
                &mut builder.events,
                1,
                source,
                fingerprint,
            );
            record_fingerprinted(
                &mut builder.states[4],
                &mut builder.events,
                4,
                0,
                fingerprint,
            );
            let (positions, fingerprints) = builder.finish(&pool);
            assert_eq!(positions.get(1), &[position(source)]);
            assert!(!positions.get_with_overcrowded(1).1);
            assert_eq!(positions.get(4), &[0]);
            let expected = if fingerprint == 0 {
                fingerprints::DISABLED
            } else {
                fingerprint
            };
            for (bucket, slot) in fingerprints.iter().enumerate() {
                assert_eq!(
                    slot.load(Ordering::Relaxed),
                    if [1, 4].contains(&bucket) {
                        expected
                    } else {
                        0
                    }
                );
            }
            for bucket in [0, 2, 3, 5, 6] {
                assert!(positions.get(bucket).is_empty());
            }
        }
    }
}

#[test]
fn large_one_pass_owner_routing_matches_independent_stable_intervals() {
    let snapshot = |records: &[Record]| {
        records
            .iter()
            .map(|record| (record.bucket, record.source(), record.fingerprint()))
            .collect::<Vec<_>>()
    };
    for buckets in [1usize, 3, 17, 256, 65535, 1 << 28] {
        for workers in [2usize, 3, 8, 17] {
            let width = buckets.div_ceil(workers).max(1);
            let owners = buckets.div_ceil(width);
            let choices = [
                0,
                (buckets - 1) as u32,
                (buckets / 2) as u32,
                buckets as u32,
                u32::MAX,
                1,
            ];
            let records: Vec<_> = (0..257)
                .map(|i| {
                    Record::new(
                        choices[(i * 17 + i / 7) % choices.len()],
                        source_code(i / 31, i % 31),
                        [0, 1, 255, 4094, 4095][i % 5],
                    )
                })
                .collect();
            for length in [0, 1, 2, 7, 8, 31, 32, 33, 256, 257] {
                let input = records[..length].to_vec();
                // Independent interval filtering: no histogram, owner key,
                // prefix offsets, or scatter code shared with route_records.
                let mut expected = Vec::new();
                for group in 0..owners {
                    let first = group * width;
                    let end = (first + width).min(buckets);
                    expected.extend(input.iter().copied().filter(|record| {
                        let bucket = record.bucket as usize;
                        bucket >= first && bucket < end
                    }));
                }
                expected.extend(
                    input
                        .iter()
                        .copied()
                        .filter(|record| record.bucket as usize >= buckets),
                );
                let actual = route_records(input, width, buckets);
                assert_eq!(snapshot(&actual), snapshot(&expected));
                for group in 0..owners {
                    let first = group * width;
                    let end = (first + width).min(buckets);
                    let lo = actual.partition_point(|record| (record.bucket as usize) < first);
                    let hi =
                        lo + actual[lo..].partition_point(|record| (record.bucket as usize) < end);
                    let independent: Vec<_> = records[..length]
                        .iter()
                        .copied()
                        .filter(|record| {
                            let bucket = record.bucket as usize;
                            bucket >= first && bucket < end
                        })
                        .collect();
                    assert_eq!(snapshot(&actual[lo..hi]), snapshot(&independent));
                }
            }
        }
    }
}

#[test]
fn large_one_pass_narrow_relative_offsets_preserve_supported_domain() {
    for length in [0, 22, 23, 31, 32, 65535] {
        assert!(unitig_length_eligible(length));
    }
    for length in [65536, 65537, usize::MAX] {
        assert!(!unitig_length_eligible(length));
    }
    for start in [0, 1, u64::from(u32::MAX) - 65535] {
        for offset in [0u64, 1, 22, 23, 31, 32, 65512, 65535] {
            assert_eq!(relative_offset(start + offset, start) as u64, offset);
        }
    }
}
