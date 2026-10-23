use super::SequenceMetadata;
use super::{BASE_CODES, EncodedBaseStore, compute_minimizer_reps_into, minimizer_count};

#[test]
fn rolling_packed_windows_match_independent_windows() {
    let seqs = vec![
        (0..193)
            .map(|i| b"ACGT"[(i * i + i / 7) % 4])
            .collect::<Vec<_>>(),
    ];
    let store = EncodedBaseStore::from_sequences(&seqs, None);
    for k in 1..=32 {
        for first in 0..=193 - k {
            let last = (first + 64 - k).min(193 - k);
            let got: Vec<_> = store.kmer_range(0, first, last, k).unwrap().collect();
            let expected: Vec<_> = (first..=last)
                .map(|pos| {
                    let expected = store.encode_kmer_at(0, pos, k).unwrap();
                    assert_eq!(store.sequence(0).unwrap().kmer_at(pos, k), Some(expected));
                    (pos, expected)
                })
                .collect();
            assert_eq!(got, expected, "k={k} first={first}");
            for flags in 0..4 {
                for query in [
                    expected[0].1,
                    expected[expected.len() / 2].1,
                    expected.last().unwrap().1,
                    u64::MAX,
                ] {
                    let reverse = query ^ 3;
                    let reference = expected.iter().find_map(|&(start, code)| {
                        if flags & 1 != 0 && code == query {
                            Some((start, false))
                        } else if flags & 2 != 0 && code == reverse {
                            Some((start, true))
                        } else {
                            None
                        }
                    });
                    assert_eq!(
                        store.kmer_range(0, first, last, k).unwrap().first_match(
                            query,
                            reverse,
                            flags & 1 != 0,
                            flags & 2 != 0
                        ),
                        reference,
                        "k={k} first={first} flags={flags}"
                    );
                }
            }
        }
    }
    assert!(store.kmer_range(0, 0, 64, 1).is_none());
    assert!(store.kmer_range(0, 4, 3, 31).is_none());
    assert!(store.kmer_range(0, 180, 190, 31).is_none());
    assert!(store.kmer_range(0, usize::MAX, usize::MAX, 31).is_none());
}

#[test]
fn packed_metadata_preserves_wide_boundaries() {
    let compact = SequenceMetadata::new(vec![0, 2], vec![63]);
    assert!(matches!(compact, SequenceMetadata::Compact(_)));
    assert_eq!(compact.get(0), Some((0, 63)));
    assert_eq!(compact.get(1), None);
    let mut unsupported = EncodedBaseStore {
        metadata: SequenceMetadata::Compact(vec![(1u64 << 31) | (32u64 << 32), 31u64 << 32]),
        values: Vec::new(),
    };
    unsupported.inline_short_words();
    assert!(matches!(unsupported.metadata, SequenceMetadata::Compact(_)));
    assert_eq!(unsupported.metadata.get(0), Some((1usize << 31, 32)));
    let Some(big) = (u32::MAX as usize).checked_add(1) else {
        return;
    };
    let wide = SequenceMetadata::new(vec![0, big, big + 1], vec![big, 1]);
    assert!(matches!(wide, SequenceMetadata::Wide(_)));
    assert_eq!(wide.get(1), Some((big, 1)));
}

#[test]
fn packed_minimizers_equal_byte_reference_for_all_g_and_word_boundaries() {
    let mut state = 43u64;
    let sequences: Vec<Vec<u8>> = (0..=129)
        .map(|len| {
            (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGTacgt"[(state & 7) as usize]
                })
                .collect()
        })
        .chain([vec![b'A'; 97], vec![b'T'; 97], b"ACGT".repeat(32)])
        .collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    for store in [
        EncodedBaseStore::from_sequences(&sequences, None),
        EncodedBaseStore::from_sequences(&sequences, Some(&pool)),
    ] {
        for (seq_id, sequence) in sequences.iter().enumerate() {
            for g in 1..=32 {
                let mut expected = vec![[0; 8]; minimizer_count(sequence.len(), g)];
                compute_minimizer_reps_into(sequence, g, &mut expected);
                for (pos, bytes) in expected.iter().enumerate() {
                    assert_eq!(
                        store.minimizer_at(seq_id, pos, g).as_ref(),
                        Some(bytes),
                        "sequence {seq_id}, position {pos}, g {g}"
                    );
                }
                assert_eq!(store.minimizer_at(seq_id, expected.len(), g), None);
                assert_eq!(store.minimizer_at(seq_id, usize::MAX, g), None);
            }
            assert_eq!(store.minimizer_at(seq_id, 0, 0), None);
            assert_eq!(store.minimizer_at(seq_id, 0, 33), None);
        }
        assert_eq!(store.minimizer_at(sequences.len(), 0, 1), None);
    }
}

#[test]
fn packed_windows_equal_original_byte_encoding_at_every_boundary() {
    let sequences = vec![
        Vec::new(),
        b"A".to_vec(),
        b"ACGTNacgtnX".to_vec(),
        vec![b'A'; 31],
        vec![b'T'; 31],
        (0..31).map(|i| b"ACGTNacgtnX"[i % 10]).collect(),
        (0..129).map(|i| b"ACGTacgtN"[i % 9]).collect(),
    ];
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .expect("thread pool");
    let mut inline = EncodedBaseStore::from_sequences(&sequences, None);
    inline.inline_short_words();
    assert!(matches!(inline.metadata, SequenceMetadata::InlineShort(_)));
    inline.inline_short_words(); // Repeated preparation is idempotent.
    for store in [
        EncodedBaseStore::from_sequences(&sequences, None),
        EncodedBaseStore::from_sequences(&sequences, Some(&pool)),
        inline,
    ] {
        for (seq_id, sequence) in sequences.iter().enumerate() {
            assert_eq!(store.len(seq_id), Some(sequence.len()));
            for start in 0..=sequence.len() + 1 {
                for k in 0..=33 {
                    let expected = if k <= 32 && start + k <= sequence.len() {
                        Some(sequence[start..start + k].iter().fold(0u64, |word, &base| {
                            (word << 2) | u64::from(BASE_CODES[base as usize] & 3)
                        }))
                    } else {
                        None
                    };
                    assert_eq!(store.encode_kmer_at(seq_id, start, k), expected);
                    assert_eq!(store.sequence(seq_id).unwrap().kmer_at(start, k), expected);
                    if k > 0
                        && let Some(expected) = expected
                    {
                        assert_eq!(
                            store.kmer_range(seq_id, start, start, k).unwrap().next(),
                            Some((start, expected))
                        );
                    }
                }
            }
            assert_eq!(store.encode_kmer_at(seq_id, usize::MAX, 1), None);
        }
        assert_eq!(store.len(sequences.len()), None);
        assert_eq!(store.encode_kmer_at(sequences.len(), 0, 0), None);
    }
}
