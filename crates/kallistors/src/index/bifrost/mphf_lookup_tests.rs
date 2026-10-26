#[test]
fn absolute_rank_words_preserve_padded_domain_coverage() {
    // Only the hash domain must have compact coverage. Serialized padding may
    // extend beyond that coverage without changing any valid query result.
    let mut bytes = Vec::new();
    for value in [64u64, 9] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    for _ in 0..9 {
        bytes.extend_from_slice(&5u64.to_le_bytes());
    }
    for value in [1u64, 17] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    let mut bits = BitVector::load(&mut std::io::Cursor::new(bytes)).unwrap();
    assert_eq!(bits.ranked_words.len(), 8);
    assert_eq!(bits._nchar, 9);
    let expected: Vec<_> = (0..64)
        .map(|pos| (bits.get(pos), bits.rank(pos), bits.get_rank(pos)))
        .collect();
    bits.promote_absolute_words();
    bits.bit_array.clear();
    bits.ranks.clear();
    for (pos, &(bit, rank, found)) in expected.iter().enumerate() {
        assert_eq!(bits.get(pos as u64), bit);
        assert_eq!(bits.rank(pos as u64), rank);
        assert_eq!(bits.get_rank(pos as u64), found);
    }
}

use super::*;

#[test]
fn bulk_bitvector_words_preserve_little_endian_and_truncated_input_errors() {
    let words = [0u64, 1, u64::MAX, 1 << 63, 0x0123456789abcdef];
    let bytes: Vec<_> = words.iter().flat_map(|word| word.to_le_bytes()).collect();
    assert_eq!(
        read_u64_array(&mut std::io::Cursor::new(&bytes), words.len()).unwrap(),
        words
    );
    assert!(
        read_u64_array(&mut std::io::Cursor::new([]), 0)
            .unwrap()
            .is_empty()
    );
    for end in 0..bytes.len() {
        assert!(read_u64_array(&mut std::io::Cursor::new(&bytes[..end]), words.len()).is_err());
    }
}

#[test]
fn reduced_domains_preserve_exact_u64_remainders_at_boundaries() {
    let mut state = 17u64;
    for domain in [
        1,
        2,
        3,
        7,
        64,
        65,
        4096,
        4160,
        8128,
        8192,
        8256,
        1 << 20,
        50_102_144,
        50_102_171,
        u32::MAX as u64,
        1u64 << 63,
        u64::MAX - 63,
        u64::MAX,
    ] {
        let reduced = StrengthReducedU64::new(domain);
        let reciprocal = ReducedDomain::new(domain);
        for value in [0, 1, domain - 1, domain, domain.saturating_add(1), u64::MAX] {
            assert_eq!(value % reduced, value % domain);
            assert_eq!(reciprocal.remainder(value), value % domain);
        }
        for _ in 0..20_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            assert_eq!(state % reduced, state % domain);
            assert_eq!(reciprocal.remainder(state), state % domain);
        }
    }
}

#[test]
fn lookup_matches_diagnostic_oracle_across_hash_levels() {
    assert_eq!(std::mem::size_of::<RankedBitWord>(), 10);
    assert_eq!(std::mem::size_of::<AbsoluteRankedBitWord>(), 16);
    assert_eq!(std::mem::align_of::<AbsoluteRankedBitWord>(), 16);
    // Independent full-width rank oracle at both1024-word base transitions,
    // including the largest possible16-bit relative rank.
    let mut words = vec![u64::MAX; 2051];
    words[1023] = 0;
    words[1024] = 5;
    words[2048] = 1;
    let mut rank = 17u64;
    let ranks: Vec<_> = words
        .chunks(8)
        .map(|chunk| {
            let base = rank;
            rank += chunk.iter().map(|v| u64::from(v.count_ones())).sum::<u64>();
            base
        })
        .collect();
    let mut bytes = Vec::new();
    for value in [words.len() as u64 * 64, words.len() as u64] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    for value in &words {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(&(ranks.len() as u64).to_le_bytes());
    for value in &ranks {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    let mut compact_words = BitVector::load(&mut std::io::Cursor::new(bytes)).unwrap();
    let probes: Vec<_> = [0usize, 1, 1022, 1023, 1024, 1025, 2047, 2048, 2050]
        .into_iter()
        .flat_map(|word| [0usize, 1, 17, 63].map(move |offset| (word * 64 + offset) as u64))
        .collect();
    let expected: Vec<_> = probes
        .iter()
        .map(|&pos| compact_words.get(pos).then(|| compact_words.rank(pos)))
        .collect();
    assert_eq!(compact_words.rank_bases.len(), 3);
    compact_words.bit_array.clear();
    compact_words.ranks.clear();
    for (&pos, expected) in probes.iter().zip(expected) {
        assert_eq!(compact_words.get_rank(pos), expected);
    }
    for base in [u64::from(u32::MAX) - 1, u64::from(u32::MAX) + 1] {
        let bytes: Vec<u8> = [64u64, 1, 0b101, 1, base]
            .into_iter()
            .flat_map(u64::to_le_bytes)
            .collect();
        let bitset = BitVector::load(&mut std::io::Cursor::new(bytes)).unwrap();
        for pos in 0..64 {
            assert_eq!(
                bitset.get_rank(pos),
                bitset.get(pos).then(|| bitset.rank(pos))
            );
        }
    }

    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("transcripts.fa");
    let mut state = 137u64;
    let mut text = String::from(">tx\n");
    for _ in 0..5000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        text.push(b"ACGT"[(state & 3) as usize] as char);
    }
    text.push('\n');
    std::fs::write(&fasta, text).unwrap();
    let path = temp.path().join("index.idx");
    crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default()).unwrap();
    let data = std::fs::read(path).unwrap();
    let mut reader = std::io::Cursor::new(data.as_slice());
    let _version = read_u64_le(&mut reader).unwrap();
    let graph_bytes = read_u64_le(&mut reader).unwrap();
    reader.set_position(16 + graph_bytes);
    let _mphf_bytes = read_u64_le(&mut reader).unwrap();
    let mphf_start = reader.position() as usize;
    let mphf = BooPhf::load(&mut reader).unwrap();
    // Select the large-table load layout without changing the serialized rank
    // oracle, hash levels, keys or expected values of this compact fixture.
    let mut large_bytes = data[mphf_start..].to_vec();
    large_bytes[20..28].copy_from_slice(&(1u64 << 28).to_le_bytes());
    let direct = BooPhf::load(&mut std::io::Cursor::new(large_bytes)).unwrap();
    assert!(direct.absolute_levels);

    let mut compact = mphf.clone();
    compact.release_redundant_rank_storage();
    for (old, new) in mphf.levels.iter().zip(&compact.levels) {
        for pos in 0..old.hash_domain {
            assert_eq!(new.bitset.get_rank(pos), old.bitset.get_rank(pos));
            assert_eq!(new.bitset.get(pos), old.bitset.get(pos));
            assert_eq!(new.bitset.rank(pos), old.bitset.rank(pos));
        }
    }
    for level in &mphf.levels {
        for pos in 0..level.hash_domain {
            assert_eq!(
                level.bitset.get_rank(pos),
                level.bitset.get(pos).then(|| level.bitset.rank(pos))
            );
        }
    }
    assert!(mphf.compact_levels);
    let mut absolute = mphf.clone();
    absolute.absolute_levels = true;
    for level in &mut absolute.levels {
        level.bitset.promote_absolute_words();
    }
    absolute.release_redundant_rank_storage();
    for (old, new) in mphf.levels.iter().zip(&absolute.levels) {
        for pos in 0..old.hash_domain {
            assert_eq!(new.bitset.get_rank(pos), old.bitset.get_rank(pos));
            assert_eq!(new.bitset.get(pos), old.bitset.get(pos));
            assert_eq!(new.bitset.rank(pos), old.bitset.rank(pos));
        }
    }
    let mut fallback = mphf.clone();
    fallback.compact_levels = false;
    for level in &mut fallback.levels {
        level.bitset.ranked_words.clear();
    }
    for batch in 0..2500 {
        let keys: [Option<[u8; 8]>; 8] = std::array::from_fn(|lane| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (lane < batch % 9).then_some(state.to_le_bytes())
        });
        let expected = keys.map(|key| key.and_then(|key| mphf.debug_lookup(&key).0));
        assert_eq!(
            mphf.lookup_batch(&MphfBatchKeys::from(&keys))
                .into_options(),
            expected
        );
        assert_eq!(
            compact
                .lookup_batch(&MphfBatchKeys::from(&keys))
                .into_options(),
            expected
        );
        assert_eq!(
            absolute
                .lookup_batch(&MphfBatchKeys::from(&keys))
                .into_options(),
            expected
        );
        assert_eq!(
            direct
                .lookup_batch(&MphfBatchKeys::from(&keys))
                .into_options(),
            expected
        );
        assert_eq!(
            fallback
                .lookup_batch(&MphfBatchKeys::from(&keys))
                .into_options(),
            expected
        );
        let wide_keys: [Option<[u8; 8]>; 16] = std::array::from_fn(|lane| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (lane < batch % 17).then_some(state.to_le_bytes())
        });
        let wide_expected = wide_keys.map(|key| key.and_then(|key| mphf.debug_lookup(&key).0));
        assert_eq!(
            mphf.lookup_batch(&MphfBatchKeys::from(&wide_keys))
                .into_options(),
            wide_expected
        );
        assert_eq!(
            compact
                .lookup_batch(&MphfBatchKeys::from(&wide_keys))
                .into_options(),
            wide_expected
        );
        assert_eq!(
            absolute
                .lookup_batch(&MphfBatchKeys::from(&wide_keys))
                .into_options(),
            wide_expected
        );
        assert_eq!(
            direct
                .lookup_batch(&MphfBatchKeys::from(&wide_keys))
                .into_options(),
            wide_expected
        );
        assert_eq!(
            fallback
                .lookup_batch(&MphfBatchKeys::from(&wide_keys))
                .into_options(),
            wide_expected
        );
        let loader_keys: [Option<[u8; 8]>; 32] = std::array::from_fn(|lane| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (lane < batch % 33).then_some(state.to_le_bytes())
        });
        let loader_expected = loader_keys.map(|key| key.and_then(|key| mphf.debug_lookup(&key).0));
        assert_eq!(
            mphf.lookup_batch(&MphfBatchKeys::from(&loader_keys))
                .into_options(),
            loader_expected
        );
        assert_eq!(
            compact
                .lookup_batch(&MphfBatchKeys::from(&loader_keys))
                .into_options(),
            loader_expected
        );
        assert_eq!(
            absolute
                .lookup_batch(&MphfBatchKeys::from(&loader_keys))
                .into_options(),
            loader_expected
        );
        assert_eq!(
            direct
                .lookup_batch(&MphfBatchKeys::from(&loader_keys))
                .into_options(),
            loader_expected
        );
        assert_eq!(
            fallback
                .lookup_batch(&MphfBatchKeys::from(&loader_keys))
                .into_options(),
            loader_expected
        );
        let widest_keys: [Option<[u8; 8]>; 64] = std::array::from_fn(|lane| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (lane < batch % 65 && (batch % 2 == 0 || lane % 3 != 1)).then_some(state.to_le_bytes())
        });
        let widest_expected = widest_keys.map(|key| key.and_then(|key| mphf.debug_lookup(&key).0));
        assert_eq!(
            mphf.lookup_batch(&MphfBatchKeys::from(&widest_keys))
                .into_options(),
            widest_expected
        );
        assert_eq!(
            compact
                .lookup_batch(&MphfBatchKeys::from(&widest_keys))
                .into_options(),
            widest_expected
        );
        assert_eq!(
            absolute
                .lookup_batch(&MphfBatchKeys::from(&widest_keys))
                .into_options(),
            widest_expected
        );
        assert_eq!(
            fallback
                .lookup_batch(&MphfBatchKeys::from(&widest_keys))
                .into_options(),
            widest_expected
        );
        assert_eq!(
            direct
                .lookup_batch(&MphfBatchKeys::from(&widest_keys))
                .into_options(),
            widest_expected
        );
        for key in keys.into_iter().flatten() {
            assert_eq!(mphf.lookup(&key), mphf.debug_lookup(&key).0);
            assert_eq!(compact.lookup(&key), mphf.debug_lookup(&key).0);
            assert_eq!(absolute.lookup(&key), mphf.debug_lookup(&key).0);
            assert_eq!(fallback.lookup(&key), mphf.debug_lookup(&key).0);
        }
    }
    let key = 1729u64.to_le_bytes();
    let mut tail = mphf;
    tail.nb_levels = 1;
    tail.levels.clear();
    tail.final_hash = HashMap::from([(key, 12)]);
    tail.lastbitsetrank = 7;
    let keys = [
        Some(key),
        None,
        Some([0; 8]),
        None,
        Some(key),
        None,
        None,
        None,
    ];
    assert_eq!(
        tail.lookup_batch(&MphfBatchKeys::from(&keys))
            .into_options(),
        [Some(19), None, None, None, Some(19), None, None, None]
    );
}

#[test]
fn direct_absolute_rank_construction_preserves_legacy_fallback_and_padding() {
    for (base, gap, rank_count) in [
        (0u64, 0u64, 3usize),
        (17, 0, 2),
        (u64::from(u32::MAX) - 2, 0, 3),
        (u64::from(u32::MAX) + 1, 0, 3),
        (100, 100_000, 3),
    ] {
        let words: Vec<u64> = (0..17).map(|i| 0x3197a5c73u64.rotate_left(i * 3)).collect();
        let mut rank = base;
        let ranks: Vec<_> = words
            .chunks(8)
            .enumerate()
            .take(rank_count)
            .map(|(i, chunk)| {
                if i != 0 {
                    rank += gap;
                }
                let first = rank;
                rank += chunk.iter().map(|v| u64::from(v.count_ones())).sum::<u64>();
                first
            })
            .collect();
        let mut bytes = Vec::new();
        for word in [12u64 * 64, words.len() as u64]
            .into_iter()
            .chain(words)
            .chain(std::iter::once(ranks.len() as u64))
            .chain(ranks)
        {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        let original = BitVector::load(&mut std::io::Cursor::new(&bytes)).unwrap();
        let direct = BitVector::load_with_layout(&mut std::io::Cursor::new(&bytes), true).unwrap();
        for position in 0..12u64 * 64 {
            assert_eq!(direct.get(position), original.get(position));
            assert_eq!(direct.rank(position), original.rank(position));
            assert_eq!(direct.get_rank(position), original.get_rank(position));
        }
        for end in [0, 7, 15, bytes.len() - 1] {
            assert!(
                BitVector::load_with_layout(&mut std::io::Cursor::new(&bytes[..end]), true)
                    .is_err()
            );
        }
    }
}

#[test]
fn masked_lookup_results_preserve_full_width_ranks_and_missing_lanes() {
    let key = 1729u64.to_le_bytes();
    let phf = BooPhf {
        _gamma: 2.0,
        nb_levels: 1,
        lastbitsetrank: 0,
        nelem: u64::MAX,
        levels: Vec::new(),
        final_hash: HashMap::from([(key, u64::MAX)]),
        compact_levels: false,
        absolute_levels: false,
    };
    let mut keys = [None; 32];
    keys[0] = Some(key);
    keys[31] = Some(key);
    keys[17] = Some([0; 8]);
    let expected = keys.map(|key| key.and_then(|value| phf.debug_lookup(&value).0));
    assert_eq!(
        phf.lookup_batch(&MphfBatchKeys::from(&keys)).into_options(),
        expected
    );
    assert_eq!(
        phf.lookup_batch(&MphfBatchKeys::from(&keys))
            .into_iter()
            .collect::<Vec<_>>(),
        expected.to_vec()
    );
    // The borrowed prefetch view must not truncate high ranks, treat0 as
    // absent, or include absent lanes; it must not consume the result either.
    for found in [0, 1, 1 << 63, 0xaaaa_aaaa_aaaa_aaaa, u64::MAX] {
        let ranks =
            std::array::from_fn::<_, 64, _>(
                |lane| if lane == 0 { 0 } else { u64::MAX - lane as u64 },
            );
        let expected: Vec<_> = ranks
            .iter()
            .enumerate()
            .filter_map(|(lane, &rank)| (found & (1u64 << lane) != 0).then_some(rank))
            .collect();
        let result = BatchLookupResult { ranks, found };
        assert_eq!(result.found_ranks().collect::<Vec<_>>(), expected);
        assert_eq!(result.into_iter().flatten().collect::<Vec<_>>(), expected);
    }
}

#[test]
fn masked_lookup_keys_preserve_zero_full_width_sparse_tails_and_reuse() {
    let mut expected = [None; 64];
    let mut keys = MphfBatchKeys::<64>::new();
    for round in 0..67 {
        keys.clear();
        expected.fill(None);
        for (lane, entry) in expected.iter_mut().enumerate().take(round.min(64)) {
            let word = if lane % 3 == 0 {
                0
            } else {
                u64::MAX - lane as u64
            };
            *entry = Some(word.to_le_bytes());
            keys.set(lane, *entry);
        }
        assert_eq!(keys.iter().collect::<Vec<_>>(), expected.to_vec());
        keys.set(63, Some(u64::MAX.to_le_bytes()));
        expected[63] = Some(u64::MAX.to_le_bytes());
        keys.set(0, Some(0u64.to_le_bytes()));
        expected[0] = Some(0u64.to_le_bytes());
        assert_eq!(keys.iter().collect::<Vec<_>>(), expected.to_vec());
        keys.set(0, None);
        expected[0] = None;
        assert_eq!(keys.iter().collect::<Vec<_>>(), expected.to_vec());
    }
    keys.clear();
    assert_eq!(keys.iter().collect::<Vec<_>>(), vec![None; 64]);
    assert_eq!(
        std::mem::align_of::<MphfBatchKeys<32>>(),
        std::mem::align_of::<u64>()
    );
    let phf = BooPhf {
        _gamma: 2.0,
        nb_levels: 1,
        lastbitsetrank: 0,
        nelem: u64::MAX,
        levels: Vec::new(),
        final_hash: HashMap::from([(0u64.to_le_bytes(), 1), (u64::MAX.to_le_bytes(), u64::MAX)]),
        compact_levels: false,
        absolute_levels: false,
    };
    keys.set(0, Some(0u64.to_le_bytes()));
    keys.set(63, Some(u64::MAX.to_le_bytes()));
    let mut wanted = [None; 64];
    wanted[0] = Some(1u64);
    wanted[63] = Some(u64::MAX);
    assert_eq!(phf.lookup_batch(&keys).into_options(), wanted);
}
