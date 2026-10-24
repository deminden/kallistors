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
    let mphf = BooPhf::load(&mut reader).unwrap();
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
    let mut fallback = mphf.clone();
    fallback.compact_levels = false;
    for level in &mut fallback.levels {
        level.bitset.ranked_words.clear();
    }
    for batch in 0..2500 {
        let keys = std::array::from_fn(|lane| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (lane < batch % 9).then_some(state.to_le_bytes())
        });
        let expected = keys.map(|key| key.and_then(|key| mphf.debug_lookup(&key).0));
        assert_eq!(mphf.lookup_batch8(&keys), expected);
        assert_eq!(compact.lookup_batch8(&keys), expected);
        assert_eq!(fallback.lookup_batch8(&keys), expected);
        for key in keys.into_iter().flatten() {
            assert_eq!(mphf.lookup(&key), mphf.debug_lookup(&key).0);
            assert_eq!(compact.lookup(&key), mphf.debug_lookup(&key).0);
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
        tail.lookup_batch8(&keys),
        [Some(19), None, None, None, Some(19), None, None, None]
    );
}
