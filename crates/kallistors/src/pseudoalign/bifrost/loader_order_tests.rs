use super::*;

#[test]
fn buffered_ec_nodes_preserve_boundaries_and_truncation_errors() {
    let mut bitmap = vec![1u8];
    bitmap.extend_from_slice(&2u32.to_le_bytes());
    bitmap.extend_from_slice(&1u32.to_le_bytes());
    bitmap.extend_from_slice(&7u32.to_le_bytes());
    let mut node = 7u32.to_le_bytes().to_vec();
    node.push(1);
    node.extend_from_slice(&0u32.to_le_bytes());
    node.extend_from_slice(&31u32.to_le_bytes());
    node.extend_from_slice(&(bitmap.len() as u64).to_le_bytes());
    node.extend_from_slice(&bitmap);
    node.extend_from_slice(&0u64.to_le_bytes());
    for capacity in [1, 5, 17, 1024] {
        let mut payload = node.clone();
        payload.push(0x5a);
        let mut reader = std::io::BufReader::with_capacity(capacity, std::io::Cursor::new(payload));
        let mut flat = FlatEcIndex::new(1);
        flat.push_buffered_node(&mut reader, node.len(), &mut Vec::new())
            .unwrap();
        assert_eq!(flat.ec(0, 0), &[1, 7]);
        assert_eq!(flat.block_bounds(0, 0), Some((0, 31)));
        let mut next = [0];
        reader.read_exact(&mut next).unwrap();
        assert_eq!(next, [0x5a]);
    }
    for cut in 0..node.len() {
        let mut reader = std::io::BufReader::with_capacity(7, std::io::Cursor::new(&node[..cut]));
        let mut flat = FlatEcIndex::new(1);
        assert!(
            flat.push_buffered_node(&mut reader, node.len(), &mut Vec::new())
                .is_err()
        );
        assert!(flat.ec_values.is_empty());
    }
}

#[test]
fn compressed_words_preserve_every_byte_and_partial_word() {
    for len in 0usize..=257 {
        let data: Vec<u8> = (0..len.div_ceil(4)).map(|i| (i * 79 + 37) as u8).collect();
        let mut values = Vec::new();
        append_compressed_words(&data, len, &mut values);
        for capacity in [1, 7, 8, 9, 64, 1024] {
            let mut payload = data.clone();
            payload.push(0x5a);
            let mut reader =
                std::io::BufReader::with_capacity(capacity, std::io::Cursor::new(payload));
            let mut actual = vec![u64::MAX];
            let mut scratch = Vec::new();
            append_buffered_compressed(&mut reader, len, &mut actual, &mut scratch).unwrap();
            assert_eq!(&actual[1..], values.as_slice());
            assert_eq!(actual[0], u64::MAX);
            let mut following = [0];
            reader.read_exact(&mut following).unwrap();
            assert_eq!(following, [0x5a]);
        }
        for cut in 0..data.len() {
            let mut reader =
                std::io::BufReader::with_capacity(7, std::io::Cursor::new(&data[..cut]));
            let mut actual = vec![u64::MAX];
            assert!(
                append_buffered_compressed(&mut reader, len, &mut actual, &mut Vec::new()).is_err()
            );
            assert_eq!(actual, vec![u64::MAX]);
        }
        let seq = EncodedSequence {
            words: &values,
            len,
        };
        let decoded = decode_compressed_sequence(&data, len);
        assert_eq!(values.len(), len.div_ceil(32));
        for (i, &base) in decoded.iter().enumerate() {
            assert_eq!(b"ACGT"[seq.kmer_at(i, 1).unwrap() as usize], base);
        }
    }
}

#[test]
fn abundant_nodes_use_node_keys_after_hash_table_reordering() {
    let first = b"ACCTCAGCCTCCCAAAATGCTAGGATTACAG";
    let second = b"CGGGTGTGGTGGCGGGCACCTGTAGTCCCAG";
    let canonical = |seq: &[u8]| {
        let (forward, reverse) = encode_kmer_pair(seq).unwrap();
        forward.min(reverse)
    };
    let base = 19;
    // The graph's hash-table serialization has the opposite order to
    // the node records, each of which supplies its own EC and head k-mer.
    let mut kmers = vec![second.to_vec(), first.to_vec()];
    let mut map = HashMap::from([(canonical(second), base), (canonical(first), base + 1)]);
    associate_abundant_node(&mut kmers, &mut map, base, base, first).unwrap();
    associate_abundant_node(&mut kmers, &mut map, base, base + 1, second).unwrap();
    assert_eq!(kmers, [first.to_vec(), second.to_vec()]);
    assert_eq!(map[&canonical(first)], base);
    assert_eq!(map[&canonical(second)], base + 1);
    let mut words: Vec<_> = [second, first]
        .iter()
        .map(|seq| encode_kmer_pair(*seq).unwrap().0 << 2)
        .collect();
    let mut word_map = HashMap::from([(canonical(second), base), (canonical(first), base + 1)]);
    associate_abundant_word(&mut words, &mut word_map, base, base, first).unwrap();
    associate_abundant_word(&mut words, &mut word_map, base, base + 1, second).unwrap();
    assert_eq!(
        words,
        [
            encode_kmer_pair(first).unwrap().0 << 2,
            encode_kmer_pair(second).unwrap().0 << 2
        ]
    );
    assert_eq!(word_map, map);
    let fallback = SelectiveFallbackKmerIndex::build(base, &[], &kmers);
    assert_eq!(fallback.lookup(canonical(first)).unwrap()[0].0, base);
    assert_eq!(fallback.lookup(canonical(second)).unwrap()[0].0, base + 1);
    let duplicates = SelectiveFallbackKmerIndex::build(
        base,
        &[first.to_vec(), second.to_vec()],
        &[first.to_vec(), first.to_vec()],
    );
    let forward = encode_kmer_pair(first).unwrap().0;
    assert_eq!(
        duplicates.lookup(canonical(first)).unwrap(),
        &[(base, forward), (base + 2, forward), (base + 3, forward)]
    );
    assert_eq!(duplicates.lookup(canonical(second)).unwrap()[0].0, base + 1);
    assert!(
        duplicates
            .lookup(canonical(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))
            .is_none()
    );

    let km = vec![first.to_vec(), second.to_vec()];
    let h = [first.to_vec(), first.to_vec()];
    let packed_km = EncodedBaseStore::from_sequences(&km, None);
    let mut h_words: Vec<u64> = h
        .iter()
        .map(|seq| encode_kmer_pair(seq).unwrap().0 << 2)
        .collect();
    h_words.reverse();
    let mut h_map = HashMap::from([(canonical(first), base + 3)]);
    associate_abundant_word(&mut h_words, &mut h_map, base + 2, base + 2, first).unwrap();
    associate_abundant_word(&mut h_words, &mut h_map, base + 2, base + 3, first).unwrap();
    assert_eq!(h_map[&canonical(first)], base + 3);
    let packed = SelectiveFallbackKmerIndex::build_packed(base, &packed_km, &h_words, 31);
    for key in [canonical(first), canonical(second)] {
        assert_eq!(packed.lookup(key), duplicates.lookup(key));
    }
    let words: Vec<u64> = km
        .iter()
        .map(|seq| encode_kmer_pair(seq).unwrap().0 << 2)
        .collect();
    let fixed = EncodedBaseStore::from_fixed_words(words, 31);
    for (id, seq) in km.iter().enumerate() {
        assert_eq!(
            decoded_graph_sequence(&[], &fixed, id).unwrap().as_ref(),
            seq
        );
    }
    assert!(
        associate_abundant_node(
            &mut kmers,
            &mut map,
            base,
            base,
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        )
        .is_err()
    );
}

#[test]
fn singleton_bucket_storage_preserves_zero_specials_and_tie_order() {
    if usize::BITS != 64 {
        return;
    }
    let mut index = MinzPositionIndex::from_counts_inline(vec![0, 1, 2, 1, 0], &[3], true);
    assert!(index.inline.is_some());
    let mut cursors = index.cursors();
    assert!(cursors.is_empty());
    index.push_at(&mut cursors, 1, 0);
    index.push_at(&mut cursors, 2, 10);
    index.push_at(&mut cursors, 3, u64::MAX);
    index.push_at(&mut cursors, 2, 20);
    index.finish_inline();
    assert_eq!(index.get(0), &[]);
    assert_eq!(index.get(1), &[0]);
    assert_eq!(index.get(2), &[10, 20]);
    assert_eq!(index.get(3), &[u64::MAX]);
    assert_eq!(index.get(4), &[]);
    for bucket in 0..5 {
        assert_eq!(index.get_with_overcrowded(bucket).1, bucket == 3);
    }
}

#[test]
fn descriptor_marker_bit_preserves_specials_and_high_singleton_bits() {
    if usize::BITS != 64 {
        return;
    }
    let counts = vec![1, 3, 1, 1];
    let mut inline = MinzPositionIndex::from_counts_inline(counts.clone(), &[2, 3], true);
    let mut csr = MinzPositionIndex::from_counts(counts);
    let mut inline_cursors = inline.cursors();
    let mut csr_cursors = csr.cursors();
    for (bucket, value) in [
        (0, (1u64 << 62) | (1u64 << 61) | 123),
        (1, 11),
        (1, u64::MAX),
        (1, 22),
        (2, 0xffff_ffff_0000_0001),
        (3, 0xffff_ffff_8000_0000),
    ] {
        inline.push_at(&mut inline_cursors, bucket, value);
        csr.push_at(&mut csr_cursors, bucket, value);
    }
    inline.finish_inline();
    for bucket in 0..4 {
        assert_eq!(
            inline.get_with_special(bucket),
            csr.get_with_special(bucket)
        );
        assert_eq!(inline.get_with_special(bucket).1, bucket != 0);
        assert_eq!(
            inline.get_with_overcrowded(bucket),
            csr.get_with_overcrowded(bucket)
        );
        assert_eq!(
            inline.get_with_overcrowded(bucket).1,
            bucket == 1 || bucket == 3
        );
    }
}

#[test]
fn ec_metadata_promotes_without_narrowing_large_ranges() {
    let mut blocks = FlatEcBlocks::new();
    blocks.push(1, 3, 0, 2);
    for (lb, ub, start, end) in [
        (4095, 4095, 0, 8191),
        (4096, 5000, 10, 12),
        (0, 4, 1 << 26, (1 << 26) + 1),
        (0, 2, 0, 8192),
    ] {
        let id = blocks.len();
        blocks.push(lb, ub, start, end);
        let b = blocks.get(id).unwrap();
        assert_eq!((b.lb, b.ub, b.ec_start, b.ec_end), (lb, ub, start, end));
    }
    let Some(high) = (u32::MAX as usize).checked_add(1) else {
        return;
    };
    blocks.push(4, 9, high, high + 3);
    assert!(matches!(blocks, FlatEcBlocks::Wide(_)));
    assert_eq!(blocks.get(0).unwrap().ec_end, 2);
    assert_eq!(blocks.get(blocks.len() - 1).unwrap().ec_start, high);
    let mut offsets = MinzOffsets::Compact(vec![0, 2]);
    offsets.push(high);
    assert!(matches!(offsets, MinzOffsets::Wide(_)));
    assert_eq!(offsets.value_at(1), Some(2));
    assert_eq!(offsets.value_at(2), Some(high));
}

#[test]
fn grouped_sequence_decode_equals_byte_oracle_for_all_tail_lengths() {
    let bytes: Vec<u8> = (0..=255).collect();
    for len in 0..=bytes.len() * 4 {
        let expected: Vec<u8> = (0..len)
            .map(|i| b"ACGT"[((bytes[i / 4] >> (2 * (i % 4))) & 3) as usize])
            .collect();
        assert_eq!(decode_compressed_sequence(&bytes, len), expected);
    }
    let mut word = 91u64;
    for _ in 0..1024 {
        word ^= word << 13;
        word ^= word >> 7;
        word ^= word << 17;
        for k in 0..=32 {
            let expected: Vec<u8> = (0..k)
                .map(|i| b"ACGT"[((word >> (62 - 2 * i)) & 3) as usize])
                .collect();
            assert_eq!(decode_kmer(&word.to_le_bytes(), k), expected);
        }
    }
}

#[test]
fn compact_prefixes_and_wide_count_buffer_reuse_preserve_offsets() {
    for counts in [vec![], vec![0], vec![2, 0, 3, 1]] {
        let mut reserved = Vec::with_capacity(counts.len() + 1);
        reserved.extend_from_slice(&counts);
        let index = MinzPositionIndex::from_counts(reserved);
        assert!(matches!(index.offsets, MinzOffsets::Compact(_)));
        assert_eq!(index.len(), counts.len());
        assert_eq!(index.values.len(), counts.iter().sum());
        for (bucket, &count) in counts.iter().enumerate() {
            assert_eq!(index.get(bucket).len(), count);
        }
    }
    if let Some(total) = (u32::MAX as usize).checked_add(1) {
        let mut counts = Vec::with_capacity(3);
        counts.extend_from_slice(&[u32::MAX as usize, 1]);
        let pointer = counts.as_ptr();
        let (offsets, sum) = MinzOffsets::from_counts(counts);
        assert_eq!(sum, total);
        assert_eq!(offsets.bounds(0), (0, u32::MAX as usize));
        assert_eq!(offsets.bounds(1), (u32::MAX as usize, total));
        match offsets {
            MinzOffsets::Wide(v) => assert_eq!(v.as_ptr(), pointer),
            MinzOffsets::Compact(_) => panic!("large offsets must retain usize width"),
        }
    }
}

#[test]
fn parallel_fill_preserves_every_bucket_in_serial_order() {
    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("transcripts.fa");
    let mut state = 19u64;
    let mut text = String::new();
    for transcript in 0..5 {
        text.push_str(&format!(">tx{transcript}\n"));
        for _ in 0..2000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            text.push(b"ACGT"[(state & 3) as usize] as char);
        }
        text.push('\n');
    }
    std::fs::write(&fasta, text).unwrap();
    let path = temp.path().join("index.idx");
    crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default()).unwrap();
    let serial = build_bifrost_index_threaded(&path, 1).unwrap();
    let parallel = build_bifrost_index_threaded(&path, 4).unwrap();
    assert_eq!(
        serial.minz_positions.offsets,
        parallel.minz_positions.offsets
    );
    assert_eq!(serial.minz_positions.values, parallel.minz_positions.values);
    let compact = build_bifrost_index_quant_threaded(&path, false, 4).unwrap();
    let packed = build_bifrost_index_quant_packed_threaded(&path, 4).unwrap();
    for (id, seq) in serial.unitigs.iter().enumerate() {
        assert!(packed.unitigs.is_empty());
        assert_eq!(packed.unitig_sequence(id).unwrap().as_ref(), seq);
    }
    assert!(packed.unitig_sequence(serial.unitigs.len()).is_none());
    assert!(compact.flat_ec.owned_ec_blocks.is_empty());
    assert_eq!(serial.minz_positions.len(), compact.minz_positions.len());
    for bucket in 0..serial.minz_positions.len() {
        assert_eq!(
            packed.minz_positions.get_with_overcrowded(bucket),
            serial.minz_positions.get_with_overcrowded(bucket)
        );
        assert_eq!(
            packed.minz_positions.get_with_special(bucket),
            serial.minz_positions.get_with_special(bucket)
        );

        assert_eq!(
            serial.minz_positions.get(bucket),
            compact.minz_positions.get(bucket)
        );
        assert_eq!(
            serial.minz_positions.get_with_overcrowded(bucket),
            compact.minz_positions.get_with_overcrowded(bucket)
        );
        assert_eq!(
            serial.minz_positions.get_with_special(bucket),
            compact.minz_positions.get_with_special(bucket)
        );
    }
    for (uid, blocks) in serial.flat_ec.owned_ec_blocks.iter().enumerate() {
        for (bid, block) in blocks.iter().enumerate() {
            assert_eq!(compact.flat_ec.ec(uid, bid), block.ec.as_ref());
            assert_eq!(
                compact.flat_ec.block_bounds(uid, bid),
                Some((block.lb, block.ub))
            );
        }
        assert_eq!(compact.flat_ec.block_bounds(uid, blocks.len()), None);
    }
    assert_eq!(compact.flat_ec.block_bounds(usize::MAX, 0), None);
    let positional = build_bifrost_index_quant_threaded(&path, true, 2).unwrap();
    assert!(!positional.flat_ec.owned_ec_blocks.is_empty());
    assert!(positional.has_strand_annotations());
    let packed_positional =
        build_bifrost_index_positional_packed_threaded(&path, 2, false).unwrap();
    assert!(packed_positional.has_strand_annotations());
    assert!(packed_positional.unitigs.is_empty());
    assert!(packed_positional.km_unitigs.is_empty());
    assert_eq!(positional.unitig_count(), packed_positional.unitig_count());
    for (expected, actual) in positional
        .flat_ec
        .owned_ec_blocks
        .iter()
        .zip(&packed_positional.flat_ec.owned_ec_blocks)
    {
        assert_eq!(expected.len(), actual.len());
        for (expected, actual) in expected.iter().zip(actual.iter()) {
            assert_eq!((expected.lb, expected.ub), (actual.lb, actual.ub));
            assert_eq!(expected.ec, actual.ec);
            for (expected, actual) in expected
                .positions
                .as_ref()
                .unwrap()
                .iter()
                .zip(actual.positions.as_ref().unwrap().iter())
            {
                let mut left = Vec::new();
                let mut right = Vec::new();
                expected.for_each(|v| left.push(v));
                actual.for_each(|v| right.push(v));
                assert_eq!(left, right);
            }
        }
    }
    for seq in &serial.unitigs {
        for read in seq.windows(100).step_by(37) {
            for strand in [
                super::super::Strand::Unstranded,
                super::super::Strand::Reverse,
            ] {
                let options = super::super::PseudoalignOptions::default();
                let expected =
                    super::super::ec_for_read_bifrost(&positional, read, strand, None, options);
                let actual = super::super::ec_for_read_bifrost(
                    &packed_positional,
                    read,
                    strand,
                    None,
                    options,
                );
                assert_eq!(format!("{expected:?}"), format!("{actual:?}"));
                if let (Some(expected), Some(actual)) = (expected, actual) {
                    for fragment in [1, 73, 2000] {
                        assert_eq!(
                            super::super::filter_ec_by_fragment(
                                &positional,
                                &expected.ec,
                                expected.best_match.unwrap(),
                                fragment
                            ),
                            super::super::filter_ec_by_fragment(
                                &packed_positional,
                                &actual.ec,
                                actual.best_match.unwrap(),
                                fragment
                            )
                        );
                    }
                }
            }
            for strand in [
                super::super::Strand::Unstranded,
                super::super::Strand::Reverse,
            ] {
                let mut options = super::super::PseudoalignOptions::default();
                options.investigation.trace_fast_path = true;
                for mutate in [false, true] {
                    let mut query = read.to_vec();
                    if mutate {
                        query[40] = b'N';
                    }
                    let live =
                        super::super::ec_for_read_bifrost(&compact, &query, strand, None, options);
                    let mut debug = super::super::ReadDebugState::default();
                    let traced = super::super::ec_for_read_bifrost(
                        &compact,
                        &query,
                        strand,
                        Some(&mut debug),
                        options,
                    );
                    assert_eq!(format!("{live:?}"), format!("{traced:?}"));
                }
                for general in [false, true] {
                    let mut debug_full = super::super::ReadDebugState::default();
                    let mut debug_compact = super::super::ReadDebugState::default();
                    let full = super::super::ec_for_read_bifrost(
                        &serial,
                        read,
                        strand,
                        general.then_some(&mut debug_full),
                        super::super::PseudoalignOptions::default(),
                    );
                    let lean = super::super::ec_for_read_bifrost(
                        &compact,
                        read,
                        strand,
                        general.then_some(&mut debug_compact),
                        super::super::PseudoalignOptions::default(),
                    );
                    assert_eq!(format!("{full:?}"), format!("{lean:?}"));
                    for query in [read.to_vec(), read.to_ascii_lowercase()] {
                        let mut debug_packed = super::super::ReadDebugState::default();
                        let got = super::super::ec_for_read_bifrost(
                            &packed,
                            &query,
                            strand,
                            general.then_some(&mut debug_packed),
                            super::super::PseudoalignOptions::default(),
                        );
                        let mut debug_expected = super::super::ReadDebugState::default();
                        let expected = super::super::ec_for_read_bifrost(
                            &serial,
                            &query,
                            strand,
                            general.then_some(&mut debug_expected),
                            super::super::PseudoalignOptions::default(),
                        );
                        assert_eq!(format!("{expected:?}"), format!("{got:?}"));
                    }
                }
            }
        }
    }
}
