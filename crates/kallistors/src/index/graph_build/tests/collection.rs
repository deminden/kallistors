//! Collection regressions.

use super::*;
#[test]
fn kmer_roundtrip() {
    let seq = b"ACGTACGTAGCCATGACTACGATCGATCGTAGCAGTACGTACGATGCATCGATGCTACGTAGCATGC";
    for k in (3..=31).step_by(2) {
        let mut observed = Vec::new();
        scan_kmer_pairs(seq, k, |pos, fwd, rev| {
            observed.push((pos, fwd, rev));
            Ok(())
        })
        .unwrap();
        let expected = seq
            .windows(k)
            .enumerate()
            .map(|(pos, bases)| {
                let (fwd, rev) = encode_kmer_pair(bases).unwrap();
                (pos, fwd, rev)
            })
            .collect::<Vec<_>>();
        assert_eq!(observed, expected);
        for invalid in 0..seq.len() {
            let mut bad = seq.to_vec();
            bad[invalid] = b'N';
            let mut count = 0;
            assert!(
                scan_kmer_pairs(&bad, k, |_, _, _| {
                    count += 1;
                    Ok(())
                })
                .is_err()
            );
            assert_eq!(count, invalid.saturating_sub(k - 1));
        }
    }

    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG";
    let word = encode_kmer(seq).expect("encode");
    assert_eq!(decode_kmer_word(word, seq.len()), seq);
    let mut random = 29u64;
    for k in 1..=32 {
        for _ in 0..100 {
            random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
            let expected = (0..k)
                .map(|i| b"ACGT"[((random >> (62 - 2 * i)) & 3) as usize])
                .collect::<Vec<_>>();
            let mut actual = vec![b'X'; 5];
            append_kmer_bases(&mut actual, random, k);
            assert_eq!(&actual[..5], b"XXXXX");
            assert_eq!(&actual[5..], expected);
        }
    }
}

#[test]
fn fast_reverse_complement_matches_reference() {
    for k in (3..=31).step_by(2) {
        let mask = kmer_mask(k);
        let seq = &b"ACGTTCAGGATCCGTAACCGTTAGCATGCTA"[..k];
        let word = encode_kmer(seq).expect("encode");
        assert_eq!(
            reverse_complement_word(word, k, mask),
            reverse_complement_word_slow(word, k)
        );
    }
}

#[test]
fn compacts_linear_transcript() {
    let tr = Transcript {
        name: "tx".into(),
        original_len: 36,
        seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
    };
    let (graph, _) = build_kmer_unitig_graph_with_report(vec![tr], 31, 23, -1, 1).expect("graph");
    assert!(graph.unitigs.len() < graph.unitigs[0].len());
    assert_eq!(graph.km_unitigs.len(), 0);
    assert_eq!(graph.nodes.len(), graph.unitigs.len());
}
