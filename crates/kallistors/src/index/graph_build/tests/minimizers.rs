//! Minimizers regressions.

use super::*;
#[test]
fn rolling_minimizer_hashes_match_byte_oracle_at_every_window() {
    let mut state = 42u64;
    let random = (0..512)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect::<Vec<_>>();
    for seq in [random, vec![b'A'; 128], b"acgt".repeat(32)] {
        for g in 3..=29 {
            let actual = minimizer_rep_hashes(&seq, g).unwrap();
            let expected = seq
                .windows(g)
                .map(|window| minimizer_rep_hash(window).unwrap())
                .collect::<Vec<_>>();
            assert_eq!(actual, expected, "g={g}");
        }
    }
    assert!(minimizer_rep_hashes(b"ACNTACGT", 3).is_err());
    assert!(minimizer_rep_hashes(b"ACGTACNT", 3).is_err());
}

#[test]
fn sliding_minimizer_selection_matches_naive_reference() {
    let mut scratch = MinimizerScratch::default();
    let mut rng = 1234567u64;
    let random = (0..128)
        .map(|_| {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            b"ACGT"[(rng & 3) as usize]
        })
        .collect::<Vec<_>>();
    let repeated = b"ACGT".repeat(32);
    let tied = vec![b'C'; 128];
    for seq in [&random, &repeated, &tied] {
        for k in (5..=31).step_by(2) {
            for g in [3, k - 2] {
                for base in [0, u64::from(u32::MAX) - 64] {
                    let mut keys = Vec::new();
                    let mut blocks = vec![Vec::new(); if base == 0 { 1 } else { 2 }];
                    insert_selected_minimizers(
                        seq,
                        k,
                        g,
                        base,
                        &mut keys,
                        &mut blocks,
                        &mut scratch,
                    )
                    .unwrap();
                    let positions = blocks.iter().enumerate().flat_map(|(block, positions)| {
                        positions
                            .iter()
                            .map(move |&pos| ((block as u64) << 32) | u64::from(pos))
                    });
                    let actual = positions
                        .zip(keys.into_iter().map(u64::to_be_bytes))
                        .collect::<Vec<_>>();
                    let mut expected = selected_minimizers_naive(seq, k, g);
                    expected.sort_unstable();
                    expected.dedup_by_key(|(pos, _)| *pos);
                    let expected = expected
                        .into_iter()
                        .map(|(pos, key)| (base + u64::from(pos), key))
                        .collect::<Vec<_>>();
                    assert_eq!(actual, expected, "k={k},g={g},base={base}");
                }
            }
        }
    }
    for position in 0..random.len() {
        let mut invalid = random.clone();
        invalid[position] = b'N';
        assert!(
            insert_selected_minimizers(
                &invalid,
                31,
                23,
                0,
                &mut Vec::new(),
                &mut [Vec::new()],
                &mut scratch
            )
            .is_err(),
            "invalid position={position}"
        );
    }
}
