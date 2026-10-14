use kallistors::ec::EcList;
use kallistors::quant::{EcCountsInput, QuantOptions, bootstrap_quantify, em_quantify};

// Preserve the dense pre-optimization update as a numerical regression oracle.
fn dense_counts(input: &EcCountsInput, lengths: &[f64], options: QuantOptions) -> Vec<f64> {
    let mut alpha = vec![1.0 / lengths.len() as f64; lengths.len()];
    let mut next = vec![0.0; lengths.len()];
    let mut final_round = false;
    for iteration in 0..options.max_iter {
        next.fill(0.0);
        for (ec, &count) in input.ec_list.classes.iter().zip(&input.counts) {
            if ec.len() == 1 {
                next[ec[0] as usize] = f64::from(count);
            }
        }
        for (ec, &count) in input.ec_list.classes.iter().zip(&input.counts) {
            if ec.len() < 2 || count == 0 {
                continue;
            }
            let count = f64::from(count);
            let mut denominator = 0.0;
            for &transcript in ec {
                let tr = transcript as usize;
                denominator += alpha[tr] * (count / lengths[tr]);
            }
            if denominator < f64::from_bits(1) {
                continue;
            }
            let scale = count / denominator;
            for &transcript in ec {
                let tr = transcript as usize;
                next[tr] += (count / lengths[tr]) * alpha[tr] * scale;
            }
        }
        let mut changes = 0;
        for i in 0..alpha.len() {
            if next[i] > 1e-2 && (next[i] - alpha[i]).abs() / next[i] > 1e-2 {
                changes += 1;
            }
            alpha[i] = next[i];
        }
        if final_round {
            break;
        }
        if changes == 0 && iteration > options.min_rounds {
            final_round = true;
            for value in &mut alpha {
                if *value < 1e-7 / 10.0 {
                    *value = 0.0;
                }
            }
        }
    }
    alpha
}

#[test]
fn sparse_em_preserves_dense_updates_and_zero_iteration_initialization() {
    let cases = [
        EcCountsInput {
            ec_list: EcList {
                classes: vec![vec![0], vec![2], vec![1, 2], vec![3, 4], vec![99]],
            },
            counts: vec![2, 3, 8, 0, 1],
        },
        EcCountsInput {
            ec_list: EcList {
                classes: vec![vec![0], vec![1, 2]],
            },
            counts: vec![0, 0],
        },
    ];
    let lengths: Vec<u32> = (0..128).map(|i| 500 + i * 7).collect();
    for input in cases {
        for max_iter in [0, 1, 2, 100] {
            let options = QuantOptions {
                max_iter,
                min_rounds: 7,
                ..QuantOptions::default()
            };
            let result = em_quantify(&input, &lengths, None, None, options).expect("quantify");
            assert_eq!(
                result.est_counts,
                dense_counts(&input, &result.eff_lengths, options)
            );
            if max_iter > 0 {
                assert_eq!(result.est_counts[127], 0.0);
            }
        }
    }
}

#[test]
fn sparse_bootstrap_updates_leave_unobserved_transcripts_zero() {
    let input = EcCountsInput {
        ec_list: EcList {
            classes: vec![vec![1], vec![1, 2]],
        },
        counts: vec![4, 8],
    };
    let options = QuantOptions {
        max_iter: 100,
        min_rounds: 7,
        ..QuantOptions::default()
    };
    let results = bootstrap_quantify(&input, &vec![800.0; 128], 3, 42, options).expect("bootstrap");
    for result in results {
        assert_eq!(result.est_counts[0], 0.0);
        assert_eq!(result.est_counts[127], 0.0);
        assert!((result.est_counts.iter().sum::<f64>() - 12.0).abs() < 1e-10);
    }
}
