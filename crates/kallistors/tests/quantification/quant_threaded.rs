use kallistors::bias::{BiasCounts, hexamer_to_int};
use kallistors::ec::EcList;
use kallistors::quant::{
    EcCountsInput, QuantOptions, em_quantify_with_fragment_histogram,
    em_quantify_with_fragment_histogram_threaded,
};

#[test]
fn threaded_em_preserves_serial_bits_with_bias_histograms_and_sparse_targets() {
    let lengths: Vec<u32> = (0..32).map(|i| 1000 + 7 * i).collect();
    let sequences: Vec<Vec<u8>> = lengths
        .iter()
        .enumerate()
        .map(|(tr, &len)| {
            let mut state = 0x9e37_79b9_7f4a_7c15u64 ^ tr as u64;
            (0..len)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect()
        })
        .collect();
    let mut classes = vec![vec![0], vec![2], vec![31]];
    let mut counts = vec![0, 3, 7];
    for i in 0..24 {
        classes.push(vec![i, (i + 3) % 24, (i + 7) % 24]);
        counts.push((i * 19) % 101);
    }
    let mut disconnected_classes = classes[..3].to_vec();
    for i in 0..24 {
        let component = i % 3;
        let local = i / 3;
        disconnected_classes.push(vec![
            i,
            ((local + 3) % 8) * 3 + component,
            ((local + 7) % 8) * 3 + component,
        ]);
    }
    let disconnected = EcCountsInput {
        ec_list: EcList {
            classes: disconnected_classes,
        },
        counts: counts.clone(),
    };
    let input = EcCountsInput {
        ec_list: EcList { classes },
        counts,
    };
    let mut histogram = vec![0; 301];
    histogram[100] = 3;
    histogram[300] = 7;
    let mut bias = BiasCounts::new();
    for (tr, seq) in sequences.iter().enumerate() {
        for window in seq.windows(6) {
            for reverse in [false, true] {
                for _ in 0..1 + tr % 5 {
                    bias.record(hexamer_to_int(window, reverse).unwrap());
                }
            }
        }
    }
    for input in [input, disconnected] {
        for hist in [None, Some(histogram.as_slice())] {
            for use_bias in [false, true] {
                let minimums: &[usize] = if use_bias {
                    &[0, 1]
                } else {
                    &[0, 1, usize::MAX]
                };
                for &min_rounds in minimums {
                    for max_iter in [0, 1, 2, 100] {
                        let options = QuantOptions {
                            max_iter,
                            min_rounds,
                            bias: use_bias,
                            ..QuantOptions::default()
                        };
                        let serial = em_quantify_with_fragment_histogram(
                            &input,
                            &lengths,
                            Some(&sequences),
                            Some(&bias),
                            hist,
                            options,
                        )
                        .unwrap();
                        assert!(
                            serial
                                .est_counts
                                .iter()
                                .chain(&serial.tpm)
                                .chain(&serial.eff_lengths)
                                .all(|v| v.is_finite() && *v >= 0.0)
                        );
                        if max_iter > 0 {
                            let expected_mass: f64 = input.counts.iter().map(|&v| v as f64).sum();
                            assert!(
                                (serial.est_counts.iter().sum::<f64>() - expected_mass).abs()
                                    < 1e-8
                            );
                            assert!((serial.tpm.iter().sum::<f64>() - 1e6).abs() < 1e-8);
                            assert!(serial.est_counts[25..31].iter().all(|&v| v == 0.0));
                        }
                        assert_eq!(
                            serial.post_bias.is_some(),
                            use_bias && max_iter > min_rounds
                        );
                        if use_bias && max_iter > min_rounds {
                            let unbiased = em_quantify_with_fragment_histogram(
                                &input,
                                &lengths,
                                Some(&sequences),
                                Some(&bias),
                                hist,
                                QuantOptions {
                                    bias: false,
                                    ..options
                                },
                            )
                            .unwrap();
                            assert_ne!(
                                serial.eff_lengths, unbiased.eff_lengths,
                                "bias fixture must exercise an actual correction"
                            );
                        }
                        for (threads, owned) in [(2, false), (8, false), (8, true)] {
                            let parallel = if owned {
                                kallistors::quant::em_quantify_owned_with_fragment_histogram_threaded(
                                kallistors::quant::EcCountsInput { ec_list: kallistors::ec::EcList { classes: input.ec_list.classes.clone() }, counts: input.counts.clone() },
                                &lengths,
                                Some(&sequences),
                                Some(&bias),
                                hist,
                                options,
                                threads,
                            )
                            .unwrap()
                            } else {
                                em_quantify_with_fragment_histogram_threaded(
                                    &input,
                                    &lengths,
                                    Some(&sequences),
                                    Some(&bias),
                                    hist,
                                    options,
                                    threads,
                                )
                                .unwrap()
                            };
                            assert_eq!(
                                serial
                                    .est_counts
                                    .iter()
                                    .map(|x| x.to_bits())
                                    .collect::<Vec<_>>(),
                                parallel
                                    .est_counts
                                    .iter()
                                    .map(|x| x.to_bits())
                                    .collect::<Vec<_>>()
                            );
                            assert_eq!(
                                serial.tpm.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                                parallel.tpm.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                            );
                            assert_eq!(serial.eff_lengths, parallel.eff_lengths);
                            assert_eq!(serial.post_bias, parallel.post_bias);
                        }
                    }
                }
            }
        }
    }
}
