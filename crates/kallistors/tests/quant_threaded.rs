use kallistors::bias::BiasCounts;
use kallistors::ec::EcList;
use kallistors::quant::{
    EcCountsInput, QuantOptions, em_quantify_with_fragment_histogram,
    em_quantify_with_fragment_histogram_threaded,
};

#[test]
fn threaded_em_preserves_serial_bits_with_bias_histograms_and_sparse_targets() {
    let lengths: Vec<u32> = (0..32).map(|i| 500 + 7 * i).collect();
    let sequences: Vec<Vec<u8>> = lengths
        .iter()
        .enumerate()
        .map(|(tr, &len)| {
            (0..len as usize)
                .map(|pos| b"ACGT"[(pos + pos / 11 + tr) % 4])
                .collect()
        })
        .collect();
    let mut classes = vec![vec![0], vec![2], vec![31]];
    let mut counts = vec![0, 3, 7];
    for i in 0..24 {
        classes.push(vec![i, (i + 3) % 24, (i + 7) % 24]);
        counts.push((i * 19) % 101);
    }
    let input = EcCountsInput {
        ec_list: EcList { classes },
        counts,
    };
    let mut histogram = vec![0; 301];
    histogram[100] = 3;
    histogram[300] = 7;
    let mut bias = BiasCounts::new();
    for hex in 0..4096 {
        bias.record(hex);
    }
    for hist in [None, Some(histogram.as_slice())] {
        for use_bias in [false, true] {
            for max_iter in [0, 1, 2, 100] {
                let options = QuantOptions {
                    max_iter,
                    min_rounds: 1,
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
                for threads in [1, 2, 8] {
                    let parallel = em_quantify_with_fragment_histogram_threaded(
                        &input,
                        &lengths,
                        Some(&sequences),
                        Some(&bias),
                        hist,
                        options,
                        threads,
                    )
                    .unwrap();
                    assert_eq!(serial.est_counts, parallel.est_counts);
                    assert_eq!(serial.tpm, parallel.tpm);
                    assert_eq!(serial.eff_lengths, parallel.eff_lengths);
                    assert_eq!(serial.post_bias, parallel.post_bias);
                }
            }
        }
    }
}
