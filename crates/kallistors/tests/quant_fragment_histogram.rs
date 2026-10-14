use kallistors::ec::EcList;
use kallistors::index::{IndexBuildOptions, build_index};
use kallistors::quant::{EcCountsInput, QuantOptions, em_quantify_with_fragment_histogram};
use std::fs;
use std::process::Command;

#[test]
fn empirical_fragment_model_uses_conditional_means_and_endpoint_observations() {
    let input = EcCountsInput {
        ec_list: EcList {
            classes: vec![vec![0], vec![1], vec![2], vec![3]],
        },
        counts: vec![3, 2, 1, 4],
    };
    let mut histogram = vec![0; 301];
    histogram[50] = 2;
    histogram[100] = 1;
    histogram[200] = 2;
    let result = em_quantify_with_fragment_histogram(
        &input,
        &[50, 100, 199, 1000],
        None,
        None,
        Some(&histogram),
        QuantOptions::default(),
    )
    .expect("empirical quantification");
    // Including the length-100 observation gives a conditional mean of 200/3.
    // The full distribution's mean is 120, including for targets beyond its support.
    let expected = [1.0, 101.0 - 200.0 / 3.0, 200.0 - 200.0 / 3.0, 881.0];
    for (actual, expected) in result.eff_lengths.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-10);
    }
    assert_eq!(result.est_counts, [3.0, 2.0, 1.0, 4.0]);
    let mass: f64 = result
        .est_counts
        .iter()
        .zip(expected)
        .map(|(count, length)| count / length)
        .sum();
    for ((count, length), tpm) in result.est_counts.iter().zip(expected).zip(result.tpm) {
        assert!((tpm - count / length / mass * 1_000_000.0).abs() < 1e-8);
    }
}

#[test]
fn empty_empirical_fragment_distribution_is_rejected() {
    let input = EcCountsInput {
        ec_list: EcList {
            classes: vec![vec![0]],
        },
        counts: vec![1],
    };
    for histogram in [&[][..], &[0, 0, 0][..]] {
        assert!(
            em_quantify_with_fragment_histogram(
                &input,
                &[1000],
                None,
                None,
                Some(histogram),
                QuantOptions::default(),
            )
            .is_err()
        );
    }
}

#[test]
fn paired_quant_command_passes_observed_histogram_to_em() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mut seed = 0x5a17_f35bu64;
    let sequence: Vec<u8> = (0..1000)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize]
        })
        .collect();
    let fasta = dir.path().join("transcripts.fa");
    fs::write(
        &fasta,
        format!(">tx\n{}\n", String::from_utf8_lossy(&sequence)),
    )
    .expect("write transcript");
    let index = dir.path().join("index.idx");
    build_index(&index, &[fasta], IndexBuildOptions::default()).expect("build index");
    let mut left = String::new();
    let mut right = String::new();
    for (pair, (start, length)) in [(100, 100), (400, 300)].into_iter().enumerate() {
        let read1 = &sequence[start..start + 50];
        let read2: Vec<u8> = sequence[start + length - 50..start + length]
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            })
            .collect();
        left.push_str(&format!(
            "@r{pair}\n{}\n+\n{}\n",
            String::from_utf8_lossy(read1),
            "I".repeat(50)
        ));
        right.push_str(&format!(
            "@r{pair}\n{}\n+\n{}\n",
            String::from_utf8_lossy(&read2),
            "I".repeat(50)
        ));
    }
    let reads1 = dir.path().join("reads1.fastq");
    let reads2 = dir.path().join("reads2.fastq");
    fs::write(&reads1, left).expect("write mate 1");
    fs::write(&reads2, right).expect("write mate 2");
    let out = dir.path().join("quant");
    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("quant")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out)
        .args(["-t", "8"])
        .arg(reads1)
        .arg(reads2)
        .output()
        .expect("run quant");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let abundance = fs::read_to_string(out.join("abundance.tsv")).expect("abundance");
    let row: Vec<&str> = abundance
        .lines()
        .nth(1)
        .expect("target row")
        .split('\t')
        .collect();
    // Two empirical fragments of 100 and 300 bases have mean 200. A fitted,
    // zero-truncated Gaussian instead shifts that mean and fails this assertion.
    assert!((row[2].parse::<f64>().expect("effective length") - 801.0).abs() < 1e-10);
    assert_eq!(row[3].parse::<f64>().expect("estimated count"), 2.0);
}
