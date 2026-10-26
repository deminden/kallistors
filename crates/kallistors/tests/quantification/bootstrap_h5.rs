use crate::helpers::{
    json_u64, read_abundance, run_kallistors, synthetic_transcripts, write_fasta, write_fastq_pairs,
};
use std::ffi::OsStr;
use std::fs;
use tempfile::TempDir;

#[test]
fn seeded_bootstraps_and_main_abundances_agree_between_hdf5_and_tsv() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("transcripts.fa");
    let index = dir.path().join("index.idx");
    let reads1 = dir.path().join("reads_1.fq");
    let reads2 = dir.path().join("reads_2.fq");
    let hdf_out = dir.path().join("hdf");
    let text_out = dir.path().join("text");
    let transcripts = synthetic_transcripts();
    write_fasta(&fasta, &transcripts).unwrap();
    write_fastq_pairs(&reads1, &reads2, &transcripts, 50, 90, 24).unwrap();
    assert!(
        run_kallistors(&[
            OsStr::new("index"),
            OsStr::new("-i"),
            index.as_os_str(),
            fasta.as_os_str(),
        ])
        .unwrap()
        .success()
    );
    for (out, plaintext) in [(&hdf_out, false), (&text_out, true)] {
        let mut args = vec![
            OsStr::new("quant"),
            OsStr::new("-t"),
            OsStr::new("2"),
            OsStr::new("-i"),
            index.as_os_str(),
            OsStr::new("-o"),
            out.as_os_str(),
            OsStr::new("-b"),
            OsStr::new("3"),
            OsStr::new("--seed"),
            OsStr::new("11"),
        ];
        if plaintext {
            args.push(OsStr::new("--plaintext"));
        }
        args.extend([reads1.as_os_str(), reads2.as_os_str()]);
        assert!(run_kallistors(&args).unwrap().success());
        let info = fs::read_to_string(out.join("run_info.json")).unwrap();
        assert_eq!(json_u64(&info, "n_bootstraps"), Some(3));
        assert_eq!(json_u64(&info, "n_processed"), Some(24));
        assert_eq!(json_u64(&info, "n_pseudoaligned"), Some(24));
    }
    assert!(!text_out.join("abundance.h5").exists());
    let h5 = hdf5_pure::File::open(hdf_out.join("abundance.h5")).unwrap();
    let ids = h5.dataset("aux/ids").unwrap().read_string().unwrap();
    assert_eq!(
        ids,
        transcripts
            .iter()
            .map(|(name, _)| name.clone())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        h5.dataset("aux/num_bootstrap").unwrap().read_i32().unwrap(),
        vec![3]
    );
    let main = read_abundance(hdf_out.join("abundance.tsv"));
    assert_eq!(main, read_abundance(text_out.join("abundance.tsv")));
    let counts = h5.dataset("est_counts").unwrap().read_f64().unwrap();
    let effective = h5.dataset("aux/eff_lengths").unwrap().read_f64().unwrap();
    let lengths = h5.dataset("aux/lengths").unwrap().read_i32().unwrap();
    assert_eq!(counts.len(), ids.len());
    for (i, name) in ids.iter().enumerate() {
        let row = main[name];
        assert_eq!(lengths[i], row.0 as i32);
        assert!((effective[i] - row.1).abs() < 1e-9);
        assert!((counts[i] - row.2).abs() < 1e-9);
    }
    assert_eq!(
        h5.dataset("aux/bias_observed").unwrap().read_i32().unwrap(),
        vec![1; 4096]
    );
    assert_eq!(
        h5.dataset("aux/bias_normalized")
            .unwrap()
            .read_f64()
            .unwrap(),
        vec![1.0; 4096]
    );
    let mut bootstraps = Vec::new();
    for i in 0..3 {
        let values = h5
            .dataset(&format!("bootstrap/bs{i}"))
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(values.len(), ids.len());
        assert!(values.iter().all(|v| v.is_finite() && *v >= 0.0));
        assert!((values.iter().sum::<f64>() - 24.0).abs() < 1e-8);
        let rows = read_abundance(text_out.join(format!("bs_abundance_{i}.tsv")));
        assert_eq!(rows.len(), ids.len());
        for (j, name) in ids.iter().enumerate() {
            assert!((values[j] - rows[name].2).abs() < 1e-9);
            assert_eq!(rows[name].0, main[name].0);
            assert_eq!(rows[name].1, main[name].1);
        }
        bootstraps.push(values);
    }
    // Fixed seed and nonuniform EC counts: this detects placeholder/copy output.
    assert!(bootstraps.windows(2).any(|pair| pair[0] != pair[1]));
}

#[test]
fn large_hdf5_vectors_keep_bounded_chunks_and_exact_tail_values() {
    use kallistors::quant::{H5QuantOutput, QuantResult, write_abundance_h5};
    let temp = TempDir::new().unwrap();
    let path = temp.path().join("large.h5");
    let n = 70_003;
    let names: Vec<_> = (0..n).map(|i| format!("transcript_{i:06}")).collect();
    let ids: Vec<_> = names.iter().map(String::as_str).collect();
    let lengths = vec![500u32; n];
    let counts: Vec<_> = (0..n).map(|i| (i % 123) as f64 * 0.125).collect();
    let result = QuantResult {
        est_counts: counts.clone(),
        tpm: vec![0.0; n],
        eff_lengths: vec![300.0; n],
        post_bias: None,
    };
    write_abundance_h5(
        &path,
        H5QuantOutput {
            target_names: &ids,
            target_lengths: &lengths,
            result: &result,
            bootstrap_est_counts: &[&counts],
            reads_processed: 100,
            fragment_length_hist: None,
            bias_observed: None,
            index_version: 13,
            kallisto_version: "0.52.0",
            call: "large-vector-regression",
            start_time: "",
        },
    )
    .unwrap();
    let h5 = hdf5_pure::File::open(&path).unwrap();
    for name in [
        "est_counts",
        "aux/ids",
        "aux/eff_lengths",
        "aux/lengths",
        "bootstrap/bs0",
    ] {
        let dataset = h5.dataset(name).unwrap();
        assert_eq!(dataset.shape().unwrap(), vec![n as u64]);
        assert_eq!(dataset.chunk_shape().unwrap(), Some(vec![16_384]));
        assert_eq!(dataset.chunks().unwrap().len(), n.div_ceil(16_384));
    }
    assert_eq!(
        h5.dataset("est_counts").unwrap().read_f64().unwrap(),
        counts
    );
    assert_eq!(
        h5.dataset("bootstrap/bs0").unwrap().read_f64().unwrap(),
        counts
    );
    assert_eq!(h5.dataset("aux/ids").unwrap().read_string().unwrap(), names);
    assert_eq!(
        h5.dataset("aux/lengths").unwrap().read_i32().unwrap(),
        vec![500; n]
    );
}

#[test]
fn parallel_bootstraps_preserve_original_serial_seed_order_and_bits() {
    use kallistors::ec::EcList;
    use kallistors::quant::{
        EcCountsInput, QuantOptions, bootstrap_quantify, bootstrap_quantify_threaded,
    };
    let input = EcCountsInput {
        ec_list: EcList {
            classes: vec![
                vec![0],
                vec![1],
                vec![0, 1],
                vec![1, 2, 3],
                vec![0, 3],
                vec![2, 3],
            ],
        },
        counts: vec![2, 3, 11, 37, 41, 0],
    };
    let lengths = [1.0, 13.0, 180.0, 512.75];
    let options = QuantOptions {
        max_iter: 512,
        min_rounds: 3,
        ..Default::default()
    };
    for seed in [0, 42, u64::MAX] {
        let expected = bootstrap_quantify(&input, &lengths, 7, seed, options).unwrap();
        for workers in [1, 2, 3, 8] {
            let actual =
                bootstrap_quantify_threaded(&input, &lengths, 7, seed, options, workers).unwrap();
            assert_eq!(actual.len(), expected.len());
            for (actual, expected) in actual.iter().zip(&expected) {
                for (a, b) in [
                    (&actual.est_counts, &expected.est_counts),
                    (&actual.eff_lengths, &expected.eff_lengths),
                    (&actual.tpm, &expected.tpm),
                ] {
                    assert_eq!(
                        a.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        b.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                    );
                }
            }
        }
    }
    assert!(
        bootstrap_quantify_threaded(&input, &lengths, 0, 42, options, 8)
            .unwrap()
            .is_empty()
    );
    let expected = bootstrap_quantify(&input, &lengths, 1, 42, options).unwrap();
    let actual = bootstrap_quantify_threaded(&input, &lengths, 1, 42, options, 8).unwrap();
    assert_eq!(actual[0].est_counts, expected[0].est_counts);
    let invalid = EcCountsInput {
        ec_list: EcList {
            classes: vec![vec![0]],
        },
        counts: vec![],
    };
    for samples in [0, 1, 7] {
        assert!(bootstrap_quantify_threaded(&invalid, &lengths, samples, 42, options, 8).is_err());
    }
}
