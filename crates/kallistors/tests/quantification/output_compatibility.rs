use crate::helpers::read_abundance;
use kallistors::quant::{QuantResult, RunInfo, write_abundance_tsv_from_parts};

#[test]
fn run_info_preserves_schema_values_and_escaped_strings() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run_info.json");
    let info = RunInfo {
        n_targets: 3,
        n_bootstraps: 2,
        n_processed: 10,
        n_pseudoaligned: 9,
        n_unique: 7,
        p_pseudoaligned: 90.0,
        p_unique: 70.0,
        kallisto_version: "test-build".into(),
        index_version: 13,
        kmer_length: 31,
        start_time: "test-time".into(),
        call: "quant \"reads\"\\path\n\r\t\u{1}é".into(),
    };
    info.write_json(&path).unwrap();
    // An independent parser checks escaping without pinning whitespace or field order.
    let json: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    for (field, expected) in [
        ("n_targets", 3),
        ("n_bootstraps", 2),
        ("n_processed", 10),
        ("n_pseudoaligned", 9),
        ("n_unique", 7),
        ("index_version", 13),
        ("k-mer length", 31),
    ] {
        assert_eq!(json[field].as_u64(), Some(expected), "{field}");
    }
    assert_eq!(json["p_pseudoaligned"].as_f64(), Some(90.0));
    assert_eq!(json["p_unique"].as_f64(), Some(70.0));
    assert_eq!(
        json["kallisto_version"].as_str(),
        Some(info.kallisto_version.as_str())
    );
    assert_eq!(json["start_time"].as_str(), Some(info.start_time.as_str()));
    assert_eq!(json["call"].as_str(), Some(info.call.as_str()));
    let allocation = kallistors::pseudoalign::MatchCacheSize::from_bytes(96)
        .unwrap()
        .allocation(2)
        .unwrap();
    let cache_path = dir.path().join("with-cache.json");
    info.write_json_with_match_cache(&cache_path, allocation)
        .unwrap();
    let with_cache: serde_json::Value =
        serde_json::from_slice(&std::fs::read(cache_path).unwrap()).unwrap();
    for (field, value) in json.as_object().unwrap() {
        assert_eq!(&with_cache[field], value, "existing field {field}");
    }
    assert_eq!(
        with_cache["match_cache"],
        serde_json::json!({
            "policy":"96B", "requested_bytes":96, "allocated_bytes":64,
            "bytes_per_worker":32, "workers":2,
        })
    );
}

#[test]
fn abundance_tsv_preserves_transcript_identity_and_numeric_values() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("abundance.tsv");
    let result = QuantResult {
        est_counts: vec![2.25, 5.75],
        tpm: vec![200_000.0, 800_000.0],
        eff_lengths: vec![50.5, 150.25],
        post_bias: None,
    };
    write_abundance_tsv_from_parts(&path, &["tx0", "tx1"], &[100, 200], &result).unwrap();
    let rows = read_abundance(path);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows["tx0"], (100, 50.5, 2.25, 200_000.0));
    assert_eq!(rows["tx1"], (200, 150.25, 5.75, 800_000.0));
}

#[test]
fn parallel_abundance_rows_match_serial_bytes_across_waves_and_short_columns() {
    let dir = tempfile::tempdir().unwrap();
    let names: Vec<_> = (0..32771).map(|i| format!("tx_{i}_é")).collect();
    let names: Vec<_> = names.iter().map(String::as_str).collect();
    let lengths = vec![100u32; names.len() - 3];
    let values = [
        -0.0,
        0.0,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::INFINITY,
        f64::NAN,
        1e-320,
        1.0 / 3.0,
    ];
    let make = |len| (0..len).map(|i| values[i % values.len()]).collect();
    let result = QuantResult {
        est_counts: make(names.len() - 2),
        tpm: make(names.len()),
        eff_lengths: make(names.len() - 1),
        post_bias: None,
    };
    let serial = dir.path().join("serial.tsv");
    write_abundance_tsv_from_parts(&serial, &names, &lengths, &result).unwrap();
    let expected = std::fs::read(serial).unwrap();
    for threads in [2, 8] {
        let path = dir.path().join(format!("parallel{threads}.tsv"));
        kallistors::quant::write_abundance_tsv_from_parts_threaded(
            &path, &names, &lengths, &result, threads,
        )
        .unwrap();
        assert_eq!(std::fs::read(path).unwrap(), expected, "threads={threads}");
    }
}
