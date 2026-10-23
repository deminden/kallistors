//! Long reads regressions.

use super::*;
#[test]
fn bus_batch_long_writes_one_flens_line_per_batch() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let a_r1 = dir.path().join("a_r1.fastq");
    let a_r2 = dir.path().join("a_r2.fastq");
    let b_r1 = dir.path().join("b_r1.fastq");
    let b_r2 = dir.path().join("b_r2.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("batch_long_flens_out");

    write_fastq(&a_r1, &[("a_cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&a_r2, &[("a_seq", &transcript[..40])]);
    write_fastq(&b_r1, &[("b_cell", b"TGCATGCATGCATGCATTTTTTTTTTTT")]);
    write_fastq(&b_r2, &[("b_seq", &transcript[..60])]);
    fs::write(
        &batch,
        format!(
            "sample_a\t{}\t{}\nsample_b\t{}\t{}\n",
            a_r1.display(),
            a_r2.display(),
            b_r1.display(),
            b_r2.display()
        ),
    )
    .expect("write batch");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("-B")
        .arg(&batch)
        .arg("--long")
        .status()
        .expect("run kallistors bus batch --long");
    assert!(status.success());

    let flens = fs::read_to_string(out_dir.join("flens.txt")).unwrap();
    let lines = flens.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2);
    let values = lines
        .iter()
        .map(|line| line.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(values, vec![9.0, 29.0]);
}

#[test]
fn bus_long_filters_reads_above_unmapped_threshold() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("long_out");

    write_fastq(
        &r1,
        &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_novel", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read", &transcript), ("seq_novel", &[b'N'; 80])],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--long")
        .arg("--threshold")
        .arg("0.5")
        .arg("--unmapped")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --long");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(2)
    );
    let ratios = fs::read_to_string(out_dir.join("unmapped_ratio.txt")).unwrap();
    let values = ratios
        .trim_end()
        .split(',')
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(values.len(), 2);
    assert!(values[0] <= 0.5);
    assert!(values[1] > 0.5);
    let flens = fs::read_to_string(out_dir.join("flens.txt")).unwrap();
    let flen_values = flens
        .split_whitespace()
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(flen_values.len(), 1);
    assert!(flen_values[0] >= 0.0);
    assert_eq!(
        fs::read(out_dir.join("index.saved")).unwrap(),
        fs::read(&index).unwrap()
    );
    let novel = fs::read_to_string(out_dir.join("novel.fastq")).unwrap();
    let novel_lines = novel.lines().collect::<Vec<_>>();
    assert_eq!(novel_lines.len(), 8);
    assert_eq!(novel_lines[0], "@unmapped");
    assert_eq!(novel_lines[1], std::str::from_utf8(&[b'N'; 80]).unwrap());
    assert_eq!(novel_lines[2], "+");
    assert_eq!(novel_lines[3], "I".repeat(80));
    assert_eq!(novel_lines[4], "@novel_disjointIntersect");
    assert_eq!(novel_lines[5], std::str::from_utf8(&[b'N'; 80]).unwrap());
    assert_eq!(novel_lines[6], "+");
    assert_eq!(novel_lines[7], "I".repeat(80));
}

#[test]
fn bus_long_computes_threshold_from_error_rate() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("long_error_rate_out");

    write_fastq(
        &r1,
        &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_novel", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read", &transcript), ("seq_novel", &[b'N'; 80])],
    );

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--long")
        .arg("-e")
        .arg("0.01")
        .arg("-P")
        .arg("ONT")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --long --error-rate");
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Using computed threshold 0.38"));

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}

#[test]
fn bus_long_computed_threshold_uses_index_k() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT";
    let (index, transcript) = build_index_with_named_transcript_and_k(&dir, "tx0", transcript, 15);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("long_error_rate_k15_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--long")
        .arg("-e")
        .arg("0.01")
        .arg("-P")
        .arg("ONT")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --long --error-rate on k15 index");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Using computed threshold 0.7"), "{stderr}");

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}

#[test]
fn bus_long_ignores_paired_flag_for_single_cdna_technologies() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("long_paired_flag_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--long")
        .arg("--paired")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --long --paired");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    assert!(out_dir.join("flens.txt").exists());
    assert!(!out_dir.join("novel.fastq").exists());
}

#[test]
fn bus_smartseq2_long_paired_uses_four_input_files_without_paired_bus() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bc1 = dir.path().join("bc1.fastq");
    let bc2 = dir.path().join("bc2.fastq");
    let read1 = dir.path().join("read1.fastq");
    let read2 = dir.path().join("read2.fastq");
    let out_dir = dir.path().join("smartseq2_long_paired_out");

    write_fastq(&bc1, &[("bc1", b"ACGT")]);
    write_fastq(&bc2, &[("bc2", b"TGCA")]);
    write_fastq(&read1, &[("read1", &transcript)]);
    write_fastq(&read2, &[("read2", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("SmartSeq2")
        .arg("--long")
        .arg("--paired")
        .arg(&bc1)
        .arg(&bc2)
        .arg(&read1)
        .arg(&read2)
        .status()
        .expect("run kallistors bus SmartSeq2 --long --paired");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"ACGTTGCA"));
    assert!(out_dir.join("flens.txt").exists());
}

#[test]
fn bus_platform_errors_clearly() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("platform_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read lane:1", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--long")
        .arg("--platform")
        .arg("bad")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --platform bad");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--platform must be PACBIO or ONT"));
}

#[test]
fn bus_ignores_long_read_knobs_when_not_long() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("short_read_long_knobs_out");

    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--threshold")
        .arg("2")
        .arg("--error-rate")
        .arg("0")
        .arg("--platform")
        .arg("not-a-platform")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus with inert long-read knobs");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}

#[test]
fn bus_long_invalid_threshold_uses_default() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("long_invalid_threshold_out");

    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--long")
        .arg("--threshold")
        .arg("2")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --long with invalid threshold");
    assert!(output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains(
            "Threshold not in (0,1). Setting default threshold for unmapped kmers to 0.8"
        )
    );

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}
