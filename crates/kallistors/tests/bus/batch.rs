//! Batch regressions.

use super::*;
#[test]
fn bus_packed_workers_preserve_order_flags_and_limits_across_batches() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let interleaved = dir.path().join("interleaved.fastq");
    // Cross batch two with skipped/ambiguous records; malformed tails beyond the limit must remain unread.
    let count = 4 * 65_536 + 77;
    let limit = 65_536 + 46;
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut combined = Vec::new();
    let mut expected = Vec::new();
    let mut prefix_ends = (0, 0, 0);
    for number in 0..count {
        let mut barcode_umi = b"ACGTACGTACGTACGTTTTTTTTTTTTT".to_vec();
        let flags = match number % 3 {
            1 => {
                barcode_umi[0] = b'N';
                1
            }
            2 => {
                barcode_umi[20] = b'N';
                17 << 8
            }
            _ => 0,
        };
        if number % 17 == 0 {
            barcode_umi.truncate(20);
        }
        let sequence = if number % 7 == 0 {
            vec![b'N'; transcript.len()]
        } else {
            transcript.clone()
        };
        let record = |seq: &[u8]| {
            format!(
                "@r{number}\n{}\n+\n{}\n",
                std::str::from_utf8(seq).unwrap(),
                "I".repeat(seq.len())
            )
            .into_bytes()
        };
        let a = record(&barcode_umi);
        let b = record(&sequence);
        left.extend_from_slice(&a);
        right.extend_from_slice(&b);
        combined.extend_from_slice(&a);
        combined.extend_from_slice(&b);
        if number + 1 == limit {
            prefix_ends = (left.len(), right.len(), combined.len());
        }
        if number < limit && number % 17 != 0 && number % 7 != 0 {
            expected.push((number, flags));
        }
    }
    fs::write(&r1, &left).unwrap();
    fs::write(&r2, &right).unwrap();
    fs::write(&interleaved, &combined).unwrap();
    let broken_r1 = dir.path().join("limited-r1.fastq");
    let broken_r2 = dir.path().join("limited-r2.fastq");
    let broken_interleaved = dir.path().join("limited-interleaved.fastq");
    let truncated = b"@unread\nACGT\n";
    fs::write(&broken_r1, [&left[..prefix_ends.0], truncated].concat()).unwrap();
    fs::write(&broken_r2, [&right[..prefix_ends.1], truncated].concat()).unwrap();
    fs::write(
        &broken_interleaved,
        [&combined[..prefix_ends.2], truncated].concat(),
    )
    .unwrap();
    drop((left, right, combined));
    for (case, interleave, numbered, broken_tail) in [
        ("paired-files", false, false, false),
        ("interleaved", true, false, false),
        ("numbered", false, true, false),
        ("trailing-truncation", false, false, true),
        ("interleaved-trailing-truncation", true, false, true),
        ("compound-tags", false, false, false),
    ] {
        let mut control = None;
        for threads in [1, 4] {
            let out = dir.path().join(format!("{case}-{threads}"));
            let mut command = bus_command(&index, &out);
            command.args([
                "-x",
                if case == "compound-tags" {
                    "0,0,8,0,8,16:0,16,22,0,22,28:1,0,0"
                } else {
                    "10xv3"
                },
                "-t",
                &threads.to_string(),
                "-N",
                &limit.to_string(),
            ]);
            if numbered {
                command.arg("--num");
            }
            if interleave {
                command.arg("--interleaved").arg(if broken_tail {
                    &broken_interleaved
                } else {
                    &interleaved
                });
            } else {
                command
                    .arg(if broken_tail { &broken_r1 } else { &r1 })
                    .arg(if broken_tail { &broken_r2 } else { &r2 });
            }
            let status =
                command_status_with_timeout(&mut command, Duration::from_secs(30)).unwrap();
            assert!(status.success(), "{case}/{threads}");
            let (_, _, records) = read_bus_records(&out.join("output.bus"));
            assert_eq!(records.len(), expected.len(), "{case}/{threads}");
            for (record, &(number, flags)) in records.iter().zip(&expected) {
                assert_eq!(
                    record.flags,
                    if numbered { number as u32 } else { flags },
                    "{case}/{threads}/{number}"
                );
            }
            let info = fs::read_to_string(out.join("run_info.json")).unwrap();
            assert_eq!(json_u64(&info, "n_processed"), Some(limit as u64));
            let files = ["output.bus", "matrix.ec", "transcripts.txt"]
                .map(|name| fs::read(out.join(name)).unwrap());
            if let Some(control) = &control {
                assert_eq!(&files, control, "{case}: threaded output order");
            } else {
                control = Some(files);
            }
        }
    }
}

#[test]
fn bus_without_technology_infers_paired_batch_rows() {
    let dir = tempfile::tempdir().expect("tempdir");
    let a_left = dir.path().join("a_left.fastq");
    let a_right = dir.path().join("a_right.fastq");
    let b_left = dir.path().join("b_left.fastq");
    let b_right = dir.path().join("b_right.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("bulk_no_technology_batch_paired_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(&a_left, &[("a/1", &transcript[..40])]);
    write_fastq(&a_right, &[("a/2", &transcript[20..60])]);
    write_fastq(&b_left, &[("b/1", &transcript[..45])]);
    write_fastq(&b_right, &[("b/2", &transcript[15..65])]);
    fs::write(
        &batch,
        format!(
            "sample_a\t{}\t{}\nsample_b\t{}\t{}\n",
            a_left.display(),
            a_right.display(),
            b_left.display(),
            b_right.display()
        ),
    )
    .expect("write batch");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("--batch")
        .arg(&batch)
        .output()
        .expect("run paired batch kallistors bus without -x");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("[bus] will try running read files supplied in batch file"),
        "{stderr}"
    );
    assert!(!stderr.contains("--paired ignored"), "{stderr}");

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].barcode, 0);
    assert_eq!(records[1].barcode, 1);
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\nsample_b\n"
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.sample.barcodes")).unwrap(),
        "AAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAC\n"
    );
    assert!(out_dir.join("flens.txt").exists());
}

#[test]
fn bus_without_technology_rejects_mixed_batch_row_widths() {
    let dir = tempfile::tempdir().expect("tempdir");
    let single = dir.path().join("single.fastq");
    let left = dir.path().join("left.fastq");
    let right = dir.path().join("right.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("bulk_no_technology_mixed_batch_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(&single, &[("single", &transcript)]);
    write_fastq(&left, &[("left", &transcript[..40])]);
    write_fastq(&right, &[("right", &transcript[20..60])]);
    fs::write(
        &batch,
        format!(
            "sample_a\t{}\nsample_b\t{}\t{}\n",
            single.display(),
            left.display(),
            right.display()
        ),
    )
    .expect("write batch");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("--batch")
        .arg(&batch)
        .arg("--paired")
        .output()
        .expect("run mixed-width batch kallistors bus without -x");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("batch file line 2 has 2 files, expected 1"),
        "{stderr}"
    );
    assert!(
        stderr.contains(
            "[bus] --paired ignored; single/paired-end is inferred from number of files supplied"
        ),
        "{stderr}"
    );
}

#[test]
fn bus_batch_barcode_prefixes_match_records_for_standard_custom_and_near_limit_layouts() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let a_r1 = dir.path().join("a_r1.fastq");
    let a_r2 = dir.path().join("a_r2.fastq");
    let b_r1 = dir.path().join("b_r1.fastq");
    let b_r2 = dir.path().join("b_r2.fastq");
    let batch = dir.path().join("batch.txt");
    write_fastq(&a_r2, &[("a_seq", &transcript)]);
    write_fastq(&b_r2, &[("b_seq", &transcript)]);
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
    .unwrap();
    for (case, technology, a_cell, b_cell, bc_len, umi_len) in [
        (
            "tenx",
            "10XV3",
            b"ACGTACGTACGTACGTTTTTTTTTTTTT".to_vec(),
            b"TGCATGCATGCATGCATTTTTTTTTTTT".to_vec(),
            16,
            12,
        ),
        (
            "custom",
            "0,0,8:0,8,20:1,0,0",
            b"ACGTACGTTTTTTTTTTTTT".to_vec(),
            b"TGCATGCATTTTTTTTTTTT".to_vec(),
            8,
            12,
        ),
        ("near_limit", "BDWTA", vec![b'A'; 60], vec![b'C'; 60], 27, 8),
    ] {
        write_fastq(&a_r1, &[("a_cell", &a_cell)]);
        write_fastq(&b_r1, &[("b_cell", &b_cell)]);
        let out = dir.path().join(case);
        let mut command = bus_command(&index, &out);
        command
            .args(["-x", technology, "--batch-barcodes", "-B"])
            .arg(&batch);
        if case == "tenx" {
            command.arg("--verbose");
        }
        let output = command.output().unwrap();
        assert_bus_success(&output, case);
        if case == "tenx" {
            let stderr = String::from_utf8_lossy(&output.stderr);
            for name in ["sample_a", "sample_b"] {
                assert!(
                    stderr.contains(&format!("will process sample {name}:")),
                    "{case}: {stderr}"
                );
            }
            for file in [&a_r1, &a_r2, &b_r1, &b_r2] {
                assert!(
                    stderr.contains(&file.display().to_string()),
                    "{case}: {stderr}"
                );
            }
        }
        let (header_bc, header_umi, records) = read_bus_records(&out.join("output.bus"));
        assert_eq!((header_bc, header_umi), (32, umi_len), "{case}");
        assert_eq!(records.len(), 2, "{case}");
        assert_eq!(
            fs::read_to_string(out.join("matrix.cells")).unwrap(),
            "sample_a\nsample_b\n",
            "{case}"
        );
        let sidecar = fs::read_to_string(out.join("matrix.sample.barcodes")).unwrap();
        let prefixes: Vec<_> = sidecar.lines().collect();
        let prefix_len = 32 - bc_len;
        assert_eq!(
            prefixes,
            vec![
                "A".repeat(prefix_len),
                format!("{}C", "A".repeat(prefix_len - 1))
            ],
            "{case}"
        );
        for (record, prefix) in records.iter().zip(prefixes) {
            assert_eq!(
                record.barcode >> (2 * bc_len),
                encode_bus_seq(prefix.as_bytes()),
                "{case}"
            );
            assert_eq!(record.ec, 0, "{case}");
            assert_eq!(record.count, 1, "{case}");
        }
        assert_eq!(
            json_u64(
                &fs::read_to_string(out.join("run_info.json")).unwrap(),
                "n_processed"
            ),
            Some(2),
            "{case}"
        );
    }
}

#[test]
fn bus_batch_resolves_relative_paths_from_batch_file_directory() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let data_dir = dir.path().join("data");
    fs::create_dir(&data_dir).expect("create data dir");
    let r1 = data_dir.join("sample_r1.fastq");
    let r2 = data_dir.join("sample_r2.fastq");
    let batch = data_dir.join("batch.txt");
    let out_dir = dir.path().join("batch_relative_out");

    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &transcript)]);
    fs::write(&batch, "sample_a\tsample_r1.fastq\tsample_r2.fastq\n").expect("write batch");

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
        .status()
        .expect("run kallistors bus relative batch");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\n"
    );
}

#[test]
fn bus_batch_rejects_wrong_number_of_files_per_row() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let batch = dir.path().join("batch_bad.txt");
    let out_dir = dir.path().join("batch_bad_out");

    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    fs::write(&batch, format!("sample_a\t{}\n", r1.display())).expect("write batch");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("-B")
        .arg(&batch)
        .output()
        .expect("run kallistors bus malformed batch");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("batch file line 1 has 1 files, expected 2")
    );
}

#[test]
fn bus_rejects_empty_batch_with_explicit_or_inferred_technology() {
    let dir = tempfile::tempdir().unwrap();
    let (index, _) = build_tiny_index(&dir);
    let batch = dir.path().join("empty.txt");
    fs::write(&batch, "\n# no samples\n\n").unwrap();
    for (case, technology) in [("explicit", Some("10XV3")), ("inferred", None)] {
        let mut command = bus_command(&index, &dir.path().join(case));
        if let Some(technology) = technology {
            command.args(["-x", technology]);
        }
        let output = command.arg("--batch").arg(&batch).output().unwrap();
        assert_bus_error(&output, case, "batch file contains no read groups");
    }
}

#[test]
fn bus_batch_barcodes_requires_batch_mode() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("batch_barcodes_without_batch_out");

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
        .arg("--batch-barcodes")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --batch-barcodes without batch");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("--batch-barcodes requires batch mode")
    );
}

#[test]
fn bus_batch_barcodes_rejects_unusable_prefix_geometry() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let batch = dir.path().join("batch.txt");
    write_fastq(&r1, &[("cell", &[b'A'; 44])]);
    write_fastq(&r2, &[("seq", &transcript)]);
    fs::write(
        &batch,
        format!("sample_a\t{}\t{}\n", r1.display(), r2.display()),
    )
    .unwrap();
    for (case, technology, error) in [
        (
            "no_room",
            "0,0,32:0,32,44:1,0,0",
            "--batch-barcodes requires barcode length shorter than 32 bases",
        ),
        (
            "unbounded",
            "0,0,0:0,8,20:1,0,0",
            "--batch-barcodes requires a bounded barcode length",
        ),
    ] {
        let out = dir.path().join(case);
        let output = bus_command(&index, &out)
            .args(["-x", technology, "--batch-barcodes", "-B"])
            .arg(&batch)
            .output()
            .unwrap();
        assert_bus_error(&output, case, error);
        assert!(
            !out.exists(),
            "{case}: invalid geometry must precede output creation"
        );
    }
}

#[test]
fn bus_batch_rejects_missing_read_file_clearly() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let batch = dir.path().join("batch.txt");
    let missing_r1 = dir.path().join("missing_r1.fastq");
    let missing_r2 = dir.path().join("missing_r2.fastq");
    let out_dir = dir.path().join("batch_missing_fastq_out");
    fs::write(
        &batch,
        format!(
            "sampleA\t{}\t{}\n",
            missing_r1.display(),
            missing_r2.display()
        ),
    )
    .expect("write batch file");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--batch")
        .arg(&batch)
        .output()
        .expect("run kallistors bus batch with missing FASTQs");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("file not found") && stderr.contains("missing_r1.fastq"),
        "{stderr}"
    );
    assert!(
        !out_dir.exists(),
        "missing batch inputs should fail before creating output directory"
    );
}

#[test]
fn bus_batch_reuses_barcode_for_repeated_sample_ids() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let a_r1 = dir.path().join("a_r1.fastq");
    let a_r2 = dir.path().join("a_r2.fastq");
    let b_r1 = dir.path().join("b_r1.fastq");
    let b_r2 = dir.path().join("b_r2.fastq");
    let c_r1 = dir.path().join("c_r1.fastq");
    let c_r2 = dir.path().join("c_r2.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("batch_repeated_ids_out");

    write_fastq(&a_r1, &[("a_cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&a_r2, &[("a_seq", &transcript)]);
    write_fastq(&b_r1, &[("b_cell", b"TGCATGCATGCATGCATTTTTTTTTTTT")]);
    write_fastq(&b_r2, &[("b_seq", &transcript)]);
    write_fastq(&c_r1, &[("c_cell", b"AAAACCCCGGGGTTTTAAAAAAAAAAAA")]);
    write_fastq(&c_r2, &[("c_seq", &transcript)]);
    fs::write(
        &batch,
        format!(
            "sample_a\t{}\t{}\nsample_b\t{}\t{}\nsample_a\t{}\t{}\n",
            a_r1.display(),
            a_r2.display(),
            b_r1.display(),
            b_r2.display(),
            c_r1.display(),
            c_r2.display()
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
        .arg("--batch-barcodes")
        .status()
        .expect("run kallistors bus batch with repeated sample IDs");
    assert!(status.success());

    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\nsample_b\nsample_a\n"
    );
    let sample_barcodes = fs::read_to_string(out_dir.join("matrix.sample.barcodes")).unwrap();
    let sample_barcodes = sample_barcodes.lines().collect::<Vec<_>>();
    assert_eq!(sample_barcodes.len(), 3);
    assert_ne!(sample_barcodes[0], sample_barcodes[1]);
    assert_eq!(sample_barcodes[0], sample_barcodes[2]);
}

#[test]
fn bus_batch_paired_writes_one_flens_line_per_batch() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_index_with_transcript(
        &dir,
        b"ACGTTGCAAGTCGATCGTACCGTTAACGGCATGTCAGTCCGATGACCTAGTGCATTCGAGTACCGATGCTAGTCGATCGT",
    );
    let a_r1 = dir.path().join("a_r1.fastq");
    let a_r2 = dir.path().join("a_r2.fastq");
    let b_r1 = dir.path().join("b_r1.fastq");
    let b_r2 = dir.path().join("b_r2.fastq");
    let batch = dir.path().join("batch.txt");

    write_fastq(&a_r1, &[("a/1", &transcript[..40]); 10001]);
    write_fastq(
        &a_r2,
        &[("a/2", reverse_complement(&transcript[20..60]).as_slice()); 10001],
    );
    write_fastq(&b_r1, &[("b/1", &transcript[..35]); 7]);
    write_fastq(
        &b_r2,
        &[("b/2", reverse_complement(&transcript[10..45]).as_slice()); 7],
    );
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

    for threads in [1, 8] {
        let out_dir = dir.path().join(format!("batch_paired_flens_out-{threads}"));
        let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
            .arg("bus")
            .args(["-t", &threads.to_string()])
            .arg("-i")
            .arg(&index)
            .arg("-o")
            .arg(&out_dir)
            .arg("-x=-1,0,0:-1,0,0:0,0,0,1,0,0")
            .arg("--paired")
            .arg("--batch")
            .arg(&batch)
            .status()
            .expect("run kallistors paired batch bus");
        assert!(status.success());

        let flens = fs::read_to_string(out_dir.join("flens.txt")).unwrap();
        let lines = flens.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 2);
        for (line, (expected_length, expected_count)) in
            lines.into_iter().zip([(60, 10001), (45, 7)])
        {
            let counts = line
                .split_whitespace()
                .map(|value| value.parse::<u32>().unwrap())
                .collect::<Vec<_>>();
            let mut expected = vec![0; kallistors::pseudoalign::MAX_FRAG_LEN as usize];
            expected[expected_length] = expected_count;
            assert_eq!(
                counts, expected,
                "batch with fragment length {expected_length}"
            );
        }
    }
}

#[test]
fn bus_batch_unmapped_writes_ratio_per_processed_read() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let a_r1 = dir.path().join("a_r1.fastq");
    let a_r2 = dir.path().join("a_r2.fastq");
    let b_r1 = dir.path().join("b_r1.fastq");
    let b_r2 = dir.path().join("b_r2.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("batch_unmapped_out");

    write_fastq(&a_r1, &[("a_cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&a_r2, &[("a_seq", &transcript)]);
    write_fastq(&b_r1, &[("b_cell", b"TGCATGCATGCATGCATTTTTTTTTTTT")]);
    write_fastq(&b_r2, &[("b_seq", &[b'N'; 80])]);
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
        .arg("--unmapped")
        .status()
        .expect("run kallistors bus batch --unmapped");
    assert!(status.success());

    let ratios = fs::read_to_string(out_dir.join("unmapped_ratio.txt")).unwrap();
    assert!(ratios.trim_end().ends_with(','));
    let values = ratios
        .trim_end()
        .split(',')
        .filter(|value| !value.is_empty())
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(values.len(), 2);
    assert!(values[0] < values[1]);
    assert_eq!(values[1], 1.0);
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\nsample_b\n"
    );
}

#[test]
fn bus_batch_num_reads_larger_than_input_returns_error_after_writing_outputs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("batch_num_reads_short_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);
    fs::write(
        &batch,
        format!("sample_a\t{}\t{}\n", r1.display(), r2.display()),
    )
    .expect("write batch");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--batch")
        .arg(&batch)
        .arg("--numReads")
        .arg("2")
        .output()
        .expect("run kallistors batch bus --numReads larger than input");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("Number of reads processed is less than --numReads: 2, returning 1")
    );
    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(1)
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\n"
    );
    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}
