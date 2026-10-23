//! Input regressions.

use super::*;
#[test]
fn bus_reads_gzipped_fastq_and_reports_verbose_inputs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq.gz");
    let r2 = dir.path().join("r2.fastq.gz");
    let out_dir = dir.path().join("bus_gz_out");
    let (index, transcript) = build_tiny_index(&dir);

    write_fastq_gz(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq_gz(&r2, &[("seq_read", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("--verbose")
        .arg("--dfk-onlist")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus with gzipped FASTQ");
    assert_bus_success(&output, "verbose input");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("will process sample 1:"), "{stderr}");

    assert!(stderr.contains(&r1.display().to_string()), "{stderr}");

    assert!(stderr.contains(&r2.display().to_string()), "{stderr}");

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 12);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].ec, 0);
}

#[test]
fn bus_reads_single_fastq_from_stdin_sentinel() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let out_dir = dir.path().join("stdin_bus_out");
    let mut fastq = Vec::new();
    fastq.extend_from_slice(b"@stdin_read\n");
    fastq.extend_from_slice(b"ACGTTTAA");
    fastq.extend_from_slice(&transcript);
    fastq.extend_from_slice(b"\n+\n");
    fastq.extend(std::iter::repeat_n(b'I', 8 + transcript.len()));
    fastq.extend_from_slice(b"\n");

    let mut child = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=0,0,4:0,4,8:0,8,0")
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .expect("spawn kallistors bus stdin");
    child
        .stdin
        .as_mut()
        .expect("child stdin")
        .write_all(&fastq)
        .expect("write FASTQ to stdin");
    let status = child.wait().expect("wait for kallistors bus stdin");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 4);
    assert_eq!(umi_len, 4);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].ec, 0);
    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(1)
    );
}

#[test]
fn bus_rejects_zero_threads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("zero_threads_out");

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
        .arg("-t")
        .arg("0")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus -t 0");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid number of threads 0"));
}

#[test]
fn bus_rejects_missing_index_file_clearly() {
    let dir = tempfile::tempdir().expect("tempdir");
    let index = dir.path().join("missing.idx");
    let read = dir.path().join("read.fastq");
    let out_dir = dir.path().join("missing_index_out");

    write_fastq(&read, &[("read", b"ACGTACGTACGTACGTACGTACGTACGTACGT")]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("MATQSEQ")
        .arg(&read)
        .output()
        .expect("run kallistors bus with missing index");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("kallisto index file not found"), "{stderr}");
    assert!(
        !out_dir.exists(),
        "missing index should fail before creating output directory"
    );
}

#[test]
fn bus_rejects_output_path_that_is_file() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_path = dir.path().join("not_a_directory");

    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &transcript)]);
    fs::write(&out_path, b"occupied").expect("write output file placeholder");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_path)
        .arg("-x")
        .arg("10XV3")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus with file output path");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("exists and is not a directory"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn bus_rejects_missing_fastq_file_clearly() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("missing_r1.fastq");
    let r2 = dir.path().join("missing_r2.fastq");
    let out_dir = dir.path().join("missing_fastq_out");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus with missing FASTQs");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("file not found") && stderr.contains("missing_r1.fastq"),
        "{stderr}"
    );
    assert!(
        !out_dir.exists(),
        "missing inputs should fail before creating output directory"
    );
}

#[test]
fn bus_num_writes_read_numbers_in_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("num_out");

    write_fastq(
        &r1,
        &[
            ("cell_read", b"NCGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_read2", b"ANGTACGTACGTACGTTTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read", &transcript), ("seq_read2", &transcript)],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--num")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --num");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].flags, 0);
    assert_eq!(records[1].flags, 1);
}

#[test]
fn bus_rejects_fastq_files_with_different_record_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    for (case, left_count, right_count) in [("first-longer", 2, 1), ("second-longer", 1, 2)] {
        let barcodes: &[(&str, &[u8])] = &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_read2", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ];
        let sequences: &[(&str, &[u8])] = &[("seq_read", &transcript), ("seq_read2", &transcript)];
        write_fastq(&r1, &barcodes[..left_count]);
        write_fastq(&r2, &sequences[..right_count]);
        for threads in [1, 4] {
            let out = dir.path().join(format!("{case}-{threads}"));
            let error = dir.path().join(format!("{case}-{threads}.stderr"));
            let mut command = bus_command(&index, &out);
            command
                .args(["-x", "10xv3", "-t", &threads.to_string()])
                .arg(&r1)
                .arg(&r2)
                .stderr(Stdio::from(fs::File::create(&error).unwrap()));
            let status =
                command_status_with_timeout(&mut command, Duration::from_secs(30)).unwrap();
            assert!(!status.success(), "{case}/{threads}");
            assert!(
                fs::read_to_string(error)
                    .unwrap()
                    .contains("FASTQ files ended at different record counts"),
                "{case}/{threads}"
            );
        }
    }
}

#[test]
fn bus_rejects_incomplete_interleaved_record_group() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let interleaved = dir.path().join("interleaved_incomplete.fastq");
    let out_dir = dir.path().join("interleaved_incomplete_out");

    write_fastq(
        &interleaved,
        &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")],
    );

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--inleaved")
        .arg(&interleaved)
        .output()
        .expect("run kallistors bus incomplete interleaved");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("interleaved FASTQ ended in the middle of a record group")
    );
}

#[test]
fn bus_num_reads_zero_means_unlimited() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("zero_out");

    write_fastq(
        &r1,
        &[
            ("cell_read_1", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_read_2", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read_1", &transcript), ("seq_read_2", &transcript)],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--numReads")
        .arg("0")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --numReads 0");
    assert!(status.success());

    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(2)
    );
    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 2);
}

#[test]
fn bus_num_reads_larger_than_input_returns_error_after_writing_outputs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("num_reads_short_out");

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
        .arg("--numReads")
        .arg("2")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --numReads larger than input");
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
    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}
