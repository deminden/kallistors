//! Layouts regressions.

use super::*;
#[test]
fn bus_cli_writes_tenx_v3_bus_outputs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("bus_out");

    let (index, transcript) = build_tiny_index(&dir);

    write_fastq(
        &r1,
        &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_read2", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read", &transcript), ("seq_read2", &transcript)],
    );

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("-N")
        .arg("1")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus");
    assert!(output.status.success());

    assert_eq!(
        fs::read_to_string(out_dir.join("transcripts.txt")).unwrap(),
        "tx0\n"
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.ec")).unwrap(),
        "0\t0\n"
    );
    let bus = fs::read(out_dir.join("output.bus")).expect("read bus");
    assert_eq!(&bus[0..4], b"BUS\0");
    assert_eq!(read_u32_le(&bus, 4), 1);
    assert_eq!(read_u32_le(&bus, 8), 16);
    assert_eq!(read_u32_le(&bus, 12), 12);
    let text_len = read_u32_le(&bus, 16) as usize;
    assert_eq!(
        std::str::from_utf8(&bus[20..20 + text_len]).unwrap(),
        "BUS file produced by kallisto"
    );
    let record_offset = 20 + text_len;
    assert_eq!(bus.len(), record_offset + 32);
    assert_ne!(read_u64_le(&bus, record_offset), 0);
    assert_ne!(read_u64_le(&bus, record_offset + 8), 0);
    assert_eq!(read_i32_le(&bus, record_offset + 16), 0);
    assert_eq!(read_u32_le(&bus, record_offset + 20), 1);
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).unwrap();
    for (field, expected) in [
        ("n_targets", 1),
        ("n_bootstraps", 0),
        ("n_processed", 1),
        ("n_pseudoaligned", 1),
        ("n_unique", 1),
        ("index_version", 13),
        ("k-mer length", 31),
    ] {
        assert_eq!(json_u64(&run_info, field), Some(expected), "{field}");
    }
    let json: serde_json::Value = serde_json::from_str(&run_info).unwrap();
    assert_eq!(json["p_pseudoaligned"].as_f64(), Some(100.0));
    assert_eq!(json["p_unique"].as_f64(), Some(100.0));
    assert!(json["start_time"].is_string());
    assert!(json["call"].as_str().unwrap().contains("bus"));
    assert!(json.get("technology").is_none(), "upstream BUS schema");
}

#[test]
fn bus_without_technology_infers_paired_interleaved_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let reads = dir.path().join("interleaved.fastq");
    let out_dir = dir.path().join("bulk_no_technology_interleaved_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(
        &reads,
        &[("left", &transcript[..40]), ("right", &transcript[20..60])],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("--interleaved")
        .arg(&reads)
        .status()
        .expect("run inferred paired interleaved kallistors bus without -x");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert!(out_dir.join("flens.txt").exists());
    assert!(out_dir.join("index.saved").exists());
}

#[test]
#[ignore = "requires original kallisto on PATH"]
fn bus_tenx_v3_matches_upstream_kallisto_on_tiny_case() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let kallistors_out = dir.path().join("kallistors_bus");
    let kallisto_out = dir.path().join("kallisto_bus");
    let upstream_index = dir.path().join("upstream.idx");
    let fasta = write_transcript_fasta(&dir, "tx0", TINY_TRANSCRIPT);

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read instrument:1", TINY_TRANSCRIPT)]);

    let status = command_status_with_timeout(
        Command::new("kallisto")
            .arg("index")
            .arg("-i")
            .arg(&upstream_index)
            .arg(&fasta),
        Duration::from_secs(30),
    )
    .expect("run kallisto index");
    assert!(status.success(), "kallisto index failed");

    let status = command_status_with_timeout(
        Command::new(env!("CARGO_BIN_EXE_kallistors"))
            .arg("bus")
            .arg("-i")
            .arg(&upstream_index)
            .arg("-o")
            .arg(&kallistors_out)
            .arg("-x")
            .arg("10XV3")
            .arg(&r1)
            .arg(&r2),
        Duration::from_secs(30),
    )
    .expect("run kallistors bus");
    assert!(status.success());

    let status = command_status_with_timeout(
        Command::new("kallisto")
            .arg("bus")
            .arg("-i")
            .arg(&upstream_index)
            .arg("-o")
            .arg(&kallisto_out)
            .arg("-x")
            .arg("10XV3")
            .arg(&r1)
            .arg(&r2),
        Duration::from_secs(30),
    )
    .expect("run kallisto bus");
    assert!(status.success(), "kallisto bus failed");

    assert_eq!(
        fs::read_to_string(kallistors_out.join("transcripts.txt")).unwrap(),
        fs::read_to_string(kallisto_out.join("transcripts.txt")).unwrap()
    );
    assert_eq!(
        fs::read_to_string(kallistors_out.join("matrix.ec")).unwrap(),
        fs::read_to_string(kallisto_out.join("matrix.ec")).unwrap()
    );
    assert_eq!(
        read_bus_records(&kallistors_out.join("output.bus")),
        read_bus_records(&kallisto_out.join("output.bus"))
    );
}

#[test]
fn bus_custom_technology_suffix_paired_consumes_extra_mate() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let r3 = dir.path().join("r3.fastq");
    let out_dir = dir.path().join("custom_suffix_paired_out");
    let (index, transcript) = build_tiny_index(&dir);

    write_fastq(&r1, &[("barcodes", b"ACGTACGTACGTACGTAAAAAAAAAAAA")]);
    write_fastq(&r2, &[("left", &transcript)]);
    write_fastq(&r3, &[("right", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("0,0,16:0,16,28:1,0,0%FORWARD%PAIRED")
        .arg(&r1)
        .arg(&r2)
        .arg(&r3)
        .status()
        .expect("run kallistors bus custom suffix paired");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 12);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].ec, 0);
}

#[test]
fn bus_custom_all_sentinel_technology_does_not_panic() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let read = dir.path().join("read.fastq");
    let out_dir = dir.path().join("custom_all_sentinel_out");

    write_fastq(&read, &[("read", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=-1,-1,-1:-1,-1,-1:-1,-1,-1")
        .arg(&read)
        .output()
        .expect("run kallistors bus with all-sentinel custom technology");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("zero reads pseudoaligned"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert!(records.is_empty());
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).unwrap();
    assert_eq!(json_u64(&run_info, "n_processed"), Some(1));
    assert_eq!(json_u64(&run_info, "n_pseudoaligned"), Some(0));
}

#[test]
fn bus_interleaved_reads_one_fastq_as_technology_groups() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let interleaved = dir.path().join("interleaved.fastq");
    let out_dir = dir.path().join("interleaved_out");

    write_fastq(
        &interleaved,
        &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("seq_read", &transcript),
        ],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--inleaved")
        .arg(&interleaved)
        .status()
        .expect("run kallistors bus interleaved");
    assert!(status.success());

    let bus = fs::read(out_dir.join("output.bus")).expect("read bus");
    let text_len = read_u32_le(&bus, 16) as usize;
    assert_eq!(bus.len(), 20 + text_len + 32);
    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(1)
    );
}

#[test]
fn bus_stranded_rejects_unpaired_multi_sequence_technology() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("unpaired_multi_stranded_out");

    write_fastq(&r1, &[("cell", b"ACGTACGT")]);
    write_fastq(&r2, &[("seq", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=-1,-1,-1:0,0,4:1,0,31,1,31,62")
        .arg("--fr-stranded")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus unpaired multi-sequence --fr-stranded");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains(
        "Strand-specific read processing is only supported for technologies with a single cDNA read file or paired-end reads"
    ));
}
