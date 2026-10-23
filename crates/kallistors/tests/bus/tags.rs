//! Tags regressions.

use super::*;
#[test]
fn bus_skips_short_barcode_and_umi_slices_without_losing_following_reads() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    write_fastq(&r2, &[("short", &transcript), ("valid", &transcript)]);
    for (case, short) in [
        ("barcode", &b"ACGTACGT"[..]),
        ("umi", &b"ACGTACGTACGTACGTTTT"[..]),
    ] {
        write_fastq(
            &r1,
            &[("short", short), ("valid", b"ACGTACGTACGTACGTTTTTTTTTTTTT")],
        );
        let out = dir.path().join(case);
        let output = bus_command(&index, &out)
            .args(["-x", "10XV3"])
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert_bus_success(&output, case);
        let (bc_len, umi_len, records) = read_bus_records(&out.join("output.bus"));
        assert_eq!((bc_len, umi_len), (16, 12), "{case}");
        assert_eq!(records.len(), 1, "{case}");
        assert_eq!(
            records[0].barcode,
            encode_bus_seq(b"ACGTACGTACGTACGT"),
            "{case}"
        );
        assert_eq!(records[0].umi, encode_bus_seq(b"TTTTTTTTTTTT"), "{case}");
        assert_eq!(records[0].ec, 0, "{case}");
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
fn bus_num_reads_counts_skipped_too_short_barcode_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("num_reads_short_barcode_out");

    write_fastq(
        &r1,
        &[
            ("short_cell", b"ACGTACGT"),
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("short_seq", &transcript), ("seq_read", &transcript)],
    );

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--numReads")
        .arg("1")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --numReads with too-short barcode");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("zero reads pseudoaligned"), "{stderr}");

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 12);
    assert!(records.is_empty());
    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(1)
    );
}

#[test]
fn bus_custom_open_ended_umi_patches_bus_header_length() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("custom_open_ended_umi_out");

    let mut bc_umi = b"ACGTACGTACGTACGT".to_vec();
    bc_umi.extend(std::iter::repeat_n(b'T', 12));
    write_fastq(&r1, &[("bc_umi", &bc_umi)]);
    write_fastq(&r2, &[("seq", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("0,0,16:0,16,0:1,0,0")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus custom open-ended UMI");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 12);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].ec, 0);
}

#[test]
fn bus_rejects_custom_barcode_and_umi_widths_above_32_bases() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    write_fastq(&r1, &[("cell", &[b'A'; 64])]);
    write_fastq(&r2, &[("seq", &transcript)]);
    for (case, technology, error) in [
        (
            "barcode",
            "0,0,33:0,33,45:1,0,0",
            "barcode length 33 exceeds BUS limit of 32 bases",
        ),
        (
            "umi",
            "0,0,16:0,16,49:1,0,0",
            "UMI length 33 exceeds BUS limit of 32 bases",
        ),
    ] {
        let out = dir.path().join(case);
        let output = bus_command(&index, &out)
            .args(["-x", technology])
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert_bus_error(&output, case, error);
    }
}

#[test]
fn bus_rejects_invalid_tag_sequences() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    write_fastq(&r1, &[("cell", &[b'A'; 56])]);
    write_fastq(&r2, &[("seq", &transcript)]);
    for (case, tag, error) in [
        (
            "too_long",
            "A".repeat(33),
            "tag sequence length 33 exceeds BUS limit of 32 bases",
        ),
        (
            "invalid_base",
            "ACGN".into(),
            "tag sequence contains invalid base N at position 4",
        ),
        ("empty", "".into(), "tag sequence cannot be empty"),
    ] {
        let out = dir.path().join(case);
        let output = bus_command(&index, &out)
            .args(["-x", "0,0,16:0,16,56:1,0,0", "--tag"])
            .arg(tag)
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert_bus_error(&output, case, error);
        assert!(
            !out.exists(),
            "{case}: validation must precede output creation"
        );
    }
}

#[test]
fn bus_unmapped_writes_ratio_for_skipped_too_short_barcode_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("unmapped_short_barcode_out");

    write_fastq(
        &r1,
        &[
            ("short_cell", b"ACGTACGT"),
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("short_seq", &[b'N'; 80]), ("seq_read", &transcript)],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--unmapped")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --unmapped with too-short barcode");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    let ratios = fs::read_to_string(out_dir.join("unmapped_ratio.txt")).unwrap();
    let values = ratios
        .trim_end()
        .split(',')
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(values.len(), 2);
    assert_eq!(values[0], 1.0);
    assert!(values[1] < values[0]);
}

#[test]
fn bus_custom_rx_umi_reads_fastq_header_comment() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("rx_bus_out");

    write_fastq(
        &r1,
        &[
            ("cell RX:Z:ACGT sample", b"TGCA"),
            ("cell_missing_rx sample", b"ACGT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read", &transcript), ("seq_missing_rx", &transcript)],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("0,0,4:RX:1,0,0")
        .arg("--pseudobam")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus custom RX");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 4);
    assert_eq!(umi_len, 4);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"TGCA"));
    assert_eq!(records[0].umi, encode_bus_seq(b"ACGT"));
    assert_eq!(
        json_u64(
            &fs::read_to_string(out_dir.join("run_info.json")).unwrap(),
            "n_processed"
        ),
        Some(2)
    );

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open pseudobam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read pseudobam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read pseudobam records");
    assert_eq!(records.len(), 1);
    assert_eq!(
        records[0].name().map(|name| name.as_ref()),
        Some(&b"seq_read"[..])
    );
}

#[test]
fn bus_smartseq3_default_tag_trims_umi() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bc1 = dir.path().join("bc1.fastq");
    let bc2 = dir.path().join("bc2.fastq");
    let r1 = dir.path().join("smartseq3_r1.fastq");
    let r2 = dir.path().join("smartseq3_r2.fastq");

    write_fastq(&bc1, &[("bc1", b"ACGT")]);
    write_fastq(&bc2, &[("bc2", b"TGCA")]);
    let mut tagged_read = b"ATTGCGCAATG".to_vec();
    tagged_read.extend_from_slice(b"TTTTTTTT");
    tagged_read.extend_from_slice(b"AAA");
    tagged_read.extend_from_slice(&transcript);
    write_fastq(&r1, &[("r1", &tagged_read)]);
    write_fastq(&r2, &[("r2", &reverse_complement(&transcript))]);

    for threads in [1, 8] {
        let out_dir = dir.path().join(format!("smartseq3_out-{threads}"));
        let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
            .arg("bus")
            .args(["-t", &threads.to_string()])
            .arg("-i")
            .arg(&index)
            .arg("-o")
            .arg(&out_dir)
            .arg("-x")
            .arg("SMARTSEQ3")
            .arg(&bc1)
            .arg(&bc2)
            .arg(&r1)
            .arg(&r2)
            .output()
            .expect("run kallistors bus smartseq3");
        assert!(output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr)
                .contains("[bus] Using ATTGCGCAATG as UMI tag sequence")
        );

        let bus = fs::read(out_dir.join("output.bus")).expect("read bus");
        assert_eq!(read_u32_le(&bus, 8), 8);
        assert_eq!(read_u32_le(&bus, 12), 8);
        let text_len = read_u32_le(&bus, 16) as usize;
        let record_offset = 20 + text_len;
        assert_eq!(bus.len(), record_offset + 32);
        assert_ne!(read_u64_le(&bus, record_offset + 8), u64::MAX);
        assert_eq!(read_i32_le(&bus, record_offset + 16), 0);
        assert_eq!(
            fs::read(out_dir.join("index.saved")).unwrap(),
            fs::read(&index).unwrap()
        );
        let flens = fs::read_to_string(out_dir.join("flens.txt")).unwrap();
        let flen_values = flens
            .split_whitespace()
            .map(|value| value.parse::<u32>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(flen_values.len(), 1000);
        assert_eq!(flen_values.iter().sum::<u32>(), 1);
    }
}

#[test]
fn bus_no_umi_technology_writes_saved_index() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let read = dir.path().join("matqseq.fastq");
    let out_dir = dir.path().join("matqseq_out");
    let mut matqseq_read = b"ACGTACGT".to_vec();
    matqseq_read.extend_from_slice(&transcript);
    write_fastq(&read, &[("matqseq", &matqseq_read)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("MATQSEQ")
        .arg(&read)
        .status()
        .expect("run kallistors bus MATQSEQ");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert_eq!(
        fs::read(out_dir.join("index.saved")).unwrap(),
        fs::read(&index).unwrap()
    );
}
