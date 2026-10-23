//! Bam input regressions.

use super::*;
#[test]
fn bus_bam_rejects_barcode_tag_longer_than_bus_encoding() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_long_barcode_out");
    let long_barcode = vec![b'A'; 33];

    write_bam(&bam, &[(&transcript, &long_barcode, b"TTTTTTTTTT")]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--bam")
        .arg(&bam)
        .output()
        .expect("run kallistors bus --bam with long barcode");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("barcode length 33 exceeds BUS limit of 32 bases"),
        "{stderr}"
    );
}

#[test]
fn bus_bam_rejects_paired_flag_in_short_read_mode() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_paired_flag_out");

    write_bam(&bam, &[(&transcript, b"ACGTACGTACGTACGT", b"TTTTTTTTTT")]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg("--paired")
        .arg(&bam)
        .output()
        .expect("run kallistors bus --bam --paired");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("Paired reads are not compatible with the specified technology")
    );
    assert!(!out_dir.join("flens.txt").exists());
}

#[test]
fn bus_bam_reads_and_normalizes_corrected_barcode_and_umi_tags() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    for (case, tag, umi, expected_umi) in [
        (
            "raw_umi",
            Tag::UMI_SEQUENCE,
            &b"TTTTTTTTTT"[..],
            &b"TTTTTTTTTT"[..],
        ),
        ("mi", Tag::UMI_ID, &b"TTTTAAAA"[..], &b"TTTTAAAA"[..]),
        (
            "mi_suffix",
            Tag::UMI_ID,
            &b"TTTTAAAA-1"[..],
            &b"TTTTAAAA"[..],
        ),
        (
            "ub",
            Tag::new(b'U', b'B'),
            &b"AAAATTTT"[..],
            &b"AAAATTTT"[..],
        ),
        (
            "ub_suffix",
            Tag::new(b'U', b'B'),
            &b"AAAATTTT-1"[..],
            &b"AAAATTTT"[..],
        ),
    ] {
        write_bam_with_corrected_barcode_and_umi(&bam, &transcript, b"ACGTACGT-1", umi, tag);
        let out = dir.path().join(case);
        let output = bus_command(&index, &out)
            .args(["-x", "10XV2", "--bam"])
            .arg(&bam)
            .output()
            .unwrap();
        assert_bus_success(&output, case);
        let (bc_len, umi_len, records) = read_bus_records(&out.join("output.bus"));
        assert_eq!((bc_len, umi_len), (8, expected_umi.len() as u32), "{case}");
        assert_eq!(records.len(), 1, "{case}");
        assert_eq!(records[0].barcode, encode_bus_seq(b"ACGTACGT"), "{case}");
        assert_eq!(records[0].umi, encode_bus_seq(expected_umi), "{case}");
        assert_eq!(records[0].ec, 0, "{case}");
    }
}

#[test]
fn bus_bam_prefers_raw_barcode_over_corrected_barcode() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("raw_and_corrected_barcode.bam");
    let out_dir = dir.path().join("bam_raw_over_corrected_barcode_out");

    write_bam_with_raw_and_corrected_barcode(
        &bam,
        &transcript,
        b"TTTTCCCC",
        b"ACGTACGT-1",
        b"AAAAAAAAAA",
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg(&bam)
        .status()
        .expect("run kallistors bus --bam raw and corrected barcode");
    assert!(status.success());

    let (bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"TTTTCCCC"));
}

#[test]
fn bus_bam_missing_tags_preserve_headers_counts_and_unmapped_ratio() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    for (case, barcode, umi, unmapped) in [
        ("missing_barcode", None, Some(&b"TTTTTTTTTT"[..]), false),
        ("missing_umi", Some(&b"ACGTACGTACGTACGT"[..]), None, false),
        (
            "missing_barcode_ratio",
            None,
            Some(&b"TTTTTTTTTT"[..]),
            true,
        ),
    ] {
        write_bam_with_optional_tags(&bam, &[(&transcript, barcode, umi)]);
        let out = dir.path().join(case);
        let mut command = bus_command(&index, &out);
        command.args(["-x", "10XV2", "--bam"]);
        if unmapped {
            command.arg("--unmapped");
        }
        let output = command.arg(&bam).output().unwrap();
        assert_bus_error(&output, case, "zero reads pseudoaligned");
        let (bc_len, umi_len, records) = read_bus_records(&out.join("output.bus"));
        assert_eq!((bc_len, umi_len), (16, 10), "{case}");
        assert!(records.is_empty(), "{case}");
        let info = fs::read_to_string(out.join("run_info.json")).unwrap();
        assert_eq!(json_u64(&info, "n_processed"), Some(1), "{case}");
        assert_eq!(json_u64(&info, "n_pseudoaligned"), Some(0), "{case}");
        if unmapped {
            let ratios = fs::read_to_string(out.join("unmapped_ratio.txt")).unwrap();
            let values: Vec<f64> = ratios
                .trim_end()
                .split(',')
                .map(|v| v.parse().unwrap())
                .collect();
            assert_eq!(values.len(), 1, "{case}");
            assert!(
                values[0] < 1.0,
                "{case}: tags must not erase the read's mapped k-mers"
            );
        }
    }
}

#[test]
fn bus_bam_long_unmapped_writes_novel_fastq_for_skipped_missing_tag_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("missing_barcode_long_unmapped.bam");
    let out_dir = dir.path().join("bam_missing_barcode_long_unmapped_out");
    let sequence = [b'N'; 80];

    write_bam_with_optional_tags(&bam, &[(&sequence, None, Some(b"TTTTTTTTTT"))]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg("--long")
        .arg("--threshold")
        .arg("0.5")
        .arg("--unmapped")
        .arg(&bam)
        .output()
        .expect("run kallistors bus --bam --long --unmapped missing barcode");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("zero reads pseudoaligned"), "{stderr}");

    let novel = fs::read_to_string(out_dir.join("novel.fastq")).unwrap();
    let novel_lines = novel.lines().collect::<Vec<_>>();
    assert_eq!(novel_lines.len(), 4);
    assert_eq!(novel_lines[0], "@skipped_missingTags");
    assert_eq!(novel_lines[1], std::str::from_utf8(&sequence).unwrap());
    assert_eq!(novel_lines[2], "+");
    assert_eq!(novel_lines[3], "!".repeat(80));
}

#[test]
fn bus_bam_skips_secondary_and_supplementary_records_before_counting_num_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_secondary_out");

    write_bam_with_flags(
        &bam,
        &[
            (
                &transcript,
                b"ACGTACGTACGTACGT",
                b"TTTTTTTTTT",
                Flags::SECONDARY,
            ),
            (
                &transcript,
                b"AAAACCCCGGGGTTTT",
                b"CCCCCCCCCC",
                Flags::SUPPLEMENTARY,
            ),
            (
                &transcript,
                b"TGCATGCATGCATGCA",
                b"AAAAAAAAAA",
                Flags::empty(),
            ),
        ],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg("-N")
        .arg("1")
        .arg(&bam)
        .status()
        .expect("run kallistors bus --bam with secondary/supplementary records");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 10);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"TGCATGCATGCATGCA"));
}

#[test]
fn bus_bam_long_writes_novel_fastq_for_filtered_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_long_out");

    write_bam(
        &bam,
        &[
            (&transcript, b"ACGTACGTACGTACGT", b"TTTTTTTTTT"),
            (&[b'N'; 80], b"TGCATGCATGCATGCA", b"AAAAAAAAAA"),
        ],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg("--long")
        .arg("--threshold")
        .arg("0.5")
        .arg("--unmapped")
        .arg(&bam)
        .status()
        .expect("run kallistors bus --bam --long");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    let ratios = fs::read_to_string(out_dir.join("unmapped_ratio.txt")).unwrap();
    assert_eq!(ratios.trim_end().split(',').count(), 2);
    let novel = fs::read_to_string(out_dir.join("novel.fastq")).unwrap();
    let novel_lines = novel.lines().collect::<Vec<_>>();
    assert_eq!(novel_lines.len(), 8);
    assert_eq!(novel_lines[0], "@unmapped");
    assert_eq!(novel_lines[1], std::str::from_utf8(&[b'N'; 80]).unwrap());
    assert_eq!(novel_lines[2], "+");
    assert_eq!(novel_lines[3], "!".repeat(80));
    assert_eq!(novel_lines[4], "@novel_disjointIntersect");
    assert_eq!(novel_lines[5], std::str::from_utf8(&[b'N'; 80]).unwrap());
    assert_eq!(novel_lines[6], "+");
    assert_eq!(novel_lines[7], "!".repeat(80));
}

#[test]
fn bus_bam_long_encodes_novel_fastq_quality_scores() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_long_quality_out");
    let sequence = [b'N'; 80];

    write_bam_with_quality(
        &bam,
        &sequence,
        b"TGCATGCATGCATGCA",
        b"AAAAAAAAAA",
        &[40; 80],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg("--long")
        .arg("--threshold")
        .arg("0.5")
        .arg("--unmapped")
        .arg(&bam)
        .status()
        .expect("run kallistors bus --bam --long with BAM qualities");
    assert!(!status.success());

    let novel = fs::read_to_string(out_dir.join("novel.fastq")).unwrap();
    let novel_lines = novel.lines().collect::<Vec<_>>();
    assert_eq!(novel_lines[0], "@unmapped");
    assert_eq!(novel_lines[2], "+");
    assert_eq!(novel_lines[3], "I".repeat(80));
}

#[test]
fn bus_bam_preserves_raw_tag_values_lengths_and_ecs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");

    let fixtures: [(&[u8], &[u8], &[u8]); 3] = [
        (&transcript, b"ACGTACGT", b"TTTTTT"),
        (&transcript, b"TGCATGCA", b"AAAAAA"),
        (&transcript, b"ACGTACGTAC", b"TTTTTTTT"),
    ];
    let input = fixtures.into_iter().cycle().take(4099).collect::<Vec<_>>();
    write_bam(&bam, &input);

    for threads in [1, 8] {
        let out_dir = dir.path().join(format!("bam_lengths_out-{threads}"));
        let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
            .arg("bus")
            .args(["-t", &threads.to_string()])
            .arg("--verbose")
            .arg("-i")
            .arg(&index)
            .arg("-o")
            .arg(&out_dir)
            .arg("-x")
            .arg("10XV3")
            .arg("--bam")
            .arg(&bam)
            .output()
            .expect("run kallistors bus --bam");
        assert_bus_success(&output, "verbose input");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("will process sample 1:"), "{stderr}");

        assert!(stderr.contains(&bam.display().to_string()), "{stderr}");

        let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
        assert_eq!(bc_len, 8);
        assert_eq!(umi_len, 6);
        assert_eq!(records.len(), 4099);
        for (record, (barcode, umi)) in records.iter().zip(
            [
                (&b"ACGTACGT"[..], &b"TTTTTT"[..]),
                (&b"TGCATGCA"[..], &b"AAAAAA"[..]),
                (&b"ACGTACGTAC"[..], &b"TTTTTTTT"[..]),
            ]
            .into_iter()
            .cycle(),
        ) {
            assert_eq!(record.barcode, encode_bus_seq(barcode));
            assert_eq!(record.umi, encode_bus_seq(umi));
            assert_eq!(record.ec, 0);
            assert_eq!(record.count, 1);
        }
    }
}

#[test]
fn bus_bam_ignores_num_and_keeps_ambiguity_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_num_out");

    write_bam(&bam, &[(&transcript, b"NCGTACGT", b"TTNTTT")]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg("--num")
        .arg(&bam)
        .output()
        .expect("run kallistors bus --bam --num");
    assert!(output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("Warning: --bam option was used, so --num option will be ignored")
    );

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].flags, 1 | (9 << 8));
}
