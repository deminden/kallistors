//! Mapping output regressions.

use super::*;
#[test]
fn bus_rejects_disjoint_transcript_hits_without_emitting_an_empty_ec() {
    let dir = tempfile::tempdir().expect("tempdir");
    let other = b"GTTCCATAGGACCTTGCGTATCAGTTCGGATACCTGACTTACCGGTATGGTCATGACCTTAGGACCTCGTACGAGTTCAG";
    let fasta = dir.path().join("transcripts.fa");
    let index = dir.path().join("transcripts.idx");
    fs::write(
        &fasta,
        format!(
            ">tx0\n{}\n>tx1\n{}\n",
            std::str::from_utf8(TINY_TRANSCRIPT).unwrap(),
            std::str::from_utf8(other).unwrap()
        ),
    )
    .unwrap();
    build_index(&index, &[fasta], IndexBuildOptions::default()).unwrap();
    let chimera = [TINY_TRANSCRIPT[..40].to_vec(), other[..40].to_vec()].concat();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    write_fastq(
        &r1,
        &[
            ("valid", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("disjoint", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(&r2, &[("valid", TINY_TRANSCRIPT), ("disjoint", &chimera)]);
    let out = dir.path().join("bus");
    let output = bus_command(&index, &out)
        .args(["-x", "10xv3"])
        .args([r1, r2])
        .output()
        .unwrap();
    assert_bus_success(&output, "disjoint transcript hits");
    let (_, _, records) = read_bus_records(&out.join("output.bus"));
    assert_eq!(
        records.len(),
        1,
        "disjoint hits must not create a BUS record"
    );
    assert_eq!(records[0].barcode, encode_bus_seq(b"ACGTACGTACGTACGT"));
    let info = fs::read_to_string(out.join("run_info.json")).unwrap();
    assert_eq!(json_u64(&info, "n_processed"), Some(2));
    assert_eq!(json_u64(&info, "n_pseudoaligned"), Some(1));
    for line in fs::read_to_string(out.join("matrix.ec")).unwrap().lines() {
        assert!(!line.split_once('\t').unwrap().1.is_empty());
    }
}

#[test]
fn bus_explicit_and_inferred_bulk_share_single_end_outputs() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let read = dir.path().join("bulk.fastq");
    write_fastq(&read, &[("bulk_read", &transcript)]);
    for (case, technology) in [("explicit", Some("Bulk")), ("inferred", None)] {
        let out = dir.path().join(case);
        let mut command = bus_command(&index, &out);
        if let Some(technology) = technology {
            command.args(["-x", technology]);
        }
        let output = command.arg(&read).output().unwrap();
        assert_bus_success(&output, case);
        let (bc_len, umi_len, records) = read_bus_records(&out.join("output.bus"));
        assert_eq!((bc_len, umi_len), (16, 1), "{case}");
        assert_eq!(records.len(), 1, "{case}");
        assert_eq!(
            (records[0].barcode, records[0].umi, records[0].ec),
            (0, 0, 0),
            "{case}"
        );
        assert_eq!(
            fs::read(out.join("index.saved")).unwrap(),
            fs::read(&index).unwrap(),
            "{case}"
        );
        assert_eq!(
            fs::read_to_string(out.join("matrix.cells")).unwrap(),
            "batch0\n",
            "{case}"
        );
        assert_eq!(
            fs::read_to_string(out.join("matrix.sample.barcodes")).unwrap(),
            "AAAAAAAAAAAAAAAA\n",
            "{case}"
        );
    }
}

#[test]
fn bus_infers_paired_direct_reads_with_or_without_paired_flag() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let left = dir.path().join("left.fastq");
    let right = dir.path().join("right.fastq");
    write_fastq(&left, &[("left", &transcript[..40])]);
    write_fastq(&right, &[("right", &transcript[20..60])]);
    for (case, paired) in [("explicit_flag", true), ("inferred", false)] {
        let out = dir.path().join(case);
        let mut command = bus_command(&index, &out);
        if paired {
            command.arg("--paired");
        }
        let output = command.arg(&left).arg(&right).output().unwrap();
        assert_bus_success(&output, case);
        let (bc_len, umi_len, records) = read_bus_records(&out.join("output.bus"));
        assert_eq!((bc_len, umi_len), (16, 1), "{case}");
        assert_eq!(records.len(), 1, "{case}");
        assert_eq!(records[0].ec, 0, "{case}");
        assert!(out.join("flens.txt").exists(), "{case}");
        assert_eq!(
            fs::read(out.join("index.saved")).unwrap(),
            fs::read(&index).unwrap(),
            "{case}"
        );
    }
}

#[test]
fn bus_unmapped_writes_ratio_per_processed_read() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("unmapped_out");

    write_fastq(
        &r1,
        &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_unmapped", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(
        &r2,
        &[("seq_read", &transcript), ("seq_unmapped", &[b'N'; 80])],
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
        .expect("run kallistors bus --unmapped");
    assert!(status.success());

    let ratios = fs::read_to_string(out_dir.join("unmapped_ratio.txt")).unwrap();
    assert!(!ratios.trim_end().ends_with(','));
    let values = ratios
        .trim_end()
        .split(',')
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(values.len(), 2);
    assert!(values[0] < values[1]);
    assert_eq!(values[1], 1.0);
}

#[test]
fn bus_zero_pseudoaligned_returns_error_after_writing_outputs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("zero_pseudoaligned_out");

    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &[b'N'; 80])]);

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
        .expect("run kallistors bus with no pseudoaligned reads");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero reads pseudoaligned"));

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 12);
    assert!(records.is_empty());
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).unwrap();
    assert_eq!(json_u64(&run_info, "n_processed"), Some(1));
    assert_eq!(json_u64(&run_info, "n_pseudoaligned"), Some(0));
}

#[test]
fn bus_unpaired_multi_sequence_query_has_only_internal_separator() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, _transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("unpaired_multi_sequence_query_out");

    write_fastq(&r1, &[("r1", b"TTTTNNNN")]);
    write_fastq(&r2, &[("r2", b"NNNN")]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=-1,-1,-1:0,0,4:0,4,8,1,0,4")
        .arg("--long")
        .arg("--threshold")
        .arg("0.5")
        .arg("--unmapped")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus unpaired multi-sequence query");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("zero reads pseudoaligned"), "{stderr}");

    let novel = fs::read_to_string(out_dir.join("novel.fastq")).unwrap();
    let novel_lines = novel.lines().collect::<Vec<_>>();
    assert_eq!(novel_lines.len(), 8);
    assert_eq!(novel_lines[0], "@unmapped");
    assert_eq!(novel_lines[1], "NNNNNNNNN");
    assert_eq!(novel_lines[2], "+");
    assert_eq!(novel_lines[3], "IIII!IIII");
    assert_eq!(novel_lines[4], "@novel_disjointIntersect");
    assert_eq!(novel_lines[5], "NNNNNNNNN");
    assert_eq!(novel_lines[6], "+");
    assert_eq!(novel_lines[7], "IIII!IIII");
}

#[test]
fn bus_list_does_not_require_index_or_output() {
    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("--list")
        .output()
        .expect("run kallistors bus --list");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    for technology in [
        "10xv1",
        "10xv2",
        "10xv3",
        "10xv4",
        "Bulk",
        "ParseV3",
        "SmartSeq2",
        "BDWTA",
        "CELSeq",
        "CELSeq2",
        "DropSeq",
        "inDropsv1",
        "inDropsv2",
        "inDropsv3",
        "MATQSEQ",
        "PETRISEQ",
        "SCRBSeq",
        "SmartSeq3",
        "SPLiT-seq",
        "STORM-seq",
        "SureCell",
        "VASA-seq",
        "Visium",
    ] {
        assert!(
            stdout.contains(technology),
            "technology list should contain {technology}: {stdout}"
        );
    }
}

#[test]
fn bus_aa_maps_coding_reads_against_aa_index() {
    let dir = tempfile::tempdir().expect("tempdir");
    let index = build_aa_index_with_transcript(&dir, b"FFFFFFFFFFFF");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("aa_bus_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(
        &r2,
        &[("seq_read", b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")],
    );

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--aa")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --aa");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].ec, 0);
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).expect("read run_info");
    assert_eq!(json_u64(&run_info, "n_frame_clashes"), Some(2));
}

#[test]
fn bus_aa_ignores_paired_flag_for_single_cdna_technologies() {
    let dir = tempfile::tempdir().expect("tempdir");
    let index = build_aa_index_with_transcript(&dir, b"FFFFFFFFFFFF");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("aa_paired_flag_ignored_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(
        &r2,
        &[("seq_read", b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")],
    );

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--aa")
        .arg("--paired")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --aa --paired");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("[bus] --paired ignored; --aa only supports single-end reads"),
        "{stderr}"
    );

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].ec, 0);
}

#[test]
fn bus_aa_rejects_paired_technologies() {
    let dir = tempfile::tempdir().expect("tempdir");
    let index = build_aa_index_with_transcript(&dir, b"FFFFFFFFFFFF");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let r3 = dir.path().join("r3.fastq");
    let r4 = dir.path().join("r4.fastq");
    let out_dir = dir.path().join("aa_paired_bus_out");

    write_fastq(&r1, &[("bc1", b"ACGTACGT")]);
    write_fastq(&r2, &[("bc2", b"TGCATGCA")]);
    write_fastq(&r3, &[("read1", b"ATTGCGCAATGTTTTTTTT")]);
    write_fastq(&r4, &[("read2", b"AAAAAAAAAAAAAAAA")]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("SMARTSEQ3")
        .arg("--aa")
        .arg(&r1)
        .arg(&r2)
        .arg(&r3)
        .arg(&r4)
        .output()
        .expect("run kallistors bus --aa SMARTSEQ3");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--aa BUS mode currently supports single cDNA read technologies"));
}

// Reverse-complement targets share every canonical k-mer. Counts alone cannot
// detect an ignored strand flag; check each barcode/UMI's transcript membership.
#[test]
fn bus_strand_flags_filter_fastq_paired_and_bam_transcripts() {
    use std::collections::BTreeMap;

    let dir = tempfile::tempdir().unwrap();
    let transcript = b"CCCTCCGCTATGTACTCCACGTCCTGTGGGAACGACCTCCTGTTGCATATGTAACACGATGATAACCCCACTGGACGTGGCACTAATTTACAAGCGCAAAGGTCCATGTAATTCCCGATACGGGGAAGTCAACTCCCCTAAGAATGTATCCAAGTCAGGCTTTGCCTCTATAAGCTATCAATCCAACCGGTCTGATCTGACACCGGTAGGTGCGGTTAAATGGTGAGCGTTTACTTACTTGTCCGGCTGCAGGTTATGCCGGCTTAGAGCGACTGGCACAGGGAAACCTAGATCCCCCAC";
    let reverse = helpers::revcomp(transcript);
    let query = &transcript[60..160];
    let reverse_query = helpers::revcomp(query);
    let fasta = dir.path().join("strand.fa");
    let index = dir.path().join("strand.idx");
    fs::write(
        &fasta,
        format!(
            ">forward\n{}\n>reverse\n{}\n",
            std::str::from_utf8(transcript).unwrap(),
            std::str::from_utf8(&reverse).unwrap()
        ),
    )
    .unwrap();
    build_index(&index, &[fasta], IndexBuildOptions::default()).unwrap();
    let tags = dir.path().join("tags.fq");
    let reads = dir.path().join("reads.fq");
    let left = dir.path().join("left.fq");
    let right = dir.path().join("right.fq");
    let bam = dir.path().join("reads.bam");
    let barcode = b"ACGTACGTACGTACGT";
    let umis = [b"AAAAAAAAAAAA", b"AAAAAAAAAAAC"];
    write_fastq(
        &tags,
        &[
            ("forward", b"ACGTACGTACGTACGTAAAAAAAAAAAA"),
            ("reverse", b"ACGTACGTACGTACGTAAAAAAAAAAAC"),
        ],
    );
    write_fastq(&reads, &[("forward", query), ("reverse", &reverse_query)]);
    let first = &transcript[50..125];
    let second = helpers::revcomp(&transcript[155..230]);
    write_fastq(&left, &[("forward", first), ("reverse", &second)]);
    write_fastq(&right, &[("forward", &second), ("reverse", first)]);
    write_bam(
        &bam,
        &[
            (query, barcode, umis[0]),
            (&reverse_query, barcode, umis[1]),
        ],
    );

    for threads in [1, 8] {
        for loader in [
            None,
            Some("--kallisto-direct-kmer"),
            Some("--kallisto-fallback"),
        ] {
            for (strand, targets) in [
                ("--unstranded", [vec![0, 1], vec![0, 1]]),
                ("--fr-stranded", [vec![0], vec![1]]),
                ("--rf-stranded", [vec![1], vec![0]]),
            ] {
                for layout in ["fastq", "paired", "bam"] {
                    let case = format!("{layout}-{strand}-{loader:?}-{threads}");
                    let out = dir.path().join(&case);
                    let mut command = bus_command(&index, &out);
                    command.args(["-t", &threads.to_string(), strand]);
                    if let Some(loader) = loader {
                        command.arg(loader);
                    }
                    match layout {
                        "fastq" => {
                            command.args(["-x", "10xv3"]).arg(&tags).arg(&reads);
                        }
                        "paired" => {
                            command
                                .args(["-x", "0,0,16:0,16,28:1,0,0,2,0,0", "--paired"])
                                .arg(&tags)
                                .arg(&left)
                                .arg(&right);
                        }
                        "bam" => {
                            command.args(["-x", "10xv3", "--bam"]).arg(&bam);
                        }
                        _ => unreachable!(),
                    }
                    assert_bus_success(&command.output().unwrap(), &case);
                    let ecs: BTreeMap<i32, Vec<u32>> = fs::read_to_string(out.join("matrix.ec"))
                        .unwrap()
                        .lines()
                        .map(|line| {
                            let (id, members) = line.split_once('\t').unwrap();
                            (
                                id.parse().unwrap(),
                                members.split(',').map(|v| v.parse().unwrap()).collect(),
                            )
                        })
                        .collect();
                    let (_, _, records) = read_bus_records(&out.join("output.bus"));
                    assert_eq!(records.len(), 2, "{case}");
                    let mut seen = [false; 2];
                    for record in records {
                        let read = umis
                            .iter()
                            .position(|umi| encode_bus_seq(*umi) == record.umi)
                            .unwrap();
                        assert!(!seen[read], "{case}: duplicated input identity");
                        seen[read] = true;
                        assert_eq!(record.barcode, encode_bus_seq(barcode), "{case}");
                        assert_eq!((record.count, record.flags), (1, 0), "{case}");
                        assert_eq!(ecs[&record.ec], targets[read], "{case}/read{read}");
                    }
                }
            }
        }
    }
}
