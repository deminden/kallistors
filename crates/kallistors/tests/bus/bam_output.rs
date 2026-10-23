//! Bam output regressions.

use super::*;
#[test]
fn bus_bulk_rejects_pseudobam_bam_and_tag_modes() {
    let dir = tempfile::tempdir().expect("tempdir");
    let read = dir.path().join("bulk.fastq");
    let bam = dir.path().join("bulk.bam");
    let (index, transcript) = build_tiny_index(&dir);

    write_fastq(&read, &[("bulk_read", &transcript)]);
    write_bam(&bam, &[(&transcript, b"ACGTACGTACGTACGT", b"TTTTTTTTTT")]);

    for (flag, input, expected) in [
        (
            "--pseudobam",
            read.as_path(),
            "Pseudobam not supported yet in this mode",
        ),
        ("--bam", bam.as_path(), "--bam not supported in this mode"),
        ("--tag", read.as_path(), "--tag not supported in this mode"),
    ] {
        let out_dir = dir.path().join(format!(
            "bulk_reject_{}",
            flag.trim_start_matches('-').replace('-', "_")
        ));
        let mut command = Command::new(env!("CARGO_BIN_EXE_kallistors"));
        command
            .arg("bus")
            .arg("-i")
            .arg(&index)
            .arg("-o")
            .arg(&out_dir)
            .arg("-x")
            .arg("Bulk");
        if flag == "--tag" {
            command.arg("--tag").arg("AAAA");
        } else {
            command.arg(flag);
        }
        let output = command
            .arg(input)
            .output()
            .expect("run kallistors bus Bulk incompatible mode");
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(expected), "{stderr}");
    }
}

#[test]
fn bus_smartseq3_missing_tag_keeps_prefix_in_pseudobam_read() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bc1 = dir.path().join("bc1.fastq");
    let bc2 = dir.path().join("bc2.fastq");
    let r1 = dir.path().join("smartseq3_no_tag_pseudobam_r1.fastq");
    let r2 = dir.path().join("smartseq3_no_tag_pseudobam_r2.fastq");
    let out_dir = dir.path().join("smartseq3_no_tag_pseudobam_out");

    write_fastq(&bc1, &[("bc1", b"ACGT")]);
    write_fastq(&bc2, &[("bc2", b"TGCA")]);
    let prefix = b"GGGGGGGGGGGGGGGGGGG";
    let mut r1_sequence = prefix.to_vec();
    r1_sequence.extend_from_slice(&transcript);
    let r1_quality = (0..r1_sequence.len())
        .map(|idx| b'!' + u8::try_from(idx % 40).unwrap())
        .collect::<Vec<_>>();
    let r2_sequence = reverse_complement(&transcript);
    write_fastq_with_qualities(&r1, &[("r1", &r1_sequence, &r1_quality)]);
    write_fastq(&r2, &[("r2", &r2_sequence)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("SMARTSEQ3")
        .arg("--pseudobam")
        .arg(&bc1)
        .arg(&bc2)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus smartseq3 missing tag --pseudobam");
    assert!(status.success());

    let (_bc_len, umi_len, bus_records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(umi_len, 8);
    assert_eq!(bus_records.len(), 1);
    assert_eq!(bus_records[0].umi, u64::MAX);

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open pseudobam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read pseudobam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read pseudobam records");
    assert_eq!(records.len(), 2);
    assert_eq!(
        records[0].sequence().iter().collect::<Vec<_>>(),
        r1_sequence
    );
    assert_eq!(
        records[0].quality_scores().iter().collect::<Vec<_>>(),
        r1_quality
    );
    assert_eq!(
        records[1].sequence().iter().collect::<Vec<_>>(),
        r2_sequence
    );
}

#[test]
fn bus_pseudobam_rejects_unpaired_multi_sequence_technology() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("unpaired_multi_pseudobam_out");

    write_fastq(&r1, &[("cell", b"ACGTACGT")]);
    write_fastq(&r2, &[("seq", &transcript)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=-1,-1,-1:0,0,4:1,0,31,1,31,62")
        .arg("--pseudobam")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus unpaired multi-sequence --pseudobam");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("BAM output is currently only supported")
    );
}

#[test]
fn bus_rejects_pseudobam_and_genomebam_together() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("conflicting_bam_output_modes_out");

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
        .arg("--pseudobam")
        .arg("--genomebam")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --pseudobam --genomebam");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--pseudobam and --genomebam are mutually exclusive"),
        "{stderr}"
    );
    assert!(!stderr.contains("--genomebam requires --gtf"), "{stderr}");
}

#[test]
fn bus_pseudobam_writes_transcriptome_bam() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("pseudobam_out");
    let qualities = (0..transcript.len())
        .map(|idx| b'!' + u8::try_from(idx % 40).unwrap())
        .collect::<Vec<_>>();

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq_with_qualities(&r2, &[("seq_read", &transcript, &qualities)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--pseudobam")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --pseudobam");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open pseudobam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let header = reader.read_header().expect("read pseudobam header");
    assert_eq!(header.reference_sequences().len(), 1);
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read pseudobam records");
    assert_eq!(records.len(), 1);
    assert_eq!(
        records[0].name().map(|name| name.as_ref()),
        Some(&b"seq_read"[..])
    );
    assert_eq!(records[0].sequence().iter().collect::<Vec<_>>(), transcript);
    assert_eq!(
        records[0].quality_scores().iter().collect::<Vec<_>>(),
        qualities
    );
    assert_eq!(
        records[0].reference_sequence_id().transpose().unwrap(),
        Some(0)
    );
    assert!(records[0].alignment_start().transpose().unwrap().is_some());
    assert!(!records[0].flags().is_unmapped());
}

#[test]
fn bus_pseudobam_writes_one_record_per_paired_mate() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("mate1.fastq");
    let r2 = dir.path().join("mate2.fastq");
    let out_dir = dir.path().join("paired_pseudobam_out");
    let mate1 = &transcript[..40];
    let mate2 = &transcript[20..60];
    let qual1 = vec![b'#'; mate1.len()];
    let qual2 = vec![b'5'; mate2.len()];

    write_fastq_with_qualities(&r1, &[("read", mate1, &qual1)]);
    write_fastq_with_qualities(&r2, &[("read", mate2, &qual2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=-1,0,0:-1,0,0:0,0,0,1,0,0")
        .arg("--paired")
        .arg("--pseudobam")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors paired bus --pseudobam");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open pseudobam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read pseudobam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read pseudobam records");
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence().iter().collect::<Vec<_>>(), mate1);
    assert_eq!(records[1].sequence().iter().collect::<Vec<_>>(), mate2);
    assert_eq!(
        records[0].quality_scores().iter().collect::<Vec<_>>(),
        qual1
    );
    assert_eq!(
        records[1].quality_scores().iter().collect::<Vec<_>>(),
        qual2
    );
    assert!(records[0].flags().is_segmented());
    assert!(records[0].flags().is_first_segment());
    assert!(records[1].flags().is_segmented());
    assert!(records[1].flags().is_last_segment());
    assert_eq!(
        records[0].mate_reference_sequence_id().transpose().unwrap(),
        records[1].reference_sequence_id().transpose().unwrap()
    );
    assert_eq!(
        records[1].mate_reference_sequence_id().transpose().unwrap(),
        records[0].reference_sequence_id().transpose().unwrap()
    );
    assert_eq!(
        records[0].mate_alignment_start().transpose().unwrap(),
        records[1].alignment_start().transpose().unwrap()
    );
    assert_eq!(
        records[1].mate_alignment_start().transpose().unwrap(),
        records[0].alignment_start().transpose().unwrap()
    );
    assert!(records[0].template_length() > 0);
    assert_eq!(records[0].template_length(), -records[1].template_length());
}

#[test]
fn bus_pseudobam_marks_unmapped_mate_without_proper_pair_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("mate1.fastq");
    let r2 = dir.path().join("mate2.fastq");
    let out_dir = dir.path().join("paired_pseudobam_one_mate_out");
    let mate1 = &transcript[..40];
    let mate2 = vec![b'N'; 40];

    write_fastq(&r1, &[("read", mate1)]);
    write_fastq(&r2, &[("read", &mate2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x=-1,0,0:-1,0,0:0,0,0,1,0,0")
        .arg("--paired")
        .arg("--pseudobam")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors paired bus --pseudobam with one unmapped mate");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open pseudobam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read pseudobam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read pseudobam records");
    assert_eq!(records.len(), 2);
    assert!(!records[0].flags().is_unmapped());
    assert!(records[0].flags().is_mate_unmapped());
    assert!(!records[0].flags().is_properly_segmented());
    assert!(records[1].flags().is_unmapped());
    assert!(!records[1].flags().is_mate_unmapped());
    assert!(!records[1].flags().is_properly_segmented());
    assert!(records[0].data().get(&Tag::ALIGNMENT_HIT_COUNT).is_some());
    assert!(records[1].data().get(&Tag::ALIGNMENT_HIT_COUNT).is_none());
    assert_eq!(
        records[1].mate_reference_sequence_id().transpose().unwrap(),
        records[0].reference_sequence_id().transpose().unwrap()
    );
    assert_eq!(
        records[1].mate_alignment_start().transpose().unwrap(),
        records[0].alignment_start().transpose().unwrap()
    );
}

#[test]
fn bus_genomebam_writes_projected_bam() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT";
    let (index, transcript) = build_index_with_named_transcript(&dir, "tx0.1", transcript);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let chromosomes = dir.path().join("chromosomes.txt");
    let out_dir = dir.path().join("genomebam_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);
    fs::write(
        &gtf,
        concat!(
            "chr1\ttest\tgene\t101\t176\t.\t+\t.\tgene_id \"gene0\"; gene_version \"2\"; gene_name \"Gene Zero\";\n",
            "chr1\ttest\texon\t101\t176\t.\t+\t.\tgene_id \"gene0\"; gene_version \"2\"; transcript_id \"tx0\"; transcript_version \"1\";\n",
        ),
    )
    .expect("write GTF");
    fs::write(&chromosomes, "chr1\t1000\n").expect("write chromosomes");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg("--chromosomes")
        .arg(&chromosomes)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --genomebam");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open genomebam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let header = reader.read_header().expect("read genomebam header");
    assert_eq!(header.reference_sequences().len(), 1);
    assert!(header.reference_sequences().contains_key(&b"chr1"[..]));
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read genomebam records");
    assert_eq!(records.len(), 1);
    assert_eq!(
        records[0].reference_sequence_id().transpose().unwrap(),
        Some(0)
    );
    assert!(records[0].alignment_start().transpose().unwrap().is_some());
    assert_eq!(records[0].sequence().iter().collect::<Vec<_>>(), transcript);
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.genelist.txt")).expect("read gene list"),
        "0\tgene0.2\tGene Zero\n"
    );
    let bai_path = out_dir.join("pseudoalignments.bam.bai");
    assert!(bai_path.exists());
    let bai = noodles_bam::bai::fs::read(bai_path).expect("read genomebam index");
    assert_eq!(bai.reference_sequences().len(), 1);
}

#[test]
fn bus_genomebam_writes_gene_list_from_exons_without_gene_feature() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT";
    let (index, transcript) = build_index_with_named_transcript(&dir, "tx0.1", transcript);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("exons_only.gtf");
    let out_dir = dir.path().join("genomebam_exons_only_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);
    fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t176\t.\t+\t.\tgene_id \"gene0\"; gene_version \"2\"; gene_name \"Gene Zero\"; transcript_id \"tx0\"; transcript_version \"1\";\n",
    )
    .expect("write exon-only GTF");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --genomebam with exon-only GTF");
    assert!(status.success());

    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.genelist.txt")).expect("read gene list"),
        "0\tgene0.2\tGene Zero\n"
    );
}

#[test]
fn bus_genomebam_writes_unmapped_record_without_gtf_projection() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT";
    let (index, transcript) = build_index_with_named_transcript(&dir, "tx0", transcript);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let out_dir = dir.path().join("genomebam_unprojected_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);
    fs::write(
        &gtf,
        "chr1\ttest\tgene\t101\t176\t.\t+\t.\tgene_id \"gene0\"; gene_name \"Gene Zero\";\n",
    )
    .expect("write GTF without transcript exon");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --genomebam without transcript projection");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open genomebam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read genomebam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read genomebam records");
    assert_eq!(records.len(), 1);
    assert!(records[0].flags().is_unmapped());
    assert!(
        records[0]
            .reference_sequence_id()
            .transpose()
            .unwrap()
            .is_none()
    );
    assert!(out_dir.join("pseudoalignments.bam.bai").exists());
}

#[test]
fn bus_genomebam_sorts_records_by_coordinate() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTTGCAACGAGTCGATCGTACGATCGATCGTTGACCGTAACTGGCATACGCTAGGTCATGCA";
    let (index, transcript) = build_index_with_transcript(&dir, transcript);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let out_dir = dir.path().join("genomebam_sorted_out");
    let late = &transcript[32..64];
    let early = &transcript[0..32];

    write_fastq(
        &r1,
        &[
            ("late_cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("early_cell", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(&r2, &[("late_read", late), ("early_read", early)]);
    fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t164\t.\t+\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n",
    )
    .expect("write GTF");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus sorted --genomebam");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open genomebam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read genomebam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read genomebam records");
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence().iter().collect::<Vec<_>>(), early);
    assert_eq!(records[1].sequence().iter().collect::<Vec<_>>(), late);
    let first_start = records[0]
        .alignment_start()
        .transpose()
        .unwrap()
        .expect("first start");
    let second_start = records[1]
        .alignment_start()
        .transpose()
        .unwrap()
        .expect("second start");
    assert!(first_start < second_start);
}

#[test]
fn bus_genomebam_projects_spliced_cigar() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCATGCATGCATGCA";
    let (index, transcript) = build_index_with_transcript(&dir, transcript);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let out_dir = dir.path().join("genomebam_spliced_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);
    fs::write(
        &gtf,
        concat!(
            "chr1\ttest\texon\t101\t132\t.\t+\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n",
            "chr1\ttest\texon\t143\t174\t.\t+\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n"
        ),
    )
    .expect("write GTF");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus spliced --genomebam");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open genomebam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read genomebam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read genomebam records");
    assert_eq!(records.len(), 1);
    let ops = records[0]
        .cigar()
        .iter()
        .collect::<Result<Vec<_>, _>>()
        .expect("read CIGAR");
    assert_eq!(ops.len(), 3);
    assert_eq!((ops[0].kind(), ops[0].len()), (Kind::Match, 32));
    assert_eq!((ops[1].kind(), ops[1].len()), (Kind::Skip, 10));
    assert_eq!((ops[2].kind(), ops[2].len()), (Kind::Match, 32));
}

#[test]
fn bus_genomebam_projects_negative_strand_spliced_cigar() {
    let dir = tempfile::tempdir().expect("tempdir");
    let transcript = b"ACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCATGCATGCATGCA";
    let (index, transcript) = build_index_with_transcript(&dir, transcript);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let out_dir = dir.path().join("genomebam_negative_spliced_out");
    let qualities = (0..transcript.len())
        .map(|idx| b'!' + u8::try_from(idx % 40).unwrap())
        .collect::<Vec<_>>();

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq_with_qualities(&r2, &[("seq_read", &transcript, &qualities)]);
    fs::write(
        &gtf,
        concat!(
            "chr1\ttest\texon\t101\t132\t.\t-\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n",
            "chr1\ttest\texon\t143\t174\t.\t-\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n"
        ),
    )
    .expect("write GTF");

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus negative-strand --genomebam");
    assert!(status.success());

    let bam = fs::File::open(out_dir.join("pseudoalignments.bam")).expect("open genomebam");
    let mut reader = noodles_bam::io::Reader::new(bam);
    let _header = reader.read_header().expect("read genomebam header");
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("read genomebam records");
    assert_eq!(records.len(), 1);
    assert!(records[0].flags().is_reverse_complemented());
    assert_eq!(
        records[0].sequence().iter().collect::<Vec<_>>(),
        reverse_complement(&transcript)
    );
    assert_eq!(
        records[0].quality_scores().iter().collect::<Vec<_>>(),
        qualities.iter().rev().copied().collect::<Vec<_>>()
    );
    let ops = records[0]
        .cigar()
        .iter()
        .collect::<Result<Vec<_>, _>>()
        .expect("read CIGAR");
    assert_eq!(ops.len(), 3);
    assert_eq!((ops[0].kind(), ops[0].len()), (Kind::Match, 32));
    assert_eq!((ops[1].kind(), ops[1].len()), (Kind::Skip, 10));
    assert_eq!((ops[2].kind(), ops[2].len()), (Kind::Match, 32));
}

#[test]
fn bus_genomebam_requires_gtf() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("genomebam_missing_gtf_out");

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
        .arg("--genomebam")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --genomebam without gtf");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--genomebam requires --gtf"));
    assert!(!stderr.contains("unexpected argument"));
}

#[test]
fn bus_genomebam_rejects_missing_gtf_file() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("missing.gtf");
    let out_dir = dir.path().join("genomebam_missing_gtf_file_out");

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
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --genomebam with missing gtf file");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("GTF file not found"));
}

#[test]
fn bus_genomebam_rejects_invalid_exons_and_inconsistent_transcripts() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &transcript)]);
    for (case, exons, error) in [
        (
            "zero_start",
            &[("chr1", 0, 75, "+")][..],
            "invalid GTF exon start: 0",
        ),
        (
            "reversed_interval",
            &[("chr1", 75, 74, "+")][..],
            "invalid GTF exon interval: 75-74",
        ),
        (
            "invalid_strand",
            &[("chr1", 1, 75, ".")][..],
            "invalid GTF exon strand: .",
        ),
        (
            "split_chromosome",
            &[("chr1", 1, 30, "+"), ("chr2", 31, 75, "+")][..],
            "GTF transcript tx0 has exons on multiple chromosomes",
        ),
        (
            "split_strand",
            &[("chr1", 1, 30, "+"), ("chr1", 31, 75, "-")][..],
            "GTF transcript tx0 has exons on multiple strands",
        ),
    ] {
        let contents: String = exons.iter().map(|(chr, start, end, strand)|
            format!("{chr}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n")).collect();
        fs::write(&gtf, contents).unwrap();
        let output = bus_command(&index, &dir.path().join(case))
            .args(["-x", "10XV3", "--genomebam", "--gtf"])
            .arg(&gtf)
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert_bus_error(&output, case, error);
    }
}

#[test]
fn bus_genomebam_rejects_missing_chromosome_file() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let chromosomes = dir.path().join("missing.chromosomes.txt");
    let out_dir = dir.path().join("genomebam_missing_chromosomes_out");

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read", &transcript)]);
    fs::write(
        &gtf,
        "chr1\ttest\texon\t1\t75\t.\t+\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n",
    )
    .expect("write GTF");

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg("--genomebam")
        .arg("--gtf")
        .arg(&gtf)
        .arg("--chromosomes")
        .arg(&chromosomes)
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("run kallistors bus --genomebam with missing chromosome file");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Chromosome file not found"));
}

#[test]
fn bus_genomebam_rejects_invalid_chromosome_tables() {
    let dir = tempfile::tempdir().unwrap();
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let gtf = dir.path().join("genes.gtf");
    let chromosomes = dir.path().join("chromosomes.txt");
    write_fastq(&r1, &[("cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq", &transcript)]);
    fs::write(
        &gtf,
        "chr1\ttest\texon\t1\t75\t.\t+\t.\tgene_id \"gene0\"; transcript_id \"tx0\";\n",
    )
    .unwrap();
    for (case, contents, error) in [
        (
            "duplicate",
            "chr1\t1000\nchr1\t2000\n",
            "duplicate chromosome name in line",
        ),
        (
            "zero_length",
            "chr1\t0\n",
            "invalid chromosome length in line",
        ),
    ] {
        fs::write(&chromosomes, contents).unwrap();
        let output = bus_command(&index, &dir.path().join(case))
            .args(["-x", "10XV3", "--genomebam", "--gtf"])
            .arg(&gtf)
            .arg("--chromosomes")
            .arg(&chromosomes)
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert_bus_error(&output, case, error);
    }
}
