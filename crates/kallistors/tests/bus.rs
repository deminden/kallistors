use std::fs;
use std::io::Write as _;
use std::process::{Command, Stdio};

use kallistors::index::{IndexBuildOptions, build_index};
use noodles_sam::alignment::{
    RecordBuf,
    io::Write as _,
    record::{Flags, cigar::op::Kind, data::field::Tag},
    record_buf::{Sequence, data::field::Value},
};

type BamRecordFixture<'a> = (&'a [u8], &'a [u8], &'a [u8], Flags);
type OptionalBamTagsFixture<'a> = (&'a [u8], Option<&'a [u8]>, Option<&'a [u8]>);

fn write_fastq(path: &std::path::Path, records: &[(&str, &[u8])]) {
    let mut out = Vec::new();
    for (name, seq) in records {
        out.extend_from_slice(b"@");
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(b"\n");
        out.extend_from_slice(seq);
        out.extend_from_slice(b"\n+\n");
        out.extend(std::iter::repeat_n(b'I', seq.len()));
        out.extend_from_slice(b"\n");
    }
    fs::write(path, out).expect("write FASTQ");
}

fn write_fastq_gz(path: &std::path::Path, records: &[(&str, &[u8])]) {
    let file = fs::File::create(path).expect("create gz FASTQ");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for (name, seq) in records {
        encoder.write_all(b"@").expect("write gz FASTQ");
        encoder.write_all(name.as_bytes()).expect("write gz FASTQ");
        encoder.write_all(b"\n").expect("write gz FASTQ");
        encoder.write_all(seq).expect("write gz FASTQ");
        encoder.write_all(b"\n+\n").expect("write gz FASTQ");
        for _ in 0..seq.len() {
            encoder.write_all(b"I").expect("write gz FASTQ");
        }
        encoder.write_all(b"\n").expect("write gz FASTQ");
    }
    encoder.finish().expect("finish gz FASTQ");
}

fn write_fastq_with_qualities(path: &std::path::Path, records: &[(&str, &[u8], &[u8])]) {
    let mut out = Vec::new();
    for (name, seq, qual) in records {
        assert_eq!(seq.len(), qual.len());
        out.extend_from_slice(b"@");
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(b"\n");
        out.extend_from_slice(seq);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(qual);
        out.extend_from_slice(b"\n");
    }
    fs::write(path, out).expect("write FASTQ");
}

fn write_bam(path: &std::path::Path, records: &[(&[u8], &[u8], &[u8])]) {
    let records = records
        .iter()
        .map(|(seq, barcode, umi)| (*seq, *barcode, *umi, Flags::empty()))
        .collect::<Vec<_>>();
    write_bam_with_flags(path, &records);
}

fn write_bam_with_corrected_barcode(
    path: &std::path::Path,
    seq: &[u8],
    barcode: &[u8],
    umi: &[u8],
) {
    let mut writer = noodles_bam::io::Writer::new(Vec::new());
    let header = noodles_sam::Header::default();
    writer.write_header(&header).expect("write BAM header");
    let barcode = std::str::from_utf8(barcode).expect("ASCII barcode");
    let umi = std::str::from_utf8(umi).expect("ASCII UMI");
    let record = RecordBuf::builder()
        .set_sequence(Sequence::from(seq))
        .set_data(
            [
                (Tag::CELL_BARCODE_ID, Value::from(barcode)),
                (Tag::UMI_SEQUENCE, Value::from(umi)),
            ]
            .into_iter()
            .collect(),
        )
        .build();
    writer
        .write_alignment_record(&header, &record)
        .expect("write BAM record");
    writer.try_finish().expect("finish BAM");
    fs::write(path, writer.get_ref().get_ref()).expect("write BAM");
}

fn write_bam_with_raw_and_corrected_barcode(
    path: &std::path::Path,
    seq: &[u8],
    raw_barcode: &[u8],
    corrected_barcode: &[u8],
    umi: &[u8],
) {
    let mut writer = noodles_bam::io::Writer::new(Vec::new());
    let header = noodles_sam::Header::default();
    writer.write_header(&header).expect("write BAM header");
    let raw_barcode = std::str::from_utf8(raw_barcode).expect("ASCII barcode");
    let corrected_barcode = std::str::from_utf8(corrected_barcode).expect("ASCII barcode");
    let umi = std::str::from_utf8(umi).expect("ASCII UMI");
    let record = RecordBuf::builder()
        .set_sequence(Sequence::from(seq))
        .set_data(
            [
                (Tag::CELL_BARCODE_SEQUENCE, Value::from(raw_barcode)),
                (Tag::CELL_BARCODE_ID, Value::from(corrected_barcode)),
                (Tag::UMI_SEQUENCE, Value::from(umi)),
            ]
            .into_iter()
            .collect(),
        )
        .build();
    writer
        .write_alignment_record(&header, &record)
        .expect("write BAM record");
    writer.try_finish().expect("finish BAM");
    fs::write(path, writer.get_ref().get_ref()).expect("write BAM");
}

fn write_bam_with_corrected_barcode_and_umi(
    path: &std::path::Path,
    seq: &[u8],
    barcode: &[u8],
    umi: &[u8],
    umi_tag: Tag,
) {
    let mut writer = noodles_bam::io::Writer::new(Vec::new());
    let header = noodles_sam::Header::default();
    writer.write_header(&header).expect("write BAM header");
    let barcode = std::str::from_utf8(barcode).expect("ASCII barcode");
    let umi = std::str::from_utf8(umi).expect("ASCII UMI");
    let record = RecordBuf::builder()
        .set_sequence(Sequence::from(seq))
        .set_data(
            [
                (Tag::CELL_BARCODE_ID, Value::from(barcode)),
                (umi_tag, Value::from(umi)),
            ]
            .into_iter()
            .collect(),
        )
        .build();
    writer
        .write_alignment_record(&header, &record)
        .expect("write BAM record");
    writer.try_finish().expect("finish BAM");
    fs::write(path, writer.get_ref().get_ref()).expect("write BAM");
}

fn write_bam_with_flags(path: &std::path::Path, records: &[BamRecordFixture<'_>]) {
    let mut writer = noodles_bam::io::Writer::new(Vec::new());
    let header = noodles_sam::Header::default();
    writer.write_header(&header).expect("write BAM header");
    for (seq, barcode, umi, flags) in records {
        let barcode = std::str::from_utf8(barcode).expect("ASCII barcode");
        let umi = std::str::from_utf8(umi).expect("ASCII UMI");
        let record = RecordBuf::builder()
            .set_flags(*flags)
            .set_sequence(Sequence::from(*seq))
            .set_data(
                [
                    (Tag::CELL_BARCODE_SEQUENCE, Value::from(barcode)),
                    (Tag::UMI_SEQUENCE, Value::from(umi)),
                ]
                .into_iter()
                .collect(),
            )
            .build();
        writer
            .write_alignment_record(&header, &record)
            .expect("write BAM record");
    }
    writer.try_finish().expect("finish BAM");
    fs::write(path, writer.get_ref().get_ref()).expect("write BAM");
}

fn write_bam_with_optional_tags(path: &std::path::Path, records: &[OptionalBamTagsFixture<'_>]) {
    let mut writer = noodles_bam::io::Writer::new(Vec::new());
    let header = noodles_sam::Header::default();
    writer.write_header(&header).expect("write BAM header");
    for (seq, barcode, umi) in records {
        let mut data = Vec::new();
        if let Some(barcode) = barcode {
            let barcode = std::str::from_utf8(barcode).expect("ASCII barcode");
            data.push((Tag::CELL_BARCODE_SEQUENCE, Value::from(barcode)));
        }
        if let Some(umi) = umi {
            let umi = std::str::from_utf8(umi).expect("ASCII UMI");
            data.push((Tag::UMI_SEQUENCE, Value::from(umi)));
        }
        let record = RecordBuf::builder()
            .set_sequence(Sequence::from(*seq))
            .set_data(data.into_iter().collect())
            .build();
        writer
            .write_alignment_record(&header, &record)
            .expect("write BAM record");
    }
    writer.try_finish().expect("finish BAM");
    fs::write(path, writer.get_ref().get_ref()).expect("write BAM");
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

fn read_u32_le(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}

fn read_i32_le(bytes: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}

fn read_u64_le(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .status()
        .map(|status| status.success())
        .unwrap_or(false)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BusRecord {
    barcode: u64,
    umi: u64,
    ec: i32,
    count: u32,
    flags: u32,
}

fn read_bus_records(path: &std::path::Path) -> (u32, u32, Vec<BusRecord>) {
    let bytes = fs::read(path).expect("read bus");
    assert_eq!(&bytes[0..4], b"BUS\0");
    let bc_len = read_u32_le(&bytes, 8);
    let umi_len = read_u32_le(&bytes, 12);
    let text_len = read_u32_le(&bytes, 16) as usize;
    let mut offset = 20 + text_len;
    let mut records = Vec::new();
    while offset < bytes.len() {
        records.push(BusRecord {
            barcode: read_u64_le(&bytes, offset),
            umi: read_u64_le(&bytes, offset + 8),
            ec: read_i32_le(&bytes, offset + 16),
            count: read_u32_le(&bytes, offset + 20),
            flags: read_u32_le(&bytes, offset + 24),
        });
        offset += 32;
    }
    (bc_len, umi_len, records)
}

fn encode_bus_seq(seq: &[u8]) -> u64 {
    seq.iter().fold(0u64, |acc, base| {
        (acc << 2)
            | match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            }
    })
}

fn build_tiny_index(dir: &tempfile::TempDir) -> (std::path::PathBuf, Vec<u8>) {
    build_index_with_transcript(
        dir,
        b"ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT",
    )
}

fn build_index_with_transcript(
    dir: &tempfile::TempDir,
    transcript: &[u8],
) -> (std::path::PathBuf, Vec<u8>) {
    build_index_with_named_transcript(dir, "tx0", transcript)
}

fn build_index_with_named_transcript(
    dir: &tempfile::TempDir,
    name: &str,
    transcript: &[u8],
) -> (std::path::PathBuf, Vec<u8>) {
    build_index_with_named_transcript_and_k(dir, name, transcript, 31)
}

fn build_index_with_named_transcript_and_k(
    dir: &tempfile::TempDir,
    name: &str,
    transcript: &[u8],
    k: usize,
) -> (std::path::PathBuf, Vec<u8>) {
    let fasta = dir.path().join("transcripts.fa");
    let index = dir.path().join("transcripts.idx");
    fs::write(
        &fasta,
        format!(">{name}\n{}\n", std::str::from_utf8(transcript).unwrap()),
    )
    .expect("write fasta");
    build_index(
        &index,
        std::slice::from_ref(&fasta),
        IndexBuildOptions {
            k,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build index");
    (index, transcript.to_vec())
}

fn build_aa_index_with_transcript(dir: &tempfile::TempDir, protein: &[u8]) -> std::path::PathBuf {
    let fasta = dir.path().join("proteins.fa");
    let index = dir.path().join("proteins.idx");
    fs::write(
        &fasta,
        format!(">tx0\n{}\n", std::str::from_utf8(protein).unwrap()),
    )
    .expect("write protein FASTA");
    build_index(
        &index,
        std::slice::from_ref(&fasta),
        IndexBuildOptions {
            k: 31,
            aa: true,
            ..IndexBuildOptions::default()
        },
    )
    .expect("build AA index");
    index
}

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
    assert!(
        String::from_utf8_lossy(&output.stderr).contains(
            "[bus] Note: Strand option was not specified; index has no strand annotations, processing as --unstranded"
        )
    );

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
    for expected in [
        "\"n_targets\": 1",
        "\"n_bootstraps\": 0",
        "\"n_processed\": 1",
        "\"n_pseudoaligned\": 1",
        "\"n_unique\": 1",
        "\"p_pseudoaligned\": 100.0",
        "\"p_unique\": 100.0",
        "\"index_version\": 13",
        "\"k-mer length\": 31",
        "\"start_time\":",
        "\"call\":",
    ] {
        assert!(
            run_info.contains(expected),
            "run_info.json should contain {expected}: {run_info}"
        );
    }
    assert!(
        !run_info.contains("\"technology\""),
        "run_info.json should match upstream BUS schema without technology: {run_info}"
    );
}

#[test]
fn bus_bulk_defaults_to_single_end_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let read = dir.path().join("bulk.fastq");
    let out_dir = dir.path().join("bulk_single_end_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(&read, &[("bulk_read", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("Bulk")
        .arg(&read)
        .status()
        .expect("run kallistors bus -x Bulk");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, 0);
    assert_eq!(records[0].umi, 0);
    assert_eq!(records[0].ec, 0);
    assert!(out_dir.join("index.saved").exists());
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "batch0\n"
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.sample.barcodes")).unwrap(),
        "AAAAAAAAAAAAAAAA\n"
    );
}

#[test]
fn bus_without_technology_defaults_to_bulk_single_end() {
    let dir = tempfile::tempdir().expect("tempdir");
    let read = dir.path().join("bulk.fastq");
    let out_dir = dir.path().join("bulk_no_technology_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(&read, &[("bulk_read", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg(&read)
        .status()
        .expect("run kallistors bus without -x");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "batch0\n"
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.sample.barcodes")).unwrap(),
        "AAAAAAAAAAAAAAAA\n"
    );
}

#[test]
fn bus_without_technology_supports_bulk_paired_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let left = dir.path().join("left.fastq");
    let right = dir.path().join("right.fastq");
    let out_dir = dir.path().join("bulk_no_technology_paired_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(&left, &[("left", &transcript[..40])]);
    write_fastq(&right, &[("right", &transcript[20..60])]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("--paired")
        .arg(&left)
        .arg(&right)
        .status()
        .expect("run paired kallistors bus without -x");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert!(out_dir.join("flens.txt").exists());
    assert!(out_dir.join("index.saved").exists());
}

#[test]
fn bus_without_technology_infers_paired_direct_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let left = dir.path().join("left.fastq");
    let right = dir.path().join("right.fastq");
    let out_dir = dir.path().join("bulk_no_technology_inferred_paired_out");

    let (index, transcript) = build_tiny_index(&dir);
    write_fastq(&left, &[("left", &transcript[..40])]);
    write_fastq(&right, &[("right", &transcript[20..60])]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg(&left)
        .arg(&right)
        .status()
        .expect("run inferred paired kallistors bus without -x");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 1);
    assert_eq!(records.len(), 1);
    assert!(out_dir.join("flens.txt").exists());
    assert!(out_dir.join("index.saved").exists());
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
fn bus_accepts_verbose_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("verbose_bus_out");
    let (index, transcript) = build_tiny_index(&dir);

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
        .arg("--verbose")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --verbose");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}

#[test]
fn bus_accepts_dfk_onlist_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("dfk_onlist_bus_out");
    let (index, transcript) = build_tiny_index(&dir);

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
        .arg("--dfk-onlist")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus --dfk-onlist");
    assert!(status.success());

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
}

#[test]
fn bus_reads_gzipped_fastq_inputs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq.gz");
    let r2 = dir.path().join("r2.fastq.gz");
    let out_dir = dir.path().join("bus_gz_out");
    let (index, transcript) = build_tiny_index(&dir);

    write_fastq_gz(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq_gz(&r2, &[("seq_read", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV3")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus with gzipped FASTQ");
    assert!(status.success());

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
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 1")
    );
}

#[test]
fn bus_tenx_v3_matches_upstream_kallisto_on_tiny_case() {
    if !kallisto_available() {
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let kallistors_out = dir.path().join("kallistors_bus");
    let kallisto_out = dir.path().join("kallisto_bus");
    let (index, transcript) = build_tiny_index(&dir);

    write_fastq(&r1, &[("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&r2, &[("seq_read instrument:1", &transcript)]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&kallistors_out)
        .arg("-x")
        .arg("10XV3")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallistors bus");
    assert!(status.success());

    let status = Command::new("kallisto")
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&kallisto_out)
        .arg("-x")
        .arg("10XV3")
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("run kallisto bus");
    if !status.success() {
        return;
    }

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
fn bus_batch_writes_cells_and_batch_barcodes() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let a_r1 = dir.path().join("a_r1.fastq");
    let a_r2 = dir.path().join("a_r2.fastq");
    let b_r1 = dir.path().join("b_r1.fastq");
    let b_r2 = dir.path().join("b_r2.fastq");
    let batch = dir.path().join("batch.txt");
    let out_dir = dir.path().join("batch_out");

    write_fastq(&a_r1, &[("a_cell", b"ACGTACGTACGTACGTTTTTTTTTTTTT")]);
    write_fastq(&a_r2, &[("a_seq", &transcript)]);
    write_fastq(&b_r1, &[("b_cell", b"TGCATGCATGCATGCATTTTTTTTTTTT")]);
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
        .expect("run kallistors bus batch");
    assert!(status.success());

    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\nsample_b\n"
    );
    let sample_barcodes = fs::read_to_string(out_dir.join("matrix.sample.barcodes")).unwrap();
    assert_eq!(sample_barcodes.lines().count(), 2);
    let bus = fs::read(out_dir.join("output.bus")).expect("read bus");
    assert_eq!(read_u32_le(&bus, 8), 32);
    assert_eq!(read_u32_le(&bus, 12), 12);
    let text_len = read_u32_le(&bus, 16) as usize;
    assert_eq!(bus.len(), 20 + text_len + 64);
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 2")
    );
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
    assert!(run_info.contains("\"n_processed\": 1"), "{run_info}");
    assert!(run_info.contains("\"n_pseudoaligned\": 0"), "{run_info}");
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
    let out_dir = dir.path().join("batch_paired_flens_out");

    write_fastq(&a_r1, &[("a/1", &transcript[..40])]);
    write_fastq(&a_r2, &[("a/2", &transcript[20..60])]);
    write_fastq(&b_r1, &[("b/1", &transcript[..35])]);
    write_fastq(&b_r2, &[("b/2", &transcript[10..45])]);
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
    for line in lines {
        let counts = line
            .split_whitespace()
            .map(|value| value.parse::<u32>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(counts.len(), kallistors::pseudoalign::MAX_FRAG_LEN as usize);
    }
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

    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert!(records.is_empty());
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).unwrap();
    assert!(run_info.contains("\"n_processed\": 1"));
    assert!(run_info.contains("\"n_pseudoaligned\": 0"));
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
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 2")
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
    assert!(novel.starts_with("@unmapped\n"));
    assert!(novel.contains("@novel_disjointIntersect\n"));
    assert!(novel.contains(std::str::from_utf8(&[b'N'; 80]).unwrap()));
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

#[test]
fn bus_bam_reads_sequence_and_barcode_umi_tags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_out");

    write_bam(
        &bam,
        &[
            (&transcript, b"ACGTACGTACGTACGT", b"TTTTTTTTTT"),
            (&transcript, b"TGCATGCATGCATGCA", b"AAAAAAAAAA"),
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
        .arg(&bam)
        .status()
        .expect("run kallistors bus --bam");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 16);
    assert_eq!(umi_len, 10);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].ec, 0);
    assert_eq!(records[1].ec, 0);
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
fn bus_bam_trims_corrected_barcode_suffix() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("corrected_barcode.bam");
    let out_dir = dir.path().join("bam_corrected_barcode_out");

    write_bam_with_corrected_barcode(&bam, &transcript, b"ACGTACGT-1", b"TTTTTTTTTT");

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
        .expect("run kallistors bus --bam corrected barcode");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(umi_len, 10);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"ACGTACGT"));
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
fn bus_bam_reads_corrected_umi_id_tag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("corrected_umi.bam");
    let out_dir = dir.path().join("bam_corrected_umi_out");

    write_bam_with_corrected_barcode_and_umi(
        &bam,
        &transcript,
        b"ACGTACGT-1",
        b"TTTTAAAA",
        Tag::UMI_ID,
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
        .expect("run kallistors bus --bam corrected UMI");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(umi_len, 8);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"ACGTACGT"));
    assert_eq!(records[0].umi, encode_bus_seq(b"TTTTAAAA"));
}

#[test]
fn bus_bam_reads_corrected_umi_barcode_tag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("corrected_umi_barcode.bam");
    let out_dir = dir.path().join("bam_corrected_umi_barcode_out");

    write_bam_with_corrected_barcode_and_umi(
        &bam,
        &transcript,
        b"ACGTACGT-1",
        b"AAAATTTT",
        Tag::new(b'U', b'B'),
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
        .expect("run kallistors bus --bam corrected UB UMI");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(umi_len, 8);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].barcode, encode_bus_seq(b"ACGTACGT"));
    assert_eq!(records[0].umi, encode_bus_seq(b"AAAATTTT"));
}

#[test]
fn bus_bam_rejects_missing_barcode_tag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("missing_barcode.bam");
    let out_dir = dir.path().join("bam_missing_barcode_out");

    write_bam_with_optional_tags(&bam, &[(&transcript, None, Some(b"TTTTTTTTTT"))]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg(&bam)
        .output()
        .expect("run kallistors bus --bam missing barcode");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("BAM record 1 is missing CR/CB barcode tag")
    );
}

#[test]
fn bus_bam_rejects_missing_umi_tag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("missing_umi.bam");
    let out_dir = dir.path().join("bam_missing_umi_out");

    write_bam_with_optional_tags(&bam, &[(&transcript, Some(b"ACGTACGTACGTACGT"), None)]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&out_dir)
        .arg("-x")
        .arg("10XV2")
        .arg("--bam")
        .arg(&bam)
        .output()
        .expect("run kallistors bus --bam missing UMI");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("BAM record 1 is missing UR/RX/MI/UB UMI tag")
    );
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
    assert!(novel.starts_with("@unmapped\n"));
    assert!(novel.contains("@novel_disjointIntersect\n"));
    assert!(novel.contains(std::str::from_utf8(&[b'N'; 80]).unwrap()));
}

#[test]
fn bus_bam_uses_aux_tag_lengths_in_header() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bam = dir.path().join("reads.bam");
    let out_dir = dir.path().join("bam_lengths_out");

    write_bam(
        &bam,
        &[
            (&transcript, b"ACGTACGT", b"TTTTTT"),
            (&transcript, b"TGCATGCA", b"AAAAAA"),
            (&transcript, b"ACGTACGTAC", b"TTTTTTTT"),
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
        .arg("--bam")
        .arg(&bam)
        .status()
        .expect("run kallistors bus --bam");
    assert!(status.success());

    let (bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(bc_len, 8);
    assert_eq!(umi_len, 6);
    assert_eq!(records.len(), 3);
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
    assert_ne!(records[0].flags, 0);
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
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 1")
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
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 1")
    );
}

#[test]
fn bus_rejects_fastq_files_with_different_record_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    let out_dir = dir.path().join("mismatched_fastq_out");

    write_fastq(
        &r1,
        &[
            ("cell_read", b"ACGTACGTACGTACGTTTTTTTTTTTTT"),
            ("cell_read2", b"TGCATGCATGCATGCATTTTTTTTTTTT"),
        ],
    );
    write_fastq(&r2, &[("seq_read", &transcript)]);

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
        .expect("run kallistors bus with mismatched FASTQs");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("FASTQ files ended at different record counts")
    );
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
fn bus_smartseq3_default_tag_trims_umi() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bc1 = dir.path().join("bc1.fastq");
    let bc2 = dir.path().join("bc2.fastq");
    let r1 = dir.path().join("smartseq3_r1.fastq");
    let r2 = dir.path().join("smartseq3_r2.fastq");
    let out_dir = dir.path().join("smartseq3_out");

    write_fastq(&bc1, &[("bc1", b"ACGT")]);
    write_fastq(&bc2, &[("bc2", b"TGCA")]);
    let mut tagged_read = b"ATTGCGCAATG".to_vec();
    tagged_read.extend_from_slice(b"TTTTTTTT");
    tagged_read.extend_from_slice(b"AAA");
    tagged_read.extend_from_slice(&transcript);
    write_fastq(&r1, &[("r1", &tagged_read)]);
    write_fastq(&r2, &[("r2", &reverse_complement(&transcript))]);

    let output = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
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

#[test]
fn bus_smartseq3_missing_tag_uses_non_umi_sentinel() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (index, transcript) = build_tiny_index(&dir);
    let bc1 = dir.path().join("bc1.fastq");
    let bc2 = dir.path().join("bc2.fastq");
    let r1 = dir.path().join("smartseq3_no_tag_r1.fastq");
    let r2 = dir.path().join("smartseq3_no_tag_r2.fastq");
    let out_dir = dir.path().join("smartseq3_no_tag_out");

    write_fastq(&bc1, &[("bc1", b"ACGT")]);
    write_fastq(&bc2, &[("bc2", b"TGCA")]);
    write_fastq(&r1, &[("r1", &transcript)]);
    write_fastq(&r2, &[("r2", &reverse_complement(&transcript))]);

    let status = Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .arg("bus")
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
        .status()
        .expect("run kallistors bus smartseq3 missing tag");
    assert!(status.success());

    let (_bc_len, umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(umi_len, 8);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].umi, u64::MAX);
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

    let (_bc_len, _umi_len, bus_records) = read_bus_records(&out_dir.join("output.bus"));
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

    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 2")
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
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 1")
    );
    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
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
    assert!(
        fs::read_to_string(out_dir.join("run_info.json"))
            .unwrap()
            .contains("\"n_processed\": 1")
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("matrix.cells")).unwrap(),
        "sample_a\n"
    );
    let (_bc_len, _umi_len, records) = read_bus_records(&out_dir.join("output.bus"));
    assert_eq!(records.len(), 1);
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
    assert!(run_info.contains("\"n_frame_clashes\": 2"));
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
