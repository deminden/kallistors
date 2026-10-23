pub mod helpers;
use helpers::{command_status_with_timeout, json_u64};
use std::fs;
use std::io::Write as _;
use std::process::{Command, Stdio};
use std::time::Duration;

use kallistors::index::{IndexBuildOptions, build_index};
use noodles_sam::alignment::{
    RecordBuf,
    io::Write as _,
    record::{Flags, cigar::op::Kind, data::field::Tag},
    record_buf::{QualityScores, Sequence, data::field::Value},
};

type BamRecordFixture<'a> = (&'a [u8], &'a [u8], &'a [u8], Flags);
type OptionalBamTagsFixture<'a> = (&'a [u8], Option<&'a [u8]>, Option<&'a [u8]>);
const TINY_TRANSCRIPT: &[u8] =
    b"ACGTGCACTGATCGTACGATCGTACGTTAGCTAGCTAGGCTAGCATCGATCGATGCTAGCTAGCTGACT";

fn bus_command(index: &std::path::Path, out: &std::path::Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_kallistors"));
    command.arg("bus").arg("-i").arg(index).arg("-o").arg(out);
    command
}

fn assert_bus_success(output: &std::process::Output, case: &str) {
    assert!(
        output.status.success(),
        "{case}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_bus_error(output: &std::process::Output, case: &str, error: &str) {
    assert!(!output.status.success(), "{case}: expected failure");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(error),
        "{case}: expected {error:?}: {stderr}"
    );
}

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

fn write_bam_with_quality(
    path: &std::path::Path,
    seq: &[u8],
    barcode: &[u8],
    umi: &[u8],
    quality: &[u8],
) {
    let mut writer = noodles_bam::io::Writer::new(Vec::new());
    let header = noodles_sam::Header::default();
    writer.write_header(&header).expect("write BAM header");
    let barcode = std::str::from_utf8(barcode).expect("ASCII barcode");
    let umi = std::str::from_utf8(umi).expect("ASCII UMI");
    let record = RecordBuf::builder()
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quality.to_vec()))
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
    build_index_with_transcript(dir, TINY_TRANSCRIPT)
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
    let fasta = write_transcript_fasta(dir, name, transcript);
    let index = dir.path().join("transcripts.idx");

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

fn write_transcript_fasta(
    dir: &tempfile::TempDir,
    name: &str,
    transcript: &[u8],
) -> std::path::PathBuf {
    let fasta = dir.path().join("transcripts.fa");
    fs::write(
        &fasta,
        format!(">{name}\n{}\n", std::str::from_utf8(transcript).unwrap()),
    )
    .expect("write fasta");
    fasta
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

#[path = "bus/bam_input.rs"]
mod bam_input;
#[path = "bus/bam_output.rs"]
mod bam_output;
#[path = "bus/batch.rs"]
mod batch;
#[path = "bus/input.rs"]
mod input;
#[path = "bus/layouts.rs"]
mod layouts;
#[path = "bus/long_reads.rs"]
mod long_reads;
#[path = "bus/mapping_output.rs"]
mod mapping_output;
#[path = "bus/tags.rs"]
mod tags;
