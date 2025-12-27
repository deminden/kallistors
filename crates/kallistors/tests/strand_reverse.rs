use std::fs::File;
use std::io::{Cursor, Read};
use std::path::PathBuf;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    Strand, build_bifrost_index, pseudoalign_single_end_bifrost_with_strand,
};

fn fixture_path(name: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests");
    path.push("fixtures");
    path.push(name);
    path
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        });
    }
    out
}

#[test]
fn reverse_strand_aligns_only_with_revcomp() {
    let index_path = fixture_path("synthetic.idx");
    let fasta_path = fixture_path("synthetic_transcripts.fa");
    if !index_path.exists() || !fasta_path.exists() {
        return;
    }

    let index = build_bifrost_index(&index_path).expect("index");

    let mut fasta = String::new();
    File::open(&fasta_path)
        .expect("fasta file")
        .read_to_string(&mut fasta)
        .expect("read fasta");
    let transcripts = parse_fasta(&fasta);
    let read_len = index.k + 10;
    let mut read = None;
    'outer: for seq in &transcripts {
        if seq.len() < read_len {
            continue;
        }
        for start in 0..=seq.len() - read_len {
            let slice = &seq[start..start + read_len];
            let rc = revcomp(slice);
            if !transcripts.iter().any(|t| contains_subseq(t, &rc)) {
                read = Some(slice.to_vec());
                break 'outer;
            }
        }
    }
    let read = read.expect("find non-palindromic read");
    let rc = revcomp(&read);
    let mut fastq = Vec::new();
    fastq.extend_from_slice(b"@rc\n");
    fastq.extend_from_slice(&rc);
    fastq.extend_from_slice(b"\n+\n");
    fastq.extend_from_slice(&vec![b'I'; rc.len()]);
    fastq.extend_from_slice(b"\n");

    let mut rc_reader = FastqReader::new(Cursor::new(fastq.clone()));
    let res_forward =
        pseudoalign_single_end_bifrost_with_strand(&index, &mut rc_reader, Strand::Forward)
            .expect("forward");
    assert_eq!(res_forward.reads_aligned, 0);

    let mut rc_reader = FastqReader::new(Cursor::new(fastq.clone()));
    let res_reverse =
        pseudoalign_single_end_bifrost_with_strand(&index, &mut rc_reader, Strand::Reverse)
            .expect("reverse");
    assert_eq!(res_reverse.reads_aligned, 1);

    let mut rc_reader = FastqReader::new(Cursor::new(fastq));
    let res_unstranded =
        pseudoalign_single_end_bifrost_with_strand(&index, &mut rc_reader, Strand::Unstranded)
            .expect("unstranded");
    assert_eq!(res_unstranded.reads_aligned, 1);
}

fn parse_fasta(input: &str) -> Vec<Vec<u8>> {
    let mut seqs = Vec::new();
    let mut current = Vec::new();
    for line in input.lines() {
        if line.starts_with('>') {
            if !current.is_empty() {
                seqs.push(current);
                current = Vec::new();
            }
            continue;
        }
        current.extend_from_slice(line.trim().as_bytes());
    }
    if !current.is_empty() {
        seqs.push(current);
    }
    seqs
}

fn contains_subseq(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}
