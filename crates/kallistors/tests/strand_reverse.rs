use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    Strand, build_bifrost_index, pseudoalign_single_end_bifrost_with_strand,
};

mod helpers;

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
    let dataset = match helpers::build_synthetic_index() {
        Some(dataset) => dataset,
        None => return,
    };
    let index = build_bifrost_index(&dataset.index_path).expect("index");
    let transcripts = dataset.transcripts;
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

fn contains_subseq(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}
