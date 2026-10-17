use crate::helpers::revcomp;
use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    Strand, build_bifrost_index, pseudoalign_single_end_bifrost_with_strand,
};

use crate::helpers;

#[test]
fn reverse_strand_aligns_only_with_revcomp() {
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index(&dataset.index_path).expect("index");
    assert_eq!(index.transcript_names.len(), dataset.transcripts.len());
    // Strand restricts lookup relative to stored unitigs. Their orientation is
    // chosen by graph construction and need not match the input transcript.
    let unitigs = index
        .unitigs
        .iter()
        .chain(&index.km_unitigs)
        .map(Vec::as_slice)
        .collect::<Vec<_>>();
    let read_len = index.k + 10;
    let mut read = None;
    'outer: for seq in &unitigs {
        if seq.len() < read_len {
            continue;
        }
        for start in 0..=seq.len() - read_len {
            let slice = &seq[start..start + read_len];
            let rc = revcomp(slice);
            if rc
                .windows(index.k)
                .all(|kmer| !unitigs.iter().any(|unitig| contains_subseq(unitig, kmer)))
            {
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
