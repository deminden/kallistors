use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    Strand, build_bifrost_index, build_kmer_ec_index, pseudoalign_paired_bifrost_with_strand,
    pseudoalign_paired_naive,
};
use std::io::Cursor;

mod helpers;

#[test]
fn paired_parity_naive_vs_bifrost() {
    let dataset = match helpers::build_synthetic_index() {
        Some(dataset) => dataset,
        None => return,
    };
    let transcripts = dataset.transcripts;
    let pairs = synthesize_paired_reads(&transcripts, 50, 90, 8);
    let (fq1, fq2) = to_fastq_pairs(&pairs);

    let naive_index = build_kmer_ec_index(&dataset.index_path).expect("naive index");
    let bifrost_index = build_bifrost_index(&dataset.index_path).expect("bifrost index");

    let mut r1 = FastqReader::new(Cursor::new(fq1.clone()));
    let mut r2 = FastqReader::new(Cursor::new(fq2.clone()));
    let naive = pseudoalign_paired_naive(&naive_index, &mut r1, &mut r2).expect("naive");

    let mut r1 = FastqReader::new(Cursor::new(fq1));
    let mut r2 = FastqReader::new(Cursor::new(fq2));
    let bifrost = pseudoalign_paired_bifrost_with_strand(
        &bifrost_index,
        &mut r1,
        &mut r2,
        Strand::Unstranded,
    )
    .expect("bifrost");

    assert_eq!(naive.reads_processed, bifrost.reads_processed);
    assert_eq!(naive.reads_aligned, bifrost.reads_aligned);
    assert_eq!(naive.ec_list.len(), bifrost.ec_list.len());
}

fn synthesize_paired_reads(
    transcripts: &[Vec<u8>],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut out = Vec::with_capacity(count);
    let mut idx = 0usize;
    let mut attempts = 0usize;
    while out.len() < count && attempts < count * 10 {
        let seq = &transcripts[idx % transcripts.len()];
        idx += 1;
        attempts += 1;
        if seq.len() < frag_len {
            continue;
        }
        let start = (idx * 7) % (seq.len() - frag_len + 1);
        let frag = &seq[start..start + frag_len];
        let read1 = frag[..read_len].to_vec();
        let read2 = revcomp(&frag[frag_len - read_len..]);
        out.push((read1, read2));
    }
    out
}

fn to_fastq_pairs(pairs: &[(Vec<u8>, Vec<u8>)]) -> (Vec<u8>, Vec<u8>) {
    let mut fq1 = Vec::new();
    let mut fq2 = Vec::new();
    for (i, (r1, r2)) in pairs.iter().enumerate() {
        fq1.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        fq1.extend_from_slice(r1);
        fq1.extend_from_slice(b"\n+\n");
        fq1.extend_from_slice(vec![b'I'; r1.len()].as_slice());
        fq1.extend_from_slice(b"\n");

        fq2.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        fq2.extend_from_slice(r2);
        fq2.extend_from_slice(b"\n+\n");
        fq2.extend_from_slice(vec![b'I'; r2.len()].as_slice());
        fq2.extend_from_slice(b"\n");
    }
    (fq1, fq2)
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
