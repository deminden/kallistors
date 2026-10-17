use crate::helpers::{synthesize_paired_reads, to_fastq_pairs};
use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    Strand, build_bifrost_index, build_kmer_ec_index, pseudoalign_paired_bifrost_with_strand,
    pseudoalign_paired_naive,
};
use std::io::Cursor;

use crate::helpers;

#[test]
fn paired_parity_naive_vs_bifrost() {
    let dataset = helpers::build_synthetic_index();
    let transcripts = dataset.transcripts;
    let pairs = synthesize_paired_reads(&transcripts, 50, 90, 64);
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
    let ec_map = |counts: &kallistors::pseudoalign::EcCounts| {
        counts
            .ec_list
            .iter()
            .cloned()
            .zip(counts.counts.iter().copied())
            .collect::<std::collections::BTreeMap<_, _>>()
    };
    assert!(naive.ec_list.iter().any(|ec| ec.len() > 1));
    assert!(naive.ec_list.iter().any(|ec| ec.len() == 1));
    assert_eq!(ec_map(&naive), ec_map(&bifrost));
}

#[test]
fn paired_allows_one_unmapped_mate() {
    let dataset = helpers::build_synthetic_index();
    let seq = &dataset.transcripts[0];
    let read_len = 50usize;
    let read1 = seq[..read_len].to_vec();
    let read2 = vec![b'N'; read_len];
    let (fq1, fq2) = to_fastq_pairs(&[(read1, read2)]);

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

    assert_eq!(naive.reads_processed, 1);
    assert_eq!(naive.reads_aligned, 1);
    assert_eq!(bifrost.reads_processed, 1);
    assert_eq!(bifrost.reads_aligned, 1);
    assert_eq!(naive.ec_list, bifrost.ec_list);
    assert_eq!(naive.counts, bifrost.counts);
}
