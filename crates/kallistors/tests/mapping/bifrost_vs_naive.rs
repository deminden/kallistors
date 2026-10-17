use std::collections::HashMap;
use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    EcCounts, build_bifrost_index, build_kmer_ec_index, pseudoalign_single_end,
    pseudoalign_single_end_bifrost,
};

use crate::helpers::{self, synthesize_reads, to_fastq};

fn normalize_ec(counts: &EcCounts) -> HashMap<Vec<u32>, u32> {
    let mut map = HashMap::new();
    for (ec, &count) in counts.ec_list.iter().zip(counts.counts.iter()) {
        let mut key = ec.clone();
        key.sort_unstable();
        map.insert(key, count);
    }
    map
}

#[test]
fn bifrost_matches_naive_ecs() {
    let dataset = helpers::build_synthetic_index();

    let reads = synthesize_reads(&dataset.transcripts, 60, 24);
    let fastq = to_fastq(&reads);

    let naive_index = build_kmer_ec_index(&dataset.index_path).expect("naive index");
    let bifrost_index = build_bifrost_index(&dataset.index_path).expect("bifrost index");

    let mut reader = FastqReader::new(Cursor::new(fastq.clone()));
    let naive_counts = pseudoalign_single_end(&naive_index, &mut reader).expect("naive counts");

    let mut reader = FastqReader::new(Cursor::new(fastq));
    let bifrost_counts =
        pseudoalign_single_end_bifrost(&bifrost_index, &mut reader).expect("bifrost counts");

    assert_eq!(naive_counts.reads_processed, 24);
    assert_eq!(naive_counts.reads_aligned, 24);
    assert!(naive_counts.ec_list.iter().any(|ec| ec.len() == 1));
    assert!(naive_counts.ec_list.iter().any(|ec| ec.len() > 1));
    assert_eq!(naive_counts.reads_processed, bifrost_counts.reads_processed);
    assert_eq!(naive_counts.reads_aligned, bifrost_counts.reads_aligned);
    assert_eq!(normalize_ec(&naive_counts), normalize_ec(&bifrost_counts));
}
