use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    EcCounts, build_bifrost_index, build_kmer_ec_index, pseudoalign_single_end,
    pseudoalign_single_end_bifrost,
};

fn fixture_path(name: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests");
    path.push("fixtures");
    path.push(name);
    path
}

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
    let index_path = fixture_path("synthetic.idx");
    let reads_path = fixture_path("synthetic_reads.fq");

    let naive_index = build_kmer_ec_index(&index_path).expect("naive index");
    let bifrost_index = build_bifrost_index(&index_path).expect("bifrost index");

    let reads_file = File::open(&reads_path).expect("reads file");
    let mut reads = FastqReader::new(BufReader::new(reads_file));
    let naive_counts = pseudoalign_single_end(&naive_index, &mut reads).expect("naive counts");

    let reads_file = File::open(&reads_path).expect("reads file");
    let mut reads = FastqReader::new(BufReader::new(reads_file));
    let bifrost_counts =
        pseudoalign_single_end_bifrost(&bifrost_index, &mut reads).expect("bifrost counts");

    assert_eq!(naive_counts.reads_processed, bifrost_counts.reads_processed);
    assert_eq!(naive_counts.reads_aligned, bifrost_counts.reads_aligned);
    assert_eq!(normalize_ec(&naive_counts), normalize_ec(&bifrost_counts));
}
