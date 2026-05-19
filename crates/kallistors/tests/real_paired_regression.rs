use std::io::Cursor;
use std::path::PathBuf;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index_with_positions_threaded,
    pseudoalign_paired_bifrost_with_options, pseudoalign_paired_bifrost_with_options_threaded,
    unique_pseudoaligned_reads,
};

const READS1: &[u8] = include_bytes!("fixtures/paired_backoff_regression_1.fixture");
const READS2: &[u8] = include_bytes!("fixtures/paired_backoff_regression_2.fixture");
const RUN_ENV: &str = "KALLISTORS_RUN_REAL_PAIRED_REGRESSION";

#[test]
fn paired_backoff_regression_reads_keep_kallisto_counts() {
    if std::env::var_os(RUN_ENV).is_none() {
        return;
    }
    let index_path = workspace_root().join("data/gencode.v49_kallisto.idx");
    if !index_path.exists() {
        return;
    }
    let index = build_bifrost_index_with_positions_threaded(&index_path, false, 2)
        .expect("load gencode index");

    let mut left = FastqReader::new(Cursor::new(READS1));
    let mut right = FastqReader::new(Cursor::new(READS2));
    let counts = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut left,
        &mut right,
        Strand::Unstranded,
        PseudoalignOptions::default(),
    )
    .expect("single-thread paired pseudoalign");
    assert_regression_counts(&counts);

    let mut left = FastqReader::new(Cursor::new(READS1));
    let mut right = FastqReader::new(Cursor::new(READS2));
    let counts = pseudoalign_paired_bifrost_with_options_threaded(
        &index,
        &mut left,
        &mut right,
        Strand::Unstranded,
        PseudoalignOptions::default(),
        2,
    )
    .expect("threaded paired pseudoalign");
    assert_regression_counts(&counts);
}

fn assert_regression_counts(counts: &kallistors::pseudoalign::EcCounts) {
    assert_eq!(counts.reads_processed, 2);
    assert_eq!(counts.reads_aligned, 1);
    assert_eq!(unique_pseudoaligned_reads(counts), 0);
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|path| path.parent())
        .expect("workspace root")
        .to_path_buf()
}
