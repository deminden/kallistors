use std::io::Cursor;
use std::path::PathBuf;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index_with_positions_threaded,
    pseudoalign_paired_bifrost_with_options, pseudoalign_paired_bifrost_with_options_threaded,
    unique_pseudoaligned_reads,
};

const READS1: &[u8] = include_bytes!("../fixtures/paired_backoff_regression_1.fixture");
const READS2: &[u8] = include_bytes!("../fixtures/paired_backoff_regression_2.fixture");

#[test]
#[ignore = "requires the downloaded GENCODE v50 full index"]
fn paired_backoff_regression_reads_keep_kallisto_counts() {
    let index_path = workspace_root().join("data/srr-small/gencode.v50.kallisto.idx");
    assert!(index_path.is_file(), "prepare the GENCODE v50 index first");
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
    assert_eq!(counts.counts, [1]);
    assert_eq!(
        counts.ec_list,
        [vec![
            240838, 240839, 240840, 240841, 240844, 240845, 240846, 240847, 240851, 240855, 240856,
            240857, 240858, 240860, 240861, 240862, 240864, 240865, 240866, 240869, 240870, 240871,
            240872, 240873, 240874, 240875, 240876, 240877, 240878, 240879, 240881, 240882, 240883,
            240884, 240885, 240886, 240887, 240888, 240889, 240890, 240891
        ]]
    );
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|path| path.parent())
        .expect("workspace root")
        .to_path_buf()
}
