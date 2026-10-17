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

// SRR10919051.8026841 and .8028299, fasterp-trimmed; expected complete ECs
// independently exported by original kallisto v0.52.0 with Pachter Lab v1
// mature + nascent FASTAs concatenated in that order. No introns are synthesized.
#[test]
#[ignore = "requires the official Pachter Lab human NAC index"]
fn paired_nac_ambiguous_jump_probes_keep_complete_upstream_ecs() {
    use kallistors::io::ReadSource;
    use kallistors::pseudoalign::{
        build_bifrost_index_quant_packed_threaded, trace_read_pair_bifrost,
    };
    use std::collections::BTreeMap;

    let path = std::env::var_os("KALLISTORS_TEST_NAC_INDEX")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            workspace_root().join("data/official-human-nac/human_nac.mature-plus-nascent.idx")
        });
    assert!(
        path.is_file(),
        "prepare the official Pachter Lab human NAC index first"
    );
    let index =
        build_bifrost_index_quant_packed_threaded(&path, 8).expect("load official NAC index");
    let reads1 = include_bytes!("../fixtures/paired_nac_ambiguous_1.fixture");
    let reads2 = include_bytes!("../fixtures/paired_nac_ambiguous_2.fixture");
    let expected = [
        vec![18719, 18720, 18721, 18722, 18732, 231083],
        vec![251846, 251888],
    ];
    for threads in [1, 8] {
        let mut left = FastqReader::new(Cursor::new(reads1));
        let mut right = FastqReader::new(Cursor::new(reads2));
        let counts = pseudoalign_paired_bifrost_with_options_threaded(
            &index,
            &mut left,
            &mut right,
            Strand::Unstranded,
            PseudoalignOptions::default(),
            threads,
        )
        .expect("map ambiguous paired reads");
        assert_eq!(counts.reads_processed, 2);
        assert_eq!(counts.reads_aligned, 2);
        let actual: BTreeMap<_, _> = counts.ec_list.into_iter().zip(counts.counts).collect();
        let oracle: BTreeMap<_, _> = expected.iter().cloned().map(|ec| (ec, 1)).collect();
        assert_eq!(actual, oracle, "threads={threads}");
    }
    let mut left = FastqReader::new(Cursor::new(reads1));
    let mut right = FastqReader::new(Cursor::new(reads2));
    for ec in expected {
        let a = left.next_record().unwrap().unwrap();
        let b = right.next_record().unwrap().unwrap();
        for fast in [false, true] {
            let mut options = PseudoalignOptions::default();
            options.investigation.trace_fast_path = fast;
            let trace =
                trace_read_pair_bifrost(&index, &a.seq, &b.seq, Strand::Unstranded, options);
            assert_eq!(
                trace.merged_ec.as_deref(),
                Some(ec.as_slice()),
                "fast={fast}"
            );
            assert!(!trace.hard_reject_pair);
        }
    }
}
