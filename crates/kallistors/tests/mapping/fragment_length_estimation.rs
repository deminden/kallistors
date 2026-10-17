use crate::helpers::{synthesize_paired_reads, to_fastq_pairs};
use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index_with_positions,
    pseudoalign_paired_bifrost_with_options,
};
use std::io::Cursor;

use crate::helpers;

#[test]
fn paired_fragment_length_estimate_matches_known() {
    let dataset = helpers::build_synthetic_index();
    let index =
        build_bifrost_index_with_positions(&dataset.index_path, true).expect("bifrost index");
    let frag_len = 90usize;
    let read_len = 50usize;
    let pairs = synthesize_paired_reads(&dataset.transcripts, read_len, frag_len, 12);
    let (fq1, fq2) = to_fastq_pairs(&pairs);

    let mut r1 = FastqReader::new(Cursor::new(fq1));
    let mut r2 = FastqReader::new(Cursor::new(fq2));
    let counts = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut r1,
        &mut r2,
        Strand::Unstranded,
        PseudoalignOptions::default(),
    )
    .expect("paired pseudoalign");

    let stats = counts.fragment_length_stats.expect("fragment length stats");
    let mean = stats.mean().expect("mean");
    assert!((mean - frag_len as f64).abs() < 1.0);
}

#[test]
fn sampling_preserves_every_ec_across_parser_chunks() {
    use kallistors::pseudoalign::pseudoalign_paired_bifrost_with_options_threaded;
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index_with_positions(&dataset.index_path, true).unwrap();
    let pairs = synthesize_paired_reads(&dataset.transcripts, 50, 90, 25_000);
    let (fq1, fq2) = to_fastq_pairs(&pairs);
    let run = |threads, limit| {
        pseudoalign_paired_bifrost_with_options_threaded(
            &index,
            &mut FastqReader::new(Cursor::new(&fq1)),
            &mut FastqReader::new(Cursor::new(&fq2)),
            Strand::Unstranded,
            PseudoalignOptions {
                fragment_length_samples: limit,
                ..Default::default()
            },
            threads,
        )
        .unwrap()
    };
    let serial = run(1, 37);
    // The sequence-only serial parser must preserve the record API's exact
    // observation order, quota, histogram and first-occurrence EC ordering.
    let record_serial = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut FastqReader::new(Cursor::new(&fq1)),
        &mut FastqReader::new(Cursor::new(&fq2)),
        Strand::Unstranded,
        PseudoalignOptions {
            fragment_length_samples: 37,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(serial.ec_list, record_serial.ec_list);
    assert_eq!(serial.counts, record_serial.counts);
    assert_eq!(
        serial.fragment_length_hist,
        record_serial.fragment_length_hist
    );
    let parallel = run(4, 37);
    let all = run(4, 0);
    let canonical = |c: &kallistors::pseudoalign::EcCounts| {
        c.ec_list
            .iter()
            .cloned()
            .zip(c.counts.iter().copied())
            .collect::<std::collections::BTreeMap<_, _>>()
    };
    assert_eq!(canonical(&serial), canonical(&all));
    assert_eq!(canonical(&parallel), canonical(&all));
    assert_eq!(serial.reads_processed, pairs.len() as u64);
    assert_eq!(parallel.reads_processed, all.reads_processed);
    let a = serial.fragment_length_stats.unwrap();
    let b = parallel.fragment_length_stats.unwrap();
    assert_eq!(a.count(), 37);
    // These 25,000 pairs fit in one upstream 8 MiB sequence buffer, even
    // though the parser returns several 12,000-pair chunks. Only one sampling
    // quota applies; parser chunks must not become concurrent sample buffers.
    assert_eq!(b.count(), 37);
    assert_eq!(
        parallel
            .fragment_length_hist
            .unwrap()
            .iter()
            .map(|&c| u64::from(c))
            .sum::<u64>(),
        b.count()
    );
    assert!((a.mean().unwrap() - b.mean().unwrap()).abs() < 1.0);
    assert!(all.fragment_length_stats.unwrap().count() > 37);
}
