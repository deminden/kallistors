use crate::helpers::{synthesize_paired_reads, synthesize_reads, to_fastq, to_fastq_pairs};
use std::collections::HashMap;
use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index, pseudoalign_paired_bifrost_with_options,
    pseudoalign_paired_bifrost_with_options_threaded, pseudoalign_single_end_bifrost_with_options,
    pseudoalign_single_end_bifrost_with_options_threaded,
};

use crate::helpers;

#[test]
fn completed_bias_sampling_preserves_ecs_samples_and_fragment_histograms() {
    use kallistors::pseudoalign::{
        StrandSpecific, build_bifrost_index_with_positions,
        pseudoalign_paired_bifrost_debug_with_options,
        pseudoalign_single_end_bifrost_debug_with_options,
    };
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index_with_positions(&dataset.index_path, true).unwrap();
    let reads = to_fastq(&synthesize_reads(&dataset.transcripts, 50, 512));
    let (left, right) = to_fastq_pairs(&synthesize_paired_reads(&dataset.transcripts, 50, 90, 512));
    for mode in [None, Some(StrandSpecific::FR), Some(StrandSpecific::RF)] {
        for limit in [0, 1, 3] {
            let options = PseudoalignOptions {
                bias: true,
                max_bias: limit,
                fragment_length_samples: 1,
                strand_specific: mode,
                ..Default::default()
            };
            let (expected, _) = pseudoalign_single_end_bifrost_debug_with_options(
                &index,
                &mut FastqReader::new(Cursor::new(&reads)),
                Strand::Unstranded,
                0,
                None,
                options,
            )
            .unwrap();
            for threads in [1, 8] {
                let got = pseudoalign_single_end_bifrost_with_options_threaded(
                    &index,
                    &mut FastqReader::new(Cursor::new(&reads)),
                    Strand::Unstranded,
                    None,
                    options,
                    threads,
                )
                .unwrap();
                assert_eq!(to_ec_map(&got), to_ec_map(&expected));
                assert_eq!(got.reads_processed, expected.reads_processed);
                assert_eq!(
                    got.bias.as_ref().unwrap().counts,
                    expected.bias.as_ref().unwrap().counts
                );
                assert_eq!(
                    got.bias.as_ref().unwrap().total,
                    expected.bias.as_ref().unwrap().total
                );
            }
            let (expected, _) = pseudoalign_paired_bifrost_debug_with_options(
                &index,
                &mut FastqReader::new(Cursor::new(&left)),
                &mut FastqReader::new(Cursor::new(&right)),
                Strand::Unstranded,
                0,
                options,
            )
            .unwrap();
            let got = pseudoalign_paired_bifrost_with_options(
                &index,
                &mut FastqReader::new(Cursor::new(&left)),
                &mut FastqReader::new(Cursor::new(&right)),
                Strand::Unstranded,
                options,
            )
            .unwrap();
            assert_eq!(to_ec_map(&got), to_ec_map(&expected));
            assert_eq!(got.fragment_length_hist, expected.fragment_length_hist);
            assert_eq!(
                got.bias.as_ref().unwrap().counts,
                expected.bias.as_ref().unwrap().counts
            );
            assert_eq!(
                got.bias.as_ref().unwrap().total,
                expected.bias.as_ref().unwrap().total
            );
            if mode.is_none() {
                assert_eq!(got.bias.unwrap().total, limit as u64);
            }
        }
    }
}

#[test]
fn threaded_single_end_matches_single_thread() {
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index(&dataset.index_path).expect("bifrost index");
    let mut reads = synthesize_reads(&dataset.transcripts, 50, 25_001);
    for read in reads.iter_mut().step_by(17) {
        read.fill(b'N');
    }
    let fq = to_fastq(&reads);

    let mut reader = FastqReader::new(Cursor::new(fq.clone()));
    let single = pseudoalign_single_end_bifrost_with_options(
        &index,
        &mut reader,
        Strand::Unstranded,
        None,
        PseudoalignOptions::default(),
    )
    .expect("single-thread");

    assert_eq!(single.reads_processed, 25_001);
    assert!(single.reads_aligned > 0 && single.reads_aligned < single.reads_processed);
    assert!(single.ec_list.iter().any(|ec| ec.len() == 1));
    assert!(single.ec_list.iter().any(|ec| ec.len() > 1));
    for threads in [0, 1, 4] {
        let mut reader = FastqReader::new(Cursor::new(&fq));
        let threaded = pseudoalign_single_end_bifrost_with_options_threaded(
            &index,
            &mut reader,
            Strand::Unstranded,
            None,
            PseudoalignOptions::default(),
            threads,
        )
        .expect("threaded");
        assert_eq!(single.reads_processed, threaded.reads_processed);
        assert_eq!(single.reads_aligned, threaded.reads_aligned);
        assert_eq!(to_ec_map(&single), to_ec_map(&threaded));
    }
}

#[test]
fn threaded_paired_matches_single_thread() {
    use kallistors::pseudoalign::{
        StrandSpecific, build_bifrost_index_with_positions,
        pseudoalign_paired_bifrost_debug_with_options,
    };
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index_with_positions(&dataset.index_path, true).unwrap();
    let mut pairs = synthesize_paired_reads(&dataset.transcripts, 50, 90, 25_001);
    for (left, right) in pairs.iter_mut().step_by(17) {
        left.fill(b'N');
        right.fill(b'N');
    }
    for (left, _) in pairs.iter_mut().step_by(31) {
        left.truncate(index.k.saturating_sub(1));
    }
    for (left, right) in pairs.iter_mut().step_by(37) {
        std::mem::swap(left, right);
    }
    let (fq1, fq2) = to_fastq_pairs(&pairs);
    // Cross serial parser batches after the first fragment observation. The
    // diagnostic oracle keeps complete scans, including strand rejection.
    for mode in [None, Some(StrandSpecific::FR), Some(StrandSpecific::RF)] {
        for fragment_length_samples in [0, 1] {
            let options = PseudoalignOptions {
                fragment_length_samples,
                strand_specific: mode,
                ..Default::default()
            };
            let (expected, _) = pseudoalign_paired_bifrost_debug_with_options(
                &index,
                &mut FastqReader::new(Cursor::new(&fq1)),
                &mut FastqReader::new(Cursor::new(&fq2)),
                Strand::Unstranded,
                0,
                options,
            )
            .unwrap();
            assert_eq!(expected.reads_processed, 25_001);
            assert!(
                expected.reads_aligned > 0 && expected.reads_aligned < expected.reads_processed
            );
            if mode.is_none() {
                assert!(expected.ec_list.iter().any(|ec| ec.len() == 1));
                assert!(expected.ec_list.iter().any(|ec| ec.len() > 1));
            }
            for threads in [1, 8] {
                let got = pseudoalign_paired_bifrost_with_options_threaded(
                    &index,
                    &mut FastqReader::new(Cursor::new(&fq1)),
                    &mut FastqReader::new(Cursor::new(&fq2)),
                    Strand::Unstranded,
                    options,
                    threads,
                )
                .unwrap();
                assert_eq!(expected.reads_processed, got.reads_processed);
                assert_eq!(expected.reads_aligned, got.reads_aligned);
                assert_eq!(to_ec_map(&expected), to_ec_map(&got));
            }
        }
    }
}

#[test]
fn threaded_input_errors_release_backpressured_producers() {
    use std::fs::{self, File};
    use std::process::{Command, Stdio};
    use std::time::Duration;

    let dataset = helpers::build_synthetic_index();
    let dir = tempfile::tempdir().unwrap();
    let long = dir.path().join("long.fq");
    let empty = dir.path().join("empty.fq");
    let log = dir.path().join("stderr");
    // Exceeds the bounded queue. A subprocess watchdog turns a deadlock into a
    // bounded failure and reaps the process instead of leaving blocked threads.
    std::fs::write(&long, to_fastq(&vec![vec![b'A'; 50]; 60_000])).unwrap();
    std::fs::write(&empty, []).unwrap();
    for (left, right) in [(&long, &empty), (&empty, &long)] {
        let status = helpers::command_status_with_timeout(
            Command::new(env!("CARGO_BIN_EXE_kallistors"))
                .arg("pseudoalign")
                .arg("--index")
                .arg(&dataset.index_path)
                .arg("--reads")
                .arg(left)
                .arg("--reads2")
                .arg(right)
                .arg("--out")
                .arg(dir.path().join("counts.tsv"))
                .args(["--threads", "4", "--fragment-length-samples", "0"])
                .stdout(Stdio::null())
                .stderr(File::create(&log).unwrap()),
            Duration::from_secs(10),
        )
        .expect("unequal mate counts must fail without blocking producer join");
        assert!(!status.success());
        assert!(
            fs::read_to_string(&log)
                .unwrap()
                .contains("paired FASTQ length mismatch")
        );
    }
    let broken = dir.path().join("broken.fq");
    let mut truncated = fs::read(&long).unwrap();
    truncated.extend_from_slice(b"@broken\nACGT\n");
    fs::write(&broken, truncated).unwrap();
    let status = helpers::command_status_with_timeout(
        Command::new(env!("CARGO_BIN_EXE_kallistors"))
            .arg("pseudoalign")
            .arg("--index")
            .arg(&dataset.index_path)
            .arg("--reads")
            .arg(&broken)
            .arg("--out")
            .arg(dir.path().join("single.tsv"))
            .args(["--threads", "4"])
            .stdout(Stdio::null())
            .stderr(File::create(&log).unwrap()),
        Duration::from_secs(10),
    )
    .expect("single-end input errors must propagate after queued batches drain");
    assert!(!status.success());
    assert!(fs::read_to_string(log).unwrap().contains("truncated FASTQ"));
}

fn to_ec_map(counts: &kallistors::pseudoalign::EcCounts) -> HashMap<Vec<u32>, u32> {
    let mut map = HashMap::new();
    for (idx, ec) in counts.ec_list.iter().enumerate() {
        let count = counts.counts.get(idx).copied().unwrap_or(0);
        map.insert(ec.clone(), count);
    }
    map
}
