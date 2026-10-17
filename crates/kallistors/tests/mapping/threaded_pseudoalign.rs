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

    let mut reader = FastqReader::new(Cursor::new(fq));
    let threaded = pseudoalign_single_end_bifrost_with_options_threaded(
        &index,
        &mut reader,
        Strand::Unstranded,
        None,
        PseudoalignOptions::default(),
        4,
    )
    .expect("threaded");

    assert_eq!(single.reads_processed, 25_001);
    assert!(single.reads_aligned > 0 && single.reads_aligned < single.reads_processed);
    assert!(single.ec_list.iter().any(|ec| ec.len() == 1));
    assert!(single.ec_list.iter().any(|ec| ec.len() > 1));
    assert_eq!(single.reads_processed, threaded.reads_processed);
    assert_eq!(single.reads_aligned, threaded.reads_aligned);
    assert_eq!(to_ec_map(&single), to_ec_map(&threaded));
}

#[test]
fn threaded_paired_matches_single_thread() {
    let dataset = helpers::build_synthetic_index();
    let index = build_bifrost_index(&dataset.index_path).expect("bifrost index");
    let mut pairs = synthesize_paired_reads(&dataset.transcripts, 50, 90, 25_001);
    for (left, right) in pairs.iter_mut().step_by(17) {
        left.fill(b'N');
        right.fill(b'N');
    }
    let options = PseudoalignOptions {
        fragment_length_samples: 0,
        ..Default::default()
    };
    let (fq1, fq2) = to_fastq_pairs(&pairs);

    let mut r1 = FastqReader::new(Cursor::new(fq1.clone()));
    let mut r2 = FastqReader::new(Cursor::new(fq2.clone()));
    let single = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut r1,
        &mut r2,
        Strand::Unstranded,
        options,
    )
    .expect("single-thread");

    let mut r1 = FastqReader::new(Cursor::new(fq1));
    let mut r2 = FastqReader::new(Cursor::new(fq2));
    let threaded = pseudoalign_paired_bifrost_with_options_threaded(
        &index,
        &mut r1,
        &mut r2,
        Strand::Unstranded,
        options,
        4,
    )
    .expect("threaded");

    assert_eq!(single.reads_processed, 25_001);
    assert!(single.reads_aligned > 0 && single.reads_aligned < single.reads_processed);
    assert!(single.ec_list.iter().any(|ec| ec.len() == 1));
    assert!(single.ec_list.iter().any(|ec| ec.len() > 1));
    assert_eq!(single.reads_processed, threaded.reads_processed);
    assert_eq!(single.reads_aligned, threaded.reads_aligned);
    assert_eq!(to_ec_map(&single), to_ec_map(&threaded));
}

#[test]
fn paired_length_error_releases_backpressured_producer() {
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
}

fn to_ec_map(counts: &kallistors::pseudoalign::EcCounts) -> HashMap<Vec<u32>, u32> {
    let mut map = HashMap::new();
    for (idx, ec) in counts.ec_list.iter().enumerate() {
        let count = counts.counts.get(idx).copied().unwrap_or(0);
        map.insert(ec.clone(), count);
    }
    map
}
