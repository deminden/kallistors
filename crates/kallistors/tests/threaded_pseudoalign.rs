use std::collections::HashMap;
use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index, pseudoalign_paired_bifrost_with_options,
    pseudoalign_paired_bifrost_with_options_threaded, pseudoalign_single_end_bifrost_with_options,
    pseudoalign_single_end_bifrost_with_options_threaded,
};

mod helpers;

#[test]
fn threaded_single_end_matches_single_thread() {
    let dataset = match helpers::build_synthetic_index() {
        Some(dataset) => dataset,
        None => return,
    };
    let index = build_bifrost_index(&dataset.index_path).expect("bifrost index");
    let reads = synthesize_reads(&dataset.transcripts, 50, 12);
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

    assert_eq!(single.reads_processed, threaded.reads_processed);
    assert_eq!(single.reads_aligned, threaded.reads_aligned);
    assert_eq!(to_ec_map(&single), to_ec_map(&threaded));
}

#[test]
fn threaded_paired_matches_single_thread() {
    let dataset = match helpers::build_synthetic_index() {
        Some(dataset) => dataset,
        None => return,
    };
    let index = build_bifrost_index(&dataset.index_path).expect("bifrost index");
    let pairs = synthesize_paired_reads(&dataset.transcripts, 50, 90, 10);
    let (fq1, fq2) = to_fastq_pairs(&pairs);

    let mut r1 = FastqReader::new(Cursor::new(fq1.clone()));
    let mut r2 = FastqReader::new(Cursor::new(fq2.clone()));
    let single = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut r1,
        &mut r2,
        Strand::Unstranded,
        PseudoalignOptions::default(),
    )
    .expect("single-thread");

    let mut r1 = FastqReader::new(Cursor::new(fq1));
    let mut r2 = FastqReader::new(Cursor::new(fq2));
    let threaded = pseudoalign_paired_bifrost_with_options_threaded(
        &index,
        &mut r1,
        &mut r2,
        Strand::Unstranded,
        PseudoalignOptions::default(),
        4,
    )
    .expect("threaded");

    assert_eq!(single.reads_processed, threaded.reads_processed);
    assert_eq!(single.reads_aligned, threaded.reads_aligned);
    assert_eq!(to_ec_map(&single), to_ec_map(&threaded));
}

fn synthesize_reads(transcripts: &[Vec<u8>], read_len: usize, count: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(count);
    let mut idx = 0usize;
    let mut attempts = 0usize;
    while out.len() < count && attempts < count * 10 {
        let seq = &transcripts[idx % transcripts.len()];
        idx += 1;
        attempts += 1;
        if seq.len() < read_len {
            continue;
        }
        let start = (idx * 7) % (seq.len() - read_len + 1);
        out.push(seq[start..start + read_len].to_vec());
    }
    out
}

fn to_fastq(reads: &[Vec<u8>]) -> Vec<u8> {
    let mut fq = Vec::new();
    for (i, read) in reads.iter().enumerate() {
        fq.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        fq.extend_from_slice(read);
        fq.extend_from_slice(b"\n+\n");
        fq.extend_from_slice(vec![b'I'; read.len()].as_slice());
        fq.extend_from_slice(b"\n");
    }
    fq
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

fn to_ec_map(counts: &kallistors::pseudoalign::EcCounts) -> HashMap<Vec<u32>, u32> {
    let mut map = HashMap::new();
    for (idx, ec) in counts.ec_list.iter().enumerate() {
        let count = counts.counts.get(idx).copied().unwrap_or(0);
        map.insert(ec.clone(), count);
    }
    map
}
