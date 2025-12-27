use std::collections::HashMap;
use std::io::Cursor;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    EcCounts, build_bifrost_index, build_kmer_ec_index, pseudoalign_single_end,
    pseudoalign_single_end_bifrost,
};

mod helpers;

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
    let dataset = match helpers::build_synthetic_index() {
        Some(dataset) => dataset,
        None => return,
    };

    let reads = synthesize_reads(&dataset.transcripts, 60, 24);
    let fastq = to_fastq(&reads);

    let naive_index = build_kmer_ec_index(&dataset.index_path).expect("naive index");
    let bifrost_index = build_bifrost_index(&dataset.index_path).expect("bifrost index");

    let mut reader = FastqReader::new(Cursor::new(fastq.clone()));
    let naive_counts = pseudoalign_single_end(&naive_index, &mut reader).expect("naive counts");

    let mut reader = FastqReader::new(Cursor::new(fastq));
    let bifrost_counts =
        pseudoalign_single_end_bifrost(&bifrost_index, &mut reader).expect("bifrost counts");

    assert_eq!(naive_counts.reads_processed, bifrost_counts.reads_processed);
    assert_eq!(naive_counts.reads_aligned, bifrost_counts.reads_aligned);
    assert_eq!(normalize_ec(&naive_counts), normalize_ec(&bifrost_counts));
}

fn synthesize_reads(transcripts: &[Vec<u8>], read_len: usize, count: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(count);
    let mut idx = 0usize;
    for _ in 0..count {
        let seq = &transcripts[idx % transcripts.len()];
        idx += 1;
        if seq.len() <= read_len {
            out.push(seq.clone());
            continue;
        }
        let start = (idx * 11) % (seq.len() - read_len);
        out.push(seq[start..start + read_len].to_vec());
    }
    out
}

fn to_fastq(reads: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, seq) in reads.iter().enumerate() {
        out.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        out.extend_from_slice(seq);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(vec![b'I'; seq.len()].as_slice());
        out.extend_from_slice(b"\n");
    }
    out
}
