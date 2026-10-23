use super::*;

fn reverse_complement_word_slow(word: u64, k: usize) -> u64 {
    let mut out = 0u64;
    for i in 0..k {
        let shift = 62 - ((i & 0x1f) << 1);
        let code = (word >> shift) & 3;
        let out_shift = 62 - (((k - 1 - i) & 0x1f) << 1);
        out |= (3 - code) << out_shift;
    }
    out
}

fn selected_minimizers_naive(seq: &[u8], k: usize, g: usize) -> Vec<(u32, [u8; 8])> {
    let Some(bounds) = bifrost_minimizer_bounds(k, g) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for kmer_start in 0..=seq.len() - k {
        let mut best_hash = u64::MAX;
        let mut selected = Vec::new();
        for rel_pos in bounds {
            let pos = kmer_start + rel_pos;
            let min_seq = &seq[pos..pos + g];
            let hash = minimizer_rep_hash(min_seq).expect("hash minimizer");
            let bytes = encode_minimizer_rep(min_seq).expect("encode minimizer");
            match hash.cmp(&best_hash) {
                std::cmp::Ordering::Less => {
                    best_hash = hash;
                    selected.clear();
                    selected.push((pos as u32, bytes));
                }
                std::cmp::Ordering::Equal => selected.push((pos as u32, bytes)),
                std::cmp::Ordering::Greater => {}
            }
        }
        out.extend(selected);
    }
    out
}

#[path = "tests/collection.rs"]
mod collection;
#[path = "tests/compaction.rs"]
mod compaction;
#[path = "tests/lookup_edges.rs"]
mod lookup_edges;
#[path = "tests/minimizers.rs"]
mod minimizers;
#[path = "tests/payloads.rs"]
mod payloads;
