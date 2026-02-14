use crate::index::bifrost::{encode_minimizer_rep, wyhash};

const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];

#[inline]
fn bifrost_neighbor_bounds(k: usize, g: usize) -> Option<(usize, usize)> {
    // Bifrost minHashKmer(..., neighbor_hash=true) scans [shift, k-g-shift] inclusive.
    // With shift=1 this is [1, k-g-1].
    let shift = 1usize;
    if k < g + shift + 1 {
        return None;
    }
    let end = k.checked_sub(g + shift)?;
    if end < shift {
        None
    } else {
        Some((shift, end))
    }
}

#[inline]
fn strict_bounds(k: usize, g: usize) -> Option<(usize, usize)> {
    if k < g + 2 {
        return None;
    }
    let start = 1usize;
    let mut end = k.saturating_sub(g + 2);
    if end < start {
        end = start;
    }
    Some((start, end))
}

pub(super) fn minimizers_for_kmer(seq: &[u8], g: usize) -> Option<Vec<([u8; 8], usize)>> {
    let (start, end) = strict_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best: Vec<([u8; 8], usize)> = Vec::new();
    for pos in start..=end {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        match best_hash {
            None => {
                best_hash = Some(h);
                best.push((bytes, pos));
            }
            Some(min) if h < min => {
                best_hash = Some(h);
                best.clear();
                best.push((bytes, pos));
            }
            Some(min) if h == min => {
                best.push((bytes, pos));
            }
            _ => {}
        }
    }
    if best.is_empty() { None } else { Some(best) }
}

pub(super) fn minimizers_ranked_for_kmer(
    seq: &[u8],
    g: usize,
    max_hashes: usize,
) -> Option<Vec<([u8; 8], usize)>> {
    if max_hashes == 0 {
        return None;
    }
    let (start, end) = strict_bounds(seq.len(), g)?;
    let mut all: Vec<(u64, [u8; 8], usize)> = Vec::new();
    for pos in start..=end {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        all.push((h, bytes, pos));
    }
    if all.is_empty() {
        return None;
    }
    all.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out = Vec::new();
    let mut seen = 0usize;
    let mut current_hash: Option<u64> = None;
    for (h, bytes, pos) in all {
        if current_hash.map(|v| v != h).unwrap_or(true) {
            seen += 1;
            if seen > max_hashes {
                break;
            }
            current_hash = Some(h);
        }
        out.push((bytes, pos));
    }
    if out.is_empty() { None } else { Some(out) }
}

type MinimizerCandidates = Vec<([u8; 8], usize)>;
type MinhashCandidates = (MinimizerCandidates, Option<([u8; 8], usize)>);

pub(super) fn minhash_primary_for_kmer(seq: &[u8], g: usize) -> Option<(u64, MinimizerCandidates)> {
    let (start, end) = bifrost_neighbor_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best: Vec<([u8; 8], usize)> = Vec::new();

    for pos in start..=end {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        match best_hash {
            None => {
                best_hash = Some(h);
                best.push((bytes, pos));
            }
            Some(min) if h < min => {
                best_hash = Some(h);
                best.clear();
                best.push((bytes, pos));
            }
            Some(min) if h == min => {
                best.push((bytes, pos));
            }
            _ => {}
        }
    }

    best_hash.map(|h| (h, best))
}

pub(super) fn minhash_next_after_hash(
    seq: &[u8],
    g: usize,
    min_hash: u64,
) -> Option<(u64, ([u8; 8], usize))> {
    let (start, end) = bifrost_neighbor_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best_rep: Option<u64> = None;
    let mut best_bytes: [u8; 8] = [0; 8];
    let mut best_pos: usize = 0;

    for pos in start..=end {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        if h <= min_hash {
            continue;
        }
        let bytes = encode_minimizer_rep(slice)?;
        let rep = u64::from_le_bytes(bytes);
        match best_hash {
            None => {
                best_hash = Some(h);
                best_rep = Some(rep);
                best_bytes = bytes;
                best_pos = pos;
            }
            Some(curr_h) if h < curr_h => {
                best_hash = Some(h);
                best_rep = Some(rep);
                best_bytes = bytes;
                best_pos = pos;
            }
            Some(curr_h) if h == curr_h => {
                // Match Bifrost minHashKmer::compute_min(min_v): on equal hash,
                // choose lexicographically smaller canonical minimizer rep.
                // If rep is identical, keep the first position encountered.
                if let Some(curr_rep) = best_rep
                    && rep < curr_rep
                {
                    best_rep = Some(rep);
                    best_bytes = bytes;
                    best_pos = pos;
                }
            }
            _ => {}
        }
    }

    best_hash.map(|h| (h, (best_bytes, best_pos)))
}

pub(super) fn minhash_candidates_for_kmer(seq: &[u8], g: usize) -> Option<MinhashCandidates> {
    let (min_hash, min_positions) = minhash_primary_for_kmer(seq, g)?;
    let next_min = minhash_next_after_hash(seq, g, min_hash).map(|(_, v)| v);
    Some((min_positions, next_min))
}

pub(super) fn minimizer_for_kmer_strict(seq: &[u8], g: usize) -> Option<([u8; 8], usize)> {
    let (start, end) = strict_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best: Option<([u8; 8], usize)> = None;
    for pos in start..=end {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        match best_hash {
            None => {
                best_hash = Some(h);
                best = Some((bytes, pos));
            }
            Some(min) if h < min => {
                best_hash = Some(h);
                best = Some((bytes, pos));
            }
            Some(min) if h == min => {
                if let Some((_, best_pos)) = best
                    && pos < best_pos
                {
                    best = Some((bytes, pos));
                }
            }
            _ => {}
        }
    }
    best
}

pub(super) fn fill_revcomp(seq: &[u8], out: &mut Vec<u8>) -> bool {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => return false,
        };
        out.push(comp);
    }
    true
}

fn rep_hash(seq: &[u8]) -> Option<u64> {
    let mut h = 0u64;
    let mut ht = 0u64;
    let len = seq.len();
    for i in 0..len {
        let b = seq[i];
        if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            return None;
        }
        h = h.rotate_left(1);
        ht = ht.rotate_left(1);
        let f = ((b & 6) >> 1) as usize;
        let rb = seq[len - 1 - i];
        let r = (((rb ^ 4) & 6) >> 1) as usize;
        h ^= REP_HASH_VALS[f];
        ht ^= REP_HASH_VALS[r];
    }
    let mut hashes = [h, ht];
    if hashes[1] < hashes[0] {
        hashes.swap(0, 1);
    }
    let mut buf = [0u8; 16];
    buf[..8].copy_from_slice(&hashes[0].to_le_bytes());
    buf[8..].copy_from_slice(&hashes[1].to_le_bytes());
    Some(wyhash(&buf, 0))
}
