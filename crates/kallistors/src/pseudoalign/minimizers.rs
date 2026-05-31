use core::range::RangeInclusive;

use crate::index::bifrost::{encode_minimizer_rep, wyhash};

const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];
const INVALID_BASE_CODE: u8 = 4;
const HASH_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'T' as usize] = 2;
    codes[b't' as usize] = 2;
    codes[b'G' as usize] = 3;
    codes[b'g' as usize] = 3;
    codes
};
const HASH_COMP_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'T' as usize] = 0;
    codes[b't' as usize] = 0;
    codes[b'G' as usize] = 1;
    codes[b'g' as usize] = 1;
    codes[b'A' as usize] = 2;
    codes[b'a' as usize] = 2;
    codes[b'C' as usize] = 3;
    codes[b'c' as usize] = 3;
    codes
};
const COMP_BASES: [u8; 256] = {
    let mut bases = [0; 256];
    bases[b'A' as usize] = b'T';
    bases[b'a' as usize] = b'T';
    bases[b'C' as usize] = b'G';
    bases[b'c' as usize] = b'G';
    bases[b'G' as usize] = b'C';
    bases[b'g' as usize] = b'C';
    bases[b'T' as usize] = b'A';
    bases[b't' as usize] = b'A';
    bases
};

#[inline]
fn bifrost_neighbor_bounds(k: usize, g: usize) -> Option<RangeInclusive<usize>> {
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
        Some(RangeInclusive {
            start: shift,
            last: end,
        })
    }
}

#[inline]
fn strict_bounds(k: usize, g: usize) -> Option<RangeInclusive<usize>> {
    if k < g + 2 {
        return None;
    }
    let start = 1usize;
    let mut end = k.saturating_sub(g + 2);
    if end < start {
        end = start;
    }
    Some(RangeInclusive { start, last: end })
}

pub(super) fn minimizers_for_kmer(seq: &[u8], g: usize) -> Option<Vec<([u8; 8], usize)>> {
    let bounds = strict_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best: Vec<([u8; 8], usize)> = Vec::new();
    for pos in bounds {
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

pub(super) fn minimizers_for_kmer_into(
    seq: &[u8],
    g: usize,
    out: &mut Vec<([u8; 8], usize)>,
) -> bool {
    let Some(bounds) = strict_bounds(seq.len(), g) else {
        out.clear();
        return false;
    };
    out.clear();
    let mut best_hash: Option<u64> = None;
    for pos in bounds {
        let slice = &seq[pos..pos + g];
        let Some(h) = rep_hash(slice) else {
            out.clear();
            return false;
        };
        let Some(bytes) = encode_minimizer_rep(slice) else {
            out.clear();
            return false;
        };
        match best_hash {
            None => {
                best_hash = Some(h);
                out.push((bytes, pos));
            }
            Some(min) if h < min => {
                best_hash = Some(h);
                out.clear();
                out.push((bytes, pos));
            }
            Some(min) if h == min => out.push((bytes, pos)),
            _ => {}
        }
    }
    !out.is_empty()
}

pub(super) fn minimizers_ranked_for_kmer(
    seq: &[u8],
    g: usize,
    max_hashes: usize,
) -> Option<Vec<([u8; 8], usize)>> {
    if max_hashes == 0 {
        return None;
    }
    let bounds = strict_bounds(seq.len(), g)?;
    let mut all: Vec<(u64, [u8; 8], usize)> = Vec::new();
    for pos in bounds {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        all.push((h, bytes, pos));
    }
    if all.is_empty() {
        return None;
    }
    all.sort_by_key(|a| a.0);
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
    let bounds = bifrost_neighbor_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best: Vec<([u8; 8], usize)> = Vec::new();

    for pos in bounds {
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
    let bounds = bifrost_neighbor_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best_rep: Option<u64> = None;
    let mut best_bytes: [u8; 8] = [0; 8];
    let mut best_pos: usize = 0;

    for pos in bounds {
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

pub(super) fn minimizer_tail_for_kmer(seq: &[u8], g: usize) -> Option<([u8; 8], usize)> {
    let bounds = bifrost_neighbor_bounds(seq.len(), g)?;
    let end = bounds.last;
    let slice = &seq[end..end + g];
    let bytes = encode_minimizer_rep(slice)?;
    Some((bytes, end))
}

pub(super) fn minimizer_for_kmer_strict(seq: &[u8], g: usize) -> Option<([u8; 8], usize)> {
    let bounds = strict_bounds(seq.len(), g)?;
    let mut best_hash: Option<u64> = None;
    let mut best: Option<([u8; 8], usize)> = None;
    for pos in bounds {
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
        let comp = COMP_BASES[b as usize];
        if comp == 0 {
            return false;
        }
        out.push(comp);
    }
    true
}

fn rep_hash(seq: &[u8]) -> Option<u64> {
    let mut h = 0u64;
    let mut ht = 0u64;
    let len = seq.len();
    for i in 0..len {
        let f = HASH_CODES[seq[i] as usize];
        if f == INVALID_BASE_CODE {
            return None;
        }
        h = h.rotate_left(1);
        ht = ht.rotate_left(1);
        let r = HASH_COMP_CODES[seq[len - 1 - i] as usize];
        if r == INVALID_BASE_CODE {
            return None;
        }
        h ^= REP_HASH_VALS[f as usize];
        ht ^= REP_HASH_VALS[r as usize];
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
