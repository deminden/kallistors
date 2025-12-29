use crate::index::bifrost::{encode_minimizer_rep, wyhash};

const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];

pub(super) fn minimizers_for_kmer(seq: &[u8], g: usize) -> Option<Vec<([u8; 8], usize)>> {
    if seq.len() < g + 2 {
        return None;
    }
    let mut best_hash: Option<u64> = None;
    let mut best: Vec<([u8; 8], usize)> = Vec::new();
    let start = 1usize;
    let mut end = seq.len().saturating_sub(g + 2);
    if end < start {
        end = start;
    }
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
    if max_hashes == 0 || seq.len() < g + 2 {
        return None;
    }
    let start = 1usize;
    let mut end = seq.len().saturating_sub(g + 2);
    if end < start {
        end = start;
    }
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

pub(super) fn minimizers_sorted_for_kmer(
    seq: &[u8],
    g: usize,
) -> Option<Vec<(u64, [u8; 8], usize)>> {
    if seq.len() < g + 2 {
        return None;
    }
    let start = 1usize;
    let mut end = seq.len().saturating_sub(g + 2);
    if end < start {
        end = start;
    }
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
    all.sort_by(|a, b| a.0.cmp(&b.0).then(a.2.cmp(&b.2)));
    Some(all)
}

type MinimizerCandidates = Vec<([u8; 8], usize)>;
type MinhashCandidates = (MinimizerCandidates, Option<([u8; 8], usize)>);

pub(super) fn minhash_candidates_for_kmer(seq: &[u8], g: usize) -> Option<MinhashCandidates> {
    if seq.len() < g + 2 {
        return None;
    }
    let k = seq.len();
    let shift = 1usize;
    let max_pos = k.saturating_sub(g + shift);
    if max_pos < shift {
        return None;
    }
    let mut all: Vec<(u64, [u8; 8], usize)> = Vec::new();
    for pos in shift..=max_pos {
        let slice = &seq[pos..pos + g];
        let h = rep_hash(slice)?;
        let bytes = encode_minimizer_rep(slice)?;
        all.push((h, bytes, pos));
    }
    if all.is_empty() {
        return None;
    }
    let min_hash = all.iter().map(|v| v.0).min()?;
    let mut min_positions: Vec<([u8; 8], usize)> = all
        .iter()
        .filter(|v| v.0 == min_hash)
        .map(|v| (v.1, v.2))
        .collect();
    min_positions.sort_by(|a, b| a.1.cmp(&b.1));

    let mut next_min: Option<(u64, [u8; 8], usize)> = None;
    for (h, bytes, pos) in all.into_iter() {
        if h <= min_hash {
            continue;
        }
        match next_min {
            None => next_min = Some((h, bytes, pos)),
            Some((best_h, best_bytes, best_pos)) => {
                if h < best_h
                    || (h == best_h
                        && (bytes < best_bytes || (bytes == best_bytes && pos < best_pos)))
                {
                    next_min = Some((h, bytes, pos));
                }
            }
        }
    }
    let next_min = next_min.map(|v| (v.1, v.2));
    Some((min_positions, next_min))
}

pub(super) fn minimizer_for_kmer_strict(seq: &[u8], g: usize) -> Option<([u8; 8], usize)> {
    if seq.len() < g + 2 {
        return None;
    }
    let start = 1usize;
    let mut end = seq.len().saturating_sub(g + 2);
    if end < start {
        end = start;
    }
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
