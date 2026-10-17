use core::range::RangeInclusive;

use crate::index::bifrost::{encode_minimizer_rep, wyhash};

const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];
// Hash groups use the same A/C/T/G REP codes as the byte oracle, but the input
// is the validated A/C/G/T two-bit word. XOR/rotation composition is exact.
const REP_HASH_GROUPS: [[u64; 256]; 5] = {
    let mut groups = [[0u64; 256]; 5];
    let mut len = 1;
    while len <= 4 {
        let mut word = 0;
        while word < 256 {
            let mut hash = 0u64;
            let mut i = 0;
            while i < len {
                let base = (word >> (2 * (len - i - 1))) & 3;
                let rep_code = base ^ (base >> 1);
                hash = hash.rotate_left(1) ^ REP_HASH_VALS[rep_code];
                i += 1;
            }
            groups[len][word] = hash;
            word += 1;
        }
        len += 1;
    }
    groups
};
// Grouped XOR/rotation composition for the common 23-base initializer.
static REP_GROUP_SIX: [u64; 4096] = {
    let mut groups = [0; 4096];
    let mut word = 0;
    while word < groups.len() {
        groups[word] =
            REP_HASH_GROUPS[4][word & 255] ^ REP_HASH_GROUPS[2][word >> 8].rotate_left(4);
        word += 1;
    }
    groups
};
static REP_GROUP_FIVE: [u64; 1024] = {
    let mut groups = [0; 1024];
    let mut word = 0;
    while word < groups.len() {
        groups[word] =
            REP_HASH_GROUPS[4][word & 255] ^ REP_HASH_GROUPS[1][word >> 8].rotate_left(4);
        word += 1;
    }
    groups
};

#[inline]
fn rep_grouped_23(word: u64) -> u64 {
    REP_GROUP_SIX[(word & 4095) as usize]
        ^ REP_GROUP_SIX[((word >> 12) & 4095) as usize].rotate_left(6)
        ^ REP_GROUP_SIX[((word >> 24) & 4095) as usize].rotate_left(12)
        ^ REP_GROUP_FIVE[((word >> 36) & 1023) as usize].rotate_left(18)
}

const REP_FORWARD: [u64; 4] = [
    REP_HASH_VALS[0],
    REP_HASH_VALS[1],
    REP_HASH_VALS[3],
    REP_HASH_VALS[2],
];
const REP_COMPLEMENT: [u64; 4] = [
    REP_HASH_VALS[2],
    REP_HASH_VALS[3],
    REP_HASH_VALS[1],
    REP_HASH_VALS[0],
];
const REP_REMOVE_FORWARD: [[u64; 4]; 33] = {
    let mut table = [[0; 4]; 33];
    let mut g = 0;
    while g <= 32 {
        let mut base = 0;
        while base < 4 {
            table[g][base] = REP_FORWARD[base].rotate_left(g as u32);
            base += 1;
        }
        g += 1;
    }
    table
};
const REP_ADD_REVERSE: [[u64; 4]; 33] = {
    let mut table = [[0; 4]; 33];
    let mut g = 1;
    while g <= 32 {
        let mut base = 0;
        while base < 4 {
            table[g][base] = REP_COMPLEMENT[base].rotate_left((g - 1) as u32);
            base += 1;
        }
        g += 1;
    }
    table
};
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

/// Reuse the exact forward/reverse words already validated by the read encoder.
pub(super) fn minimizers_for_encoded_kmer_into(
    seq: &[u8],
    g: usize,
    forward: u64,
    reverse: u64,
    out: &mut Vec<([u8; 8], usize)>,
) -> bool {
    if seq.len() > 32 {
        return minimizers_for_kmer_into(seq, g, out);
    }
    if seq.len() == 31 && g == 23 {
        return minimizers_31_23(forward, reverse, out);
    }
    out.clear();
    let Some(bounds) = strict_bounds(seq.len(), g) else {
        return false;
    };
    let mask = if g == 0 { 0 } else { u64::MAX >> (64 - 2 * g) };
    let mut best_hash = None;
    let mut rolling = RepHashWindow::from_encoded(
        (forward >> (2 * (seq.len() - g - bounds.start))) & mask,
        (reverse >> (2 * bounds.start)) & mask,
        g,
    );
    let first = bounds.start;
    for pos in bounds {
        if pos > first {
            let old = ((forward >> (2 * (seq.len() - pos))) & 3) as usize;
            let new = ((forward >> (2 * (seq.len() - pos - g))) & 3) as usize;
            rolling = rolling.advance_encoded(old, new, g);
        }
        let hash = rolling.finish();
        if best_hash.is_some_and(|best| hash > best) {
            continue;
        }
        let rep = if g == 0 {
            0
        } else {
            let f = (forward >> (2 * (seq.len() - g - pos))) & mask;
            let r = (reverse >> (2 * pos)) & mask;
            f.min(r) << (64 - 2 * g)
        };
        if best_hash != Some(hash) {
            out.clear();
            best_hash = Some(hash);
        }
        out.push((rep.to_le_bytes(), pos));
    }
    !out.is_empty()
}

fn minimizers_31_23(forward: u64, reverse: u64, out: &mut Vec<([u8; 8], usize)>) -> bool {
    const MASK: u64 = (1u64 << 46) - 1;
    let mut rolling =
        RepHashWindow::from_encoded((forward >> 14) & MASK, (reverse >> 2) & MASK, 23);
    let mut hashes = [0u64; 6];
    for (offset, hash) in hashes.iter_mut().enumerate() {
        let pos = offset + 1;
        if offset != 0 {
            let old = ((forward >> (2 * (31 - pos))) & 3) as usize;
            let new = ((forward >> (2 * (31 - pos - 23))) & 3) as usize;
            rolling = rolling.advance_encoded(old, new, 23);
        }
        *hash = rolling.finish();
    }
    let best = hashes.iter().copied().min().expect("six windows");
    out.clear();
    for (offset, hash) in hashes.into_iter().enumerate() {
        if hash == best {
            let pos = offset + 1;
            let f = (forward >> (2 * (8 - pos))) & MASK;
            let r = (reverse >> (2 * pos)) & MASK;
            out.push(((f.min(r) << 18).to_le_bytes(), pos));
        }
    }
    true
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

#[derive(Clone, Copy)]
struct RepHashWindow {
    forward: u64,
    reverse: u64,
}

impl RepHashWindow {
    #[inline]
    fn from_encoded(mut fword: u64, mut rword: u64, g: usize) -> Self {
        if g == 23 {
            return Self {
                forward: rep_grouped_23(fword),
                reverse: rep_grouped_23(rword),
            };
        }
        let mut forward = 0u64;
        let mut reverse = 0u64;
        for chunk in 0..g / 4 {
            forward ^= REP_HASH_GROUPS[4][(fword & 255) as usize].rotate_left((chunk * 4) as u32);
            reverse ^= REP_HASH_GROUPS[4][(rword & 255) as usize].rotate_left((chunk * 4) as u32);
            fword >>= 8;
            rword >>= 8;
        }
        let tail = g % 4;
        if tail != 0 {
            let rotation = (g - tail) as u32;
            forward ^= REP_HASH_GROUPS[tail][fword as usize].rotate_left(rotation);
            reverse ^= REP_HASH_GROUPS[tail][rword as usize].rotate_left(rotation);
        }
        Self { forward, reverse }
    }

    #[cfg(test)]
    fn new(seq: &[u8]) -> Option<Self> {
        let mut forward = 0u64;
        let mut reverse = 0u64;
        for (i, &base) in seq.iter().enumerate() {
            let f = HASH_CODES[base as usize];
            let r = HASH_COMP_CODES[seq[seq.len() - i - 1] as usize];
            if f == INVALID_BASE_CODE || r == INVALID_BASE_CODE {
                return None;
            }
            forward = forward.rotate_left(1) ^ REP_HASH_VALS[f as usize];
            reverse = reverse.rotate_left(1) ^ REP_HASH_VALS[r as usize];
        }
        Some(Self { forward, reverse })
    }

    #[inline]
    fn advance_encoded(self, old: usize, new: usize, g: usize) -> Self {
        if g == 0 {
            return self;
        }
        Self {
            forward: self.forward.rotate_left(1) ^ REP_REMOVE_FORWARD[g][old] ^ REP_FORWARD[new],
            reverse: self.reverse.rotate_right(1)
                ^ REP_COMPLEMENT[old].rotate_right(1)
                ^ REP_ADD_REVERSE[g][new],
        }
    }

    #[cfg(test)]
    fn advance(self, old: u8, new: u8, g: usize) -> Option<Self> {
        if g == 0 {
            return Some(self);
        }
        let old_f = HASH_CODES[old as usize];
        let new_f = HASH_CODES[new as usize];
        let old_r = HASH_COMP_CODES[old as usize];
        let new_r = HASH_COMP_CODES[new as usize];
        if [old_f, new_f, old_r, new_r].contains(&INVALID_BASE_CODE) {
            return None;
        }
        Some(Self {
            forward: self.forward.rotate_left(1)
                ^ REP_HASH_VALS[old_f as usize].rotate_left(g as u32)
                ^ REP_HASH_VALS[new_f as usize],
            reverse: self.reverse.rotate_right(1)
                ^ REP_HASH_VALS[old_r as usize].rotate_right(1)
                ^ REP_HASH_VALS[new_r as usize].rotate_left((g - 1) as u32),
        })
    }

    fn finish(self) -> u64 {
        let mut buf = [0u8; 16];
        buf[..8].copy_from_slice(&self.forward.min(self.reverse).to_le_bytes());
        buf[8..].copy_from_slice(&self.forward.max(self.reverse).to_le_bytes());
        wyhash(&buf, 0)
    }
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

#[cfg(test)]
mod encoded_tests {
    use super::*;

    #[test]
    fn common_window_specialization_keeps_every_homopolymer_tie() {
        let seq = [b'A'; 31];
        let (forward, reverse) = crate::pseudoalign::ec::encode_kmer_pair(&seq).unwrap();
        let mut got = Vec::new();
        assert!(minimizers_for_encoded_kmer_into(
            &seq, 23, forward, reverse, &mut got
        ));
        assert_eq!(got, minimizers_for_kmer(&seq, 23).unwrap());
        assert_eq!(
            got.iter().map(|entry| entry.1).collect::<Vec<_>>(),
            vec![1, 2, 3, 4, 5, 6]
        );
    }

    #[test]
    fn rolling_rep_hash_equals_byte_oracle_for_every_window() {
        let seq: Vec<u8> = (0..161)
            .map(|i| b"ACGTacgt"[(i * i + 3 * i + i / 7) % 8])
            .collect();
        for g in 0..=32 {
            let mut window = RepHashWindow::new(&seq[..g]).expect("valid bases");
            for pos in 0..=seq.len() - g {
                if pos > 0 {
                    window = window
                        .advance(seq[pos - 1], seq[pos + g - 1], g)
                        .expect("valid slide");
                }
                assert_eq!(
                    Some(window.finish()),
                    rep_hash(&seq[pos..pos + g]),
                    "g={g} pos={pos}"
                );
            }
        }
        assert!(RepHashWindow::new(b"ACNT").is_none());
        assert!(
            RepHashWindow::new(b"ACGT")
                .unwrap()
                .advance(b'A', b'N', 4)
                .is_none()
        );
    }

    #[test]
    fn encoded_candidates_preserve_hashes_positions_and_ties() {
        let mut state = 7u64;
        for k in 2..=32 {
            for case in 0..80 {
                let seq: Vec<u8> = (0..k)
                    .map(|_| {
                        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                        if case < 8 {
                            b"ACGTacgt"[case]
                        } else {
                            b"ACGTacgt"[(state >> 32) as usize % 8]
                        }
                    })
                    .collect();
                let mut f = 0u64;
                let mut r = 0u64;
                for (i, b) in seq.iter().enumerate() {
                    let code = match b.to_ascii_uppercase() {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => unreachable!(),
                    };
                    f = (f << 2) | code;
                    r |= (3 - code) << (2 * i);
                }
                for g in 0..=k {
                    let mut old = Vec::new();
                    let mut new = Vec::new();
                    assert_eq!(
                        minimizers_for_kmer_into(&seq, g, &mut old),
                        minimizers_for_encoded_kmer_into(&seq, g, f, r, &mut new)
                    );
                    assert_eq!(old, new, "k={k}, g={g}, case={case}");
                }
            }
        }
    }
}
