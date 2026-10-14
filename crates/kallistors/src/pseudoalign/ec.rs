use super::KmerEcIndex;

const INVALID_BASE_CODE: u8 = 4;
const BASE_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'g' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes[b't' as usize] = 3;
    codes
};

pub(crate) fn ec_has_offlist(ec: &[u32], onlist: &[bool]) -> bool {
    for &t in ec {
        let idx = t as usize;
        if idx >= onlist.len() || !onlist[idx] {
            return true;
        }
    }
    false
}

pub(crate) fn filter_shades(ec: &[u32], shade_sequences: &[bool], out: &mut Vec<u32>) {
    out.clear();
    for &tr in ec {
        let idx = tr as usize;
        if idx < shade_sequences.len() && shade_sequences[idx] {
            continue;
        }
        out.push(tr);
    }
}

pub(crate) fn filter_shades_in_place(ec: &mut Vec<u32>, shade_sequences: &[bool]) {
    ec.retain(|&tr| {
        let idx = tr as usize;
        idx >= shade_sequences.len() || !shade_sequences[idx]
    });
}

pub(crate) fn merge_sorted_unique_vec(base: &mut Vec<u32>, extra: &[u32]) {
    if extra.is_empty() {
        return;
    }
    if base.is_empty() {
        base.extend_from_slice(extra);
        return;
    }
    let mut merged = Vec::with_capacity(base.len() + extra.len());
    let mut i = 0;
    let mut j = 0;
    while i < base.len() && j < extra.len() {
        if base[i] == extra[j] {
            merged.push(base[i]);
            i += 1;
            j += 1;
        } else if base[i] < extra[j] {
            merged.push(base[i]);
            i += 1;
        } else {
            merged.push(extra[j]);
            j += 1;
        }
    }
    if i < base.len() {
        merged.extend_from_slice(&base[i..]);
    }
    if j < extra.len() {
        merged.extend_from_slice(&extra[j..]);
    }
    merged.dedup();
    *base = merged;
}

pub(crate) fn intersect_sorted(a: &[u32], b: &[u32], out: &mut Vec<u32>) {
    out.clear();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            out.push(a[i]);
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
}

pub(crate) fn merge_sorted_unique(a: &mut Vec<u32>, b: &[u32]) {
    if a.is_empty() {
        a.extend_from_slice(b);
        return;
    }
    let mut merged = Vec::with_capacity(a.len() + b.len());
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            merged.push(a[i]);
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            merged.push(a[i]);
            i += 1;
        } else {
            merged.push(b[j]);
            j += 1;
        }
    }
    if i < a.len() {
        merged.extend_from_slice(&a[i..]);
    }
    if j < b.len() {
        merged.extend_from_slice(&b[j..]);
    }
    merged.dedup();
    *a = merged;
}

pub(crate) fn encode_kmer_pair(seq: &[u8]) -> Option<(u64, u64)> {
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut rev_shift = 0u32;
    for &b in seq {
        let v = BASE_CODES[b as usize];
        if v == INVALID_BASE_CODE {
            return None;
        }
        let v = u64::from(v);
        fwd = (fwd << 2) | v;
        rev |= (3 - v) << rev_shift;
        rev_shift += 2;
    }
    Some((fwd, rev))
}

pub(crate) fn lookup_ec(index: &KmerEcIndex, code: u64) -> Option<&[u32]> {
    let idx = index.mphf.try_hash(&code)? as usize;
    if index.keys_by_index.get(idx) == Some(&code) {
        return Some(&index.ecs[idx]);
    }
    None
}

/// Reuse exact words only when the scanner advances one base; jumps re-encode.
pub(crate) struct KmerPairCursor {
    k: usize,
    mask: u64,
    previous: Option<(usize, u64, u64)>,
}

impl KmerPairCursor {
    pub(crate) fn new(k: usize) -> Self {
        Self {
            k,
            mask: if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            },
            previous: None,
        }
    }

    #[inline]
    pub(crate) fn encode(&mut self, kmer: &[u8], pos: usize) -> Option<(u64, u64)> {
        let codes = if self.k > 0
            && kmer.len() == self.k
            && let Some((previous_pos, forward, reverse)) = self.previous
            && previous_pos.checked_add(1) == Some(pos)
        {
            let code = BASE_CODES[kmer[self.k - 1] as usize];
            if code == INVALID_BASE_CODE {
                None
            } else {
                let code = u64::from(code);
                Some((
                    ((forward << 2) | code) & self.mask,
                    (reverse >> 2) | ((3 - code) << (2 * (self.k - 1))),
                ))
            }
        } else {
            encode_kmer_pair(kmer)
        };
        self.previous = codes.map(|(forward, reverse)| (pos, forward, reverse));
        codes
    }
}

#[cfg(test)]
mod rolling_read_tests {
    use super::*;

    #[test]
    fn adjacent_words_preserve_byte_encoding_with_jumps_and_invalid_bases() {
        let seq: Vec<u8> = (0..180)
            .map(|i| b"ACGTacgtN"[(i * i + 3 * i + i / 11) % 9])
            .collect();
        let valid: Vec<u8> = (0..180)
            .map(|i| b"ACGTacgt"[(i * i + 3 * i + i / 11) % 8])
            .collect();
        for seq in [&seq, &valid] {
            for k in 1..=32 {
                let mut cursor = KmerPairCursor::new(k);
                let positions = (0..=seq.len() - k).chain([19, 20, 21, 70, 1, 2, 3]);
                for pos in positions {
                    let kmer = &seq[pos..pos + k];
                    assert_eq!(
                        cursor.encode(kmer, pos),
                        encode_kmer_pair(kmer),
                        "k={k} pos={pos}"
                    );
                }
            }
        }
    }
}
