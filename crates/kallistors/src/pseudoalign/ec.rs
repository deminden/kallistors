use super::KmerEcIndex;

use crate::dna::{BASE_CODES, INVALID_BASE_CODE};

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
    if a == b {
        out.extend_from_slice(a);
        return;
    }
    if let [value] = a {
        if b.binary_search(value).is_ok() {
            out.push(*value);
        }
        return;
    }
    if let [value] = b {
        if a.binary_search(value).is_ok() {
            out.push(*value);
        }
        return;
    }
    let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if !small.is_empty() && large.len() / small.len() >= 16 {
        let mut cursor = 0;
        for &value in small {
            cursor += large[cursor..].partition_point(|&candidate| candidate < value);
            if large.get(cursor) == Some(&value) {
                out.push(value);
                cursor += 1;
            }
        }
        return;
    }
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

// Both inputs are sorted multisets. Writes never overtake the unread prefix.
#[inline(always)]
pub(crate) fn intersect_sorted_in_place(a: &mut Vec<u32>, b: &[u32]) {
    if a.is_empty() || b.is_empty() || a[0] > b[b.len() - 1] || b[0] > a[a.len() - 1] {
        a.clear();
        return;
    }
    if a.as_slice() == b {
        return;
    }
    if a.len() == 1 {
        if b.binary_search(&a[0]).is_err() {
            a.clear();
        }
        return;
    }
    if b.len() == 1 {
        if a.binary_search(&b[0]).is_ok() {
            a[0] = b[0];
            a.truncate(1);
        } else {
            a.clear();
        }
        return;
    }
    let mut written = 0;
    if a.len() <= b.len() / 16 {
        let mut cursor = 0;
        for i in 0..a.len() {
            let value = a[i];
            cursor += b[cursor..].partition_point(|&candidate| candidate < value);
            if b.get(cursor) == Some(&value) {
                a[written] = value;
                written += 1;
                cursor += 1;
            }
        }
    } else if b.len() <= a.len() / 16 {
        let mut cursor = 0;
        for &value in b {
            cursor += a[cursor..].partition_point(|&candidate| candidate < value);
            if a.get(cursor) == Some(&value) {
                a[written] = value;
                written += 1;
                cursor += 1;
            }
        }
    } else {
        let (mut i, mut j) = (0, 0);
        let common = a.len().min(b.len());
        if common >= 32 {
            // Four position-wise equal values are exactly the next four
            // multiset matches. Their output already occupies this prefix;
            // leave it untouched, then resume the original merge at the first
            // unequal block. Complete equality was handled above as before.
            while common - i >= 4 && a[i..i + 4] == b[i..i + 4] {
                i += 4;
            }
            j = i;
            written = i;
        }
        while i < a.len() && j < b.len() {
            if a[i] == b[j] {
                let value = a[i];
                // SAFETY: written=i at loop entry (including a preserved prefix). An equal pair increments
                // both; only advancing i can increase their gap. Therefore
                // written<=i<a.len() throughout this branch. No allocation
                // or truncation occurs in the loop, and value is read before
                // the mutable slot borrow, including the written==i case.
                unsafe {
                    *a.get_unchecked_mut(written) = value;
                }
                written += 1;
                i += 1;
                j += 1;
            } else if a[i] < b[j] {
                i += 1;
            } else {
                j += 1;
            }
        }
    }
    a.truncate(written);
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

pub(crate) use crate::dna::encode_kmer_pair;

#[cfg(test)]
fn encode_kmer_pair_byte_oracle(seq: &[u8]) -> Option<(u64, u64)> {
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

/// Match Bifrost's iterator advance: a requested probe skips windows containing
/// ambiguous bases and may reach the end instead of producing a lookup miss.
#[inline]
pub(crate) fn valid_kmer_start_at_or_after(
    seq: &[u8],
    mut start: usize,
    k: usize,
) -> Option<usize> {
    if k == 0 {
        return None;
    }
    loop {
        let window = seq.get(start..start.checked_add(k)?)?;
        let Some(invalid) = window
            .iter()
            .rposition(|&base| BASE_CODES[base as usize] == INVALID_BASE_CODE)
        else {
            return Some(start);
        };
        start += invalid + 1;
    }
}

#[inline]
pub(crate) fn encoded_kmer_at_or_after(
    seq: &[u8],
    start: usize,
    k: usize,
) -> Option<(usize, (u64, u64))> {
    let window = seq.get(start..start.checked_add(k)?)?;
    if let Some(codes) = encode_kmer_pair(window) {
        return Some((start, codes));
    }
    let valid = valid_kmer_start_at_or_after(seq, start, k)?;
    encode_kmer_pair(&seq[valid..valid + k]).map(|codes| (valid, codes))
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
    fn jump_probes_skip_ambiguous_windows_like_the_upstream_iterator() {
        let mut sequences = vec![b"ACGTacgtACGT".to_vec(), b"NNNNNN".to_vec()];
        let mut read = b"ACGT".repeat(32);
        read[58] = b'N';
        read[95] = b'n';
        sequences.push(read);
        for sequence in sequences {
            for k in [1, 3, 31, 32, 47] {
                for start in 0..=sequence.len() + 1 {
                    let expected = sequence.get(start..).and_then(|suffix| {
                        suffix
                            .windows(k)
                            .position(|window| {
                                window.iter().all(|base| {
                                    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                                })
                            })
                            .map(|offset| start + offset)
                    });
                    assert_eq!(valid_kmer_start_at_or_after(&sequence, start, k), expected);
                    if k <= 32 {
                        let encoded = encoded_kmer_at_or_after(&sequence, start, k);
                        assert_eq!(encoded.map(|(position, _)| position), expected);
                        if let Some((position, codes)) = encoded {
                            assert_eq!(
                                Some(codes),
                                encode_kmer_pair_byte_oracle(&sequence[position..position + k])
                            );
                        }
                    }
                }
            }
        }
        assert_eq!(valid_kmer_start_at_or_after(b"ACGT", 0, 0), None);
        assert_eq!(valid_kmer_start_at_or_after(b"ACGT", usize::MAX, 31), None);
    }

    #[test]
    fn sorted_intersection_preserves_multisets_and_reused_output() {
        fn check(a: &[u32], b: &[u32], out: &mut Vec<u32>) {
            // Independent multiplicity oracle; shares no cursor, write
            // frontier, binary/galloping search or vector comparison code.
            let mut remaining = std::collections::BTreeMap::<u32, usize>::new();
            for &value in b {
                *remaining.entry(value).or_default() += 1;
            }
            let mut expected = Vec::new();
            for &value in a {
                if let Some(count) = remaining.get_mut(&value)
                    && *count > 0
                {
                    expected.push(value);
                    *count -= 1;
                }
            }
            intersect_sorted(a, b, out);
            assert_eq!(*out, expected);
            let mut inplace = a.to_vec();
            intersect_sorted_in_place(&mut inplace, b);
            assert_eq!(inplace, expected);
            let mut same = a.to_vec();
            intersect_sorted_in_place(&mut same, a);
            assert_eq!(same, a);
            intersect_sorted(a, a, out);
            assert_eq!(*out, a);
        }
        let mut out = vec![u32::MAX; 100];
        // Preserve all existing varied-ratio/duplicate/reused-buffer cases.
        for n in [0, 1, 2, 7, 15, 16, 17, 31, 67, 257] {
            for m in [0, 1, 2, 7, 15, 16, 17, 31, 67, 257] {
                for salt in [0, 1, 17] {
                    let mut a: Vec<_> = (0..n).map(|i| (i * 13 + salt) % 37).collect();
                    let mut b: Vec<_> = (0..m).map(|i| (i * 7 + salt) % 43).collect();
                    a.sort_unstable();
                    b.sort_unstable();
                    check(&a, &b, &mut out);
                }
            }
        }
        // Exercise equal prefix/self writes, then compact by one or more
        // positions across block/tail boundaries, with duplicates and high IDs.
        for length in [3usize, 4, 5, 15, 16, 17, 31, 32, 33, 63, 64, 65, 257, 4097] {
            for multiplicity in [1usize, 2, 3] {
                let a: Vec<_> = (0..length).map(|i| (i / multiplicity) as u32).collect();
                for removed in [0, 1, length / 2, length - 1] {
                    let mut b = a.clone();
                    b.remove(removed);
                    b.push(u32::MAX);
                    check(&a, &b, &mut out);
                    check(&b, &a, &mut out);
                }
                let mut b = a.clone();
                b.extend([u32::MAX - 1, u32::MAX, u32::MAX]);
                check(&a, &b, &mut out);
                check(&b, &a, &mut out);
            }
        }
    }

    #[test]
    fn grouped_read_encoding_matches_independent_byte_oracle() {
        for k in 0..=32 {
            let mut seq: Vec<u8> = (0..k).map(|i| b"ACGTacgt"[(i * 7 + k) % 8]).collect();
            assert_eq!(encode_kmer_pair(&seq), encode_kmer_pair_byte_oracle(&seq));
            for pos in 0..k {
                let original = seq[pos];
                for byte in 0..=255u8 {
                    seq[pos] = byte;
                    assert_eq!(
                        encode_kmer_pair(&seq),
                        encode_kmer_pair_byte_oracle(&seq),
                        "k={k} pos={pos} byte={byte}"
                    );
                }
                seq[pos] = original;
            }
        }
        let mut state = 197u64;
        for k in 0..=32 {
            for _ in 0..2000 {
                let seq: Vec<u8> = (0..k)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        b"ACGTacgtN"[state as usize % 9]
                    })
                    .collect();
                assert_eq!(encode_kmer_pair(&seq), encode_kmer_pair_byte_oracle(&seq));
            }
        }
    }

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
