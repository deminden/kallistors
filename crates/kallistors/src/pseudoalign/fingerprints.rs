//! Conservative bucket filtering; false positives retain ordinary exact matching.

use std::sync::atomic::{AtomicU16, Ordering};

pub(crate) const DISABLED: u16 = (1 << 12) - 1;

#[inline]
pub(crate) fn word_fingerprint(word: u64) -> u16 {
    (word.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 52) as u16
}

pub(crate) fn observe(slot: &AtomicU16, word: u64) {
    let fingerprint = word_fingerprint(word);
    if fingerprint == 0 || fingerprint == DISABLED {
        slot.store(DISABLED, Ordering::Relaxed);
        return;
    }
    loop {
        let current = slot.load(Ordering::Relaxed);
        if current == DISABLED || current == fingerprint {
            return;
        }
        if current != 0 {
            slot.store(DISABLED, Ordering::Relaxed);
            return;
        }
        if slot
            .compare_exchange(0, fingerprint, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
        {
            return;
        }
    }
}

#[derive(Default)]
struct QueryFingerprints {
    first: [u16; 8],
    last: u16,
}
impl QueryFingerprints {
    fn new(forward: u64, reverse: u64) -> Self {
        const MASK: u64 = (1u64 << 46) - 1;
        let mut result = Self::default();
        for position in 0..=8 {
            let f = (forward >> (2 * (8 - position))) & MASK;
            let r = (reverse >> (2 * position)) & MASK;
            let fingerprint = word_fingerprint(f.min(r) << 18);
            if position < 8 {
                result.first[position] = fingerprint;
            } else {
                result.last = fingerprint;
            }
        }
        result
    }
    #[inline]
    fn advance(&mut self, last: u16) {
        #[cfg(target_arch = "aarch64")]
        {
            use std::arch::aarch64::*;
            // SAFETY: both vector memory accesses cover exactly first[0..8].
            // NEON is baseline on AArch64. The lane after old first[7] is
            // supplied from a duplicated local scalar, never an overread.
            unsafe {
                let prior = vld1q_u16(self.first.as_ptr());
                let shifted = vextq_u16::<1>(prior, vdupq_n_u16(self.last));
                vst1q_u16(self.first.as_mut_ptr(), shifted);
            }
        }
        #[cfg(not(target_arch = "aarch64"))]
        {
            self.first.copy_within(1..8, 0);
            self.first[7] = self.last;
        }
        self.last = last;
    }

    #[inline]
    fn contains(&self, fingerprint: u16) -> bool {
        #[cfg(target_arch = "aarch64")]
        {
            use std::arch::aarch64::*;
            // SAFETY: first is exactly eight initialized u16 values; NEON
            // permits this unaligned bounded load on every AArch64 target.
            let any = unsafe {
                vmaxvq_u16(vceqq_u16(
                    vld1q_u16(self.first.as_ptr()),
                    vdupq_n_u16(fingerprint),
                ))
            };
            any != 0 || self.last == fingerprint
        }
        #[cfg(not(target_arch = "aarch64"))]
        {
            self.first.contains(&fingerprint) || self.last == fingerprint
        }
    }
}
#[derive(Default)]
pub(crate) struct QueryFingerprintCache {
    key: Option<u64>,
    fingerprints: QueryFingerprints,
}
impl QueryFingerprintCache {
    #[inline]
    pub(crate) fn contains(&mut self, forward: u64, reverse: u64, fingerprint: u16) -> bool {
        if self.key != Some(forward) {
            const OVERLAP: u64 = (1u64 << 60) - 1;
            if self
                .key
                .is_some_and(|previous| previous & OVERLAP == forward >> 2)
            {
                const MASK: u64 = (1u64 << 46) - 1;
                self.fingerprints.advance(word_fingerprint(
                    (forward & MASK).min((reverse >> 16) & MASK) << 18,
                ));
            } else {
                self.fingerprints = QueryFingerprints::new(forward, reverse);
            }
            self.key = Some(forward);
        }
        self.fingerprints.contains(fingerprint)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_fingerprint_advance_matches_independent_scalar_lanes() {
        for base in [0u16, 1, 255, 4094, 4095, u16::MAX - 8] {
            for prior_last in [0u16, 1, 4095, u16::MAX] {
                let mut values = QueryFingerprints {
                    first: std::array::from_fn(|i| base.wrapping_add(i as u16)),
                    last: prior_last,
                };
                for new_last in [0u16, 1, 4095, u16::MAX, 7, 31, 0, 4095, 9, 10] {
                    let old = values.first;
                    let expected = [
                        old[1],
                        old[2],
                        old[3],
                        old[4],
                        old[5],
                        old[6],
                        old[7],
                        values.last,
                    ];
                    values.advance(new_last);
                    assert_eq!(values.first, expected);
                    assert_eq!(values.last, new_last);
                }
            }
        }
    }

    #[test]
    fn every_anchor_and_fingerprint_membership_match_independent_oracles() {
        let sequence: Vec<_> = (0..191).map(|i| b"ACGT"[(i * i + i / 7) % 4]).collect();
        let mut rolling = QueryFingerprintCache::default();
        for reversed in [false, true] {
            for position in (0..=160).chain([19, 20, 21, 1, 2, 3, 159, 160]) {
                let mut query = sequence[position..position + 31].to_vec();
                if reversed {
                    query.reverse();
                    for byte in &mut query {
                        *byte = match *byte {
                            b'A' => b'T',
                            b'C' => b'G',
                            b'G' => b'C',
                            _ => b'A',
                        };
                    }
                }
                let (forward, reverse) = crate::dna::encode_kmer_pair(&query).unwrap();
                let expected: Vec<_> = query
                    .windows(23)
                    .map(|window| {
                        word_fingerprint(u64::from_le_bytes(
                            crate::index::bifrost::encode_minimizer_rep(window).unwrap(),
                        ))
                    })
                    .collect();
                for fingerprint in 0..=DISABLED {
                    assert_eq!(
                        rolling.contains(forward, reverse, fingerprint),
                        expected.contains(&fingerprint),
                        "position={position} reversed={reversed}"
                    );
                }
            }
        }
        let mut state = 17u64;
        let mut cache = QueryFingerprintCache::default();
        for _ in 0..1000 {
            let seq: Vec<_> = (0..31)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect();
            let (forward, reverse) = crate::dna::encode_kmer_pair(&seq).unwrap();
            let expected: Vec<_> = seq
                .windows(23)
                .map(|window| {
                    word_fingerprint(u64::from_le_bytes(
                        crate::index::bifrost::encode_minimizer_rep(window).unwrap(),
                    ))
                })
                .collect();
            // Rolling cases above exhaust the complete fingerprint domain.
            // Random words cover every real anchor plus absent/boundary values.
            for fingerprint in expected.iter().copied().chain([
                0,
                1,
                255,
                256,
                DISABLED - 1,
                DISABLED,
                (state as u16) & DISABLED,
                ((state >> 12) as u16) & DISABLED,
            ]) {
                assert_eq!(
                    cache.contains(forward, reverse, fingerprint),
                    expected.contains(&fingerprint)
                );
                assert_eq!(
                    cache.contains(reverse, forward, fingerprint),
                    expected.contains(&fingerprint)
                );
            }
        }
    }

    #[test]
    fn conflicting_and_special_buckets_stay_disabled_under_parallel_updates() {
        for bucket in [0, 1, u32::MAX] {
            for word in [0, 1, 17, 123456789, u64::MAX] {
                let tagged = TaggedBucket::new(bucket, word);
                assert_eq!(tagged.bucket(), bucket);
                assert_eq!(tagged.fingerprint(), word_fingerprint(word));
                let decoded =
                    (0..<TaggedBucket as rdst::RadixKey>::LEVELS).fold(0u64, |v, level| {
                        v | (u64::from(rdst::RadixKey::get_level(&tagged, level)) << (8 * level))
                    });
                assert_eq!(decoded, u64::from(bucket));
            }
        }
        let words: Vec<_> = (1..10000u64)
            .map(|word| word << 18)
            .filter(|&word| word_fingerprint(word) != 0 && word_fingerprint(word) != DISABLED)
            .take(2)
            .collect();
        assert_ne!(word_fingerprint(words[0]), word_fingerprint(words[1]));
        let slot = AtomicU16::new(0);
        std::thread::scope(|scope| {
            for &word in &words {
                let slot = &slot;
                scope.spawn(move || observe(slot, word));
            }
        });
        assert_eq!(slot.load(Ordering::Relaxed), DISABLED);
        observe(&slot, words[0]);
        assert_eq!(slot.load(Ordering::Relaxed), DISABLED);
        let special = AtomicU16::new(DISABLED);
        observe(&special, words[1]);
        assert_eq!(special.load(Ordering::Relaxed), DISABLED);
    }
}

// Fingerprint construction travels with the existing sorted count batches.
// Ordinary loaders retain their original compact u32 scratch representation.
pub(crate) trait BucketCountKey: rdst::RadixKey + Send + Sync + Copy {
    const FINGERPRINT: bool;
    fn new(bucket: u32, word: u64) -> Self;
    fn bucket(self) -> u32;
    fn fingerprint(self) -> u16;
}
impl BucketCountKey for u32 {
    const FINGERPRINT: bool = false;
    fn new(bucket: u32, _: u64) -> Self {
        bucket
    }
    fn bucket(self) -> u32 {
        self
    }
    fn fingerprint(self) -> u16 {
        0
    }
}
#[derive(Clone, Copy)]
pub(crate) struct TaggedBucket(u64);
impl rdst::RadixKey for TaggedBucket {
    // Only bucket identity is required by count grouping. Fingerprint
    // collisions are checked across the complete group by its owner.
    const LEVELS: usize = 4;
    fn get_level(&self, level: usize) -> u8 {
        (self.0 >> (16 + level * 8)) as u8
    }
}
impl BucketCountKey for TaggedBucket {
    const FINGERPRINT: bool = true;
    fn new(bucket: u32, word: u64) -> Self {
        Self((u64::from(bucket) << 16) | u64::from(word_fingerprint(word)))
    }
    fn bucket(self) -> u32 {
        (self.0 >> 16) as u32
    }
    fn fingerprint(self) -> u16 {
        self.0 as u16
    }
}
// The caller gives every bucket exactly one updating owner. Other buckets
// may be merged in parallel; simultaneous updates cannot target this slot.
pub(crate) fn observe_serial(slot: &AtomicU16, fingerprint: u16) {
    let old = slot.load(Ordering::Relaxed);
    if fingerprint == 0 || fingerprint == DISABLED || (old != 0 && old != fingerprint) {
        slot.store(DISABLED, Ordering::Relaxed);
    } else if old == 0 {
        slot.store(fingerprint, Ordering::Relaxed);
    }
}
