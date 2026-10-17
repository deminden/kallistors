//! Conservative bucket filtering; false positives retain ordinary exact matching.

use std::sync::atomic::{AtomicU8, Ordering};

pub(crate) const DISABLED: u8 = u8::MAX;

#[inline]
pub(crate) fn word_fingerprint(word: u64) -> u8 {
    (word.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 56) as u8
}

pub(crate) fn observe(slot: &AtomicU8, word: u64) {
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
    first: u64,
    last: u8,
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
                result.first |= u64::from(fingerprint) << (8 * position);
            } else {
                result.last = fingerprint;
            }
        }
        result
    }

    #[inline]
    fn contains(&self, fingerprint: u8) -> bool {
        let different = self.first ^ (u64::from(fingerprint) * 0x0101_0101_0101_0101);
        let equal = !((different & 0x7f7f_7f7f_7f7f_7f7f).wrapping_add(0x7f7f_7f7f_7f7f_7f7f)
            | different)
            & 0x8080_8080_8080_8080;
        equal != 0 || self.last == fingerprint
    }
}

#[derive(Default)]
pub(crate) struct QueryFingerprintCache {
    key: Option<u64>,
    fingerprints: QueryFingerprints,
}

impl QueryFingerprintCache {
    #[inline]
    pub(crate) fn contains(&mut self, forward: u64, reverse: u64, fingerprint: u8) -> bool {
        if self.key != Some(forward) {
            self.fingerprints = QueryFingerprints::new(forward, reverse);
            self.key = Some(forward);
        }
        self.fingerprints.contains(fingerprint)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_anchor_and_byte_membership_match_independent_oracles() {
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
            for fingerprint in 0..=u8::MAX {
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
        let words: Vec<_> = (1..10000u64)
            .map(|word| word << 18)
            .filter(|&word| word_fingerprint(word) != 0 && word_fingerprint(word) != DISABLED)
            .take(2)
            .collect();
        assert_ne!(word_fingerprint(words[0]), word_fingerprint(words[1]));
        let slot = AtomicU8::new(0);
        std::thread::scope(|scope| {
            for &word in &words {
                let slot = &slot;
                scope.spawn(move || observe(slot, word));
            }
        });
        assert_eq!(slot.load(Ordering::Relaxed), DISABLED);
        observe(&slot, words[0]);
        assert_eq!(slot.load(Ordering::Relaxed), DISABLED);
        let special = AtomicU8::new(DISABLED);
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
    fn fingerprint(self) -> u8;
}
impl BucketCountKey for u32 {
    const FINGERPRINT: bool = false;
    fn new(bucket: u32, _: u64) -> Self {
        bucket
    }
    fn bucket(self) -> u32 {
        self
    }
    fn fingerprint(self) -> u8 {
        0
    }
}
#[derive(Clone, Copy)]
pub(crate) struct TaggedBucket(u64);
impl rdst::RadixKey for TaggedBucket {
    const LEVELS: usize = 5;
    fn get_level(&self, level: usize) -> u8 {
        (self.0 >> (level * 8)) as u8
    }
}
impl BucketCountKey for TaggedBucket {
    const FINGERPRINT: bool = true;
    fn new(bucket: u32, word: u64) -> Self {
        Self((u64::from(bucket) << 8) | u64::from(word_fingerprint(word)))
    }
    fn bucket(self) -> u32 {
        (self.0 >> 8) as u32
    }
    fn fingerprint(self) -> u8 {
        self.0 as u8
    }
}
// Called by the serial count merger after every parallel chunk has completed.
// No concurrent fingerprint writes exist during this phase.
pub(crate) fn observe_serial(slot: &AtomicU8, fingerprint: u8) {
    let old = slot.load(Ordering::Relaxed);
    if fingerprint == 0 || fingerprint == DISABLED || (old != 0 && old != fingerprint) {
        slot.store(DISABLED, Ordering::Relaxed);
    } else if old == 0 {
        slot.store(fingerprint, Ordering::Relaxed);
    }
}
