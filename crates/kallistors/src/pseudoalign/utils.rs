use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use crate::bias::BiasCounts;

use super::EcCounts;

// EC keys contain internally generated numeric transcript IDs; HashMap still
// compares the complete vectors, including when hashes collide.
#[derive(Default)]
pub(crate) struct EcKeyHasher(u64);

impl EcKeyHasher {
    #[inline]
    fn add(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

impl Hasher for EcKeyHasher {
    #[inline]
    fn finish(&self) -> u64 {
        let hash = (self.0 ^ (self.0 >> 33)).wrapping_mul(0xff51_afd7_ed55_8ccd);
        hash ^ (hash >> 33)
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let (chunks, rest) = bytes.as_chunks::<8>();
        for &chunk in chunks {
            self.add(u64::from_ne_bytes(chunk));
        }
        if !rest.is_empty() {
            let mut word = [0; 8];
            word[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_ne_bytes(word));
        }
    }

    #[inline]
    fn write_usize(&mut self, value: usize) {
        self.add(value as u64);
    }
}

pub(crate) type EcCountMap = HashMap<Vec<u32>, usize, BuildHasherDefault<EcKeyHasher>>;

pub(crate) fn add_ec_count(target: &mut EcCounts, target_map: &mut EcCountMap, ec: Vec<u32>) {
    let ec_id = match target_map.get(&ec) {
        Some(id) => *id,
        None => {
            let id = target.ec_list.len();
            target_map.insert(ec.clone(), id);
            target.ec_list.push(ec);
            target.counts.push(0);
            id
        }
    };
    target.counts[ec_id] = target.counts[ec_id].saturating_add(1);
}

pub(crate) fn merge_ec_counts(target: &mut EcCounts, target_map: &mut EcCountMap, other: EcCounts) {
    target.reads_processed = target.reads_processed.saturating_add(other.reads_processed);
    target.reads_aligned = target.reads_aligned.saturating_add(other.reads_aligned);
    for (idx, ec) in other.ec_list.into_iter().enumerate() {
        let count = other.counts.get(idx).copied().unwrap_or(0);
        if let Some(&id) = target_map.get(&ec) {
            if let Some(slot) = target.counts.get_mut(id) {
                *slot = slot.saturating_add(count);
            }
        } else {
            let id = target.ec_list.len();
            target_map.insert(ec.clone(), id);
            target.ec_list.push(ec);
            target.counts.push(count);
        }
    }
    match (&mut target.bias, other.bias) {
        (Some(target_bias), Some(other_bias)) => merge_bias_counts(target_bias, &other_bias),
        (None, Some(other_bias)) => {
            let mut bias = BiasCounts::new();
            merge_bias_counts(&mut bias, &other_bias);
            target.bias = Some(bias);
        }
        _ => {}
    }
    match (
        &mut target.fragment_length_stats,
        other.fragment_length_stats,
    ) {
        (Some(target_stats), Some(other_stats)) => target_stats.merge(&other_stats),
        (None, Some(other_stats)) => target.fragment_length_stats = Some(other_stats),
        _ => {}
    }
    match (&mut target.fragment_length_hist, other.fragment_length_hist) {
        (Some(target_hist), Some(other_hist)) => merge_fragment_hist(target_hist, &other_hist),
        (None, Some(other_hist)) => target.fragment_length_hist = Some(other_hist),
        _ => {}
    }
}

fn merge_bias_counts(target: &mut BiasCounts, other: &BiasCounts) {
    let len = target.counts.len().min(other.counts.len());
    for i in 0..len {
        target.counts[i] = target.counts[i].saturating_add(other.counts[i]);
    }
    target.total = target.total.saturating_add(other.total);
}

fn merge_fragment_hist(target: &mut [u32], other: &[u32]) {
    let len = target.len().min(other.len());
    for i in 0..len {
        target[i] = target[i].saturating_add(other[i]);
    }
}

#[cfg(test)]
mod count_map_tests {
    use super::*;

    #[test]
    fn numeric_map_keeps_complete_ec_keys_and_counts() {
        let mut fast = EcCountMap::default();
        let mut standard = HashMap::new();
        for i in 0..20_000u32 {
            let key = vec![0, i % 193, (i * 17) % 197, (i * 31) % 199];
            *fast.entry(key.clone()).or_insert(0) += 1;
            *standard.entry(key).or_insert(0) += 1;
        }
        assert_eq!(fast.len(), standard.len());
        for (key, count) in standard {
            assert_eq!(fast.get(&key), Some(&count));
        }
        for key in [vec![], vec![0], vec![0, 0], vec![0, 0, 0]] {
            fast.insert(key.clone(), key.len());
        }
        assert_eq!(fast.get(&vec![0, 0]), Some(&2));
    }
}
