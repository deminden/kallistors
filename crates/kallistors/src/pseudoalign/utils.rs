use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::sync::{Arc, Mutex};

use crate::bias::BiasCounts;

use super::EcCounts;

// EC keys contain internally generated numeric transcript IDs. Complete vector
// checks below establish identity, including when fingerprints collide.
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

#[derive(Default)]
struct EcFingerprintHasher(u64);
impl Hasher for EcFingerprintHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(byte);
        }
    }
    fn write_u64(&mut self, value: u64) {
        self.0 = value;
    }
    fn write_u32(&mut self, value: u32) {
        // Mix into both halves so SwissTable's high-bit control tags retain entropy.
        self.0 = u64::from(value).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    }
}

pub(crate) struct SharedEcKey {
    ec: Box<[u32]>,
    previous: Option<Arc<SharedEcKey>>,
}
pub(crate) type SharedEcKeys = Vec<Arc<SharedEcKey>>;
type SharedShard = Mutex<HashMap<u64, Arc<SharedEcKey>, BuildHasherDefault<EcFingerprintHasher>>>;
pub(crate) struct SharedEcPool {
    shards: [SharedShard; 64],
}
impl Default for SharedEcPool {
    fn default() -> Self {
        Self {
            shards: std::array::from_fn(|_| Mutex::new(HashMap::default())),
        }
    }
}
impl SharedEcPool {
    fn intern(&self, hash: u64, ec: Vec<u32>) -> Arc<SharedEcKey> {
        let mut shard = self.shards[(hash >> 58) as usize]
            .lock()
            .expect("EC pool poisoned");
        let mut node = shard.get(&hash);
        while let Some(key) = node {
            if key.ec.as_ref() == ec.as_slice() {
                return Arc::clone(key);
            }
            node = key.previous.as_ref();
        }
        let previous = shard.get(&hash).cloned();
        let key = Arc::new(SharedEcKey {
            ec: ec.into_boxed_slice(),
            previous,
        });
        shard.insert(hash, Arc::clone(&key));
        key
    }
    pub(crate) fn clear(&self) {
        for shard in &self.shards {
            shard.lock().expect("EC pool poisoned").clear();
        }
    }
}
struct SharedWorker {
    pool: Arc<SharedEcPool>,
    keys: SharedEcKeys,
}

// The public EC list owns each key once. Fingerprints select candidates;
// complete transcript vectors establish identity, including hash collisions.
#[derive(Default)]
pub(crate) struct EcCountMap {
    entries: HashMap<u32, u32, BuildHasherDefault<EcFingerprintHasher>>,
    wide_heads: HashMap<u32, usize>,
    collisions: HashMap<usize, usize>,
    shared: Option<SharedWorker>,
}
impl EcCountMap {
    pub(crate) fn shared(pool: Arc<SharedEcPool>) -> Self {
        Self {
            shared: Some(SharedWorker {
                pool,
                keys: Vec::new(),
            }),
            ..Self::default()
        }
    }
    pub(crate) fn into_shared_keys(self) -> SharedEcKeys {
        self.shared.expect("shared worker EC map").keys
    }
    fn store_new(&mut self, ec: Vec<u32>, classes: &mut Vec<Vec<u32>>) {
        if let Some(worker) = self.shared.as_mut() {
            worker
                .keys
                .push(worker.pool.intern(Self::fingerprint(&ec), ec));
        } else {
            classes.push(ec);
        }
    }
    fn folded(hash: u64) -> u32 {
        hash as u32 ^ (hash >> 32) as u32
    }
    fn head_id(raw: u32, hash: u32, wide: &HashMap<u32, usize>) -> usize {
        if raw == u32::MAX {
            wide[&hash]
        } else {
            raw as usize
        }
    }
    fn encode_head(id: usize, hash: u32, wide: &mut HashMap<u32, usize>) -> u32 {
        if id >= u32::MAX as usize {
            wide.insert(hash, id);
            u32::MAX
        } else {
            if !wide.is_empty() {
                wide.remove(&hash);
            }
            id as u32
        }
    }
    fn fingerprint(ec: &[u32]) -> u64 {
        let mut hash = EcKeyHasher::default();
        ec.hash(&mut hash);
        hash.finish()
    }
    // The caller appends the complete key when the returned `new` flag is true.
    // A single hash/table probe covers both the lookup and insertion paths.
    fn get_or_insert(&mut self, ec: &[u32], classes: &[Vec<u32>]) -> (usize, bool) {
        self.get_or_insert_hashed(Self::fingerprint(ec), ec, classes)
    }
    fn get_or_insert_hashed(
        &mut self,
        hash: u64,
        ec: &[u32],
        classes: &[Vec<u32>],
    ) -> (usize, bool) {
        let keys = self.shared.as_ref().map(|w| &w.keys);
        let next = keys.map_or(classes.len(), Vec::len);
        let hash = Self::folded(hash);
        match self.entries.entry(hash) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(Self::encode_head(next, hash, &mut self.wide_heads));
            }
            std::collections::hash_map::Entry::Occupied(mut entry) => {
                let head = Self::head_id(*entry.get(), hash, &self.wide_heads);
                let mut id = head;
                loop {
                    if keys.map_or_else(|| classes[id].as_slice(), |k| k[id].ec.as_ref()) == ec {
                        return (id, false);
                    }
                    match self.collisions.get(&id) {
                        Some(&previous) => id = previous,
                        None => break,
                    }
                }
                self.collisions.insert(next, head);
                entry.insert(Self::encode_head(next, hash, &mut self.wide_heads));
            }
        }
        (next, true)
    }
    pub(crate) fn get(&self, ec: &[u32], classes: &[Vec<u32>]) -> Option<usize> {
        self.get_hashed(Self::fingerprint(ec), ec, classes)
    }
    fn get_hashed(&self, hash: u64, ec: &[u32], classes: &[Vec<u32>]) -> Option<usize> {
        let keys = self.shared.as_ref().map(|w| &w.keys);
        let hash = Self::folded(hash);
        let mut id = Self::head_id(*self.entries.get(&hash)?, hash, &self.wide_heads);
        loop {
            if keys.map_or_else(|| classes[id].as_slice(), |k| k[id].ec.as_ref()) == ec {
                return Some(id);
            }
            id = *self.collisions.get(&id)?;
        }
    }
    pub(crate) fn insert(&mut self, ec: &[u32], id: usize) {
        self.insert_hashed(Self::fingerprint(ec), id);
    }
    fn insert_hashed(&mut self, hash: u64, id: usize) {
        let hash = Self::folded(hash);
        let previous = self
            .entries
            .get(&hash)
            .map(|&raw| Self::head_id(raw, hash, &self.wide_heads));
        let raw = Self::encode_head(id, hash, &mut self.wide_heads);
        self.entries.insert(hash, raw);
        if let Some(previous) = previous {
            self.collisions.insert(id, previous);
        }
    }
}

pub(crate) fn add_ec_count(target: &mut EcCounts, target_map: &mut EcCountMap, ec: Vec<u32>) {
    let (ec_id, new) = target_map.get_or_insert(&ec, &target.ec_list);
    if new {
        target_map.store_new(ec, &mut target.ec_list);
        target.counts.push(0);
    }
    target.counts[ec_id] = target.counts[ec_id].saturating_add(1);
}

// Reusable paired-EC scratch remains with the worker; only new keys are copied.
pub(crate) fn add_ec_count_slice(target: &mut EcCounts, target_map: &mut EcCountMap, ec: &[u32]) {
    let (ec_id, new) = target_map.get_or_insert(ec, &target.ec_list);
    if new {
        target_map.store_new(ec.to_vec(), &mut target.ec_list);
        target.counts.push(0);
    }
    target.counts[ec_id] = target.counts[ec_id].saturating_add(1);
}

pub(crate) fn merge_ec_counts(target: &mut EcCounts, target_map: &mut EcCountMap, other: EcCounts) {
    target.reads_processed = target.reads_processed.saturating_add(other.reads_processed);
    target.reads_aligned = target.reads_aligned.saturating_add(other.reads_aligned);
    for (idx, ec) in other.ec_list.into_iter().enumerate() {
        let count = other.counts.get(idx).copied().unwrap_or(0);
        let (id, new) = target_map.get_or_insert(&ec, &target.ec_list);
        if new {
            target_map.store_new(ec, &mut target.ec_list);
            target.counts.push(count);
        } else if let Some(slot) = target.counts.get_mut(id) {
            *slot = slot.saturating_add(count);
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

pub(crate) fn merge_shared_ec_counts(
    target: &mut EcCounts,
    target_map: &mut EcCountMap,
    other: EcCounts,
    keys: SharedEcKeys,
) {
    debug_assert!(other.ec_list.is_empty());
    for (key, &count) in keys.into_iter().zip(&other.counts) {
        let (id, new) = target_map.get_or_insert(&key.ec, &target.ec_list);
        if new {
            let ec =
                Arc::try_unwrap(key).map_or_else(|key| key.ec.to_vec(), |key| key.ec.into_vec());
            target.ec_list.push(ec);
            target.counts.push(count);
        } else {
            target.counts[id] = target.counts[id].saturating_add(count);
        }
    }
    merge_ec_counts(target, target_map, other);
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
    fn shared_workers_keep_complete_collision_keys_and_merge_counts() {
        let pool = Arc::new(SharedEcPool::default());
        let first = pool.intern(123, vec![1, 2]);
        let different = pool.intern(123, vec![2, 1]);
        assert!(!Arc::ptr_eq(&first, &different));
        assert!(Arc::ptr_eq(&first, &pool.intern(123, vec![1, 2])));
        let mut target = EcCounts::new(false);
        let mut map = EcCountMap::default();
        for _ in 0..3 {
            let mut worker = EcCounts::new(false);
            let mut worker_map = EcCountMap::shared(Arc::clone(&pool));
            for key in [&[1, 2][..], &[2, 1][..], &[1, 2][..], &[]] {
                add_ec_count_slice(&mut worker, &mut worker_map, key);
            }
            assert!(worker.ec_list.is_empty());
            worker.reads_processed = 4;
            merge_shared_ec_counts(&mut target, &mut map, worker, worker_map.into_shared_keys());
        }
        pool.clear();
        assert_eq!(first.ec.as_ref(), [1, 2]);
        assert_eq!(different.ec.as_ref(), [2, 1]);
        assert_eq!(target.reads_processed, 12);
        for (key, count) in [(vec![1, 2], 6), (vec![2, 1], 3), (vec![], 3)] {
            assert_eq!(
                target.counts[map.get(&key, &target.ec_list).unwrap()],
                count
            );
        }
    }

    #[test]
    fn numeric_map_keeps_complete_ec_keys_and_counts() {
        let mut fast = EcCountMap::default();
        let mut target = EcCounts::new(false);
        let mut standard = HashMap::new();
        for i in 0..20_000u32 {
            let key = vec![0, i % 193, (i * 17) % 197, (i * 31) % 199];
            add_ec_count_slice(&mut target, &mut fast, &key);
            *standard.entry(key).or_insert(0u32) += 1;
        }
        assert_eq!(target.ec_list.len(), standard.len());
        for (key, count) in standard {
            let id = fast.get(&key, &target.ec_list).unwrap();
            assert_eq!(target.counts[id], count);
        }
        let keys = vec![vec![], vec![0], vec![0, 0], vec![0, 0, 0]];
        let mut forced = EcCountMap::default();
        for id in 0..keys.len() {
            assert_eq!(
                forced.get_or_insert_hashed(123, &keys[id], &keys[..id]),
                (id, true)
            );
        }
        for (id, key) in keys.iter().enumerate() {
            assert_eq!(forced.get_or_insert_hashed(123, key, &keys), (id, false));
        }
        for (id, key) in keys.iter().enumerate() {
            assert_eq!(forced.get_hashed(123, key, &keys), Some(id));
        }
        assert_eq!(forced.get_hashed(123, &[1], &keys), None);
        assert_eq!(forced.get_hashed(124, &[0], &keys), None);
        // Distinct full fingerprints can fold to the same key; vectors still decide.
        assert_eq!(
            forced.get_or_insert_hashed(123u64 << 32, &keys[0], &keys),
            (0, false)
        );
        let mut wide = HashMap::new();
        for id in [0, (u32::MAX - 1) as usize, u32::MAX as usize, usize::MAX, 1] {
            let raw = EcCountMap::encode_head(id, 7, &mut wide);
            assert_eq!(EcCountMap::head_id(raw, 7, &wide), id);
        }
    }
}
