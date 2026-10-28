//! Bounded exact reuse of cDNA answers across barcode/UMI records and lanes.
use kallistors::pseudoalign::{BifrostIndex, PseudoalignOptions, Strand};
use std::borrow::Cow;
use std::sync::Mutex;
mod key;
use key::{QueryKey, StoredKey};

const SHARDS: usize = 64;
const SLOTS: usize = 8192;
const WAYS: usize = 3;
struct Entry<'index> {
    hash: usize,
    sequence: StoredKey,
    ec: Option<Cow<'index, [u32]>>,
}
struct Bucket<'index> {
    entries: [Option<Entry<'index>>; WAYS],
    // Three two-bit way IDs in least-to-most recently used order.
    replace: usize,
}
impl Default for Bucket<'_> {
    fn default() -> Self {
        Self {
            entries: Default::default(),
            replace: 0b10_01_00,
        }
    }
}
impl Bucket<'_> {
    #[inline]
    fn touch(&mut self, way: usize) {
        let least = self.replace & 3;
        let middle = (self.replace >> 2) & 3;
        let most = self.replace >> 4;
        if way == least {
            self.replace = middle | (most << 2) | (way << 4);
        } else if way == middle {
            self.replace = least | (most << 2) | (way << 4);
        }
    }
}
pub(in crate::commands::bus) struct QueryCache<'index> {
    index: &'index BifrostIndex,
    options: PseudoalignOptions,
    shards: Vec<Mutex<Vec<Bucket<'index>>>>,
}
impl<'index> QueryCache<'index> {
    pub(in crate::commands::bus) fn new(
        index: &'index BifrostIndex,
        options: PseudoalignOptions,
    ) -> Self {
        Self {
            index,
            options,
            shards: (0..SHARDS).map(|_| Mutex::new(Vec::new())).collect(),
        }
    }
    pub(in crate::commands::bus) fn lookup(&self, sequence: &[u8]) -> Option<Cow<'index, [u32]>> {
        if sequence.len() > 256 {
            return self.map(sequence);
        }
        let mut compressed = [0; 64];
        let key = QueryKey::new(sequence, &mut compressed);
        let hash = key.hash();
        let shard = &self.shards[hash & (SHARDS - 1)];
        let slot = (hash >> SHARDS.trailing_zeros()) & (SLOTS - 1);
        {
            let mut entries = shard.lock().expect("BUS sequence cache worker panicked");
            if let Some(bucket) = entries.get_mut(slot)
                && let Some(way) = bucket.entries.iter().position(|entry| {
                    entry
                        .as_ref()
                        .is_some_and(|entry| entry.hash == hash && key.matches(&entry.sequence))
                })
            {
                bucket.touch(way);
                return bucket.entries[way]
                    .as_ref()
                    .expect("matched cache entry")
                    .ec
                    .clone();
            }
        }
        // Mapping never holds a shared lock. Racing misses may do duplicate work.
        let ec = self.map(sequence);
        if cacheable_ec(&ec) {
            let mut entries = shard.lock().expect("BUS sequence cache worker panicked");
            if entries.is_empty() {
                entries.resize_with(SLOTS, Bucket::default);
            }
            let bucket = &mut entries[slot];
            let way = bucket.replace & 3;
            bucket.touch(way);
            if let Some(entry) = &mut bucket.entries[way] {
                entry.hash = hash;
                key.replace(&mut entry.sequence);
                entry.ec.clone_from(&ec);
            } else {
                bucket.entries[way] = Some(Entry {
                    hash,
                    sequence: key.stored(),
                    ec: ec.clone(),
                });
            }
        }
        ec
    }
    fn map(&self, sequence: &[u8]) -> Option<Cow<'index, [u32]>> {
        // The cache belongs to one run with a fixed strand policy. The counting
        // API applies that policy before either positive or negative answers
        // enter the cache, without constructing a diagnostic read trace.
        kallistors::pseudoalign::ec_for_sequence_bifrost_counting(
            self.index,
            sequence,
            Strand::Unstranded,
            self.options,
        )
    }
}

fn cacheable_ec(ec: &Option<Cow<'_, [u32]>>) -> bool {
    // Borrowed answers cost only a slice into the live index. Owned intersections
    // need a size cap so a few highly ambiguous reads cannot inflate every slot.
    match ec {
        None | Some(Cow::Borrowed(_)) => true,
        Some(Cow::Owned(values)) => values.len() <= 256,
    }
}

#[cfg(test)]
fn query_hash(sequence: &[u8]) -> usize {
    QueryKey::new(sequence, &mut [0; 64]).hash()
}
#[cfg(test)]
fn stored_key(sequence: &[u8]) -> StoredKey {
    QueryKey::new(sequence, &mut [0; 64]).stored()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn full_sequence_check_survives_collisions_and_replacements() {
        let mut bucket = Bucket::default();
        let mut expected = vec![0, 1, 2];
        for way in [0, 1, 2, 0, 2, 1, 1, 2, 0] {
            expected.retain(|&old| old != way);
            expected.push(way);
            bucket.touch(way);
            assert_eq!(
                [
                    bucket.replace & 3,
                    (bucket.replace >> 2) & 3,
                    bucket.replace >> 4
                ],
                expected.as_slice()
            );
        }
        let large: Vec<u32> = (0..257).collect();
        assert!(cacheable_ec(&Some(Cow::Borrowed(&large))));
        assert!(!cacheable_ec(&Some(Cow::Owned(large))));
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("reference.fa");
        let path = dir.path().join("reference.idx");
        let mut state = 81u64;
        let mut generate = || -> Vec<u8> {
            (0..300)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect()
        };
        let a = generate();
        let b = generate();
        std::fs::write(
            &fasta,
            format!(
                ">a\n{}\n>b\n{}\n",
                String::from_utf8_lossy(&a),
                String::from_utf8_lossy(&b)
            ),
        )
        .unwrap();
        kallistors::index::build_index(&path, &[fasta], Default::default()).unwrap();
        let index = kallistors::pseudoalign::build_bifrost_index_threaded(&path, 1).unwrap();
        let cache = QueryCache::new(&index, Default::default());
        let a = &a[70..161];
        let b = &b[70..161];
        assert_eq!(cache.map(a).unwrap().as_ref(), &[0]);
        assert_eq!(cache.map(b).unwrap().as_ref(), &[1]);
        let hash = query_hash(a);
        let shard = &cache.shards[hash & (SHARDS - 1)];
        let slot = (hash >> SHARDS.trailing_zeros()) & (SLOTS - 1);
        shard.lock().unwrap().resize_with(SLOTS, Bucket::default);
        // Force foreign positive/negative entries into this query's slot.
        let negative = [b'N'; 100];
        for (sequence, ec) in [(b, Some(Cow::Borrowed(&[1][..]))), (&negative[..], None)] {
            for entry in &mut shard.lock().unwrap()[slot].entries {
                *entry = Some(Entry {
                    hash,
                    sequence: stored_key(sequence),
                    ec: ec.clone(),
                });
            }
            assert_eq!(cache.lookup(a).unwrap().as_ref(), &[0]);
        }
        std::thread::scope(|scope| {
            scope.spawn(|| {
                for _ in 0..1000 {
                    for entry in &mut shard.lock().unwrap()[slot].entries {
                        *entry = Some(Entry {
                            hash,
                            sequence: stored_key(b),
                            ec: Some(Cow::Borrowed(&[1])),
                        });
                    }
                }
            });
            for _ in 0..1000 {
                assert_eq!(cache.lookup(a).unwrap().as_ref(), &[0]);
            }
        });
        assert!(cache.lookup(&[b'N'; 100]).is_none());
        assert!(cache.lookup(&[b'N'; 100]).is_none());
    }
}
