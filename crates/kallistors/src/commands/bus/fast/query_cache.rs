//! Bounded exact reuse of cDNA answers across barcode/UMI records and lanes.
use kallistors::pseudoalign::{BifrostIndex, PseudoalignOptions, Strand};
use std::borrow::Cow;
use std::sync::Mutex;

const SHARDS: usize = 64;
const SLOTS: usize = 8192;
const WAYS: usize = 2;
struct Entry<'index> {
    hash: usize,
    sequence: Vec<u8>,
    ec: Option<Cow<'index, [u32]>>,
}
#[derive(Default)]
struct Bucket<'index> {
    entries: [Option<Entry<'index>>; WAYS],
    replace: usize,
}
pub(super) struct QueryCache<'index> {
    index: &'index BifrostIndex,
    options: PseudoalignOptions,
    shards: Vec<Mutex<Vec<Bucket<'index>>>>,
}
impl<'index> QueryCache<'index> {
    pub(super) fn new(index: &'index BifrostIndex, options: PseudoalignOptions) -> Self {
        Self {
            index,
            options,
            shards: (0..SHARDS).map(|_| Mutex::new(Vec::new())).collect(),
        }
    }
    pub(super) fn lookup(&self, sequence: &[u8]) -> Option<Cow<'index, [u32]>> {
        if sequence.len() > 256 {
            return self.map(sequence);
        }
        let hash = query_hash(sequence);
        let shard = &self.shards[hash & (SHARDS - 1)];
        let slot = (hash >> SHARDS.trailing_zeros()) & (SLOTS - 1);
        {
            let mut entries = shard.lock().expect("BUS sequence cache worker panicked");
            if let Some(bucket) = entries.get_mut(slot)
                && let Some(way) = bucket.entries.iter().position(|entry| {
                    entry
                        .as_ref()
                        .is_some_and(|entry| entry.hash == hash && entry.sequence == sequence)
                })
            {
                bucket.replace = way ^ 1;
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
            let way = bucket.replace;
            bucket.replace ^= 1;
            if let Some(entry) = &mut bucket.entries[way] {
                entry.hash = hash;
                entry.sequence.clear();
                entry.sequence.extend_from_slice(sequence);
                entry.ec.clone_from(&ec);
            } else {
                bucket.entries[way] = Some(Entry {
                    hash,
                    sequence: sequence.to_vec(),
                    ec: ec.clone(),
                });
            }
        }
        ec
    }
    fn map(&self, sequence: &[u8]) -> Option<Cow<'index, [u32]>> {
        kallistors::pseudoalign::ec_for_sequence_bifrost_counting(
            self.index,
            sequence,
            Strand::Unstranded,
            self.options,
        )
    }
}

fn cacheable_ec(ec: &Option<Cow<'_, [u32]>>) -> bool {
    match ec {
        None | Some(Cow::Borrowed(_)) => true,
        Some(Cow::Owned(values)) => values.len() <= 256,
    }
}

fn query_hash(sequence: &[u8]) -> usize {
    let mut hash = sequence.len() as u64;
    let (words, tail) = sequence.as_chunks::<8>();
    for word in words {
        hash = (hash ^ u64::from_le_bytes(*word))
            .wrapping_mul(0x9e3779b97f4a7c15)
            .rotate_left(27);
    }
    for &base in tail {
        hash = (hash ^ u64::from(base)).wrapping_mul(0x100000001b3);
    }
    (hash ^ (hash >> 32)) as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn full_sequence_check_survives_collisions_and_replacements() {
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
        let a = &a[70..170];
        let b = &b[70..170];
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
                    sequence: sequence.to_vec(),
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
                            sequence: b.to_vec(),
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
