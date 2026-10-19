//! Exact prefix-relative keys for large-graph topology and transcript mapping.
use super::super::workers::IndexWorkers;
use crate::Result;
use rayon::prelude::*;

use super::{KmerDirectory, KmerLookup, LINEAR_KMER_LOOKUP_MAX, prefetch_rank_delta};

pub(super) trait MappingKeys: Sync {
    fn len(&self) -> usize;
    fn word(&self, lookup: &KmerLookup, index: usize) -> u64;
    fn retire_word_hints(&mut self) {}
    fn raw_words(&self) -> Option<&[u64]> {
        None
    }
    fn window<'a>(
        &'a self,
        lookup: &KmerLookup,
        first: usize,
        count: usize,
        buffer: &'a mut [u64],
    ) -> &'a [u64];

    fn find(&self, lookup: &KmerLookup, word: u64) -> Option<usize>;
    fn find_many_count<const N: usize>(
        &self,
        lookup: &KmerLookup,
        words: [u64; N],
        count: usize,
    ) -> [Option<usize>; N];
    fn find_many_ranks<const N: usize>(
        &self,
        lookup: &KmerLookup,
        words: [u64; N],
        count: usize,
    ) -> [u32; N] {
        self.find_many_count(lookup, words, count)
            .map(|hit| hit.map_or(u32::MAX, |index| index as u32))
    }
}

impl MappingKeys for [u64] {
    fn raw_words(&self) -> Option<&[u64]> {
        Some(self)
    }
    #[inline]
    fn window<'a>(
        &'a self,
        _lookup: &KmerLookup,
        first: usize,
        count: usize,
        _buffer: &'a mut [u64],
    ) -> &'a [u64] {
        &self[first..first + count]
    }

    #[inline]
    fn word(&self, _lookup: &KmerLookup, index: usize) -> u64 {
        self[index]
    }
    fn len(&self) -> usize {
        self.len()
    }
    fn find(&self, lookup: &KmerLookup, word: u64) -> Option<usize> {
        lookup.find(self, word)
    }
    fn find_many_count<const N: usize>(
        &self,
        lookup: &KmerLookup,
        words: [u64; N],
        count: usize,
    ) -> [Option<usize>; N] {
        lookup.find_many_count(self, words, count)
    }
}

impl MappingKeys for Vec<u64> {
    fn raw_words(&self) -> Option<&[u64]> {
        Some(self)
    }
    #[inline]
    fn window<'a>(
        &'a self,
        _lookup: &KmerLookup,
        first: usize,
        count: usize,
        _buffer: &'a mut [u64],
    ) -> &'a [u64] {
        &self[first..first + count]
    }

    #[inline]
    fn word(&self, _lookup: &KmerLookup, index: usize) -> u64 {
        self[index]
    }
    fn len(&self) -> usize {
        self.len()
    }
    fn find(&self, lookup: &KmerLookup, word: u64) -> Option<usize> {
        lookup.find(self, word)
    }
    fn find_many_count<const N: usize>(
        &self,
        lookup: &KmerLookup,
        words: [u64; N],
        count: usize,
    ) -> [Option<usize>; N] {
        lookup.find_many_count(self, words, count)
    }
}

pub(super) struct PrefixKeys {
    // Retain the original allocation and its allocator alignment. Only the
    // initialized packed prefix remains live after unused pages are retired.
    words: Vec<u64>,
    count: usize,
    mask: u64,
    // One exact bucket hint per 64 ranks bounds full-key reconstruction.
    // Empty directory buckets are resolved by rank, never treated as matches.
    prefixes: Vec<u32>,
}

impl PrefixKeys {
    pub(super) fn supports(count: usize, lookup: &KmerLookup) -> bool {
        let bits = lookup.shift.saturating_sub(lookup.padding);
        count > (1 << 28)
            && matches!(lookup.directory, KmerDirectory::Ranks(_))
            && (1..=40).contains(&bits)
    }

    pub(super) fn for_topology(
        words: Vec<u64>,
        lookup: &KmerLookup,
        threads: usize,
    ) -> Result<Self> {
        assert!(
            matches!(lookup.directory, KmerDirectory::Ranks(_))
                && (1..=40).contains(&lookup.shift.saturating_sub(lookup.padding))
        );
        let prefixes = words
            .chunks(64)
            .map(|chunk| (chunk[0] >> lookup.shift) as u32)
            .collect();
        let mut packed = Self::pack(words, lookup, threads)?;
        packed.prefixes = prefixes;
        Ok(packed)
    }
    fn pack(mut words: Vec<u64>, lookup: &KmerLookup, threads: usize) -> Result<Self> {
        let count = words.len();
        let mask = (1u64 << (lookup.shift - lookup.padding)) - 1;
        // SAFETY: u64 permits every bit pattern. This exclusive byte view
        // covers initialized elements only; no other view survives the loop.
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(words.as_mut_ptr().cast::<u8>(), count * 8) };
        let prefix = count.min(65536);
        for index in 0..prefix {
            // The source word is read completely before its shorter result is
            // written. Five-byte destinations never overwrite a future word.
            let word = u64::from_ne_bytes(bytes[index * 8..index * 8 + 8].try_into().unwrap());
            let suffix = ((word >> lookup.padding) & mask).to_le_bytes();
            bytes[index * 5..index * 5 + 5].copy_from_slice(&suffix[..5]);
        }
        let pack_wave = |source: &[u8], dest: &mut [u8]| {
            for (source, dest) in source
                .as_chunks::<8>()
                .0
                .iter()
                .zip(dest.as_chunks_mut::<5>().0)
            {
                let word = u64::from_ne_bytes(*source);
                let suffix = ((word >> lookup.padding) & mask).to_le_bytes();
                dest.copy_from_slice(&suffix[..5]);
            }
        };
        let pack = || {
            let mut start = prefix;
            while start < count {
                // All output is below the first unread source word. The two
                // borrowed regions are disjoint; joined waves preserve this
                // invariant while the initialized packed prefix grows.
                let end = (start * 8 / 5).min(count);
                let (written, unread) = bytes.split_at_mut(start * 8);
                let dest = &mut written[start * 5..end * 5];
                let source = &unread[..(end - start) * 8];
                if threads > 1 {
                    source
                        .par_chunks(16384 * 8)
                        .zip(dest.par_chunks_mut(16384 * 5))
                        .for_each(|(source, dest)| pack_wave(source, dest));
                } else {
                    pack_wave(source, dest);
                }
                start = end;
            }
        };
        if threads > 1 && count > prefix {
            IndexWorkers::new(threads)?.install(pack);
        } else {
            let mut pack = pack;
            pack();
        }
        // Retain three initialized padding bytes after the final key. They
        // were already live in the original u64 allocation; this adds at most
        // one retained word and permits a single bounded eight-byte load.
        words.truncate(if count == 0 {
            0
        } else {
            (count * 5 + 3).div_ceil(8)
        });
        super::super::scratch_pages::discard_unused(&mut words, count);
        Ok(Self {
            words,
            count,
            mask,
            prefixes: Vec::new(),
        })
    }

    #[inline]
    fn prefix_at(&self, lookup: &KmerLookup, index: usize) -> usize {
        let KmerDirectory::Ranks(offsets) = &lookup.directory else {
            unreachable!("exact prefix rank directory")
        };
        let group = index >> 6;
        let mut first = self.prefixes[group] as usize;
        let mut end = self
            .prefixes
            .get(group + 1)
            .map_or(1usize << (64 - lookup.shift), |&prefix| prefix as usize + 1);
        while first + 1 < end {
            let mid = first + (end - first) / 2;
            if offsets.get(mid) <= index {
                first = mid;
            } else {
                end = mid;
            }
        }
        first
    }
    #[inline]
    fn packed_bytes(&self) -> &[u8] {
        // SAFETY: u8 accepts every bit pattern. Every retained u64, including
        // the final padding, was initialized before in-place packing; the
        // shared view stays within that live allocation and cannot mutate it.
        unsafe {
            std::slice::from_raw_parts(self.words.as_ptr().cast::<u8>(), self.words.len() * 8)
        }
    }

    #[inline(always)]
    fn suffix(&self, index: usize) -> u64 {
        debug_assert!(index < self.count);
        let first = index * 5;
        let bytes = &self.packed_bytes()[first..first + 8];
        u64::from_le_bytes(bytes.try_into().expect("padded key load")) & ((1 << 40) - 1)
    }

    #[inline]
    fn prefetch(&self, index: usize) {
        prefetch_rank_delta(&self.packed_bytes()[index * 5]);
    }

    #[inline]
    fn search(&self, word: u64, lookup: &KmerLookup, start: usize, end: usize) -> Option<usize> {
        if word & lookup.padding_mask != 0 {
            return None;
        }
        let suffix = (word >> lookup.padding) & self.mask;
        let mut size = end - start;
        if size <= LINEAR_KMER_LOOKUP_MAX {
            return (start..end).find(|&index| self.suffix(index) == suffix);
        }
        // The same fixed-depth lower-half rule as slice binary search, now
        // using one safe unaligned load instead of separate four/one-byte loads.
        let mut base = start;
        while size > 1 {
            let half = size / 2;
            let mid = base + half;
            base = std::hint::select_unpredictable(self.suffix(mid) > suffix, base, mid);
            size -= half;
        }
        (self.suffix(base) == suffix).then_some(base)
    }
}

impl MappingKeys for PrefixKeys {
    fn retire_word_hints(&mut self) {
        super::super::scratch_pages::retire(std::mem::take(&mut self.prefixes));
    }

    #[inline]
    fn word(&self, lookup: &KmerLookup, index: usize) -> u64 {
        ((self.prefix_at(lookup, index) as u64) << lookup.shift)
            | (self.suffix(index) << lookup.padding)
    }
    fn window<'a>(
        &'a self,
        lookup: &KmerLookup,
        first: usize,
        count: usize,
        buffer: &'a mut [u64],
    ) -> &'a [u64] {
        let buffer = &mut buffer[..count];
        if count == 0 {
            return buffer;
        }
        let KmerDirectory::Ranks(offsets) = &lookup.directory else {
            unreachable!("exact prefix rank directory")
        };
        let mut prefix = self.prefix_at(lookup, first);
        let mut end = offsets.get(prefix + 1);
        for (local, output) in buffer.iter_mut().enumerate() {
            let index = first + local;
            while index >= end {
                prefix += 1;
                end = offsets.get(prefix + 1);
            }
            *output = ((prefix as u64) << lookup.shift) | (self.suffix(index) << lookup.padding);
        }
        buffer
    }

    fn len(&self) -> usize {
        self.count
    }
    fn find(&self, lookup: &KmerLookup, word: u64) -> Option<usize> {
        let KmerDirectory::Ranks(offsets) = &lookup.directory else {
            unreachable!("prefix keys require exact rank ranges")
        };
        let (start, end) = offsets.range((word >> lookup.shift) as usize);
        self.search(word, lookup, start, end)
    }
    fn find_many_count<const N: usize>(
        &self,
        lookup: &KmerLookup,
        words: [u64; N],
        count: usize,
    ) -> [Option<usize>; N] {
        self.find_many_ranks(lookup, words, count)
            .map(|rank| (rank != u32::MAX).then_some(rank as usize))
    }
    fn find_many_ranks<const N: usize>(
        &self,
        lookup: &KmerLookup,
        words: [u64; N],
        count: usize,
    ) -> [u32; N] {
        debug_assert!(count <= N);
        let KmerDirectory::Ranks(offsets) = &lookup.directory else {
            unreachable!("prefix keys require exact rank ranges")
        };
        for &word in &words[..count] {
            offsets.prefetch((word >> lookup.shift) as usize);
        }
        // Neighbor probes often alternate between two exact prefix ranges.
        // Cache only their complete boundaries; every suffix is still searched
        // independently. Narrow bounds cover the directory's complete u32 domain.
        let mut recent = (usize::MAX, (0u32, 0u32));
        let mut other = recent;
        let mut bounds = [(0u32, 0u32); N];
        for i in 0..count {
            let prefix = (words[i] >> lookup.shift) as usize;
            let range = if recent.0 == prefix {
                recent.1
            } else {
                let range = if other.0 == prefix {
                    other.1
                } else {
                    let (start, end) = offsets.range(prefix);
                    (start as u32, end as u32)
                };
                other = recent;
                recent = (prefix, range);
                range
            };
            bounds[i] = range;
        }
        for &(start, end) in &bounds[..count] {
            let (start, end) = (start as usize, end as usize);
            if start < end {
                let at = if end - start > LINEAR_KMER_LOOKUP_MAX {
                    start + (end - start) / 2
                } else {
                    start
                };
                self.prefetch(at);
            }
        }
        let mut found = [u32::MAX; N];
        for i in 0..count {
            found[i] = self
                .search(words[i], lookup, bounds[i].0 as usize, bounds[i].1 as usize)
                .map_or(u32::MAX, |index| index as u32);
        }
        found
    }
}

#[cfg(test)]
mod tests {
    use super::super::RankDirectory;
    use super::*;

    #[test]
    fn prefix_keys_preserve_exact_membership_padding_and_partial_batches() {
        let original = (0..131077u64)
            .map(|i| (i << 12) | (3u64 << 52))
            .collect::<Vec<_>>();
        let mut lookup = KmerLookup::with_threads(&original, 1, 26).unwrap();
        lookup.directory = KmerDirectory::Ranks(RankDirectory::new(
            (0..=(1usize << (64 - lookup.shift)))
                .map(|bucket| {
                    original.partition_point(|&word| (word >> lookup.shift) < bucket as u64) as u32
                })
                .collect(),
            false,
        ));
        for threads in [1, 8] {
            let mut packed = PrefixKeys::for_topology(original.clone(), &lookup, threads).unwrap();
            for (i, &word) in original.iter().enumerate() {
                assert_eq!(packed.suffix(i), i as u64);
                assert_eq!(packed.word(&lookup, i), word);
            }
            for first in (0..original.len()).step_by(1024) {
                let mut buffer = [0; 1024];
                let count = (original.len() - first).min(buffer.len());
                assert_eq!(
                    packed.window(&lookup, first, count, &mut buffer),
                    &original[first..first + count]
                );
            }
            for index in [0usize, 65535, 65536, 104856, 104857, 131076] {
                assert_eq!(packed.find(&lookup, original[index]), Some(index));
            }
            let batch =
                [0usize, 1, 65535, 65536, 131075, 131076, 104856, 104857].map(|i| original[i]);
            for count in 0..=batch.len() {
                let mut expected = [None; 8];
                for i in 0..count {
                    expected[i] = original.binary_search(&batch[i]).ok();
                }
                assert_eq!(packed.find_many_count(&lookup, batch, count), expected);
            }
            packed.retire_word_hints();
            for index in [0usize, 65535, 65536, 104856, 104857, 131076] {
                assert_eq!(packed.find(&lookup, original[index]), Some(index));
            }
        }
        for k in [7usize, 10, 14, 22, 26] {
            for count in [0usize, 1, 2, 7, 8, 9, 15, 257] {
                let padding = 64 - 2 * k;
                let bits = 2 * k - 12;
                let mask = (1u64 << bits) - 1;
                let mut words = (0..count)
                    .map(|i| {
                        let prefix = (i % 7) as u64;
                        let suffix = ((i as u64).wrapping_mul(6364136223846793005) >> 7) & mask;
                        ((prefix << bits) | suffix) << padding
                    })
                    .collect::<Vec<_>>();
                words.sort_unstable();
                words.dedup();
                let original = words.clone();
                let mut lookup = KmerLookup::with_threads(&words, 1, k).unwrap();
                let offsets = (0..=(1usize << (64 - lookup.shift)))
                    .map(|bucket| {
                        words.partition_point(|&word| (word >> lookup.shift) < bucket as u64) as u32
                    })
                    .collect();
                lookup.directory = KmerDirectory::Ranks(RankDirectory::new(offsets, false));
                let packed = PrefixKeys::for_topology(words, &lookup, 1).unwrap();
                assert_eq!(packed.len(), original.len());
                assert_eq!(
                    packed.words.len(),
                    if original.is_empty() {
                        0
                    } else {
                        (original.len() * 5 + 3).div_ceil(8)
                    }
                );
                for (index, &word) in original.iter().enumerate() {
                    assert_eq!(packed.suffix(index), (word >> padding) & mask);
                    assert_eq!(packed.word(&lookup, index), word);
                }
                let mut buffer = vec![0; original.len()];
                assert_eq!(
                    packed.window(&lookup, 0, original.len(), &mut buffer),
                    original.as_slice()
                );
                let mut queries = original.clone();
                queries.extend(original.iter().map(|&word| word | 1));
                queries.extend(original.iter().map(|&word| word ^ (1 << padding)));
                queries.extend([0, u64::MAX, 8u64 << 52]);
                // Alternate distant prefixes, reuse both cached ranges and
                // evict them with the following pair across partial batches.
                queries
                    .extend((0..original.len()).flat_map(|i| {
                        [original[i], original[original.len() - 1 - i], original[i]]
                    }));
                for group in queries.chunks(4) {
                    let mut batch = [0; 4];
                    batch[..group.len()].copy_from_slice(group);
                    for used in 0..=group.len() {
                        let mut expected = [None; 4];
                        for at in 0..used {
                            expected[at] = original.binary_search(&batch[at]).ok();
                        }
                        assert_eq!(
                            packed.find_many_count(&lookup, batch, used),
                            expected,
                            "k={k}, count={count}, used={used}"
                        );
                    }
                    for &word in group {
                        assert_eq!(
                            packed.find(&lookup, word),
                            original.binary_search(&word).ok()
                        );
                    }
                }
            }
        }
    }
}
