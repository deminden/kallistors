//! Bounded conservative word lookup with exact ordered candidate verification.
use super::BifrostIndex;
use crate::pseudoalign::{PseudoalignOptions, core::fast_path_enabled};
mod directory;
use directory::WordDirectory;
use rayon::prelude::*;

const MINIMUM_BUCKET: usize = 8;
// Share conservative word selectors across ordered postings. The transient
// directory has a fixed payload ceiling; small indexes allocate only the
// entries they can use. Quantification and BUS share the exact verifier.
const MAX_BYTES: usize = 48 * 1024 * 1024;
// At most nine word/ordinal pairs per posting. Charge temporary vectors,
// both directory representations during conversion, span/ordinal storage,
// geometric capacity rounding and a possible shrink copy conservatively.
// This bound also governs parallel microbatches and oversized fallbacks.
const CONSTRUCTION_BYTES_PER_POSTING: usize = 768;
const TAG_MASK: u64 = (1 << 24) - 1;
const OFFSET_MASK: u64 = (1 << 20) - 1;

pub(crate) struct CachedBucket {
    source_addr: usize,
    source_len: usize,
    // Tags select a conservative, ordered superset. Every possible hit still
    // passes exact reference-word verification before it can affect mapping.
    // One cache-line stream carries the tag and both posting offsets. Tags
    // remain conservative; the original full-word verifier decides every hit.
    directory: WordDirectory,
    ordinals: Box<[u16]>,
}

struct WordSpan {
    word: u32,
    first: u32,
    end: u32,
}

impl BifrostIndex {
    /// Prepare a bounded, transient accelerator for supported read workloads.
    /// The serialized index and exact mapping semantics are unchanged.
    #[doc(hidden)]
    pub fn prepare_single_cell_candidate_cache(&mut self) {
        self.prepare_candidate_selectors(1);
    }

    /// Prepare the same bounded selectors with the requested worker limit.
    #[doc(hidden)]
    pub fn prepare_single_cell_candidate_cache_with_threads(&mut self, threads: usize) {
        self.prepare_candidate_selectors(threads);
    }

    /// Prepare bounded exact selectors when the requested mapping can use them.
    /// Bulk and BUS callers share the scanner's eligibility policy. Diagnostic
    /// callers that require the general scanner should leave the index unprepared.
    #[doc(hidden)]
    pub fn prepare_read_mapping(&mut self, options: PseudoalignOptions, threads: usize) {
        if fast_path_enabled(self, false, options) {
            self.prepare_candidate_selectors(threads);
        }
    }

    fn prepare_candidate_selectors(&mut self, threads: usize) {
        if self.k == 31 && self.g == 23 {
            self.encoded_unitigs.unitig_bases.inline_short_words();
        }
        // Huge nascent references have many more repeated minimizer buckets.
        // Keep the ordinary ceiling for smaller references; all entries remain
        // conservative selectors followed by the same exact verifier.
        let budget = if self.mphf.size() >= (1 << 28) {
            256 * 1024 * 1024
        } else {
            MAX_BYTES
        };
        self.prepare_candidate_cache_threaded(MINIMUM_BUCKET, budget, threads);
    }

    #[cfg(test)]
    fn prepare_candidate_cache(&mut self, minimum: usize, budget: usize) {
        self.prepare_candidate_cache_threaded(minimum, budget, 1);
    }

    fn prepare_candidate_cache_threaded(&mut self, minimum: usize, budget: usize, threads: usize) {
        if self.k != 31 || self.g != 23 || !self.large_candidate_cache.is_empty() {
            return;
        }
        let mut stored = Vec::new();
        let mut bytes = 0usize;
        let mut candidates = Vec::new();
        self.minz_positions
            .for_each_large_bucket(minimum, |bucket, positions| {
                candidates.push((bucket, positions.len()));
            });
        // Admit the largest buckets first so a bounded directory retains its
        // coverage when the reference also contains many medium buckets.
        candidates.sort_unstable_by_key(|(_, positions)| std::cmp::Reverse(*positions));
        // Charge the exact prospective sparse-directory size on each
        // admission. Eligible but uncached buckets require no address slots.
        // The serial and parallel paths use exactly the same builder and
        // largest-first admission. Speculative work is bounded independently
        // of the retained payload; results are always admitted in input order.
        let pool = (threads > 1)
            .then(|| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .ok()
            })
            .flatten();
        let mut begin = 0;
        while begin < candidates.len() {
            let Some(available) = candidate_available(stored.len(), bytes, budget) else {
                break;
            };
            // All remaining candidates have at least `minimum` postings.
            // This is the exact lower bound checked by the existing builder;
            // available capacity cannot increase after this point. Finish the
            // current result batch normally, but never build a futile next one.
            if available < minimum_candidate_bytes(minimum) {
                break;
            }
            let mut end = begin + 1;
            if pool.is_some() {
                let mut scratch = candidates[begin]
                    .1
                    .saturating_mul(CONSTRUCTION_BYTES_PER_POSTING);
                while end < candidates.len() && end - begin < 256 {
                    let next = candidates[end]
                        .1
                        .saturating_mul(CONSTRUCTION_BYTES_PER_POSTING);
                    if scratch.saturating_add(next) > 8 * 1024 * 1024 {
                        break;
                    }
                    scratch += next;
                    end += 1;
                }
            }
            let work = &candidates[begin..end];
            let results: Vec<_> = if let Some(pool) = &pool {
                pool.install(|| {
                    work.par_iter()
                        .map(|&(bucket, _)| self.build_candidate_bucket(bucket, available))
                        .collect()
                })
            } else {
                work.iter()
                    .map(|&(bucket, _)| self.build_candidate_bucket(bucket, available))
                    .collect()
            };
            for result in results {
                let Some(available) = candidate_available(stored.len(), bytes, budget) else {
                    break;
                };
                if let Some((bucket, size)) = result
                    && size <= available
                {
                    bytes += size;
                    stored.push(bucket);
                }
            }
            begin = end;
        }
        drop(pool);
        crate::index::scratch_pages::retire(candidates);
        stored.sort_unstable_by_key(|bucket| bucket.source_addr);
        if stored.is_empty() || stored.len() >= u32::MAX as usize {
            return;
        }
        let Some(slot_count) = candidate_directory_slots(stored.len()) else {
            return;
        };
        let mut slots = vec![0u32; slot_count];
        for (ordinal, bucket) in stored.iter().enumerate() {
            let mut slot = bucket_address_hash(bucket.source_addr) & (slot_count - 1);
            while slots[slot] != 0 {
                slot = (slot + 1) & (slot_count - 1);
            }
            slots[slot] = ordinal as u32 + 1;
        }
        self.candidate_cache_slots = slots;
        self.large_candidate_cache = stored;
    }

    fn build_candidate_bucket(
        &self,
        bucket: usize,
        available: usize,
    ) -> Option<(CachedBucket, usize)> {
        let positions = self.minz_positions.get(bucket);
        // Bound temporary construction independently of retained data.
        // Oversized or mixed buckets always keep the ordinary verifier.
        if positions.len() >= u16::MAX as usize
            || positions
                .len()
                .saturating_mul(CONSTRUCTION_BYTES_PER_POSTING)
                > 8 * 1024 * 1024
            || self.minz_positions.get_with_overcrowded(bucket).1
        {
            return None;
        }
        // Every entry contributes at least one ordered word posting.
        // Avoid construction only when its provable lower bound cannot fit.
        let minimum_size = minimum_candidate_bytes(positions.len());
        if minimum_size > available {
            return None;
        }
        // At most nine distinct windows per posting. For <=4096
        // postings, a 15-bit span index covers every possible multi-list.
        let tiny = (8..16).contains(&positions.len());
        let dense = (129..=4096).contains(&positions.len());
        let mut words = Vec::with_capacity(positions.len() * 9);
        for (ordinal, &position) in positions.iter().enumerate() {
            // Markers and dedicated KM nodes preserve their complete path.
            let (uid, _, km) = crate::index::bifrost::decode_pos_id(position);
            if (position >> 32) as u32 == u32::MAX || km {
                return None;
            }
            let sequence = self.encoded_unitigs.unitig_bases.sequence(uid as usize)?;
            // Every strict start lies within this complete relaxed window
            // when the query minimizer offsets are bounded by k-g. Store
            // conservative word tags; postings keep their
            // original priority and still pass the established verifier.
            let maximum = sequence.len.checked_sub(31)?;
            let (_, relative, _) = crate::index::bifrost::decode_pos_id(position);
            let first = (relative as usize).saturating_sub(8).min(maximum);
            let last = (relative as usize).min(maximum);
            for start in first..=last {
                let code = sequence.kmer_at(start, 31)?;
                let tag = candidate_tag(code);
                let tag = if tiny {
                    tag >> 13
                } else if dense {
                    tag >> 5
                } else {
                    tag
                };
                words.push((tag, ordinal as u32));
            }
        }
        words.sort_unstable();
        words.dedup();
        let mut keys: Vec<WordSpan> = Vec::new();
        let mut ordinals = Vec::with_capacity(words.len());
        for &(word, ordinal) in &words {
            let first = ordinals.len() as u32;
            if let Some(key) = keys.last_mut()
                && key.word == word
            {
                key.end = first + 1;
            } else {
                keys.push(WordSpan {
                    word,
                    first,
                    end: first + 1,
                });
            }
            ordinals.push(ordinal as u16);
        }
        crate::index::scratch_pages::retire(words);
        // A singleton ordinal fits in the existing offset field. Zero
        // count denotes that inline ordinal; real lists contain >=2 IDs.
        // Compact remaining lists in place without changing their order.
        let mut written = 0;
        for key in &mut keys {
            let first = key.first as usize;
            let end = key.end as usize;
            if end - first == 1 {
                key.first = u32::from(ordinals[first]);
                key.end = 0;
            } else {
                ordinals.copy_within(first..end, written);
                key.first = written as u32;
                key.end = (end - first) as u32;
                written += end - first;
            }
        }
        ordinals.truncate(written);
        // Unsupported new small buckets keep their original verifier.
        // All pre-existing16+ bucket layouts retain their original width.
        if tiny && !WordDirectory::tiny_layout(&keys) {
            return None;
        }
        // More precise tags use the descriptor bits guaranteed unused by
        // <=4096-source buckets. An unsupported span shape retains verification.
        if dense && !WordDirectory::dense19_layout(&keys) {
            return None;
        }
        let masked_tiny = tiny && WordDirectory::tiny_masks_supported(&keys, &ordinals);
        let narrow = !tiny && !dense && WordDirectory::narrow_layout(&keys, positions.len());
        let directory_bytes = if masked_tiny {
            WordDirectory::tiny_mask_payload_bytes(&keys)
        } else if tiny {
            WordDirectory::tiny_payload_bytes(&keys)
        } else {
            WordDirectory::payload_bytes(&keys, narrow || dense, dense || !narrow)
        };
        let size = directory_bytes
            + if masked_tiny {
                0
            } else {
                ordinals.len() * std::mem::size_of::<u16>()
            }
            + std::mem::size_of::<CachedBucket>();
        if size > available || ordinals.len() as u64 > OFFSET_MASK {
            return None;
        }
        let directory = if masked_tiny {
            let directory = WordDirectory::build_tiny_masks(&keys, &ordinals)
                .expect("validated unique tiny ordinals");
            crate::index::scratch_pages::retire(ordinals);
            ordinals = Vec::new();
            directory
        } else if tiny {
            WordDirectory::build_tiny(&keys)
        } else if dense {
            WordDirectory::build_dense19(&keys)?
        } else {
            WordDirectory::build(&keys, narrow)
        };
        Some((
            CachedBucket {
                source_addr: positions.as_ptr() as usize,
                source_len: positions.len(),
                directory,
                ordinals: ordinals.into_boxed_slice(),
            },
            size,
        ))
    }

    #[inline(always)]
    pub(crate) fn has_large_candidate_cache(&self) -> bool {
        !self.large_candidate_cache.is_empty()
    }

    #[inline]
    pub(crate) fn large_candidates_for(&self, positions: &[u64]) -> Option<&CachedBucket> {
        if positions.len() < MINIMUM_BUCKET {
            return None;
        }
        // Match immutable allocation identity and length. Hash collisions only
        // advance the directory probe; every slot compares the exact address.
        // Integer addresses are never reconstructed or dereferenced.
        let mask = self.candidate_cache_slots.len().checked_sub(1)?;
        let address = positions.as_ptr() as usize;
        let mut slot = bucket_address_hash(address) & mask;
        loop {
            let ordinal = self.candidate_cache_slots[slot].checked_sub(1)?;
            let bucket = self.large_candidate_cache.get(ordinal as usize)?;
            if bucket.source_addr == address {
                return (bucket.source_len == positions.len()).then_some(bucket);
            }
            slot = (slot + 1) & mask;
        }
    }
}

#[inline]
fn candidate_tag(word: u64) -> u32 {
    // Multiplication before taking the high bits avoids ignoring the first
    // bases of related31-mers. Tags remain only conservative selectors: all
    // collisions retain their complete ordered posting union and verifier.
    let folded = word as u32 ^ (word >> 32) as u32;
    folded.wrapping_mul(0x9e37_79b9) >> 8
}

// Share the builder's own refusal bound with its outer admission loop.
// Saturation also handles unsupported diagnostic thresholds without overflow.
fn minimum_candidate_bytes(postings: usize) -> usize {
    postings
        .saturating_mul(std::mem::size_of::<u16>())
        .saturating_add(std::mem::size_of::<u32>())
        .saturating_add(std::mem::size_of::<CachedBucket>())
}

fn candidate_available(stored: usize, bytes: usize, budget: usize) -> Option<usize> {
    let slots = candidate_directory_slots(stored.checked_add(1)?)?;
    budget
        .checked_sub(bytes)?
        .checked_sub(slots.checked_mul(std::mem::size_of::<u32>())?)
}

fn candidate_directory_slots(buckets: usize) -> Option<usize> {
    buckets.checked_mul(4)?.checked_next_power_of_two()
}

#[inline]
fn bucket_address_hash(address: usize) -> usize {
    let word = (address as u64) >> 3;
    let mixed = (word ^ (word >> 23)).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    (mixed ^ (mixed >> 32)) as usize
}

impl CachedBucket {
    fn word_positions(&self, word: u64) -> CandidateOrdinals<'_> {
        let Some(entry) = self.directory.entry(u64::from(candidate_tag(word))) else {
            return CandidateOrdinals::empty();
        };
        let first = ((entry >> 24) & OFFSET_MASK) as usize;
        let count = (entry >> 44) as usize;
        if count == 0 {
            CandidateOrdinals {
                first: first as u16,
                rest: &[],
            }
        } else {
            CandidateOrdinals::from_slice(&self.ordinals[first..first + count])
        }
    }

    // A full pair of masks is resolved before entering the verifier loop.
    // Other layouts keep their existing concrete list iterator and do not
    // dispatch an extra enum variant for every selected candidate.
    #[inline]
    pub(crate) fn matching_mask_candidates<'a>(
        &'a self,
        positions: &'a [u64],
        forward: u64,
        reverse: u64,
        allow_forward: bool,
        allow_reverse: bool,
    ) -> Option<impl Iterator<Item = u64> + 'a> {
        if !matches!(self.directory, WordDirectory::TinyMasks { .. }) {
            return None;
        }
        let mask = self.directory.tiny_mask(
            u64::from(candidate_tag(forward)),
            u64::from(candidate_tag(reverse)),
            allow_forward,
            allow_reverse,
        )?;
        debug_assert_eq!(positions.len(), self.source_len);
        debug_assert_eq!(positions.as_ptr() as usize, self.source_addr);
        Some(MaskMatchingCandidates {
            entries: positions,
            mask,
        })
    }

    // The existing independent ASCII/posting oracle still exercises every
    // layout. Production dispatch above keeps both iterators monomorphic.
    #[cfg(test)]
    pub(crate) fn matching_candidates<'a>(
        &'a self,
        positions: &'a [u64],
        forward: u64,
        reverse: u64,
        allow_forward: bool,
        allow_reverse: bool,
    ) -> impl Iterator<Item = u64> + 'a {
        if let Some(mask) = self.directory.tiny_mask(
            u64::from(candidate_tag(forward)),
            u64::from(candidate_tag(reverse)),
            allow_forward,
            allow_reverse,
        ) {
            TestCandidates::Mask(MaskMatchingCandidates {
                entries: positions,
                mask,
            })
        } else {
            TestCandidates::Lists(self.matching_candidates_lists(
                positions,
                forward,
                reverse,
                allow_forward,
                allow_reverse,
            ))
        }
    }

    pub(crate) fn matching_candidates_lists<'a>(
        &'a self,
        positions: &'a [u64],
        forward: u64,
        reverse: u64,
        allow_forward: bool,
        allow_reverse: bool,
    ) -> MatchingCandidates<'a> {
        debug_assert_eq!(positions.len(), self.source_len);
        debug_assert_eq!(positions.as_ptr() as usize, self.source_addr);
        MatchingCandidates {
            forward: if allow_forward {
                self.word_positions(forward)
            } else {
                CandidateOrdinals::empty()
            },
            reverse: if allow_reverse {
                self.word_positions(reverse)
            } else {
                CandidateOrdinals::empty()
            },
            entries: positions,
        }
    }
}

struct MaskMatchingCandidates<'a> {
    entries: &'a [u64],
    mask: u16,
}
impl Iterator for MaskMatchingCandidates<'_> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<u64> {
        if self.mask == 0 {
            return None;
        }
        let ordinal = self.mask.trailing_zeros() as usize;
        self.mask &= self.mask - 1;
        self.entries.get(ordinal).copied()
    }
}
#[cfg(test)]
enum TestCandidates<'a> {
    Mask(MaskMatchingCandidates<'a>),
    Lists(MatchingCandidates<'a>),
}
#[cfg(test)]
impl Iterator for TestCandidates<'_> {
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
        match self {
            Self::Mask(iter) => iter.next(),
            Self::Lists(iter) => iter.next(),
        }
    }
}

pub(crate) struct MatchingCandidates<'a> {
    forward: CandidateOrdinals<'a>,
    reverse: CandidateOrdinals<'a>,
    entries: &'a [u64],
}
impl<'a> Iterator for MatchingCandidates<'a> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let a = self.forward.first;
        let b = self.reverse.first;
        let ordinal = a.min(b);
        if ordinal == u16::MAX {
            return None;
        }
        if a == ordinal {
            self.forward.advance();
        }
        if b == ordinal {
            self.reverse.advance();
        }
        self.entries.get(ordinal as usize).copied()
    }
}

// The first ordinal is already in a register for both inline and list cases.
// All remaining IDs borrow the immutable, original-order posting arena.
struct CandidateOrdinals<'a> {
    first: u16,
    rest: &'a [u16],
}
impl<'a> CandidateOrdinals<'a> {
    fn empty() -> Self {
        Self {
            first: u16::MAX,
            rest: &[],
        }
    }
    fn from_slice(values: &'a [u16]) -> Self {
        if let Some((&first, rest)) = values.split_first() {
            Self { first, rest }
        } else {
            Self::empty()
        }
    }
    fn advance(&mut self) {
        *self = Self::from_slice(self.rest);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn bucket_directory_preserves_identity_through_collisions() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("reference.fa");
        let path = dir.path().join("reference.idx");
        std::fs::write(&fasta, ">tx\nCGGGTGTGGTGGCGGGCACCTGTAGTCCCAG\n").unwrap();
        crate::index::build_index(&path, &[fasta], Default::default()).unwrap();
        let mut index = crate::pseudoalign::build_bifrost_index(&path).unwrap();
        // Pigeonhole128 distinct allocations into16 slots to guarantee at
        // least eight colliding addresses, independently of allocator layout.
        let mut groups: [Vec<Vec<u64>>; 16] = std::array::from_fn(|_| Vec::new());
        for i in 0..128u64 {
            let sample = vec![i; 65];
            let slot = super::bucket_address_hash(sample.as_ptr() as usize) & 15;
            groups[slot].push(sample);
        }
        let samples: Vec<_> = groups
            .into_iter()
            .max_by_key(Vec::len)
            .unwrap()
            .into_iter()
            .take(8)
            .collect();
        assert_eq!(samples.len(), 8);
        index.candidate_cache_slots = vec![0; 16];
        index.large_candidate_cache = samples
            .iter()
            .map(|sample| super::CachedBucket {
                source_addr: sample.as_ptr() as usize,
                source_len: sample.len(),
                directory: super::WordDirectory::default(),
                ordinals: Box::default(),
            })
            .collect();
        for (ordinal, bucket) in index.large_candidate_cache.iter().enumerate() {
            let mut slot = super::bucket_address_hash(bucket.source_addr) & 15;
            while index.candidate_cache_slots[slot] != 0 {
                slot = (slot + 1) & 15;
            }
            index.candidate_cache_slots[slot] = ordinal as u32 + 1;
        }
        for sample in &samples {
            assert!(index.large_candidates_for(sample).is_some());
            assert!(index.large_candidates_for(&sample[..64]).is_none());
        }
        assert!(index.large_candidates_for(&[0; 65]).is_none());
        index.large_candidate_cache.clear();
        assert!(index.large_candidates_for(&samples[0]).is_none());
    }

    #[test]
    fn cached_words_preserve_posting_order_and_memory_budget() {
        for len in [128u32, 257] {
            let keys: Vec<_> = (0..len)
                .map(|i| super::WordSpan {
                    word: if i + 1 == len {
                        super::TAG_MASK as u32
                    } else {
                        i * 7
                    },
                    first: if len == 257 && i == 127 {
                        super::OFFSET_MASK as u32
                    } else {
                        i % 128
                    },
                    end: if (len == 128 && i == 127) || (len == 257 && i < 128) {
                        2
                    } else {
                        0
                    },
                })
                .collect();
            assert!(super::WordDirectory::narrow_layout(&keys, 128));
            assert!(!super::WordDirectory::narrow_layout(&keys, 129));
            let narrow = super::WordDirectory::build(&keys, true);
            let wide = super::WordDirectory::build(&keys, false);
            for query in (0..=2000).chain([super::TAG_MASK]) {
                let expected = keys
                    .iter()
                    .find(|key| u64::from(key.word) == query)
                    .map(|key| {
                        u64::from(key.word)
                            | (u64::from(key.first) << 24)
                            | (u64::from(key.end) << 44)
                    });
                assert_eq!(narrow.entry(query), expected);
                assert_eq!(wide.entry(query), expected);
            }
        }
        // Inline and list heads share the same sorted, duplicate-free merge.
        let positions: Vec<_> = (0..u16::MAX).map(u64::from).collect();
        for left in [vec![], vec![0], vec![65534], vec![0, 5, 65534]] {
            for right in [vec![], vec![0], vec![65534], vec![1, 5, 65534]] {
                let mut expected = left
                    .iter()
                    .chain(&right)
                    .map(|&i| u64::from(i))
                    .collect::<Vec<_>>();
                expected.sort_unstable();
                expected.dedup();
                let actual: Vec<_> = super::MatchingCandidates {
                    forward: super::CandidateOrdinals::from_slice(&left),
                    reverse: super::CandidateOrdinals::from_slice(&right),
                    entries: &positions,
                }
                .collect();
                assert_eq!(actual, expected);
            }
        }
        for first in [0u64, 1, 65535, super::OFFSET_MASK] {
            for end in [first, super::OFFSET_MASK] {
                let entry = super::TAG_MASK | (first << 24) | (end << 44);
                assert_eq!(entry & super::TAG_MASK, super::TAG_MASK);
                assert_eq!((entry >> 24) & super::OFFSET_MASK, first);
                assert_eq!(entry >> 44, end);
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("reference.fa");
        let path = dir.path().join("reference.idx");
        std::fs::write(&fasta, ">short\nCGGGTGTGGTGGCGGGCACCTGTAGTCCCAG\n>long\nACGTCGATCGTAGCATCGATCGTACGATCGATGCTAGCTAGCTAGC\n").unwrap();
        crate::index::build_index(&path, &[fasta], Default::default()).unwrap();
        let mut index = crate::pseudoalign::build_bifrost_index(&path).unwrap();
        index.prepare_candidate_cache(1, 0);
        assert!(!index.has_large_candidate_cache());
        index.prepare_candidate_cache(1, 1 << 20);
        assert!(index.has_large_candidate_cache());
        let mut groups: Vec<Vec<u64>> = (0..index.minz_positions.len())
            .map(|i| index.minz_positions.get(i).to_vec())
            .collect();
        let ordinary = groups
            .iter()
            .flatten()
            .copied()
            .find(|p| (*p >> 32) as u32 != u32::MAX && *p & 0x8000_0000 == 0)
            .unwrap();
        groups.push(vec![ordinary; 8]);
        groups.push(vec![ordinary; 15]);
        groups.push(vec![ordinary; 257]);
        index.large_candidate_cache.clear();
        index.candidate_cache_slots.clear();
        index.minz_positions =
            super::super::MinzPositionIndex::from_counts(groups.iter().map(Vec::len).collect());
        index
            .minz_positions
            .values
            .copy_from_slice(&groups.concat());
        index.prepare_candidate_cache(1, 1 << 20);
        assert!(index.large_candidate_cache.iter().any(|b| matches!(
            b.directory,
            super::WordDirectory::Dense19 { .. } | super::WordDirectory::HashedDense19 { .. }
        )));
        assert!(index.large_candidate_cache.iter().any(|b| matches!(
            b.directory,
            super::WordDirectory::Tiny11 { .. } | super::WordDirectory::TinyMasks { .. }
        )));
        fn snapshot(index: &super::BifrostIndex) -> Vec<(usize, usize, Vec<u16>, String)> {
            index
                .large_candidate_cache
                .iter()
                .map(|bucket| {
                    (
                        bucket.source_addr,
                        bucket.source_len,
                        bucket.ordinals.to_vec(),
                        format!("{:?}", bucket.directory),
                    )
                })
                .collect()
        }
        for budget in [0, 128, 512, 4096, 1 << 20] {
            index.large_candidate_cache.clear();
            index.candidate_cache_slots.clear();
            index.prepare_candidate_cache_threaded(1, budget, 1);
            let serial = snapshot(&index);
            // Check charged payload against actual boxed-array lengths,
            // independently of the prospective constructor accounting.
            let charged = index.candidate_cache_slots.len() * std::mem::size_of::<u32>()
                + index
                    .large_candidate_cache
                    .iter()
                    .map(|bucket| {
                        std::mem::size_of::<super::CachedBucket>()
                            + bucket.directory.retained_payload_bytes()
                            + bucket.ordinals.len() * std::mem::size_of::<u16>()
                    })
                    .sum::<usize>();
            assert!(charged <= budget, "charged payload {charged} > {budget}");
            let serial_slots = index.candidate_cache_slots.clone();
            index.large_candidate_cache.clear();
            index.candidate_cache_slots.clear();
            index.prepare_candidate_cache_threaded(1, budget, 4);
            assert_eq!(snapshot(&index), serial);
            assert_eq!(index.candidate_cache_slots, serial_slots);
        }
        let mut short = 0;
        let mut long = 0;
        for bucket in &index.large_candidate_cache {
            let positions = (0..index.minz_positions.len())
                .map(|i| index.minz_positions.get(i))
                .find(|p| p.as_ptr() as usize == bucket.source_addr)
                .unwrap();
            assert_eq!(positions.len(), bucket.source_len);
            for &position in positions {
                let (uid, _, km) = crate::index::bifrost::decode_pos_id(position);
                assert!(!km);
                let sequence = index.unitig_sequence(uid as usize).unwrap();
                if sequence.len() == 31 {
                    short += 1;
                } else {
                    long += 1;
                }
            }
        }
        for bucket in &index.large_candidate_cache {
            // Derive queries from independent ASCII reference windows, not the
            // directory being tested, so omitted keys are covered as well.
            let mut queries = std::collections::BTreeSet::new();
            let positions = (0..index.minz_positions.len())
                .map(|i| index.minz_positions.get(i))
                .find(|p| p.as_ptr() as usize == bucket.source_addr)
                .unwrap();
            for &position in positions {
                let (uid, _, _) = crate::index::bifrost::decode_pos_id(position);
                let sequence = index.unitig_sequence(uid as usize).unwrap();
                for window in sequence.windows(31) {
                    let word = crate::dna::encode_kmer_pair(window).unwrap().0;
                    // Include a guaranteed collision, which must retain the
                    // original candidate for the exact downstream verifier.
                    // Exercise every bit, including the high folded byte
                    // that a low24-bit truncation could not distinguish.
                    queries.extend((0..62).map(|bit| word ^ (1u64 << bit)));
                    // Invert the odd multiplicative mix to force a collision
                    // only in the 16-bit selector, with unequal full24 tags.
                    let folded = word as u32 ^ (word >> 32) as u32;
                    let other =
                        (folded.wrapping_mul(0x9e37_79b9) ^ 0x100).wrapping_mul(0x144c_bc89);
                    let collision = (word & !0xffff_ffff) | u64::from(other ^ (word >> 32) as u32);
                    assert_ne!(super::candidate_tag(word), super::candidate_tag(collision));
                    assert_eq!(
                        super::candidate_tag(word) >> 8,
                        super::candidate_tag(collision) >> 8
                    );
                    queries.insert(collision);
                    let other =
                        (folded.wrapping_mul(0x9e37_79b9) ^ 0x1000).wrapping_mul(0x144c_bc89);
                    let collision = (word & !0xffff_ffff) | u64::from(other ^ (word >> 32) as u32);
                    assert_ne!(super::candidate_tag(word), super::candidate_tag(collision));
                    assert_eq!(
                        super::candidate_tag(word) >> 5,
                        super::candidate_tag(collision) >> 5
                    );
                    queries.insert(collision);
                    let other =
                        (folded.wrapping_mul(0x9e37_79b9) ^ 0x10000).wrapping_mul(0x144c_bc89);
                    let collision = (word & !0xffff_ffff) | u64::from(other ^ (word >> 32) as u32);
                    assert_ne!(super::candidate_tag(word), super::candidate_tag(collision));
                    assert_eq!(
                        super::candidate_tag(word) >> 13,
                        super::candidate_tag(collision) >> 13
                    );
                    queries.insert(collision);
                    queries.extend([
                        word,
                        word ^ 1,
                        word ^ (1u64 << 60),
                        word ^ 1 ^ (1u64 << 32),
                        word ^ (1u64 << 24),
                    ]);
                }
            }
            let selector = |word| {
                let tag = super::candidate_tag(word);
                if matches!(
                    bucket.directory,
                    super::WordDirectory::Tiny11 { .. } | super::WordDirectory::TinyMasks { .. }
                ) {
                    tag >> 13
                } else if matches!(
                    bucket.directory,
                    super::WordDirectory::Dense19 { .. }
                        | super::WordDirectory::HashedDense19 { .. }
                ) {
                    tag >> 5
                } else {
                    tag
                }
            };
            for query in queries {
                for flags in 0..4 {
                    let reverse = crate::dna::reverse_complement_word(query, 31);
                    let expected: Vec<_> = positions
                        .iter()
                        .copied()
                        .filter(|&position| {
                            let (uid, relative, _) = crate::index::bifrost::decode_pos_id(position);
                            let sequence = index.unitig_sequence(uid as usize).unwrap();
                            let maximum = sequence.len() - 31;
                            let first = (relative as usize).saturating_sub(8).min(maximum);
                            let last = (relative as usize).min(maximum);
                            sequence.windows(31).enumerate().any(|(start, word)| {
                                if start < first || start > last {
                                    return false;
                                }
                                let code = crate::dna::encode_kmer_pair(word).unwrap().0;
                                (flags & 1 != 0 && selector(code) == selector(query))
                                    || (flags & 2 != 0 && selector(code) == selector(reverse))
                            })
                        })
                        .collect();
                    let actual: Vec<_> = bucket
                        .matching_candidates(
                            positions,
                            query,
                            reverse,
                            flags & 1 != 0,
                            flags & 2 != 0,
                        )
                        .collect();
                    assert_eq!(actual, expected);
                }
            }
        }
        assert!(short > 0 && long > 0);
        let ordinary = (0..index.minz_positions.len())
            .map(|i| index.minz_positions.get(i))
            .find(|p| p.as_ptr() as usize == index.large_candidate_cache[0].source_addr)
            .unwrap()[0];
        index.large_candidate_cache.clear();
        index.minz_positions = super::super::MinzPositionIndex::from_counts(vec![2, 2, 1]);
        let invalid = ((index.unitig_count() as u64 + 1) << 32) | 1;
        index.minz_positions.values.copy_from_slice(&[
            ordinary,
            1u64 << 31,
            ordinary,
            (u64::from(u32::MAX) << 32) | 0x8000_0000,
            invalid,
        ]);
        index.prepare_candidate_cache(1, 1 << 20);
        assert!(!index.has_large_candidate_cache());
    }
}
