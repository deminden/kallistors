//! Temporary, conservative membership rejection for large graph construction.
use super::{IndexWorkers, KmerLookup, MappingKeys, prefetch_rank_delta};
use crate::{Error, Result};
use rayon::prelude::*;

struct Block {
    words: Box<[u64]>,
}
pub(super) struct NegativeFilter {
    blocks: Vec<Block>,
    shift: usize,
}

impl NegativeFilter {
    #[cfg(test)]
    pub(super) fn build(entries: &[u64], k: usize, threads: usize) -> Result<Self> {
        let lookup = KmerLookup::with_threads(entries, threads, k)?;
        Self::build_mapped(entries, &lookup, k, threads)
    }

    pub(super) fn build_mapped(
        entries: &(impl MappingKeys + ?Sized),
        lookup: &KmerLookup,
        k: usize,
        threads: usize,
    ) -> Result<Self> {
        let bits = 2 * k.min(6);
        let shift = 64 - bits;
        let boundaries = (0..=(1usize << bits))
            .map(|prefix| {
                let mut first = 0;
                let mut end = entries.len();
                while first < end {
                    let mid = first + (end - first) / 2;
                    if entries.word(lookup, mid) >> shift < prefix as u64 {
                        first = mid + 1;
                    } else {
                        end = mid;
                    }
                }
                first
            })
            .collect::<Vec<_>>();
        let mut blocks = boundaries
            .windows(2)
            .map(|bounds| {
                let count = bounds[1] - bounds[0];
                if count == 0 {
                    return Ok(Block {
                        // A zero word makes empty prefixes exact misses too,
                        // without a conditional in the query loop.
                        words: vec![0].into_boxed_slice(),
                    });
                }
                let bits = count
                    .checked_mul(12)
                    .and_then(usize::checked_next_power_of_two)
                    .ok_or_else(|| Error::InvalidFormat("membership filter size overflow".into()))?
                    .max(64);
                Ok(Block {
                    words: vec![0; bits / 64].into_boxed_slice(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let fill = |(prefix, block): (usize, &mut Block)| {
            let end = boundaries[prefix + 1];
            let mask = block.words.len() * 64 - 1;
            let mut buffer = [0u64; 1024];
            for first in (boundaries[prefix]..end).step_by(buffer.len()) {
                let count = (end - first).min(buffer.len());
                for &word in entries.window(lookup, first, count, &mut buffer) {
                    let hash = Self::hash(word);
                    let bit = hash as usize & mask;
                    block.words[bit >> 6] |= Self::probe_mask(hash);
                }
            }
        };
        if threads > 1 {
            IndexWorkers::new(threads)?
                .install(|| blocks.par_iter_mut().enumerate().for_each(fill));
        } else {
            blocks.iter_mut().enumerate().for_each(fill);
        }
        Ok(Self { blocks, shift })
    }

    #[inline]
    fn hash(word: u64) -> u64 {
        let folded = (word ^ (word >> 33)).wrapping_mul(0xff51_afd7_ed55_8ccd);
        folded ^ (folded >> 33)
    }

    #[inline]
    fn probe_mask(hash: u64) -> u64 {
        // Three hash-derived bit positions share one loaded word. Requiring
        // all of them reduces false admissions without another random read
        // or a larger filter; every admitted key is still checked exactly.
        (1 << (hash & 63)) | (1 << ((hash >> 32) & 63)) | (1 << ((hash >> 48) & 63))
    }

    // False is an exact absence proof. Every true result still needs the
    // ordinary full-key lookup; no filter bit is ever accepted as membership.
    pub(super) fn admit<const N: usize>(&self, words: &[u64; N], count: usize) -> [bool; N] {
        debug_assert!(count <= N);
        let mut probes = [(&self.blocks[0].words[0], 0u64); N];
        for i in 0..count {
            let block = &self.blocks[(words[i] >> self.shift) as usize];
            let hash = Self::hash(words[i]);
            // Each block has a positive power-of-two word count. Deriving the
            // bit mask from its slice length halves the hot metadata while
            // retaining independently allocated, reclaimable blocks.
            let bit = hash as usize & (block.words.len() * 64 - 1);
            let slot = &block.words[bit >> 6];
            prefetch_rank_delta(slot);
            probes[i] = (slot, Self::probe_mask(hash));
        }
        let mut admitted = [false; N];
        for (i, (word, mask)) in probes.into_iter().enumerate().take(count) {
            admitted[i] = *word & mask == mask;
        }
        admitted
    }

    pub(super) fn retire(self) {
        for block in self.blocks {
            super::super::scratch_pages::retire(block.words.into_vec());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_membership_survives_prefixes_empty_blocks_and_partial_batches() {
        for k in [3, 15, 31] {
            for count in [0, 1, 3, 129, 10001] {
                let mut random = 19u64;
                let mut entries = (0..count)
                    .map(|_| {
                        random ^= random << 13;
                        random ^= random >> 7;
                        random ^= random << 17;
                        random & super::super::kmer_mask(k)
                    })
                    .collect::<Vec<_>>();
                if count > 0 {
                    entries.extend([0, super::super::kmer_mask(k)]);
                }
                entries.sort_unstable();
                entries.dedup();
                for threads in [1, 8] {
                    let filter = NegativeFilter::build(&entries, k, threads).unwrap();
                    for group in entries.chunks(8) {
                        let mut words = [0; 8];
                        words[..group.len()].copy_from_slice(group);
                        for used in 0..=group.len() {
                            let admitted = filter.admit(&words, used);
                            assert!(admitted[..used].iter().all(|&present| present));
                            assert!(admitted[used..].iter().all(|&present| !present));
                        }
                    }
                    if entries.is_empty() {
                        assert_eq!(filter.admit(&[0, u64::MAX], 2), [false; 2]);
                    }
                    filter.retire();
                }
            }
        }
    }
}
