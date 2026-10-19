//! Exact landmark locations for reference-transcript unitig walks.
use super::{
    AnchorMap, IndexWorkers, KmerLoc, PackedSequence, ReferenceKmerLocations, UnitigPath,
    concatenate_shards, kmer_mask, reverse_complement_word, scan_kmer_pairs,
};
use crate::{Error, Result};
use rayon::prelude::*;
use rdst::RadixSort;

const MAX_PREFIX_BITS: usize = 24;
const EMPTY: u64 = u64::MAX;

#[derive(Clone, Copy)]
struct Record {
    word: u64,
    location: KmerLoc,
}
impl rdst::RadixKey for Record {
    const LEVELS: usize = 8;
    fn get_level(&self, level: usize) -> u8 {
        (self.word >> (level * 8)) as u8
    }
}

pub(super) struct SparseLocations {
    blocks: Vec<Box<[Record]>>,
    shift: usize,
}

pub(super) fn worthwhile(kmers: usize, paths: usize, transcripts: usize) -> bool {
    let upper_records = 2 * paths as u64 + transcripts as u64;
    let bytes =
        upper_records * std::mem::size_of::<Record>() as u64 + ((1u64 << MAX_PREFIX_BITS) + 1) * 4;
    // Keep small references unchanged. Reserve four copies for construction
    // scratch against the dense four-byte table; endpoint-dense graphs retain
    // their original representation.
    kmers > 1 << 28 && 4 * bytes < 4 * kmers as u64
}

impl SparseLocations {
    // A reference walk starts at its transcript's first k-mer. Each exact
    // unitig jump then ends the transcript or enters an oriented endpoint of
    // another unitig. These landmarks cover every location the walk queries;
    // transcript starts inside unitigs must be included as well.
    pub(super) fn build(
        paths: &[UnitigPath],
        transcripts: &[PackedSequence],
        k: usize,
        threads: usize,
    ) -> Result<Self> {
        let mask = kmer_mask(k);
        let mut starts = AnchorMap::with_capacity_and_hasher(transcripts.len(), Default::default());
        for transcript in transcripts
            .iter()
            .filter(|transcript| transcript.len() >= k)
        {
            let word = transcript.word(0, k);
            starts.insert(word.min(reverse_complement_word(word, k, mask)), ());
        }
        let collect = |paths: &[UnitigPath], first: usize| -> Result<Vec<Record>> {
            let mut records = Vec::with_capacity(paths.len() * 2);
            for (relative, path) in paths.iter().enumerate() {
                let unitig = first + relative;
                let last = path.seq.len().checked_sub(k).ok_or_else(|| {
                    Error::InvalidFormat("joined unitig is shorter than k".into())
                })?;
                KmerLoc::new(unitig, last, true)?;
                scan_kmer_pairs(&path.seq, k, |offset, forward, reverse| {
                    let canonical = forward.min(reverse);
                    if offset == 0 || offset == last || starts.contains_key(&canonical) {
                        records.push(Record {
                            word: canonical,
                            location: KmerLoc::new(unitig, offset, forward <= reverse)?,
                        });
                    }
                    Ok(())
                })?;
            }
            Ok(records)
        };
        let workers = IndexWorkers::new(threads)?;
        let shards = if threads > 1 {
            workers.install(|| {
                paths
                    .par_chunks(64)
                    .enumerate()
                    .map(|(group, paths)| collect(paths, group * 64))
                    .collect::<Result<Vec<_>>>()
            })?
        } else {
            paths
                .chunks(64)
                .enumerate()
                .map(|(group, paths)| collect(paths, group * 64))
                .collect::<Result<Vec<_>>>()?
        };
        drop(starts);
        let mut records = concatenate_shards(shards);
        workers.collect_unused(&records);
        if threads > 1 {
            workers.install(|| records.radix_sort_unstable());
        } else {
            records.sort_unstable_by_key(|record| record.word);
        }
        let written = records.len();
        let mut conflict = false;
        records.dedup_by(|a, b| {
            if a.word == b.word {
                conflict |= a.location.0 != b.location.0;
                true
            } else {
                false
            }
        });
        if conflict || records.len() > u32::MAX as usize {
            return Err(Error::InvalidFormat(
                "inconsistent sparse unitig locations".into(),
            ));
        }
        super::super::scratch_pages::discard_unused(&mut records, written);
        // Sorted full words give each prefix one exclusive table owner.
        // Power-of-two rounding at <=3/4 load needs fewer than three table
        // records per landmark; the construction policy reserves four copies.
        let target = records.len().div_ceil(256).max(2);
        let bits = (usize::BITS as usize - target.saturating_sub(1).leading_zeros() as usize)
            .clamp(1, (2 * k).min(16));
        let shift = 64 - bits;
        let bounds = (0..=(1usize << bits))
            .map(|prefix| records.partition_point(|record| record.word >> shift < prefix as u64))
            .collect::<Vec<_>>();
        let build = |bounds: &[usize]| -> Result<Box<[Record]>> {
            let source = &records[bounds[0]..bounds[1]];
            if source.is_empty() {
                return Ok(Box::default());
            }
            let capacity = source
                .len()
                .checked_mul(4)
                .map(|n| n.div_ceil(3))
                .and_then(usize::checked_next_power_of_two)
                .ok_or_else(|| Error::InvalidFormat("landmark table size overflow".into()))?;
            let mut table = vec![
                Record {
                    word: EMPTY,
                    location: KmerLoc::default()
                };
                capacity
            ];
            for &record in source {
                debug_assert_ne!(record.word, EMPTY);
                let mut slot = hash(record.word) & (capacity - 1);
                while table[slot].word != EMPTY {
                    slot = (slot + 1) & (capacity - 1);
                }
                table[slot] = record;
            }
            Ok(table.into_boxed_slice())
        };
        let blocks = if threads > 1 {
            workers.install(|| bounds.par_windows(2).map(build).collect::<Result<Vec<_>>>())?
        } else {
            bounds.windows(2).map(build).collect::<Result<Vec<_>>>()?
        };
        super::super::scratch_pages::retire(records);
        workers.collect_unused(&blocks);
        Ok(Self { blocks, shift })
    }
}

#[inline]
fn hash(word: u64) -> usize {
    let folded = (word ^ (word >> 33)).wrapping_mul(0xff51_afd7_ed55_8ccd);
    (folded ^ (folded >> 33)) as usize
}

impl ReferenceKmerLocations for SparseLocations {
    #[inline]
    fn find(&self, word: u64) -> Option<KmerLoc> {
        let table = &self.blocks[(word >> self.shift) as usize];
        if table.is_empty() {
            return None;
        }
        let mask = table.len() - 1;
        let mut slot = hash(word) & mask;
        loop {
            let record = table[slot];
            if record.word == EMPTY {
                return None;
            }
            if record.word == word {
                return Some(record.location);
            }
            slot = (slot + 1) & mask;
        }
    }
}
