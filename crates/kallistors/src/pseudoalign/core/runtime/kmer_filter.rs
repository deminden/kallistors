//! A no-false-negative membership filter over every k-mer of the graph.
//!
//! Every successful fast-path match requires the read k-mer (or its reverse
//! complement) to equal a graph k-mer: a unitig window, a KM node or an
//! abundant k-mer. A canonical k-mer the filter rejects therefore cannot
//! match, so the probe returns `None` without the minimizer, MPHF, bucket and
//! unitig loads of a full lookup. Most such probes are k-mers that contain a
//! sequencing error. A false positive only falls through to the full lookup.

use crate::pseudoalign::BifrostIndex;

const BITS_PER_KEY: usize = 12;
const PROBES: u32 = 6;

#[derive(Clone, Copy, Default)]
#[repr(C, align(64))]
struct Block([u64; 8]);

/// Partitions of the filter written by one thread each during the build;
/// each is a contiguous slice small enough to stay in L2.
const PARTITIONS: usize = 1024;
/// Unitig rounds, bounding the partition buffers to about a quarter of the
/// k-mers at a time.
const ROUNDS: usize = 4;

pub(crate) struct KmerFilter {
    blocks: Box<[Block]>,
    per_partition: usize,
    k: usize,
}

#[inline(always)]
fn hashes(canonical: u64) -> (u64, u64) {
    let mut h = canonical.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    h ^= h >> 29;
    let g = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    (h, g ^ (g >> 31))
}

/// High hash bits choose the partition, low bits the block inside it, so the
/// build can split work by partition without any division.
#[inline(always)]
fn partition_of(h: u64) -> usize {
    (((h >> 32) * PARTITIONS as u64) >> 32) as usize
}

#[inline(always)]
fn filter_slot(h: u64, per_partition: usize) -> usize {
    partition_of(h) * per_partition + (((h & 0xffff_ffff) * per_partition as u64) >> 32) as usize
}

impl KmerFilter {
    /// Build from all k-mers of `index` with `threads` workers. `None` when
    /// the graph uses features this filter does not enumerate.
    pub(crate) fn build(index: &BifrostIndex, threads: usize) -> Option<Self> {
        if index.use_shade {
            return None;
        }
        let abundant: Vec<u64> = index
            .h_kmer_map
            .as_ref()
            .map_or_else(Vec::new, |map| map.keys().copied().collect());
        Self::build_from(&index.encoded_unitigs, &abundant, index.k, threads)
    }

    /// Build from the graph's unitig and KM sequences plus the canonical
    /// abundant k-mers, before the rest of the index exists.
    pub(crate) fn build_from(
        encoded: &crate::pseudoalign::bifrost::EncodedUnitigs,
        abundant: &[u64],
        k: usize,
        threads: usize,
    ) -> Option<Self> {
        if !(1..=32).contains(&k) {
            return None;
        }
        let unitigs = &encoded.unitig_bases;
        let km = &encoded.km_bases;
        let count = unitigs.count();
        let km_count = km.count();
        let mut windows = 0usize;
        for id in 0..count {
            windows += unitigs.len(id)?.saturating_sub(k - 1);
        }
        windows += km_count + abundant.len();
        let per_partition = (windows * BITS_PER_KEY).div_ceil(512 * PARTITIONS).max(1);
        let blocks = per_partition * PARTITIONS;
        let mut filter = Self {
            blocks: vec![Block::default(); blocks].into_boxed_slice(),
            per_partition,
            k,
        };
        let threads = threads.max(1);
        let mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let mut extra: Vec<u64> = (0..km_count)
            .filter_map(|id| km.encode_kmer_at(id, 0, k))
            .map(|word| word.min(crate::dna::reverse_complement_word(word, k)))
            .collect();
        extra.extend_from_slice(abundant);
        let round_len = count.div_ceil(ROUNDS).max(1);
        // Per-worker partition buckets, reused by every round so their pages
        // are faulted in once. Sized for an even share of one round's k-mers.
        let reserve = (windows / ROUNDS / threads / PARTITIONS) * 5 / 4 + 16;
        let mut buckets: Vec<Vec<Vec<u64>>> = (0..threads)
            .map(|_| {
                (0..PARTITIONS)
                    .map(|_| Vec::with_capacity(reserve))
                    .collect()
            })
            .collect();
        for round_start in (0..count).step_by(round_len) {
            let round_end = (round_start + round_len).min(count);
            let step = (round_end - round_start).div_ceil(threads).max(1);
            let filter_ref = &filter;
            let extra = (round_start == 0).then_some(extra.as_slice());
            // Phase 1: every worker sorts its unitigs' canonical k-mers into
            // partition buckets with sequential appends.
            std::thread::scope(|scope| {
                for (worker, buckets) in buckets.iter_mut().enumerate() {
                    let first = round_start + worker * step;
                    scope.spawn(move || {
                        for bucket in buckets.iter_mut() {
                            bucket.clear();
                        }
                        let mut push = |canonical: u64| {
                            buckets[filter_ref.partition(canonical)].push(canonical);
                        };
                        if worker == 0 {
                            for &canonical in extra.unwrap_or(&[]) {
                                push(canonical);
                            }
                        }
                        for id in first..(first + step).min(round_end) {
                            let Some(sequence) = unitigs.sequence(id) else {
                                continue;
                            };
                            if sequence.len < k {
                                continue;
                            }
                            let (mut forward, mut reverse) = (0u64, 0u64);
                            let mut seen = 0usize;
                            // Bases are packed 32 per word, first base in the top bits.
                            for (index, &word) in sequence.words().iter().enumerate() {
                                let bases = (sequence.len - index * 32).min(32);
                                for offset in 0..bases {
                                    let base = (word >> (62 - 2 * offset)) & 3;
                                    forward = ((forward << 2) | base) & mask;
                                    reverse = (reverse >> 2) | ((3 - base) << (2 * (k - 1)));
                                    seen += 1;
                                    if seen >= k {
                                        push(forward.min(reverse));
                                    }
                                }
                            }
                        }
                    });
                }
            });
            // Phase 2: each partition's slice is written by exactly one thread.
            let slices: Vec<&mut [Block]> = filter.blocks.chunks_mut(per_partition).collect();
            let mut by_thread: Vec<Vec<(usize, &mut [Block])>> =
                (0..threads).map(|_| Vec::new()).collect();
            for (partition, slice) in slices.into_iter().enumerate() {
                by_thread[partition % threads].push((partition, slice));
            }
            let buckets = &buckets;
            std::thread::scope(|scope| {
                for work in by_thread {
                    scope.spawn(move || {
                        for (partition, slice) in work {
                            let base = partition * per_partition;
                            for worker in buckets {
                                for &canonical in &worker[partition] {
                                    let (h, g) = hashes(canonical);
                                    let block = &mut slice[filter_slot(h, per_partition) - base];
                                    for probe in 0..PROBES {
                                        let bit = (g >> (9 * probe)) & 511;
                                        block.0[(bit >> 6) as usize] |= 1 << (bit & 63);
                                    }
                                }
                            }
                        }
                    });
                }
            });
        }
        let _ = &mut filter;
        Some(filter)
    }

    #[inline(always)]
    fn partition(&self, canonical: u64) -> usize {
        partition_of(hashes(canonical).0)
    }

    #[inline(always)]
    fn block(&self, h: u64) -> &Block {
        &self.blocks[filter_slot(h, self.per_partition)]
    }

    /// `false` proves no graph k-mer equals `canonical` = min(forward, rc).
    #[inline(always)]
    pub(crate) fn may_contain(&self, canonical: u64) -> bool {
        let (h, g) = hashes(canonical);
        let block = self.block(h);
        let mut all = true;
        for probe in 0..PROBES {
            let bit = (g >> (9 * probe)) & 511;
            all &= block.0[(bit >> 6) as usize] & (1 << (bit & 63)) != 0;
        }
        all
    }

    #[inline(always)]
    pub(crate) fn block_pointer(&self, canonical: u64) -> *const u8 {
        std::ptr::from_ref(self.block(hashes(canonical).0)).cast()
    }

    pub(crate) fn k(&self) -> usize {
        self.k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_bases(state: &mut u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                *state ^= *state << 13;
                *state ^= *state >> 7;
                *state ^= *state << 17;
                b"ACGT"[(*state >> 32) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn every_graph_kmer_passes_and_most_absent_kmers_are_rejected() {
        let k = 31;
        let mut state = 0x1234_5678_9abc_def0u64;
        // Unitigs shorter than, equal to and spanning many packed words.
        let unitigs: Vec<Vec<u8>> = [30, 31, 32, 33, 63, 64, 65, 200, 1000]
            .iter()
            .flat_map(|&len| (0..20).map(move |_| len))
            .map(|len| random_bases(&mut state, len))
            .collect();
        let km: Vec<Vec<u8>> = (0..50).map(|_| random_bases(&mut state, k)).collect();
        let abundant_seq = random_bases(&mut state, k);
        let (forward, reverse) = crate::dna::encode_kmer_pair(&abundant_seq).unwrap();
        let encoded = crate::pseudoalign::bifrost::encoded_unitigs_for_test(&unitigs, &km, 23);
        for threads in [1, 3, 8] {
            let filter = KmerFilter::build_from(&encoded, &[forward.min(reverse)], k, threads)
                .expect("filter");
            for sequence in unitigs
                .iter()
                .chain(&km)
                .chain(std::iter::once(&abundant_seq))
            {
                for window in sequence.windows(k) {
                    let (forward, reverse) = crate::dna::encode_kmer_pair(window).unwrap();
                    assert!(filter.may_contain(forward.min(reverse)));
                    // Both orientations map to the same canonical word.
                    let (rc_forward, rc_reverse) = (reverse, forward);
                    assert!(filter.may_contain(rc_forward.min(rc_reverse)));
                }
            }
            let absent = (0..20_000)
                .filter(|_| {
                    let probe = random_bases(&mut state, k);
                    let (forward, reverse) = crate::dna::encode_kmer_pair(&probe).unwrap();
                    filter.may_contain(forward.min(reverse))
                })
                .count();
            assert!(
                absent < 20_000 / 20,
                "false positive rate too high: {absent}"
            );
        }
    }
}
