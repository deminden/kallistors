//! Exact large-graph edges with adaptive prefix-relative rank storage.
use super::{
    CachedEdges, EDGE_HAS_NON_SELF_SHIFT, EDGE_OUT_CODE_SHIFT, EDGE_OUT_DEGREE_SHIFT,
    EDGE_STATE_BITS, KmerLookup, MappingKeys, cached_edges_from_found, edge_neighbor_words,
    kmer_mask, reverse_complement_word,
};
use super::{KmerDirectory, prefetch_kmer, prefetch_rank_delta};
use rayon::prelude::*;

pub(super) struct RelativeEdges<const WIDTH: usize> {
    pub(super) words: Vec<[u8; WIDTH]>,
    branches: Box<[[u32; 4]]>,
    starts: Box<[u32]>,
    prefix_shift: usize,
    k: usize,
    padding: usize,
}
impl<const WIDTH: usize> RelativeEdges<WIDTH> {
    const RANK_BITS: usize = {
        assert!(WIDTH >= 4 && WIDTH <= 8);
        if WIDTH <= 6 {
            4 * WIDTH - 2
        } else {
            4 * WIDTH - 4
        }
    };
    const EDGE_BITS: usize = 4 * WIDTH;
    const RANK_MASK: u32 = (1 << Self::RANK_BITS) - 1;
    const LIMIT: u32 = (1 << Self::RANK_BITS) - 5;
    pub(super) fn prefixes(entries: &[u64], lookup: &KmerLookup, k: usize) -> Option<Box<[u32]>> {
        // Choose the smallest exact directory that fits the rank field. The
        // ordinary eight-byte fallback keeps its established three-base layout.
        let choices: &[usize] = if WIDTH <= 6 {
            &[7, 8, 9]
        } else if WIDTH == 7 {
            &[5, 6, 7]
        } else {
            &[3]
        };
        for &bases in choices {
            let bits = 2 * k.min(bases);
            let starts = (0..=(1usize << bits))
                .map(|prefix| {
                    if let KmerDirectory::Ranks(offsets) = &lookup.directory
                        && lookup.shift <= 64 - bits
                    {
                        offsets.get(prefix << (64 - bits - lookup.shift)) as u32
                    } else {
                        entries.partition_point(|word| (word >> (64 - bits)) < prefix as u64) as u32
                    }
                })
                .collect::<Vec<_>>()
                .into_boxed_slice();
            if Self::fits(&starts) {
                return Some(starts);
            }
        }
        None
    }
    fn fits(starts: &[u32]) -> bool {
        starts
            .windows(2)
            .all(|pair| pair[1] - pair[0] <= Self::LIMIT)
    }
    pub(super) fn build(
        entries: &(impl MappingKeys + ?Sized),
        lookup: &KmerLookup,
        k: usize,
        starts: Box<[u32]>,
        parallel: bool,
    ) -> Self {
        Self::build_plan(
            entries,
            lookup,
            k,
            super::adaptive_prefixes::PrefixPlan {
                roots: starts,
                branches: Box::default(),
            },
            parallel,
        )
    }
    pub(super) fn adaptive_prefixes(
        entries: &[u64],
        lookup: &KmerLookup,
        k: usize,
    ) -> super::adaptive_prefixes::PrefixPlan {
        super::adaptive_prefixes::PrefixPlan::new(
            entries,
            lookup,
            k,
            if WIDTH == 4 { 10 } else { 8 },
            Self::LIMIT,
        )
    }
    pub(super) fn build_plan(
        entries: &(impl MappingKeys + ?Sized),
        lookup: &KmerLookup,
        k: usize,
        plan: super::adaptive_prefixes::PrefixPlan,
        parallel: bool,
    ) -> Self {
        Self::build_plan_filtered(entries, lookup, k, plan, parallel, None)
    }
    pub(super) fn build_filtered(
        entries: &(impl MappingKeys + ?Sized),
        lookup: &KmerLookup,
        k: usize,
        plan: super::adaptive_prefixes::PrefixPlan,
        filter: &super::negative_filter::NegativeFilter,
        parallel: bool,
    ) -> Self {
        Self::build_plan_filtered(entries, lookup, k, plan, parallel, Some(filter))
    }
    fn build_plan_filtered(
        entries: &(impl MappingKeys + ?Sized),
        lookup: &KmerLookup,
        k: usize,
        plan: super::adaptive_prefixes::PrefixPlan,
        parallel: bool,
        negative: Option<&super::negative_filter::NegativeFilter>,
    ) -> Self {
        let super::adaptive_prefixes::PrefixPlan {
            roots: starts,
            branches,
        } = plan;
        let padding = 64 - 2 * k;
        let prefix_shift = 64 - (starts.len() - 1).ilog2() as usize;
        let mask = kmer_mask(k);
        let encode = |cached: CachedEdges, word: u64| {
            let reverse = reverse_complement_word(word, k, mask);
            let mut value = 0u64;
            for orientation in 0..2 {
                let bits = (cached.packed >> (orientation * EDGE_STATE_BITS))
                    & ((1 << EDGE_STATE_BITS) - 1);
                let degree = (bits >> EDGE_OUT_DEGREE_SHIFT) & 7;
                let code = (bits >> EDGE_OUT_CODE_SHIFT) & 3;
                let edge = if degree == 1 {
                    let forward = cached.packed & (1 << (2 * EDGE_STATE_BITS + orientation)) != 0;
                    let oriented = if orientation == 0 { word } else { reverse };
                    let prefix = Self::target_prefix(
                        oriented,
                        u64::from(code),
                        forward,
                        k,
                        padding,
                        prefix_shift,
                    );
                    let relative = cached.next[orientation] - starts[prefix];
                    debug_assert!(relative < Self::LIMIT);
                    relative
                        | if WIDTH <= 6 {
                            0
                        } else {
                            u32::from(forward) << (Self::RANK_BITS + 2)
                        }
                } else {
                    if WIDTH <= 6 {
                        debug_assert_eq!(
                            (bits >> EDGE_HAS_NON_SELF_SHIFT) & 1,
                            u32::from(degree >= 2)
                        );
                        Self::LIMIT + degree
                    } else {
                        (Self::LIMIT + degree)
                            | (((bits >> EDGE_HAS_NON_SELF_SHIFT) & 1) << (Self::RANK_BITS + 2))
                    }
                };
                value |=
                    u64::from(edge | (code << Self::RANK_BITS)) << (orientation * Self::EDGE_BITS);
            }
            value
        };
        let mut words = Vec::<[u8; WIDTH]>::with_capacity(entries.len());
        let fill = |(chunk, slots): (usize, &mut [std::mem::MaybeUninit<[u8; WIDTH]>])| {
            let mut source_words = [0u64; 1024];
            let source = entries.window(lookup, chunk * 1024, slots.len(), &mut source_words);
            let mut queries = [0u64; 256];
            for (source, slots) in source.chunks(32).zip(slots.chunks_mut(32)) {
                for (i, &word) in source.iter().enumerate() {
                    queries[i * 8..i * 8 + 8].copy_from_slice(&edge_neighbor_words(word, k, mask));
                }
                let found = if let Some(negative) = negative {
                    let admitted = negative.admit(&queries, source.len() * 8);
                    let mut points = [0u64; 256];
                    let mut destinations = [0u16; 256];
                    let mut count = 0;
                    for i in 0..source.len() * 8 {
                        if admitted[i] {
                            points[count] = queries[i];
                            destinations[count] = i as u16;
                            count += 1;
                        }
                    }
                    let hits = entries.find_many_ranks(lookup, points, count);
                    let mut found = [u32::MAX; 256];
                    for i in 0..count {
                        found[destinations[i] as usize] = hits[i];
                    }
                    found
                } else if WIDTH < 8
                    && let Some(raw) = entries.raw_words()
                {
                    Self::find_neighbors(raw, lookup, source, k, &queries)
                        .map(|hit| hit.map_or(u32::MAX, |index| index as u32))
                } else {
                    entries.find_many_ranks(lookup, queries, source.len() * 8)
                };
                for (i, (&word, slot)) in source.iter().zip(slots).enumerate() {
                    let encoded = if WIDTH <= 6 {
                        let neighbors = &queries[i * 8..i * 8 + 8];
                        let found = &found[i * 8..i * 8 + 8];
                        let mut value = 0u64;
                        for orientation in 0..2 {
                            let mut presence = 0u32;
                            for code in 0..4 {
                                presence |=
                                    u32::from(found[2 * code + orientation] != u32::MAX) << code;
                            }
                            let degree = presence.count_ones();
                            let code = if presence == 0 {
                                0
                            } else if orientation == 0 {
                                31 - presence.leading_zeros()
                            } else {
                                3 - presence.trailing_zeros()
                            };
                            let edge = if degree == 1 {
                                let slot = if orientation == 0 {
                                    2 * code as usize
                                } else {
                                    2 * (3 - code) as usize + 1
                                };
                                let start = Self::target_start(
                                    &starts,
                                    &branches,
                                    prefix_shift,
                                    neighbors[slot],
                                );
                                let relative = found[slot] - start;
                                debug_assert!(relative < Self::LIMIT);
                                relative
                            } else {
                                Self::LIMIT + degree
                            };
                            value |= u64::from(edge | (code << Self::RANK_BITS))
                                << (orientation * Self::EDGE_BITS);
                        }
                        value
                    } else {
                        encode(
                            cached_edges_from_found(
                                word,
                                k,
                                mask,
                                queries[i * 8..i * 8 + 8]
                                    .try_into()
                                    .expect("eight neighbors"),
                                std::array::from_fn(|slot| {
                                    let rank = found[i * 8 + slot];
                                    (rank != u32::MAX).then_some(rank as usize)
                                }),
                            ),
                            word,
                        )
                    };
                    slot.write(
                        encoded.to_le_bytes()[..WIDTH]
                            .try_into()
                            .expect("edge width"),
                    );
                }
            }
        };
        let slots = &mut words.spare_capacity_mut()[..entries.len()];
        if parallel {
            slots.par_chunks_mut(1024).enumerate().for_each(fill);
        } else {
            slots.chunks_mut(1024).enumerate().for_each(fill);
        }
        // SAFETY: disjoint slots cover all keys and each is initialized exactly
        // once. Both paths join before publication; panic leaves length zero.
        unsafe {
            words.set_len(entries.len());
        }
        Self {
            words,
            branches,
            starts,
            prefix_shift,
            k,
            padding,
        }
    }

    fn find_neighbors(
        entries: &[u64],
        lookup: &KmerLookup,
        source: &[u64],
        k: usize,
        neighbors: &[u64; 256],
    ) -> [Option<usize>; 256] {
        let padding = 64 - 2 * k;
        let KmerDirectory::Ranks(offsets) = &lookup.directory else {
            return lookup.find_many_count(entries, *neighbors, source.len() * 8);
        };
        if lookup.shift < padding + 2 {
            return lookup.find_many_count(entries, *neighbors, source.len() * 8);
        }
        debug_assert!(source.len() <= 32);
        let mut families = [(0u64, 0u8, 0usize, 0usize); 64];
        let mut family_count = 0;
        let mut points = [0u64; 256];
        let mut destinations = [0u16; 256];
        let mut point_count = 0;
        let mask = kmer_mask(k);
        for (i, &word) in source.iter().enumerate() {
            let reverse = reverse_complement_word(word, k, mask);
            for (orientation, oriented) in [word, reverse].into_iter().enumerate() {
                let base = super::shift_append_word(oriented, 0, k, mask);
                let mut family_codes = 0u8;
                for code in 0..4usize {
                    let slot = if orientation == 0 {
                        2 * code
                    } else {
                        2 * (3 - code) + 1
                    };
                    let append = base | ((code as u64) << padding);
                    family_codes |= u8::from(neighbors[i * 8 + slot] == append) << code;
                }
                // Single requests retain the ordinary point lookup. Families
                // are useful only when at least two complete keys share work.
                if family_codes.count_ones() < 2 {
                    family_codes = 0;
                }
                if family_codes != 0 {
                    families[family_count] = (base, family_codes, i, orientation);
                    family_count += 1;
                }
                for code in 0..4usize {
                    if family_codes & (1 << code) != 0 {
                        continue;
                    }
                    let slot = i * 8
                        + if orientation == 0 {
                            2 * code
                        } else {
                            2 * (3 - code) + 1
                        };
                    points[point_count] = neighbors[slot];
                    destinations[point_count] = slot as u16;
                    point_count += 1;
                }
            }
        }
        let mut bounds = [(0usize, 0usize); 64];
        if !offsets.deltas.is_empty() {
            for &(base, _, _, _) in &families[..family_count] {
                let bucket = (base >> lookup.shift) as usize;
                prefetch_rank_delta(&offsets.deltas[(bucket >> 1) * 3]);
                prefetch_rank_delta(&offsets.bases[bucket >> 8]);
            }
        }
        for (i, &(base, _, _, _)) in families[..family_count].iter().enumerate() {
            bounds[i] = offsets.range((base >> lookup.shift) as usize);
            let (start, end) = bounds[i];
            if start < end {
                prefetch_kmer(
                    &entries[if end - start <= 8 {
                        start
                    } else {
                        start + (end - start) / 2
                    }],
                );
            }
        }
        // Ordinary requests give prefetched family bounds and keys time to
        // arrive. Every hit still requires complete numeric key equality.
        let point_found = lookup.find_many_count(entries, points, point_count);
        let mut found = [None; 256];
        for i in 0..point_count {
            found[destinations[i] as usize] = point_found[i];
        }
        let code_mask = 3u64 << padding;
        for (i, &(base, codes, source, orientation)) in families[..family_count].iter().enumerate()
        {
            let (start, end) = bounds[i];
            let first = if end - start <= 8 {
                start
            } else {
                start + entries[start..end].partition_point(|&key| key < base)
            };
            let upper = base | code_mask;
            for (relative, &key) in entries[first..end].iter().enumerate() {
                if key > upper {
                    break;
                }
                if key & !code_mask != base {
                    continue;
                }
                let code = ((key >> padding) & 3) as usize;
                if codes & (1 << code) != 0 {
                    let slot = source * 8
                        + if orientation == 0 {
                            2 * code
                        } else {
                            2 * (3 - code) + 1
                        };
                    found[slot] = Some(first + relative);
                }
            }
        }
        found
    }

    #[inline(always)]
    fn target_prefix(
        oriented: u64,
        code: u64,
        forward: bool,
        k: usize,
        padding: usize,
        prefix_shift: usize,
    ) -> usize {
        if WIDTH == 8 {
            if forward {
                (((oriented << 2) | (code << padding)) >> 58) as usize
            } else {
                let tail = !(oriented >> padding) & 15;
                (((3 - code) << 4) | ((tail & 3) << 2) | (tail >> 2)) as usize
            }
        } else {
            let next = super::shift_append_word(oriented, code, k, kmer_mask(k));
            let canonical = if forward {
                next
            } else {
                reverse_complement_word(next, k, kmer_mask(k))
            };
            (canonical >> prefix_shift) as usize
        }
    }
    #[inline]
    fn target_start(starts: &[u32], branches: &[[u32; 4]], shift: usize, word: u64) -> u32 {
        if WIDTH <= 5 {
            super::adaptive_prefixes::start(starts, branches, shift, word)
        } else {
            starts[(word >> shift) as usize]
        }
    }
    pub(super) fn scratch_bytes(&self) -> usize {
        self.words.len() * WIDTH + self.starts.len() * 4 + self.branches.len() * 16
    }
    #[inline(always)]
    fn edge(&self, index: usize, orientation: usize) -> u32 {
        let pair = &self.words[index];
        if WIDTH == 8 {
            (u64::from_le_bytes(pair.as_slice().try_into().expect("eight-byte edge"))
                >> (orientation * 32)) as u32
        } else if WIDTH == 4 {
            let first = orientation * 2;
            u32::from(u16::from_le_bytes([pair[first], pair[first + 1]]))
        } else if WIDTH <= 6 {
            let first = orientation * (WIDTH / 2);
            let shift = if WIDTH == 5 { orientation * 4 } else { 0 };
            (u32::from_le_bytes([pair[first], pair[first + 1], pair[first + 2], 0]) >> shift)
                & ((1 << Self::EDGE_BITS) - 1)
        } else {
            let first = orientation * 3;
            (u32::from_le_bytes(pair[first..first + 4].try_into().unwrap()) >> (orientation * 4))
                & 0x0fff_ffff
        }
    }
    #[inline]
    fn degree(edge: u32) -> u32 {
        let rank = edge & Self::RANK_MASK;
        if rank < Self::LIMIT {
            1
        } else {
            rank - Self::LIMIT
        }
    }
    #[inline]
    pub(super) fn is_internal(&self, index: usize) -> bool {
        (self.edge(index, 0) & Self::RANK_MASK) < Self::LIMIT
            && (self.edge(index, 1) & Self::RANK_MASK) < Self::LIMIT
    }
    #[inline]
    pub(super) fn incoming_degree(&self, index: usize, forward: bool) -> u8 {
        Self::degree(self.edge(index, usize::from(forward))) as u8
    }
    #[inline(always)]
    pub(super) fn outgoing(
        &self,
        index: usize,
        forward: bool,
        oriented: u64,
    ) -> Option<(usize, bool, u64)> {
        let edge = self.edge(index, usize::from(!forward));
        let relative = edge & Self::RANK_MASK;
        if relative >= Self::LIMIT {
            return None;
        }
        let code = u64::from((edge >> Self::RANK_BITS) & 3);
        let (start, forward) = if WIDTH <= 6 {
            let next = super::shift_append_word(oriented, code, self.k, kmer_mask(self.k));
            let reverse = reverse_complement_word(next, self.k, kmer_mask(self.k));
            let canonical = next.min(reverse);
            let start =
                Self::target_start(&self.starts, &self.branches, self.prefix_shift, canonical);
            if WIDTH >= 6 {
                let prefix = (canonical >> self.prefix_shift) as usize;
                debug_assert!(relative < self.starts[prefix + 1] - start);
            }
            (start, next <= reverse)
        } else {
            let forward = edge & (1 << (Self::RANK_BITS + 2)) != 0;
            let prefix = Self::target_prefix(
                oriented,
                code,
                forward,
                self.k,
                self.padding,
                self.prefix_shift,
            );
            debug_assert!(relative < self.starts[prefix + 1] - self.starts[prefix]);
            (self.starts[prefix], forward)
        };
        Some(((start + relative) as usize, forward, code))
    }

    #[inline]
    pub(super) fn prefetch_incoming(&self, index: usize, forward: bool) {
        let pair = &self.words[index];
        let offset = usize::from(forward) * (WIDTH / 2);
        prefetch_rank_delta(&pair[offset]);
    }
    #[inline]
    pub(super) fn continuation(
        &self,
        index: usize,
        forward: bool,
        oriented: u64,
    ) -> Option<(usize, bool, u64)> {
        self.outgoing(index, forward, oriented)
            .filter(|&(target, forward, _)| self.incoming_degree(target, forward) == 1)
    }
    pub(super) fn state_bits(&self, index: usize, forward: bool, canonical: u64) -> u32 {
        let edge = self.edge(index, usize::from(!forward));
        let degree = Self::degree(edge);
        let non_self = if degree == 1 {
            let oriented = if forward {
                canonical
            } else {
                reverse_complement_word(canonical, self.k, kmer_mask(self.k))
            };
            self.outgoing(index, forward, oriented)
                .expect("single edge")
                .0
                != index
        } else if WIDTH <= 6 {
            // k >= 3: at most one appended base returns to the same canonical
            // k-mer. A homopolymer self edge and a reverse self edge cannot
            // coexist, so every degree >= 2 has a non-self neighbor.
            degree >= 2
        } else {
            edge & (1 << (Self::RANK_BITS + 2)) != 0
        };
        u32::from(self.incoming_degree(index, forward))
            | (degree << EDGE_OUT_DEGREE_SHIFT)
            | (((edge >> Self::RANK_BITS) & 3) << EDGE_OUT_CODE_SHIFT)
            | (u32::from(non_self) << EDGE_HAS_NON_SELF_SHIFT)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn relative_ranks_preserve_prefixes_sentinels_and_fallback_boundary() {
        fn check<const WIDTH: usize>() {
            type Table<const W: usize> = RelativeEdges<W>;
            let limit = Table::<WIDTH>::LIMIT;
            let bits = Table::<WIDTH>::RANK_BITS;
            let prefix_bases = if WIDTH == 4 {
                9
            } else if WIDTH == 5 {
                8
            } else if WIDTH == 6 {
                9
            } else if WIDTH == 7 {
                7
            } else {
                3
            };
            for size in [limit - 1, limit, limit + 1] {
                let mut starts = [size; 65];
                starts[0] = 0;
                assert_eq!(Table::<WIDTH>::fits(&starts), size <= limit);
            }
            for degree in [0, 2, 3, 4] {
                assert_eq!(Table::<WIDTH>::degree(limit + degree), degree);
            }
            for rank in [0, limit - 1] {
                assert_eq!(
                    Table::<WIDTH>::degree(
                        rank | (3 << bits) | if WIDTH <= 6 { 0 } else { 1 << (bits + 2) }
                    ),
                    1
                );
            }
            for k in (3..=31).step_by(2) {
                let padding = 64 - 2 * k;
                let prefix_shift = 64 - 2 * k.min(prefix_bases);
                for first in 0..64u64 {
                    for last in 0..4u64 {
                        for code in 0..4u64 {
                            let oriented = ((first << 58) | (last << padding)) & kmer_mask(k);
                            let next =
                                super::super::shift_append_word(oriented, code, k, kmer_mask(k));
                            for forward in [true, false] {
                                let canonical = if forward {
                                    next
                                } else {
                                    reverse_complement_word(next, k, kmer_mask(k))
                                };
                                assert_eq!(
                                    Table::<WIDTH>::target_prefix(
                                        oriented,
                                        code,
                                        forward,
                                        k,
                                        padding,
                                        prefix_shift
                                    ),
                                    (canonical >> prefix_shift) as usize
                                );
                            }
                        }
                    }
                }
            }
            let mut starts = vec![limit; 1 << (2 * prefix_bases)];
            starts.insert(0, 0);
            let encoded = u64::from((limit - 1) | if WIDTH <= 6 { 0 } else { 1 << (bits + 2) });
            let words = vec![encoded.to_le_bytes()[..WIDTH].try_into().unwrap()];
            let edges = Table::<WIDTH> {
                words,
                branches: Box::default(),
                starts: starts.into_boxed_slice(),
                k: 31,
                padding: 2,
                prefix_shift: 64 - 2 * prefix_bases,
            };
            assert_eq!(
                edges.outgoing(0, true, 0),
                Some(((limit - 1) as usize, true, 0))
            );
        }
        check::<8>();
        check::<7>();
        check::<6>();
        check::<{ super::super::ADAPTIVE_EDGE_WIDTH }>();
    }
}
