//! Compact edge storage, topology discovery and finalization.

use super::*;

#[derive(Clone, Copy)]
pub(super) struct OrientedNeighbor {
    pub(super) oriented: u64,
    pub(super) index: usize,
    pub(super) forward: bool,
}

pub(super) struct GraphTopology<'a, K: MappingKeys + ?Sized = [u64]> {
    pub(super) entries: &'a K,
    pub(super) lookup: &'a KmerLookup,
    pub(super) edge_states: &'a EdgeStateCache,
    pub(super) k: usize,
    pub(super) mask: u64,
    pub(super) last_shift: usize,
}

pub(super) struct CachedEdges {
    pub(super) packed: u32,
    pub(super) next: [u32; 2],
}

// Large graphs need 31 target bits. Non-singleton edges use unused ranks
// above the graph size for their exact degree; bit31 then stores non-self.
// Keeping both base codes separately preserves the cached graph walk.
#[repr(C, packed)]
pub(super) struct RankedEdges {
    pub(super) links: u64,
    pub(super) codes: AtomicU8,
}

impl RankedEdges {
    pub(super) fn new(cached: CachedEdges, count: usize) -> Self {
        debug_assert!(count <= (1 << 31) - 5);
        let mut links = 0u64;
        let mut codes = 0u8;
        for orientation in 0..2 {
            let bits =
                (cached.packed >> (orientation * EDGE_STATE_BITS)) & ((1 << EDGE_STATE_BITS) - 1);
            let degree = (bits >> EDGE_OUT_DEGREE_SHIFT) & 7;
            let edge = if degree == 1 {
                debug_assert!(cached.next[orientation] < count as u32);
                cached.next[orientation]
                    | (((cached.packed >> (2 * EDGE_STATE_BITS + orientation)) & 1) << 31)
            } else {
                (count as u32 + degree) | (((bits >> EDGE_HAS_NON_SELF_SHIFT) & 1) << 31)
            };
            links |= u64::from(edge) << (orientation * 32);
            codes |= (((bits >> EDGE_OUT_CODE_SHIFT) & 3) as u8) << (orientation * 2);
        }
        Self {
            links,
            codes: AtomicU8::new(codes),
        }
    }

    pub(super) fn build(
        entries: &[u64],
        lookup: &KmerLookup,
        k: usize,
        parallel: bool,
    ) -> Vec<Self> {
        let mut records = Vec::<Self>::with_capacity(entries.len());
        let fill = |(words, slots): (&[u64], &mut [std::mem::MaybeUninit<Self>])| {
            let mask = kmer_mask(k);
            let mut queries = [0u64; 256];
            for (words, slots) in words.chunks(32).zip(slots.chunks_mut(32)) {
                for (i, &word) in words.iter().enumerate() {
                    queries[i * 8..i * 8 + 8].copy_from_slice(&edge_neighbor_words(word, k, mask));
                }
                let found = lookup.find_many_count(entries, queries, words.len() * 8);
                for (i, (&word, slot)) in words.iter().zip(slots).enumerate() {
                    slot.write(Self::new(
                        cached_edges_from_found(
                            word,
                            k,
                            mask,
                            queries[i * 8..i * 8 + 8]
                                .try_into()
                                .expect("eight neighbors"),
                            found[i * 8..i * 8 + 8].try_into().expect("eight results"),
                        ),
                        entries.len(),
                    ));
                }
            }
        };
        let destinations = &mut records.spare_capacity_mut()[..entries.len()];
        if parallel {
            entries
                .par_chunks(1024)
                .zip(destinations.par_chunks_mut(1024))
                .for_each(fill);
        } else {
            entries
                .chunks(1024)
                .zip(destinations.chunks_mut(1024))
                .for_each(fill);
        }
        // SAFETY: each input key has one disjoint output slot. Both iteration
        // paths join after initializing every slot; on panic the length is zero
        // and the records contain no resource-owning fields.
        unsafe {
            records.set_len(entries.len());
        }
        Self::mark_continuations(&records, parallel);
        records
    }

    pub(super) fn mark_continuations(records: &[Self], parallel: bool) {
        // Only spare flag bits change. Atomic byte access permits neighboring
        // records to be inspected concurrently without copying the edge table.
        let batch = |batch: &[Self]| {
            for record in batch {
                for orientation in 0..2 {
                    let target = (record.edge(orientation) & 0x7fff_ffff) as usize;
                    if target < records.len() {
                        prefetch_rank_delta(&records[target].codes);
                    }
                }
            }
            for record in batch {
                let mut flags = 0;
                for orientation in 0..2 {
                    let edge = record.edge(orientation);
                    let target = (edge & 0x7fff_ffff) as usize;
                    if target < records.len() {
                        let incoming = records[target].edge(usize::from(edge & (1 << 31) != 0));
                        flags |= u8::from((incoming & 0x7fff_ffff) < records.len() as u32)
                            << (4 + orientation);
                    }
                }
                record.codes.store(
                    record.codes.load(Ordering::Relaxed) | flags,
                    Ordering::Relaxed,
                );
            }
        };
        if parallel {
            records.par_chunks(32).for_each(batch);
        } else {
            records.chunks(32).for_each(batch);
        }
    }

    #[inline]
    pub(super) fn edge(&self, orientation: usize) -> u32 {
        (self.links >> (orientation * 32)) as u32
    }

    #[inline]
    pub(super) fn degree(edge: u32, count: usize) -> u32 {
        let rank = edge & 0x7fff_ffff;
        if rank < count as u32 {
            1
        } else {
            rank - count as u32
        }
    }
}

pub(super) enum EdgeStorage {
    Compact(Vec<u64>),
    Follow(Vec<u64>),
    Ranked(Vec<RankedEdges>),
    Relative(RelativeEdges),
    Relative7(PackedRelativeEdges),
    Relative6(DenseRelativeEdges),
    Adaptive(AdaptiveRelativeEdges),
    Wide(Vec<CachedEdges>),
}

pub(super) struct EdgeStateCache {
    pub(super) packed: EdgeStorage,
}

// Align like the final u64 table so final collection can reuse its allocation.
#[repr(C, align(8))]
pub(super) struct AtomicEdgePair(pub(super) [AtomicU32; 2]);

// Reuse the atomic allocation when its size and alignment match immutable
// u64 words; other allocator layouts retain the ordinary collection fallback. All discovery workers have joined
// before exclusive finalization; no atomics or aliases survive this conversion.
pub(super) fn finish_atomic_edge_records(
    mut records: Vec<AtomicEdgePair>,
    parallel: bool,
) -> EdgeStorage {
    const {
        assert!(std::mem::size_of::<AtomicEdgePair>() == std::mem::size_of::<u64>());
    }
    let count = records.len();
    let follow = count < (1 << 28) - 4;
    let finish = |pair: &mut AtomicEdgePair| {
        for edge in &mut pair.0 {
            let value = edge.get_mut();
            if *value & (1 << 31) == 0 {
                *value = if follow {
                    (count as u32 + (*value & 7))
                        | (*value & (3 << 28))
                        | (u32::from(*value & 8 != 0) << 30)
                } else {
                    *value & !0xf0
                };
            }
        }
    };
    if parallel {
        records.par_iter_mut().for_each(finish);
    } else {
        records.iter_mut().for_each(finish);
    }
    if follow {
        // Rank fields are now immutable: only bit31 changes in this pass.
        // Atomic reads permit concurrent neighbor inspection without a second
        // snapshot. All target ranks and incoming degrees remain exact.
        let mark = |pair: &AtomicEdgePair| {
            for edge in &pair.0 {
                let value = edge.load(Ordering::Relaxed);
                let target = (value & 0x0fff_ffff) as usize;
                if target < count {
                    let incoming = records[target].0[usize::from(value & (1 << 30) != 0)]
                        .load(Ordering::Relaxed);
                    let flag = u32::from((incoming & 0x0fff_ffff) < count as u32) << 31;
                    edge.store((value & !(1 << 31)) | flag, Ordering::Relaxed);
                }
            }
        };
        let batch = |pairs: &[AtomicEdgePair]| {
            // Neighbor ranks are stable. Hint independent target lines before
            // resolving their flags, bounding scratch and preserving ownership.
            for pair in pairs {
                for edge in &pair.0 {
                    let value = edge.load(Ordering::Relaxed);
                    let target = (value & 0x0fff_ffff) as usize;
                    if target < count {
                        prefetch_rank_delta(
                            &records[target].0[usize::from(value & (1 << 30) != 0)],
                        );
                    }
                }
            }
            pairs.iter().for_each(mark);
        };
        if parallel {
            records.par_chunks(16).for_each(batch);
        } else {
            records.chunks(16).for_each(batch);
        }
    }
    let words = if std::mem::align_of::<AtomicEdgePair>() != std::mem::align_of::<u64>() {
        records
            .into_iter()
            .map(|pair| compact_pair_word(pair.0.map(AtomicU32::into_inner)))
            .collect()
    } else {
        if cfg!(target_endian = "big") {
            let swap = |pair: &mut AtomicEdgePair| pair.0.swap(0, 1);
            if parallel {
                records.par_iter_mut().for_each(swap);
            } else {
                records.iter_mut().for_each(swap);
            }
        }
        let mut records = std::mem::ManuallyDrop::new(records);
        // SAFETY: discovery, normalization and atomic flag workers have joined;
        // no aliases survive. Size is asserted, alignment checked, and native
        // byte order normalized. Both element types permit all bit patterns and
        // have no drop glue. Transfer the allocation exactly once.
        unsafe {
            Vec::from_raw_parts(
                records.as_mut_ptr().cast::<u64>(),
                records.len(),
                records.capacity(),
            )
        }
    };
    if follow {
        EdgeStorage::Follow(words)
    } else {
        EdgeStorage::Compact(words)
    }
}

#[inline]
pub(super) fn follow_degree(record: u32, count: usize) -> u32 {
    let rank = record & 0x0fff_ffff;
    if rank < count as u32 {
        1
    } else {
        rank - count as u32
    }
}

impl EdgeStateCache {
    pub(super) fn scratch_bytes(&self) -> usize {
        match &self.packed {
            EdgeStorage::Relative6(records) => records.scratch_bytes(),
            EdgeStorage::Adaptive(records) => records.scratch_bytes(),
            EdgeStorage::Compact(records) | EdgeStorage::Follow(records) => records.len() * 8,
            EdgeStorage::Ranked(records) => records.len() * std::mem::size_of::<RankedEdges>(),
            EdgeStorage::Relative(records) => records.scratch_bytes(),
            EdgeStorage::Relative7(records) => records.scratch_bytes(),
            EdgeStorage::Wide(records) => records.len() * std::mem::size_of::<CachedEdges>(),
        }
    }

    pub(super) fn retire(self) {
        match self.packed {
            EdgeStorage::Relative6(records) => super::super::scratch_pages::retire(records.words),
            EdgeStorage::Adaptive(records) => super::super::scratch_pages::retire(records.words),
            EdgeStorage::Relative7(records) => super::super::scratch_pages::retire(records.words),
            EdgeStorage::Relative(records) => super::super::scratch_pages::retire(records.words),
            EdgeStorage::Compact(records) | EdgeStorage::Follow(records) => {
                super::super::scratch_pages::retire(records)
            }
            EdgeStorage::Ranked(records) => super::super::scratch_pages::retire(records),
            EdgeStorage::Wide(records) => super::super::scratch_pages::retire(records),
        }
    }

    pub(super) fn build(
        entries: &[u64],
        lookup: &KmerLookup,
        k: usize,
        threads: usize,
    ) -> Result<Self> {
        let mask = kmer_mask(k);
        let collect = |parallel| {
            if entries.len() <= (1 << 28) {
                // A serial zero-value collection can use the allocator's
                // zeroed pages without an extra parallel first-touch pass.
                let result = (0..entries.len())
                    .map(|_| AtomicEdgePair([AtomicU32::new(0), AtomicU32::new(0)]))
                    .collect::<Vec<_>>();
                let visit = |(chunk, words): (usize, &[u64])| {
                    let mut scratch = EdgeQueryScratch::default();
                    for (batch, words) in words.chunks(32).enumerate() {
                        cache_symmetric_edges(
                            chunk * 4096 + batch * 32,
                            words,
                            entries,
                            &result,
                            lookup,
                            k,
                            mask,
                            &mut scratch,
                        );
                    }
                };
                if parallel {
                    entries.par_chunks(4096).enumerate().for_each(visit);
                } else {
                    entries.chunks(4096).enumerate().for_each(visit);
                }
                finish_atomic_edge_records(result, parallel)
            } else if entries.len() <= (1 << 31) - 5 {
                if let Some(starts) = DenseRelativeEdges::prefixes(entries, lookup, k) {
                    EdgeStorage::Relative6(DenseRelativeEdges::build(
                        entries, lookup, k, starts, parallel,
                    ))
                } else if let Some(starts) = PackedRelativeEdges::prefixes(entries, lookup, k) {
                    EdgeStorage::Relative7(PackedRelativeEdges::build(
                        entries, lookup, k, starts, parallel,
                    ))
                } else if let Some(starts) = RelativeEdges::prefixes(entries, lookup, k) {
                    EdgeStorage::Relative(RelativeEdges::build(
                        entries, lookup, k, starts, parallel,
                    ))
                } else {
                    EdgeStorage::Ranked(RankedEdges::build(entries, lookup, k, parallel))
                }
            } else {
                let encode = |&word| cache_edge_states(word, entries, lookup, k, mask);
                EdgeStorage::Wide(if parallel {
                    entries.par_iter().map(encode).collect()
                } else {
                    entries.iter().map(encode).collect()
                })
            }
        };
        let packed = if threads > 1 && entries.len() >= PARALLEL_KMER_COLLECTION_MIN {
            let pool = IndexWorkers::new(threads)?;
            pool.install(|| collect(true))
        } else {
            collect(false)
        };
        Ok(Self { packed })
    }

    #[inline]
    pub(super) fn is_internal(&self, idx: usize) -> bool {
        match &self.packed {
            EdgeStorage::Relative6(records) => records.is_internal(idx),
            EdgeStorage::Adaptive(records) => records.is_internal(idx),
            EdgeStorage::Relative7(records) => records.is_internal(idx),
            EdgeStorage::Relative(records) => records.is_internal(idx),
            EdgeStorage::Ranked(entries) => {
                let links = entries[idx].links;
                ((links as u32) & 0x7fff_ffff) < entries.len() as u32
                    && (((links >> 32) as u32) & 0x7fff_ffff) < entries.len() as u32
            }
            EdgeStorage::Follow(entries) => {
                let value = entries[idx];
                ((value as u32) & 0x0fff_ffff) < entries.len() as u32
                    && (((value >> 32) as u32) & 0x0fff_ffff) < entries.len() as u32
            }
            EdgeStorage::Wide(entries) => entries[idx].packed & 0x3f == 9,
            EdgeStorage::Compact(entries) => {
                (entries[idx] as u32) & ((entries[idx] >> 32) as u32) & (1 << 31) != 0
            }
        }
    }

    pub(super) fn state_bits(&self, idx: usize, forward: bool, canonical: u64) -> u32 {
        let orientation = usize::from(!forward);
        match &self.packed {
            EdgeStorage::Relative6(records) => records.state_bits(idx, forward, canonical),
            EdgeStorage::Adaptive(records) => records.state_bits(idx, forward, canonical),
            EdgeStorage::Relative7(records) => records.state_bits(idx, forward, canonical),
            EdgeStorage::Relative(records) => records.state_bits(idx, forward, canonical),
            EdgeStorage::Ranked(entries) => {
                let record = &entries[idx];
                let edge = record.edge(orientation);
                let rank = edge & 0x7fff_ffff;
                let non_self = if rank < entries.len() as u32 {
                    rank as usize != idx
                } else {
                    edge & (1 << 31) != 0
                };
                RankedEdges::degree(record.edge(1 - orientation), entries.len())
                    | (RankedEdges::degree(edge, entries.len()) << EDGE_OUT_DEGREE_SHIFT)
                    | (u32::from((record.codes.load(Ordering::Relaxed) >> (orientation * 2)) & 3)
                        << EDGE_OUT_CODE_SHIFT)
                    | (u32::from(non_self) << EDGE_HAS_NON_SELF_SHIFT)
            }
            EdgeStorage::Follow(entries) => {
                let edge = (entries[idx] >> (orientation * 32)) as u32;
                let incoming = (entries[idx] >> ((1 - orientation) * 32)) as u32;
                let rank = edge & 0x0fff_ffff;
                let non_self = if rank < entries.len() as u32 {
                    rank as usize != idx
                } else {
                    edge & (1 << 30) != 0
                };
                follow_degree(incoming, entries.len())
                    | (follow_degree(edge, entries.len()) << EDGE_OUT_DEGREE_SHIFT)
                    | (((edge >> 28) & 3) << EDGE_OUT_CODE_SHIFT)
                    | (u32::from(non_self) << EDGE_HAS_NON_SELF_SHIFT)
            }
            EdgeStorage::Wide(entries) => {
                (entries[idx].packed >> (orientation * EDGE_STATE_BITS))
                    & ((1 << EDGE_STATE_BITS) - 1)
            }
            EdgeStorage::Compact(entries) => {
                let edge = (entries[idx] >> (orientation * 32)) as u32;
                let degree = |value: u32| if value & (1 << 31) != 0 { 1 } else { value & 7 };
                let non_self = if edge & (1 << 31) != 0 {
                    (edge & 0x0fff_ffff) as usize != idx
                } else {
                    edge & 8 != 0
                };
                degree((entries[idx] >> ((1 - orientation) * 32)) as u32)
                    | (degree(edge) << EDGE_OUT_DEGREE_SHIFT)
                    | (((edge >> 28) & 3) << EDGE_OUT_CODE_SHIFT)
                    | (u32::from(non_self) << EDGE_HAS_NON_SELF_SHIFT)
            }
        }
    }

    #[inline]
    pub(super) fn incoming_degree(&self, idx: usize, forward: bool) -> u8 {
        let orientation = usize::from(!forward);
        match &self.packed {
            EdgeStorage::Relative6(records) => records.incoming_degree(idx, forward),
            EdgeStorage::Adaptive(records) => records.incoming_degree(idx, forward),
            EdgeStorage::Relative7(records) => records.incoming_degree(idx, forward),
            EdgeStorage::Relative(records) => records.incoming_degree(idx, forward),
            EdgeStorage::Ranked(entries) => {
                RankedEdges::degree(entries[idx].edge(1 - orientation), entries.len()) as u8
            }
            EdgeStorage::Follow(entries) => follow_degree(
                (entries[idx] >> ((1 - orientation) * 32)) as u32,
                entries.len(),
            ) as u8,
            EdgeStorage::Wide(entries) => {
                ((entries[idx].packed >> (orientation * EDGE_STATE_BITS)) & 7) as u8
            }
            EdgeStorage::Compact(entries) => {
                let incoming = (entries[idx] >> ((1 - orientation) * 32)) as u32;
                if incoming & (1 << 31) != 0 {
                    1
                } else {
                    (incoming & 7) as u8
                }
            }
        }
    }

    #[inline]
    pub(super) fn outgoing(
        &self,
        idx: usize,
        forward: bool,
        oriented: u64,
    ) -> Option<(usize, bool, u64)> {
        let orientation = usize::from(!forward);
        match &self.packed {
            EdgeStorage::Relative6(records) => records.outgoing(idx, forward, oriented),
            EdgeStorage::Adaptive(records) => records.outgoing(idx, forward, oriented),
            EdgeStorage::Relative7(records) => records.outgoing(idx, forward, oriented),
            EdgeStorage::Relative(records) => records.outgoing(idx, forward, oriented),
            EdgeStorage::Ranked(entries) => {
                let record = &entries[idx];
                let edge = record.edge(orientation);
                let rank = edge & 0x7fff_ffff;
                (rank < entries.len() as u32).then_some((
                    rank as usize,
                    edge & (1 << 31) != 0,
                    u64::from((record.codes.load(Ordering::Relaxed) >> (orientation * 2)) & 3),
                ))
            }
            EdgeStorage::Follow(entries) => {
                let edge = (entries[idx] >> (orientation * 32)) as u32;
                let rank = edge & 0x0fff_ffff;
                (rank < entries.len() as u32).then_some((
                    rank as usize,
                    edge & (1 << 30) != 0,
                    u64::from((edge >> 28) & 3),
                ))
            }
            EdgeStorage::Wide(entries) => {
                let cached = &entries[idx];
                let bits = cached.packed >> (orientation * EDGE_STATE_BITS);
                if (bits >> EDGE_OUT_DEGREE_SHIFT) & 7 != 1 {
                    return None;
                }
                Some((
                    cached.next[orientation] as usize,
                    cached.packed & (1 << (2 * EDGE_STATE_BITS + orientation)) != 0,
                    u64::from((bits >> EDGE_OUT_CODE_SHIFT) & 3),
                ))
            }
            EdgeStorage::Compact(entries) => {
                let edge = (entries[idx] >> (orientation * 32)) as u32;
                (edge & (1 << 31) != 0).then_some((
                    (edge & 0x0fff_ffff) as usize,
                    edge & (1 << 30) != 0,
                    u64::from((edge >> 28) & 3),
                ))
            }
        }
    }

    #[inline]
    pub(super) fn outgoing_continuation(
        &self,
        idx: usize,
        forward: bool,
        oriented: u64,
    ) -> Option<(usize, bool, u64)> {
        if let EdgeStorage::Relative6(records) = &self.packed {
            return records.continuation(idx, forward, oriented);
        }
        if let EdgeStorage::Adaptive(records) = &self.packed {
            return records.continuation(idx, forward, oriented);
        }
        if let EdgeStorage::Relative7(records) = &self.packed {
            return records.continuation(idx, forward, oriented);
        }
        if let EdgeStorage::Relative(records) = &self.packed {
            return records.continuation(idx, forward, oriented);
        }
        if let EdgeStorage::Ranked(entries) = &self.packed {
            let record = &entries[idx];
            let orientation = usize::from(!forward);
            let codes = record.codes.load(Ordering::Relaxed);
            let edge = record.edge(orientation);
            return (codes & (1 << (4 + orientation)) != 0).then_some((
                (edge & 0x7fff_ffff) as usize,
                edge & (1 << 31) != 0,
                u64::from((codes >> (orientation * 2)) & 3),
            ));
        }
        if let EdgeStorage::Follow(entries) = &self.packed {
            let edge = (entries[idx] >> (usize::from(!forward) * 32)) as u32;
            return (edge & (1 << 31) != 0).then_some((
                (edge & 0x0fff_ffff) as usize,
                edge & (1 << 30) != 0,
                u64::from((edge >> 28) & 3),
            ));
        }
        self.outgoing(idx, forward, oriented)
            .filter(|&(next, orientation, _)| self.incoming_degree(next, orientation) == 1)
    }

    #[cfg(test)]
    pub(super) fn neighbor(&self, idx: usize, forward: bool, oriented: u64) -> (usize, bool) {
        let orientation = usize::from(!forward);
        match &self.packed {
            EdgeStorage::Adaptive(records) => {
                let (index, forward, _) = records
                    .outgoing(idx, forward, oriented)
                    .expect("single edge");
                (index, forward)
            }
            EdgeStorage::Relative6(records) => {
                let (index, forward, _) = records
                    .outgoing(idx, forward, oriented)
                    .expect("single edge");
                (index, forward)
            }
            EdgeStorage::Relative7(records) => {
                let (index, forward, _) = records
                    .outgoing(idx, forward, oriented)
                    .expect("single edge");
                (index, forward)
            }
            EdgeStorage::Relative(records) => {
                let (index, forward, _) = records
                    .outgoing(idx, forward, oriented)
                    .expect("single edge");
                (index, forward)
            }
            EdgeStorage::Ranked(entries) => {
                let edge = entries[idx].edge(orientation);
                ((edge & 0x7fff_ffff) as usize, edge & (1 << 31) != 0)
            }
            EdgeStorage::Follow(entries) => {
                let edge = (entries[idx] >> (orientation * 32)) as u32;
                ((edge & 0x0fff_ffff) as usize, edge & (1 << 30) != 0)
            }
            EdgeStorage::Wide(entries) => (
                entries[idx].next[orientation] as usize,
                entries[idx].packed & (1 << (2 * EDGE_STATE_BITS + orientation)) != 0,
            ),
            EdgeStorage::Compact(entries) => {
                let edge = (entries[idx] >> (orientation * 32)) as u32;
                ((edge & 0x0fff_ffff) as usize, edge & (1 << 30) != 0)
            }
        }
    }
}

#[cfg(test)]
pub(super) fn compact_cached_edges(cached: CachedEdges) -> [u32; 2] {
    std::array::from_fn(|orientation| {
        let bits =
            (cached.packed >> (orientation * EDGE_STATE_BITS)) & ((1 << EDGE_STATE_BITS) - 1);
        let degree = (bits >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK;
        let code = (bits >> EDGE_OUT_CODE_SHIFT) & 3;
        if degree == 1 {
            debug_assert!(cached.next[orientation] < (1 << 28));
            (1 << 31)
                | cached.next[orientation]
                | (code << 28)
                | (((cached.packed >> (2 * EDGE_STATE_BITS + orientation)) & 1) << 30)
        } else {
            degree | (((bits >> EDGE_HAS_NON_SELF_SHIFT) & 1) << 3) | (code << 28)
        }
    })
}

pub(super) fn cache_edge_states(
    word: u64,
    entries: &[u64],
    lookup: &KmerLookup,
    k: usize,
    mask: u64,
) -> CachedEdges {
    let neighbors = edge_neighbor_words(word, k, mask);
    let found = lookup.find_many(entries, neighbors);
    cached_edges_from_found(word, k, mask, neighbors, found)
}

pub(super) fn compact_pair_word(pair: [u32; 2]) -> u64 {
    u64::from(pair[0]) | (u64::from(pair[1]) << 32)
}

#[inline(always)]
pub(super) fn prefetch_graph_edge(record: &AtomicEdgePair) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: The reference is an in-bounds live atomic record. The store
    // prefetch hint changes no Rust memory or ordering of atomic operations.
    unsafe {
        core::arch::asm!("prfm pstl1keep, [{address}]",address=in(reg) record as *const AtomicEdgePair,options(nostack,readonly,preserves_flags));
    }
    #[cfg(not(target_arch = "aarch64"))]
    let _ = record;
}

#[inline(always)]
pub(super) fn accumulate_graph_edge(
    record: &AtomicU32,
    owner: usize,
    target: usize,
    forward: bool,
    base: u32,
) {
    debug_assert!(target < (1 << 28));
    let first = (1u32 << 31) | (target as u32) | (base << 28) | (u32::from(forward) << 30);
    // Most orientations have one edge. Publish it directly instead of
    // loading the old record and decoding it before the first CAS.
    let mut old = match record.compare_exchange(0, first, Ordering::Relaxed, Ordering::Relaxed) {
        Ok(_) => return,
        Err(actual) => actual,
    };
    loop {
        let single = old & (1 << 31) != 0;
        let seen = if single {
            1 << ((old >> 28) & 3)
        } else {
            (old >> 4) & 15
        };
        let combined = seen | (1 << base);
        if combined == seen {
            return;
        }
        let old_non_self = if single {
            (old & 0x0fff_ffff) as usize != owner
        } else {
            old & 8 != 0
        };
        let non_self = old_non_self || target != owner;
        let merged = combined.count_ones()
            | (u32::from(non_self) << 3)
            | (combined << 4)
            | ((31 - combined.leading_zeros()) << 28);
        match record.compare_exchange_weak(old, merged, Ordering::Relaxed, Ordering::Relaxed) {
            Ok(_) => return,
            Err(actual) => old = actual,
        }
    }
}

pub(super) struct EdgeQueryScratch {
    pub(super) words: [u64; 256],
    pub(super) metadata: [u16; 256],
}
impl Default for EdgeQueryScratch {
    fn default() -> Self {
        Self {
            words: [0; 256],
            metadata: [0; 256],
        }
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(super) fn cache_symmetric_edges(
    first: usize,
    words: &[u64],
    entries: &[u64],
    records: &[AtomicEdgePair],
    lookup: &KmerLookup,
    k: usize,
    mask: u64,
    scratch: &mut EdgeQueryScratch,
) {
    let padding = 64 - 2 * k;
    // Only active slots are read. Keep bounded worker buffers across adjacent
    // batches rather than clearing all inactive entries for every source group.
    let queries = &mut scratch.words;
    let metadata = &mut scratch.metadata;
    let mut count = 0;
    for (relative, &word) in words.iter().enumerate() {
        let reverse = reverse_complement_word(word, k, mask);
        for orientation in 0..2 {
            let (source, reverse_source) = if orientation == 0 {
                (word, reverse)
            } else {
                (reverse, word)
            };
            // Compare the first base before constructing full (k+1)-words.
            // Smaller append codes are strictly canonical; only the boundary
            // code needs the remaining lexicographic tie comparison.
            let boundary = 3 - (source >> 62) as u32;
            for base in 0..=boundary {
                if base == boundary {
                    let edge = source | (u64::from(base) << (padding - 2));
                    let reverse_edge = (u64::from(3 - base) << 62) | (reverse_source >> 2);
                    if edge > reverse_edge {
                        continue;
                    }
                }
                let next = shift_append_word(source, u64::from(base), k, mask);
                let reverse_next = shift_prepend_word(reverse_source, u64::from(3 - base), mask);
                queries[count] = next.min(reverse_next);
                metadata[count] = (base
                    | (orientation << 2)
                    | ((relative as u32) << 3)
                    | (u32::from(next <= reverse_next) << 8))
                    as u16;
                count += 1;
            }
        }
    }
    let found = lookup.find_many_count(entries, *queries, count);
    for target in found[..count].iter().flatten() {
        prefetch_graph_edge(&records[*target]);
    }
    for i in 0..count {
        let Some(target) = found[i] else {
            continue;
        };
        let info = metadata[i];
        let relative = ((info >> 3) & 31) as usize;
        let orientation = ((info >> 2) & 1) as usize;
        let owner = first + relative;
        let target_forward = info & 256 != 0;
        accumulate_graph_edge(
            &records[owner].0[orientation],
            owner,
            target,
            target_forward,
            u32::from(info & 3),
        );
        let reverse_base = if orientation == 0 {
            3 - (words[relative] >> 62) as u32
        } else {
            ((words[relative] >> padding) & 3) as u32
        };
        accumulate_graph_edge(
            &records[target].0[usize::from(target_forward)],
            target,
            owner,
            orientation != 0,
            reverse_base,
        );
    }
}

pub(super) fn edge_neighbor_words(word: u64, k: usize, mask: u64) -> [u64; 8] {
    let reverse = reverse_complement_word(word, k, mask);
    let neighbors: [u64; 8] = std::array::from_fn(|i| {
        let code = (i / 2) as u64;
        if i & 1 == 0 {
            shift_append_word(word, code, k, mask).min(shift_prepend_word(reverse, 3 - code, mask))
        } else {
            shift_prepend_word(word, code, mask).min(shift_append_word(reverse, 3 - code, k, mask))
        }
    });
    neighbors
}

pub(super) fn cached_edges_from_found(
    word: u64,
    k: usize,
    mask: u64,
    neighbors: [u64; 8],
    found: [Option<usize>; 8],
) -> CachedEdges {
    // Decode membership as four-bit masks before selecting neighbour metadata.
    let mut outgoing = 0u32;
    let mut incoming = 0u32;
    let mut non_self_out = 0u32;
    let mut non_self_in = 0u32;
    for code in 0..4 {
        outgoing |= u32::from(found[2 * code].is_some()) << code;
        incoming |= u32::from(found[2 * code + 1].is_some()) << code;
        non_self_out |= u32::from(neighbors[2 * code] != word) << code;
        non_self_in |= u32::from(neighbors[2 * code + 1] != word) << code;
    }
    let highest = |bits: u32| {
        if bits == 0 {
            0
        } else {
            31 - bits.leading_zeros()
        }
    };
    let out_code = highest(outgoing);
    let prev_code = highest(incoming);
    let reverse_code = if incoming == 0 {
        0
    } else {
        3 - incoming.trailing_zeros()
    };
    let next_indices = [
        found[2 * out_code as usize].unwrap_or(0) as u32,
        found[2 * prev_code as usize + 1].unwrap_or(0) as u32,
    ];
    let next = shift_append_word(word, u64::from(out_code), k, mask);
    let prev = shift_prepend_word(word, u64::from(prev_code), mask);
    let orientations = u32::from(outgoing != 0 && next == neighbors[2 * out_code as usize])
        | (u32::from(incoming != 0 && prev != neighbors[2 * prev_code as usize + 1]) << 1);
    let out_degree = outgoing.count_ones();
    let in_degree = incoming.count_ones();
    let forward = in_degree
        | (out_degree << EDGE_OUT_DEGREE_SHIFT)
        | (out_code << EDGE_OUT_CODE_SHIFT)
        | (u32::from(outgoing & non_self_out != 0) << EDGE_HAS_NON_SELF_SHIFT);
    let reverse = out_degree
        | (in_degree << EDGE_OUT_DEGREE_SHIFT)
        | (reverse_code << EDGE_OUT_CODE_SHIFT)
        | (u32::from(incoming & non_self_in != 0) << EDGE_HAS_NON_SELF_SHIFT);
    CachedEdges {
        packed: forward | (reverse << EDGE_STATE_BITS) | (orientations << (2 * EDGE_STATE_BITS)),
        next: next_indices,
    }
}

#[cfg(test)]
pub(super) fn pack_edge_states(
    word: u64,
    entries: &[u64],
    lookup: &KmerLookup,
    k: usize,
    mask: u64,
) -> u32 {
    cache_edge_states(word, entries, lookup, k, mask).packed & ((1 << (2 * EDGE_STATE_BITS)) - 1)
}

impl<'a, K: MappingKeys + ?Sized> GraphTopology<'a, K> {
    pub(super) fn new(
        entries: &'a K,
        lookup: &'a KmerLookup,
        edge_states: &'a EdgeStateCache,
        k: usize,
    ) -> Self {
        let mask = kmer_mask(k);
        let last_shift = 62 - ((k - 1) << 1);
        Self {
            entries,
            lookup,
            edge_states,
            k,
            mask,
            last_shift,
        }
    }

    #[inline]
    pub(super) fn outgoing_neighbor(
        &self,
        oriented: u64,
        idx: usize,
        forward: bool,
    ) -> Option<OrientedNeighbor> {
        let (index, forward, code) = self.edge_states.outgoing(idx, forward, oriented)?;
        Some(OrientedNeighbor {
            oriented: shift_append_word(oriented, code, self.k, self.mask),
            index,
            forward,
        })
    }

    #[inline]
    pub(super) fn incoming_degree(&self, idx: usize, forward: bool) -> u8 {
        self.edge_states.incoming_degree(idx, forward)
    }

    pub(super) fn canonical(&self, oriented: u64) -> u64 {
        oriented.min(self.reverse_complement(oriented))
    }

    pub(super) fn reverse_complement(&self, word: u64) -> u64 {
        reverse_complement_word(word, self.k, self.mask)
    }

    pub(super) fn last_base_code(&self, word: u64) -> u64 {
        (word >> self.last_shift) & 3
    }
}

pub(super) fn reverse_complement_word(word: u64, k: usize, mask: u64) -> u64 {
    #[cfg(target_arch = "aarch64")]
    {
        // Native bit reversal plus an adjacent-bit swap reverses two-bit bases.
        let bits = (!word).reverse_bits();
        let pairs = ((bits & 0x5555_5555_5555_5555) << 1) | ((bits >> 1) & 0x5555_5555_5555_5555);
        (pairs << (64 - 2 * k)) & mask
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        let mut word = (!word) & mask;
        word = ((word & 0x3333_3333_3333_3333) << 2) | ((word >> 2) & 0x3333_3333_3333_3333);
        word = ((word & 0x0f0f_0f0f_0f0f_0f0f) << 4) | ((word >> 4) & 0x0f0f_0f0f_0f0f_0f0f);
        word = ((word & 0x00ff_00ff_00ff_00ff) << 8) | ((word >> 8) & 0x00ff_00ff_00ff_00ff);
        word = ((word & 0x0000_ffff_0000_ffff) << 16) | ((word >> 16) & 0x0000_ffff_0000_ffff);
        word = word.rotate_left(32);
        (word << (64 - (2 * k))) & mask
    }
}
