mod adaptive_prefixes;
mod keys;
mod negative_filter;
mod relative_edges;
mod sparse_locations;
use keys::{MappingKeys, PrefixKeys};
type RelativeEdges = relative_edges::RelativeEdges<8>;
type PackedRelativeEdges = relative_edges::RelativeEdges<7>;
type DenseRelativeEdges = relative_edges::RelativeEdges<6>;
const ADAPTIVE_EDGE_WIDTH: usize = 4;
type AdaptiveRelativeEdges = relative_edges::RelativeEdges<ADAPTIVE_EDGE_WIDTH>;

use super::block_positions::BlockPositions;
use super::packed_transcripts::PackedSequence;
use super::transcript_runs::{TRInfo, TrInfoShard, TranscriptLocations};
use super::workers::IndexWorkers;
use core::range::RangeInclusive;
#[cfg(test)]
use std::collections::HashSet;
use std::collections::{BTreeMap, HashMap};
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::atomic::{AtomicU8, AtomicU32, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use rayon::prelude::*;
use rdst::RadixSort;

use super::bifrost::{encode_minimizer_rep, wyhash};
use super::builder::Transcript;
use crate::minimizer_hash::{REP_HASH_VALS, grouped_state};
use crate::{Error, Result};

// Anchors are bounded candidate maps, never membership proofs. The full
// u64 key and subsequent byte comparison remain authoritative on collisions.
#[derive(Default)]
struct AnchorHasher(u64);
impl Hasher for AnchorHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        self.0 = wyhash(bytes, self.0);
    }
    fn write_u64(&mut self, word: u64) {
        self.0 = wyhash(&word.to_le_bytes(), 0);
    }
}
type AnchorMap<V> = HashMap<u64, V, BuildHasherDefault<AnchorHasher>>;

const INVALID_BASE_CODE: u8 = 4;
const PARALLEL_KMER_COLLECTION_MIN: usize = 1_000_000;
const PARALLEL_TRINFO_COLLECTION_MIN: usize = 1_000_000;
const SWEEP_BLOCK_MIN_BREAKPOINTS: usize = 32;
const MIN_KMER_LOOKUP_PREFIX_BITS: usize = 12;
const MAX_KMER_LOOKUP_PREFIX_BITS: usize = 27;
const TARGET_KMERS_PER_LOOKUP_BUCKET: usize = 1;
const LINEAR_KMER_LOOKUP_MAX: usize = 8;
const EDGE_STATE_BITS: usize = 9;
const EDGE_IN_DEGREE_MASK: u32 = 0x7;
const EDGE_OUT_DEGREE_SHIFT: usize = 3;
const EDGE_OUT_CODE_SHIFT: usize = 6;
const EDGE_HAS_NON_SELF_SHIFT: usize = 8;
const BASE_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'g' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes[b't' as usize] = 3;
    codes
};
const HASH_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'T' as usize] = 2;
    codes[b't' as usize] = 2;
    codes[b'G' as usize] = 3;
    codes[b'g' as usize] = 3;
    codes
};
const HASH_COMP_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'T' as usize] = 0;
    codes[b't' as usize] = 0;
    codes[b'G' as usize] = 1;
    codes[b'g' as usize] = 1;
    codes[b'A' as usize] = 2;
    codes[b'a' as usize] = 2;
    codes[b'C' as usize] = 3;
    codes[b'c' as usize] = 3;
    codes
};

type TrInfoShards = Vec<TrInfoShard>;

type MinimizerIndexParts = (Vec<u64>, Vec<Vec<u32>>);

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct GraphBuildReport {
    pub(super) graph_build: Duration,
    pub(super) minimizer_index: Duration,
    pub(super) ec_build: Duration,
    pub(super) edge_scratch_bytes: usize,
}

#[derive(Clone)]
struct UnitigPath {
    seq: Box<[u8]>,
    ends: [u32; 2],
}

#[derive(Clone, Copy, Default)]
struct KmerLoc(u64);

trait ReferenceKmerLocations: Sync {
    fn find(&self, word: u64) -> Option<KmerLoc>;
}
struct DenseReferenceLocations<'a, K: MappingKeys + ?Sized> {
    keys: &'a K,
    lookup: &'a KmerLookup,
    locations: &'a KmerLocations,
}
impl<K: MappingKeys + ?Sized> ReferenceKmerLocations for DenseReferenceLocations<'_, K> {
    #[inline]
    fn find(&self, word: u64) -> Option<KmerLoc> {
        self.keys
            .find(self.lookup, word)
            .map(|rank| self.locations.get(rank))
    }
}

enum KmerLocations {
    Compact {
        positions: Vec<u32>,
        starts: Vec<u32>,
        buckets: Vec<u32>,
    },
    Wide(Vec<KmerLoc>),
    GlobalCompact {
        positions: Vec<u32>,
        starts: Vec<u32>,
        buckets: Vec<u32>,
    },
}

impl KmerLocations {
    #[inline]
    fn get(&self, index: usize) -> KmerLoc {
        match self {
            Self::Wide(locations) => locations[index],
            Self::Compact {
                positions,
                starts,
                buckets,
            } => {
                let encoded = positions[index];
                if encoded & (1 << 31) == 0 {
                    return KmerLoc((u64::from(encoded >> 9) << 32) | u64::from(encoded & 511));
                }
                let position = (encoded & !(1 << 31)) >> 1;
                let mut unitig = buckets[(position >> 6) as usize] as usize;
                while starts[unitig + 1] <= position {
                    unitig += 1;
                }
                KmerLoc(
                    ((unitig as u64) << 32)
                        | (u64::from(position - starts[unitig]) << 1)
                        | u64::from(encoded & 1),
                )
            }
            Self::GlobalCompact {
                positions,
                starts,
                buckets,
            } => {
                let encoded = positions[index];
                let position = encoded >> 1;
                let mut unitig = buckets[(position >> 6) as usize] as usize;
                while starts[unitig + 1] <= position {
                    unitig += 1;
                }
                KmerLoc(
                    ((unitig as u64) << 32)
                        | (u64::from(position - starts[unitig]) << 1)
                        | u64::from(encoded & 1),
                )
            }
        }
    }
}

#[derive(Clone, Copy)]
struct TrInfoEvent {
    coord: u32,
    is_start: bool,
    trid: u32,
    pos: u32,
}

struct KmerStore {
    sorted: Vec<u64>,
}

#[derive(Clone, Copy)]
struct OrientedNeighbor {
    oriented: u64,
    index: usize,
    forward: bool,
}

// Canonical claims bound the path count by the already-validated u32
// k-mer count. Supported k<32 leaves a low word bit for orientation.
#[repr(C, packed)]
#[derive(Clone, Copy)]
struct OrientedStart {
    word: u64,
    unitig: u32,
}
impl OrientedStart {
    fn word(self) -> u64 {
        self.word & !1
    }
    fn unitig_id(self) -> usize {
        self.unitig as usize
    }
    fn forward(self) -> bool {
        self.word & 1 != 0
    }
}
impl rdst::RadixKey for OrientedStart {
    const LEVELS: usize = 8;
    fn get_level(&self, level: usize) -> u8 {
        (self.word() >> (level * 8)) as u8
    }
}

// Bound exact oriented-start searches to a numeric prefix range. Complete
// words and duplicate checks still establish membership and uniqueness.
struct OrientedStartLookup<'a> {
    starts: &'a [OrientedStart],
    offsets: Vec<u32>,
    shift: u32,
}
impl<'a> OrientedStartLookup<'a> {
    fn new(starts: &'a [OrientedStart], threads: usize) -> Self {
        if starts.len() < 1024 || starts.len() > u32::MAX as usize {
            return Self {
                starts,
                offsets: Vec::new(),
                shift: 0,
            };
        }
        let target = starts.len().div_ceil(8);
        let bits = (usize::BITS - target.saturating_sub(1).leading_zeros()).min(20);
        let shift = 64 - bits;
        let mut offsets = vec![0; (1usize << bits) + 1];
        let fill = |(chunk, dest): (usize, &mut [u32])| {
            let first = chunk * 4096;
            let mut at = starts.partition_point(|start| (start.word() >> shift) < first as u64);
            for (local, slot) in dest.iter_mut().enumerate() {
                let bucket = (first + local) as u64;
                while at < starts.len() && (starts[at].word() >> shift) < bucket {
                    at += 1;
                }
                *slot = at as u32;
            }
        };
        if threads > 1 {
            IndexWorkers::new(threads)
                .expect("existing index worker budget")
                .install(|| offsets.par_chunks_mut(4096).enumerate().for_each(fill));
        } else {
            offsets.chunks_mut(4096).enumerate().for_each(fill);
        }
        Self {
            starts,
            offsets,
            shift,
        }
    }
    fn find(&self, word: u64) -> Option<OrientedStart> {
        let starts = if self.offsets.is_empty() {
            self.starts
        } else {
            let bucket = (word >> self.shift) as usize;
            &self.starts[self.offsets[bucket] as usize..self.offsets[bucket + 1] as usize]
        };
        find_unique_oriented_start(starts, word)
    }
}

enum JoinLinks {
    Narrow(Vec<u32>),
    Wide(Vec<usize>),
}

impl JoinLinks {
    fn retire(self) {
        match self {
            Self::Narrow(values) => super::scratch_pages::retire(values),
            Self::Wide(values) => super::scratch_pages::retire(values),
        }
    }

    fn new(count: usize) -> Self {
        if count <= u32::MAX as usize {
            Self::Narrow(vec![u32::MAX; count])
        } else {
            Self::Wide(vec![usize::MAX; count])
        }
    }
    #[inline]
    fn get(&self, state: usize) -> Option<usize> {
        match self {
            Self::Narrow(values) => (values[state] != u32::MAX).then_some(values[state] as usize),
            Self::Wide(values) => (values[state] != usize::MAX).then_some(values[state]),
        }
    }
    #[inline]
    fn set(&mut self, state: usize, next: usize) {
        match self {
            Self::Narrow(values) => values[state] = next as u32,
            Self::Wide(values) => values[state] = next,
        }
    }
}

enum JoinCandidates {
    Narrow(Vec<[u32; 2]>),
    Wide(Vec<[usize; 2]>),
}

struct GraphTopology<'a, K: MappingKeys + ?Sized = [u64]> {
    entries: &'a K,
    lookup: &'a KmerLookup,
    edge_states: &'a EdgeStateCache,
    k: usize,
    mask: u64,
    last_shift: usize,
}

struct KmerLookup {
    directory: KmerDirectory,
    shift: usize,
    padding: usize,
    suffix_bits: usize,
    suffix_mask: u64,
    rank_bits: usize,
    rank_mask: u64,
    padding_mask: u64,
    pair_bits: usize,
    pair_mask: u64,
    triple_bits: usize,
    triple_mask: u64,
    quad_bits: usize,
    quad_mask: u64,
}

enum KmerDirectory {
    Ranks(RankDirectory),
    Inline(Vec<u64>),
}

struct RankDirectory {
    delta_bits: usize,
    bases: Vec<u32>,
    deltas: Vec<u8>,
    wide: Vec<u32>,
}

impl RankDirectory {
    fn new(offsets: Vec<u32>, parallel: bool) -> Self {
        if offsets.last().copied().unwrap_or(0) >= (1 << 31) {
            return Self {
                delta_bits: 12,
                bases: Vec::new(),
                deltas: Vec::new(),
                wide: offsets,
            };
        }
        if offsets.last().copied().unwrap_or(0) > 1 << 28 {
            Self::packed::<16>(offsets, parallel)
        } else {
            Self::packed::<12>(offsets, parallel)
        }
    }
    fn packed<const BITS: usize>(offsets: Vec<u32>, parallel: bool) -> Self {
        let mut bases = Vec::with_capacity(offsets.len().div_ceil(256));
        let mut deltas = vec![0u8; offsets.len().div_ceil(2) * (BITS / 4) + 1];
        let mut wide = Vec::new();
        for chunk in offsets.chunks(256) {
            let base = chunk[0];
            if chunk[chunk.len() - 1] - base < (1 << BITS) {
                bases.push(base);
            } else {
                bases.push((1 << 31) | wide.len() as u32);
                wide.extend_from_slice(chunk);
            }
        }
        let fill = |(dest, source): (&mut [u8], &[u32])| {
            let base = source[0];
            if source[source.len() - 1] - base < (1 << BITS) {
                for (dest, pair) in dest.chunks_mut(BITS / 4).zip(source.chunks(2)) {
                    let first = pair[0] - base;
                    let second = pair.get(1).map_or(0, |value| value - base);
                    let packed = (first | (second << BITS)).to_le_bytes();
                    dest.copy_from_slice(&packed[..BITS / 4]);
                }
            }
        };
        if parallel {
            deltas
                .par_chunks_mut(256 * BITS / 8)
                .zip(offsets.par_chunks(256))
                .for_each(fill);
        } else {
            deltas
                .chunks_mut(256 * BITS / 8)
                .zip(offsets.chunks(256))
                .for_each(fill);
        }
        // The raw prefix table has been copied exactly into compact blocks.
        // Retire its initialized pages before allocating the edge cache.
        super::scratch_pages::retire(offsets);
        Self {
            delta_bits: BITS,
            bases,
            deltas,
            wide,
        }
    }
    #[inline]
    fn prefetch(&self, bucket: usize) {
        if !self.deltas.is_empty() {
            let first = if self.delta_bits == 16 {
                bucket * 2
            } else {
                (bucket >> 1) * 3
            };
            prefetch_rank_delta(&self.deltas[first]);
            prefetch_rank_delta(&self.bases[bucket >> 8]);
        }
    }
    #[inline]
    fn get(&self, bucket: usize) -> usize {
        if self.bases.is_empty() {
            return self.wide[bucket] as usize;
        }
        let base = self.bases[bucket >> 8];
        if base & (1 << 31) != 0 {
            self.wide[(base & !(1 << 31)) as usize + (bucket & 255)] as usize
        } else {
            let first = if self.delta_bits == 16 {
                bucket * 2
            } else {
                (bucket >> 1) * 3 + (bucket & 1)
            };
            let value = u16::from_le_bytes(
                self.deltas[first..first + 2]
                    .try_into()
                    .expect("packed rank pair"),
            );
            let value = if self.delta_bits == 16 {
                value
            } else {
                (value >> ((bucket & 1) * 4)) & 4095
            };
            (base + u32::from(value)) as usize
        }
    }
    #[inline(always)]
    fn range(&self, bucket: usize) -> (usize, usize) {
        if self.bases.is_empty() {
            let pair = &self.wide[bucket..bucket + 2];
            return (pair[0] as usize, pair[1] as usize);
        }
        if bucket & 255 == 255 {
            return (self.get(bucket), self.get(bucket + 1));
        }
        let base = self.bases[bucket >> 8];
        if base & (1 << 31) != 0 {
            let index = (base & !(1 << 31)) as usize + (bucket & 255);
            let pair = &self.wide[index..index + 2];
            (pair[0] as usize, pair[1] as usize)
        } else {
            let first = if self.delta_bits == 16 {
                bucket * 2
            } else {
                (bucket >> 1) * 3 + (bucket & 1)
            };
            let packed = u32::from_le_bytes(
                self.deltas[first..first + 4]
                    .try_into()
                    .expect("packed rank range"),
            );
            let shift = if self.delta_bits == 16 {
                0
            } else {
                (bucket & 1) * 4
            };
            let mask = (1 << self.delta_bits) - 1;
            (
                (base + ((packed >> shift) & mask)) as usize,
                (base + ((packed >> (self.delta_bits + shift)) & mask)) as usize,
            )
        }
    }
}

struct CachedEdges {
    packed: u32,
    next: [u32; 2],
}

// Large graphs need 31 target bits. Non-singleton edges use unused ranks
// above the graph size for their exact degree; bit31 then stores non-self.
// Keeping both base codes separately preserves the cached graph walk.
#[repr(C, packed)]
struct RankedEdges {
    links: u64,
    codes: AtomicU8,
}

impl RankedEdges {
    fn new(cached: CachedEdges, count: usize) -> Self {
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

    fn build(entries: &[u64], lookup: &KmerLookup, k: usize, parallel: bool) -> Vec<Self> {
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

    fn mark_continuations(records: &[Self], parallel: bool) {
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
    fn edge(&self, orientation: usize) -> u32 {
        (self.links >> (orientation * 32)) as u32
    }

    #[inline]
    fn degree(edge: u32, count: usize) -> u32 {
        let rank = edge & 0x7fff_ffff;
        if rank < count as u32 {
            1
        } else {
            rank - count as u32
        }
    }
}

enum EdgeStorage {
    Compact(Vec<u64>),
    Follow(Vec<u64>),
    Ranked(Vec<RankedEdges>),
    Relative(RelativeEdges),
    Relative7(PackedRelativeEdges),
    Relative6(DenseRelativeEdges),
    Adaptive(AdaptiveRelativeEdges),
    Wide(Vec<CachedEdges>),
}

struct EdgeStateCache {
    packed: EdgeStorage,
}

struct VisitedKmers {
    words: Vec<u64>,
}

pub(super) struct BuildGraph {
    pub(super) k: usize,
    pub(super) g: usize,
    pub(super) unitigs: Vec<Vec<u8>>,
    pub(super) km_unitigs: Vec<Vec<u8>>,
    pub(super) nodes: Vec<NodePayload>,
    pub(super) minimizer_keys: Vec<u64>,
    pub(super) unitig_bitmap_blocks: Vec<Vec<u32>>,
    pub(super) km_bitmap_blocks: Vec<Vec<u32>>,
    pub(super) transcript_names: Vec<String>,
    pub(super) transcript_lengths: Vec<u32>,
}

pub(super) struct NodePayload {
    pub(super) id: u32,
    pub(super) blocks: Vec<NodeBlock>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct NodeBlock {
    pub(super) lb: u32,
    pub(super) ub: u32,
    pub(super) ec: Vec<u32>,
    pub(super) positions: BlockPositions,
}

pub(super) fn build_kmer_unitig_graph_with_report(
    mut transcripts: Vec<Transcript>,
    k: usize,
    g: usize,
    ec_max_size: i32,
    threads: usize,
) -> Result<(BuildGraph, GraphBuildReport)> {
    let graph_start = Instant::now();
    // Build topology first; EC payloads are keyed to the final compacted unitig coordinates.
    // Keep identifier bytes in one owned allocation during topology building,
    // allowing retired sequence/name size-class pages to be reclaimed together.
    let mut identifiers = String::with_capacity(transcripts.iter().map(|tr| tr.name.len()).sum());
    let mut name_ends = Vec::with_capacity(transcripts.len() + 1);
    name_ends.push(0);
    for transcript in &mut transcripts {
        identifiers.push_str(&transcript.name);
        name_ends.push(identifiers.len());
        transcript.name = String::new();
    }
    let (kmers, mut packed) = collect_kmers_and_pack(&mut transcripts, k, threads)?;
    // Collection consumed every ASCII sequence and identifiers already live
    // in their arena. Keep only original lengths through topology building.
    let transcript_lengths = transcripts.iter().map(|tr| tr.original_len).collect();
    super::scratch_pages::retire(transcripts);
    let lookup = KmerLookup::with_threads(&kmers.sorted, threads, k)?;
    let early_prefixes = if PrefixKeys::supports(kmers.sorted.len(), &lookup)
        && kmers.sorted.len() <= (1usize << 31) - 5
    {
        Some(AdaptiveRelativeEdges::adaptive_prefixes(
            &kmers.sorted,
            &lookup,
            k,
        ))
    } else {
        None
    };
    let (edge_scratch_bytes, mapping) = if let Some(starts) = early_prefixes {
        let keys = PrefixKeys::for_topology(kmers.sorted, &lookup, threads)?;
        // Reclaim completed collection/packing scratch on its owners before
        // allocating the temporary membership blocks; live keys stay intact.
        IndexWorkers::new(threads)?.collect_unused(&packed);
        let negative = negative_filter::NegativeFilter::build_mapped(&keys, &lookup, k, threads)?;
        let build = || {
            AdaptiveRelativeEdges::build_filtered(&keys, &lookup, k, starts, &negative, threads > 1)
        };
        let records = if threads > 1 {
            IndexWorkers::new(threads)?.install(build)
        } else {
            build()
        };
        negative.retire();
        IndexWorkers::new(threads)?.collect_unused(&records.words);
        let edges = EdgeStateCache {
            packed: EdgeStorage::Adaptive(records),
        };
        let bytes = edges.scratch_bytes();
        (
            bytes,
            finish_graph_mapping(keys, edges, &lookup, &mut packed, k, threads, graph_start)?,
        )
    } else {
        let edges = EdgeStateCache::build(&kmers.sorted, &lookup, k, threads)?;
        let bytes = edges.scratch_bytes();
        let mapping = if PrefixKeys::supports(kmers.sorted.len(), &lookup) {
            let keys = PrefixKeys::for_topology(kmers.sorted, &lookup, threads)?;
            finish_graph_mapping(keys, edges, &lookup, &mut packed, k, threads, graph_start)?
        } else {
            finish_graph_mapping(
                kmers.sorted,
                edges,
                &lookup,
                &mut packed,
                k,
                threads,
                graph_start,
            )?
        };
        (bytes, mapping)
    };
    let (paths, trinfo_shards, graph_build, ec_start) = mapping;
    // Keys and locations are retired by the mapping helper before EC payloads.
    drop(packed);
    drop(lookup);
    // Mapping is complete. Move identifiers and lengths out, releasing sequence
    // buffers before EC payloads and minimizer records are allocated.
    let copy_name = |bounds: &[usize]| identifiers[bounds[0]..bounds[1]].to_owned();
    let transcript_names = if threads > 1 {
        name_ends.par_windows(2).map(copy_name).collect()
    } else {
        name_ends.windows(2).map(copy_name).collect()
    };
    super::scratch_pages::retire(identifiers.into_bytes());
    super::scratch_pages::retire(name_ends);
    let counts = paths
        .iter()
        .map(|path| (path.seq.len() - k + 1) as u32)
        .collect::<Vec<_>>();
    let trinfos = TranscriptLocations::from_shards(trinfo_shards, &counts);
    drop(counts);
    let ec_threshold = ec_threshold(ec_max_size);
    let nodes = build_node_payloads(&paths, trinfos, ec_threshold, k, threads)?;
    let ec_build = ec_start.elapsed();

    let unitigs: Vec<Vec<u8>> = paths.into_iter().map(|path| path.seq.into_vec()).collect();
    let minimizer_start = Instant::now();
    let (minimizer_keys, unitig_bitmap_blocks) =
        build_unitig_minimizer_index(&unitigs, k, g, threads)?;
    let minimizer_index = minimizer_start.elapsed();

    Ok((
        BuildGraph {
            k,
            g,
            unitigs,
            km_unitigs: Vec::new(),
            nodes,
            minimizer_keys,
            unitig_bitmap_blocks,
            km_bitmap_blocks: vec![Vec::new()],
            transcript_names,
            transcript_lengths,
        },
        GraphBuildReport {
            graph_build,
            minimizer_index,
            ec_build,
            edge_scratch_bytes,
        },
    ))
}

fn finish_graph_mapping<K: MappingKeys>(
    mut keys: K,
    edge_states: EdgeStateCache,
    lookup: &KmerLookup,
    packed: &mut [PackedSequence],
    k: usize,
    threads: usize,
    graph_start: Instant,
) -> Result<(Vec<UnitigPath>, TrInfoShards, Duration, Instant)> {
    let paths = {
        let topology = GraphTopology::new(&keys, lookup, &edge_states, k);
        compact_unitigs(&keys, &topology, threads)?
    };
    // Join probes already carry complete endpoint words and ranks. Only
    // exact suffix membership remains, so rank-to-word hints can retire now.
    keys.retire_word_hints();
    let topology = GraphTopology::new(&keys, lookup, &edge_states, k);
    // Collect only retired worker allocations while the immutable keys remain live.
    IndexWorkers::new(threads)?.collect_unused(&paths);
    let queries = collect_join_queries(&paths, &topology, threads)?;

    edge_states.retire();
    let paths = join_unitig_paths_with_queries(paths, queries, k, threads)?;
    let (trinfos, graph_build, ec_start) =
        map_transcript_runs(keys, packed, &paths, lookup, k, threads, graph_start)?;
    Ok((paths, trinfos, graph_build, ec_start))
}

fn map_transcript_runs<K: MappingKeys>(
    keys: K,
    packed: &mut [PackedSequence],
    paths: &[UnitigPath],
    lookup: &KmerLookup,
    k: usize,
    threads: usize,
    graph_start: Instant,
) -> Result<(TrInfoShards, Duration, Instant)> {
    if sparse_locations::worthwhile(keys.len(), paths.len(), packed.len()) {
        // Full-word landmarks replace both the key directory search and the
        // dense rank-to-location table, so graph-sized keys can retire now.
        drop(keys);
        let locations = sparse_locations::SparseLocations::build(paths, packed, k, threads)?;
        let graph_build = graph_start.elapsed();
        let ec_start = Instant::now();
        let shards = build_trinfos(packed, k, &locations, paths, threads)?;
        return Ok((shards, graph_build, ec_start));
    }
    let kmer_locs = build_kmer_locs(paths, &keys, lookup, k, threads)?;
    let graph_build = graph_start.elapsed();
    let ec_start = Instant::now();
    let locations = DenseReferenceLocations {
        keys: &keys,
        lookup,
        locations: &kmer_locs,
    };
    let shards = build_trinfos(packed, k, &locations, paths, threads)?;
    Ok((shards, graph_build, ec_start))
}

// Keep the private word oracle independent of destructive scratch retirement.
#[cfg(test)]
fn collect_kmers(transcripts: &[Transcript], k: usize, threads: usize) -> Result<KmerStore> {
    collect_kmers_and_pack(&mut transcripts.to_vec(), k, threads).map(|(words, _)| words)
}

fn collect_kmers_and_pack(
    transcripts: &mut [Transcript],
    k: usize,
    threads: usize,
) -> Result<(KmerStore, Vec<PackedSequence>)> {
    let total_positions = transcripts
        .iter()
        .map(|tr| kmer_position_count(tr.seq.len(), k))
        .sum::<usize>();
    if threads > 1 && transcripts.len() > 1 && total_positions >= PARALLEL_KMER_COLLECTION_MIN {
        let pool = IndexWorkers::new(threads)?;
        let mut sorted = Vec::with_capacity(total_positions);
        struct Scratch {
            words: Vec<u64>,
            written: usize,
            packed: Vec<(usize, Vec<PackedSequence>)>,
        }
        let buffers = {
            let collected = std::sync::Mutex::new(&mut sorted);
            // Moving disjoint chunks out of this iterator is safe; workers do
            // all packing only after their group's borrowed anchors are gone.
            let groups = std::sync::Mutex::new(
                bounded_transcript_groups(transcripts, k, 1 << 24).enumerate(),
            );
            pool.install(|| {
                rayon::broadcast(|_| -> Result<Scratch> {
                    let mut words = Vec::new();
                    let mut written = 0usize;
                    let mut packed = Vec::new();
                    let mut prefixes = AnchorMap::default();
                    let mut suffixes = AnchorMap::default();
                    loop {
                        let Some((group, chunk)) = groups.lock().unwrap().next() else {
                            break;
                        };
                        let count = chunk
                            .iter()
                            .map(|tr| kmer_position_count(tr.seq.len(), k))
                            .sum();
                        if words.capacity() < count {
                            super::scratch_pages::discard_unused(&mut words, written);
                            drop(std::mem::take(&mut words));
                            words = Vec::with_capacity(count);
                            written = 0;
                        }
                        prefixes.clear();
                        suffixes.clear();
                        // These values borrow this group. A fresh map bounds their
                        // lifetime before its ASCII sequences are packed and freed.
                        let mut interior = AnchorMap::default();
                        for i in 0..chunk.len() {
                            let sequence = anchored_transcript_region(
                                chunk,
                                i,
                                k,
                                &mut prefixes,
                                &mut suffixes,
                            );
                            insert_interior_unique_kmers::<4>(
                                &mut words,
                                sequence,
                                k,
                                i,
                                &mut interior,
                            )?;
                        }
                        drop(interior);
                        let group_packed = pack_transcript_sequences(chunk, 1)?;
                        if words.len()
                            >= total_positions
                                .div_ceil(threads)
                                .max(PARALLEL_KMER_COLLECTION_MIN)
                        {
                            words.radix_sort_unstable();
                        } else {
                            words
                                .radix_sort_builder()
                                .with_parallel(false)
                                .with_single_threaded_tuner()
                                .sort();
                        }
                        written = written.max(words.len());
                        words.dedup();
                        collected.lock().unwrap().extend_from_slice(&words);
                        packed.push((group, group_packed));
                        words.clear();
                    }
                    Ok(Scratch {
                        words,
                        written,
                        packed,
                    })
                })
            })
        };
        let mut groups = Vec::new();
        for buffer in buffers {
            let Scratch {
                mut words,
                written,
                packed,
            } = buffer?;
            super::scratch_pages::discard_unused(&mut words, written);
            groups.extend(packed);
        }
        groups.sort_unstable_by_key(|(group, _)| *group);
        let packed = concatenate_shards(groups.into_iter().map(|(_, packed)| packed).collect());
        // All original ASCII groups and their borrowed maps are gone now.
        // Reclaim on their owners before the global sorting peak, once only.
        pool.collect_unused(&sorted);
        pool.install(|| sorted.radix_sort_unstable());
        let written = sorted.len();
        sorted.dedup();
        super::scratch_pages::discard_unused(&mut sorted, written);
        return Ok((KmerStore { sorted }, packed));
    }
    let mut sorted = Vec::with_capacity(total_positions);
    for transcript in transcripts.iter() {
        insert_transcript_kmers(&mut sorted, &transcript.seq, k)?;
    }
    sorted.sort_unstable();
    sorted.dedup();
    sorted.shrink_to_fit();
    let packed = pack_transcript_sequences(transcripts, threads)?;
    if sorted.len() >= PARALLEL_KMER_COLLECTION_MIN {
        IndexWorkers::new(threads)?.collect_unused(&sorted);
    }
    Ok((KmerStore { sorted }, packed))
}

#[cfg(test)]
fn novel_transcript_region<'a>(seq: &'a [u8], previous: &[u8], k: usize) -> &'a [u8] {
    &seq[novel_transcript_bounds(seq, previous, k)]
}

#[cfg(test)]
fn novel_transcript_history<'a>(seq: &'a [u8], previous: &[Transcript], k: usize) -> &'a [u8] {
    if seq.len() < k {
        return seq;
    }
    let mut first = 0;
    let mut end = seq.len();
    for prior in previous.iter().rev().take(2) {
        let bounds = novel_transcript_bounds(seq, &prior.seq, k);
        first = first.max(bounds.start);
        end = end.min(bounds.end);
        if first + k > end {
            return &seq[..0];
        }
    }
    &seq[first..end]
}

fn anchored_transcript_region<'a>(
    chunk: &'a [Transcript],
    i: usize,
    k: usize,
    prefixes: &mut AnchorMap<usize>,
    suffixes: &mut AnchorMap<usize>,
) -> &'a [u8] {
    let seq = &chunk[i].seq;
    if seq.len() < k {
        return seq;
    }
    let prefix = if seq.len() >= 32 {
        crate::dna::encode_kmer_pair(&seq[..32]).map(|pair| pair.0)
    } else {
        None
    };
    let suffix = if seq.len() >= 32 {
        crate::dna::encode_kmer_pair(&seq[seq.len() - 32..]).map(|pair| pair.0)
    } else {
        None
    };
    let candidates = [
        i.checked_sub(1),
        prefix.and_then(|key| prefixes.get(&key).copied()),
        suffix.and_then(|key| suffixes.get(&key).copied()),
    ];
    let mut first = 0;
    let mut end = seq.len();
    for (slot, id) in candidates.iter().copied().enumerate() {
        let Some(id) = id else {
            continue;
        };
        if candidates[..slot].contains(&Some(id)) {
            continue;
        }
        let bounds = novel_transcript_bounds(seq, &chunk[id].seq, k);
        first = first.max(bounds.start);
        end = end.min(bounds.end);
    }
    if let Some(key) = prefix {
        prefixes.insert(key, i);
    }
    if let Some(key) = suffix {
        suffixes.insert(key, i);
    }
    if first + k > end {
        &seq[..0]
    } else {
        &seq[first..end]
    }
}

fn shared_prefix_bytes(a: &[u8], b: &[u8]) -> usize {
    let length = a.len().min(b.len());
    let mut matched = 0;
    while matched + 32 <= length && a[matched..matched + 32] == b[matched..matched + 32] {
        matched += 32;
    }
    while matched + 8 <= length {
        let difference = u64::from_le_bytes(a[matched..matched + 8].try_into().unwrap())
            ^ u64::from_le_bytes(b[matched..matched + 8].try_into().unwrap());
        if difference != 0 {
            return matched + difference.trailing_zeros() as usize / 8;
        }
        matched += 8;
    }
    matched
        + a[matched..]
            .iter()
            .zip(&b[matched..])
            .take_while(|(a, b)| a == b)
            .count()
}

fn shared_suffix_bytes(a: &[u8], b: &[u8]) -> usize {
    let length = a.len().min(b.len());
    let mut matched = 0;
    while matched + 32 <= length
        && a[a.len() - matched - 32..a.len() - matched]
            == b[b.len() - matched - 32..b.len() - matched]
    {
        matched += 32;
    }
    while matched + 8 <= length {
        let difference = u64::from_le_bytes(
            a[a.len() - matched - 8..a.len() - matched]
                .try_into()
                .unwrap(),
        ) ^ u64::from_le_bytes(
            b[b.len() - matched - 8..b.len() - matched]
                .try_into()
                .unwrap(),
        );
        if difference != 0 {
            return matched + difference.leading_zeros() as usize / 8;
        }
        matched += 8;
    }
    matched
        + a[..a.len() - matched]
            .iter()
            .rev()
            .zip(b[..b.len() - matched].iter().rev())
            .take_while(|(a, b)| a == b)
            .count()
}

fn novel_transcript_bounds(seq: &[u8], previous: &[u8], k: usize) -> std::ops::Range<usize> {
    if seq.len() < k || previous.len() < k {
        return 0..seq.len();
    }
    let prefix = shared_prefix_bytes(seq, previous);
    if prefix == seq.len() {
        return 0..0;
    }
    let suffix = shared_suffix_bytes(seq, previous);
    let first = prefix.saturating_sub(k - 1);
    let end = kmer_position_count(seq.len(), k).saturating_sub(suffix.saturating_sub(k - 1));
    if first >= end {
        0..0
    } else {
        first..end + k - 1
    }
}

// Anchors choose only candidates. Complete byte equality proves each skipped
// interval; every seam window remains and every transcript retains its metadata.
fn insert_interior_unique_kmers<'a, const ANCHOR_BITS: u32>(
    entries: &mut Vec<u64>,
    seq: &'a [u8],
    k: usize,
    transcript: usize,
    anchors: &mut AnchorMap<(usize, &'a [u8])>,
) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    if seq
        .chunks(32)
        .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none())
    {
        return Err(Error::InvalidFormat("invalid DNA base".into()));
    }
    let mask = kmer_mask(k);
    let shift = 64 - 2 * k;
    let (mut forward, mut reverse) = encode_kmer_pair(&seq[..k])?;
    let mut pos = 0;
    let anchor_mask = (1u64 << ANCHOR_BITS) - 1;
    while pos <= seq.len() - k {
        let hash = (forward ^ (forward >> 27)).wrapping_mul(0x9e3779b97f4a7c15) >> 32;
        if hash & anchor_mask == 0 {
            let (prior_id, prior) = *anchors.entry(forward).or_insert((transcript, &seq[pos..]));
            if prior_id < transcript {
                let current = &seq[pos..];
                let matched = shared_prefix_bytes(current, prior);
                if matched >= k {
                    pos += matched - k + 1;
                    if pos > seq.len() - k {
                        break;
                    }
                    (forward, reverse) = encode_kmer_pair(&seq[pos..pos + k])?;
                    continue;
                }
            }
        }
        entries.push(forward.min(reverse));
        pos += 1;
        if pos > seq.len() - k {
            break;
        }
        let base = seq[pos + k - 1];
        let code = u64::from(((base >> 1) ^ (base >> 2)) & 3);
        forward = ((forward << 2) & mask) | (code << shift);
        reverse = ((reverse >> 2) | ((3 - code) << 62)) & mask;
    }
    Ok(())
}

fn insert_transcript_kmers(entries: &mut Vec<u64>, seq: &[u8], k: usize) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    let mask = kmer_mask(k);
    let shift = 64 - 2 * k;
    let (mut fwd, mut rev) = encode_kmer_pair(&seq[..k])?;
    entries.push(fwd.min(rev));
    // Transcript sequences have already been normalized by the FASTA loader.
    // Validate the tail once, allowing the hot rolling loop to use bit arithmetic.
    if seq[k..]
        .chunks(32)
        .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none())
    {
        return Err(Error::InvalidFormat("invalid DNA base".into()));
    }
    entries.extend(seq[k..].iter().map(|&base| {
        let code = u64::from(((base >> 1) ^ (base >> 2)) & 3);
        fwd = ((fwd << 2) & mask) | (code << shift);
        rev = ((rev >> 2) | ((3 - code) << 62)) & mask;
        fwd.min(rev)
    }));

    Ok(())
}

// Keep ordinary groups at their established target count while bounding raw
// scratch for long transcripts. A single oversized transcript remains intact,
// preserving every seam window; the bound is budget plus one transcript.
fn bounded_transcript_groups(
    transcripts: &mut [Transcript],
    k: usize,
    position_budget: usize,
) -> impl Iterator<Item = &mut [Transcript]> {
    let mut targets = 0usize;
    let mut positions = 0usize;
    transcripts.split_inclusive_mut(move |transcript| {
        targets += 1;
        positions = positions.saturating_add(kmer_position_count(transcript.seq.len(), k));
        if targets >= 2048 || positions >= position_budget {
            targets = 0;
            positions = 0;
            true
        } else {
            false
        }
    })
}

fn kmer_position_count(seq_len: usize, k: usize) -> usize {
    if seq_len >= k { seq_len - k + 1 } else { 0 }
}

fn compact_unitigs(
    entries: &(impl MappingKeys + ?Sized),
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    threads: usize,
) -> Result<Vec<UnitigPath>> {
    if threads > 1 && entries.len() >= PARALLEL_KMER_COLLECTION_MIN {
        return compact_unitigs_parallel(entries, topology, threads);
    }
    let mut visited = VisitedKmers::new(entries.len());
    let mut out = Vec::new();

    // Branch starts produce maximal paths; the second pass catches isolated cycles.
    for idx in 0..entries.len() {
        let word = entries.word(topology.lookup, idx);
        if visited.contains_index(idx) {
            continue;
        }
        let Some(start) = choose_branch_start_orientation(word, idx, topology) else {
            continue;
        };
        let path = walk_unitig(start, idx, topology, &mut visited)?;
        if path.seq.len() >= topology.k {
            out.push(path);
        }
    }

    finish_unitig_cycles(entries, topology, &mut visited, &mut out)?;

    Ok(out)
}

fn finish_unitig_cycles(
    entries: &(impl MappingKeys + ?Sized),
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &mut VisitedKmers,
    out: &mut Vec<UnitigPath>,
) -> Result<()> {
    for word in 0..visited.words.len() {
        loop {
            let remaining = !visited.words[word];
            if remaining == 0 {
                break;
            }
            let idx = (word << 6) + remaining.trailing_zeros() as usize;
            if idx >= entries.len() {
                break;
            }
            let start =
                choose_any_start_orientation(entries.word(topology.lookup, idx), idx, topology);
            out.push(walk_unitig(start, idx, topology, visited)?);
        }
    }
    Ok(())
}

trait PathStep: Copy + Send + Sync {
    fn pack(index: usize, base: u64) -> Self;
    fn index(self) -> usize;
    fn base(self) -> u64;
}
impl PathStep for u32 {
    #[inline]
    fn pack(index: usize, base: u64) -> Self {
        debug_assert!(index < 1 << 30);
        ((index as u32) << 2) | (base as u32)
    }
    #[inline]
    fn index(self) -> usize {
        (self >> 2) as usize
    }
    #[inline]
    fn base(self) -> u64 {
        u64::from(self & 3)
    }
}
impl PathStep for u64 {
    #[inline]
    fn pack(index: usize, base: u64) -> Self {
        ((index as u64) << 2) | base
    }
    #[inline]
    fn index(self) -> usize {
        (self >> 2) as usize
    }
    #[inline]
    fn base(self) -> u64 {
        self & 3
    }
}

struct PreparedPath<Step: PathStep> {
    start: u64,
    steps: Vec<Step>,
}

fn compact_unitigs_parallel(
    entries: &(impl MappingKeys + ?Sized),
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    threads: usize,
) -> Result<Vec<UnitigPath>> {
    if entries.len() <= 1 << 30 {
        compact_unitigs_parallel_steps::<u32>(entries, topology, threads)
    } else {
        compact_unitigs_parallel_steps::<u64>(entries, topology, threads)
    }
}

// The exact combined length prevents geometric growth and repeated prefix
// copies while large ordered worker shards are moved into their final buffer.
fn concatenate_shards<T>(shards: Vec<Vec<T>>) -> Vec<T> {
    let mut combined = Vec::with_capacity(shards.iter().map(Vec::len).sum());
    for shard in shards {
        combined.extend(shard);
    }
    combined
}

// Ordered bounded metadata avoids geometric reallocations retaining large
// prefixes while path sequences are live. Final assembly retires each moved
// chunk's unused pages before filling the next one.
#[derive(Default)]
struct PathChunks {
    chunks: Vec<Vec<UnitigPath>>,
    len: usize,
}
impl PathChunks {
    fn push(&mut self, path: UnitigPath) {
        if self
            .chunks
            .last()
            .is_none_or(|chunk| chunk.len() == chunk.capacity())
        {
            self.chunks.push(Vec::with_capacity(32768));
        }
        self.chunks.last_mut().unwrap().push(path);
        self.len += 1;
    }
    fn extend(&mut self, paths: impl IntoIterator<Item = UnitigPath>) {
        for path in paths {
            self.push(path);
        }
    }
    fn finish(self) -> Vec<UnitigPath> {
        let mut paths = Vec::with_capacity(self.len);
        for mut chunk in self.chunks {
            let written = chunk.len();
            paths.append(&mut chunk);
            super::scratch_pages::discard_unused(&mut chunk, written);
        }
        paths
    }
}

fn compact_unitigs_parallel_steps<Step: PathStep>(
    entries: &(impl MappingKeys + ?Sized),
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    threads: usize,
) -> Result<Vec<UnitigPath>> {
    let pool = IndexWorkers::new(threads)?;
    let starts = pool.install(|| {
        (0..entries.len().div_ceil(65536))
            .into_par_iter()
            .map(|chunk| {
                (chunk * 65536..((chunk + 1) * 65536).min(entries.len()))
                    .filter_map(|idx| {
                        if topology.edge_states.is_internal(idx) {
                            return None;
                        }
                        let canonical = entries.word(topology.lookup, idx);
                        choose_branch_start_orientation(canonical, idx, topology)
                            .map(|word| Step::pack(idx, u64::from(word == canonical)))
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    });
    let starts = concatenate_shards(starts);
    let mut visited = VisitedKmers::new(entries.len());
    let mut out = PathChunks::default();
    let mut buffers = (0..(1024 * threads).min(starts.len()))
        .map(|_| Vec::<Step>::with_capacity(64))
        .collect::<Vec<_>>();
    for batch in starts.chunks(1024 * threads) {
        // Keep only rank and orientation globally; materialize exact start
        // words in the bounded batch needed by the path workers.
        let seeds = batch
            .iter()
            .map(|&seed| {
                let index = seed.index();
                let canonical = entries.word(topology.lookup, index);
                let word = if seed.base() != 0 {
                    canonical
                } else {
                    reverse_complement_word(canonical, topology.k, topology.mask)
                };
                (index, word)
            })
            .collect::<Vec<_>>();
        // Workers see only paths committed by earlier batches. They never claim
        // words; commit in original start order resolves every overlap exactly
        // as the serial walk, including reverse-orientation collisions.
        let paths = pool.install(|| {
            seeds
                .par_chunks(32)
                .zip(buffers.par_chunks_mut(32))
                .map(|(seeds, buffers)| {
                    prepare_path_group::<Step, false>(seeds, buffers, topology, &visited)
                })
                .collect::<Vec<_>>()
        });
        let mut paths = paths;
        // Only the claimed prefix is materialized or reused. Truncate in place
        // after ordered claims rather than retaining a parallel count array.
        for path in paths.iter_mut().flatten().flatten() {
            let count = claim_prepared_path(path, &mut visited);
            path.steps.truncate(count);
        }
        let sequences = pool.install(|| {
            paths
                .par_iter()
                .map(|group| {
                    group
                        .iter()
                        .filter_map(|path| {
                            path.as_ref()
                                .filter(|path| !path.steps.is_empty())
                                .map(|path| {
                                    materialize_prepared_path(path, topology.k, path.steps.len())
                                })
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        });
        out.extend(sequences.into_iter().flatten());
        for (slot, path) in paths.into_iter().flatten().enumerate() {
            let Some(path) = path else {
                continue;
            };
            let mut steps = path.steps;
            if steps.capacity() <= 512 {
                steps.clear();
                buffers[slot] = steps;
            }
        }
    }

    // Unbranched cycles keep the serial ordered rule and orientation choice.
    super::scratch_pages::retire(starts);
    super::scratch_pages::retire(buffers);
    finish_unitig_cycles_parallel::<Step>(
        entries,
        topology,
        &mut visited,
        &mut out,
        &pool,
        threads,
    )?;
    drop(visited);
    Ok(out.finish())
}

fn claim_prepared_path<Step: PathStep>(
    path: &PreparedPath<Step>,
    visited: &mut VisitedKmers,
) -> usize {
    if visited.contains_index(path.steps[0].index()) {
        return 0;
    }
    path.steps
        .iter()
        .take_while(|&&step| visited.insert_index(step.index()))
        .count()
}
fn materialize_prepared_path<Step: PathStep>(
    path: &PreparedPath<Step>,
    k: usize,
    count: usize,
) -> UnitigPath {
    let mut seq = Vec::with_capacity(k + count - 1);
    append_kmer_bases(&mut seq, path.start, k);
    seq.extend(
        path.steps[1..count]
            .iter()
            .map(|step| code_to_base(step.base())),
    );
    UnitigPath {
        seq: seq.into_boxed_slice(),
        ends: [
            path.steps[0].index() as u32,
            path.steps[count - 1].index() as u32,
        ],
    }
}
fn commit_prepared_path<Step: PathStep>(
    path: &PreparedPath<Step>,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &mut VisitedKmers,
) -> Option<UnitigPath> {
    let count = claim_prepared_path(path, visited);
    (count != 0).then(|| materialize_prepared_path(path, topology.k, count))
}

fn finish_unitig_cycles_parallel<Step: PathStep>(
    entries: &(impl MappingKeys + ?Sized),
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &mut VisitedKmers,
    out: &mut PathChunks,
    pool: &IndexWorkers,
    threads: usize,
) -> Result<()> {
    let batch = if entries.len() > 1 << 28 {
        512 * threads
    } else {
        64 * threads
    };
    finish_unitig_cycles_batched::<Step>(entries, topology, visited, out, pool, batch)
}

fn finish_unitig_cycles_batched<Step: PathStep>(
    entries: &(impl MappingKeys + ?Sized),
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &mut VisitedKmers,
    out: &mut PathChunks,
    pool: &IndexWorkers,
    batch: usize,
) -> Result<()> {
    let mut cursor = 0;
    let mut buffers = (0..batch)
        .map(|_| Vec::<Step>::with_capacity(64))
        .collect::<Vec<_>>();
    while cursor < entries.len() {
        let mut seeds = Vec::with_capacity(batch);
        // Select the next unvisited entries in the same ascending order as the
        // serial sweep. Earlier commits can only remove later candidates.
        while cursor < entries.len() && seeds.len() < batch {
            let word = cursor >> 6;
            let remaining = !visited.words[word] & (u64::MAX << (cursor & 63));
            if remaining == 0 {
                cursor = (word + 1) << 6;
                continue;
            }
            let idx = (word << 6) + remaining.trailing_zeros() as usize;
            if idx >= entries.len() {
                cursor = entries.len();
                break;
            }
            cursor = idx + 1;
            seeds.push((
                idx,
                choose_any_start_orientation(entries.word(topology.lookup, idx), idx, topology),
            ));
        }
        let prepared = pool.install(|| {
            seeds
                .par_chunks(32)
                .zip(buffers.par_chunks_mut(32))
                .map(|(seeds, buffers)| {
                    prepare_path_group::<Step, true>(seeds, buffers, topology, visited)
                })
                .collect::<Vec<_>>()
        });
        for (slot, ((idx, start), path)) in seeds
            .into_iter()
            .zip(prepared.into_iter().flatten())
            .enumerate()
        {
            if visited.contains_index(idx) {
                continue;
            }
            if let Some(path) = path {
                if let Some(unitig) = commit_prepared_path(&path, topology, visited) {
                    out.push(unitig);
                }
                let mut steps = path.steps;
                if steps.capacity() <= 512 {
                    steps.clear();
                    buffers[slot] = steps;
                }
            } else {
                // Bound speculative work and memory for exceptionally long
                // components; complete these with the original serial walk.
                out.push(walk_unitig(start, idx, topology, visited)?);
            }
        }
    }
    Ok(())
}

struct PreparingPath<Step: PathStep> {
    path: PreparedPath<Step>,
    start_idx: usize,
    slot: usize,
    index: usize,
    forward: bool,
    base: u64,
}

fn prepare_path_group<Step: PathStep, const RESIDUAL: bool>(
    seeds: &[(usize, u64)],
    buffers: &mut [Vec<Step>],
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &VisitedKmers,
) -> [Option<PreparedPath<Step>>; 32] {
    if let EdgeStorage::Follow(follow) = &topology.edge_states.packed {
        prepare_path_group_inner::<Step, RESIDUAL, true, false, 8>(
            seeds, buffers, topology, visited, follow, None,
        )
    } else if let EdgeStorage::Relative(relative) = &topology.edge_states.packed {
        prepare_path_group_inner::<Step, RESIDUAL, false, true, 8>(
            seeds,
            buffers,
            topology,
            visited,
            &[],
            Some(relative),
        )
    } else if let EdgeStorage::Relative7(relative) = &topology.edge_states.packed {
        prepare_path_group_inner::<Step, RESIDUAL, false, true, 7>(
            seeds,
            buffers,
            topology,
            visited,
            &[],
            Some(relative),
        )
    } else if let EdgeStorage::Adaptive(relative) = &topology.edge_states.packed {
        prepare_path_group_inner::<Step, RESIDUAL, false, true, ADAPTIVE_EDGE_WIDTH>(
            seeds,
            buffers,
            topology,
            visited,
            &[],
            Some(relative),
        )
    } else if let EdgeStorage::Relative6(relative) = &topology.edge_states.packed {
        prepare_path_group_inner::<Step, RESIDUAL, false, true, 6>(
            seeds,
            buffers,
            topology,
            visited,
            &[],
            Some(relative),
        )
    } else {
        prepare_path_group_inner::<Step, RESIDUAL, false, false, 8>(
            seeds,
            buffers,
            topology,
            visited,
            &[],
            None,
        )
    }
}

fn prepare_path_group_inner<
    Step: PathStep,
    const RESIDUAL: bool,
    const FOLLOW: bool,
    const RELATIVE: bool,
    const WIDTH: usize,
>(
    seeds: &[(usize, u64)],
    buffers: &mut [Vec<Step>],
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &VisitedKmers,
    follow: &[u64],
    relative: Option<&relative_edges::RelativeEdges<WIDTH>>,
) -> [Option<PreparedPath<Step>>; 32] {
    let mut output = std::array::from_fn(|_| None);
    let mut pending: [Option<PreparingPath<Step>>; 4] = std::array::from_fn(|_| None);
    let mut cursor = 0;
    let mut oriented = [0u64; 4];
    loop {
        for (lane_id, lane) in pending.iter_mut().enumerate() {
            while lane.is_none() && cursor < seeds.len() {
                let slot = cursor;
                cursor += 1;
                let (idx, start) = seeds[slot];
                if visited.contains_index(idx) {
                    continue;
                }
                if RELATIVE {
                    oriented[lane_id] = start;
                }
                *lane = Some(PreparingPath {
                    path: PreparedPath {
                        start,
                        steps: std::mem::take(&mut buffers[slot]),
                    },
                    start_idx: idx,
                    slot,
                    index: idx,
                    forward: if RELATIVE && WIDTH <= 6 {
                        start <= reverse_complement_word(start, topology.k, topology.mask)
                    } else {
                        start == topology.entries.word(topology.lookup, idx)
                    },
                    base: topology.last_base_code(start),
                });
            }
        }
        if pending.iter().all(Option::is_none) {
            break;
        }
        // Interleave independent streams so outgoing and incoming metadata loads
        // from several paths can be in flight together on each worker.
        let mut next: [Option<(usize, bool, u64)>; 4] = std::array::from_fn(|i| {
            let cur = pending[i].as_mut()?;
            cur.path.steps.push(Step::pack(cur.index, cur.base));
            if FOLLOW {
                let edge = (follow[cur.index] >> (usize::from(!cur.forward) * 32)) as u32;
                (edge & (1 << 31) != 0).then_some((
                    (edge & 0x0fff_ffff) as usize,
                    edge & (1 << 30) != 0,
                    u64::from((edge >> 28) & 3),
                ))
            } else if RELATIVE {
                let relative = relative.expect("relative walk");
                if WIDTH < 8 {
                    let next = relative.outgoing(cur.index, cur.forward, oriented[i]);
                    if let Some((index, forward, _)) = next {
                        relative.prefetch_incoming(index, forward);
                        if RESIDUAL {
                            prefetch_rank_delta(&visited.words[index >> 6]);
                        }
                    }
                    next
                } else {
                    relative.continuation(cur.index, cur.forward, oriented[i])
                }
            } else {
                topology
                    .edge_states
                    .outgoing_continuation(cur.index, cur.forward, 0)
            }
        });
        if RELATIVE && WIDTH < 8 {
            let relative = relative.expect("relative walk");
            for next in &mut next {
                if let Some((index, forward, _)) = *next
                    && relative.incoming_degree(index, forward) != 1
                {
                    *next = None;
                }
            }
        }
        let advance: [bool; 4] = std::array::from_fn(|i| match (pending[i].as_ref(), next[i]) {
            (Some(cur), Some((idx, _, _))) => {
                idx != cur.start_idx && (!RESIDUAL || !visited.contains_index(idx))
            }
            _ => false,
        });
        for i in 0..pending.len() {
            if RESIDUAL
                && advance[i]
                && pending[i]
                    .as_ref()
                    .is_some_and(|cur| cur.path.steps.len() >= 8192)
            {
                pending[i] = None;
            } else if advance[i] {
                let (idx, forward, base) = next[i].unwrap();
                let cur = pending[i].as_mut().unwrap();
                if RELATIVE {
                    oriented[i] = shift_append_word(oriented[i], base, topology.k, topology.mask);
                }
                cur.index = idx;
                cur.forward = forward;
                cur.base = base;
                debug_assert!(cur.path.steps.len() <= 2 * topology.entries.len());
            } else if let Some(cur) = pending[i].take() {
                output[cur.slot] = Some(cur.path);
            }
        }
    }
    output
}

fn choose_branch_start_orientation(
    word: u64,
    idx: usize,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
) -> Option<u64> {
    if topology.edge_states.is_internal(idx) {
        return None;
    }
    let forward_bits = topology.edge_states.state_bits(idx, true, word);
    let reverse_bits = topology.edge_states.state_bits(idx, false, word);
    let starts = |bits: u32| {
        (bits & EDGE_IN_DEGREE_MASK) != 1
            && ((bits >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK) > 0
    };
    let non_self = |bits: u32| bits & (1 << EDGE_HAS_NON_SELF_SHIFT) != 0;
    // Most unvisited words are internal to a path. Their degrees reject them
    // without computing either oriented sequence word.
    match (starts(forward_bits), starts(reverse_bits)) {
        (true, true) if non_self(forward_bits) => Some(word),
        (true, true) if non_self(reverse_bits) => Some(topology.reverse_complement(word)),
        (true, _) => Some(word),
        (_, true) => Some(topology.reverse_complement(word)),
        _ => None,
    }
}

fn choose_any_start_orientation(
    word: u64,
    idx: usize,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
) -> u64 {
    let forward = topology.edge_states.state_bits(idx, true, word);
    let reverse = topology.edge_states.state_bits(idx, false, word);
    if forward & (1 << EDGE_HAS_NON_SELF_SHIFT) != 0 {
        word
    } else if reverse & (1 << EDGE_HAS_NON_SELF_SHIFT) != 0 {
        topology.reverse_complement(word)
    } else if (forward >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK > 0
        || (reverse >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK == 0
    {
        word
    } else {
        topology.reverse_complement(word)
    }
}

fn walk_unitig(
    start: u64,
    start_idx: usize,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &mut VisitedKmers,
) -> Result<UnitigPath> {
    let mut seq = decode_kmer_word(start, topology.k);
    let mut index = start_idx;
    let mut oriented = start;
    let mut forward = start == topology.entries.word(topology.lookup, start_idx);
    loop {
        if !visited.insert_index(index) {
            break;
        }
        let Some((next, next_forward, base)) = topology
            .edge_states
            .outgoing_continuation(index, forward, oriented)
        else {
            break;
        };
        if visited.contains_index(next) {
            break;
        }
        seq.push(code_to_base(base));
        oriented = shift_append_word(oriented, base, topology.k, topology.mask);
        index = next;
        forward = next_forward;
    }
    Ok(UnitigPath {
        seq: seq.into_boxed_slice(),
        ends: [start_idx as u32, index as u32],
    })
}

#[cfg(test)]
fn join_unitig_paths(
    paths: Vec<UnitigPath>,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    threads: usize,
) -> Result<Vec<UnitigPath>> {
    let queries = collect_join_queries(&paths, topology, threads)?;
    join_unitig_paths_with_queries(paths, queries, topology.k, threads)
}

#[repr(C)]
struct JoinQuery {
    word: [u8; 8],
    source: [u8; 5],
}
impl JoinQuery {
    fn new(source: usize, word: u64) -> Result<Self> {
        let source = source as u64;
        if source >= 1 << 40 {
            return Err(Error::InvalidFormat(
                "join source exceeds packed range".into(),
            ));
        }
        Ok(Self {
            word: word.to_ne_bytes(),
            source: source.to_le_bytes()[..5].try_into().unwrap(),
        })
    }
    fn word(&self) -> u64 {
        u64::from_ne_bytes(self.word)
    }
    fn source(&self) -> usize {
        let mut bytes = [0; 8];
        bytes[..5].copy_from_slice(&self.source);
        u64::from_le_bytes(bytes) as usize
    }
}

fn collect_join_queries(
    paths: &[UnitigPath],
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    threads: usize,
) -> Result<Vec<Vec<JoinQuery>>> {
    let collect = |(chunk, paths): (usize, &[UnitigPath])| {
        let mut out = Vec::new();
        for (relative, path) in paths.iter().enumerate() {
            for forward in [true, false] {
                let end = oriented_end_word(&path.seq, forward, topology.k, topology.mask)?;
                let canonical = topology.canonical(end);
                let rank = path.ends[usize::from(forward)];
                let index = if rank == u32::MAX {
                    topology
                        .entries
                        .find(topology.lookup, canonical)
                        .expect("known endpoint")
                } else {
                    let index = rank as usize;
                    debug_assert_eq!(
                        topology.entries.find(topology.lookup, canonical),
                        Some(index)
                    );
                    index
                };
                if let Some(next) = topology.outgoing_neighbor(end, index, end == canonical)
                    && topology.incoming_degree(next.index, next.forward) == 1
                {
                    out.push(JoinQuery::new(
                        state_id(chunk * 512 + relative, forward),
                        next.oriented,
                    )?);
                }
            }
        }
        Ok::<_, Error>(out)
    };
    if threads > 1 && paths.len() >= 1024 {
        IndexWorkers::new(threads)?
            .install(|| paths.par_chunks(512).enumerate().map(collect).collect())
    } else {
        paths.chunks(512).enumerate().map(collect).collect()
    }
}

fn join_unitig_paths_with_queries(
    mut paths: Vec<UnitigPath>,
    queries: Vec<Vec<JoinQuery>>,
    k: usize,
    threads: usize,
) -> Result<Vec<UnitigPath>> {
    if queries.iter().all(Vec::is_empty) {
        return Ok(paths);
    }
    let active = vec![true; paths.len()];
    // Stitch only unique k-1 overlaps left behind after the primary unitig walk.
    let starts = build_oriented_starts(&paths, &active, k, threads)?;
    let candidates = resolve_join_queries(paths.len(), queries, &starts, threads)?;
    super::scratch_pages::retire(starts);
    super::scratch_pages::retire(active);
    let (next, prev) = build_join_edges(paths.len(), candidates);
    let mut used = vec![false; paths.len()];
    let mut out = Vec::with_capacity(paths.len());

    for unitig_id in 0..paths.len() {
        if used[unitig_id] {
            continue;
        }
        let start = rewind_join_start(state_id(unitig_id, true), &prev);
        let seq = joined_sequence_from_state(start, &mut paths, &next, &mut used, k);
        out.push(UnitigPath {
            seq: seq.into_boxed_slice(),
            ends: [u32::MAX; 2],
        });
    }

    super::scratch_pages::retire(paths);
    next.retire();
    prev.retire();
    super::scratch_pages::retire(used);
    Ok(out)
}

fn build_oriented_starts(
    paths: &[UnitigPath],
    active: &[bool],
    k: usize,
    threads: usize,
) -> Result<Vec<OrientedStart>> {
    if threads > 1 && paths.len() >= 1024 {
        let pool = IndexWorkers::new(threads)?;
        return pool.install(|| {
            let mask = kmer_mask(k);
            let counts = paths
                .chunks(512)
                .enumerate()
                .map(|(chunk, paths)| {
                    let first = chunk * 512;
                    active[first..first + paths.len()]
                        .iter()
                        .filter(|&&active| active)
                        .count()
                        * 2
                })
                .collect::<Vec<_>>();
            // Count first and fill disjoint exact slices directly. No full set
            // of shard allocations overlaps with the final sortable start table.
            // Workers initialize the exact spare range directly, avoiding a
            // full serial zero-fill and allowing parallel first touch.
            let total = counts.iter().sum();
            let mut starts = Vec::<OrientedStart>::with_capacity(total);
            let mut remaining = &mut starts.spare_capacity_mut()[..total];
            let mut destinations = Vec::with_capacity(counts.len());
            for count in counts {
                let (destination, tail) = remaining.split_at_mut(count);
                destinations.push(destination);
                remaining = tail;
            }
            destinations
                .into_par_iter()
                .enumerate()
                .try_for_each(|(chunk, destination)| {
                    let first = chunk * 512;
                    let (slots, remainder) = destination.as_chunks_mut::<2>();
                    debug_assert!(remainder.is_empty());
                    let mut slots = slots.iter_mut();
                    for (relative, path) in paths[first..(first + 512).min(paths.len())]
                        .iter()
                        .enumerate()
                    {
                        let unitig_id = first + relative;
                        if active[unitig_id] {
                            let first = encode_kmer(&path.seq[..k])?;
                            let last = encode_kmer(&path.seq[path.seq.len() - k..])?;
                            let slots = slots.next().expect("counted oriented start");
                            slots[0].write(OrientedStart {
                                word: first | 1,
                                unitig: unitig_id as u32,
                            });
                            slots[1].write(OrientedStart {
                                word: reverse_complement_word(last, k, mask),
                                unitig: unitig_id as u32,
                            });
                        }
                    }
                    Ok::<_, Error>(())
                })?;
            // SAFETY: counts assign two slots to every active path; immutable
            // masks and disjoint slices ensure each slot is written once. A
            // successful try_for_each joins every worker before publishing len.
            // On error/panic len stays zero; OrientedStart is Copy with no Drop.
            unsafe {
                starts.set_len(total);
            }
            starts.radix_sort_builder().with_low_mem_tuner().sort();
            pool.collect_unused(&starts);
            Ok(starts)
        });
    }

    let mut starts = Vec::with_capacity(active.iter().filter(|&&is_active| is_active).count() * 2);
    let mask = kmer_mask(k);
    for (unitig_id, path) in paths.iter().enumerate() {
        if !active[unitig_id] {
            continue;
        }
        let first = encode_kmer(&path.seq[..k])?;
        let last = encode_kmer(&path.seq[path.seq.len() - k..])?;
        starts.push(OrientedStart {
            word: first | 1,
            unitig: unitig_id as u32,
        });
        starts.push(OrientedStart {
            word: reverse_complement_word(last, k, mask),
            unitig: unitig_id as u32,
        });
    }
    starts
        .radix_sort_builder()
        .with_parallel(false)
        .with_single_threaded_tuner()
        .sort();
    Ok(starts)
}

fn resolve_join_queries(
    count: usize,
    queries: Vec<Vec<JoinQuery>>,
    starts: &[OrientedStart],
    threads: usize,
) -> Result<JoinCandidates> {
    if count <= (u32::MAX as usize) / 2 {
        resolve_join_queries_typed(queries, starts, threads, |state| state as u32)
            .map(JoinCandidates::Narrow)
    } else {
        resolve_join_queries_typed(queries, starts, threads, std::convert::identity)
            .map(JoinCandidates::Wide)
    }
}
fn resolve_join_queries_typed<State: Copy + Send + Sync>(
    queries: Vec<Vec<JoinQuery>>,
    starts: &[OrientedStart],
    threads: usize,
    encode: impl Fn(usize) -> State + Send + Sync,
) -> Result<Vec<[State; 2]>> {
    let lookup = OrientedStartLookup::new(starts, threads);
    let resolve = |queries: Vec<JoinQuery>| {
        queries
            .into_iter()
            .filter_map(|query| {
                let target = lookup.find(query.word())?;
                (target.unitig_id() != query.source() / 2).then(|| {
                    [
                        encode(query.source()),
                        encode(state_id(target.unitig_id(), target.forward())),
                    ]
                })
            })
            .collect::<Vec<_>>()
    };
    let shards = if threads > 1 {
        IndexWorkers::new(threads)?.install(|| queries.into_par_iter().map(resolve).collect())
    } else {
        queries.into_iter().map(resolve).collect()
    };
    Ok(concatenate_shards(shards))
}

fn build_join_edges(unitig_count: usize, candidates: JoinCandidates) -> (JoinLinks, JoinLinks) {
    let mut next = JoinLinks::new(unitig_count * 2);
    let mut prev = JoinLinks::new(unitig_count * 2);
    let mut insert = |source, target| {
        if next.get(source).is_none() && prev.get(target).is_none() {
            next.set(source, target);
            prev.set(target, source);
        }
    };
    match candidates {
        JoinCandidates::Narrow(values) => {
            for &[source, target] in &values {
                insert(source as usize, target as usize);
            }
            super::scratch_pages::retire(values);
        }
        JoinCandidates::Wide(values) => {
            for &[source, target] in &values {
                insert(source, target);
            }
            super::scratch_pages::retire(values);
        }
    }
    (next, prev)
}

// Valid states occupy0..2*unitig_count. An allocated state array cannot
// reach usize::MAX, leaving an exact out-of-domain missing-link sentinel.
#[inline]
fn join_state(states: &JoinLinks, state: usize) -> Option<usize> {
    states.get(state)
}

fn rewind_join_start(mut state: usize, prev: &JoinLinks) -> usize {
    let start = state;
    let mut slow = state;
    let mut fast = state;
    while let Some(next_slow) = join_state(prev, slow) {
        let Some(next_fast) = join_state(prev, fast).and_then(|value| join_state(prev, value))
        else {
            break;
        };
        slow = next_slow;
        fast = next_fast;
        if slow == fast {
            // Join edges have at most one successor and one predecessor,
            // so a predecessor cycle cannot have a tail entering it.
            return start;
        }
    }
    while let Some(previous) = join_state(prev, state) {
        state = previous;
    }

    state
}

fn joined_sequence_from_state(
    start: usize,
    paths: &mut [UnitigPath],
    next: &JoinLinks,
    used: &mut [bool],
    k: usize,
) -> Vec<u8> {
    let mut state = start;
    let (unitig_id, forward) = state_parts(state);
    let mut seq = std::mem::take(&mut paths[unitig_id].seq).into_vec();
    if !forward {
        seq.reverse();
        for base in &mut seq {
            *base = match *base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
        }
    }

    loop {
        let (unitig_id, _) = state_parts(state);
        used[unitig_id] = true;
        let Some(next_state) = join_state(next, state) else {
            break;
        };
        let (next_unitig_id, next_forward) = state_parts(next_state);
        if used[next_unitig_id] {
            break;
        }
        let next_seq = std::mem::take(&mut paths[next_unitig_id].seq);
        if next_forward {
            seq.extend_from_slice(&next_seq[k - 1..]);
        } else {
            seq.extend(
                next_seq[..next_seq.len() - (k - 1)]
                    .iter()
                    .rev()
                    .map(|&base| code_to_base(3 - u64::from(BASE_CODES[base as usize]))),
            );
        }
        state = next_state;
    }

    seq
}

fn state_id(unitig_id: usize, forward: bool) -> usize {
    unitig_id * 2 + usize::from(!forward)
}

fn state_parts(state: usize) -> (usize, bool) {
    (state / 2, state.is_multiple_of(2))
}

fn find_unique_oriented_start(starts: &[OrientedStart], word: u64) -> Option<OrientedStart> {
    if starts.len() <= 8 {
        let index = starts.iter().position(|start| start.word() == word)?;
        return (index + 1 == starts.len() || starts[index + 1].word() != word)
            .then_some(starts[index]);
    }
    let mut idx = starts
        .binary_search_by_key(&word, |start| start.word())
        .ok()?;
    while idx > 0 && starts[idx - 1].word() == word {
        idx -= 1;
    }
    if idx + 1 < starts.len() && starts[idx + 1].word() == word {
        None
    } else {
        Some(starts[idx])
    }
}

fn oriented_end_word(seq: &[u8], forward: bool, k: usize, mask: u64) -> Result<u64> {
    if forward {
        encode_kmer(&seq[seq.len() - k..])
    } else {
        Ok(reverse_complement_word(encode_kmer(&seq[..k])?, k, mask))
    }
}

fn build_kmer_locs<K: MappingKeys + ?Sized>(
    paths: &[UnitigPath],
    sorted_kmers: &K,
    lookup: &KmerLookup,
    k: usize,
    threads: usize,
) -> Result<KmerLocations> {
    let mut starts = Vec::with_capacity(paths.len() + 1);
    let mut positions = 0u64;
    for path in paths {
        starts.push(positions);
        positions += (path.seq.len() - k + 1) as u64;
    }
    starts.push(positions);
    if positions <= u64::from(u32::MAX >> 2) && paths.len() <= (1 << 22) {
        let starts = starts
            .into_iter()
            .map(|value| value as u32)
            .collect::<Vec<_>>();
        let locations = (0..sorted_kmers.len())
            .map(|_| AtomicU32::new(0))
            .collect::<Vec<_>>();
        let store = |(group, paths): (usize, &[UnitigPath])| {
            store_location_group::<64, _, _, _>(
                paths,
                group * 64,
                k,
                sorted_kmers,
                lookup,
                |index| prefetch_location(&locations[index]),
                |index, loc| {
                    locations[index].store(
                        if loc.offset() < 256 {
                            ((loc.unitig_id() as u32) << 9)
                                | ((loc.offset() as u32) << 1)
                                | u32::from(loc.unitig_is_canonical())
                        } else {
                            (1 << 31)
                                | ((starts[loc.unitig_id()] + loc.offset() as u32) << 1)
                                | u32::from(loc.unitig_is_canonical())
                        },
                        Ordering::Relaxed,
                    );
                },
            )
        };
        if threads > 1 && paths.len() > 1 && sorted_kmers.len() >= PARALLEL_TRINFO_COLLECTION_MIN {
            let pool = IndexWorkers::new(threads)?;
            pool.install(|| paths.par_chunks(64).enumerate().try_for_each(store))?;
        } else {
            for (unitig, path) in paths.iter().enumerate() {
                store_unitig_locations(
                    path,
                    unitig,
                    k,
                    sorted_kmers,
                    lookup,
                    |index| prefetch_location(&locations[index]),
                    |index, loc| {
                        locations[index].store(
                            if loc.offset() < 256 {
                                ((unitig as u32) << 9)
                                    | ((loc.offset() as u32) << 1)
                                    | u32::from(loc.unitig_is_canonical())
                            } else {
                                (1 << 31)
                                    | ((starts[unitig] + loc.offset() as u32) << 1)
                                    | u32::from(loc.unitig_is_canonical())
                            },
                            Ordering::Relaxed,
                        )
                    },
                )?;
            }
        }
        let mut buckets = Vec::with_capacity((positions as usize).div_ceil(64));
        let mut unitig = 0;
        for base in (0..positions as u32).step_by(64) {
            while starts[unitig + 1] <= base {
                unitig += 1;
            }
            buckets.push(unitig as u32);
        }
        return Ok(KmerLocations::Compact {
            positions: locations.into_iter().map(AtomicU32::into_inner).collect(),
            starts,
            buckets,
        });
    }
    // Large graphs can outgrow the direct unitig/offset encoding while still
    // fitting a global position and orientation in one word. Keeping this
    // scratch table narrow avoids random writes to an additional eight-byte
    // table for every k-mer. Smaller graphs retain their direct lookup path.
    if positions <= u64::from(u32::MAX >> 1) {
        let starts = starts
            .into_iter()
            .map(|value| value as u32)
            .collect::<Vec<_>>();
        let locations = (0..sorted_kmers.len())
            .map(|_| AtomicU32::new(0))
            .collect::<Vec<_>>();
        let store = |(group, paths): (usize, &[UnitigPath])| {
            store_location_group::<64, _, _, _>(
                paths,
                group * 64,
                k,
                sorted_kmers,
                lookup,
                |index| prefetch_location(&locations[index]),
                |index, loc| {
                    let position = starts[loc.unitig_id()] + loc.offset() as u32;
                    locations[index].store(
                        (position << 1) | u32::from(loc.unitig_is_canonical()),
                        Ordering::Relaxed,
                    );
                },
            )
        };
        if threads > 1 && paths.len() > 1 && sorted_kmers.len() >= PARALLEL_TRINFO_COLLECTION_MIN {
            IndexWorkers::new(threads)?
                .install(|| paths.par_chunks(64).enumerate().try_for_each(store))?;
        } else {
            for (unitig, path) in paths.iter().enumerate() {
                store_unitig_locations(
                    path,
                    unitig,
                    k,
                    sorted_kmers,
                    lookup,
                    |index| prefetch_location(&locations[index]),
                    |index, loc| {
                        let position = starts[unitig] + loc.offset() as u32;
                        locations[index].store(
                            (position << 1) | u32::from(loc.unitig_is_canonical()),
                            Ordering::Relaxed,
                        );
                    },
                )?;
            }
        }
        let mut buckets = Vec::with_capacity((positions as usize).div_ceil(64));
        let mut unitig = 0;
        for base in (0..positions as u32).step_by(64) {
            while starts[unitig + 1] <= base {
                unitig += 1;
            }
            buckets.push(unitig as u32);
        }
        return Ok(KmerLocations::GlobalCompact {
            positions: locations.into_iter().map(AtomicU32::into_inner).collect(),
            starts,
            buckets,
        });
    }
    if threads > 1 && paths.len() > 1 && sorted_kmers.len() >= PARALLEL_TRINFO_COLLECTION_MIN {
        let out = (0..sorted_kmers.len())
            .map(|_| AtomicU64::new(0))
            .collect::<Vec<_>>();
        let pool = IndexWorkers::new(threads)?;
        pool.install(|| {
            paths
                .par_chunks(64)
                .enumerate()
                .try_for_each(|(group, paths)| {
                    store_location_group::<64, _, _, _>(
                        paths,
                        group * 64,
                        k,
                        sorted_kmers,
                        lookup,
                        |index| prefetch_location(&out[index]),
                        |idx, loc| out[idx].store(loc.0, Ordering::Relaxed),
                    )
                })
        })?;
        return Ok(KmerLocations::Wide(
            out.into_iter()
                .map(|loc| KmerLoc(loc.into_inner()))
                .collect(),
        ));
    }
    let mut out = vec![KmerLoc::default(); sorted_kmers.len()];
    for (unitig_id, path) in paths.iter().enumerate() {
        store_unitig_locations(
            path,
            unitig_id,
            k,
            sorted_kmers,
            lookup,
            |_| {},
            |idx, loc| out[idx] = loc,
        )?;
    }
    Ok(KmerLocations::Wide(out))
}

#[inline(always)]
fn prefetch_location<T>(slot: &T) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: a live in-bounds destination supplies a cache hint only. No
    // memory is written and no alias or atomic ordering is changed.
    unsafe {
        core::arch::asm!("prfm pstl1keep, [{address}]", address=in(reg) slot as *const T, options(nostack,readonly,preserves_flags));
    }
    #[cfg(not(target_arch = "aarch64"))]
    let _ = slot;
}

#[allow(clippy::too_many_arguments)]
fn store_unitig_locations<K: MappingKeys + ?Sized, P, F>(
    path: &UnitigPath,
    unitig_id: usize,
    k: usize,
    sorted_kmers: &K,
    lookup: &KmerLookup,
    prefetch: P,
    mut emit: F,
) -> Result<()>
where
    P: Fn(usize),
    F: FnMut(usize, KmerLoc),
{
    KmerLoc::new(unitig_id, path.seq.len() - k, true)?;
    let mask = kmer_mask(k);
    let mut flush = |batch: [u64; 16], count: usize, first: usize| -> Result<()> {
        let found = sorted_kmers.find_many_count(lookup, batch.map(|word| word & mask), count);
        for &index in found[..count].iter().flatten() {
            prefetch(index);
        }
        for (relative, index) in found.into_iter().enumerate().take(count) {
            let idx = index.ok_or_else(|| {
                Error::InvalidFormat("joined unitig contains unknown k-mer".into())
            })?;
            let loc = KmerLoc(
                ((unitig_id as u64) << 32)
                    | (((first + relative) as u64) << 1)
                    | (batch[relative] & 1),
            );
            emit(idx, loc);
        }
        Ok(())
    };
    let mut batch = [0; 16];
    let mut count = 0;
    let mut first = 0;
    scan_kmer_pairs(&path.seq, k, |_, fwd, rev| {
        batch[count] = fwd.min(rev) | u64::from(fwd <= rev);
        count += 1;
        if count == batch.len() {
            flush(batch, count, first)?;
            first += count;
            count = 0;
        }
        Ok(())
    })?;
    if count != 0 {
        flush(batch, count, first)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn store_location_group<const N: usize, K: MappingKeys + ?Sized, P, F>(
    paths: &[UnitigPath],
    first_unitig: usize,
    k: usize,
    sorted_kmers: &K,
    lookup: &KmerLookup,
    prefetch: P,
    mut emit: F,
) -> Result<()>
where
    P: Fn(usize),
    F: FnMut(usize, KmerLoc),
{
    let mut flush = |keys: [u64; N], locations: &[KmerLoc; N], count: usize| -> Result<()> {
        let found = sorted_kmers.find_many_count(lookup, keys, count);
        for &index in found[..count].iter().flatten() {
            prefetch(index);
        }
        for (slot, index) in found.into_iter().enumerate().take(count) {
            let index = index.ok_or_else(|| {
                Error::InvalidFormat("joined unitig contains unknown k-mer".into())
            })?;
            emit(index, locations[slot]);
        }
        Ok(())
    };
    let mut keys = [0; N];
    let mut locations = [KmerLoc::default(); N];
    let mut count = 0;
    for (relative, path) in paths.iter().enumerate() {
        let unitig_id = first_unitig + relative;
        KmerLoc::new(unitig_id, path.seq.len() - k, true)?;
        scan_kmer_pairs(&path.seq, k, |offset, forward, reverse| {
            keys[count] = forward.min(reverse);
            locations[count] = KmerLoc(
                ((unitig_id as u64) << 32) | ((offset as u64) << 1) | u64::from(forward <= reverse),
            );
            count += 1;
            if count == N {
                flush(keys, &locations, count)?;
                count = 0;
            }
            Ok(())
        })?;
    }
    if count != 0 {
        flush(keys, &locations, count)?;
    }
    Ok(())
}

fn pack_transcript_sequences(
    transcripts: &mut [Transcript],
    threads: usize,
) -> Result<Vec<PackedSequence>> {
    let pack = |group: &mut [Transcript]| {
        let packed =
            PackedSequence::group(group.iter().map(|transcript| transcript.seq.as_slice()))?;
        for transcript in group {
            transcript.seq = Vec::new();
        }
        Ok::<_, Error>(packed)
    };
    let shards = if threads > 1 {
        IndexWorkers::new(threads)?.install(|| {
            transcripts
                .par_chunks_mut(1024)
                .map(pack)
                .collect::<Result<Vec<_>>>()
        })?
    } else {
        transcripts
            .chunks_mut(1024)
            .map(pack)
            .collect::<Result<Vec<_>>>()?
    };
    Ok(concatenate_shards(shards))
}

fn build_trinfos<L: ReferenceKmerLocations + ?Sized>(
    transcripts: &mut [PackedSequence],
    k: usize,
    locations: &L,
    paths: &[UnitigPath],
    threads: usize,
) -> Result<TrInfoShards> {
    let total_positions = transcripts
        .iter()
        .map(|tr| kmer_position_count(tr.len(), k))
        .sum::<usize>();
    if threads > 1
        && transcripts.len() > 1
        && total_positions >= PARALLEL_TRINFO_COLLECTION_MIN
        && !paths.is_empty()
    {
        let pool = IndexWorkers::new(threads)?;
        return pool.install(|| {
            transcripts
                .par_iter_mut()
                .enumerate()
                .map(|(tr_id, transcript)| {
                    let records =
                        collect_transcript_trinfos(tr_id as u32, transcript, k, locations, paths)?;
                    transcript.clear();
                    Ok(records)
                })
                .collect()
        });
    }
    let mut records = TrInfoShard::new(paths.len() <= (1 << 31));
    for (tr_id, transcript) in transcripts.iter_mut().enumerate() {
        scan_transcript_trinfos(
            tr_id as u32,
            transcript,
            k,
            locations,
            paths,
            |unitig_id, info| {
                records.push(unitig_id, info, (paths[unitig_id].seq.len() - k + 1) as u32)
            },
        )?;
        transcript.clear();
    }
    Ok(vec![records])
}

fn collect_transcript_trinfos<L: ReferenceKmerLocations + ?Sized>(
    trid: u32,
    transcript: &PackedSequence,
    k: usize,
    locations: &L,
    paths: &[UnitigPath],
) -> Result<TrInfoShard> {
    let mut out = TrInfoShard::new(paths.len() <= (1 << 31));
    scan_transcript_trinfos(trid, transcript, k, locations, paths, |unitig_id, info| {
        out.push(unitig_id, info, (paths[unitig_id].seq.len() - k + 1) as u32);
    })?;
    Ok(out)
}

// The mapping context is shared by serial and parallel transcript collection.
fn scan_transcript_trinfos<L: ReferenceKmerLocations + ?Sized, F>(
    trid: u32,
    transcript: &PackedSequence,
    k: usize,
    locations: &L,
    paths: &[UnitigPath],
    mut emit: F,
) -> Result<()>
where
    F: FnMut(usize, TRInfo),
{
    let seq = transcript;
    if seq.len() < k {
        return Ok(());
    }
    let mask = kmer_mask(k);
    let mut fwd = seq.word(0, k);
    let mut rev = reverse_complement_word(fwd, k, mask);
    let mut pos = 0;
    loop {
        let mut advance = 1;
        if let Some(loc) = locations.find(fwd.min(rev)) {
            let same_strand = (fwd <= rev) == loc.unitig_is_canonical();
            let unitig = &paths[loc.unitig_id()].seq;
            let offset = loc.offset();
            let remaining = seq.len() - pos - k;
            // Every interior edge of these compacted and joined unitigs has
            // exactly one successor. The matched kmer and strand therefore
            // prove the entire continuation up to the unitig/transcript end.
            // Branches end unitigs; no sampling or unchecked ambiguous jump.
            let matched = if same_strand {
                remaining.min(unitig.len() - offset - k)
            } else {
                remaining.min(offset)
            };
            #[cfg(test)]
            {
                let checked = (0..matched)
                    .take_while(|&i| {
                        let expected = if same_strand {
                            unitig[offset + k + i]
                        } else {
                            code_to_base(3 - u64::from(BASE_CODES[unitig[offset - 1 - i] as usize]))
                        };
                        seq.base(pos + k + i) == expected
                    })
                    .count();
                assert_eq!(
                    matched, checked,
                    "unitig interior must have an exact continuation"
                );
            }

            emit(
                loc.unitig_id(),
                TRInfo {
                    trid,
                    pos: pos as u32 | if same_strand { 0 } else { 0x8000_0000 },
                    start: (if same_strand {
                        offset
                    } else {
                        offset - matched
                    }) as u32,
                    stop: (if same_strand {
                        offset + matched + 1
                    } else {
                        offset + 1
                    }) as u32,
                },
            );
            advance += matched;
        }
        pos += advance;
        if pos > seq.len() - k {
            break;
        }
        if advance == 1 {
            let code = seq.code(pos + k - 1);
            fwd = shift_append_word(fwd, code, k, mask);
            rev = shift_prepend_word(rev, 3 - code, mask);
        } else {
            fwd = seq.word(pos, k);
            rev = reverse_complement_word(fwd, k, mask);
        }
    }
    Ok(())
}

fn build_node_payloads(
    paths: &[UnitigPath],
    mut trinfos: TranscriptLocations,
    ec_threshold: usize,
    k: usize,
    threads: usize,
) -> Result<Vec<NodePayload>> {
    if trinfos.offsets.len() != paths.len() + 1 {
        return Err(Error::InvalidFormat(
            "invalid transcript-location table".into(),
        ));
    }
    let mut tail = trinfos.records.as_mut_slice();
    let mut slices = Vec::with_capacity(paths.len());
    for pair in trinfos.offsets.windows(2) {
        let (current, remaining) = tail.split_at_mut(pair[1] - pair[0]);
        slices.push(current);
        tail = remaining;
    }
    if threads > 1 && paths.len() >= 1024 {
        let pool = IndexWorkers::new(threads)?;
        return Ok(pool.install(|| {
            paths
                .par_iter()
                .zip(slices.into_par_iter())
                .enumerate()
                .map(|(id, (path, infos))| NodePayload {
                    id: id as u32,
                    blocks: build_blocks_for_unitig(
                        infos,
                        (path.seq.len().saturating_sub(k) + 1) as u32,
                        ec_threshold,
                    ),
                })
                .collect()
        }));
    }
    Ok(paths
        .iter()
        .zip(slices)
        .enumerate()
        .map(|(id, (path, infos))| NodePayload {
            id: id as u32,
            blocks: build_blocks_for_unitig(
                infos,
                (path.seq.len().saturating_sub(k) + 1) as u32,
                ec_threshold,
            ),
        })
        .collect())
}

fn build_blocks_for_unitig(
    trinfos: &mut [TRInfo],
    kmer_count: u32,
    ec_threshold: usize,
) -> Vec<NodeBlock> {
    if trinfos.is_empty() {
        return vec![NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec: Vec::new(),
            positions: BlockPositions::default(),
        }];
    }

    if kmer_count > 0
        && trinfos
            .iter()
            .all(|info| info.start == 0 && info.stop == kmer_count)
    {
        // A complete-unitig EC needs one block. Group exact positions directly,
        // reserving each transcript's actual group size rather than a tree.
        trinfos.sort_unstable_by_key(|info| (info.trid, info.pos));
        let count = 1 + trinfos
            .windows(2)
            .filter(|pair| pair[0].trid != pair[1].trid)
            .count();
        let capacity = if count <= ec_threshold { count } else { 0 };
        let (mut ec, mut positions) = (
            Vec::with_capacity(capacity),
            BlockPositions::with_capacity(capacity, if capacity > 0 { trinfos.len() } else { 0 }),
        );
        if count <= ec_threshold {
            let mut first = 0;
            while first < trinfos.len() {
                let end = first
                    + trinfos[first..].partition_point(|info| info.trid == trinfos[first].trid);
                ec.push(trinfos[first].trid);
                positions.push(trinfos[first..end].iter().map(|info| info.pos));
                first = end;
            }
        }
        return vec![NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec,
            positions,
        }];
    }

    let mut breakpoints = Vec::with_capacity(trinfos.len() * 2 + 2);
    breakpoints.push(0);
    breakpoints.push(kmer_count);
    for info in trinfos.iter() {
        breakpoints.push(info.start);
        breakpoints.push(info.stop);
    }
    breakpoints.sort_unstable();
    breakpoints.dedup();

    // Small EC breakpoint sets are cheaper to scan directly; large ones use a sweep-line.
    if breakpoints.len() <= SWEEP_BLOCK_MIN_BREAKPOINTS {
        return build_blocks_for_unitig_by_scan(trinfos, &breakpoints, kmer_count, ec_threshold);
    }

    let mut events = Vec::with_capacity(trinfos.len() * 2);
    for info in trinfos.iter() {
        events.push(TrInfoEvent {
            coord: info.start,
            is_start: true,
            trid: info.trid,
            pos: info.pos,
        });
        events.push(TrInfoEvent {
            coord: info.stop,
            is_start: false,
            trid: info.trid,
            pos: info.pos,
        });
    }
    events.sort_unstable_by_key(|event| event.coord);

    let mut blocks = Vec::new();
    let mut event_idx = 0usize;
    let mut active: BTreeMap<u32, BTreeMap<u32, u32>> = BTreeMap::new();
    for win in breakpoints.windows(2) {
        let lb = win[0];
        let ub = win[1];
        let event_start = event_idx;
        while event_idx < events.len() && events[event_idx].coord == lb {
            event_idx += 1;
        }
        for event in events[event_start..event_idx]
            .iter()
            .filter(|event| !event.is_start)
        {
            remove_active_trinfo(&mut active, event.trid, event.pos);
        }
        for event in events[event_start..event_idx]
            .iter()
            .filter(|event| event.is_start)
        {
            *active
                .entry(event.trid)
                .or_default()
                .entry(event.pos)
                .or_default() += 1;
        }
        if lb == ub {
            continue;
        }
        let mut ec = Vec::new();
        let mut positions = BlockPositions::default();
        if active.len() <= ec_threshold {
            ec.reserve(active.len());
            positions = BlockPositions::with_capacity(
                active.len(),
                active.values().map(BTreeMap::len).sum(),
            );
            for (&trid, values) in &active {
                ec.push(trid);
                positions.push(values.keys().copied());
            }
        }
        push_or_merge_block(
            &mut blocks,
            NodeBlock {
                lb,
                ub,
                ec,
                positions,
            },
        );
    }

    if blocks.is_empty() {
        blocks.push(NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec: Vec::new(),
            positions: BlockPositions::default(),
        });
    }
    blocks
}

fn build_blocks_for_unitig_by_scan(
    trinfos: &mut [TRInfo],
    breakpoints: &[u32],
    kmer_count: u32,
    ec_threshold: usize,
) -> Vec<NodeBlock> {
    trinfos.sort_unstable_by_key(|info| (info.trid, info.pos, info.start, info.stop));
    let mut blocks = Vec::new();
    for win in breakpoints.windows(2) {
        let lb = win[0];
        let ub = win[1];
        if lb == ub {
            continue;
        }
        // Sorted TRInfo already groups transcript IDs and positions. A second
        // tree and per-transcript sort would reconstruct the same order.
        let mut ec = Vec::new();
        let mut positions = BlockPositions::default();
        for info in trinfos.iter() {
            if info.start <= lb && info.stop >= ub {
                if ec.last() != Some(&info.trid) {
                    ec.push(info.trid);
                    positions.start_group();
                }
                positions.push_last_unique(info.pos);
            }
        }
        if ec.len() > ec_threshold {
            ec.clear();
            positions.clear();
        }
        push_or_merge_block(
            &mut blocks,
            NodeBlock {
                lb,
                ub,
                ec,
                positions,
            },
        );
    }

    if blocks.is_empty() {
        blocks.push(NodeBlock {
            lb: 0,
            ub: kmer_count,
            ec: Vec::new(),
            positions: BlockPositions::default(),
        });
    }
    blocks
}

fn remove_active_trinfo(active: &mut BTreeMap<u32, BTreeMap<u32, u32>>, trid: u32, pos: u32) {
    let Some(pos_counts) = active.get_mut(&trid) else {
        return;
    };
    let Some(count) = pos_counts.get_mut(&pos) else {
        return;
    };
    *count -= 1;
    if *count == 0 {
        pos_counts.remove(&pos);
    }
    if pos_counts.is_empty() {
        active.remove(&trid);
    }
}

fn push_or_merge_block(blocks: &mut Vec<NodeBlock>, block: NodeBlock) {
    if let Some(last) = blocks.last_mut()
        && last.ub == block.lb
        && last.ec == block.ec
        && last.positions == block.positions
    {
        last.ub = block.ub;
        return;
    }
    blocks.push(block);
}

fn build_unitig_minimizer_index(
    unitigs: &[Vec<u8>],
    k: usize,
    g: usize,
    threads: usize,
) -> Result<MinimizerIndexParts> {
    let mut keys = Vec::new();
    let total_positions = unitigs.iter().try_fold(0usize, |acc, seq| {
        acc.checked_add(seq.len().saturating_sub(g) + 1)
            .ok_or_else(|| Error::InvalidFormat("too many minimizer positions".into()))
    })?;
    let bitmap_count = (total_positions as u64 >> 32) as usize + 1;
    let mut bitmap_blocks = vec![Vec::new(); bitmap_count.max(1)];
    if threads > 1 && total_positions >= PARALLEL_KMER_COLLECTION_MIN {
        let mut bases = Vec::with_capacity(unitigs.len());
        let mut global = 0u64;
        for seq in unitigs {
            bases.push(global);
            global += minimizer_position_count(seq.len(), g) as u64;
        }
        let pool = IndexWorkers::new(threads)?;
        // Bounded groups amortize allocation without retaining one output per unitig.
        for first in (0..unitigs.len()).step_by(4096 * threads) {
            let end = (first + 4096 * threads).min(unitigs.len());
            let shards = pool.install(|| {
                unitigs[first..end]
                    .par_chunks(512)
                    .enumerate()
                    .map(|(chunk, sequences)| {
                        let mut local_keys = Vec::new();
                        let mut local_blocks = vec![Vec::new(); bitmap_count.max(1)];
                        let mut scratch = MinimizerScratch::default();
                        for (relative, seq) in sequences.iter().enumerate() {
                            insert_selected_minimizers(
                                seq,
                                k,
                                g,
                                bases[first + chunk * 512 + relative],
                                &mut local_keys,
                                &mut local_blocks,
                                &mut scratch,
                            )?;
                        }
                        Ok::<_, Error>((local_keys, local_blocks))
                    })
                    .collect::<Result<Vec<_>>>()
            })?;
            for (local_keys, local_blocks) in shards {
                keys.extend(local_keys);
                for (dest, source) in bitmap_blocks.iter_mut().zip(local_blocks) {
                    dest.extend(source);
                }
            }
        }
    } else {
        let mut global = 0u64;
        let mut scratch = MinimizerScratch::default();
        for seq in unitigs {
            insert_selected_minimizers(
                seq,
                k,
                g,
                global,
                &mut keys,
                &mut bitmap_blocks,
                &mut scratch,
            )?;
            global += minimizer_position_count(seq.len(), g) as u64;
        }
    }

    // Keep exact numeric keys throughout collection and sorting, avoiding a
    // second full-size allocation solely to change element alignment.
    if threads > 1 && keys.len() >= PARALLEL_KMER_COLLECTION_MIN {
        let pool = IndexWorkers::new(threads)?;
        pool.install(|| keys.radix_sort_unstable());
    } else {
        keys.sort_unstable();
    }
    let written = keys.len();
    keys.dedup();
    super::scratch_pages::discard_unused(&mut keys, written);
    Ok((keys, bitmap_blocks))
}

#[derive(Default)]
struct MinimizerScratch {
    window: [(usize, u64); 32],
}

fn insert_selected_minimizers(
    seq: &[u8],
    k: usize,
    g: usize,
    global_base: u64,
    keys: &mut Vec<u64>,
    bitmap_blocks: &mut [Vec<u32>],
    scratch: &mut MinimizerScratch,
) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    let Some(bounds) = bifrost_minimizer_bounds(k, g) else {
        return Ok(());
    };
    if seq
        .chunks(32)
        .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none())
    {
        return Err(Error::InvalidFormat("invalid minimizer".into()));
    }
    let window = &mut scratch.window;
    let mut head = 0usize;
    let mut length = 0usize;
    debug_assert!(bounds.last - bounds.start < 32);
    let (mut forward, mut reverse) =
        minimizer_rep_hash_words(&seq[bounds.start..bounds.start + g])?;
    let remove_forward = REP_HASH_VALS.map(|value| value.rotate_left(g as u32));
    let remove_reverse = REP_HASH_VALS.map(|value| value.rotate_right(1));
    let add_reverse = REP_HASH_VALS.map(|value| value.rotate_left((g - 1) as u32));
    let mut next_pos = bounds.start;
    let mut started = false;
    let mut last_selected = None;
    let mut previous_minimum = None;
    for kmer_start in 0..=seq.len() - k {
        let window_start = kmer_start + bounds.start;
        let window_end = kmer_start + bounds.last;
        while length != 0 && window[head].0 < window_start {
            head = (head + 1) & 31;
            length -= 1;
        }
        while next_pos <= window_end {
            if started {
                let old = seq[next_pos - 1] as usize;
                let new = seq[next_pos + g - 1] as usize;
                forward = forward.rotate_left(1)
                    ^ remove_forward[HASH_CODES[old] as usize]
                    ^ REP_HASH_VALS[HASH_CODES[new] as usize];
                reverse = reverse.rotate_right(1)
                    ^ remove_reverse[HASH_COMP_CODES[old] as usize]
                    ^ add_reverse[HASH_COMP_CODES[new] as usize];
            } else {
                started = true;
            }
            let hash = finish_minimizer_rep_hash(forward, reverse);
            while length != 0 && window[(head + length - 1) & 31].1 > hash {
                length -= 1;
            }
            window[(head + length) & 31] = (next_pos, hash);
            length += 1;
            next_pos += 1;
        }
        if length == 0 {
            continue;
        }
        let minimum = window[head].1;
        // With an unchanged minimum, all old ties were emitted already. Only
        // the one newly admitted position can add another tied minimum.
        let first = if previous_minimum == Some(minimum) {
            length - 1
        } else {
            0
        };
        previous_minimum = Some(minimum);
        for slot in first..length {
            let (pos, hash) = window[(head + slot) & 31];
            if hash != minimum {
                break;
            }
            if last_selected.is_none_or(|last| pos > last) {
                let key = encode_minimizer_rep(&seq[pos..pos + g])
                    .ok_or_else(|| Error::InvalidFormat("invalid minimizer".into()))?;
                add_minimizer_position(global_base + pos as u64, key, keys, bitmap_blocks);
                last_selected = Some(pos);
            }
        }
    }
    Ok(())
}

fn add_minimizer_position(
    pos: u64,
    bytes: [u8; 8],
    keys: &mut Vec<u64>,
    bitmap_blocks: &mut [Vec<u32>],
) {
    keys.push(u64::from_be_bytes(bytes));
    bitmap_blocks[(pos >> 32) as usize].push(pos as u32);
}

fn minimizer_position_count(seq_len: usize, g: usize) -> usize {
    if g == 0 || seq_len < g {
        0
    } else {
        seq_len - g + 1
    }
}

fn bifrost_minimizer_bounds(k: usize, g: usize) -> Option<RangeInclusive<usize>> {
    let shift = 1usize;
    if k < g + shift + 1 {
        return None;
    }
    let end = k.checked_sub(g + shift)?;
    if end < shift {
        None
    } else {
        Some(RangeInclusive {
            start: shift,
            last: end,
        })
    }
}

fn minimizer_rep_hash_words(seq: &[u8]) -> Result<(u64, u64)> {
    if (1..=32).contains(&seq.len()) {
        let (forward, reverse) = crate::dna::encode_kmer_pair(seq)
            .ok_or_else(|| Error::InvalidFormat("invalid minimizer".into()))?;
        return Ok((
            grouped_state(forward, seq.len()),
            grouped_state(reverse, seq.len()),
        ));
    }
    let mut forward = 0u64;
    let mut reverse = 0u64;
    for i in 0..seq.len() {
        let f = HASH_CODES[seq[i] as usize];
        if f == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid minimizer".into()));
        }
        forward = forward.rotate_left(1) ^ REP_HASH_VALS[f as usize];

        let r = HASH_COMP_CODES[seq[seq.len() - 1 - i] as usize];
        if r == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid minimizer".into()));
        }
        reverse = reverse.rotate_left(1) ^ REP_HASH_VALS[r as usize];
    }
    Ok((forward, reverse))
}

fn finish_minimizer_rep_hash(forward: u64, reverse: u64) -> u64 {
    let mut buf = [0u8; 16];
    buf[..8].copy_from_slice(&forward.min(reverse).to_le_bytes());
    buf[8..].copy_from_slice(&forward.max(reverse).to_le_bytes());
    wyhash(&buf, 0)
}

#[cfg(test)]
fn minimizer_rep_hashes_into(seq: &[u8], g: usize, hashes: &mut Vec<u64>) -> Result<()> {
    let (mut forward, mut reverse) = minimizer_rep_hash_words(&seq[..g])?;
    let remove_forward = REP_HASH_VALS.map(|value| value.rotate_left(g as u32));
    let remove_reverse = REP_HASH_VALS.map(|value| value.rotate_right(1));
    let add_reverse = REP_HASH_VALS.map(|value| value.rotate_left((g - 1) as u32));
    hashes.clear();
    hashes.reserve(seq.len() - g + 1);
    hashes.push(finish_minimizer_rep_hash(forward, reverse));
    for pos in g..seq.len() {
        let old = seq[pos - g] as usize;
        let new = seq[pos] as usize;
        let new_forward = HASH_CODES[new];
        if new_forward == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid minimizer".into()));
        }
        // Remove the outgoing base and add the incoming base in both orientations.
        // This is the same REP hash as re-encoding the complete overlapping window.
        forward = forward.rotate_left(1)
            ^ remove_forward[HASH_CODES[old] as usize]
            ^ REP_HASH_VALS[new_forward as usize];
        reverse = reverse.rotate_right(1)
            ^ remove_reverse[HASH_COMP_CODES[old] as usize]
            ^ add_reverse[HASH_COMP_CODES[new] as usize];
        hashes.push(finish_minimizer_rep_hash(forward, reverse));
    }
    Ok(())
}

#[cfg(test)]
fn minimizer_rep_hashes(seq: &[u8], g: usize) -> Result<Vec<u64>> {
    let mut hashes = Vec::new();
    minimizer_rep_hashes_into(seq, g, &mut hashes)?;
    Ok(hashes)
}

#[cfg(test)]
fn minimizer_rep_hash(seq: &[u8]) -> Result<u64> {
    let (forward, reverse) = minimizer_rep_hash_words(seq)?;
    Ok(finish_minimizer_rep_hash(forward, reverse))
}

fn ec_threshold(ec_max_size: i32) -> usize {
    if ec_max_size > 0 {
        ec_max_size as usize
    } else {
        usize::MAX
    }
}

pub(super) fn encode_kmer(seq: &[u8]) -> Result<u64> {
    let (fwd, _) = encode_kmer_pair(seq)?;
    Ok(fwd)
}

fn encode_kmer_pair(seq: &[u8]) -> Result<(u64, u64)> {
    if (1..=32).contains(&seq.len()) {
        let (forward, reverse) = crate::dna::encode_kmer_pair(seq)
            .ok_or_else(|| Error::InvalidFormat("invalid DNA base".into()))?;
        let shift = 64 - 2 * seq.len();
        return Ok((forward << shift, reverse << shift));
    }

    let mut fwd = 0u64;
    let mut rev = 0u64;
    let last = seq.len().saturating_sub(1);
    for (i, &b) in seq.iter().enumerate() {
        let code = BASE_CODES[b as usize];
        if code == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid DNA base".into()));
        }
        let code = u64::from(code);
        fwd |= code << (62 - ((i & 0x1f) << 1));
        rev |= (3 - code) << (62 - (((last - i) & 0x1f) << 1));
    }
    Ok((fwd, rev))
}

fn scan_kmer_pairs<F>(seq: &[u8], k: usize, mut f: F) -> Result<()>
where
    F: FnMut(usize, u64, u64) -> Result<()>,
{
    if seq.len() < k {
        return Ok(());
    }
    let mask = kmer_mask(k);
    let (mut fwd, mut rev) = encode_kmer_pair(&seq[..k])?;
    f(0, fwd, rev)?;
    let mut pos = 1;
    for chunk in seq[k..].chunks(32) {
        if crate::dna::encode_kmer_pair(chunk).is_some() {
            for &base in chunk {
                let code = u64::from(((base >> 1) ^ (base >> 2)) & 3);
                fwd = shift_append_word(fwd, code, k, mask);
                rev = shift_prepend_word(rev, 3 - code, mask);
                f(pos, fwd, rev)?;
                pos += 1;
            }
        } else {
            // Preserve the first invalid-base error and preceding callbacks.
            for &base in chunk {
                let code = base_code(base)?;
                fwd = shift_append_word(fwd, code, k, mask);
                rev = shift_prepend_word(rev, 3 - code, mask);
                f(pos, fwd, rev)?;
                pos += 1;
            }
        }
    }
    Ok(())
}

fn base_code(base: u8) -> Result<u64> {
    let code = BASE_CODES[base as usize];
    if code == INVALID_BASE_CODE {
        return Err(Error::InvalidFormat("invalid DNA base".into()));
    }
    Ok(u64::from(code))
}

fn code_to_base(code: u64) -> u8 {
    match code {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

fn kmer_mask(k: usize) -> u64 {
    u64::MAX << (64 - (2 * k))
}

fn shift_append_word(word: u64, code: u64, k: usize, mask: u64) -> u64 {
    let last_shift = 62 - ((k - 1) << 1);
    ((word << 2) & mask) | (code << last_shift)
}

fn shift_prepend_word(word: u64, code: u64, mask: u64) -> u64 {
    ((word >> 2) | (code << 62)) & mask
}

fn append_kmer_bases(out: &mut Vec<u8>, word: u64, k: usize) {
    #[cfg(all(target_arch = "aarch64", target_endian = "little"))]
    if (16..=32).contains(&k) {
        use std::arch::aarch64::*;
        out.reserve(k);
        let start = out.len();
        // SAFETY: NEON is baseline on AArch64. Both 16-byte writes lie in the
        // reserved k-byte destination; any overlap encodes identical bases.
        // u8 has no invalid values or drop glue, and length is published only
        // after the complete destination is initialized. Input tables are 16B.
        unsafe {
            let indexes = vld1q_u8([3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0].as_ptr());
            let shifts =
                vld1q_s8([-6, -4, -2, 0, -6, -4, -2, 0, -6, -4, -2, 0, -6, -4, -2, 0].as_ptr());
            let alphabet = vld1q_u8(b"ACGTACGTACGTACGT".as_ptr());
            let decode = |packed| {
                let repeated = vqtbl1q_u8(vreinterpretq_u8_u32(vdupq_n_u32(packed)), indexes);
                vqtbl1q_u8(
                    alphabet,
                    vandq_u8(vshlq_u8(repeated, shifts), vdupq_n_u8(3)),
                )
            };
            vst1q_u8(out.as_mut_ptr().add(start), decode((word >> 32) as u32));
            if k > 16 {
                vst1q_u8(
                    out.as_mut_ptr().add(start + k - 16),
                    decode((word << (2 * (k - 16)) >> 32) as u32),
                );
            }
            out.set_len(start + k);
        }
        return;
    }
    out.extend((0..k).map(|i| code_to_base((word >> (62 - 2 * (i & 31))) & 3)));
}

fn decode_kmer_word(word: u64, k: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(k);
    append_kmer_bases(&mut out, word, k);
    out
}

#[cfg(test)]
fn reverse_complement_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

impl KmerLoc {
    fn new(unitig_id: usize, offset: usize, unitig_is_canonical: bool) -> Result<Self> {
        let unitig_id = u32::try_from(unitig_id)
            .map_err(|_| Error::InvalidFormat("too many unitigs for v13 index".into()))?;
        if offset >= (1usize << 31) {
            return Err(Error::InvalidFormat(
                "unitig is too long for packed v13 index metadata".into(),
            ));
        }
        Ok(Self(
            ((unitig_id as u64) << 32) | ((offset as u64) << 1) | u64::from(unitig_is_canonical),
        ))
    }

    fn unitig_id(self) -> usize {
        (self.0 >> 32) as usize
    }

    fn offset(self) -> usize {
        ((self.0 >> 1) & 0x7fff_ffff) as usize
    }

    fn unitig_is_canonical(self) -> bool {
        self.0 & 1 != 0
    }
}

impl VisitedKmers {
    fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
        }
    }

    fn insert_index(&mut self, idx: usize) -> bool {
        let word = &mut self.words[idx >> 6];
        let mask = 1u64 << (idx & 63);
        let previous = *word;
        // Re-setting an existing bit is inert. Publish the same exact claim
        // without a conditional store in the ordered hot loop.
        *word = previous | mask;
        previous & mask == 0
    }

    fn contains_index(&self, idx: usize) -> bool {
        self.words[idx >> 6] & (1u64 << (idx & 63)) != 0
    }
}

#[inline(always)]
fn prefetch_rank_delta<T>(delta: &T) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: the live in-bounds reference supplies a cache hint only.
    unsafe {
        core::arch::asm!("prfm pldl1keep, [{address}]",address=in(reg) delta as *const T,options(nostack,readonly,preserves_flags));
    }
    #[cfg(not(target_arch = "aarch64"))]
    let _ = delta;
}

#[inline(always)]
fn prefetch_kmer(word: &u64) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: The reference identifies an in-bounds live key. PRFM is a
    // nonblocking cache hint; it neither reads a Rust value nor changes memory.
    unsafe {
        core::arch::asm!("prfm pldl1keep, [{address}]",address=in(reg) word as *const u64,options(nostack,readonly,preserves_flags));
    }
    #[cfg(not(target_arch = "aarch64"))]
    let _ = word;
}

impl KmerLookup {
    #[cfg(test)]
    fn new(entries: &[u64]) -> Result<Self> {
        Self::with_threads(entries, 1, 32)
    }

    fn with_threads(entries: &[u64], threads: usize, k: usize) -> Result<Self> {
        if entries.len() > u32::MAX as usize {
            return Err(Error::InvalidFormat(
                "too many distinct k-mers for compact lookup".into(),
            ));
        }
        let target = entries.len().div_ceil(2).max(1);
        let bits = (usize::BITS as usize - target.saturating_sub(1).leading_zeros() as usize)
            .clamp(MIN_KMER_LOOKUP_PREFIX_BITS, MAX_KMER_LOOKUP_PREFIX_BITS);
        let suffix_bits = (2 * k).saturating_sub(bits);
        let padding = 64 - 2 * k;
        let rank_bits =
            usize::BITS as usize - entries.len().saturating_sub(1).leading_zeros() as usize;
        // Keep the directory within512MiB and leave the descriptor's high bit
        // free to distinguish a complete singleton key from a rank range.
        let inline = rank_bits + suffix_bits <= 63
            && entries.len() < (1 << 27)
            && (1usize << bits) <= 1 << 26;
        let prefix = if inline {
            bits
        } else {
            kmer_lookup_prefix_bits(entries.len())
        };
        let shift = 64 - prefix;
        let count = 1usize << prefix;
        let suffix_mask = if suffix_bits == 0 {
            0
        } else {
            (1u64 << suffix_bits) - 1
        };
        let pair_bits = (61 - rank_bits) / 2;
        let pair_mask = (1u64 << pair_bits) - 1;
        let triple_bits = (60 - rank_bits) / 3;
        let triple_mask = (1u64 << triple_bits) - 1;
        let quad_bits = (59 - rank_bits) / 4;
        let quad_mask = (1u64 << quad_bits) - 1;
        let pool = if threads > 1 && entries.len() >= PARALLEL_KMER_COLLECTION_MIN {
            Some(IndexWorkers::new(threads)?)
        } else {
            None
        };
        let directory = if inline {
            let mut slots = vec![0u64; count];
            let fill = |(chunk, dest): (usize, &mut [u64])| {
                let first = chunk << 12;
                let mut index = entries.partition_point(|&word| (word >> shift) < first as u64);
                for (local, slot) in dest.iter_mut().enumerate() {
                    let start = index;
                    let end_bucket = (first + local + 1) as u64;
                    while index < entries.len() && (entries[index] >> shift) < end_bucket {
                        index += 1;
                    }
                    *slot = if index == start {
                        0
                    } else if index == start + 1 {
                        (1u64 << 63)
                            | ((start as u64) << suffix_bits)
                            | ((entries[start] >> padding) & suffix_mask)
                    } else if index == start + 2 {
                        let first = (entries[start] >> padding) & pair_mask;
                        let second = (entries[start + 1] >> padding) & pair_mask;
                        (1u64 << 61)
                            | start as u64
                            | (first << rank_bits)
                            | (second << (rank_bits + pair_bits))
                    } else if index == start + 3 {
                        let mut descriptor = (1u64 << 60) | start as u64;
                        for slot in 0..3 {
                            descriptor |= ((entries[start + slot] >> padding) & triple_mask)
                                << (rank_bits + slot * triple_bits);
                        }
                        descriptor
                    } else if index == start + 4 {
                        let mut descriptor = (1u64 << 59) | start as u64;
                        for slot in 0..4 {
                            descriptor |= ((entries[start + slot] >> padding) & quad_mask)
                                << (rank_bits + slot * quad_bits);
                        }
                        descriptor
                    } else {
                        ((index as u64) << 32) | start as u64
                    };
                }
            };
            if let Some(pool) = &pool {
                pool.install(|| slots.par_chunks_mut(1 << 12).enumerate().for_each(fill));
            } else {
                slots.chunks_mut(1 << 12).enumerate().for_each(fill);
            }
            KmerDirectory::Inline(slots)
        } else {
            let mut offsets = vec![0u32; count + 1];
            let fill = |(chunk, dest): (usize, &mut [u32])| {
                let first = chunk << 18;
                let mut index = entries.partition_point(|&word| (word >> shift) < first as u64);
                for (local, slot) in dest.iter_mut().enumerate() {
                    let bucket = (first + local) as u64;
                    while index < entries.len() && (entries[index] >> shift) < bucket {
                        index += 1;
                    }
                    *slot = index as u32;
                }
            };
            if let Some(pool) = &pool {
                pool.install(|| offsets.par_chunks_mut(1 << 18).enumerate().for_each(fill));
            } else {
                offsets.chunks_mut(1 << 18).enumerate().for_each(fill);
            }
            KmerDirectory::Ranks(if let Some(pool) = &pool {
                pool.install(|| RankDirectory::new(offsets, true))
            } else {
                RankDirectory::new(offsets, false)
            })
        };
        Ok(Self {
            directory,
            shift,
            padding,
            suffix_bits,
            suffix_mask,
            rank_bits,
            rank_mask: (1u64 << rank_bits) - 1,
            padding_mask: (1u64 << padding) - 1,
            pair_bits,
            pair_mask,
            triple_bits,
            triple_mask,
            quad_bits,
            quad_mask,
        })
    }

    #[inline]
    fn search(entries: &[u64], word: u64, start: usize, end: usize) -> Option<usize> {
        let slice = &entries[start..end];
        if slice.len() == 1 {
            return (slice[0] == word).then_some(start);
        }
        if slice.len() == 2 {
            return if slice[0] == word {
                Some(start)
            } else {
                (slice[1] == word).then_some(start + 1)
            };
        }
        if slice.len() <= LINEAR_KMER_LOOKUP_MAX {
            slice
                .iter()
                .position(|&entry| entry == word)
                .map(|idx| start + idx)
        } else {
            slice.binary_search(&word).ok().map(|idx| start + idx)
        }
    }

    #[inline]
    fn guarded_result<const N: usize>(
        &self,
        entries: &[u64],
        word: u64,
        descriptor: u64,
        bits: usize,
        mask: u64,
    ) -> Option<usize> {
        let start = (descriptor & self.rank_mask) as usize;
        let guard = (word >> self.padding) & mask;
        for slot in 0..N {
            if guard == (descriptor >> (self.rank_bits + slot * bits)) & mask
                && entries[start + slot] == word
            {
                return Some(start + slot);
            }
        }
        None
    }

    #[inline(always)]
    fn inline_result(&self, entries: &[u64], word: u64, descriptor: u64) -> Option<usize> {
        if descriptor == 0 {
            return None;
        }
        match descriptor.leading_zeros() {
            0 => {
                if word & self.padding_mask != 0
                    || (word >> self.padding) & self.suffix_mask != descriptor & self.suffix_mask
                {
                    return None;
                }
                Some(((descriptor & !(1 << 63)) >> self.suffix_bits) as usize)
            }
            2 => {
                let start = (descriptor & self.rank_mask) as usize;
                let fingerprint = (word >> self.padding) & self.pair_mask;
                if fingerprint == (descriptor >> self.rank_bits) & self.pair_mask
                    && entries[start] == word
                {
                    return Some(start);
                }
                if fingerprint == (descriptor >> (self.rank_bits + self.pair_bits)) & self.pair_mask
                    && entries[start + 1] == word
                {
                    return Some(start + 1);
                }
                None
            }
            3 => self.guarded_result::<3>(
                entries,
                word,
                descriptor,
                self.triple_bits,
                self.triple_mask,
            ),
            4 => {
                self.guarded_result::<4>(entries, word, descriptor, self.quad_bits, self.quad_mask)
            }
            _ => Self::search(
                entries,
                word,
                (descriptor as u32) as usize,
                (descriptor >> 32) as usize,
            ),
        }
    }

    #[inline]
    fn find_many<const N: usize>(&self, entries: &[u64], words: [u64; N]) -> [Option<usize>; N] {
        self.find_many_count(entries, words, N)
    }

    #[inline(always)]
    fn find_many_count<const N: usize>(
        &self,
        entries: &[u64],
        words: [u64; N],
        count: usize,
    ) -> [Option<usize>; N] {
        debug_assert!(count <= N);
        let mut found = [None; N];
        match &self.directory {
            KmerDirectory::Inline(slots) => {
                let mut descriptors = [0u64; N];
                for i in 0..count {
                    descriptors[i] = slots[(words[i] >> self.shift) as usize];
                }
                for i in 0..count {
                    found[i] = self.inline_result(entries, words[i], descriptors[i]);
                }
            }
            KmerDirectory::Ranks(offsets) => {
                for &word in &words[..count] {
                    offsets.prefetch((word >> self.shift) as usize);
                }
                let mut bounds = [(0usize, 0usize); N];
                for i in 0..count {
                    let bucket = (words[i] >> self.shift) as usize;
                    bounds[i] = offsets.range(bucket);
                }
                // Hint independent exact key loads before searching each range.
                // Empty buckets issue no hint; large buckets start at midpoint.
                for &(first, end) in &bounds[..count] {
                    if first < end {
                        let index = if end - first > LINEAR_KMER_LOOKUP_MAX {
                            first + (end - first) / 2
                        } else {
                            first
                        };
                        prefetch_kmer(&entries[index]);
                    }
                }
                for i in 0..count {
                    found[i] = Self::search(entries, words[i], bounds[i].0, bounds[i].1);
                }
            }
        }
        found
    }

    fn find(&self, entries: &[u64], word: u64) -> Option<usize> {
        let bucket = (word >> self.shift) as usize;
        match &self.directory {
            KmerDirectory::Inline(slots) => self.inline_result(entries, word, slots[bucket]),
            KmerDirectory::Ranks(offsets) => {
                let (start, end) = offsets.range(bucket);
                Self::search(entries, word, start, end)
            }
        }
    }
}

fn kmer_lookup_prefix_bits(entry_count: usize) -> usize {
    if entry_count == 0 {
        return MIN_KMER_LOOKUP_PREFIX_BITS;
    }
    let target_buckets = entry_count.div_ceil(TARGET_KMERS_PER_LOOKUP_BUCKET);
    let bits = usize::BITS as usize - target_buckets.saturating_sub(1).leading_zeros() as usize;
    bits.clamp(MIN_KMER_LOOKUP_PREFIX_BITS, MAX_KMER_LOOKUP_PREFIX_BITS)
}

// Align like the final u64 table so final collection can reuse its allocation.
#[repr(C, align(8))]
struct AtomicEdgePair([AtomicU32; 2]);

// Reuse the atomic allocation when its size and alignment match immutable
// u64 words; other allocator layouts retain the ordinary collection fallback. All discovery workers have joined
// before exclusive finalization; no atomics or aliases survive this conversion.
fn finish_atomic_edge_records(mut records: Vec<AtomicEdgePair>, parallel: bool) -> EdgeStorage {
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
fn follow_degree(record: u32, count: usize) -> u32 {
    let rank = record & 0x0fff_ffff;
    if rank < count as u32 {
        1
    } else {
        rank - count as u32
    }
}

impl EdgeStateCache {
    fn scratch_bytes(&self) -> usize {
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

    fn retire(self) {
        match self.packed {
            EdgeStorage::Relative6(records) => super::scratch_pages::retire(records.words),
            EdgeStorage::Adaptive(records) => super::scratch_pages::retire(records.words),
            EdgeStorage::Relative7(records) => super::scratch_pages::retire(records.words),
            EdgeStorage::Relative(records) => super::scratch_pages::retire(records.words),
            EdgeStorage::Compact(records) | EdgeStorage::Follow(records) => {
                super::scratch_pages::retire(records)
            }
            EdgeStorage::Ranked(records) => super::scratch_pages::retire(records),
            EdgeStorage::Wide(records) => super::scratch_pages::retire(records),
        }
    }

    fn build(entries: &[u64], lookup: &KmerLookup, k: usize, threads: usize) -> Result<Self> {
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
    fn is_internal(&self, idx: usize) -> bool {
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

    fn state_bits(&self, idx: usize, forward: bool, canonical: u64) -> u32 {
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
    fn incoming_degree(&self, idx: usize, forward: bool) -> u8 {
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
    fn outgoing(&self, idx: usize, forward: bool, oriented: u64) -> Option<(usize, bool, u64)> {
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
    fn outgoing_continuation(
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
    fn neighbor(&self, idx: usize, forward: bool, oriented: u64) -> (usize, bool) {
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
fn compact_cached_edges(cached: CachedEdges) -> [u32; 2] {
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

fn cache_edge_states(
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

fn compact_pair_word(pair: [u32; 2]) -> u64 {
    u64::from(pair[0]) | (u64::from(pair[1]) << 32)
}

#[inline(always)]
fn prefetch_graph_edge(record: &AtomicEdgePair) {
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
fn accumulate_graph_edge(
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

struct EdgeQueryScratch {
    words: [u64; 256],
    metadata: [u16; 256],
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
fn cache_symmetric_edges(
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

fn edge_neighbor_words(word: u64, k: usize, mask: u64) -> [u64; 8] {
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

fn cached_edges_from_found(
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
fn pack_edge_states(word: u64, entries: &[u64], lookup: &KmerLookup, k: usize, mask: u64) -> u32 {
    cache_edge_states(word, entries, lookup, k, mask).packed & ((1 << (2 * EDGE_STATE_BITS)) - 1)
}

impl<'a, K: MappingKeys + ?Sized> GraphTopology<'a, K> {
    fn new(
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
    fn outgoing_neighbor(
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
    fn incoming_degree(&self, idx: usize, forward: bool) -> u8 {
        self.edge_states.incoming_degree(idx, forward)
    }

    fn canonical(&self, oriented: u64) -> u64 {
        oriented.min(self.reverse_complement(oriented))
    }

    fn reverse_complement(&self, word: u64) -> u64 {
        reverse_complement_word(word, self.k, self.mask)
    }

    fn last_base_code(&self, word: u64) -> u64 {
        (word >> self.last_shift) & 3
    }
}

fn reverse_complement_word(word: u64, k: usize, mask: u64) -> u64 {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn reverse_complement_word_slow(word: u64, k: usize) -> u64 {
        let mut out = 0u64;
        for i in 0..k {
            let shift = 62 - ((i & 0x1f) << 1);
            let code = (word >> shift) & 3;
            let out_shift = 62 - (((k - 1 - i) & 0x1f) << 1);
            out |= (3 - code) << out_shift;
        }
        out
    }

    #[test]
    fn oriented_start_prefix_ranges_preserve_unique_duplicate_and_missing_words() {
        let mut records = Vec::new();
        for i in 0..2048u64 {
            for copy in 0..=i % 3 {
                records.push(OrientedStart {
                    word: (i << 52) | (copy & 1),
                    unitig: (i * 3 + copy) as u32,
                });
            }
        }
        let identity = |v: Option<OrientedStart>| v.map(|v| (v.word(), v.unitig_id(), v.forward()));
        let mut expected = HashMap::<u64, Option<OrientedStart>>::new();
        for &record in &records {
            expected
                .entry(record.word())
                .and_modify(|entry| *entry = None)
                .or_insert(Some(record));
        }
        for threads in [1, 8] {
            let lookup = OrientedStartLookup::new(&records, threads);
            for i in 0..4096u64 {
                for word in [i << 52, (i << 52) | (1 << 40)] {
                    assert_eq!(
                        identity(lookup.find(word)),
                        identity(expected.get(&word).copied().flatten())
                    );
                }
            }
        }
        for records in [&records[..0], &records[..1], &records[..17]] {
            let lookup = OrientedStartLookup::new(records, 1);
            for word in [0, 1 << 52, u64::MAX] {
                assert_eq!(
                    identity(lookup.find(word)),
                    identity(find_unique_oriented_start(records, word))
                );
            }
        }
    }

    #[test]
    fn shared_transcript_ends_preserve_every_seam_window() {
        // Exact differing-byte positions are the proof boundary, including
        // shorter slices and every 8/32-byte seam in either direction.
        for length in 0usize..97 {
            let a = (0..length).map(|i| (i * 37) as u8).collect::<Vec<_>>();
            for other_length in [length, length.saturating_sub(3)] {
                for change in 0..=other_length {
                    let mut b = a[..other_length].to_vec();
                    if change < b.len() {
                        b[change] ^= 0x80;
                    }
                    assert_eq!(
                        shared_prefix_bytes(&a, &b),
                        a.iter().zip(&b).take_while(|(a, b)| a == b).count()
                    );
                    assert_eq!(
                        shared_suffix_bytes(&a, &b),
                        a.iter()
                            .rev()
                            .zip(b.iter().rev())
                            .take_while(|(a, b)| a == b)
                            .count()
                    );
                }
            }
        }
        let mut seed = 17u64;
        let mut sequence = |length| {
            (0..length)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    code_to_base(seed >> 62)
                })
                .collect::<Vec<_>>()
        };
        let interior = sequence(4096);
        let mut first = sequence(89);
        first.extend_from_slice(&interior);
        first.extend(sequence(47));
        let mut second = sequence(137);
        second.extend_from_slice(&interior);
        second.extend(sequence(23));
        let mut third = second.clone();
        third[700] = code_to_base((u64::from(BASE_CODES[third[700] as usize]) + 1) & 3);
        for k in (3..=31).step_by(2) {
            let transcripts = [first.as_slice(), second.as_slice(), third.as_slice()];
            let expected = transcripts
                .iter()
                .flat_map(|seq| seq.windows(k))
                .map(|window| {
                    let (a, b) = encode_kmer_pair(window).unwrap();
                    a.min(b)
                })
                .collect::<HashSet<_>>();
            let mut anchors = AnchorMap::default();
            let mut sparse = Vec::new();
            let mut every_anchor = AnchorMap::default();
            let mut dense = Vec::new();
            for (id, seq) in transcripts.into_iter().enumerate() {
                insert_interior_unique_kmers::<4>(&mut sparse, seq, k, id, &mut anchors).unwrap();
                insert_interior_unique_kmers::<0>(&mut dense, seq, k, id, &mut every_anchor)
                    .unwrap();
            }
            assert_eq!(
                sparse.iter().copied().collect::<HashSet<_>>(),
                expected,
                "sparse interior k={k}"
            );
            assert_eq!(
                dense.iter().copied().collect::<HashSet<_>>(),
                expected,
                "dense interior k={k}"
            );
            assert!(sparse.len() < first.len() + second.len() + third.len() - 3 * (k - 1));
        }
        let previous = b"ACGTTGCATGTCAGTACCTGATCGATGCTAGCATCGTACTGATGCGTACGTAGCTAGCATGCTAGT";
        for k in (5..=31).step_by(2) {
            for prefix in [0, k - 1, k, k + 1, 40] {
                for suffix in [0, k - 1, k, k + 1, 40] {
                    let mut current = previous[..prefix].to_vec();
                    current.extend_from_slice(b"GTACCGATTCGA");
                    current.extend_from_slice(&previous[previous.len() - suffix..]);
                    let words = |seq: &[u8]| {
                        seq.windows(k)
                            .map(|word| {
                                let (forward, reverse) = encode_kmer_pair(word).unwrap();
                                forward.min(reverse)
                            })
                            .collect::<HashSet<_>>()
                    };
                    let mut expected = words(previous);
                    expected.extend(words(&current));
                    let mut actual = words(previous);
                    actual.extend(words(novel_transcript_region(&current, previous, k)));
                    assert_eq!(actual, expected, "k={k}, prefix={prefix}, suffix={suffix}");
                    let unrelated = vec![b'C'; previous.len()];
                    let history = [
                        previous.to_vec(),
                        unrelated,
                        current.clone(),
                        previous.to_vec(),
                    ]
                    .into_iter()
                    .enumerate()
                    .map(|(id, seq)| Transcript {
                        name: id.to_string(),
                        original_len: seq.len() as u32,
                        seq,
                    })
                    .collect::<Vec<_>>();
                    let mut anchored = HashSet::new();
                    let mut prefixes = AnchorMap::default();
                    let mut suffixes = AnchorMap::default();
                    for i in 0..history.len() {
                        anchored.extend(words(anchored_transcript_region(
                            &history,
                            i,
                            k,
                            &mut prefixes,
                            &mut suffixes,
                        )));
                    }
                    let all = history
                        .iter()
                        .flat_map(|transcript| words(&transcript.seq))
                        .collect::<HashSet<_>>();
                    assert_eq!(anchored, all, "anchor seam k={k}");
                    let mut observed = HashSet::new();
                    let mut complete = HashSet::new();
                    for (i, transcript) in history.iter().enumerate() {
                        observed.extend(words(novel_transcript_history(
                            &transcript.seq,
                            &history[..i],
                            k,
                        )));
                        complete.extend(words(&transcript.seq));
                    }
                    assert_eq!(
                        observed, complete,
                        "history seam k={k},prefix={prefix},suffix={suffix}"
                    );
                }
            }
            assert!(novel_transcript_region(previous, previous, k).is_empty());
        }
        let mut transcripts = (0..3073)
            .map(|id| {
                let mut seq = b"ACGT".repeat(96);
                seq[id % 384] = b"TGCA"[id & 3];
                Transcript {
                    name: id.to_string(),
                    original_len: seq.len() as u32,
                    seq,
                }
            })
            .collect::<Vec<_>>();
        for seq in [Vec::new(), b"ACG".to_vec()] {
            transcripts.push(Transcript {
                name: transcripts.len().to_string(),
                original_len: seq.len() as u32,
                seq,
            });
        }
        let expected = collect_kmers(&transcripts, 31, 1).unwrap().sorted;
        for budget in [64, 65536, usize::MAX] {
            let mut owned = transcripts.clone();
            let mut actual = Vec::new();
            let mut identifiers = Vec::new();
            for group in bounded_transcript_groups(&mut owned, 31, budget) {
                assert!(!group.is_empty() && group.len() <= 2048);
                let before_last = group[..group.len() - 1]
                    .iter()
                    .map(|tr| kmer_position_count(tr.seq.len(), 31))
                    .sum::<usize>();
                assert!(before_last < budget);
                let mut prefixes = AnchorMap::default();
                let mut suffixes = AnchorMap::default();
                let mut anchors = AnchorMap::default();
                for (id, transcript) in group.iter().enumerate() {
                    identifiers.push(transcript.name.clone());
                    let sequence =
                        anchored_transcript_region(group, id, 31, &mut prefixes, &mut suffixes);
                    insert_interior_unique_kmers::<4>(&mut actual, sequence, 31, id, &mut anchors)
                        .unwrap();
                }
            }
            actual.sort_unstable();
            actual.dedup();
            assert_eq!(actual, expected, "bounded groups budget={budget}");
            assert_eq!(
                identifiers,
                transcripts
                    .iter()
                    .map(|tr| tr.name.clone())
                    .collect::<Vec<_>>()
            );
        }

        let mut owned = transcripts.clone();
        let (actual, packed) = collect_kmers_and_pack(&mut owned, 31, 8).unwrap();
        assert_eq!(actual.sorted, expected);
        assert_eq!(packed.len(), transcripts.len());
        for ((packed, source), retired) in packed.iter().zip(&transcripts).zip(&owned) {
            assert_eq!(packed.len(), source.seq.len());
            assert_eq!(retired.name, source.name);
            assert_eq!(retired.original_len, source.original_len);
            assert!(retired.seq.is_empty());
            for (index, &base) in source.seq.iter().enumerate() {
                assert_eq!(packed.base(index), base);
            }
        }
    }

    #[test]
    fn kmer_roundtrip() {
        let seq = b"ACGTACGTAGCCATGACTACGATCGATCGTAGCAGTACGTACGATGCATCGATGCTACGTAGCATGC";
        for k in (3..=31).step_by(2) {
            let mut observed = Vec::new();
            scan_kmer_pairs(seq, k, |pos, fwd, rev| {
                observed.push((pos, fwd, rev));
                Ok(())
            })
            .unwrap();
            let expected = seq
                .windows(k)
                .enumerate()
                .map(|(pos, bases)| {
                    let (fwd, rev) = encode_kmer_pair(bases).unwrap();
                    (pos, fwd, rev)
                })
                .collect::<Vec<_>>();
            assert_eq!(observed, expected);
            for invalid in 0..seq.len() {
                let mut bad = seq.to_vec();
                bad[invalid] = b'N';
                let mut count = 0;
                assert!(
                    scan_kmer_pairs(&bad, k, |_, _, _| {
                        count += 1;
                        Ok(())
                    })
                    .is_err()
                );
                assert_eq!(count, invalid.saturating_sub(k - 1));
            }
        }

        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let word = encode_kmer(seq).expect("encode");
        assert_eq!(decode_kmer_word(word, seq.len()), seq);
        let mut random = 29u64;
        for k in 1..=32 {
            for _ in 0..100 {
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let expected = (0..k)
                    .map(|i| b"ACGT"[((random >> (62 - 2 * i)) & 3) as usize])
                    .collect::<Vec<_>>();
                let mut actual = vec![b'X'; 5];
                append_kmer_bases(&mut actual, random, k);
                assert_eq!(&actual[..5], b"XXXXX");
                assert_eq!(&actual[5..], expected);
            }
        }
    }

    #[test]
    fn fast_reverse_complement_matches_reference() {
        for k in (3..=31).step_by(2) {
            let mask = kmer_mask(k);
            let seq = &b"ACGTTCAGGATCCGTAACCGTTAGCATGCTA"[..k];
            let word = encode_kmer(seq).expect("encode");
            assert_eq!(
                reverse_complement_word(word, k, mask),
                reverse_complement_word_slow(word, k)
            );
        }
    }

    #[test]
    fn bidirectional_edge_states_match_separate_orientation_oracle() {
        for parallel in [false, true] {
            for bias in [0, 1 << 28] {
                let offsets = (0u32..1537)
                    .map(|i| {
                        bias + (i / 256) * 100_000
                            + (i % 256) * [1, 20, 300][(i / 256) as usize % 3]
                    })
                    .collect::<Vec<_>>();
                let ranks = RankDirectory::new(offsets.clone(), parallel);
                assert_eq!(ranks.delta_bits, if bias == 0 { 12 } else { 16 });
                for (i, &expected) in offsets.iter().enumerate() {
                    assert_eq!(ranks.get(i), expected as usize);
                    if i + 1 < offsets.len() {
                        assert_eq!(ranks.range(i), (expected as usize, offsets[i + 1] as usize));
                    }
                }
            }
            let wide = [0, (1 << 31) - 1, 1 << 31, u32::MAX];
            let ranks = RankDirectory::new(wide.to_vec(), parallel);
            for (i, expected) in wide.into_iter().enumerate() {
                assert_eq!(ranks.get(i), expected as usize);
            }
        }

        for k in [3, 5, 7] {
            let mask = kmer_mask(k);
            for value in 0..1u64 << (2 * k) {
                let word = value << (64 - 2 * k);
                let reverse = reverse_complement_word(word, k, mask);
                let self_edges = (0..4)
                    .filter(|&base| {
                        let next = shift_append_word(word, base, k, mask);
                        next == word || next == reverse
                    })
                    .count();
                assert!(
                    self_edges <= 1,
                    "multiple self successors k={k} word={word}"
                );
            }
        }
        for length in 5..=12 {
            let mut entries = (0..length - 1)
                .map(|value| (value as u64) << 2)
                .collect::<Vec<_>>();
            entries.push(1u64 << 36);
            let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
            for query in [0, 1, 4, 8, 44, 48, 1u64 << 36, 1u64 << 37] {
                assert_eq!(
                    lookup.find(&entries, query),
                    entries.binary_search(&query).ok()
                );
                assert_eq!(
                    lookup.find_many(&entries, [query; 4]),
                    [entries.binary_search(&query).ok(); 4]
                );
            }
        }

        let entries = [0, 4, 8, 1u64 << 16];
        let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
        for query in [0, 1, 4, 8, 12, 1u64 << 16, 1u64 << 17] {
            assert_eq!(
                lookup.find(&entries, query),
                entries.binary_search(&query).ok()
            );
        }

        let entries = [0, 4, 1u64 << 21];
        let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
        for query in [0, 1, 4, 8, 1u64 << 21, 1u64 << 22] {
            assert_eq!(
                lookup.find(&entries, query),
                entries.binary_search(&query).ok()
            );
        }

        // Equal guard bits still require complete-key equality, including padding.
        for entries in [vec![0, 1u64 << 32], vec![0, 4]] {
            let lookup = KmerLookup::with_threads(&entries, 1, 31).unwrap();
            for query in [0, 1, 4, 8, 1u64 << 32, 1u64 << 33] {
                assert_eq!(
                    lookup.find(&entries, query),
                    entries.binary_search(&query).ok()
                );
                assert_eq!(
                    lookup.find_many(&entries, [query; 4]),
                    [entries.binary_search(&query).ok(); 4]
                );
            }
        }

        for next in [0, (1 << 28) - 1] {
            let cached = CachedEdges {
                packed: (1 << EDGE_OUT_DEGREE_SHIFT)
                    | (1 << (EDGE_STATE_BITS + EDGE_OUT_DEGREE_SHIFT))
                    | (3 << EDGE_OUT_CODE_SHIFT)
                    | (2 << (EDGE_STATE_BITS + EDGE_OUT_CODE_SHIFT))
                    | (1 << (2 * EDGE_STATE_BITS)),
                next: [next, next],
            };
            let compact = compact_cached_edges(cached);
            let cache = EdgeStateCache {
                packed: EdgeStorage::Compact(vec![compact_pair_word(compact)]),
            };
            assert_eq!(cache.neighbor(0, true, 0), (next as usize, true));
            assert_eq!(cache.neighbor(0, false, 0), (next as usize, false));
        }
        let wide = EdgeStateCache {
            packed: EdgeStorage::Wide(vec![CachedEdges {
                packed: 1 << (2 * EDGE_STATE_BITS + 1),
                next: [u32::MAX, 1 << 28],
            }]),
        };
        assert_eq!(wide.neighbor(0, true, 0), (u32::MAX as usize, false));
        assert_eq!(wide.neighbor(0, false, 0), (1 << 28, true));

        assert_eq!(std::mem::size_of::<RankedEdges>(), 9);
        for count in [(1usize << 28) + 1, (1usize << 31) - 5] {
            for degrees in [[0u32, 4], [1, 1], [2, 3]] {
                for flags in 0..4u32 {
                    let cached = CachedEdges {
                        packed: (degrees[0] << EDGE_OUT_DEGREE_SHIFT)
                            | (degrees[1] << (EDGE_STATE_BITS + EDGE_OUT_DEGREE_SHIFT))
                            | (3 << EDGE_OUT_CODE_SHIFT)
                            | (2 << (EDGE_STATE_BITS + EDGE_OUT_CODE_SHIFT))
                            | (1 << EDGE_HAS_NON_SELF_SHIFT)
                            | (flags << (2 * EDGE_STATE_BITS)),
                        next: [0, count as u32 - 1],
                    };
                    let ranked = RankedEdges::new(cached, count);
                    for (orientation, &degree) in degrees.iter().enumerate() {
                        let edge = ranked.edge(orientation);
                        assert_eq!(RankedEdges::degree(edge, count), degree);
                        assert_eq!(
                            (ranked.codes.load(Ordering::Relaxed) >> (2 * orientation)) & 3,
                            3 - orientation as u8
                        );
                        if degree == 1 {
                            assert_eq!(
                                edge & 0x7fff_ffff,
                                if orientation == 0 {
                                    0
                                } else {
                                    count as u32 - 1
                                }
                            );
                            assert_eq!(edge >> 31, (flags >> orientation) & 1);
                        } else {
                            assert_eq!(edge & 0x7fff_ffff, count as u32 + degree);
                            assert_eq!(edge >> 31, u32::from(orientation == 0));
                        }
                    }
                }
            }
        }

        for k in (3..=31).step_by(2) {
            let mask = kmer_mask(k);
            let mut entries = vec![0];
            let mut state = 42u64;
            for _ in 0..512 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let word = state & mask;
                entries.push(word.min(reverse_complement_word_slow(word, k)));
            }
            let center = entries[1];
            for code in 0..4 {
                for word in [
                    shift_append_word(center, code, k, mask),
                    shift_prepend_word(center, code, mask),
                ] {
                    entries.push(word.min(reverse_complement_word_slow(word, k)));
                }
            }
            entries.sort_unstable();
            entries.dedup();
            let membership = entries.iter().copied().collect::<HashSet<_>>();
            let lookup = KmerLookup::with_threads(&entries, 1, k).unwrap();
            let symmetric_cache = EdgeStateCache::build(&entries, &lookup, k, 1).unwrap();
            let ranked_cache = EdgeStateCache {
                packed: EdgeStorage::Ranked(RankedEdges::build(&entries, &lookup, k, k % 4 == 1)),
            };
            let relative_cache = EdgeStateCache {
                packed: EdgeStorage::Relative(RelativeEdges::build(
                    &entries,
                    &lookup,
                    k,
                    RelativeEdges::prefixes(&entries, &lookup, k).unwrap(),
                    k % 4 == 1,
                )),
            };
            let mut fallback = KmerLookup::new(&entries).unwrap();
            fallback.directory = KmerDirectory::Ranks(RankDirectory::new(
                (0..=(1usize << (64 - fallback.shift)))
                    .map(|bucket| {
                        entries.partition_point(|&word| (word >> fallback.shift) < bucket as u64)
                            as u32
                    })
                    .collect(),
                false,
            ));
            let packed_cache = EdgeStateCache {
                packed: EdgeStorage::Relative7(PackedRelativeEdges::build(
                    &entries,
                    &fallback,
                    k,
                    PackedRelativeEdges::prefixes(&entries, &lookup, k).unwrap(),
                    true,
                )),
            };
            let dense_cache = EdgeStateCache {
                packed: EdgeStorage::Relative6(DenseRelativeEdges::build(
                    &entries,
                    &fallback,
                    k,
                    DenseRelativeEdges::prefixes(&entries, &lookup, k).unwrap(),
                    true,
                )),
            };
            let adaptive_cache = EdgeStateCache {
                packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_plan(
                    &entries,
                    &fallback,
                    k,
                    adaptive_prefixes::PrefixPlan::new(&entries, &fallback, k, 1, 4),
                    true,
                )),
            };
            let mut superset = entries.clone();
            for &word in &entries {
                superset.extend(edge_neighbor_words(word, k, mask));
            }
            superset.sort_unstable();
            superset.dedup();
            let filtered_caches = [&entries, &superset].map(|filter_keys| {
                let filter = negative_filter::NegativeFilter::build(filter_keys, k, 1).unwrap();
                let cache = EdgeStateCache {
                    packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_filtered(
                        &entries,
                        &fallback,
                        k,
                        adaptive_prefixes::PrefixPlan::new(&entries, &fallback, k, 1, 4),
                        &filter,
                        true,
                    )),
                };
                filter.retire();
                cache
            });
            for &word in &entries {
                let queries = [
                    word,
                    word | 1,
                    word ^ (1u64 << (64 - 2 * k)),
                    word ^ (1u64 << 63),
                ];
                let expected = queries.map(|query| entries.binary_search(&query).ok());
                for count in 0..=queries.len() {
                    let mut partial = expected;
                    partial[count..].fill(None);
                    assert_eq!(lookup.find_many_count(&entries, queries, count), partial);
                    assert_eq!(fallback.find_many_count(&entries, queries, count), partial);
                }

                assert_eq!(
                    lookup.find_many(&entries, queries),
                    expected,
                    "inline k={k}"
                );
                assert_eq!(
                    fallback.find_many(&entries, queries),
                    expected,
                    "rank k={k}"
                );
                for (query, expected) in queries.into_iter().zip(expected) {
                    assert_eq!(lookup.find(&entries, query), expected);
                    assert_eq!(fallback.find(&entries, query), expected);
                }
            }
            for &word in &entries {
                let oracle = |oriented| {
                    let mut incoming = 0;
                    let mut outgoing = 0;
                    let mut out_code = 0;
                    let mut non_self = false;
                    let mut self_edges = 0;
                    for code in 0..4 {
                        let next = shift_append_word(oriented, code, k, mask);
                        let next = next.min(reverse_complement_word_slow(next, k));
                        self_edges += usize::from(next == word);
                        if membership.contains(&next) {
                            outgoing += 1;
                            out_code = code as u32;
                            non_self |= next != word;
                        }
                        let prev = shift_prepend_word(oriented, code, mask);
                        let prev = prev.min(reverse_complement_word_slow(prev, k));
                        incoming += u32::from(membership.contains(&prev));
                    }
                    assert!(self_edges <= 1, "canonical self-edge bound k={k}");
                    incoming
                        | (outgoing << EDGE_OUT_DEGREE_SHIFT)
                        | (out_code << EDGE_OUT_CODE_SHIFT)
                        | (u32::from(non_self) << EDGE_HAS_NON_SELF_SHIFT)
                };
                let cached = cache_edge_states(word, &entries, &lookup, k, mask);
                let idx = lookup.find(&entries, word).unwrap();
                let compact =
                    compact_cached_edges(cache_edge_states(word, &entries, &lookup, k, mask));
                let mut compact_entries = vec![0u64; entries.len()];
                compact_entries[idx] = compact_pair_word(compact);
                let compact_cache = EdgeStateCache {
                    packed: EdgeStorage::Compact(compact_entries),
                };
                let mut wide_entries = (0..entries.len())
                    .map(|_| CachedEdges {
                        packed: 0,
                        next: [0; 2],
                    })
                    .collect::<Vec<_>>();
                wide_entries[idx] = cache_edge_states(word, &entries, &lookup, k, mask);
                let wide_cache = EdgeStateCache {
                    packed: EdgeStorage::Wide(wide_entries),
                };
                for forward in [true, false] {
                    let oriented = if forward {
                        word
                    } else {
                        reverse_complement_word_slow(word, k)
                    };
                    assert_eq!(
                        compact_cache.state_bits(idx, forward, word),
                        wide_cache.state_bits(idx, forward, word)
                    );
                    let bits = wide_cache.state_bits(idx, forward, word);
                    assert_eq!(ranked_cache.state_bits(idx, forward, word), bits);
                    assert_eq!(relative_cache.state_bits(idx, forward, word), bits);
                    assert_eq!(packed_cache.state_bits(idx, forward, word), bits);
                    assert_eq!(dense_cache.state_bits(idx, forward, word), bits);
                    assert_eq!(adaptive_cache.state_bits(idx, forward, word), bits);
                    for filtered in &filtered_caches {
                        assert_eq!(filtered.state_bits(idx, forward, word), bits);
                        assert_eq!(
                            filtered.outgoing(idx, forward, oriented),
                            wide_cache.outgoing(idx, forward, oriented)
                        );
                    }

                    assert_eq!(
                        adaptive_cache.outgoing(idx, forward, oriented),
                        wide_cache.outgoing(idx, forward, oriented)
                    );

                    assert_eq!(
                        dense_cache.outgoing(idx, forward, oriented),
                        relative_cache.outgoing(idx, forward, oriented)
                    );
                    assert_eq!(
                        packed_cache.outgoing(idx, forward, oriented),
                        relative_cache.outgoing(idx, forward, oriented)
                    );
                    assert_eq!(relative_cache.is_internal(idx), bits & 0x3f == 9);
                    assert_eq!(
                        relative_cache.incoming_degree(idx, forward),
                        (bits & 7) as u8
                    );
                    assert_eq!(
                        relative_cache.outgoing_continuation(idx, forward, oriented),
                        symmetric_cache.outgoing_continuation(idx, forward, oriented)
                    );
                    assert_eq!(ranked_cache.incoming_degree(idx, forward), (bits & 7) as u8);
                    assert_eq!(ranked_cache.is_internal(idx), bits & 0x3f == 9);
                    assert_eq!(
                        ranked_cache.outgoing_continuation(idx, forward, oriented),
                        symmetric_cache.outgoing_continuation(idx, forward, oriented)
                    );
                    assert_eq!(
                        symmetric_cache.state_bits(idx, forward, word),
                        bits,
                        "symmetric k={k}"
                    );
                    if (bits >> EDGE_OUT_DEGREE_SHIFT) & 7 == 1 {
                        assert_eq!(
                            symmetric_cache.neighbor(idx, forward, oriented),
                            wide_cache.neighbor(idx, forward, oriented)
                        );
                    }

                    assert_eq!(wide_cache.is_internal(idx), bits & 0x3f == 9);
                    assert_eq!(compact_cache.is_internal(idx), bits & 0x3f == 9);

                    assert_eq!(
                        compact_cache.incoming_degree(idx, forward),
                        (bits & 7) as u8
                    );
                    assert_eq!(wide_cache.incoming_degree(idx, forward), (bits & 7) as u8);
                    let expected = if (bits >> EDGE_OUT_DEGREE_SHIFT) & 7 == 1 {
                        let (index, orientation) = wide_cache.neighbor(idx, forward, oriented);
                        Some((
                            index,
                            orientation,
                            u64::from((bits >> EDGE_OUT_CODE_SHIFT) & 3),
                        ))
                    } else {
                        None
                    };
                    assert_eq!(wide_cache.outgoing(idx, forward, oriented), expected);
                    assert_eq!(ranked_cache.outgoing(idx, forward, oriented), expected);
                    assert_eq!(relative_cache.outgoing(idx, forward, oriented), expected);
                    assert_eq!(compact_cache.outgoing(idx, forward, oriented), expected);
                    assert_eq!(symmetric_cache.outgoing(idx, forward, oriented), expected);
                    assert_eq!(
                        symmetric_cache.outgoing_continuation(idx, forward, oriented),
                        expected.filter(|&(target, orientation, _)| symmetric_cache
                            .incoming_degree(target, orientation)
                            == 1)
                    );

                    if ((wide_cache.state_bits(idx, forward, word) >> EDGE_OUT_DEGREE_SHIFT)
                        & EDGE_IN_DEGREE_MASK)
                        == 1
                    {
                        assert_eq!(
                            compact_cache.neighbor(idx, forward, oriented),
                            wide_cache.neighbor(idx, forward, oriented)
                        );
                    }
                }

                for (orientation, oriented) in [word, reverse_complement_word_slow(word, k)]
                    .into_iter()
                    .enumerate()
                {
                    let bits = (cached.packed >> (orientation * EDGE_STATE_BITS))
                        & ((1 << EDGE_STATE_BITS) - 1);
                    if (bits >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK == 1 {
                        let code = (bits >> EDGE_OUT_CODE_SHIFT) & 3;
                        let next = shift_append_word(oriented, u64::from(code), k, mask);
                        let canonical = next.min(reverse_complement_word_slow(next, k));
                        assert_eq!(entries[cached.next[orientation] as usize], canonical);
                        assert_eq!(
                            cached.packed & (1 << (2 * EDGE_STATE_BITS + orientation)) != 0,
                            next == canonical
                        );
                    }
                }
                let expected = oracle(word)
                    | (oracle(reverse_complement_word_slow(word, k)) << EDGE_STATE_BITS);
                assert_eq!(
                    pack_edge_states(word, &entries, &lookup, k, mask),
                    expected,
                    "k={k}, word={word}"
                );
            }
        }
    }

    #[test]
    fn rolling_minimizer_hashes_match_byte_oracle_at_every_window() {
        let mut state = 42u64;
        let random = (0..512)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect::<Vec<_>>();
        for seq in [random, vec![b'A'; 128], b"acgt".repeat(32)] {
            for g in 3..=29 {
                let actual = minimizer_rep_hashes(&seq, g).unwrap();
                let expected = seq
                    .windows(g)
                    .map(|window| minimizer_rep_hash(window).unwrap())
                    .collect::<Vec<_>>();
                assert_eq!(actual, expected, "g={g}");
            }
        }
        assert!(minimizer_rep_hashes(b"ACNTACGT", 3).is_err());
        assert!(minimizer_rep_hashes(b"ACGTACNT", 3).is_err());
    }

    #[test]
    fn block_sweep_matches_scan_builder() {
        let row = |trid| TRInfo {
            trid,
            start: 0,
            stop: 3,
            ..TRInfo::default()
        };
        let counts = [3u32, 1, 5, 1];
        for compact in [false, true] {
            let source = vec![vec![(2usize, row(1)), (0, row(2))], vec![(2, row(3))]];
            let shards = source
                .into_iter()
                .map(|records| {
                    let mut shard = TrInfoShard::new(compact);
                    for (unitig, info) in records {
                        shard.push(unitig, info, counts[unitig]);
                    }
                    shard
                })
                .collect();
            let flat = TranscriptLocations::from_shards(shards, &counts);
            assert_eq!(flat.offsets, [0, 1, 1, 3, 3]);
            assert_eq!(flat.records, [row(2), row(1), row(3)]);
        }
        let mut wide = TrInfoShard::new(false);
        wide.push(u32::MAX as usize, row(7), 3);
        assert_eq!(
            wide.iter(&[]).collect::<Vec<_>>(),
            [(u32::MAX as usize, row(7))]
        );

        let mut grouped = BlockPositions::default();
        grouped.push([]);
        grouped.push([1, 1, 2, 2]);
        grouped.push([2, 3]);
        assert_eq!(
            grouped.iter().collect::<Vec<_>>(),
            [vec![], vec![1, 2], vec![2, 3]]
        );
        grouped.clear();
        assert_eq!(grouped.len(), 0);
        grouped.push([u32::MAX]);
        assert_eq!(grouped.iter().collect::<Vec<_>>(), [vec![u32::MAX]]);
        let kmer_count = 64;
        let mut infos = (0..48)
            .map(|i| TRInfo {
                trid: (i % 7) as u32,
                pos: (i * 11) as u32,
                start: i as u32,
                stop: (i + 5) as u32,
            })
            .collect::<Vec<_>>();
        let mut breakpoints = vec![0, kmer_count];
        for info in &infos {
            breakpoints.push(info.start);
            breakpoints.push(info.stop);
        }
        breakpoints.sort_unstable();
        breakpoints.dedup();

        let mut scan_infos = infos.clone();
        let expected =
            build_blocks_for_unitig_by_scan(&mut scan_infos, &breakpoints, kmer_count, usize::MAX);
        let actual = build_blocks_for_unitig(&mut infos, kmer_count, usize::MAX);
        assert_eq!(actual, expected);
    }

    #[test]
    fn sliding_minimizer_selection_matches_naive_reference() {
        let mut scratch = MinimizerScratch::default();
        let mut rng = 1234567u64;
        let random = (0..128)
            .map(|_| {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                b"ACGT"[(rng & 3) as usize]
            })
            .collect::<Vec<_>>();
        let repeated = b"ACGT".repeat(32);
        let tied = vec![b'C'; 128];
        for seq in [&random, &repeated, &tied] {
            for k in (5..=31).step_by(2) {
                for g in [3, k - 2] {
                    for base in [0, u64::from(u32::MAX) - 64] {
                        let mut keys = Vec::new();
                        let mut blocks = vec![Vec::new(); if base == 0 { 1 } else { 2 }];
                        insert_selected_minimizers(
                            seq,
                            k,
                            g,
                            base,
                            &mut keys,
                            &mut blocks,
                            &mut scratch,
                        )
                        .unwrap();
                        let positions = blocks.iter().enumerate().flat_map(|(block, positions)| {
                            positions
                                .iter()
                                .map(move |&pos| ((block as u64) << 32) | u64::from(pos))
                        });
                        let actual = positions
                            .zip(keys.into_iter().map(u64::to_be_bytes))
                            .collect::<Vec<_>>();
                        let mut expected = selected_minimizers_naive(seq, k, g);
                        expected.sort_unstable();
                        expected.dedup_by_key(|(pos, _)| *pos);
                        let expected = expected
                            .into_iter()
                            .map(|(pos, key)| (base + u64::from(pos), key))
                            .collect::<Vec<_>>();
                        assert_eq!(actual, expected, "k={k},g={g},base={base}");
                    }
                }
            }
        }
        for position in 0..random.len() {
            let mut invalid = random.clone();
            invalid[position] = b'N';
            assert!(
                insert_selected_minimizers(
                    &invalid,
                    31,
                    23,
                    0,
                    &mut Vec::new(),
                    &mut [Vec::new()],
                    &mut scratch
                )
                .is_err(),
                "invalid position={position}"
            );
        }
    }

    fn selected_minimizers_naive(seq: &[u8], k: usize, g: usize) -> Vec<(u32, [u8; 8])> {
        let Some(bounds) = bifrost_minimizer_bounds(k, g) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for kmer_start in 0..=seq.len() - k {
            let mut best_hash = u64::MAX;
            let mut selected = Vec::new();
            for rel_pos in bounds {
                let pos = kmer_start + rel_pos;
                let min_seq = &seq[pos..pos + g];
                let hash = minimizer_rep_hash(min_seq).expect("hash minimizer");
                let bytes = encode_minimizer_rep(min_seq).expect("encode minimizer");
                match hash.cmp(&best_hash) {
                    std::cmp::Ordering::Less => {
                        best_hash = hash;
                        selected.clear();
                        selected.push((pos as u32, bytes));
                    }
                    std::cmp::Ordering::Equal => selected.push((pos as u32, bytes)),
                    std::cmp::Ordering::Greater => {}
                }
            }
            out.extend(selected);
        }
        out
    }

    #[test]
    fn parallel_path_preparation_preserves_serial_claim_order() {
        assert_eq!(std::mem::size_of::<JoinQuery>(), 13);
        for source in [0u64, (1 << 33) - 1] {
            if let Ok(source) = usize::try_from(source) {
                let query = JoinQuery::new(source, u64::MAX).unwrap();
                assert_eq!(query.source(), source);
                assert_eq!(query.word(), u64::MAX);
            }
        }
        if let Ok(overflow) = usize::try_from(1u64 << 40) {
            assert!(JoinQuery::new(overflow, 0).is_err());
        }

        assert!(sparse_locations::worthwhile(1 << 30, 1 << 20, 1 << 20));
        assert!(!sparse_locations::worthwhile(1 << 30, 1 << 29, 1));

        assert_eq!(join_state(&JoinLinks::Wide(vec![usize::MAX]), 0), None);
        assert_eq!(join_state(&JoinLinks::Wide(vec![0]), 0), Some(0));
        let paths = (0..1200)
            .map(|id| UnitigPath {
                ends: [u32::MAX; 2],
                seq: if id & 1 == 0 {
                    Box::from(&b"ACGTTGCATGCA"[..])
                } else {
                    Box::from(&b"TGCATGCAACGT"[..])
                },
            })
            .collect::<Vec<_>>();
        for mode in 0..3 {
            let active = (0..paths.len())
                .map(|id| match mode {
                    0 => true,
                    1 => false,
                    _ => id % 3 != 0,
                })
                .collect::<Vec<_>>();
            let mut serial = build_oriented_starts(&paths, &active, 5, 1)
                .unwrap()
                .into_iter()
                .map(|entry| (entry.word(), entry.unitig_id(), entry.forward()))
                .collect::<Vec<_>>();
            let mut parallel = build_oriented_starts(&paths, &active, 5, 8)
                .unwrap()
                .into_iter()
                .map(|entry| (entry.word(), entry.unitig_id(), entry.forward()))
                .collect::<Vec<_>>();
            serial.sort_unstable();
            parallel.sort_unstable();
            assert_eq!(parallel, serial);
        }
        let mut invalid_paths = paths
            .iter()
            .map(|path| UnitigPath {
                ends: [u32::MAX; 2],
                seq: path.seq.clone(),
            })
            .collect::<Vec<_>>();
        invalid_paths[713].seq[0] = b'N';
        assert!(
            build_oriented_starts(&invalid_paths, &vec![true; invalid_paths.len()], 5, 8).is_err()
        );
        let wide_start = OrientedStart {
            word: 0x9abc_def0_1234_5671,
            unitig: u32::MAX,
        };
        assert_eq!(wide_start.unitig_id(), u32::MAX as usize);
        assert!(wide_start.forward());
        assert_eq!(wide_start.word(), 0x9abc_def0_1234_5670);

        assert_eq!(join_state(&JoinLinks::Narrow(vec![u32::MAX]), 0), None);
        assert_eq!(join_state(&JoinLinks::Narrow(vec![0]), 0), Some(0));
        assert_eq!(
            join_state(&JoinLinks::Narrow(vec![u32::MAX - 1]), 0),
            Some(u32::MAX as usize - 1)
        );
        assert_eq!(
            join_state(&JoinLinks::Wide(vec![usize::MAX - 1]), 0),
            Some(usize::MAX - 1)
        );

        for index in [0, 1, (1 << 30) - 1] {
            for base in 0..4 {
                let narrow = <u32 as PathStep>::pack(index, base);
                let wide = <u64 as PathStep>::pack(index, base);
                assert_eq!(narrow.index(), index);
                assert_eq!(wide.index(), index);
                assert_eq!(narrow.base(), base);
                assert_eq!(wide.base(), base);
            }
        }
        let wide = <u64 as PathStep>::pack(u32::MAX as usize, 3);
        assert_eq!(wide.index(), u32::MAX as usize);

        let mut cycle_rng = 987654321u64;
        let mut ring = (0..10000)
            .map(|_| {
                cycle_rng ^= cycle_rng << 13;
                cycle_rng ^= cycle_rng >> 7;
                cycle_rng ^= cycle_rng << 17;
                b"ACGT"[(cycle_rng & 3) as usize]
            })
            .collect::<Vec<_>>();
        ring.extend_from_within(..30);
        let cycle = vec![Transcript {
            name: "bounded-cycle".into(),
            original_len: ring.len() as u32,
            seq: ring,
        }];
        let words = collect_kmers(&cycle, 31, 1).unwrap();
        let lookup = KmerLookup::new(&words.sorted).unwrap();
        let edges = EdgeStateCache::build(&words.sorted, &lookup, 31, 1).unwrap();
        let topology = GraphTopology::new(&words.sorted, &lookup, &edges, 31);
        let start = choose_any_start_orientation(words.sorted[0], 0, &topology);
        assert!(
            prepare_path_group::<u64, true>(
                &[(0, start)],
                &mut [Vec::new()],
                &topology,
                &VisitedKmers::new(words.sorted.len())
            )[0]
            .is_none()
        );
        let serial = compact_unitigs(&words.sorted, &topology, 1)
            .unwrap()
            .into_iter()
            .map(|path| path.seq)
            .collect::<Vec<_>>();
        let parallel = compact_unitigs_parallel(&words.sorted, &topology, 8)
            .unwrap()
            .into_iter()
            .map(|path| path.seq)
            .collect::<Vec<_>>();
        assert_eq!(parallel, serial, "bounded cycle fallback");
        let wide = compact_unitigs_parallel_steps::<u64>(&words.sorted, &topology, 8)
            .unwrap()
            .into_iter()
            .map(|path| path.seq)
            .collect::<Vec<_>>();
        assert_eq!(wide, serial, "wide bounded cycle fallback");
        for batch in [512, 4096] {
            let mut visited = VisitedKmers::new(words.sorted.len());
            let mut paths = PathChunks::default();
            let workers = IndexWorkers::new(8).unwrap();
            finish_unitig_cycles_batched::<u64>(
                &words.sorted,
                &topology,
                &mut visited,
                &mut paths,
                &workers,
                batch,
            )
            .unwrap();
            let actual = paths
                .finish()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            assert_eq!(actual, serial, "ordered cycle batch {batch}");
        }

        let relative_edges = EdgeStateCache {
            packed: EdgeStorage::Relative(RelativeEdges::build(
                &words.sorted,
                &lookup,
                31,
                RelativeEdges::prefixes(&words.sorted, &lookup, 31).unwrap(),
                true,
            )),
        };
        let relative_topology = GraphTopology::new(&words.sorted, &lookup, &relative_edges, 31);
        let relative = compact_unitigs_parallel_steps::<u64>(&words.sorted, &relative_topology, 8)
            .unwrap()
            .into_iter()
            .map(|path| path.seq)
            .collect::<Vec<_>>();
        assert_eq!(relative, serial, "relative bounded cycle fallback");

        let mut rng = 42u64;
        let mut transcripts = Vec::new();
        for id in 0..1200 {
            let seq = (0..96)
                .map(|_| {
                    rng ^= rng << 13;
                    rng ^= rng >> 7;
                    rng ^= rng << 17;
                    b"ACGT"[(rng & 3) as usize]
                })
                .collect::<Vec<_>>();
            transcripts.push(Transcript {
                name: id.to_string(),
                original_len: seq.len() as u32,
                seq,
            });
        }
        let first = transcripts[0].seq.clone();
        for (id, seq) in [
            first[7..83].to_vec(),
            reverse_complement_seq(&first[11..89]),
        ]
        .into_iter()
        .enumerate()
        {
            transcripts.push(Transcript {
                name: format!("contained-{id}"),
                original_len: seq.len() as u32,
                seq,
            });
        }

        transcripts.push(Transcript {
            name: "reverse".into(),
            original_len: first.len() as u32,
            seq: reverse_complement_seq(&first),
        });
        transcripts.push(Transcript {
            name: "cycle".into(),
            original_len: 256,
            seq: b"ACGT".repeat(64),
        });
        for k in [5, 15, 31] {
            let words = collect_kmers(&transcripts, k, 1).unwrap();
            let lookup = KmerLookup::new(&words.sorted).unwrap();
            let edges = EdgeStateCache::build(&words.sorted, &lookup, k, 1).unwrap();
            let topology = GraphTopology::new(&words.sorted, &lookup, &edges, k);
            let expected = compact_unitigs(&words.sorted, &topology, 1)
                .unwrap()
                .into_iter()
                .map(|path| path.seq)
                .collect::<Vec<_>>();
            let relative_edges = EdgeStateCache {
                packed: EdgeStorage::Relative(RelativeEdges::build(
                    &words.sorted,
                    &lookup,
                    k,
                    RelativeEdges::prefixes(&words.sorted, &lookup, k).unwrap(),
                    true,
                )),
            };
            let packed_edges = EdgeStateCache {
                packed: EdgeStorage::Relative7(PackedRelativeEdges::build(
                    &words.sorted,
                    &lookup,
                    k,
                    PackedRelativeEdges::prefixes(&words.sorted, &lookup, k).unwrap(),
                    true,
                )),
            };
            let dense_edges = EdgeStateCache {
                packed: EdgeStorage::Relative6(DenseRelativeEdges::build(
                    &words.sorted,
                    &lookup,
                    k,
                    DenseRelativeEdges::prefixes(&words.sorted, &lookup, k).unwrap(),
                    true,
                )),
            };
            let relative_topology = GraphTopology::new(&words.sorted, &lookup, &relative_edges, k);
            let packed_topology = GraphTopology::new(&words.sorted, &lookup, &packed_edges, k);
            let dense_topology = GraphTopology::new(&words.sorted, &lookup, &dense_edges, k);
            let adaptive_edges = EdgeStateCache {
                packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_plan(
                    &words.sorted,
                    &lookup,
                    k,
                    adaptive_prefixes::PrefixPlan::new(&words.sorted, &lookup, k, 1, 15),
                    true,
                )),
            };
            let adaptive_topology = GraphTopology::new(&words.sorted, &lookup, &adaptive_edges, k);

            if k >= 15 {
                // Force the large exact key geometry on this serial-order oracle.
                // Sparse and empty prefix ranges must reconstruct full words too.
                let mut rank_lookup = KmerLookup::with_threads(&words.sorted, 1, k).unwrap();
                let bits = (2 * k).saturating_sub(40).max(12);
                rank_lookup.shift = 64 - bits;
                let mut cursor = 0;
                let offsets = (0..=(1usize << bits))
                    .map(|bucket| {
                        while cursor < words.sorted.len()
                            && words.sorted[cursor] >> rank_lookup.shift < bucket as u64
                        {
                            cursor += 1;
                        }
                        cursor as u32
                    })
                    .collect();
                rank_lookup.directory = KmerDirectory::Ranks(RankDirectory::new(offsets, false));
                let mut keys =
                    PrefixKeys::for_topology(words.sorted.clone(), &rank_lookup, 1).unwrap();
                let filter =
                    negative_filter::NegativeFilter::build_mapped(&keys, &rank_lookup, k, 8)
                        .unwrap();
                let packed_edges = EdgeStateCache {
                    packed: EdgeStorage::Adaptive(AdaptiveRelativeEdges::build_filtered(
                        &keys,
                        &rank_lookup,
                        k,
                        adaptive_prefixes::PrefixPlan::new(&words.sorted, &rank_lookup, k, 1, 128),
                        &filter,
                        true,
                    )),
                };
                filter.retire();
                let packed_topology = GraphTopology::new(&keys, &rank_lookup, &packed_edges, k);
                for threads in [1, 8] {
                    let paths = if threads == 1 {
                        compact_unitigs(&keys, &packed_topology, threads)
                    } else {
                        compact_unitigs_parallel(&keys, &packed_topology, threads)
                    }
                    .unwrap();
                    assert_eq!(
                        paths.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                        expected.iter().collect::<Vec<_>>(),
                        "packed-key paths k={k}, threads={threads}"
                    );
                    let joined = join_unitig_paths(paths, &packed_topology, threads).unwrap();
                    let raw_paths = compact_unitigs(&words.sorted, &topology, 1).unwrap();
                    let raw_joined = join_unitig_paths(raw_paths, &topology, 1).unwrap();
                    assert_eq!(
                        joined.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                        raw_joined.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                        "packed-key joins k={k}, threads={threads}"
                    );
                }
                let raw_paths = compact_unitigs(&words.sorted, &topology, 1).unwrap();
                let normalize = |groups: Vec<Vec<JoinQuery>>| {
                    groups
                        .into_iter()
                        .flatten()
                        .map(|query| (query.source(), query.word()))
                        .collect::<Vec<_>>()
                };
                let expected_probes =
                    normalize(collect_join_queries(&raw_paths, &packed_topology, 1).unwrap());
                keys.retire_word_hints();
                let retired_topology = GraphTopology::new(&keys, &rank_lookup, &packed_edges, k);
                let actual_probes =
                    normalize(collect_join_queries(&raw_paths, &retired_topology, 1).unwrap());
                assert_eq!(
                    actual_probes, expected_probes,
                    "joins after hint retirement"
                );
            }

            for threads in [2, 8] {
                let actual = compact_unitigs_parallel(&words.sorted, &topology, threads)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                assert_eq!(actual, expected, "k={k}, threads={threads}");
                let relative_paths = compact_unitigs_parallel_steps::<u64>(
                    &words.sorted,
                    &relative_topology,
                    threads,
                )
                .unwrap();
                assert_eq!(
                    relative_paths
                        .iter()
                        .map(|path| &path.seq)
                        .collect::<Vec<_>>(),
                    expected.iter().collect::<Vec<_>>(),
                    "relative paths k={k}, threads={threads}"
                );
                let packed_paths =
                    compact_unitigs_parallel_steps::<u64>(&words.sorted, &packed_topology, threads)
                        .unwrap();
                assert_eq!(
                    packed_paths
                        .iter()
                        .map(|path| &path.seq)
                        .collect::<Vec<_>>(),
                    expected.iter().collect::<Vec<_>>(),
                    "packed paths k={k}, threads={threads}"
                );
                let dense_paths =
                    compact_unitigs_parallel_steps::<u64>(&words.sorted, &dense_topology, threads)
                        .unwrap();
                assert_eq!(
                    dense_paths.iter().map(|path| &path.seq).collect::<Vec<_>>(),
                    expected.iter().collect::<Vec<_>>(),
                    "dense paths k={k}, threads={threads}"
                );
                let adaptive_paths = compact_unitigs_parallel_steps::<u64>(
                    &words.sorted,
                    &adaptive_topology,
                    threads,
                )
                .unwrap();
                assert_eq!(
                    adaptive_paths
                        .iter()
                        .map(|path| &path.seq)
                        .collect::<Vec<_>>(),
                    expected.iter().collect::<Vec<_>>(),
                    "adaptive paths k={k}, threads={threads}"
                );
                let paths = expected
                    .iter()
                    .map(|seq| UnitigPath {
                        seq: seq.clone(),
                        ends: [u32::MAX; 2],
                    })
                    .collect::<Vec<_>>();
                let serial_join = join_unitig_paths(paths.clone(), &topology, 1)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                let adaptive_join = join_unitig_paths(adaptive_paths, &adaptive_topology, threads)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                assert_eq!(
                    adaptive_join, serial_join,
                    "adaptive join k={k}, threads={threads}"
                );
                let relative_join = join_unitig_paths(relative_paths, &relative_topology, threads)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                assert_eq!(
                    relative_join, serial_join,
                    "relative join k={k}, threads={threads}"
                );
                let packed_join = join_unitig_paths(packed_paths, &packed_topology, threads)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                assert_eq!(
                    packed_join, serial_join,
                    "packed join k={k}, threads={threads}"
                );
                let dense_join = join_unitig_paths(dense_paths, &dense_topology, threads)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                assert_eq!(
                    dense_join, serial_join,
                    "dense join k={k}, threads={threads}"
                );
                let parallel_join = join_unitig_paths(paths, &topology, threads)
                    .unwrap()
                    .into_iter()
                    .map(|path| path.seq)
                    .collect::<Vec<_>>();
                assert_eq!(parallel_join, serial_join, "join k={k}, threads={threads}");
                let joined = serial_join
                    .into_iter()
                    .map(|seq| UnitigPath {
                        seq,
                        ends: [u32::MAX; 2],
                    })
                    .collect::<Vec<_>>();
                // Exercise cross-unitig batching directly even below the public
                // parallel-size threshold, including group and tail seams.
                if threads == 2 {
                    let mut expected_locations = Vec::new();
                    for (unitig, path) in joined.iter().enumerate() {
                        store_unitig_locations(
                            path,
                            unitig,
                            k,
                            &words.sorted,
                            &lookup,
                            |_| {},
                            |rank, loc| expected_locations.push((rank, loc.0)),
                        )
                        .unwrap();
                    }
                    for group_size in [1, 63, 64, 65] {
                        let mut actual_locations = Vec::new();
                        for (group, paths) in joined.chunks(group_size).enumerate() {
                            store_location_group::<64, _, _, _>(
                                paths,
                                group * group_size,
                                k,
                                &words.sorted,
                                &lookup,
                                |_| {},
                                |rank, loc| actual_locations.push((rank, loc.0)),
                            )
                            .unwrap();
                        }
                        assert_eq!(
                            actual_locations, expected_locations,
                            "location groups k={k}, size={group_size}"
                        );
                    }
                }
                let locs = build_kmer_locs(&joined, &words.sorted, &lookup, k, threads).unwrap();
                let mut inputs = transcripts.clone();
                let mut packed = pack_transcript_sequences(&mut inputs, threads).unwrap();
                let sparse =
                    sparse_locations::SparseLocations::build(&joined, &packed, k, threads).unwrap();
                assert!(
                    sparse.find(u64::MAX).is_none(),
                    "empty slot is not membership"
                );
                let dense = DenseReferenceLocations {
                    keys: &words.sorted,
                    lookup: &lookup,
                    locations: &locs,
                };
                let dense_shards = build_trinfos(&mut packed, k, &dense, &joined, threads).unwrap();
                let mut inputs = transcripts.clone();
                let mut packed = pack_transcript_sequences(&mut inputs, threads).unwrap();
                let sparse_shards =
                    build_trinfos(&mut packed, k, &sparse, &joined, threads).unwrap();
                let counts = joined
                    .iter()
                    .map(|path| (path.seq.len() - k + 1) as u32)
                    .collect::<Vec<_>>();
                let dense_rows = dense_shards
                    .iter()
                    .flat_map(|shard| shard.iter(&counts))
                    .collect::<Vec<_>>();
                let sparse_rows = sparse_shards
                    .iter()
                    .flat_map(|shard| shard.iter(&counts))
                    .collect::<Vec<_>>();
                assert_eq!(
                    sparse_rows, dense_rows,
                    "sparse exact runs k={k}, threads={threads}"
                );
                for path in &joined {
                    let word = encode_kmer(&path.seq[..k]).unwrap();
                    assert!(
                        sparse.find(word | 1).is_none(),
                        "complete-key padding check"
                    );
                }
            }
        }
    }

    #[test]
    fn compacts_linear_transcript() {
        let tr = Transcript {
            name: "tx".into(),
            original_len: 36,
            seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        };
        let (graph, _) =
            build_kmer_unitig_graph_with_report(vec![tr], 31, 23, -1, 1).expect("graph");
        assert!(graph.unitigs.len() < graph.unitigs[0].len());
        assert_eq!(graph.km_unitigs.len(), 0);
        assert_eq!(graph.nodes.len(), graph.unitigs.len());
    }

    #[test]
    fn joins_cyclic_repeat_unitigs() {
        let seq = b"ACGT".repeat(20);
        let tr = Transcript {
            name: "tx".into(),
            original_len: seq.len() as u32,
            seq,
        };
        let (graph, _) =
            build_kmer_unitig_graph_with_report(vec![tr], 31, 23, -1, 1).expect("graph");
        assert_eq!(graph.unitigs.len(), 1);
        assert_eq!(graph.km_unitigs.len(), 0);
        assert_eq!(graph.nodes.len(), graph.unitigs.len());
    }
}
