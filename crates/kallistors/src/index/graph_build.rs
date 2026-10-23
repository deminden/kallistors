use compaction::UnitigPath;
#[cfg(test)]
use edges::RankedEdges;
use edges::{CachedEdges, EdgeStateCache, EdgeStorage, GraphTopology};
#[cfg(test)]
use joins::{JoinLinks, OrientedStart, OrientedStartLookup};
use locations::DenseReferenceLocations;
#[cfg(test)]
use lookup::RankDirectory;
use lookup::{KmerDirectory, KmerLoc, KmerLookup, VisitedKmers};
mod collection;
mod compaction;
mod edges;
mod joins;
mod locations;
mod lookup;
mod minimizers;
mod payloads;

#[cfg(test)]
use collection::{
    anchored_transcript_region, bounded_transcript_groups, collect_kmers,
    insert_interior_unique_kmers, novel_transcript_history, novel_transcript_region,
    shared_prefix_bytes, shared_suffix_bytes,
};
use collection::{collect_kmers_and_pack, kmer_position_count};
#[cfg(test)]
use compaction::{
    PathChunks, PathStep, choose_any_start_orientation, compact_unitigs_parallel,
    compact_unitigs_parallel_steps, finish_unitig_cycles_batched, prepare_path_group,
};
use compaction::{compact_unitigs, concatenate_shards};
#[cfg(test)]
use edges::{cache_edge_states, compact_cached_edges, compact_pair_word, pack_edge_states};
use edges::{cached_edges_from_found, edge_neighbor_words, reverse_complement_word};
#[cfg(test)]
use joins::find_unique_oriented_start;
#[cfg(test)]
use joins::{JoinQuery, build_oriented_starts, join_state, join_unitig_paths};
use joins::{collect_join_queries, join_unitig_paths_with_queries};
use locations::{build_kmer_locs, pack_transcript_sequences};
#[cfg(test)]
use locations::{store_location_group, store_unitig_locations};
use lookup::{prefetch_kmer, prefetch_rank_delta};
use minimizers::build_unitig_minimizer_index;
#[cfg(test)]
use minimizers::{
    MinimizerScratch, bifrost_minimizer_bounds, insert_selected_minimizers, minimizer_rep_hash,
    minimizer_rep_hashes,
};
#[cfg(test)]
use payloads::{build_blocks_for_unitig, build_blocks_for_unitig_by_scan};
use payloads::{build_node_payloads, build_trinfos};

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

trait ReferenceKmerLocations: Sync {
    fn find(&self, word: u64) -> Option<KmerLoc>;
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

#[cfg(test)]
#[path = "graph_build/tests.rs"]
mod tests;
