use core::range::RangeInclusive;
use std::collections::{BTreeMap, HashSet, VecDeque};
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use super::bifrost::{encode_minimizer_rep, wyhash};
use super::builder::Transcript;
use crate::{Error, Result};

const INVALID_BASE_CODE: u8 = 4;
const PARALLEL_KMER_COLLECTION_MIN: usize = 1_000_000;
const PARALLEL_TRINFO_COLLECTION_MIN: usize = 1_000_000;
const SWEEP_BLOCK_MIN_BREAKPOINTS: usize = 32;
const MIN_KMER_LOOKUP_PREFIX_BITS: usize = 12;
const MAX_KMER_LOOKUP_PREFIX_BITS: usize = 27;
const TARGET_KMERS_PER_LOOKUP_BUCKET: usize = 2;
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
const REP_HASH_VALS: [u64; 4] = [
    2053695854357871005,
    5073395517033431291,
    10060236952204337488,
    7783083932390163561,
];

type MinimizerIndexParts = (Vec<[u8; 8]>, Vec<Vec<u32>>);

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct GraphBuildReport {
    pub(super) graph_build: Duration,
    pub(super) minimizer_index: Duration,
    pub(super) ec_build: Duration,
}

#[derive(Clone)]
struct UnitigPath {
    seq: Vec<u8>,
}

#[derive(Clone, Copy, Default)]
struct KmerLoc(u64);

#[derive(Clone, Copy)]
struct TRInfo {
    trid: u32,
    pos: u32,
    start: u32,
    stop: u32,
}

struct Run {
    unitig_id: usize,
    tr_start_pos: u32,
    start: usize,
    stop: usize,
    last_offset: usize,
    same_strand: bool,
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
    canonical: u64,
    oriented: u64,
}

#[derive(Clone, Copy)]
struct OrientedStart {
    word: u64,
    unitig_id: usize,
    forward: bool,
}

#[derive(Clone, Copy)]
struct JoinCandidate {
    source_id: usize,
    source_forward: bool,
    target_id: usize,
    target_forward: bool,
}

#[derive(Clone, Copy, Default)]
struct EdgeState {
    in_degree: u8,
    out_degree: u8,
    out_neighbor: Option<OrientedNeighbor>,
    has_non_self_out: bool,
}

struct GraphTopology<'a> {
    entries: &'a [u64],
    lookup: &'a KmerLookup,
    edge_states: &'a EdgeStateCache,
    k: usize,
    mask: u64,
    last_shift: usize,
}

struct KmerLookup {
    offsets: Vec<u32>,
    shift: usize,
}

struct EdgeStateCache {
    packed: Vec<u32>,
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
    pub(super) minimizer_keys: Vec<[u8; 8]>,
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
    pub(super) positions: Vec<Vec<u32>>,
}

pub(super) fn build_kmer_unitig_graph_with_report(
    transcripts: &[Transcript],
    k: usize,
    g: usize,
    ec_max_size: i32,
    threads: usize,
) -> Result<(BuildGraph, GraphBuildReport)> {
    let graph_start = Instant::now();
    // Build topology first; EC payloads are keyed to the final compacted unitig coordinates.
    let kmers = collect_kmers(transcripts, k, threads)?;
    let lookup = KmerLookup::new(&kmers.sorted)?;
    let paths = {
        let edge_states = EdgeStateCache::build(&kmers.sorted, &lookup, k, threads)?;
        let topology = GraphTopology::new(&kmers.sorted, &lookup, &edge_states, k);
        let paths = compact_unitigs(&kmers.sorted, &lookup, &topology)?;
        join_unitig_paths(paths, &topology)?
    };
    let kmer_locs = build_kmer_locs(&paths, &kmers.sorted, &lookup, k, threads)?;
    let graph_build = graph_start.elapsed();

    let ec_start = Instant::now();
    let trinfos = build_trinfos(
        transcripts,
        k,
        &kmers.sorted,
        &lookup,
        &kmer_locs,
        paths.len(),
        threads,
    )?;
    let ec_threshold = ec_threshold(ec_max_size);
    let nodes = build_node_payloads(&paths, trinfos, ec_threshold, k);
    let ec_build = ec_start.elapsed();

    let unitigs: Vec<Vec<u8>> = paths.into_iter().map(|path| path.seq).collect();
    let minimizer_start = Instant::now();
    let (minimizer_keys, unitig_bitmap_blocks) = build_unitig_minimizer_index(&unitigs, k, g)?;
    let minimizer_index = minimizer_start.elapsed();

    let transcript_names = transcripts.iter().map(|tr| tr.name.clone()).collect();
    let transcript_lengths = transcripts.iter().map(|tr| tr.original_len).collect();
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
        },
    ))
}

fn collect_kmers(transcripts: &[Transcript], k: usize, threads: usize) -> Result<KmerStore> {
    let counts = transcripts
        .iter()
        .map(|tr| kmer_position_count(tr.seq.len(), k))
        .collect::<Vec<_>>();
    let total_positions = counts.iter().sum::<usize>();
    if threads > 1 && transcripts.len() > 1 && total_positions >= PARALLEL_KMER_COLLECTION_MIN {
        let mut offsets = Vec::with_capacity(transcripts.len() + 1);
        offsets.push(0usize);
        for &count in &counts {
            offsets.push(offsets.last().copied().unwrap_or(0) + count);
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| {
                Error::InvalidFormat(format!("failed to create index worker pool: {err}"))
            })?;
        let mut entries = Vec::with_capacity(total_positions);
        entries.resize_with(total_positions, MaybeUninit::uninit);
        let entries_ptr = entries.as_mut_ptr() as usize;
        pool.install(|| {
            transcripts
                .par_iter()
                .enumerate()
                .try_for_each(|(idx, transcript)| {
                    let start = offsets[idx];
                    let len = offsets[idx + 1] - start;
                    // The prefix offsets partition the preallocated k-mer buffer, so each
                    // transcript writes a distinct range.
                    let out = unsafe {
                        std::slice::from_raw_parts_mut(
                            (entries_ptr as *mut MaybeUninit<u64>).add(start),
                            len,
                        )
                    };
                    fill_transcript_kmers_uninit(out, transcript, k)
                })
        })?;
        let sorted = pool.install(|| {
            let mut entries = assume_init_u64_vec(entries);
            entries.par_sort_unstable();
            entries.dedup();
            entries.shrink_to_fit();
            Ok::<Vec<u64>, Error>(entries)
        })?;
        return Ok(KmerStore { sorted });
    }

    let mut sorted = Vec::with_capacity(total_positions);
    for transcript in transcripts {
        insert_transcript_kmers(&mut sorted, transcript, k)?;
    }

    sorted.sort_unstable();
    sorted.dedup();
    sorted.shrink_to_fit();
    Ok(KmerStore { sorted })
}

fn insert_transcript_kmers(
    entries: &mut Vec<u64>,
    transcript: &Transcript,
    k: usize,
) -> Result<()> {
    if transcript.seq.len() < k {
        return Ok(());
    }
    scan_kmer_pairs(&transcript.seq, k, |_, fwd, rev| {
        entries.push(fwd.min(rev));
        Ok(())
    })
}

fn fill_transcript_kmers_uninit(
    out: &mut [MaybeUninit<u64>],
    transcript: &Transcript,
    k: usize,
) -> Result<()> {
    if transcript.seq.len() < k {
        return Ok(());
    }
    let mut idx = 0usize;
    scan_kmer_pairs(&transcript.seq, k, |_, fwd, rev| {
        out[idx].write(fwd.min(rev));
        idx += 1;
        Ok(())
    })?;
    debug_assert_eq!(idx, out.len());
    Ok(())
}

fn assume_init_u64_vec(mut entries: Vec<MaybeUninit<u64>>) -> Vec<u64> {
    let ptr = entries.as_mut_ptr() as *mut u64;
    let len = entries.len();
    let cap = entries.capacity();
    std::mem::forget(entries);
    unsafe { Vec::from_raw_parts(ptr, len, cap) }
}

fn kmer_position_count(seq_len: usize, k: usize) -> usize {
    if seq_len >= k { seq_len - k + 1 } else { 0 }
}

fn compact_unitigs(
    entries: &[u64],
    lookup: &KmerLookup,
    topology: &GraphTopology<'_>,
) -> Result<Vec<UnitigPath>> {
    let mut visited = VisitedKmers::new(entries.len());
    let mut out = Vec::new();

    // Branch starts produce maximal paths; the second pass catches isolated cycles.
    for &word in entries {
        if visited.contains_word(entries, lookup, word)? {
            continue;
        }
        let Some(start) = choose_branch_start_orientation(word, topology) else {
            continue;
        };
        let path = walk_unitig(start, topology, entries, lookup, &mut visited)?;
        if path.seq.len() >= topology.k {
            out.push(path);
        }
    }

    for &word in entries {
        if visited.contains_word(entries, lookup, word)? {
            continue;
        }
        let start = choose_any_start_orientation(word, topology);
        let path = walk_unitig(start, topology, entries, lookup, &mut visited)?;
        if path.seq.len() >= topology.k {
            out.push(path);
        }
    }

    Ok(out)
}

fn choose_branch_start_orientation(word: u64, topology: &GraphTopology<'_>) -> Option<u64> {
    let reverse = topology.reverse_complement(word);
    let forward_state = topology.degree_state(word);
    let reverse_state = topology.degree_state(reverse);

    let forward_is_start = forward_state.in_degree != 1 && forward_state.out_degree > 0;
    let reverse_is_start = reverse_state.in_degree != 1 && reverse_state.out_degree > 0;
    match (forward_is_start, reverse_is_start) {
        (true, true) if has_non_self_out(forward_state) => Some(word),
        (true, true) if has_non_self_out(reverse_state) => Some(reverse),
        (true, _) => Some(word),
        (_, true) => Some(reverse),
        _ => None,
    }
}

fn choose_any_start_orientation(word: u64, topology: &GraphTopology<'_>) -> u64 {
    let reverse = topology.reverse_complement(word);
    let forward_state = topology.degree_state(word);
    let reverse_state = topology.degree_state(reverse);

    if has_non_self_out(forward_state) {
        word
    } else if has_non_self_out(reverse_state) {
        reverse
    } else if forward_state.out_degree > 0 || reverse_state.out_degree == 0 {
        word
    } else {
        reverse
    }
}

fn has_non_self_out(state: EdgeState) -> bool {
    state.has_non_self_out
}

fn walk_unitig(
    start: u64,
    topology: &GraphTopology<'_>,
    entries: &[u64],
    lookup: &KmerLookup,
    visited: &mut VisitedKmers,
) -> Result<UnitigPath> {
    let mut seq = decode_kmer_word(start, topology.k);
    let mut cur = start;

    loop {
        let word = topology.canonical(cur);
        let Some(_word_idx) = visited.insert_word(entries, lookup, word)? else {
            break;
        };

        let state = topology.state(cur);
        if state.out_degree != 1 {
            break;
        }
        let next = state.out_neighbor.expect("single outgoing neighbor");
        if visited.contains_word(entries, lookup, next.canonical)?
            || topology.in_degree(next.oriented) != 1
        {
            break;
        }
        seq.push(code_to_base(topology.last_base_code(next.oriented)));
        cur = next.oriented;
    }

    Ok(UnitigPath { seq })
}

fn join_unitig_paths(
    paths: Vec<UnitigPath>,
    topology: &GraphTopology<'_>,
) -> Result<Vec<UnitigPath>> {
    let active = vec![true; paths.len()];
    // Stitch only unique k-1 overlaps left behind after the primary unitig walk.
    let starts = build_oriented_starts(&paths, &active, topology.k)?;
    let candidates = build_join_candidates(&paths, &active, &starts, topology)?;
    let (next, prev) = build_join_edges(paths.len(), candidates);
    let mut used = vec![false; paths.len()];
    let mut out = Vec::new();

    for unitig_id in 0..paths.len() {
        if used[unitig_id] {
            continue;
        }
        let start = rewind_join_start(state_id(unitig_id, true), &prev);
        let seq = joined_sequence_from_state(start, &paths, &next, &mut used, topology.k);
        out.push(UnitigPath { seq });
    }

    Ok(out)
}

fn build_oriented_starts(
    paths: &[UnitigPath],
    active: &[bool],
    k: usize,
) -> Result<Vec<OrientedStart>> {
    let mut starts = Vec::with_capacity(active.iter().filter(|&&is_active| is_active).count() * 2);
    let mask = kmer_mask(k);
    for (unitig_id, path) in paths.iter().enumerate() {
        if !active[unitig_id] {
            continue;
        }
        let first = encode_kmer(&path.seq[..k])?;
        let last = encode_kmer(&path.seq[path.seq.len() - k..])?;
        starts.push(OrientedStart {
            word: first,
            unitig_id,
            forward: true,
        });
        starts.push(OrientedStart {
            word: reverse_complement_word(last, k, mask),
            unitig_id,
            forward: false,
        });
    }
    starts.sort_unstable_by_key(|start| start.word);
    Ok(starts)
}

fn build_join_candidates(
    paths: &[UnitigPath],
    active: &[bool],
    starts: &[OrientedStart],
    topology: &GraphTopology<'_>,
) -> Result<Vec<JoinCandidate>> {
    let mut out = Vec::new();
    for (source_id, path) in paths.iter().enumerate() {
        if !active[source_id] {
            continue;
        }
        for source_forward in [true, false] {
            let end = oriented_end_word(&path.seq, source_forward, topology.k, topology.mask)?;
            let Some(target) = join_target(end, starts, topology) else {
                continue;
            };
            if target.unitig_id != source_id {
                out.push(JoinCandidate {
                    source_id,
                    source_forward,
                    target_id: target.unitig_id,
                    target_forward: target.forward,
                });
            }
        }
    }
    out.sort_unstable_by_key(|candidate| {
        (
            candidate.source_id,
            !candidate.source_forward,
            candidate.target_id,
            !candidate.target_forward,
        )
    });
    Ok(out)
}

fn build_join_edges(
    unitig_count: usize,
    candidates: Vec<JoinCandidate>,
) -> (Vec<Option<usize>>, Vec<Option<usize>>) {
    let mut next = vec![None; unitig_count * 2];
    let mut prev = vec![None; unitig_count * 2];
    for candidate in candidates {
        let source = state_id(candidate.source_id, candidate.source_forward);
        let target = state_id(candidate.target_id, candidate.target_forward);
        if next[source].is_some() || prev[target].is_some() {
            continue;
        }
        next[source] = Some(target);
        prev[target] = Some(source);
    }
    (next, prev)
}

fn rewind_join_start(mut state: usize, prev: &[Option<usize>]) -> usize {
    let mut seen = HashSet::new();
    while let Some(prev_state) = prev[state] {
        if !seen.insert(prev_state) {
            break;
        }
        state = prev_state;
    }
    state
}

fn joined_sequence_from_state(
    start: usize,
    paths: &[UnitigPath],
    next: &[Option<usize>],
    used: &mut [bool],
    k: usize,
) -> Vec<u8> {
    let mut state = start;
    let (unitig_id, forward) = state_parts(state);
    let mut seq = oriented_sequence(&paths[unitig_id].seq, forward);

    loop {
        let (unitig_id, _) = state_parts(state);
        used[unitig_id] = true;
        let Some(next_state) = next[state] else {
            break;
        };
        let (next_unitig_id, next_forward) = state_parts(next_state);
        if used[next_unitig_id] {
            break;
        }
        let next_seq = oriented_sequence(&paths[next_unitig_id].seq, next_forward);
        seq.extend_from_slice(&next_seq[k - 1..]);
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

fn join_target(
    end: u64,
    starts: &[OrientedStart],
    topology: &GraphTopology<'_>,
) -> Option<OrientedStart> {
    let state = topology.state(end);
    if state.out_degree != 1 {
        return None;
    }
    let next = state.out_neighbor?;
    if topology.in_degree(next.oriented) != 1 {
        return None;
    }
    find_unique_oriented_start(starts, next.oriented)
}

fn find_unique_oriented_start(starts: &[OrientedStart], word: u64) -> Option<OrientedStart> {
    let mut idx = starts
        .binary_search_by_key(&word, |start| start.word)
        .ok()?;
    while idx > 0 && starts[idx - 1].word == word {
        idx -= 1;
    }
    if idx + 1 < starts.len() && starts[idx + 1].word == word {
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

fn oriented_sequence(seq: &[u8], forward: bool) -> Vec<u8> {
    if forward {
        seq.to_vec()
    } else {
        reverse_complement_seq(seq)
    }
}

fn build_kmer_locs(
    paths: &[UnitigPath],
    sorted_kmers: &[u64],
    lookup: &KmerLookup,
    k: usize,
    threads: usize,
) -> Result<Vec<KmerLoc>> {
    if threads > 1 && paths.len() > 1 && sorted_kmers.len() >= PARALLEL_TRINFO_COLLECTION_MIN {
        let out = (0..sorted_kmers.len())
            .map(|_| AtomicU64::new(0))
            .collect::<Vec<_>>();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| {
                Error::InvalidFormat(format!("failed to create index worker pool: {err}"))
            })?;
        pool.install(|| {
            paths
                .par_iter()
                .enumerate()
                .try_for_each(|(unitig_id, path)| {
                    scan_kmer_pairs(&path.seq, k, |offset, fwd, rev| {
                        let canonical = fwd.min(rev);
                        let Some(loc_idx) = lookup.find(sorted_kmers, canonical) else {
                            return Err(Error::InvalidFormat(
                                "joined unitig contains unknown k-mer".into(),
                            ));
                        };
                        let loc = KmerLoc::new(unitig_id, offset, fwd <= rev)?;
                        out[loc_idx].store(loc.0, Ordering::Relaxed);
                        Ok(())
                    })
                })
        })?;
        return Ok(out
            .into_iter()
            .map(|loc| KmerLoc(loc.into_inner()))
            .collect());
    }

    let mut out = vec![KmerLoc::default(); sorted_kmers.len()];
    for (unitig_id, path) in paths.iter().enumerate() {
        scan_kmer_pairs(&path.seq, k, |offset, fwd, rev| {
            let canonical = fwd.min(rev);
            let Some(loc_idx) = lookup.find(sorted_kmers, canonical) else {
                return Err(Error::InvalidFormat(
                    "joined unitig contains unknown k-mer".into(),
                ));
            };
            out[loc_idx] = KmerLoc::new(unitig_id, offset, fwd <= rev)?;
            Ok(())
        })?;
    }
    Ok(out)
}

fn build_trinfos(
    transcripts: &[Transcript],
    k: usize,
    sorted_kmers: &[u64],
    lookup: &KmerLookup,
    kmer_locs: &[KmerLoc],
    unitig_count: usize,
    threads: usize,
) -> Result<Vec<Vec<TRInfo>>> {
    let mut trinfos = vec![Vec::new(); unitig_count];
    let total_positions = transcripts
        .iter()
        .map(|tr| kmer_position_count(tr.seq.len(), k))
        .sum::<usize>();

    if threads > 1
        && transcripts.len() > 1
        && total_positions >= PARALLEL_TRINFO_COLLECTION_MIN
        && unitig_count > 0
    {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| {
                Error::InvalidFormat(format!("failed to create index worker pool: {err}"))
            })?;
        let shards = pool.install(|| {
            transcripts
                .par_iter()
                .enumerate()
                .map(|(tr_id, transcript)| {
                    collect_transcript_trinfos(
                        tr_id as u32,
                        transcript,
                        k,
                        sorted_kmers,
                        lookup,
                        kmer_locs,
                    )
                })
                .collect::<Result<Vec<_>>>()
        })?;
        for shard in shards {
            for (unitig_id, info) in shard {
                trinfos[unitig_id].push(info);
            }
        }
        return Ok(trinfos);
    }

    for (tr_id, transcript) in transcripts.iter().enumerate() {
        scan_transcript_trinfos(
            tr_id as u32,
            transcript,
            k,
            sorted_kmers,
            lookup,
            kmer_locs,
            |unitig_id, info| {
                trinfos[unitig_id].push(info);
            },
        )?;
    }

    Ok(trinfos)
}

fn collect_transcript_trinfos(
    trid: u32,
    transcript: &Transcript,
    k: usize,
    sorted_kmers: &[u64],
    lookup: &KmerLookup,
    kmer_locs: &[KmerLoc],
) -> Result<Vec<(usize, TRInfo)>> {
    let mut out = Vec::new();
    scan_transcript_trinfos(
        trid,
        transcript,
        k,
        sorted_kmers,
        lookup,
        kmer_locs,
        |unitig_id, info| {
            out.push((unitig_id, info));
        },
    )?;
    Ok(out)
}

fn scan_transcript_trinfos<F>(
    trid: u32,
    transcript: &Transcript,
    k: usize,
    sorted_kmers: &[u64],
    lookup: &KmerLookup,
    kmer_locs: &[KmerLoc],
    mut emit: F,
) -> Result<()>
where
    F: FnMut(usize, TRInfo),
{
    if transcript.seq.len() < k {
        return Ok(());
    }
    let mut run: Option<Run> = None;
    scan_kmer_pairs(&transcript.seq, k, |pos, fwd, rev| {
        let canonical = fwd.min(rev);
        let Some(loc_idx) = lookup.find(sorted_kmers, canonical) else {
            flush_run(trid, &mut run, &mut emit);
            return Ok(());
        };
        let loc = kmer_locs[loc_idx];
        let transcript_is_canonical = fwd <= rev;
        let same_strand = transcript_is_canonical == loc.unitig_is_canonical();
        let unitig_id = loc.unitig_id();
        let offset = loc.offset();
        if let Some(active) = run.as_mut()
            && active.unitig_id == unitig_id
            && active.same_strand == same_strand
        {
            let extends = if same_strand {
                offset == active.last_offset + 1
            } else {
                offset + 1 == active.last_offset
            };
            if extends {
                if same_strand {
                    active.stop = offset + 1;
                } else {
                    active.start = offset;
                }
                active.last_offset = offset;
                return Ok(());
            }
        }
        flush_run(trid, &mut run, &mut emit);
        run = Some(Run {
            unitig_id,
            tr_start_pos: pos as u32,
            start: offset,
            stop: offset + 1,
            last_offset: offset,
            same_strand,
        });
        Ok(())
    })?;
    flush_run(trid, &mut run, &mut emit);
    Ok(())
}

fn flush_run<F>(trid: u32, run: &mut Option<Run>, emit: &mut F)
where
    F: FnMut(usize, TRInfo),
{
    let Some(run) = run.take() else {
        return;
    };
    let pos = if run.same_strand {
        run.tr_start_pos
    } else {
        run.tr_start_pos | 0x8000_0000
    };
    emit(
        run.unitig_id,
        TRInfo {
            trid,
            pos,
            start: run.start as u32,
            stop: run.stop as u32,
        },
    );
}

fn build_node_payloads(
    paths: &[UnitigPath],
    mut trinfos: Vec<Vec<TRInfo>>,
    ec_threshold: usize,
    k: usize,
) -> Vec<NodePayload> {
    let mut nodes = Vec::with_capacity(paths.len());
    for (id, path) in paths.iter().enumerate() {
        let kmer_count = path.seq.len().saturating_sub(k) + 1;
        let infos = trinfos
            .get_mut(id)
            .map(Vec::as_mut_slice)
            .unwrap_or(&mut []);
        let blocks = build_blocks_for_unitig(infos, kmer_count as u32, ec_threshold);
        nodes.push(NodePayload {
            id: id as u32,
            blocks,
        });
    }
    nodes
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
            positions: Vec::new(),
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
        let (ec, positions) = if active.len() <= ec_threshold {
            active
                .iter()
                .map(|(&tr, pos)| (tr, pos.keys().copied().collect()))
                .unzip()
        } else {
            (Vec::new(), Vec::new())
        };
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
            positions: Vec::new(),
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
    trinfos.sort_by_key(|info| (info.trid, info.pos, info.start, info.stop));
    let mut blocks = Vec::new();
    for win in breakpoints.windows(2) {
        let lb = win[0];
        let ub = win[1];
        if lb == ub {
            continue;
        }
        let mut by_tr: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        for info in trinfos.iter() {
            if info.start <= lb && info.stop >= ub {
                by_tr.entry(info.trid).or_default().push(info.pos);
            }
        }
        let (ec, positions) = if by_tr.len() <= ec_threshold {
            by_tr
                .into_iter()
                .map(|(tr, mut pos)| {
                    pos.sort_unstable();
                    pos.dedup();
                    (tr, pos)
                })
                .unzip()
        } else {
            (Vec::new(), Vec::new())
        };
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
            positions: Vec::new(),
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
) -> Result<MinimizerIndexParts> {
    let mut keys = Vec::new();
    let total_positions = unitigs.iter().try_fold(0usize, |acc, seq| {
        acc.checked_add(seq.len().saturating_sub(g) + 1)
            .ok_or_else(|| Error::InvalidFormat("too many minimizer positions".into()))
    })?;
    let bitmap_count = (total_positions as u64 >> 32) as usize + 1;
    let mut bitmap_blocks = vec![Vec::new(); bitmap_count.max(1)];
    let mut global = 0u64;

    for seq in unitigs {
        insert_selected_minimizers(seq, k, g, global, &mut keys, &mut bitmap_blocks)?;
        global += minimizer_position_count(seq.len(), g) as u64;
    }

    for block in &mut bitmap_blocks {
        block.sort_unstable();
        block.dedup();
    }

    keys.sort_unstable();
    keys.dedup();
    Ok((keys, bitmap_blocks))
}

fn insert_selected_minimizers(
    seq: &[u8],
    k: usize,
    g: usize,
    global_base: u64,
    keys: &mut Vec<[u8; 8]>,
    bitmap_blocks: &mut [Vec<u32>],
) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    let Some(bounds) = bifrost_minimizer_bounds(k, g) else {
        return Ok(());
    };
    let mut hashes = Vec::with_capacity(seq.len() - g + 1);
    for pos in 0..=seq.len() - g {
        hashes.push(minimizer_rep_hash(&seq[pos..pos + g])?);
    }

    // Keep the window minimum with a monotonic queue instead of rescanning every k-mer.
    let mut selected_positions = Vec::new();
    let mut window = VecDeque::with_capacity(bounds.last - bounds.start + 1);
    let mut next_pos = bounds.start;
    for kmer_start in 0..=seq.len() - k {
        let window_start = kmer_start + bounds.start;
        let window_end = kmer_start + bounds.last;
        while window.front().is_some_and(|&pos| pos < window_start) {
            window.pop_front();
        }
        while next_pos <= window_end {
            let hash = hashes[next_pos];
            while window.back().is_some_and(|&pos| hashes[pos] > hash) {
                window.pop_back();
            }
            window.push_back(next_pos);
            next_pos += 1;
        }
        let Some(&first) = window.front() else {
            continue;
        };
        let min_hash = hashes[first];
        selected_positions.extend(
            window
                .iter()
                .copied()
                .take_while(|&pos| hashes[pos] == min_hash),
        );
    }

    selected_positions.sort_unstable();
    selected_positions.dedup();
    for pos in selected_positions {
        let bytes = encode_minimizer_rep(&seq[pos..pos + g])
            .ok_or_else(|| Error::InvalidFormat("invalid minimizer".into()))?;
        add_minimizer_position(global_base + pos as u64, bytes, keys, bitmap_blocks);
    }
    Ok(())
}

fn add_minimizer_position(
    pos: u64,
    bytes: [u8; 8],
    keys: &mut Vec<[u8; 8]>,
    bitmap_blocks: &mut [Vec<u32>],
) {
    keys.push(bytes);
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

fn minimizer_rep_hash(seq: &[u8]) -> Result<u64> {
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
    let mut hashes = [forward, reverse];
    if hashes[1] < hashes[0] {
        hashes.swap(0, 1);
    }
    let mut buf = [0u8; 16];
    buf[..8].copy_from_slice(&hashes[0].to_le_bytes());
    buf[8..].copy_from_slice(&hashes[1].to_le_bytes());
    Ok(wyhash(&buf, 0))
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
    for pos in 1..=seq.len() - k {
        let code = base_code(seq[pos + k - 1])?;
        fwd = shift_append_word(fwd, code, k, mask);
        rev = shift_prepend_word(rev, 3 - code, mask);
        f(pos, fwd, rev)?;
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

fn decode_kmer_word(word: u64, k: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(k);
    for i in 0..k {
        let shift = 62 - ((i & 0x1f) << 1);
        let code = (word >> shift) & 3;
        out.push(match code {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}

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

    fn contains_word(&self, entries: &[u64], lookup: &KmerLookup, word: u64) -> Result<bool> {
        let idx = lookup
            .find(entries, word)
            .ok_or_else(|| Error::InvalidFormat("graph traversal found unknown k-mer".into()))?;
        Ok(self.contains_index(idx))
    }

    fn insert_word(
        &mut self,
        entries: &[u64],
        lookup: &KmerLookup,
        word: u64,
    ) -> Result<Option<usize>> {
        let idx = lookup
            .find(entries, word)
            .ok_or_else(|| Error::InvalidFormat("graph traversal found unknown k-mer".into()))?;
        if self.contains_index(idx) {
            return Ok(None);
        }
        self.words[idx >> 6] |= 1u64 << (idx & 63);
        Ok(Some(idx))
    }

    fn contains_index(&self, idx: usize) -> bool {
        self.words[idx >> 6] & (1u64 << (idx & 63)) != 0
    }
}

impl KmerLookup {
    fn new(entries: &[u64]) -> Result<Self> {
        if entries.len() > u32::MAX as usize {
            return Err(Error::InvalidFormat(
                "too many distinct k-mers for compact lookup".into(),
            ));
        }
        // Prefix buckets keep the table compact and make small buckets linear-scan friendly.
        let prefix_bits = kmer_lookup_prefix_bits(entries.len());
        let bucket_count = 1usize << prefix_bits;
        let shift = 64 - prefix_bits;
        let mut offsets = vec![0u32; bucket_count + 1];
        for &word in entries {
            offsets[((word >> shift) as usize) + 1] += 1;
        }
        for idx in 1..offsets.len() {
            offsets[idx] += offsets[idx - 1];
        }
        Ok(Self { offsets, shift })
    }

    fn find(&self, entries: &[u64], word: u64) -> Option<usize> {
        let bucket = (word >> self.shift) as usize;
        let start = self.offsets[bucket] as usize;
        let end = self.offsets[bucket + 1] as usize;
        let slice = &entries[start..end];
        if slice.len() <= LINEAR_KMER_LOOKUP_MAX {
            slice
                .iter()
                .position(|&entry| entry == word)
                .map(|idx| start + idx)
        } else {
            slice.binary_search(&word).ok().map(|idx| start + idx)
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

impl EdgeStateCache {
    fn build(entries: &[u64], lookup: &KmerLookup, k: usize, threads: usize) -> Result<Self> {
        let mask = kmer_mask(k);
        if threads > 1 && entries.len() >= PARALLEL_KMER_COLLECTION_MIN {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|err| {
                    Error::InvalidFormat(format!("failed to create index worker pool: {err}"))
                })?;
            let packed = pool.install(|| {
                entries
                    .par_iter()
                    .map(|&word| pack_edge_states(word, entries, lookup, k, mask))
                    .collect()
            });
            return Ok(Self { packed });
        }

        Ok(Self {
            packed: entries
                .iter()
                .map(|&word| pack_edge_states(word, entries, lookup, k, mask))
                .collect(),
        })
    }

    fn state_bits(&self, idx: usize, forward: bool) -> u32 {
        let shift = if forward { 0 } else { EDGE_STATE_BITS };
        (self.packed[idx] >> shift) & ((1u32 << EDGE_STATE_BITS) - 1)
    }
}

fn pack_edge_states(word: u64, entries: &[u64], lookup: &KmerLookup, k: usize, mask: u64) -> u32 {
    let forward = pack_edge_state(word, word, entries, lookup, k, mask);
    let reverse_word = reverse_complement_word(word, k, mask);
    let reverse = pack_edge_state(reverse_word, word, entries, lookup, k, mask);
    forward | (reverse << EDGE_STATE_BITS)
}

fn pack_edge_state(
    oriented: u64,
    canonical_word: u64,
    entries: &[u64],
    lookup: &KmerLookup,
    k: usize,
    mask: u64,
) -> u32 {
    let mut in_degree = 0u32;
    let mut out_degree = 0u32;
    let mut out_code = 0u32;
    let mut has_non_self_out = false;

    for code in 0..4 {
        let next = shift_append_word(oriented, code, k, mask);
        let next_canonical = next.min(reverse_complement_word(next, k, mask));
        if lookup.find(entries, next_canonical).is_some() {
            out_degree += 1;
            out_code = code as u32;
            has_non_self_out |= next_canonical != canonical_word;
        }

        let prev = shift_prepend_word(oriented, code, mask);
        let prev_canonical = prev.min(reverse_complement_word(prev, k, mask));
        if lookup.find(entries, prev_canonical).is_some() {
            in_degree += 1;
        }
    }

    in_degree
        | (out_degree << EDGE_OUT_DEGREE_SHIFT)
        | (out_code << EDGE_OUT_CODE_SHIFT)
        | (u32::from(has_non_self_out) << EDGE_HAS_NON_SELF_SHIFT)
}

impl<'a> GraphTopology<'a> {
    fn new(
        entries: &'a [u64],
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

    fn state(&self, oriented: u64) -> EdgeState {
        self.state_from_bits(oriented, true)
    }

    fn degree_state(&self, oriented: u64) -> EdgeState {
        self.state_from_bits(oriented, false)
    }

    fn in_degree(&self, oriented: u64) -> u8 {
        self.state_from_bits(oriented, false).in_degree
    }

    fn state_from_bits(&self, oriented: u64, include_neighbor: bool) -> EdgeState {
        let canonical = self.canonical(oriented);
        let idx = self
            .lookup
            .find(self.entries, canonical)
            .expect("state lookup for known k-mer");
        let forward = oriented == canonical;
        let bits = self.edge_states.state_bits(idx, forward);
        let in_degree = (bits & EDGE_IN_DEGREE_MASK) as u8;
        let out_degree = ((bits >> EDGE_OUT_DEGREE_SHIFT) & EDGE_IN_DEGREE_MASK) as u8;
        let out_neighbor = if include_neighbor && out_degree == 1 {
            let code = u64::from((bits >> EDGE_OUT_CODE_SHIFT) & 0x3);
            let oriented = shift_append_word(oriented, code, self.k, self.mask);
            let canonical = self.canonical(oriented);
            Some(OrientedNeighbor {
                canonical,
                oriented,
            })
        } else {
            None
        };
        EdgeState {
            in_degree,
            out_degree,
            out_neighbor,
            has_non_self_out: bits & (1u32 << EDGE_HAS_NON_SELF_SHIFT) != 0,
        }
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
    let mut word = (!word) & mask;
    word = ((word & 0x3333_3333_3333_3333) << 2) | ((word >> 2) & 0x3333_3333_3333_3333);
    word = ((word & 0x0f0f_0f0f_0f0f_0f0f) << 4) | ((word >> 4) & 0x0f0f_0f0f_0f0f_0f0f);
    word = ((word & 0x00ff_00ff_00ff_00ff) << 8) | ((word >> 8) & 0x00ff_00ff_00ff_00ff);
    word = ((word & 0x0000_ffff_0000_ffff) << 16) | ((word >> 16) & 0x0000_ffff_0000_ffff);
    word = word.rotate_left(32);
    (word << (64 - (2 * k))) & mask
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
    fn kmer_roundtrip() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let word = encode_kmer(seq).expect("encode");
        assert_eq!(decode_kmer_word(word, seq.len()), seq);
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
    fn block_sweep_matches_scan_builder() {
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
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let k = 31;
        let g = 23;
        let mut keys = Vec::new();
        let mut bitmap_blocks = vec![Vec::new()];
        insert_selected_minimizers(seq, k, g, 0, &mut keys, &mut bitmap_blocks)
            .expect("sliding minimizers");

        let mut actual = bitmap_blocks[0]
            .iter()
            .copied()
            .zip(keys)
            .collect::<Vec<_>>();
        actual.sort_unstable();

        let mut expected = selected_minimizers_naive(seq, k, g);
        expected.sort_unstable();
        expected.dedup_by_key(|(pos, _)| *pos);
        assert_eq!(actual, expected);
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
    fn compacts_linear_transcript() {
        let tr = Transcript {
            name: "tx".into(),
            original_len: 36,
            seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        };
        let (graph, _) = build_kmer_unitig_graph_with_report(&[tr], 31, 23, -1, 1).expect("graph");
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
        let (graph, _) = build_kmer_unitig_graph_with_report(&[tr], 31, 23, -1, 1).expect("graph");
        assert_eq!(graph.unitigs.len(), 1);
        assert_eq!(graph.km_unitigs.len(), 0);
        assert_eq!(graph.nodes.len(), graph.unitigs.len());
    }
}
