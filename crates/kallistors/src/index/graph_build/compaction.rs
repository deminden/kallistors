//! Ordered unitig compaction, endpoint claims and graph walks.

use super::*;

#[derive(Clone)]
pub(super) struct UnitigPath {
    pub(super) seq: Box<[u8]>,
    pub(super) ends: [u32; 2],
}

pub(super) fn compact_unitigs(
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

pub(super) fn finish_unitig_cycles(
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

pub(super) trait PathStep: Copy + Send + Sync {
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

pub(super) struct PreparedPath<Step: PathStep> {
    pub(super) start: u64,
    pub(super) steps: Vec<Step>,
}

pub(super) fn compact_unitigs_parallel(
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
pub(super) fn concatenate_shards<T>(shards: Vec<Vec<T>>) -> Vec<T> {
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
pub(super) struct PathChunks {
    pub(super) chunks: Vec<Vec<UnitigPath>>,
    pub(super) len: usize,
}
impl PathChunks {
    pub(super) fn push(&mut self, path: UnitigPath) {
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
    pub(super) fn extend(&mut self, paths: impl IntoIterator<Item = UnitigPath>) {
        for path in paths {
            self.push(path);
        }
    }
    pub(super) fn finish(self) -> Vec<UnitigPath> {
        let mut paths = Vec::with_capacity(self.len);
        for mut chunk in self.chunks {
            let written = chunk.len();
            paths.append(&mut chunk);
            super::super::scratch_pages::discard_unused(&mut chunk, written);
        }
        paths
    }
}

pub(super) fn compact_unitigs_parallel_steps<Step: PathStep>(
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
    super::super::scratch_pages::retire(starts);
    super::super::scratch_pages::retire(buffers);
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

pub(super) fn claim_prepared_path<Step: PathStep>(
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
pub(super) fn materialize_prepared_path<Step: PathStep>(
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
pub(super) fn commit_prepared_path<Step: PathStep>(
    path: &PreparedPath<Step>,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    visited: &mut VisitedKmers,
) -> Option<UnitigPath> {
    let count = claim_prepared_path(path, visited);
    (count != 0).then(|| materialize_prepared_path(path, topology.k, count))
}

pub(super) fn finish_unitig_cycles_parallel<Step: PathStep>(
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

pub(super) fn finish_unitig_cycles_batched<Step: PathStep>(
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

pub(super) struct PreparingPath<Step: PathStep> {
    pub(super) path: PreparedPath<Step>,
    pub(super) start_idx: usize,
    pub(super) slot: usize,
    pub(super) index: usize,
    pub(super) forward: bool,
    pub(super) base: u64,
}

pub(super) fn prepare_path_group<Step: PathStep, const RESIDUAL: bool>(
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

pub(super) fn prepare_path_group_inner<
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

pub(super) fn choose_branch_start_orientation(
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

pub(super) fn choose_any_start_orientation(
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

pub(super) fn walk_unitig(
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
