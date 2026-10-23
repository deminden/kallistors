//! Unitig path joining and exact overlap validation.

use super::*;

// Canonical claims bound the path count by the already-validated u32
// k-mer count. Supported k<32 leaves a low word bit for orientation.
#[repr(C, packed)]
#[derive(Clone, Copy)]
pub(super) struct OrientedStart {
    pub(super) word: u64,
    pub(super) unitig: u32,
}
impl OrientedStart {
    pub(super) fn word(self) -> u64 {
        self.word & !1
    }
    pub(super) fn unitig_id(self) -> usize {
        self.unitig as usize
    }
    pub(super) fn forward(self) -> bool {
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
pub(super) struct OrientedStartLookup<'a> {
    pub(super) starts: &'a [OrientedStart],
    pub(super) offsets: Vec<u32>,
    pub(super) shift: u32,
}
impl<'a> OrientedStartLookup<'a> {
    pub(super) fn new(starts: &'a [OrientedStart], threads: usize) -> Self {
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
    pub(super) fn find(&self, word: u64) -> Option<OrientedStart> {
        let starts = if self.offsets.is_empty() {
            self.starts
        } else {
            let bucket = (word >> self.shift) as usize;
            &self.starts[self.offsets[bucket] as usize..self.offsets[bucket + 1] as usize]
        };
        find_unique_oriented_start(starts, word)
    }
}

pub(super) enum JoinLinks {
    Narrow(Vec<u32>),
    Wide(Vec<usize>),
}

impl JoinLinks {
    pub(super) fn retire(self) {
        match self {
            Self::Narrow(values) => super::super::scratch_pages::retire(values),
            Self::Wide(values) => super::super::scratch_pages::retire(values),
        }
    }

    pub(super) fn new(count: usize) -> Self {
        if count <= u32::MAX as usize {
            Self::Narrow(vec![u32::MAX; count])
        } else {
            Self::Wide(vec![usize::MAX; count])
        }
    }
    #[inline]
    pub(super) fn get(&self, state: usize) -> Option<usize> {
        match self {
            Self::Narrow(values) => (values[state] != u32::MAX).then_some(values[state] as usize),
            Self::Wide(values) => (values[state] != usize::MAX).then_some(values[state]),
        }
    }
    #[inline]
    pub(super) fn set(&mut self, state: usize, next: usize) {
        match self {
            Self::Narrow(values) => values[state] = next as u32,
            Self::Wide(values) => values[state] = next,
        }
    }
}

pub(super) enum JoinCandidates {
    Narrow(Vec<[u32; 2]>),
    Wide(Vec<[usize; 2]>),
}

#[cfg(test)]
pub(super) fn join_unitig_paths(
    paths: Vec<UnitigPath>,
    topology: &GraphTopology<'_, impl MappingKeys + ?Sized>,
    threads: usize,
) -> Result<Vec<UnitigPath>> {
    let queries = collect_join_queries(&paths, topology, threads)?;
    join_unitig_paths_with_queries(paths, queries, topology.k, threads)
}

#[repr(C)]
pub(super) struct JoinQuery {
    pub(super) word: [u8; 8],
    pub(super) source: [u8; 5],
}
impl JoinQuery {
    pub(super) fn new(source: usize, word: u64) -> Result<Self> {
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
    pub(super) fn word(&self) -> u64 {
        u64::from_ne_bytes(self.word)
    }
    pub(super) fn source(&self) -> usize {
        let mut bytes = [0; 8];
        bytes[..5].copy_from_slice(&self.source);
        u64::from_le_bytes(bytes) as usize
    }
}

pub(super) fn collect_join_queries(
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

pub(super) fn join_unitig_paths_with_queries(
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
    super::super::scratch_pages::retire(starts);
    super::super::scratch_pages::retire(active);
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

    super::super::scratch_pages::retire(paths);
    next.retire();
    prev.retire();
    super::super::scratch_pages::retire(used);
    Ok(out)
}

pub(super) fn build_oriented_starts(
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

pub(super) fn resolve_join_queries(
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
pub(super) fn resolve_join_queries_typed<State: Copy + Send + Sync>(
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

pub(super) fn build_join_edges(
    unitig_count: usize,
    candidates: JoinCandidates,
) -> (JoinLinks, JoinLinks) {
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
            super::super::scratch_pages::retire(values);
        }
        JoinCandidates::Wide(values) => {
            for &[source, target] in &values {
                insert(source, target);
            }
            super::super::scratch_pages::retire(values);
        }
    }
    (next, prev)
}

// Valid states occupy0..2*unitig_count. An allocated state array cannot
// reach usize::MAX, leaving an exact out-of-domain missing-link sentinel.
#[inline]
pub(super) fn join_state(states: &JoinLinks, state: usize) -> Option<usize> {
    states.get(state)
}

pub(super) fn rewind_join_start(mut state: usize, prev: &JoinLinks) -> usize {
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

pub(super) fn joined_sequence_from_state(
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

pub(super) fn state_id(unitig_id: usize, forward: bool) -> usize {
    unitig_id * 2 + usize::from(!forward)
}

pub(super) fn state_parts(state: usize) -> (usize, bool) {
    (state / 2, state.is_multiple_of(2))
}

pub(super) fn find_unique_oriented_start(
    starts: &[OrientedStart],
    word: u64,
) -> Option<OrientedStart> {
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

pub(super) fn oriented_end_word(seq: &[u8], forward: bool, k: usize, mask: u64) -> Result<u64> {
    if forward {
        encode_kmer(&seq[seq.len() - k..])
    } else {
        Ok(reverse_complement_word(encode_kmer(&seq[..k])?, k, mask))
    }
}
