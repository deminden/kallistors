//! K-mer locations and packed transcript sequences.

use super::*;

pub(super) struct DenseReferenceLocations<'a, K: MappingKeys + ?Sized> {
    pub(super) keys: &'a K,
    pub(super) lookup: &'a KmerLookup,
    pub(super) locations: &'a KmerLocations,
}
impl<K: MappingKeys + ?Sized> ReferenceKmerLocations for DenseReferenceLocations<'_, K> {
    #[inline]
    fn find(&self, word: u64) -> Option<KmerLoc> {
        self.keys
            .find(self.lookup, word)
            .map(|rank| self.locations.get(rank))
    }
}

pub(super) enum KmerLocations {
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
    pub(super) fn get(&self, index: usize) -> KmerLoc {
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

pub(super) fn build_kmer_locs<K: MappingKeys + ?Sized>(
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
pub(super) fn prefetch_location<T>(slot: &T) {
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
pub(super) fn store_unitig_locations<K: MappingKeys + ?Sized, P, F>(
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
pub(super) fn store_location_group<const N: usize, K: MappingKeys + ?Sized, P, F>(
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

pub(super) fn pack_transcript_sequences(
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
