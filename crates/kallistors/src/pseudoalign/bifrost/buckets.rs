//! Minimizer bucket counts, fingerprints and position population.

use super::loader_probe::CachedUnitigMinimizers;
use super::*;
mod count_merge;
use count_merge::merge_count_chunks;
mod cached_counts;
use crate::index::bifrost::MphfBatchKeys;

pub(super) fn precompute_encoded_unitigs(
    unitigs: &[Vec<u8>],
    km_unitigs: &[Vec<u8>],
    g: usize,
    pool: Option<&ThreadPool>,
) -> EncodedUnitigs {
    EncodedUnitigs {
        g,
        unitig_bases: EncodedBaseStore::from_sequences(unitigs, pool),
        km_bases: EncodedBaseStore::from_sequences(km_unitigs, pool),
    }
}

#[derive(Clone, Copy)]
pub(super) struct UnitigBlockState {
    pub(super) unitig_id: usize,
    pub(super) total_unitig_len: u64,
    pub(super) current_unitig_len: u64,
}

pub(super) fn advance_unitig_state(
    state: &mut UnitigBlockState,
    encoded_unitigs: &EncodedUnitigs,
    pos: u64,
) {
    while pos >= state.total_unitig_len + state.current_unitig_len {
        state.unitig_id += 1;
        state.total_unitig_len += state.current_unitig_len;
        state.current_unitig_len = encoded_unitigs.unitig_minimizer_count(state.unitig_id) as u64;
    }
}

pub(super) fn precompute_unitig_block_states(
    bmp_unitigs: &[BitmapPositions],
    encoded_unitigs: &EncodedUnitigs,
) -> Vec<UnitigBlockState> {
    let mut out = Vec::with_capacity(bmp_unitigs.len());
    if encoded_unitigs.unitig_bases.metadata.is_empty() {
        return out;
    }

    let mut state = UnitigBlockState {
        unitig_id: 0,
        total_unitig_len: 0,
        current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
    };
    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        advance_unitig_state(&mut state, encoded_unitigs, base);
        out.push(state);
        if let Some(last_pos) = bitmap.last() {
            advance_unitig_state(&mut state, encoded_unitigs, base + last_pos as u64);
        }
    }
    out
}

pub(super) fn populate_unitig_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_unitigs: &[BitmapPositions],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
    bucket_cache: Option<&mut Vec<u32>>,
    fingerprints: &[AtomicU16],
) {
    if fingerprints.is_empty() {
        populate_unitig_counts_typed::<u32>(
            minz_counts,
            bmp_unitigs,
            encoded_unitigs,
            mphf,
            pool,
            bucket_cache,
            fingerprints,
        );
    } else {
        populate_unitig_counts_typed::<TaggedBucket>(
            minz_counts,
            bmp_unitigs,
            encoded_unitigs,
            mphf,
            pool,
            bucket_cache,
            fingerprints,
        );
    }
}

pub(super) fn populate_unitig_counts_typed<K: BucketCountKey>(
    minz_counts: &mut [usize],
    bmp_unitigs: &[BitmapPositions],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
    bucket_cache: Option<&mut Vec<u32>>,
    fingerprints: &[AtomicU16],
) {
    // Larger reference tables benefit from more independent outstanding probes.
    // Select once per load; ordinary references retain eight-key batches.
    if mphf.size() >= (1 << 28) {
        populate_unitig_counts_batched::<K, 64, true>(
            minz_counts,
            bmp_unitigs,
            encoded_unitigs,
            mphf,
            pool,
            bucket_cache,
            fingerprints,
        );
    } else {
        populate_unitig_counts_batched::<K, 8, false>(
            minz_counts,
            bmp_unitigs,
            encoded_unitigs,
            mphf,
            pool,
            bucket_cache,
            fingerprints,
        );
    }
}

fn populate_unitig_counts_batched<K: BucketCountKey, const BATCH: usize, const VIEWS: bool>(
    minz_counts: &mut [usize],
    bmp_unitigs: &[BitmapPositions],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
    mut bucket_cache: Option<&mut Vec<u32>>,
    fingerprints: &[AtomicU16],
) {
    if bmp_unitigs.is_empty() || encoded_unitigs.unitig_bases.metadata.is_empty() {
        return;
    }

    if mphf.size() >= (1 << 28)
        && bmp_unitigs.len() == 1
        && let Some(pool) = pool
        && pool.current_num_threads() > 1
        && let Some(cache) = bucket_cache.as_mut()
    {
        cached_counts::populate::<K, BATCH, VIEWS>(
            minz_counts,
            &bmp_unitigs[0],
            encoded_unitigs,
            mphf,
            pool,
            cache,
            fingerprints,
        );
        return;
    }

    let Some(pool) = pool else {
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        let mut probe = CachedUnitigMinimizers::new(encoded_unitigs);
        let mut keys = MphfBatchKeys::<BATCH>::new();
        let mut used = 0;
        let mut flush = |keys: &MphfBatchKeys<BATCH>, used: usize| {
            let buckets = mphf.lookup_batch(keys);
            // Lookup is independent, but cache entries, fingerprints and counts
            // retain the exact serialized bitmap order, including absent keys.
            for (bucket, key) in buckets.into_iter().zip(keys.iter()).take(used) {
                if let (Some(bucket), Some(bytes)) = (bucket, key)
                    && let Some(slot) = fingerprints.get(bucket as usize)
                {
                    super::super::fingerprints::observe(slot, u64::from_le_bytes(bytes));
                }
                if let Some(cache) = bucket_cache.as_mut() {
                    cache.push(bucket.map_or(u32::MAX, |bucket| bucket as u32));
                }
                if let Some(bucket) = bucket {
                    minz_counts[bucket as usize] += 1;
                }
            }
        };
        for (i, bitmap) in bmp_unitigs.iter().enumerate() {
            let base = (i as u64) << 32;
            for pos_bmp in bitmap.iter() {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                keys.set(
                    used,
                    if VIEWS {
                        probe.at(state.unitig_id, rel)
                    } else {
                        encoded_unitigs.unitig_minimizer_at(state.unitig_id, rel)
                    },
                );
                used += 1;
                if used == BATCH {
                    flush(&keys, used);
                    keys.clear();
                    used = 0;
                }
            }
        }
        if used != 0 {
            flush(&keys, used);
        }
        return;
    };

    let states = precompute_unitig_block_states(bmp_unitigs, encoded_unitigs);
    let mut batch_pos_ids = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut flush_batch = |batch_pos_ids: &mut Vec<u64>| {
        if batch_pos_ids.is_empty() {
            return;
        }
        if let Some(cache) = bucket_cache.as_mut() {
            let chunk_len = (batch_pos_ids.len() / (pool.current_num_threads() * 4).max(1))
                .max(LOAD_MIN_PAR_CHUNK);
            let old_len = cache.len();
            cache.reserve(batch_pos_ids.len());
            let destinations = &mut cache.spare_capacity_mut()[..batch_pos_ids.len()];
            let chunks: Vec<Vec<K>> = pool.install(|| {
                destinations
                    .par_chunks_mut(chunk_len)
                    .zip(batch_pos_ids.par_chunks(chunk_len))
                    .map(|(destinations, positions)| {
                        let mut probe = CachedUnitigMinimizers::new(encoded_unitigs);
                        let mut indices = Vec::with_capacity(positions.len());
                        for (destinations, positions) in
                            destinations.chunks_mut(BATCH).zip(positions.chunks(BATCH))
                        {
                            let mut keys = MphfBatchKeys::<BATCH>::new();
                            for (lane, &pos_id) in positions.iter().enumerate() {
                                keys.set(
                                    lane,
                                    if VIEWS {
                                        probe.at((pos_id >> 32) as usize, pos_id as u32 as usize)
                                    } else {
                                        encoded_unitigs.unitig_minimizer_at(
                                            (pos_id >> 32) as usize,
                                            pos_id as u32 as usize,
                                        )
                                    },
                                );
                            }
                            let buckets = mphf.lookup_batch(&keys);
                            for (destination, (bucket, key)) in destinations
                                .iter_mut()
                                .zip(buckets.into_iter().zip(keys.iter()))
                            {
                                destination.write(bucket.map_or(u32::MAX, |bucket| bucket as u32));
                                if let (Some(bucket), Some(key)) = (bucket, key) {
                                    indices.push(K::new(bucket as u32, u64::from_le_bytes(key)));
                                }
                            }
                        }
                        indices
                            .radix_sort_builder()
                            .with_parallel(false)
                            .with_single_threaded_tuner()
                            .sort();
                        indices
                    })
                    .collect()
            });
            // SAFETY: the indexed parallel zip writes every destination exactly once.
            // A panic unwinds before this publication and leaves the old length intact.
            unsafe {
                cache.set_len(old_len + batch_pos_ids.len());
            }
            merge_count_chunks(minz_counts, fingerprints, &chunks, pool);
        } else {
            let grouped = collect_unitig_batch_counts(batch_pos_ids, encoded_unitigs, mphf, pool);
            merge_count_batches(minz_counts, grouped);
        }
        batch_pos_ids.clear();
    };

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let mut state = states[i];
        let base = (i as u64) << 32;
        for pos_bmp in bitmap.iter() {
            let pos = base + pos_bmp as u64;
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            let rel = (pos - state.total_unitig_len) as usize;
            batch_pos_ids.push(((state.unitig_id as u64) << 32) | rel as u64);
            if batch_pos_ids.len() >= LOAD_BATCH_POSITIONS {
                flush_batch(&mut batch_pos_ids);
            }
        }
    }
    flush_batch(&mut batch_pos_ids);
}

pub(super) fn populate_unitig_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_unitigs: &[BitmapPositions],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
    bucket_cache: Option<&[u32]>,
) {
    if bmp_unitigs.is_empty() || encoded_unitigs.unitig_bases.metadata.is_empty() {
        return;
    }

    if let Some(cache) = bucket_cache {
        debug_assert_eq!(bmp_unitigs.len(), 1);
        debug_assert_eq!(cache.len(), bmp_unitigs[0].len());
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        for (_at, (pos, &bucket)) in bmp_unitigs[0].iter().zip(cache).enumerate() {
            #[cfg(target_arch = "aarch64")]
            {
                if let Some(&future) = cache.get(_at + 64)
                    && future != u32::MAX
                {
                    minz_positions.prefetch_inline_descriptor(future as usize);
                }
                if let Some(&future) = cache.get(_at + 16)
                    && future != u32::MAX
                {
                    minz_positions.prefetch_inline_overflow(future as usize);
                }
            }
            let pos = u64::from(pos);
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            if bucket != u32::MAX {
                let rel = pos - state.total_unitig_len;
                minz_positions.push_at(
                    minz_cursors,
                    bucket as usize,
                    ((state.unitig_id as u64) << 32) | rel,
                );
            }
        }
        return;
    }
    let Some(pool) = pool else {
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        for (i, bitmap) in bmp_unitigs.iter().enumerate() {
            let base = (i as u64) << 32;
            for pos_bmp in bitmap.iter() {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizer_at(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(&bytes)
                {
                    let pos_id = ((state.unitig_id as u64) << 32) | rel as u64;
                    minz_positions.push_at(minz_cursors, idx as usize, pos_id);
                }
            }
        }
        return;
    };

    if bmp_unitigs.len() == 1 {
        // A single bitmap otherwise gives the outer parallel loop only one task.
        // Resolve MPHF buckets concurrently, then append in the original order.
        // Bucket order affects which exact matching candidate is visited first.
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        let mut batch = Vec::with_capacity(LOAD_BATCH_POSITIONS);
        let mut flush = |batch: &mut Vec<u64>| {
            let chunk_len =
                (batch.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
            let events: Vec<Vec<(usize, u64)>> = pool.install(|| {
                batch
                    .par_chunks(chunk_len)
                    .map(|chunk| {
                        chunk
                            .iter()
                            .filter_map(|&pos_id| {
                                let uid = (pos_id >> 32) as usize;
                                let rel = (pos_id & 0xffff_ffff) as usize;
                                let bytes = encoded_unitigs.unitig_minimizer_at(uid, rel)?;
                                Some((mphf.lookup(&bytes)? as usize, pos_id))
                            })
                            .collect()
                    })
                    .collect()
            });
            for chunk in events {
                for (bucket, pos_id) in chunk {
                    minz_positions.push_at(minz_cursors, bucket, pos_id);
                }
            }
            batch.clear();
        };
        for pos in bmp_unitigs[0].iter() {
            let pos = u64::from(pos);
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            let rel = pos - state.total_unitig_len;
            batch.push(((state.unitig_id as u64) << 32) | rel);
            if batch.len() >= LOAD_BATCH_POSITIONS {
                flush(&mut batch);
            }
        }
        flush(&mut batch);
        return;
    }

    let states = precompute_unitig_block_states(bmp_unitigs, encoded_unitigs);
    let cursors = minz_cursors
        .iter()
        .copied()
        .map(AtomicUsize::new)
        .collect::<Vec<_>>();
    let values_ptr = minz_positions.values.as_mut_ptr() as usize;
    pool.install(|| {
        bmp_unitigs.par_iter().enumerate().for_each(|(i, bitmap)| {
            let mut state = states[i];
            let base = (i as u64) << 32;
            for pos_bmp in bitmap.iter() {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizer_at(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(&bytes)
                {
                    let pos_id = ((state.unitig_id as u64) << 32) | rel as u64;
                    let slot = cursors[idx as usize].fetch_add(1, Ordering::Relaxed);
                    // SAFETY: the precomputed bucket counts make each fetched slot unique,
                    // so parallel workers never write the same position.
                    unsafe {
                        (values_ptr as *mut u64).add(slot).write(pos_id);
                    }
                }
            }
        });
    });
    for (dst, src) in minz_cursors.iter_mut().zip(cursors) {
        *dst = src.into_inner();
    }
}

pub(super) fn populate_kmer_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_km: &[BitmapPositions],
    ctx: &KmerLoadContext<'_>,
    pool: Option<&ThreadPool>,
    fingerprints: &[AtomicU16],
) {
    let Some(pool) = pool else {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for pos_bmp in bitmap.iter() {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                    && let Some(idx) = ctx.mphf.lookup(&bytes)
                {
                    minz_counts[idx as usize] += 1;
                    if let Some(slot) = fingerprints.get(idx as usize) {
                        slot.store(super::super::fingerprints::DISABLED, Ordering::Relaxed);
                    }
                }
            }
        }
        return;
    };

    let mut batch_pos_ids = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut flush_batch = |batch_pos_ids: &mut Vec<u64>| {
        if batch_pos_ids.is_empty() {
            return;
        }
        let grouped = collect_kmer_batch_counts(batch_pos_ids, ctx, pool);
        if !fingerprints.is_empty() {
            for chunk in &grouped {
                for &idx in chunk {
                    fingerprints[idx]
                        .store(super::super::fingerprints::DISABLED, Ordering::Relaxed);
                }
            }
        }
        merge_count_batches(minz_counts, grouped);
        batch_pos_ids.clear();
    };

    for (i, bitmap) in bmp_km.iter().enumerate() {
        let base = (i as u64) << 32;
        for pos_bmp in bitmap.iter() {
            let pos = base + pos_bmp as u64;
            let km_id = (pos / ctx.km_glen as u64) as usize;
            let km_pos = (pos % ctx.km_glen as u64) as usize;
            batch_pos_ids.push(((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64);
            if batch_pos_ids.len() >= LOAD_BATCH_POSITIONS {
                flush_batch(&mut batch_pos_ids);
            }
        }
    }
    flush_batch(&mut batch_pos_ids);
}

pub(super) fn populate_kmer_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_km: &[BitmapPositions],
    ctx: &KmerLoadContext<'_>,
    pool: Option<&ThreadPool>,
) {
    // A single bitmap supplies only one outer Rayon task, so atomics do not
    // provide parallelism and duplicate the entire reference-sized cursor table.
    let Some(pool) = pool.filter(|_| bmp_km.len() > 1) else {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for pos_bmp in bitmap.iter() {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                    && let Some(idx) = ctx.mphf.lookup(&bytes)
                {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64;
                    minz_positions.push_at(minz_cursors, idx as usize, pos_id);
                }
            }
        }
        return;
    };

    let cursors = minz_cursors
        .iter()
        .copied()
        .map(AtomicUsize::new)
        .collect::<Vec<_>>();
    let values_ptr = minz_positions.values.as_mut_ptr() as usize;
    pool.install(|| {
        bmp_km.par_iter().enumerate().for_each(|(i, bitmap)| {
            let base = (i as u64) << 32;
            for pos_bmp in bitmap.iter() {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                    && let Some(idx) = ctx.mphf.lookup(&bytes)
                {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64;
                    let slot = cursors[idx as usize].fetch_add(1, Ordering::Relaxed);
                    // SAFETY: the precomputed bucket counts make each fetched slot unique,
                    // so parallel workers never write the same position.
                    unsafe {
                        (values_ptr as *mut u64).add(slot).write(pos_id);
                    }
                }
            }
        });
    });
    for (dst, src) in minz_cursors.iter_mut().zip(cursors) {
        *dst = src.into_inner();
    }
}

pub(super) fn collect_unitig_batch_counts(
    pos_ids: &[u64],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: &ThreadPool,
) -> Vec<Vec<usize>> {
    let chunk_len =
        (pos_ids.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    pool.install(|| {
        pos_ids
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut indices = Vec::with_capacity(chunk.len());
                for &pos_id in chunk {
                    let unitig_id = (pos_id >> 32) as usize;
                    let rel = (pos_id & 0xffff_ffff) as usize;
                    if let Some(bytes) = encoded_unitigs.unitig_minimizer_at(unitig_id, rel)
                        && let Some(idx) = mphf.lookup(&bytes)
                    {
                        indices.push(idx as usize);
                    }
                }
                indices
                    .radix_sort_builder()
                    .with_parallel(false)
                    .with_single_threaded_tuner()
                    .sort();
                indices
            })
            .collect()
    })
}

pub(super) fn collect_kmer_batch_counts(
    pos_ids: &[u64],
    ctx: &KmerLoadContext<'_>,
    pool: &ThreadPool,
) -> Vec<Vec<usize>> {
    let chunk_len =
        (pos_ids.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    pool.install(|| {
        pos_ids
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut indices = Vec::with_capacity(chunk.len());
                for &pos_id in chunk {
                    let km_id = (pos_id >> 32) as usize;
                    let km_pos = (pos_id & 0x7fff_ffff) as usize;
                    if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                        && let Some(idx) = ctx.mphf.lookup(&bytes)
                    {
                        indices.push(idx as usize);
                    }
                }
                indices
                    .radix_sort_builder()
                    .with_parallel(false)
                    .with_single_threaded_tuner()
                    .sort();
                indices
            })
            .collect()
    })
}

pub(super) fn merge_count_batches(minz_counts: &mut [usize], grouped: Vec<Vec<usize>>) {
    for indices in grouped {
        let mut iter = indices.into_iter();
        let Some(mut idx) = iter.next() else {
            continue;
        };
        let mut count = 1usize;
        for next_idx in iter {
            if next_idx == idx {
                count += 1;
            } else {
                minz_counts[idx] += count;
                idx = next_idx;
                count = 1;
            }
        }
        minz_counts[idx] += count;
    }
}

pub(super) fn build_kmer_pos_index(index: &BifrostIndex) -> KmerPosIndex {
    let mut map: KmerPosIndex = HashMap::new();
    let k = index.k;
    for (unitig_id, seq) in index.unitigs.iter().enumerate() {
        if seq.len() < k {
            continue;
        }
        for pos in 0..=seq.len() - k {
            if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + k]) {
                map.entry(fwd).or_default().push((unitig_id, pos, false));
                if rev != fwd {
                    map.entry(rev).or_default().push((unitig_id, pos, true));
                }
            }
        }
    }
    let base_km = index.unitig_count();
    for (km_id, seq) in index.km_unitigs.iter().enumerate() {
        if seq.len() != k {
            continue;
        }
        if let Some((fwd, rev)) = encode_kmer_pair(seq) {
            let unitig_id = base_km + km_id;
            map.entry(fwd).or_default().push((unitig_id, 0, false));
            if rev != fwd {
                map.entry(rev).or_default().push((unitig_id, 0, true));
            }
        }
    }
    let base_h = index.unitig_count() + index.km_unitig_count();
    for (hk_id, seq) in index.h_kmers.iter().enumerate() {
        if seq.len() != k {
            continue;
        }
        if let Some((fwd, rev)) = encode_kmer_pair(seq) {
            let unitig_id = base_h + hk_id;
            map.entry(fwd).or_default().push((unitig_id, 0, false));
            if rev != fwd {
                map.entry(rev).or_default().push((unitig_id, 0, true));
            }
        }
    }
    map
}
