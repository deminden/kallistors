//! Parallel reconstruction of large-reference cached positions in place.
use super::*;

pub(super) fn populate<K: BucketCountKey, const BATCH: usize, const VIEWS: bool>(
    minz_counts: &mut [usize],
    bitmap: &BitmapPositions,
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: &ThreadPool,
    cache: &mut Vec<u32>,
    fingerprints: &[AtomicU16],
) {
    populate_sized::<K, BATCH, VIEWS, LOAD_BATCH_POSITIONS, LOAD_MIN_PAR_CHUNK>(
        minz_counts,
        bitmap,
        encoded_unitigs,
        mphf,
        pool,
        cache,
        fingerprints,
    );
}

fn populate_sized<
    K: BucketCountKey,
    const BATCH: usize,
    const VIEWS: bool,
    const POSITION_BATCH: usize,
    const MIN_CHUNK: usize,
>(
    minz_counts: &mut [usize],
    bitmap: &BitmapPositions,
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: &ThreadPool,
    cache: &mut Vec<u32>,
    fingerprints: &[AtomicU16],
) {
    assert!(POSITION_BATCH.is_power_of_two() && MIN_CHUNK > 0);
    if bitmap.len() == 0 || encoded_unitigs.unitig_bases.metadata.is_empty() {
        return;
    }
    let old_len = cache.len();
    cache.extend(bitmap.iter());
    let positions = &mut cache[old_len..];
    let chunk_len = (POSITION_BATCH / pool.current_num_threads().saturating_mul(4).max(1))
        .max(MIN_CHUNK)
        .next_power_of_two()
        .min(POSITION_BATCH);
    // Only chunk boundaries are resolved on the caller. Every individual
    // position is reconstructed by its worker, before that cache slot changes.
    let mut state = UnitigBlockState {
        unitig_id: 0,
        total_unitig_len: 0,
        current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
    };
    let mut states = Vec::with_capacity(positions.len().div_ceil(chunk_len));
    for chunk in positions.chunks(chunk_len) {
        advance_unitig_state(&mut state, encoded_unitigs, u64::from(chunk[0]));
        states.push(state);
    }
    for (batch_number, batch) in positions.chunks_mut(POSITION_BATCH).enumerate() {
        let first = batch_number * (POSITION_BATCH / chunk_len);
        let initial = &states[first..first + batch.len().div_ceil(chunk_len)];
        let chunks: Vec<Vec<K>> = pool.install(|| {
            batch
                .par_chunks_mut(chunk_len)
                .zip(initial.par_iter())
                .map(|(entries, &initial)| {
                    let mut state = initial;

                    let mut probe = CachedUnitigMinimizers::new(encoded_unitigs);
                    let mut indices = Vec::with_capacity(entries.len());
                    for destinations in entries.chunks_mut(BATCH) {
                        let mut keys = MphfBatchKeys::<BATCH>::new();
                        for (lane, &raw) in destinations.iter().enumerate() {
                            let position = u64::from(raw);
                            advance_unitig_state(&mut state, encoded_unitigs, position);
                            let relative = (position - state.total_unitig_len) as usize;
                            keys.set(
                                lane,
                                if VIEWS {
                                    probe.at(state.unitig_id, relative)
                                } else {
                                    encoded_unitigs.unitig_minimizer_at(state.unitig_id, relative)
                                },
                            );
                        }
                        let buckets = mphf.lookup_batch(&keys);
                        for (destination, (bucket, key)) in destinations
                            .iter_mut()
                            .zip(buckets.into_iter().zip(keys.iter()))
                        {
                            *destination = bucket.map_or(u32::MAX, |bucket| bucket as u32);
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
        merge_count_chunks(minz_counts, fingerprints, &chunks, pool);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check<K: BucketCountKey, const BATCH: usize>(
        mphf: &BooPhf,
        sequences: &[Vec<u8>],
        g: usize,
        layout: usize,
        threads: usize,
    ) {
        let mut encoded = super::super::precompute_encoded_unitigs(sequences, &[], g, None);
        if layout == 1 {
            let SequenceMetadata::Compact(records) = &encoded.unitig_bases.metadata else {
                panic!("small fixture");
            };
            encoded.unitig_bases.metadata = SequenceMetadata::Wide(
                records
                    .iter()
                    .map(|&word| (word as u32 as usize, (word >> 32) as usize))
                    .collect(),
            );
        } else if layout == 2 {
            encoded.unitig_bases.inline_short_words();
        }
        let mut selected = Vec::new();
        let mut expected_cache = vec![u32::MAX, 17];
        let mut expected_counts = vec![0usize; mphf.size() as usize];
        let fresh = || {
            if K::FINGERPRINT {
                (0..mphf.size())
                    .map(|i| {
                        AtomicU16::new(if i == 1 {
                            crate::pseudoalign::fingerprints::DISABLED
                        } else {
                            0
                        })
                    })
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            }
        };
        let expected_fp = fresh();
        let actual_fp = fresh();
        let mut global = 0u32;
        // Independent ASCII windows and diagnostic full-width MPHF rank oracle.
        // Neither packed metadata nor the worker's position reconstruction is
        // used to construct the expected stream or count/fingerprint state.
        for sequence in sequences {
            if sequence.len() < g {
                continue;
            }
            for query in sequence.windows(g) {
                if global % 5 != 2 {
                    selected.push(global);
                    let key = crate::index::bifrost::encode_minimizer_rep(query).unwrap();
                    let bucket = mphf.debug_lookup(&key).0;
                    expected_cache.push(bucket.map_or(u32::MAX, |value| value as u32));
                    if let Some(bucket) = bucket {
                        expected_counts[bucket as usize] += 1;
                        if K::FINGERPRINT {
                            crate::pseudoalign::fingerprints::observe(
                                &expected_fp[bucket as usize],
                                u64::from_le_bytes(key),
                            );
                        }
                    }
                }
                global += 1;
            }
        }
        assert!(
            selected.len() > 512,
            "fixture must span multiple position batches"
        );
        let bitmap = BitmapPositions::new(selected, true);
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut actual_counts = vec![0usize; mphf.size() as usize];
        let mut actual_cache = vec![u32::MAX, 17];
        populate_sized::<K, BATCH, true, 512, 4>(
            &mut actual_counts,
            &bitmap,
            &encoded,
            mphf,
            &pool,
            &mut actual_cache,
            &actual_fp,
        );
        assert_eq!(
            actual_cache, expected_cache,
            "g={g} layout={layout} threads={threads}"
        );
        assert_eq!(actual_counts, expected_counts);
        for (actual, expected) in actual_fp.iter().zip(&expected_fp) {
            assert_eq!(
                actual.load(Ordering::Relaxed),
                expected.load(Ordering::Relaxed)
            );
        }
    }

    #[test]
    fn inplace_cached_counts_preserve_ascii_order_chunk_boundaries_and_fingerprints() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("transcripts.fa");
        let mut random = 137u64;
        let mut sequence = |len| {
            (0..len)
                .map(|_| {
                    random ^= random << 13;
                    random ^= random >> 7;
                    random ^= random << 17;
                    b"ACGT"[(random & 3) as usize]
                })
                .collect::<Vec<_>>()
        };
        let source = sequence(5000);
        std::fs::write(
            &fasta,
            format!(">t\n{}\n", std::str::from_utf8(&source).unwrap()),
        )
        .unwrap();
        let index = directory.path().join("index.idx");
        crate::index::build_index(&index, &[fasta], crate::index::IndexBuildOptions::default())
            .unwrap();
        let bytes = std::fs::read(index).unwrap();
        let graph = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
        let mphf = BooPhf::load(&mut std::io::Cursor::new(&bytes[16 + graph + 8..])).unwrap();
        let sequences: Vec<_> = [0, 22, 31, 65, 450, 1100, 23, 31]
            .into_iter()
            .map(sequence)
            .collect();
        for (g, layout, threads) in [(23, 0, 2), (23, 1, 3), (23, 2, 8), (1, 0, 3), (32, 1, 2)] {
            check::<u32, 8>(&mphf, &sequences, g, layout, threads);
            check::<TaggedBucket, 64>(&mphf, &sequences, g, layout, threads);
        }
    }
}
