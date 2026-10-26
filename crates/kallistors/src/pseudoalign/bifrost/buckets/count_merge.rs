//! Exclusive range ownership for exact minimizer count construction.
use super::*;

pub(super) fn merge_count_chunks<K: BucketCountKey>(
    counts: &mut [usize],
    fingerprints: &[AtomicU16],
    chunks: &[Vec<K>],
    pool: &ThreadPool,
) {
    if counts.len() >= (1 << 28)
        && pool.current_num_threads() > 1
        && chunks.iter().map(Vec::len).sum::<usize>() >= LOAD_MIN_PAR_CHUNK
    {
        merge_count_shards(counts, fingerprints, chunks, pool);
    } else {
        for chunk in chunks {
            merge_count_range(counts, fingerprints, chunk, 0);
        }
    }
}

fn merge_count_shards<K: BucketCountKey>(
    counts: &mut [usize],
    fingerprints: &[AtomicU16],
    chunks: &[Vec<K>],
    pool: &ThreadPool,
) {
    let width = counts.len().div_ceil(pool.current_num_threads()).max(1);
    pool.install(|| {
        counts
            .par_chunks_mut(width)
            .enumerate()
            .for_each(|(owner, counts)| {
                let first = owner * width;
                let end = first + counts.len();
                // Every source chunk is already sorted by bucket. An owner's
                // slices are disjoint, and source-chunk order remains unchanged.
                for chunk in chunks {
                    let lo = chunk.partition_point(|record| (record.bucket() as usize) < first);
                    let hi =
                        lo + chunk[lo..].partition_point(|record| (record.bucket() as usize) < end);
                    merge_count_range(counts, fingerprints, &chunk[lo..hi], first);
                }
            });
    });
}

fn merge_count_range<K: BucketCountKey>(
    minz_counts: &mut [usize],
    fingerprints: &[AtomicU16],
    chunk: &[K],
    first_bucket: usize,
) {
    let mut at = 0;
    while at < chunk.len() {
        let bucket = chunk[at].bucket();
        let first = chunk[at].fingerprint();
        let mut fingerprint_difference = 0u16;
        let mut end = at + 1;
        while end < chunk.len() && chunk[end].bucket() == bucket {
            if K::FINGERPRINT {
                fingerprint_difference |= chunk[end].fingerprint() ^ first;
            }
            end += 1;
        }
        if K::FINGERPRINT {
            // Bucket order is sufficient. Check every member while finding
            // the group boundary, including conflicts between equal endpoints.
            let value = if fingerprint_difference == 0 {
                first
            } else {
                crate::pseudoalign::fingerprints::DISABLED
            };
            crate::pseudoalign::fingerprints::observe_serial(&fingerprints[bucket as usize], value);
        }
        minz_counts[bucket as usize - first_bucket] += end - at;
        at = end;
    }
}

#[cfg(test)]
mod owned_count_tests {
    use super::*;

    fn check<K: BucketCountKey>(buckets: usize, threads: usize) {
        let mut random = 173u64;
        let mut expected = (0..buckets).map(|i| i * 3 + 7).collect::<Vec<_>>();
        let mut actual = expected.clone();
        let fresh = || {
            if K::FINGERPRINT {
                (0..buckets)
                    .map(|i| {
                        AtomicU16::new(if i == 5 {
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
        let mut chunks = Vec::new();
        for size in [0, 1, 512, 17, 0, 129, 1025, 3] {
            let mut records = Vec::new();
            for _ in 0..size {
                random ^= random << 13;
                random ^= random >> 7;
                random ^= random << 17;
                let bucket = random as usize % buckets;
                let word = random.rotate_left(23);
                expected[bucket] += 1;
                if K::FINGERPRINT {
                    crate::pseudoalign::fingerprints::observe(&expected_fp[bucket], word);
                }
                records.push(K::new(bucket as u32, word));
                // Repeated exact records and conflicting fingerprints can
                // cross either source-chunk or owner-range boundaries.
                if random & 7 == 0 {
                    expected[bucket] += 1;
                    records.push(K::new(bucket as u32, word));
                }
            }
            records
                .radix_sort_builder()
                .with_parallel(false)
                .with_single_threaded_tuner()
                .sort();
            chunks.push(records);
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        merge_count_shards(&mut actual, &actual_fp, &chunks, &pool);
        assert_eq!(actual, expected, "buckets={buckets} workers={threads}");
        for (a, b) in actual_fp.iter().zip(&expected_fp) {
            assert_eq!(a.load(Ordering::Relaxed), b.load(Ordering::Relaxed));
        }
    }

    #[test]
    fn exclusive_count_ranges_match_independent_record_updates_and_fingerprint_collisions() {
        for buckets in [1, 17, 113, 257] {
            for threads in [2, 3, 8] {
                check::<u32>(buckets, threads);
                check::<TaggedBucket>(buckets, threads);
            }
        }
        let pool = ThreadPoolBuilder::new().num_threads(3).build().unwrap();
        merge_count_shards::<u32>(&mut [], &[], &[], &pool);
    }
}

#[cfg(test)]
#[test]
fn bucket_only_groups_detect_hidden_fingerprint_conflicts_and_special_words() {
    let first = 17u64;
    let second = 31u64;
    assert_ne!(
        crate::pseudoalign::fingerprints::word_fingerprint(first),
        crate::pseudoalign::fingerprints::word_fingerprint(second)
    );
    for middle in [first, second, 0, u64::MAX] {
        let words = [first, middle, first];
        // This is bucket-sorted but deliberately not fingerprint-sorted. An
        // endpoint-only check would miss the middle conflicting/zero value.
        let records = words.map(|word| TaggedBucket::new(0, word));
        let expected = AtomicU16::new(0);
        for word in words {
            crate::pseudoalign::fingerprints::observe(&expected, word);
        }
        let actual = [AtomicU16::new(0)];
        let mut counts = [11usize];
        merge_count_range(&mut counts, &actual, &records, 0);
        assert_eq!(counts, [14]);
        assert_eq!(
            actual[0].load(Ordering::Relaxed),
            expected.load(Ordering::Relaxed)
        );
    }
}
