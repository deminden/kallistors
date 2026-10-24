use super::*;

#[test]
fn worker_budget_bounds_and_restores_exact_match_cache() {
    let length = || MPHF_LOOKUP_CACHE.with(|cache| cache.borrow().matches.entries.len());
    assert_eq!(length(), MATCH_CACHE_DEFAULT_ENTRIES);
    {
        let _outer = MatchCacheBudgetScope::enter(1 << 10).unwrap();
        assert_eq!(length(), 1 << 10);
        write_fast_match_cache(42, 7, None);
        assert!(matches!(read_fast_match_cache(42, 7), Some(None)));
        {
            let _inner = MatchCacheBudgetScope::enter(1 << 8).unwrap();
            assert_eq!(length(), 1 << 8);
            assert!(read_fast_match_cache(42, 7).is_none());
            reset_thread_local_caches();
            assert_eq!(length(), 1 << 8);
        }
        assert_eq!(length(), 1 << 10);
        reset_thread_local_caches();
        assert_eq!(length(), 1 << 10);
    }
    assert_eq!(length(), MATCH_CACHE_DEFAULT_ENTRIES);
    let result = std::panic::catch_unwind(|| {
        let _scope = MatchCacheBudgetScope::enter(1 << 10).unwrap();
        panic!("scope unwind oracle");
    });
    assert!(result.is_err());
    assert_eq!(length(), MATCH_CACHE_DEFAULT_ENTRIES);
    assert_eq!(
        MATCH_CACHE_ENTRIES.with(Cell::get),
        MATCH_CACHE_DEFAULT_ENTRIES
    );
}

#[test]
fn compact_match_entries_keep_all_bits_and_reject_wide_coordinates() {
    assert_eq!(std::mem::size_of::<FastKmerMatchCacheEntry>(), 32);
    let mut cache = FastKmerMatchCache::new();
    for bits in 0..8u8 {
        let value = FastKmerMatch {
            unitig_id: u32::MAX as usize,
            start: 71,
            block_idx: 19,
            used_revcomp: bits & 1 != 0,
            forward_strand: bits & 2 != 0,
            matched_relaxed: bits & 4 != 0,
        };
        cache.put(123, 7, Some(value));
        let got = cache.get(123, 7).unwrap().unwrap();
        assert_eq!(
            (got.unitig_id, got.start, got.block_idx),
            (value.unitig_id, value.start, value.block_idx)
        );
        assert_eq!(
            (got.used_revcomp, got.forward_strand, got.matched_relaxed),
            (
                value.used_revcomp,
                value.forward_strand,
                value.matched_relaxed
            )
        );
        assert!(cache.get(123, 6).is_none());
    }
    let collision = (124..u64::MAX)
        .find(|&key| cache.slot(key) == cache.slot(123))
        .unwrap();
    cache.put(collision, 255, None);
    assert!(cache.get(123, 7).unwrap().is_some());
    assert!(matches!(cache.get(collision, 255), Some(None)));
    assert!(cache.get(collision, 15).is_none());
    cache.put(123, 7, None);
    assert!(matches!(cache.get(123, 7), Some(None)));
    if let Some(wide) = (u32::MAX as usize).checked_add(1) {
        for which in 0..3 {
            let mut value = FastKmerMatch::default();
            match which {
                0 => value.unitig_id = wide,
                1 => value.start = wide,
                _ => value.block_idx = wide,
            }
            cache.put(123, 7, Some(value));
            assert!(cache.get(123, 7).is_none());
        }
    }
}

#[test]
fn cache_entries_never_cross_index_or_encoding_context() {
    reset_thread_local_caches();
    ensure_cache_context_key((100, 31, 23));
    write_fast_match_cache(
        42,
        15,
        Some(FastKmerMatch {
            unitig_id: 7,
            ..FastKmerMatch::default()
        }),
    );
    MPHF_LOOKUP_CACHE.with(|cache| cache.borrow_mut().minimizers.put(42, &[([0; 8], 3)]));
    BLOCK_LOOKUP_CACHE.with(|cache| cache.borrow_mut().put(7, 3, 0, 0, 10));
    ensure_cache_context_key((100, 31, 23));
    assert_eq!(
        MPHF_LOOKUP_CACHE.with(|cache| cache.borrow().minimizers.get(42)),
        Some(vec![([0; 8], 3)])
    );
    assert_eq!(
        read_fast_match_cache(42, 15).flatten().unwrap().unitig_id,
        7
    );
    ensure_cache_context_key((101, 31, 23));
    assert!(read_fast_match_cache(42, 15).is_none());
    assert!(MPHF_LOOKUP_CACHE.with(|cache| cache.borrow().minimizers.get(42).is_none()));
    assert!(BLOCK_LOOKUP_CACHE.with(|cache| cache.borrow().get(7, 3).is_none()));
    write_fast_match_cache(42, 15, None);
    ensure_cache_context_key((101, 31, 21));
    assert!(read_fast_match_cache(42, 15).is_none());
}

#[test]
fn candidate_cache_never_truncates_a_tie_list() {
    assert_eq!(std::mem::size_of::<MinimizerCandidateCacheEntry>(), 16);
    let mut cache = MinimizerCandidateCache::new();
    assert!(cache.entries.is_empty());
    assert!(cache.get(0).is_none());
    for word in [0, 0x8000_0000_0000_0000u64, 0xffff_ffff_ffff_ff00] {
        for position in [0, 1, 31, 127] {
            cache.put(40, &[(word.to_le_bytes(), position)]);
            assert_eq!(cache.get(40), Some(vec![(word.to_le_bytes(), position)]));
        }
    }
    for (word, position) in [(0u64, 128), (1, 31), (255, 127), (0, usize::MAX)] {
        cache.put(40, &[(word.to_le_bytes(), position)]);
        assert!(cache.get(40).is_none());
    }
    cache.put(41, &[]);
    assert_eq!(cache.get(41), Some(vec![]));
    cache.put(42, &[([1; 8], usize::MAX); 8]);
    cache.put(42, &[([1; 8], 3); 9]);
    assert!(cache.get(42).is_none());
    cache.put(43, &[([1; 8], usize::MAX); 8]);
    assert_eq!(cache.get(43), Some(vec![([1; 8], usize::MAX); 8]));
    let single = (0x1234_5600u64.to_le_bytes(), 5);
    cache.put(43, &[single]);
    assert_eq!(cache.get(43), Some(vec![single]));
    let collision = (44..u64::MAX)
        .find(|&key| cache.slot(key) == cache.slot(43))
        .unwrap();
    cache.put(collision, &[(0x6789_ab00u64.to_le_bytes(), 7)]);
    assert!(cache.get(43).is_none());
    assert_eq!(
        cache.get(collision),
        Some(vec![(0x6789_ab00u64.to_le_bytes(), 7)])
    );
}
