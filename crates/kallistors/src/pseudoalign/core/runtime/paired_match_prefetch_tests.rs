use super::*;
#[test]
fn hints_do_not_insert_evict_or_change_exact_match_answers() {
    let mut cache = MphfLookupCache::default();
    let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACG";
    let (key, _) = crate::dna::encode_kmer_pair(sequence).unwrap();
    let value = FastKmerMatch {
        unitig_id: 13,
        start: 23,
        block_idx: 8,
        used_revcomp: true,
        forward_strand: false,
        matched_relaxed: true,
    };
    cache.match_put(key, 15, Some(value));
    cache.match_put(key ^ 1, 14, None);
    let before = cache.matches.entries.as_ptr();
    let length = cache.matches.entries.len();
    for sequence in [
        sequence.as_slice(),
        b"ACGT".as_slice(),
        b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".as_slice(),
        b"acgtacgtacgtacgtacgtacgtacgtacg".as_slice(),
    ] {
        cache.prefetch_initial_match_words(sequence, sequence);
    }
    let got = cache.match_get(key, 15).flatten().unwrap();
    assert_eq!(
        (
            got.unitig_id,
            got.start,
            got.block_idx,
            got.used_revcomp,
            got.forward_strand,
            got.matched_relaxed
        ),
        (13, 23, 8, true, false, true)
    );
    assert!(matches!(cache.match_get(key ^ 1, 14), Some(None)));
    assert!(cache.match_get(key, 14).is_none());
    assert_eq!(before, cache.matches.entries.as_ptr());
    assert_eq!(length, cache.matches.entries.len());
    assert!(!cache.prefetch_match_capacity());
}
