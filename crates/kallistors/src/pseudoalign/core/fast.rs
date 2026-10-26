mod matching;
mod pair_cache;
use pair_cache::PairData;
pub(crate) use pair_cache::ScopedEcPairs;
mod rescue;
mod scanner;

use matching::{all_ascii_uppercase, match_kmer_at_pos_fast_with_codes};
use scanner::ec_for_read_bifrost_fast_with_cache;

use super::super::BorrowedReadEc;
use super::super::debug::{
    DroppedHitTrace, JumpDecisionTrace, MinimizerCandidateTrace, ReadDebugState,
    format_positions_sample,
};
use super::super::ec::{
    KmerPairCursor, encode_kmer_pair, encoded_kmer_at_or_after, intersect_sorted_in_place,
};
use super::super::minimizers::minimizer_tail_for_encoded;
use super::super::{BifrostIndex, Hit, MatchInfo, PseudoalignOptions, ReadEc, Strand};
use super::matcher::{
    jump_distance_for_match, match_kmer_at_pos, match_kmer_direct,
    match_unitig_candidate_encoded_typed, special_unitig_for_kmer,
};
use super::runtime::{
    DISABLE_MATCH_CACHE, DISABLE_MINIMIZER_CACHE, FastKmerMatch, MphfLookupCache,
    block_index_for_position_fast, ec_slice, env_flag, minimizer_bucket_has_overcrowded_marker,
    minimizer_bucket_idx, minimizer_bucket_with_fingerprint,
    minimizer_candidates_cached_into_with_code, minimizer_candidates_single_cell, minimizer_hex,
    onlist_cardinality,
};

// The lazy final-hit fold shares this prefix. Keep an immutable graph slice
// until an actually different constraint requires writable intersection storage.
struct OnlineIntersection<'index, 'scratch> {
    initial: Option<&'index [u32]>,
    owned: &'scratch mut Vec<u32>,
    pairs: Option<&'scratch mut PairData>,
}
impl std::ops::Deref for OnlineIntersection<'_, '_> {
    type Target = [u32];
    fn deref(&self) -> &[u32] {
        self.initial.unwrap_or(self.owned)
    }
}
impl<'index> OnlineIntersection<'index, '_> {
    fn initialize(&mut self, ec: &'index [u32]) {
        self.initial = Some(ec);
    }
    #[inline(always)]
    fn intersect(&mut self, ec: &'index [u32]) {
        let initial = self.initial;
        let key = self.pairs.as_ref().and_then(|_| {
            self.owned
                .is_empty()
                .then_some(initial)
                .flatten()
                .and_then(|a| pair_cache::Key::eligible(a, ec))
        });
        if let (Some(key), Some(initial), Some(pairs)) = (key, initial, self.pairs.as_deref())
            && pairs.apply(key, initial, ec, &mut self.initial, self.owned)
        {
            return;
        }
        self.intersect_original(ec);
        if let (Some(key), Some(initial), Some(pairs)) = (key, initial, self.pairs.as_deref_mut()) {
            pairs.remember(key, initial, ec, self.initial, self.owned);
        }
    }
    #[inline(always)]
    fn intersect_original(&mut self, ec: &'index [u32]) {
        if let Some(initial) = self.initial {
            if std::ptr::eq(initial, ec) || initial == ec {
                return;
            }
            // Nonoverlapping transcript ranges prove rejection before a copy.
            if initial.is_empty()
                || ec.is_empty()
                || initial[0] > ec[ec.len() - 1]
                || ec[0] > initial[initial.len() - 1]
            {
                self.initial = Some(&initial[..0]);
                return;
            }
            let (shorter, longer) = if initial.len() <= ec.len() {
                (initial, ec)
            } else {
                (ec, initial)
            };
            // Many adjacent transcript annotations add/remove alternatives. A
            // subset is already the exact intersection and can remain borrowed.
            if shorter.len() != longer.len() && sorted_subset(shorter, longer) {
                self.initial = Some(shorter);
                return;
            }
            self.owned.extend_from_slice(shorter);
            self.initial = None;
            intersect_sorted_in_place(self.owned, longer);
        } else {
            intersect_sorted_in_place(self.owned, ec);
        }
    }
}

fn sorted_subset(shorter: &[u32], longer: &[u32]) -> bool {
    if shorter.is_empty() {
        return true;
    }
    if longer.is_empty()
        || shorter[0] < longer[0]
        || shorter[shorter.len() - 1] > longer[longer.len() - 1]
    {
        return false;
    }
    let mut cursor = 0;
    // Positive shorter length is established above. This is exactly
    // floor(longer/shorter)>=16, using a constant shift instead of division.
    if shorter.len() <= longer.len() / 16 {
        for &value in shorter {
            if longer.get(cursor) == Some(&value) {
                cursor += 1;
                continue;
            }
            let remaining = &longer[cursor..];
            let mut upper = 1usize;
            while upper < remaining.len() && remaining[upper] < value {
                upper = upper.saturating_mul(2);
            }
            let lower = (upper >> 1).min(remaining.len());
            let upper = upper.min(remaining.len());
            cursor +=
                lower + remaining[lower..upper].partition_point(|&candidate| candidate < value);
            if longer.get(cursor) != Some(&value) {
                return false;
            }
            cursor += 1;
        }
        return true;
    }
    // Near-sized ECs often retain long identical runs around an internal
    // insertion. Consume both equal blocks once; a mismatch resumes the same
    // ordered scalar membership traversal, including duplicate multiplicity.
    // Sparse galloping above and the dense general fallback below are intact.
    if shorter.len() >= 8
        && shorter.len() <= longer.len()
        && longer.len() - shorter.len() < shorter.len()
    {
        let (mut remaining, mut available) = (shorter, longer);
        while !remaining.is_empty() {
            if remaining.len() >= 4 && available.len() >= 4 && remaining[..4] == available[..4] {
                remaining = &remaining[4..];
                available = &available[4..];
                continue;
            }
            let value = remaining[0];
            while available
                .first()
                .is_some_and(|&candidate| candidate < value)
            {
                available = &available[1..];
            }
            if available.first() != Some(&value) {
                return false;
            }
            remaining = &remaining[1..];
            available = &available[1..];
        }
        return true;
    }
    for &value in shorter {
        while cursor < longer.len() && longer[cursor] < value {
            cursor += 1;
        }
        if longer.get(cursor) != Some(&value) {
            return false;
        }
        cursor += 1;
    }
    true
}

#[derive(Default)]
pub(super) struct FastScannerScratch {
    min_candidates: Vec<([u8; 8], usize)>,
    jump_candidates: Vec<([u8; 8], usize)>,
    online_intersection: Vec<u32>,
    pairs: Option<Box<PairData>>,
}

impl FastScannerScratch {
    fn clear(&mut self) {
        self.min_candidates.clear();
        self.jump_candidates.clear();
        self.online_intersection.clear();
    }
}

pub(crate) fn ec_for_read_bifrost_fast(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    super::runtime::with_mphf_lookup_cache(|cache| {
        // Reuse one heap owner instead of moving three Vec headers twice.
        // Default recovery preserves usability after an unwinding scanner.
        let mut scratch = cache.scan_scratch.take().unwrap_or_default();
        scratch.clear();
        let result = if dbg.is_some() {
            ec_for_read_bifrost_fast_with_cache::<true, false, true>(
                index,
                seq,
                strand,
                dbg,
                options,
                cache,
                &mut scratch,
            )
        } else {
            ec_for_read_bifrost_fast_with_cache::<false, false, true>(
                index,
                seq,
                strand,
                dbg,
                options,
                cache,
                &mut scratch,
            )
        };
        cache.scan_scratch = Some(scratch);
        result.map(BorrowedReadEc::into_owned)
    })
}

#[inline(always)]
pub(super) fn ec_for_read_bifrost_fast_counting_borrowed<'index, const PLACEMENT: bool>(
    index: &'index BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<BorrowedReadEc<'index>> {
    super::runtime::with_mphf_lookup_cache(|cache| {
        // Reuse one heap owner instead of moving three Vec headers twice.
        // Default recovery preserves usability after an unwinding scanner.
        let mut scratch = cache.scan_scratch.take().unwrap_or_default();
        scratch.clear();
        let result = ec_for_read_bifrost_fast_with_cache::<false, true, PLACEMENT>(
            index,
            seq,
            strand,
            None,
            options,
            cache,
            &mut scratch,
        );
        cache.scan_scratch = Some(scratch);
        result
    })
}

// Scope existing per-thread state over an ordered sequence batch. The same
// scanner handles every read; reset/diagnostic/alternate paths keep old calls.
pub(crate) fn try_single_counting_batch<'index>(
    index: &'index BifrostIndex,
    batch: &crate::io::PackedSeqBatch,
    strand: Strand,
    options: PseudoalignOptions,
    pairs: &mut ScopedEcPairs<'index>,
    mut count: impl FnMut(Option<BorrowedReadEc<'index>>),
) -> bool {
    if index.k == 0
        || super::super::reset_all_caches_per_read()
        || options.investigation.any_enabled()
        || !super::runtime::fast_path_enabled(index, false, options)
    {
        return false;
    }
    super::runtime::ensure_cache_context(index);
    super::runtime::with_mphf_lookup_cache(|cache| {
        let mut scratch = cache.scan_scratch.take().unwrap_or_default();
        let pair_scope_entered = pairs.enter(index, &mut scratch.pairs);
        for record in batch.records() {
            if record.seq.len() < index.k {
                count(None);
                continue;
            }
            scratch.clear();
            count(ec_for_read_bifrost_fast_with_cache::<false, true, true>(
                index,
                record.seq,
                strand,
                None,
                options,
                cache,
                &mut scratch,
            ));
        }
        if pair_scope_entered {
            pairs.leave(&mut scratch.pairs);
        } else {
            scratch.pairs = None;
        }
        cache.scan_scratch = Some(scratch);
    });
    true
}

// Use the established read scanner with one cache/scratch borrow per ordered
// batch. Filtering mate one before scanning mate two preserves hard rejection
// and the existing cache access order. The caller has exhausted sampling.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_paired_counting_batch<'index>(
    index: &'index BifrostIndex,
    left: &crate::io::PackedSeqBatch,
    right: &crate::io::PackedSeqBatch,
    range: std::ops::Range<usize>,
    strand: Strand,
    options: PseudoalignOptions,
    pairs: &mut ScopedEcPairs<'index>,
    count: impl FnMut(Option<BorrowedReadEc<'index>>, Option<BorrowedReadEc<'index>>),
) -> bool {
    if index.k == 0
        || super::super::reset_all_caches_per_read()
        || options.investigation.any_enabled()
        || !super::runtime::fast_path_enabled(index, false, options)
    {
        return false;
    }
    super::runtime::ensure_cache_context(index);
    if options.strand_specific.is_some() {
        paired_counting_batch_with_cache::<true>(
            index, left, right, range, strand, options, pairs, count,
        );
    } else {
        paired_counting_batch_with_cache::<false>(
            index, left, right, range, strand, options, pairs, count,
        );
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn paired_counting_batch_with_cache<'index, const PLACEMENT: bool>(
    index: &'index BifrostIndex,
    left: &crate::io::PackedSeqBatch,
    right: &crate::io::PackedSeqBatch,
    range: std::ops::Range<usize>,
    strand: Strand,
    options: PseudoalignOptions,
    pairs: &mut ScopedEcPairs<'index>,
    mut count: impl FnMut(Option<BorrowedReadEc<'index>>, Option<BorrowedReadEc<'index>>),
) {
    super::runtime::with_mphf_lookup_cache(|cache| {
        let mut scratch = cache.scan_scratch.take().unwrap_or_default();
        let pair_scope_entered = pairs.enter(index, &mut scratch.pairs);
        let prefetch = index.k == 31
            && index.g == 23
            && cache.prefetch_match_capacity()
            && !super::runtime::env_flag(
                &super::runtime::DISABLE_MATCH_CACHE,
                "KALLISTORS_DISABLE_MATCH_CACHE",
            );
        for idx in range {
            // SAFETY: the caller validated both batch lengths and this range.
            let left_sequence = unsafe { left.seq_unchecked(idx) };
            let right_sequence = unsafe { right.seq_unchecked(idx) };
            // Pure encoding and read-only hardware hints. Actual matcher gets,
            // puts, filters, rejected mates and callback order remain unchanged.
            if prefetch {
                cache.prefetch_initial_match_words(left_sequence, right_sequence);
            }
            let mut scan = |seq: &[u8]| {
                if seq.len() < index.k {
                    return None;
                }
                scratch.clear();
                ec_for_read_bifrost_fast_with_cache::<false, true, PLACEMENT>(
                    index,
                    seq,
                    strand,
                    None,
                    options,
                    cache,
                    &mut scratch,
                )
            };
            let a = scan(left_sequence).map(|read| {
                if let Some(mode) = options.strand_specific {
                    read.filter_strand(index, mode, true)
                } else {
                    read
                }
            });
            if a.as_ref().is_some_and(|read| read.read.hard_reject_pair) {
                count(a, None);
                continue;
            }
            let b = scan(right_sequence).map(|read| {
                if let Some(mode) = options.strand_specific {
                    read.filter_strand(index, mode, false)
                } else {
                    read
                }
            });
            count(a, b);
        }
        if pair_scope_entered {
            pairs.leave(&mut scratch.pairs);
        } else {
            scratch.pairs = None;
        }
        cache.scan_scratch = Some(scratch);
    });
}

#[cfg(test)]
mod uppercase_tests {
    use super::all_ascii_uppercase;
    #[test]
    fn uppercase_vectors_match_scalar_for_all_bytes_and_record_tails() {
        for len in [
            0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 95, 100, 127, 128, 129, 257,
        ] {
            let mut input = vec![b'N'; len];
            assert!(all_ascii_uppercase(&input));
            for position in 0..len {
                for value in 0..=u8::MAX {
                    input[position] = value;
                    assert_eq!(
                        all_ascii_uppercase(&input),
                        input.iter().all(u8::is_ascii_uppercase),
                        "{len} {position} {value}"
                    );
                }
                input[position] = b'N';
            }
        }
    }
}

#[cfg(test)]
mod borrowed_prefix_tests {
    use super::*;
    #[test]
    fn borrowed_prefix_materializes_only_when_a_constraint_changes() {
        let original = [1, 1, 3, 7];
        let same = [1, 1, 3, 7];
        let mut scratch = Vec::new();
        let mut online = OnlineIntersection {
            initial: None,
            owned: &mut scratch,
            pairs: None,
        };
        online.initialize(&original);
        online.intersect(&same);
        assert!(online.initial.is_some());
        assert_eq!(&*online, &original);
        online.intersect(&[1, 3, 3, 8]);
        assert!(online.initial.is_none());
        assert_eq!(&*online, &[1, 3]);
        online.intersect(&[3]);
        assert_eq!(&*online, &[3]);
        online.intersect(&[]);
        assert!(online.is_empty());
        let mut multisets: Vec<Vec<u32>> = (0..81)
            .map(|mut code| {
                let mut values = Vec::new();
                for value in 0..4 {
                    for _ in 0..code % 3 {
                        values.push(value);
                    }
                    code /= 3;
                }
                values
            })
            .collect();
        multisets.extend([(0..256).collect(), vec![15, 15], vec![255], vec![254, 256]]);
        // Exercise equal blocks before/after internal gaps and scalar tails.
        // The existing independent intersection oracle also checks ownership.
        multisets.push((0..64).collect());
        for removed in [0, 1, 3, 4, 7, 8, 15, 31, 63] {
            multisets.push((0..64).filter(|&value| value != removed).collect());
        }
        let duplicates: Vec<_> = (0..32).flat_map(|value| [value, value]).collect();
        multisets.push(duplicates.clone());
        for removed in [1, 7, 31, 62] {
            let mut values = duplicates.clone();
            values.remove(removed);
            multisets.push(values);
        }
        for left in &multisets {
            for right in &multisets {
                let mut expected = Vec::new();
                super::super::super::ec::intersect_sorted(left, right, &mut expected);
                let mut scratch = Vec::new();
                let mut intersection = OnlineIntersection {
                    initial: Some(left),
                    owned: &mut scratch,
                    pairs: None,
                };
                intersection.intersect(right);
                assert_eq!(&*intersection, &expected);
                if expected == *left || expected == *right {
                    assert!(intersection.initial.is_some());
                }
            }
        }
    }
}

#[cfg(test)]
mod skewed_subset_tests {
    use super::sorted_subset;

    #[test]
    fn skewed_subset_matches_independent_multiset_oracle() {
        for length in [0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 255, 256, 257] {
            let longer: Vec<u32> = (0..length).map(|i| ((i / 2) * 3) as u32).collect();
            let mut cases = vec![vec![], vec![0], vec![1], vec![u32::MAX], vec![0, 0, 0]];
            for count in 0..=8 {
                let count = count.min(length);
                cases.push(longer[..count].to_vec());
                cases.push(longer[length - count..].to_vec());
            }
            for stride in [15, 16, 17, 32] {
                cases.push(longer.iter().step_by(stride).copied().take(8).collect());
            }
            let missing: Vec<_> = cases
                .iter()
                .filter(|case| case.len() >= 3)
                .map(|case| {
                    let mut case = case.clone();
                    let middle = case.len() / 2;
                    case[middle] += 1;
                    case.sort_unstable();
                    case
                })
                .collect();
            cases.extend(missing);
            for shorter in cases {
                // Independent multiplicity oracle, with no cursor/gallop or
                // borrowed subset helper shared with the implementation.
                let mut remaining = longer.clone();
                let expected = shorter.iter().all(|value| {
                    if let Some(at) = remaining.iter().position(|item| item == value) {
                        remaining.remove(at);
                        true
                    } else {
                        false
                    }
                });
                assert_eq!(
                    sorted_subset(&shorter, &longer),
                    expected,
                    "shorter={shorter:?} longer={longer:?}"
                );
            }
        }
    }
}
