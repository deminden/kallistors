mod matching;
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
        if let Some(initial) = self.initial {
            if initial == ec {
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
    if longer.len() / shorter.len() >= 16 {
        for &value in shorter {
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
        for left in &multisets {
            for right in &multisets {
                let mut expected = Vec::new();
                super::super::super::ec::intersect_sorted(left, right, &mut expected);
                let mut scratch = Vec::new();
                let mut intersection = OnlineIntersection {
                    initial: Some(left),
                    owned: &mut scratch,
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
