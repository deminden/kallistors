//! The final all-hit EC fold for counting after placement sampling has finished.

use super::super::{BifrostIndex, Hit};
use super::runtime::ec_slice;

#[derive(Default)]
pub(super) struct FinalHitFold {
    pub(super) ec: Vec<u32>,
    pub(super) rejected: bool,
    found_nonempty: bool,
    last: Option<(usize, usize)>,
}

impl FinalHitFold {
    pub(super) fn from_online(ec: &[u32], last: Option<(usize, usize)>) -> Self {
        Self {
            ec: ec.to_vec(),
            rejected: false,
            found_nonempty: !ec.is_empty(),
            last,
        }
    }

    #[inline(always)]
    pub(super) fn push(&mut self, index: &BifrostIndex, hit: Hit) {
        // Empty annotations contribute no constraint; an empty intersection of
        // nonempty annotations is terminal and must never restart the fold.
        if self.rejected {
            return;
        }
        let key = (hit.unitig_id, hit.block_idx);
        let ec = ec_slice(index, hit.unitig_id, hit.block_idx);
        if !self.found_nonempty {
            if !ec.is_empty() {
                self.ec.extend_from_slice(ec);
                self.last = Some(key);
                self.found_nonempty = true;
            }
            return;
        }
        if self.last == Some(key) || ec.is_empty() {
            return;
        }
        super::super::ec::intersect_sorted_in_place(&mut self.ec, ec);
        self.rejected = self.ec.is_empty();
        self.last = Some(key);
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::{PseudoalignOptions, ReadEc, Strand};
    use super::super::read::{ec_for_read_bifrost, ec_for_read_bifrost_counting};

    #[test]
    fn counting_fold_matches_complete_hit_collection_after_sampling() {
        let dir = tempfile::TempDir::new().unwrap();
        let fasta = dir.path().join("transcripts.fa");
        let path = dir.path().join("index.idx");
        let mut state = 19u64;
        let mut make = || -> Vec<u8> {
            (0..400)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect()
        };
        let a = make();
        let mut b = make();
        b[..200].copy_from_slice(&a[..200]);
        let c = make();
        std::fs::write(
            &fasta,
            format!(
                ">a\n{}\n>b\n{}\n>c\n{}\n",
                String::from_utf8_lossy(&a),
                String::from_utf8_lossy(&b),
                String::from_utf8_lossy(&c)
            ),
        )
        .unwrap();
        crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default())
            .unwrap();
        for mut index in [
            crate::pseudoalign::build_bifrost_index_threaded(&path, 2).unwrap(),
            crate::pseudoalign::build_bifrost_index_positional_packed_threaded(&path, 2, false)
                .unwrap(),
        ] {
            for mask in [
                Some(vec![true; 3]),
                Some(vec![false, true, true]),
                Some(vec![false; 3]),
                None,
            ] {
                index.onlist = mask;
                let borrowed = super::super::read::ec_for_read_bifrost_counting_borrowed(
                    &index,
                    &a[31..151],
                    Strand::Unstranded,
                    PseudoalignOptions::default(),
                );
                let full = ec_for_read_bifrost(
                    &index,
                    &a[31..151],
                    Strand::Unstranded,
                    None,
                    PseudoalignOptions::default(),
                );
                assert_eq!(
                    borrowed.as_ref().map(|v| v.ec()),
                    full.as_ref().map(|v| v.ec.as_slice())
                );
                if index.onlist.as_ref().is_none_or(|m| m.iter().all(|&x| x)) {
                    assert!(borrowed.as_ref().unwrap().borrowed.is_some());
                }
            }
            index.onlist = Some(vec![true; 3]);
            let mut queries = Vec::new();
            for seq in [&a, &b, &c] {
                for window in seq.windows(120).step_by(37) {
                    queries.push(window.to_vec());
                    let mut altered = window.to_vec();
                    altered[45] = if altered[45] == b'A' { b'C' } else { b'A' };
                    queries.push(altered);
                    let mut ambiguous = window.to_vec();
                    ambiguous[50] = b'N';
                    queries.push(ambiguous);
                }
            }
            queries.push([&a[..70], &c[..70]].concat());
            queries.push(a[..90].to_ascii_lowercase());
            queries.push(a[..20].to_vec());
            let key = |value: Option<ReadEc>| {
                value.map(|value| (value.ec, value.hard_reject_pair, value.had_offlist))
            };
            for read in queries {
                for strand in [Strand::Unstranded, Strand::Forward, Strand::Reverse] {
                    super::super::runtime::reset_thread_local_caches();
                    let full = ec_for_read_bifrost(
                        &index,
                        &read,
                        strand,
                        None,
                        PseudoalignOptions::default(),
                    );
                    super::super::runtime::reset_thread_local_caches();
                    let counted = ec_for_read_bifrost_counting(
                        &index,
                        &read,
                        strand,
                        PseudoalignOptions::default(),
                    );
                    let borrowed = super::super::read::ec_for_read_bifrost_counting_borrowed(
                        &index,
                        &read,
                        strand,
                        PseudoalignOptions::default(),
                    )
                    .map(super::super::super::BorrowedReadEc::into_owned);
                    let placement_key = |value: &ReadEc| {
                        (
                            value.ec.clone(),
                            value.hard_reject_pair,
                            value.had_offlist,
                            value
                                .best_match
                                .map(|m| (m.unitig_id, m.unitig_pos, m.read_pos, m.used_revcomp)),
                            value
                                .first_hit
                                .as_ref()
                                .map(|h| (h.unitig_id, h.block_idx, h.read_pos, h.used_revcomp)),
                        )
                    };
                    let placement = super::super::read::ec_for_read_bifrost_single(
                        &index,
                        &read,
                        strand,
                        PseudoalignOptions::default(),
                    );
                    // Rejected ECs can stop before producing unused placement metadata.
                    if full.as_ref().is_some_and(|v| !v.ec.is_empty()) {
                        assert_eq!(
                            placement.as_ref().map(placement_key),
                            full.as_ref().map(placement_key)
                        );
                    } else {
                        assert_eq!(
                            placement.as_ref().map(|v| &v.ec),
                            full.as_ref().map(|v| &v.ec)
                        );
                    }
                    if let Some(value) = full.as_ref()
                        && let Some(best) = value.best_match
                    {
                        for length in [0, 100, 200, 800] {
                            let expected = super::super::filters::filter_ec_by_fragment(
                                &index, &value.ec, best, length,
                            );
                            for ec in [
                                std::borrow::Cow::Borrowed(value.ec.as_slice()),
                                std::borrow::Cow::Owned(value.ec.clone()),
                            ] {
                                let filtered = super::super::filters::filter_ec_by_fragment_cow(
                                    &index, ec, best, length,
                                );
                                assert_eq!(filtered.as_ref(), expected.as_slice());
                            }
                            if let Some(hit) = value.first_hit.as_ref() {
                                for mode in [
                                    crate::pseudoalign::StrandSpecific::FR,
                                    crate::pseudoalign::StrandSpecific::RF,
                                ] {
                                    for ec in [
                                        std::borrow::Cow::Borrowed(value.ec.as_slice()),
                                        std::borrow::Cow::Owned(value.ec.clone()),
                                        std::borrow::Cow::Borrowed(
                                            index.flat_ec.ec(hit.unitig_id, hit.block_idx),
                                        ),
                                    ] {
                                        let fragment =
                                            super::super::filters::filter_ec_by_fragment_cow(
                                                &index,
                                                ec.clone(),
                                                best,
                                                length,
                                            );
                                        let oracle = super::super::filters::filter_ec_by_strand_cow(
                                            &index, fragment, hit, mode, true,
                                        );
                                        let combined = super::super::filters::filter_ec_by_fragment_and_strand_cow(&index, ec, best, length, hit, mode);
                                        assert_eq!(combined.as_ref(), oracle.as_ref());
                                    }
                                }
                            }
                            let borrowed = super::super::filters::filter_ec_by_fragment_cow(
                                &index,
                                std::borrow::Cow::Borrowed(&value.ec),
                                best,
                                length,
                            );
                            assert_eq!(
                                matches!(borrowed, std::borrow::Cow::Borrowed(_)),
                                expected == value.ec
                            );
                        }
                    }
                    if let Some(value) = full.as_ref()
                        && let Some(hit) = value.first_hit.as_ref()
                    {
                        let block = &index.flat_ec.owned_ec_blocks[hit.unitig_id][hit.block_idx];
                        let mut missing = value.ec.clone();
                        missing.push(u32::MAX);
                        for candidates in [block.ec.as_ref(), value.ec.as_slice(), &missing] {
                            for mode in [
                                crate::pseudoalign::StrandSpecific::FR,
                                crate::pseudoalign::StrandSpecific::RF,
                            ] {
                                for first in [false, true] {
                                    let oracle = ReadEc {
                                        ec: candidates.to_vec(),
                                        best_match: value.best_match,
                                        first_hit: Some(hit.clone()),
                                        hits: Vec::new(),
                                        had_offlist: false,
                                        shade_union: Vec::new(),
                                        hard_reject_pair: false,
                                    };
                                    let expected = super::super::filters::apply_strand_filter(
                                        &index, &oracle, mode, first, false, &mut None, b"",
                                    )
                                    .unwrap_or_else(|| candidates.to_vec());
                                    for ec in [
                                        std::borrow::Cow::Borrowed(candidates),
                                        std::borrow::Cow::Owned(candidates.to_vec()),
                                    ] {
                                        let got = super::super::filters::filter_ec_by_strand_cow(
                                            &index, ec, hit, mode, first,
                                        );
                                        assert_eq!(got.as_ref(), expected);
                                    }
                                }
                            }
                        }
                    }
                    let expected = key(full);
                    assert_eq!(key(borrowed), expected);
                    assert_eq!(
                        key(counted),
                        expected,
                        "{} {strand:?}",
                        String::from_utf8_lossy(&read)
                    );
                }
            }
        }
    }
}
