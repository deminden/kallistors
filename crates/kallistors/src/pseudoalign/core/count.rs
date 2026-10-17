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
    pub(super) fn push(&mut self, index: &BifrostIndex, hit: Hit) {
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
        let index = crate::pseudoalign::build_bifrost_index_threaded(&path, 2).unwrap();
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
                let full =
                    ec_for_read_bifrost(&index, &read, strand, None, PseudoalignOptions::default());
                super::super::runtime::reset_thread_local_caches();
                let counted = ec_for_read_bifrost_counting(
                    &index,
                    &read,
                    strand,
                    PseudoalignOptions::default(),
                );
                assert_eq!(
                    key(counted),
                    key(full),
                    "{} {strand:?}",
                    String::from_utf8_lossy(&read)
                );
            }
        }
    }
}
