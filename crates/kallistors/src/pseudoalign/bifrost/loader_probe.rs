//! One immutable unitig view per ordered loader chunk, with no retained table.

use super::encoded::{EncodedSequence, EncodedUnitigs};

pub(super) struct CachedUnitigMinimizers<'a> {
    unitigs: &'a EncodedUnitigs,
    last: usize,
    sequence: Option<EncodedSequence<'a>>,
}

impl<'a> CachedUnitigMinimizers<'a> {
    pub(super) fn new(unitigs: &'a EncodedUnitigs) -> Self {
        Self {
            unitigs,
            last: usize::MAX,
            sequence: None,
        }
    }

    #[inline]
    pub(super) fn at(&mut self, unitig: usize, position: usize) -> Option<[u8; 8]> {
        self.at_g(unitig, position, self.unitigs.g)
    }

    // The large ordinary loader proves g=23 before selecting this entry point.
    // Keep the same metadata view and checked window helper, with a constant
    // width so bounds, shifts and complement masks need no runtime g dispatch.
    #[inline]
    pub(super) fn at23(&mut self, unitig: usize, position: usize) -> Option<[u8; 8]> {
        self.at_g(unitig, position, 23)
    }

    #[inline(always)]
    fn at_g(&mut self, unitig: usize, position: usize, g: usize) -> Option<[u8; 8]> {
        if unitig != self.last {
            self.sequence = self.unitigs.unitig_bases.sequence(unitig);
            self.last = unitig;
        }
        self.sequence.as_ref()?.minimizer_at(position, g)
    }
}

#[cfg(test)]
mod tests {
    use super::super::{SequenceMetadata, precompute_encoded_unitigs};
    use super::*;

    #[test]
    fn cached_loader_views_preserve_metadata_transitions_and_query_tails() {
        let sequences: Vec<Vec<u8>> = [0, 31, 32, 33, 65]
            .into_iter()
            .map(|n| (0..n).map(|i| b"ACGT"[(i * 13 + n) % 4]).collect())
            .collect();
        for g in [0, 1, 23, 31, 32, 33] {
            for layout in 0..3 {
                let mut encoded = precompute_encoded_unitigs(&sequences, &[], g, None);
                if layout == 1 {
                    let SequenceMetadata::Compact(records) = &encoded.unitig_bases.metadata else {
                        panic!("small fixture");
                    };
                    encoded.unitig_bases.metadata = SequenceMetadata::Wide(
                        records
                            .iter()
                            .map(|&v| (v as u32 as usize, (v >> 32) as usize))
                            .collect(),
                    );
                } else if layout == 2 {
                    encoded.unitig_bases.inline_short_words();
                }
                let mut probe = CachedUnitigMinimizers::new(&encoded);
                for id in [usize::MAX, 0, 1, 1, 2, 4, 4, 3, 0, 7, 7, 1, usize::MAX] {
                    let len = sequences.get(id).map_or(0, Vec::len);
                    for position in (0..len + 2).chain(std::iter::once(usize::MAX)) {
                        if g == 23 {
                            assert_eq!(
                                probe.at23(id, position),
                                encoded.unitig_minimizer_at(id, position),
                                "fixed23 layout={layout} unitig={id} position={position}"
                            );
                        }
                        assert_eq!(
                            probe.at(id, position),
                            encoded.unitig_minimizer_at(id, position),
                            "layout={layout} g={g} unitig={id} position={position}"
                        );
                    }
                }
            }
        }
    }
}
