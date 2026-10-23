//! Pair ECs and fragment lengths without diagnostic traces or placement output.

use super::{BifrostIndex, PseudoalignOptions, ReadEc, Strand};

pub struct PairedCountingResult {
    pub ec: Vec<u32>,
    pub fragment_length: Option<i64>,
}

/// Preserve BUS pair rejection and complete fragment-length observations.
/// Unlike bulk fragment sampling, every eligible pair returns its length.
pub fn ec_for_pair_bifrost_counting(
    index: &BifrostIndex,
    left: &[u8],
    right: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<PairedCountingResult> {
    let scan = |sequence: &[u8], first| {
        let mut read = super::ec_for_read_bifrost_single(index, sequence, strand, options)?;
        if let Some(mode) = options.strand_specific {
            let comprehensive = options.do_union || options.no_jump || index.use_shade;
            if let Some(filtered) =
                super::apply_strand_filter(index, &read, mode, first, comprehensive, &mut None, b"")
            {
                read.hard_reject_pair |= filtered.is_empty();
                read.ec = filtered;
            }
        }
        Some(read)
    };
    let mut a: Option<ReadEc> = scan(left, true);
    let mut b: Option<ReadEc> = scan(right, false);
    if a.as_ref().is_some_and(|v| v.hard_reject_pair)
        || b.as_ref().is_some_and(|v| v.hard_reject_pair)
    {
        return None;
    }
    if a.as_ref().is_some_and(|v| v.ec.is_empty()) {
        a = None;
    }
    if b.as_ref().is_some_and(|v| v.ec.is_empty()) {
        b = None;
    }
    let mut ec = Vec::new();
    match (&a, &b) {
        (None, None) => return None,
        (Some(read), None) | (None, Some(read)) => ec.extend_from_slice(&read.ec),
        (Some(a), Some(b)) => super::intersect_sorted(&a.ec, &b.ec, &mut ec),
    }
    if ec.is_empty() {
        return None;
    }
    let fragment_length = if ec.len() == 1 {
        a.as_ref()
            .zip(b.as_ref())
            .and_then(|(a, b)| super::paired::estimate_fragment_length_for_pair(index, ec[0], a, b))
    } else {
        None
    };
    Some(PairedCountingResult {
        ec,
        fragment_length,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn counting_matches_complete_traces_with_strands_positions_and_rejections() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("reference.fa");
        let path = dir.path().join("reference.idx");
        let mut state = 181u64;
        let mut generate = || {
            (0..400)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect::<Vec<_>>()
        };
        let a = generate();
        let mut b = generate();
        b[..100].copy_from_slice(&a[..100]);
        let reverse = crate::util::reverse_complement(&b);
        for shade in [false, true] {
            let shade_record = if shade {
                format!(">a_shade_alt\n{}\n", String::from_utf8_lossy(&a))
            } else {
                String::new()
            };
            std::fs::write(
                &fasta,
                format!(
                    ">a\n{}\n>b\n{}\n>reverse\n{}\n{}",
                    String::from_utf8_lossy(&a),
                    String::from_utf8_lossy(&b),
                    String::from_utf8_lossy(&reverse),
                    shade_record
                ),
            )
            .unwrap();
            crate::index::build_index(&path, std::slice::from_ref(&fasta), Default::default())
                .unwrap();
            let reads = [
                a[25..105].to_vec(),
                a[125..205].to_vec(),
                b[125..205].to_vec(),
                crate::util::reverse_complement(&a[180..260]),
                b[..20].to_vec(),
                vec![b'N'; 100],
                [&a[20..85], &b[210..275]].concat(),
            ];
            for index in [
                super::super::build_bifrost_index_with_positions(&path, true).unwrap(),
                super::super::build_bifrost_index_positional_packed_threaded(&path, 2, false)
                    .unwrap(),
            ] {
                assert!(index.has_strand_annotations());
                for mode in [
                    None,
                    Some(super::super::StrandSpecific::FR),
                    Some(super::super::StrandSpecific::RF),
                ] {
                    for policy in 0..3 {
                        let options = PseudoalignOptions {
                            strand_specific: mode,
                            do_union: policy == 1,
                            no_jump: policy == 2,
                            ..Default::default()
                        };
                        for left in &reads {
                            let reference = super::super::trace_read_bifrost(
                                &index,
                                left,
                                Strand::Unstranded,
                                None,
                                options,
                            )
                            .ec_after_strand_filter
                            .filter(|ec| !ec.is_empty());
                            let counted = super::super::ec_for_sequence_bifrost_counting(
                                &index,
                                left,
                                Strand::Unstranded,
                                options,
                            );
                            assert_eq!(counted.as_deref(), reference.as_deref());
                            for right in &reads {
                                let reference = super::super::trace_read_pair_bifrost(
                                    &index,
                                    left,
                                    right,
                                    Strand::Unstranded,
                                    options,
                                );
                                let counted = ec_for_pair_bifrost_counting(
                                    &index,
                                    left,
                                    right,
                                    Strand::Unstranded,
                                    options,
                                );
                                assert_eq!(
                                    counted.as_ref().map(|v| &v.ec),
                                    reference.merged_ec.as_ref()
                                );
                                assert_eq!(
                                    counted.and_then(|v| v.fragment_length),
                                    reference.estimated_fragment_length
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
