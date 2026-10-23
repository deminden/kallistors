//! BUS query assembly, alignment and skipped-read diagnostics.

use super::*;

pub(super) fn record_unmapped_ratio_for_skipped_read(
    enabled: bool,
    ratios: &mut Vec<f64>,
    index: &kallistors::pseudoalign::BifrostIndex,
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> Result<()> {
    if enabled {
        let query = bus_query_sequence(batch, spec, false)?;
        ratios.push(kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
            index, &query, options,
        ));
    }
    Ok(())
}

pub(super) fn record_bam_outputs_for_skipped_read(
    args: &BusArgs,
    ratios: &mut Vec<f64>,
    novel_writer: &mut Option<BufWriter<File>>,
    index: &kallistors::pseudoalign::BifrostIndex,
    sequence: &[u8],
    quality: &[u8],
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> Result<()> {
    if args.unmapped {
        ratios.push(kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
            index, sequence, options,
        ));
    }
    if args.long
        && args.unmapped
        && let Some(writer) = novel_writer.as_mut()
    {
        write_novel_read(writer, "skipped_missingTags", sequence, quality)?;
    }
    Ok(())
}

pub(super) fn pseudoalign_bus_read(
    index: &kallistors::pseudoalign::BifrostIndex,
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    options: kallistors::pseudoalign::PseudoalignOptions,
    tag_present: bool,
    aa: bool,
) -> Result<BusAlignment> {
    if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices");
        }
        let left = extract_sequence_slice(batch, spec.seq[0], spec, tag_present)?;
        let right = extract_sequence_slice(batch, spec.seq[1], spec, tag_present)?;
        let trace = kallistors::pseudoalign::trace_read_pair_bifrost(
            index,
            &left,
            &right,
            kallistors::pseudoalign::Strand::Unstranded,
            options,
        );
        Ok(BusAlignment {
            ec: trace.merged_ec.filter(|ec| !ec.is_empty()),
            fragment_length: trace.estimated_fragment_length,
            placement: trace.placement,
            mate_placements: vec![trace.left.placement, trace.right.placement],
            mate_has_ec: vec![
                trace
                    .left
                    .ec_after_strand_filter
                    .as_ref()
                    .is_some_and(|ec| !ec.is_empty()),
                trace
                    .right
                    .ec_after_strand_filter
                    .as_ref()
                    .is_some_and(|ec| !ec.is_empty()),
            ],
            frame_clashes: 0,
        })
    } else {
        let query = concat_sequence_slices(batch, spec, tag_present)?;
        if aa {
            return Ok(aa_alignment_for_query(index, &query, options));
        }
        let trace = kallistors::pseudoalign::trace_read_bifrost(
            index,
            &query,
            kallistors::pseudoalign::Strand::Unstranded,
            None,
            options,
        );
        Ok(BusAlignment {
            // Traces retain an empty EC to explain a rejected intersection.
            // BUS records require at least one compatible transcript.
            ec: trace.ec_after_strand_filter.filter(|ec| !ec.is_empty()),
            fragment_length: None,
            placement: trace.placement,
            mate_placements: Vec::new(),
            mate_has_ec: Vec::new(),
            frame_clashes: 0,
        })
    }
}

pub(super) fn aa_alignment_for_query(
    index: &kallistors::pseudoalign::BifrostIndex,
    query: &[u8],
    options: kallistors::pseudoalign::PseudoalignOptions,
) -> BusAlignment {
    let mut best_ec: Option<Vec<u32>> = None;
    let mut best_cardinality = usize::MAX;
    let mut frame_clashes = 0u64;
    let mut frame = Vec::new();
    for seq in [query.to_vec(), kallistors::util::reverse_complement(query)] {
        for offset in 0..3 {
            if offset >= seq.len() {
                continue;
            }
            frame.clear();
            frame.extend_from_slice(&kallistors::util::nucleotide_to_comma_free(&seq[offset..]));
            if frame.len() < index.k {
                continue;
            }
            if let Some(ec) = kallistors::pseudoalign::ec_for_sequence_bifrost(
                index,
                &frame,
                kallistors::pseudoalign::Strand::Unstranded,
                options,
            ) {
                match ec.len().cmp(&best_cardinality) {
                    std::cmp::Ordering::Less => {
                        best_cardinality = ec.len();
                        best_ec = Some(ec);
                    }
                    std::cmp::Ordering::Equal => {
                        frame_clashes = frame_clashes.saturating_add(1);
                    }
                    std::cmp::Ordering::Greater => {}
                }
            }
        }
    }
    BusAlignment {
        ec: best_ec,
        fragment_length: None,
        placement: None,
        mate_placements: Vec::new(),
        mate_has_ec: Vec::new(),
        frame_clashes,
    }
}

pub(super) fn bus_query_sequence(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices");
        }
        let mut query = extract_sequence_slice(batch, spec.seq[0], spec, tag_present)?;
        query.push(b'N');
        query.extend_from_slice(&extract_sequence_slice(
            batch,
            spec.seq[1],
            spec,
            tag_present,
        )?);
        Ok(query)
    } else {
        concat_sequence_slices(batch, spec, tag_present)
    }
}

pub(super) fn bus_query_quality(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag_present: bool,
) -> Result<Vec<u8>> {
    if spec.paired {
        if spec.seq.len() != 2 {
            bail!("paired BUS technology requires exactly two sequence slices");
        }
        let mut quality = extract_quality_slice(batch, spec.seq[0], spec, tag_present)?;
        quality.push(b'!');
        quality.extend_from_slice(&extract_quality_slice(
            batch,
            spec.seq[1],
            spec,
            tag_present,
        )?);
        Ok(quality)
    } else {
        concat_quality_slices(batch, spec, tag_present)
    }
}
