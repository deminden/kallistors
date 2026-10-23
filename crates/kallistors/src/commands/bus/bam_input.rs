//! BAM record ingestion and barcode/UMI tag decoding.

use super::*;

pub(super) fn process_bam_records(
    args: &BusArgs,
    index: &kallistors::pseudoalign::BifrostIndex,
    processing: BusProcessingConfig,
    state: &mut BusWriteState<'_, impl Write>,
    unmapped_ratios: &mut Vec<f64>,
    novel_writer: &mut Option<BufWriter<File>>,
) -> Result<(u64, u64, u64, u64)> {
    if counting::supports(args) {
        return bam_counting::process(args, index, processing, state);
    }
    let bam_path = args
        .reads
        .first()
        .ok_or_else(|| anyhow!("--bam expects exactly one BAM file"))?;
    let file = File::open(bam_path)?;
    let mut reader = noodles_bam::io::Reader::new(file);
    reader
        .read_header()
        .map_err(|err| anyhow!("failed to read BAM header: {err}"))?;

    let mut records = 0u64;
    let mut aligned = 0u64;
    let mut unique = 0u64;
    let mut frame_clashes = 0u64;

    for result in reader.records() {
        let record = result.map_err(|err| anyhow!("failed to read BAM record: {err}"))?;
        if record.flags().is_secondary() || record.flags().is_supplementary() {
            continue;
        }
        if let Some(max_reads) = processing.max_reads
            && records >= max_reads
        {
            break;
        }
        records += 1;

        let sequence = record.sequence().iter().collect::<Vec<_>>();
        let quality = if args.long && args.unmapped {
            record
                .quality_scores()
                .iter()
                .map(|score| score.saturating_add(b'!'))
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let quality = if quality.len() == sequence.len() {
            quality
        } else if args.long && args.unmapped {
            vec![b'!'; sequence.len()]
        } else {
            quality
        };
        let Some(barcode_seq) = bam_string_tag(&record, *b"CR")
            .or_else(|| bam_string_tag(&record, *b"CB").map(normalize_corrected_barcode))
        else {
            record_bam_outputs_for_skipped_read(
                args,
                unmapped_ratios,
                novel_writer,
                index,
                &sequence,
                &quality,
                processing.options,
            )?;
            continue;
        };
        let Some(umi_seq) = bam_string_tag(&record, *b"UR")
            .or_else(|| bam_string_tag(&record, *b"RX"))
            .or_else(|| bam_string_tag(&record, *b"MI").map(normalize_corrected_barcode))
            .or_else(|| bam_string_tag(&record, *b"UB").map(normalize_corrected_barcode))
        else {
            record_bam_outputs_for_skipped_read(
                args,
                unmapped_ratios,
                novel_writer,
                index,
                &sequence,
                &quality,
                processing.options,
            )?;
            continue;
        };
        let too_many_empty_kmers = if args.long || args.unmapped {
            let ratio = kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
                index,
                &sequence,
                processing.options,
            );
            if args.unmapped {
                unmapped_ratios.push(ratio);
            }
            args.long && ratio > processing.threshold
        } else {
            false
        };

        record_length(state.bc_len_hist, barcode_seq.len());
        record_length(state.umi_len_hist, umi_seq.len());

        let mut barcode_flag = 0;
        let barcode = string_to_binary(&barcode_seq, &mut barcode_flag, "barcode")?;
        let mut umi_flag = 0;
        let umi = string_to_binary(&umi_seq, &mut umi_flag, "UMI")?;

        let alignment = if args.aa {
            aa_alignment_for_query(index, &sequence, processing.options)
        } else {
            BusAlignment {
                ec: kallistors::pseudoalign::ec_for_sequence_bifrost(
                    index,
                    &sequence,
                    kallistors::pseudoalign::Strand::Unstranded,
                    processing.options,
                ),
                fragment_length: None,
                placement: None,
                mate_placements: Vec::new(),
                mate_has_ec: Vec::new(),
                frame_clashes: 0,
            }
        };
        let disjoint_intersect = alignment.ec.is_none();
        let novel_long_read = args.long && (too_many_empty_kmers || disjoint_intersect);
        if args.long
            && args.unmapped
            && disjoint_intersect
            && let Some(writer) = novel_writer.as_mut()
        {
            write_novel_read(writer, "unmapped", &sequence, &quality)?;
        }
        if novel_long_read && let Some(writer) = novel_writer.as_mut() {
            let label = if disjoint_intersect {
                "novel_disjointIntersect"
            } else {
                "novel_tooManyEmptyKmers"
            };
            write_novel_read(writer, label, &sequence, &quality)?;
        }
        frame_clashes = frame_clashes.saturating_add(alignment.frame_clashes);
        if args.long
            && !novel_long_read
            && let Some(ec) = alignment.ec.as_deref()
        {
            record_long_lengths(
                state.long_len_sums,
                state.long_len_counts,
                ec,
                sequence.len(),
            );
        }
        let ec = if novel_long_read {
            -1
        } else {
            alignment
                .ec
                .map(|ec| intern_ec(state.ec_map, state.ec_list, ec))
                .unwrap_or(-1)
        };

        if ec >= 0 {
            aligned += 1;
            if is_unique_ec(state.ec_list, ec) {
                unique += 1;
            }
            let flags = barcode_flag | (umi_flag << 8);
            write_bus_record(
                state.writer,
                &BusRecord {
                    barcode,
                    umi,
                    ec,
                    count: 1,
                    flags,
                },
            )?;
        }
    }

    Ok((records, aligned, unique, frame_clashes))
}

pub(super) fn bam_string_tag(record: &noodles_bam::Record, tag: [u8; 2]) -> Option<Vec<u8>> {
    match record.data().get(&tag)?.ok()? {
        SamDataValue::String(value) | SamDataValue::Hex(value) => {
            let bytes: &[u8] = value.as_ref();
            Some(bytes.to_vec())
        }
        _ => None,
    }
}
