mod alignment;
mod bam_counting;
mod bam_input;
mod configuration;
mod counting;
mod normalized;
mod tags;

use alignment::{
    aa_alignment_for_query, bus_query_quality, bus_query_sequence, pseudoalign_bus_read,
    record_bam_outputs_for_skipped_read, record_unmapped_ratio_for_skipped_read,
};
use bam_input::process_bam_records;
use configuration::{
    BusConfiguration, bus_pseudoalign_options, bus_strand_specific, load_index,
    long_read_threshold, prepare, print_technology_list, report_default_strand,
    report_missing_strand_annotations, report_threshold_choice, user_specified_strand,
};
use tags::{
    barcode_sequence, configure_tag, normalize_corrected_barcode, string_to_binary, umi_value,
};

use anyhow::{Result, anyhow, bail};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::time::Instant;

use noodles_sam::alignment::record::data::field::Value as SamDataValue;

mod fast;
mod genome;
mod input;
mod output;
mod pseudobam;
mod technology;

use genome::{GenomeModel, write_gene_list};

use pseudobam::{PseudoBamBusRead, PseudoBamWriter, write_pseudobam_reads};

use input::{
    next_interleaved_batch, next_parallel_batch, read_groups, report_verbose_read_groups,
    should_infer_bulk_paired, validate_read_files,
};

use output::{
    BusRecord, RunInfo, fake_barcode, write_batch_long_flens, write_batch_paired_flens,
    write_bus_header, write_bus_record, write_cells, write_long_flens, write_matrix_ec,
    write_novel_read, write_paired_flens, write_run_info, write_sample_barcodes,
    write_unmapped_ratios,
};

#[cfg(test)]
use technology::concat_slices;

use technology::{
    SliceSpec, TechnologySpec, add_unpaired_mate_slice, concat_quality_slices,
    concat_sequence_slices, concat_slices_optional, extract_quality_slice, extract_sequence_slice,
    extract_slice_optional, mark_technology_paired, technology_base_upper, technology_spec,
    total_len,
};

#[derive(Default)]
pub struct BusArgs {
    pub index: Option<PathBuf>,
    pub out_dir: Option<PathBuf>,
    pub technology: Option<String>,
    pub reads: Vec<PathBuf>,
    pub threads: usize,
    pub list: bool,
    pub max_reads: Option<u64>,
    pub batch: Option<PathBuf>,
    pub bam: bool,
    pub pseudobam: bool,
    pub genomebam: bool,
    pub gtf: Option<PathBuf>,
    pub chromosomes: Option<PathBuf>,
    pub batch_barcodes: bool,
    pub interleaved: bool,
    pub tag_sequence: Option<String>,
    pub long: bool,
    pub platform: Option<String>,
    pub error_rate: Option<f64>,
    pub threshold: Option<f64>,
    pub unmapped: bool,
    pub aa: bool,
    pub num: bool,
    pub paired: bool,
    pub unstranded: bool,
    pub fr_stranded: bool,
    pub rf_stranded: bool,
    pub do_union: bool,
    pub no_jump: bool,
    pub verbose: bool,
    pub dfk_onlist: bool,
    pub kallisto_enum: bool,
    pub kallisto_strict: bool,
    pub kallisto_local_fallback: bool,
    pub kallisto_fallback: bool,
    pub discard_special_only: bool,
    pub skip_overcrowded_minimizer: bool,
    pub kallisto_direct_kmer: bool,
    pub kallisto_bifrost_find: bool,
    pub kallisto_sparse_hits: bool,
}

#[derive(Debug, Clone)]
struct TagConfig {
    sequence: Vec<u8>,
    binary: u64,
}

struct ConfiguredTag {
    config: TagConfig,
    inferred_default: bool,
}

struct UmiValue {
    binary: u64,
    flag: u32,
    len: usize,
    seq: Vec<u8>,
    tag_present: bool,
    disable_strand_specific: bool,
}

struct BusWriteState<'a, W: Write> {
    writer: &'a mut W,
    ec_map: &'a mut HashMap<Vec<u32>, i32>,
    ec_list: &'a mut Vec<Vec<u32>>,
    bc_len_hist: &'a mut [u64; 33],
    umi_len_hist: &'a mut [u64; 33],
    long_len_sums: &'a mut [u64],
    long_len_counts: &'a mut [u64],
    paired_flens: &'a mut [u32],
    batch_paired_flens: &'a mut [Vec<u32>],
}

#[derive(Debug, Clone, Copy)]
struct BusProcessingConfig {
    options: kallistors::pseudoalign::PseudoalignOptions,
    threshold: f64,
    max_reads: Option<u64>,
}

struct BusAlignment {
    ec: Option<Vec<u32>>,
    fragment_length: Option<i64>,
    placement: Option<kallistors::pseudoalign::PseudoalignPlacement>,
    mate_placements: Vec<Option<kallistors::pseudoalign::PseudoalignPlacement>>,
    mate_has_ec: Vec<bool>,
    frame_clashes: u64,
}

pub fn run(args: BusArgs) -> Result<()> {
    if args.list {
        print_technology_list();
        return Ok(());
    }
    let BusConfiguration {
        index_path,
        out_dir,
        spec,
        tag_config,
        read_groups,
        write_pseudobam,
        synthetic_bc,
        fixed_bc_len,
        fixed_umi_len,
        batch_barcode_prefix_len,
        bc_len,
    } = prepare(&args)?;

    let start = Instant::now();
    let mut index = load_index(&args, &spec, index_path)?;
    let threshold = if args.long {
        let choice = long_read_threshold(args.threshold, args.error_rate, index.k);
        if let Some(report) = choice.report.as_ref() {
            report_threshold_choice(report);
        }
        choice.value
    } else {
        0.8
    };
    let strand_specific = bus_strand_specific(&args, &spec, &index);
    if strand_specific.is_some() && !strand_specific_compatible(&spec) {
        bail!(
            "Strand-specific read processing is only supported for technologies with a single cDNA read file or paired-end reads"
        );
    }
    if !user_specified_strand(&args)
        && let Some(strand) = strand_specific
        && spec.default_strand == Some(strand)
    {
        report_default_strand(strand);
    } else if !user_specified_strand(&args)
        && strand_specific.is_none()
        && spec.default_strand.is_some()
        && !index.has_strand_annotations()
    {
        report_missing_strand_annotations();
    }
    let options = bus_pseudoalign_options(&args, strand_specific);
    let processing = BusProcessingConfig {
        options,
        threshold,
        max_reads: args.max_reads.filter(|&value| value != 0),
    };

    fs::create_dir_all(out_dir)?;
    let mut bus_writer = BufWriter::new(File::create(out_dir.join("output.bus"))?);
    let genome_model = if args.genomebam {
        Some(GenomeModel::from_gtf(
            args.gtf.as_deref().expect("validated"),
            args.chromosomes.as_deref(),
            &index,
        )?)
    } else {
        None
    };
    if let Some(model) = genome_model.as_ref() {
        write_gene_list(&out_dir.join("matrix.genelist.txt"), model)?;
    }
    let mut pseudo_bam = if write_pseudobam {
        Some(PseudoBamWriter::create(out_dir, &index, genome_model)?)
    } else {
        None
    };
    let umi_len = if args.bam {
        0
    } else if spec.bulk_like() {
        1
    } else {
        fixed_umi_len.unwrap_or(0)
    };
    let header_bc_len = if args.bam { 0 } else { bc_len };
    write_bus_header(&mut bus_writer, header_bc_len, umi_len)?;

    let mut ec_map: HashMap<Vec<u32>, i32> = HashMap::new();
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut records = 0u64;
    let mut aligned = 0u64;
    let mut unique = 0u64;
    let mut bc_len_hist = [0u64; 33];
    let mut umi_len_hist = [0u64; 33];
    let mut unmapped_ratios = Vec::new();
    let long_targets = if args.long {
        index.transcript_names.len()
    } else {
        0
    };
    let mut long_len_sums = vec![0u64; long_targets];
    let mut long_len_counts = vec![0u64; long_targets];
    let mut paired_flens = vec![0u32; kallistors::pseudoalign::MAX_FRAG_LEN as usize];
    let batch_count = read_groups
        .iter()
        .map(|group| group.batch_index)
        .max()
        .map_or(0, |idx| idx + 1);
    let mut batch_long_len_sums = vec![
        vec![0u64; long_targets];
        if args.long && args.batch.is_some() {
            batch_count
        } else {
            0
        }
    ];
    let mut batch_long_len_counts = vec![
        vec![0u64; long_targets];
        if args.long && args.batch.is_some() {
            batch_count
        } else {
            0
        }
    ];
    let mut batch_paired_flens =
        vec![vec![0u32; kallistors::pseudoalign::MAX_FRAG_LEN as usize]; batch_count];
    let mut frame_clashes = 0u64;
    let mut novel_writer = if args.long && args.unmapped {
        Some(BufWriter::new(File::create(out_dir.join("novel.fastq"))?))
    } else {
        None
    };

    if args.bam {
        let mut state = BusWriteState {
            writer: &mut bus_writer,
            ec_map: &mut ec_map,
            ec_list: &mut ec_list,
            bc_len_hist: &mut bc_len_hist,
            umi_len_hist: &mut umi_len_hist,
            long_len_sums: &mut long_len_sums,
            long_len_counts: &mut long_len_counts,
            paired_flens: &mut paired_flens,
            batch_paired_flens: &mut batch_paired_flens,
        };
        let bam_records = process_bam_records(
            &args,
            &index,
            processing,
            &mut state,
            &mut unmapped_ratios,
            &mut novel_writer,
        )?;
        records = bam_records.0;
        aligned = bam_records.1;
        unique = bam_records.2;
        frame_clashes = bam_records.3;
    } else if fast::supports(&args, &spec) {
        // Ordinary BUS counting applies strand filtering in the counting API;
        // it does not need the serial diagnostic/placement path below.
        (records, aligned, unique) = fast::process(
            &args,
            &spec,
            &index,
            &read_groups,
            processing,
            &mut BusWriteState {
                writer: &mut bus_writer,
                ec_map: &mut ec_map,
                ec_list: &mut ec_list,
                bc_len_hist: &mut bc_len_hist,
                umi_len_hist: &mut umi_len_hist,
                long_len_sums: &mut long_len_sums,
                long_len_counts: &mut long_len_counts,
                paired_flens: &mut paired_flens,
                batch_paired_flens: &mut batch_paired_flens,
            },
        )?;
    } else if normalized::supports(&args) {
        (records, aligned, unique) = normalized::process(
            &args,
            &spec,
            &index,
            &read_groups,
            processing,
            tag_config.as_ref().map(|tag| &tag.config),
            &mut BusWriteState {
                writer: &mut bus_writer,
                ec_map: &mut ec_map,
                ec_list: &mut ec_list,
                bc_len_hist: &mut bc_len_hist,
                umi_len_hist: &mut umi_len_hist,
                long_len_sums: &mut long_len_sums,
                long_len_counts: &mut long_len_counts,
                paired_flens: &mut paired_flens,
                batch_paired_flens: &mut batch_paired_flens,
            },
        )?;
    } else {
        'groups: for group in &read_groups {
            let mut readers = group
                .files
                .iter()
                .map(|path| kallistors::io::open_fastq_reader(path))
                .collect::<kallistors::Result<Vec<_>>>()
                .map_err(|err| anyhow!("failed to open FASTQ: {err}"))?;

            loop {
                let batch = if group.interleaved {
                    next_interleaved_batch(&mut readers, spec.nfiles)?
                } else {
                    next_parallel_batch(&mut readers)?
                };
                if batch.is_empty() {
                    break;
                }
                if let Some(max_reads) = processing.max_reads
                    && records >= max_reads
                {
                    break 'groups;
                }
                records += 1;

                let Some(barcode_seq) =
                    barcode_sequence(&batch, &spec, group.batch_index, args.batch_barcodes)?
                else {
                    record_unmapped_ratio_for_skipped_read(
                        args.unmapped,
                        &mut unmapped_ratios,
                        &index,
                        &batch,
                        &spec,
                        options,
                    )?;
                    continue;
                };
                let Some(umi_value) =
                    umi_value(&batch, &spec, tag_config.as_ref().map(|tag| &tag.config))?
                else {
                    record_unmapped_ratio_for_skipped_read(
                        args.unmapped,
                        &mut unmapped_ratios,
                        &index,
                        &batch,
                        &spec,
                        options,
                    )?;
                    continue;
                };
                let umi_seq_len = umi_value.len;
                record_length(&mut bc_len_hist, barcode_seq.len());
                record_length(&mut umi_len_hist, umi_seq_len);
                let mut barcode_flag = 0;
                let barcode = string_to_binary(&barcode_seq, &mut barcode_flag, "barcode")?;

                let query = if args.long || args.unmapped {
                    Some(bus_query_sequence(&batch, &spec, umi_value.tag_present)?)
                } else {
                    None
                };
                let query_quality = if args.long || args.unmapped {
                    Some(bus_query_quality(&batch, &spec, umi_value.tag_present)?)
                } else {
                    None
                };
                let too_many_empty_kmers = if let Some(query) = query.as_deref() {
                    let ratio = kallistors::pseudoalign::unmapped_kmer_ratio_bifrost(
                        &index, query, options,
                    );
                    if args.unmapped {
                        unmapped_ratios.push(ratio);
                    }
                    args.long && ratio > threshold
                } else {
                    false
                };
                let mut read_options = options;
                if umi_value.disable_strand_specific {
                    read_options.strand_specific = None;
                }
                let alignment = pseudoalign_bus_read(
                    &index,
                    &batch,
                    &spec,
                    read_options,
                    umi_value.tag_present,
                    args.aa,
                )?;
                let disjoint_intersect = alignment.ec.is_none();
                let novel_long_read = args.long && (too_many_empty_kmers || disjoint_intersect);
                if args.long
                    && args.unmapped
                    && disjoint_intersect
                    && let (Some(writer), Some(query), Some(quality)) = (
                        novel_writer.as_mut(),
                        query.as_deref(),
                        query_quality.as_deref(),
                    )
                {
                    write_novel_read(writer, "unmapped", query, quality)?;
                }
                if novel_long_read
                    && let (Some(writer), Some(query), Some(quality)) = (
                        novel_writer.as_mut(),
                        query.as_deref(),
                        query_quality.as_deref(),
                    )
                {
                    let label = if disjoint_intersect {
                        "novel_disjointIntersect"
                    } else {
                        "novel_tooManyEmptyKmers"
                    };
                    write_novel_read(writer, label, query, quality)?;
                }
                if args.long
                    && !novel_long_read
                    && let Some(ec) = alignment.ec.as_deref()
                {
                    record_long_lengths(
                        &mut long_len_sums,
                        &mut long_len_counts,
                        ec,
                        query.as_ref().map_or(0, Vec::len),
                    );
                    if let (Some(batch_sums), Some(batch_counts)) = (
                        batch_long_len_sums.get_mut(group.batch_index),
                        batch_long_len_counts.get_mut(group.batch_index),
                    ) {
                        record_long_lengths(
                            batch_sums,
                            batch_counts,
                            ec,
                            query.as_ref().map_or(0, Vec::len),
                        );
                    }
                }
                if !args.long
                    && spec.paired
                    && let (Some(ec), Some(fragment_length)) =
                        (alignment.ec.as_deref(), alignment.fragment_length)
                {
                    record_paired_fragment_length(&mut paired_flens, ec, fragment_length);
                    if let Some(batch_flens) = batch_paired_flens.get_mut(group.batch_index) {
                        record_paired_fragment_length(batch_flens, ec, fragment_length);
                    }
                }
                let pseudo_bam_placement = alignment.placement;
                let pseudo_bam_mate_placements = alignment.mate_placements.clone();
                frame_clashes = frame_clashes.saturating_add(alignment.frame_clashes);
                let ec = if novel_long_read {
                    -1
                } else {
                    alignment
                        .ec
                        .map(|ec| intern_ec(&mut ec_map, &mut ec_list, ec))
                        .unwrap_or(-1)
                };
                if ec >= 0 {
                    aligned += 1;
                    if is_unique_ec(&ec_list, ec) {
                        unique += 1;
                    }
                    let flags = if args.num {
                        u32::try_from(records - 1).unwrap_or(u32::MAX)
                    } else {
                        barcode_flag | (umi_value.flag << 8)
                    };
                    write_bus_record(
                        &mut bus_writer,
                        &BusRecord {
                            barcode,
                            umi: umi_value.binary,
                            ec,
                            count: 1,
                            flags,
                        },
                    )?;
                    if let Some(writer) = pseudo_bam.as_mut() {
                        write_pseudobam_reads(
                            writer,
                            PseudoBamBusRead {
                                read_id: records - 1,
                                ec_id: ec,
                                ec_list: &ec_list,
                                batch: &batch,
                                spec: &spec,
                                tag_present: umi_value.tag_present,
                                barcode: &barcode_seq,
                                umi: &umi_value.seq,
                                placement: pseudo_bam_placement,
                                mate_placements: &pseudo_bam_mate_placements,
                                mate_has_ec: &alignment.mate_has_ec,
                            },
                        )?;
                    }
                }
            }
        }
    }
    bus_writer.flush()?;
    let mut bus_file = bus_writer.into_inner()?;
    if header_bc_len == 0 {
        let observed_bc_len = modal_length(&bc_len_hist).or(fixed_bc_len);
        if let Some(observed_bc_len) = observed_bc_len {
            bus_file.seek(SeekFrom::Start(8))?;
            bus_file.write_all(&observed_bc_len.to_le_bytes())?;
        }
    }
    if umi_len == 0 {
        let observed_umi_len = modal_length(&umi_len_hist).or(fixed_umi_len);
        if let Some(observed_umi_len) = observed_umi_len {
            bus_file.seek(SeekFrom::Start(12))?;
            bus_file.write_all(&observed_umi_len.to_le_bytes())?;
        }
    }

    write_matrix_ec(&out_dir.join("matrix.ec"), &ec_list)?;
    if args.unmapped {
        write_unmapped_ratios(
            &out_dir.join("unmapped_ratio.txt"),
            &unmapped_ratios,
            args.batch.is_some(),
        )?;
    }
    if args.long {
        if args.batch.is_some() {
            write_batch_long_flens(
                &out_dir.join("flens.txt"),
                &index.transcript_lengths,
                &batch_long_len_sums,
                &batch_long_len_counts,
                index.k,
            )?;
        } else {
            write_long_flens(
                &out_dir.join("flens.txt"),
                &index.transcript_lengths,
                &long_len_sums,
                &long_len_counts,
                index.k,
            )?;
        }
    } else if spec.paired {
        if args.batch.is_some() {
            write_batch_paired_flens(&out_dir.join("flens.txt"), &batch_paired_flens)?;
        } else {
            write_paired_flens(&out_dir.join("flens.txt"), &paired_flens)?;
        }
    }
    if args.long || spec.paired || no_umi(&spec) {
        fs::copy(index_path, out_dir.join("index.saved"))?;
    }
    if let Some(mut writer) = novel_writer {
        writer.flush()?;
    }
    if let Some(writer) = pseudo_bam {
        writer.finish()?;
    }
    if args.batch.is_some() || synthetic_bc {
        write_cells(
            &out_dir.join("matrix.cells"),
            read_groups
                .iter()
                .enumerate()
                .map(|(idx, group)| (idx, group.id.as_deref())),
        )?;
        if synthetic_bc || args.batch_barcodes {
            let sample_barcode_len = batch_barcode_prefix_len.unwrap_or(16);
            write_sample_barcodes(
                &out_dir.join("matrix.sample.barcodes"),
                read_groups.iter().map(|group| group.batch_index as u64),
                sample_barcode_len,
            )?;
        }
    }
    let mut transcripts_writer = BufWriter::new(File::create(out_dir.join("transcripts.txt"))?);
    let targets = index
        .take_transcript_names()
        .write_lines(&mut transcripts_writer, index.onlist.as_deref())?;
    transcripts_writer.flush()?;
    write_run_info(
        &out_dir.join("run_info.json"),
        RunInfo {
            processed: records,
            aligned,
            unique,
            targets,
            k: index.k,
            frame_clashes: args.aa.then_some(frame_clashes),
        },
    )?;

    println!(
        "processed: {records}, aligned: {aligned}, ecs: {}, time: {:.2?}",
        ec_list.len(),
        start.elapsed()
    );
    if let Some(max_reads) = processing.max_reads
        && records < max_reads
    {
        eprintln!(
            "Note: Number of reads processed is less than --numReads: {max_reads}, returning 1"
        );
        bail!("Number of reads processed is less than --numReads: {max_reads}, returning 1");
    }
    if aligned == 0 {
        bail!("zero reads pseudoaligned");
    }
    Ok(())
}

fn is_unique_ec(ec_list: &[Vec<u32>], ec: i32) -> bool {
    usize::try_from(ec)
        .ok()
        .and_then(|idx| ec_list.get(idx))
        .is_some_and(|ec| ec.len() == 1)
}

fn no_umi(spec: &TechnologySpec) -> bool {
    spec.umi.first().is_some_and(|spec| spec.file.is_none())
}

fn strand_specific_compatible(spec: &TechnologySpec) -> bool {
    spec.seq.len() == 1 || (spec.seq.len() == 2 && spec.paired)
}

fn pseudobam_compatible(spec: &TechnologySpec) -> bool {
    strand_specific_compatible(spec)
}

fn record_length(hist: &mut [u64; 33], len: usize) {
    hist[len.min(32)] = hist[len.min(32)].saturating_add(1);
}

fn modal_length(hist: &[u64; 33]) -> Option<u32> {
    hist.iter()
        .enumerate()
        .max_by_key(|&(len, count)| (*count, std::cmp::Reverse(len)))
        .and_then(|(len, &count)| (count > 0).then_some(len as u32))
}

fn intern_ec(
    ec_map: &mut HashMap<Vec<u32>, i32>,
    ec_list: &mut Vec<Vec<u32>>,
    ec: Vec<u32>,
) -> i32 {
    if let Some(id) = ec_map.get(&ec) {
        *id
    } else {
        let id = i32::try_from(ec_list.len()).unwrap_or(i32::MAX);
        ec_map.insert(ec.clone(), id);
        ec_list.push(ec);
        id
    }
}

fn record_long_lengths(sums: &mut [u64], counts: &mut [u64], ec: &[u32], read_len: usize) {
    if ec.len() != 1 {
        return;
    }
    let Ok(transcript_id) = usize::try_from(ec[0]) else {
        return;
    };
    let Some(sum) = sums.get_mut(transcript_id) else {
        return;
    };
    *sum = sum.saturating_add(read_len as u64);
    if let Some(count) = counts.get_mut(transcript_id) {
        *count = count.saturating_add(1);
    }
}

fn record_paired_fragment_length(flens: &mut [u32], ec: &[u32], fragment_length: i64) {
    if ec.len() != 1 {
        return;
    }
    let Ok(idx) = usize::try_from(fragment_length) else {
        return;
    };
    let Some(bin) = flens.get_mut(idx) else {
        return;
    };
    *bin = bin.saturating_add(1);
}

#[cfg(test)]
mod tests;
