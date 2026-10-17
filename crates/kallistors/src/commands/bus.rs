use anyhow::{Result, anyhow, bail};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::time::Instant;

use noodles_sam::alignment::record::data::field::Value as SamDataValue;

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
    write_novel_read, write_paired_flens, write_run_info, write_sample_barcodes, write_transcripts,
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

struct ThresholdChoice {
    value: f64,
    report: Option<ThresholdReport>,
}

enum ThresholdReport {
    Computed(f64),
    InvalidSupplied,
    InvalidComputed(f64),
}

fn long_read_threshold(
    threshold: Option<f64>,
    error_rate: Option<f64>,
    k: usize,
) -> ThresholdChoice {
    if let Some(threshold) = threshold {
        if 0.0 < threshold && threshold < 1.0 {
            return ThresholdChoice {
                value: threshold,
                report: None,
            };
        }
        return ThresholdChoice {
            value: 0.8,
            report: Some(ThresholdReport::InvalidSupplied),
        };
    }
    if let Some(error_rate) = error_rate {
        let computed = (1.0 / error_rate - 2.0 * k as f64) * error_rate;
        if 0.0 < computed && computed < 1.0 {
            return ThresholdChoice {
                value: computed,
                report: Some(ThresholdReport::Computed(computed)),
            };
        }
        return ThresholdChoice {
            value: 0.8,
            report: Some(ThresholdReport::InvalidComputed(computed)),
        };
    }
    ThresholdChoice {
        value: 0.8,
        report: None,
    }
}

fn report_threshold_choice(report: &ThresholdReport) {
    match report {
        ThresholdReport::Computed(threshold) => {
            eprintln!("Using computed threshold {threshold}");
        }
        ThresholdReport::InvalidSupplied => {
            eprintln!(
                "Threshold not in (0,1). Setting default threshold for unmapped kmers to 0.8"
            );
        }
        ThresholdReport::InvalidComputed(threshold) => {
            eprintln!(
                "Supplied and computed threshold are invalid, using default value of 0.8 (computed: {threshold})"
            );
        }
    }
}

fn bus_pseudoalign_options(
    args: &BusArgs,
    strand_specific: Option<kallistors::pseudoalign::StrandSpecific>,
) -> kallistors::pseudoalign::PseudoalignOptions {
    kallistors::pseudoalign::PseudoalignOptions {
        kallisto_enum: args.kallisto_enum,
        kallisto_strict: args.kallisto_strict,
        kallisto_local_fallback: args.kallisto_local_fallback,
        kallisto_fallback: args.kallisto_fallback,
        discard_special_only: args.discard_special_only,
        skip_overcrowded_minimizer: args.skip_overcrowded_minimizer,
        kallisto_direct_kmer: args.kallisto_direct_kmer,
        kallisto_bifrost_find: args.kallisto_bifrost_find,
        kallisto_sparse_hits: args.kallisto_sparse_hits,
        strand_specific,
        do_union: args.do_union,
        no_jump: args.no_jump,
        dfk_onlist: args.dfk_onlist || args.aa,
        investigation: super::investigation_options_from_env(),
        ..kallistors::pseudoalign::PseudoalignOptions::default()
    }
}

fn bus_strand_specific(
    args: &BusArgs,
    spec: &TechnologySpec,
    index: &kallistors::pseudoalign::BifrostIndex,
) -> Option<kallistors::pseudoalign::StrandSpecific> {
    if args.unstranded {
        None
    } else if args.fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if args.rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else if index.has_strand_annotations() {
        spec.default_strand
    } else {
        None
    }
}

fn user_specified_strand(args: &BusArgs) -> bool {
    args.unstranded || args.fr_stranded || args.rf_stranded
}

fn report_default_strand(strand: kallistors::pseudoalign::StrandSpecific) {
    let option = match strand {
        kallistors::pseudoalign::StrandSpecific::FR => "--fr-stranded",
        kallistors::pseudoalign::StrandSpecific::RF => "--rf-stranded",
    };
    eprintln!(
        "[bus] Note: Strand option was not specified; setting it to {option} for specified technology"
    );
}

fn report_missing_strand_annotations() {
    eprintln!(
        "[bus] Note: Strand option was not specified; index has no strand annotations, processing as --unstranded"
    );
}

pub fn run(args: BusArgs) -> Result<()> {
    if args.list {
        print_technology_list();
        return Ok(());
    }
    let technology = args.technology.as_deref().unwrap_or("Bulk");
    let index_path = args
        .index
        .as_deref()
        .ok_or_else(|| anyhow!("kallisto index file missing"))?;
    if !index_path.exists() {
        bail!("kallisto index file not found {}", index_path.display());
    }
    let out_dir = args
        .out_dir
        .as_deref()
        .ok_or_else(|| anyhow!("need to specify output directory"))?;
    if out_dir.exists() && !out_dir.is_dir() {
        bail!("file {} exists and is not a directory", out_dir.display());
    }
    if args.batch.is_none() && args.reads.is_empty() {
        if args.technology.is_none() {
            bail!(
                "the technology must be specified via -x, use \"bulk\" for regular RNA-seq reads"
            );
        }
        bail!("missing read files");
    }
    if args.threads == 0 {
        bail!("invalid number of threads 0");
    }
    let write_pseudobam = args.pseudobam || args.genomebam;
    if args.pseudobam && args.genomebam {
        bail!("--pseudobam and --genomebam are mutually exclusive");
    }
    if args.genomebam && args.gtf.is_none() {
        bail!("--genomebam requires --gtf");
    }
    if args.gtf.is_some() && !args.genomebam {
        bail!("--gtf is only valid with --genomebam");
    }
    if args.chromosomes.is_some() && !args.genomebam {
        bail!("--chromosomes is only valid with --genomebam");
    }
    if let Some(gtf) = args.gtf.as_deref()
        && !gtf.exists()
    {
        bail!("GTF file not found: {}", gtf.display());
    }
    if let Some(chromosomes) = args.chromosomes.as_deref()
        && !chromosomes.exists()
    {
        bail!("Chromosome file not found: {}", chromosomes.display());
    }
    if args.batch.is_some() && !args.reads.is_empty() {
        bail!("cannot specify batch mode and supply read files");
    }
    if let Some(batch) = args.batch.as_deref()
        && !batch.exists()
    {
        bail!("file not found {}", batch.display());
    }
    if args.batch.is_some() {
        eprintln!("[bus] will try running read files supplied in batch file");
        if args.technology.is_none() && args.paired {
            eprintln!(
                "[bus] --paired ignored; single/paired-end is inferred from number of files supplied"
            );
        }
    }
    if args.bam && args.batch.is_some() {
        bail!("--bam cannot be combined with batch mode");
    }
    if args.bam && args.interleaved {
        bail!("--bam cannot be combined with interleaved FASTQ input");
    }
    if args.bam && args.reads.len() != 1 {
        bail!("--bam expects exactly one BAM file");
    }
    if args.bam && args.num {
        eprintln!("Warning: --bam option was used, so --num option will be ignored");
    }
    if args.bam && args.batch_barcodes {
        bail!("--batch-barcodes requires batch mode and cannot be used with --bam");
    }
    if args.bam && args.paired && !args.long && !args.aa {
        bail!("Paired reads are not compatible with the specified technology");
    }
    if args.batch_barcodes && args.batch.is_none() {
        bail!("--batch-barcodes requires batch mode");
    }
    if args.interleaved && args.batch.is_some() {
        bail!("interleaved input cannot be specified with a batch file");
    }
    if args.interleaved && args.reads.len() != 1 {
        bail!("interleaved input expects exactly one FASTQ file");
    }
    if args.fr_stranded && args.rf_stranded {
        bail!("--fr-stranded and --rf-stranded are mutually exclusive");
    }
    if args.unstranded && (args.fr_stranded || args.rf_stranded) {
        bail!("--unstranded cannot be combined with --fr-stranded or --rf-stranded");
    }
    if args.long
        && let Some(error_rate) = args.error_rate
        && error_rate <= 0.0
    {
        bail!("--error-rate must be greater than 0");
    }
    if args.long
        && let Some(platform) = args.platform.as_deref()
        && !platform.eq_ignore_ascii_case("PACBIO")
        && !platform.eq_ignore_ascii_case("ONT")
    {
        bail!("--platform must be PACBIO or ONT");
    }
    if (args.do_union || args.no_jump) && (args.long || args.aa) {
        bail!("--union and --no-jump are not compatible with --long or --aa");
    }
    let mut spec = technology_spec(technology)?;
    let technology_base = technology_base_upper(technology);
    if technology_base == "BULK" {
        if write_pseudobam {
            bail!("Pseudobam not supported yet in this mode");
        }
        if args.bam {
            bail!("--bam not supported in this mode");
        }
        if args.tag_sequence.is_some() {
            bail!("--tag not supported in this mode");
        }
    }
    let tag_config = configure_tag(technology, args.tag_sequence.as_deref(), &mut spec)?;
    if let Some(tag) = tag_config.as_ref()
        && tag.inferred_default
    {
        eprintln!(
            "[bus] Using {} as UMI tag sequence",
            String::from_utf8_lossy(&tag.config.sequence)
        );
    }
    if args.bam && tag_config.is_some() {
        bail!("--tag is not supported with --bam input");
    }
    if write_pseudobam && args.bam {
        bail!("--pseudobam/--genomebam is only supported for FASTQ BUS input");
    }
    if args.technology.is_none()
        && !args.aa
        && !args.long
        && !spec.paired
        && should_infer_bulk_paired(&args)?
    {
        mark_technology_paired(&mut spec, "Bulk input")?;
    }
    if args.paired && !args.aa && !args.long && !spec.paired {
        mark_technology_paired(&mut spec, "--paired")?;
    } else if args.paired && args.long && !args.aa && technology_base == "SMARTSEQ2" {
        add_unpaired_mate_slice(&mut spec, "--paired")?;
    }
    if args.aa && args.paired && !spec.paired {
        eprintln!("[bus] --paired ignored; --aa only supports single-end reads");
    }
    if args.aa && spec.paired {
        bail!("--aa BUS mode currently supports single cDNA read technologies");
    }
    if write_pseudobam && !pseudobam_compatible(&spec) {
        bail!(
            "BAM output is currently only supported for technologies with a single cDNA read file or paired-end reads"
        );
    }
    let read_groups = read_groups(&args, spec.nfiles)?;
    if read_groups.is_empty() {
        bail!("missing read files");
    }
    validate_read_files(&args, &read_groups)?;
    if args.verbose {
        report_verbose_read_groups(&read_groups);
    }
    let synthetic_bc = spec.bc.first().is_some_and(|v| v.file.is_none());
    let fixed_bc_len = total_len(&spec.bc);
    if let Some(fixed_bc_len) = fixed_bc_len
        && fixed_bc_len > 32
    {
        bail!("barcode length {fixed_bc_len} exceeds BUS limit of 32 bases");
    }
    let fixed_umi_len = total_len(&spec.umi);
    if let Some(fixed_umi_len) = fixed_umi_len
        && fixed_umi_len > 32
    {
        bail!("UMI length {fixed_umi_len} exceeds BUS limit of 32 bases");
    }
    let batch_barcode_prefix_len = if args.batch_barcodes {
        let Some(fixed_bc_len) = fixed_bc_len else {
            bail!("--batch-barcodes requires a bounded barcode length");
        };
        if fixed_bc_len >= 32 {
            bail!("--batch-barcodes requires barcode length shorter than 32 bases");
        }
        Some(32 - fixed_bc_len)
    } else {
        None
    };
    if let Some(prefix_len) = batch_barcode_prefix_len
        && prefix_len == 0
    {
        bail!("--batch-barcodes requires barcode length shorter than 32 bases");
    }
    let bc_len = if args.batch_barcodes && fixed_bc_len.is_some_and(|v| v > 0) {
        32
    } else if synthetic_bc {
        16
    } else {
        fixed_bc_len.unwrap_or(0)
    };
    if !args.bam
        && args.batch.is_none()
        && !args.interleaved
        && !args.reads.len().is_multiple_of(spec.nfiles)
    {
        bail!(
            "technology {} expects files in groups of {}, got {}",
            technology,
            spec.nfiles,
            args.reads.len()
        );
    }

    let start = Instant::now();
    let index = if args.kallisto_direct_kmer {
        kallistors::pseudoalign::build_bifrost_index_with_kmer_pos(index_path, false)
    } else if args.kallisto_fallback {
        kallistors::pseudoalign::build_bifrost_index_with_kmer(index_path)
    } else {
        kallistors::pseudoalign::build_bifrost_index_with_positions_threaded(
            index_path,
            false,
            args.threads.max(1),
        )
    }
    .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
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
    let mut long_len_sums = vec![0u64; index.transcript_names.len()];
    let mut long_len_counts = vec![0u64; index.transcript_names.len()];
    let mut paired_flens = vec![0u32; kallistors::pseudoalign::MAX_FRAG_LEN as usize];
    let batch_count = read_groups
        .iter()
        .map(|group| group.batch_index)
        .max()
        .map_or(0, |idx| idx + 1);
    let mut batch_long_len_sums = vec![vec![0u64; index.transcript_names.len()]; batch_count];
    let mut batch_long_len_counts = vec![vec![0u64; index.transcript_names.len()]; batch_count];
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
    write_transcripts(
        &out_dir.join("transcripts.txt"),
        &index.transcript_names,
        index.onlist.as_deref(),
    )?;
    write_run_info(
        &out_dir.join("run_info.json"),
        RunInfo {
            processed: records,
            aligned,
            unique,
            targets: onlist_target_count(&index.transcript_names, index.onlist.as_deref()),
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

fn process_bam_records(
    args: &BusArgs,
    index: &kallistors::pseudoalign::BifrostIndex,
    processing: BusProcessingConfig,
    state: &mut BusWriteState<'_, impl Write>,
    unmapped_ratios: &mut Vec<f64>,
    novel_writer: &mut Option<BufWriter<File>>,
) -> Result<(u64, u64, u64, u64)> {
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
        let quality = record
            .quality_scores()
            .iter()
            .map(|score| score.saturating_add(b'!'))
            .collect::<Vec<_>>();
        let quality = if quality.len() == sequence.len() {
            quality
        } else {
            vec![b'!'; sequence.len()]
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

fn bam_string_tag(record: &noodles_bam::Record, tag: [u8; 2]) -> Option<Vec<u8>> {
    match record.data().get(&tag)?.ok()? {
        SamDataValue::String(value) | SamDataValue::Hex(value) => {
            let bytes: &[u8] = value.as_ref();
            Some(bytes.to_vec())
        }
        _ => None,
    }
}

fn normalize_corrected_barcode(mut barcode: Vec<u8>) -> Vec<u8> {
    if let Some(pos) = barcode.iter().position(|base| *base == b'-') {
        barcode.truncate(pos);
    }
    barcode
}

fn print_technology_list() {
    println!("List of supported single-cell technologies");
    println!();
    println!("short name       description");
    println!("----------       -----------");
    println!("10xv1            10x version 1 chemistry");
    println!("10xv2            10x version 2 chemistry");
    println!("10xv3            10x version 3 chemistry");
    println!("10xv4            10x version 4 chemistry");
    println!("Bulk             Bulk RNA-seq");
    println!("ParseV3          Parse Evercode V3");
    println!("SmartSeq2        Smart-seq2 (multiplexed)");
    println!("BDWTA            BD Rhapsody WTA");
    println!("CELSeq           CEL-Seq");
    println!("CELSeq2          CEL-Seq version 2");
    println!("DropSeq          DropSeq");
    println!("inDropsv1        inDrops version 1 chemistry");
    println!("inDropsv2        inDrops version 2 chemistry");
    println!("inDropsv3        inDrops version 3 chemistry");
    println!("MATQSEQ          MATQ-SEQ");
    println!("PETRISEQ         PETRI-SEQ");
    println!("SCRBSeq          SCRB-Seq");
    println!("SmartSeq3        Smart-seq3");
    println!("SPLiT-seq        SPLiT-seq");
    println!("STORM-seq        STORM-seq");
    println!("SureCell         SureCell for ddSEQ");
    println!("VASA-seq         VASA-seq");
    println!("Visium           10x Visium Spatial Transcriptomics");
}

fn barcode_sequence(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    batch_index: usize,
    batch_barcodes: bool,
) -> Result<Option<Vec<u8>>> {
    let extracted = if spec.bc.first().is_some_and(|v| v.file.is_none()) {
        fake_barcode(batch_index as u64, 16)
    } else {
        let Some(barcode) = concat_slices_optional(batch, &spec.bc, false)? else {
            return Ok(None);
        };
        barcode
    };
    if batch_barcodes && spec.bc.first().is_none_or(|v| v.file.is_some()) {
        if extracted.len() >= 32 {
            bail!("--batch-barcodes requires barcode length shorter than 32 bases");
        }
        let prefix_len = 32usize.saturating_sub(extracted.len()).min(32);
        let mut prefixed = fake_barcode(batch_index as u64, prefix_len);
        prefixed.extend_from_slice(&extracted);
        Ok(Some(prefixed))
    } else {
        Ok(Some(extracted))
    }
}

fn configure_tag(
    technology: &str,
    tag_sequence: Option<&str>,
    spec: &mut TechnologySpec,
) -> Result<Option<ConfiguredTag>> {
    let inferred_default = tag_sequence.is_none() && technology.eq_ignore_ascii_case("SMARTSEQ3");
    let tag = tag_sequence
        .map(str::to_string)
        .or_else(|| inferred_default.then(|| "ATTGCGCAATG".to_string()));
    let Some(tag) = tag else {
        return Ok(None);
    };
    validate_tag_sequence(&tag)?;
    let Some(first_umi) = spec.umi.first_mut() else {
        bail!("tag sequence requires a UMI slice");
    };
    let Some(stop) = first_umi.stop else {
        bail!("tag sequence requires a bounded UMI slice");
    };
    first_umi.start = first_umi.start.saturating_add(tag.len());
    if first_umi.start >= stop {
        bail!("tag sequence must be shorter than UMI sequence");
    }
    spec.tag_len = Some(tag.len());
    let mut flag = 0;
    Ok(Some(ConfiguredTag {
        config: TagConfig {
            binary: string_to_binary(tag.as_bytes(), &mut flag, "tag sequence")?,
            sequence: tag.into_bytes(),
        },
        inferred_default,
    }))
}

fn validate_tag_sequence(tag: &str) -> Result<()> {
    if tag.is_empty() {
        bail!("tag sequence cannot be empty");
    }
    for (idx, base) in tag.bytes().enumerate() {
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            bail!(
                "tag sequence contains invalid base {} at position {}",
                char::from(base),
                idx + 1
            );
        }
    }
    Ok(())
}

fn umi_value(
    batch: &[kallistors::io::FastqRecord],
    spec: &TechnologySpec,
    tag: Option<&TagConfig>,
) -> Result<Option<UmiValue>> {
    if spec.bulk_like() {
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(b"A", &mut flag, "UMI")?,
            flag,
            len: 1,
            seq: b"A".to_vec(),
            tag_present: false,
            disable_strand_specific: false,
        }));
    }
    if spec.keep_fastq_comments {
        let Some(umi_seq) = rx_umi_from_header(
            batch
                .first()
                .ok_or_else(|| anyhow!("RX UMI extraction requires at least one FASTQ record"))?,
        ) else {
            return Ok(None);
        };
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(&umi_seq, &mut flag, "UMI")?,
            flag,
            len: umi_seq.len(),
            seq: umi_seq,
            tag_present: false,
            disable_strand_specific: false,
        }));
    }

    let Some(tag) = tag else {
        let Some(umi_seq) = concat_slices_optional(batch, &spec.umi, false)? else {
            return Ok(None);
        };
        let mut flag = 0;
        return Ok(Some(UmiValue {
            binary: string_to_binary(&umi_seq, &mut flag, "UMI")?,
            flag,
            len: umi_seq.len(),
            seq: umi_seq,
            tag_present: false,
            disable_strand_specific: false,
        }));
    };

    let Some(first_umi) = spec.umi.first() else {
        bail!("tag sequence requires a UMI slice");
    };
    let Some(mut tagged_umi) = extract_tagged_umi_optional(batch, *first_umi, tag.sequence.len())?
    else {
        return Ok(None);
    };
    for spec in spec.umi.iter().skip(1) {
        let Some(slice) = extract_slice_optional(batch, *spec)? else {
            return Ok(None);
        };
        tagged_umi.extend_from_slice(&slice);
    }
    let mut flag = 0;
    let tagged_binary = string_to_binary(&tagged_umi, &mut flag, "tagged UMI")?;
    let tag_len = tag.sequence.len();
    let hamming_threshold = usize::from(tag_len > 5);
    let tag_present = tagged_umi.len() >= tag_len
        && hamming(
            tag.binary,
            tagged_binary >> (2 * (tagged_umi.len() - tag_len)),
            tag_len,
        ) <= hamming_threshold;
    if tag_present {
        let umi_len = tagged_umi.len() - tag_len;
        let binary = if umi_len == 0 {
            0
        } else {
            tagged_binary & !(u64::MAX << (2 * umi_len))
        };
        Ok(Some(UmiValue {
            binary,
            flag,
            len: umi_len,
            seq: tagged_umi[tag_len..].to_vec(),
            tag_present: true,
            disable_strand_specific: false,
        }))
    } else {
        Ok(Some(UmiValue {
            binary: u64::MAX,
            flag: 0,
            len: 0,
            seq: Vec::new(),
            tag_present: false,
            disable_strand_specific: true,
        }))
    }
}

fn extract_tagged_umi_optional(
    records: &[kallistors::io::FastqRecord],
    spec: SliceSpec,
    tag_len: usize,
) -> Result<Option<Vec<u8>>> {
    let mut tagged_spec = spec;
    tagged_spec.start = tagged_spec.start.saturating_sub(tag_len);
    extract_slice_optional(records, tagged_spec)
}

fn rx_umi_from_header(record: &kallistors::io::FastqRecord) -> Option<Vec<u8>> {
    let header = record.header.strip_prefix(b"@").unwrap_or(&record.header);
    let pos = header.windows(5).position(|window| window == b"RX:Z:")?;
    let start = pos + 5;
    let end = header[start..]
        .iter()
        .position(|base| matches!(base, b' ' | b'\t'))
        .map_or(header.len(), |offset| start + offset);
    (end > start).then(|| header[start..end].to_vec())
}

fn record_unmapped_ratio_for_skipped_read(
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

fn record_bam_outputs_for_skipped_read(
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

fn pseudoalign_bus_read(
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
            ec: trace.merged_ec,
            fragment_length: trace.estimated_fragment_length,
            placement: trace.placement,
            mate_placements: vec![trace.left.placement, trace.right.placement],
            mate_has_ec: vec![
                trace.left.ec_after_strand_filter.is_some(),
                trace.right.ec_after_strand_filter.is_some(),
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
            ec: trace.ec_after_strand_filter,
            fragment_length: None,
            placement: trace.placement,
            mate_placements: Vec::new(),
            mate_has_ec: Vec::new(),
            frame_clashes: 0,
        })
    }
}

fn aa_alignment_for_query(
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

fn bus_query_sequence(
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

fn bus_query_quality(
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

fn onlist_target_count(names: &[String], onlist: Option<&[bool]>) -> usize {
    match onlist {
        Some(onlist) => names
            .iter()
            .enumerate()
            .filter(|(idx, _)| onlist.get(*idx).copied().unwrap_or(false))
            .count(),
        None => names.len(),
    }
}

fn string_to_binary(seq: &[u8], flag: &mut u32, label: &str) -> Result<u64> {
    if seq.len() > 32 {
        bail!("{label} length {} exceeds BUS limit of 32 bases", seq.len());
    }
    *flag = 0;
    let mut result = 0u64;
    let mut num_n = 0u32;
    let mut pos_n = 0u32;
    for (idx, base) in seq.iter().enumerate() {
        let x = (base & 4) >> 1;
        if (base & 3) == 2 {
            if num_n == 0 {
                pos_n = u32::try_from(idx).unwrap_or(0);
            }
            num_n += 1;
        }
        result <<= 2;
        result |= u64::from(x + ((x ^ (base & 2)) >> 1));
    }
    if num_n > 0 {
        *flag = (num_n.min(3) & 3) | ((pos_n & 31) << 2);
    }
    Ok(result)
}

fn hamming(left: u64, right: u64, len: usize) -> usize {
    let diff = left ^ right;
    let mut distance = 0usize;
    for idx in 0..len {
        let shift = 2 * (len - idx - 1);
        if ((diff >> shift) & 0x03) != 0 {
            distance += 1;
        }
    }
    distance
}

#[cfg(test)]
mod tests;
