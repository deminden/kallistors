//! BUS option validation, technology setup and index selection.

use super::*;

pub(super) struct ThresholdChoice {
    pub(super) value: f64,
    pub(super) report: Option<ThresholdReport>,
}

pub(super) enum ThresholdReport {
    Computed(f64),
    InvalidSupplied,
    InvalidComputed(f64),
}

pub(super) fn long_read_threshold(
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

pub(super) fn report_threshold_choice(report: &ThresholdReport) {
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

pub(super) fn bus_pseudoalign_options(
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
        investigation: super::super::investigation_options_from_env(),
        ..kallistors::pseudoalign::PseudoalignOptions::default()
    }
}

pub(super) fn bus_strand_specific(
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

pub(super) fn user_specified_strand(args: &BusArgs) -> bool {
    args.unstranded || args.fr_stranded || args.rf_stranded
}

pub(super) fn report_default_strand(strand: kallistors::pseudoalign::StrandSpecific) {
    let option = match strand {
        kallistors::pseudoalign::StrandSpecific::FR => "--fr-stranded",
        kallistors::pseudoalign::StrandSpecific::RF => "--rf-stranded",
    };
    eprintln!(
        "[bus] Note: Strand option was not specified; setting it to {option} for specified technology"
    );
}

pub(super) fn report_missing_strand_annotations() {
    eprintln!(
        "[bus] Note: Strand option was not specified; index has no strand annotations, processing as --unstranded"
    );
}

pub(super) struct BusConfiguration<'a> {
    pub(super) index_path: &'a std::path::Path,
    pub(super) out_dir: &'a std::path::Path,
    pub(super) spec: TechnologySpec,
    pub(super) tag_config: Option<ConfiguredTag>,
    pub(super) read_groups: Vec<input::ReadGroup>,
    pub(super) write_pseudobam: bool,
    pub(super) synthetic_bc: bool,
    pub(super) fixed_bc_len: Option<u32>,
    pub(super) fixed_umi_len: Option<u32>,
    pub(super) batch_barcode_prefix_len: Option<u32>,
    pub(super) bc_len: u32,
}

pub(super) fn prepare(args: &BusArgs) -> Result<BusConfiguration<'_>> {
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
        && should_infer_bulk_paired(args)?
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
    let read_groups = read_groups(args, spec.nfiles)?;
    if read_groups.is_empty() {
        bail!("missing read files");
    }
    validate_read_files(args, &read_groups)?;
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

    Ok(BusConfiguration {
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
    })
}

pub(super) fn load_index(
    args: &BusArgs,
    index_path: &std::path::Path,
) -> Result<kallistors::pseudoalign::BifrostIndex> {
    // Strand filtering needs transcript orientation even without fragment positions.
    let load_strand_annotations = args.fr_stranded || args.rf_stranded;
    let strand = if args.fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if args.rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else {
        None
    };
    let options = bus_pseudoalign_options(args, strand);
    let mut index = if args.kallisto_direct_kmer {
        kallistors::pseudoalign::build_bifrost_index_with_kmer_pos(
            index_path,
            load_strand_annotations,
        )
    } else if args.kallisto_fallback {
        kallistors::pseudoalign::build_bifrost_index_with_positions_and_kmer(
            index_path,
            load_strand_annotations,
        )
    } else if counting::supports(args) {
        super::super::build_quant_index(index_path, false, args.threads.max(1), options).and_then(
            |index| {
                if index.use_shade {
                    // Shade filtering needs the full annotation representation.
                    drop(index);
                    kallistors::pseudoalign::build_bifrost_index_with_positions_threaded(
                        index_path,
                        load_strand_annotations,
                        args.threads.max(1),
                    )
                } else {
                    Ok(index)
                }
            },
        )
    } else {
        kallistors::pseudoalign::build_bifrost_index_with_positions_threaded(
            index_path,
            load_strand_annotations,
            args.threads.max(1),
        )
    }
    .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
    // All ordinary counting layouts share the exact scanner, including
    // normalized FASTQ and BAM. Input layout does not change selector semantics.
    if counting::supports(args) {
        index.prepare_read_mapping(options, args.threads.max(1));
    }
    Ok(index)
}

pub(super) fn print_technology_list() {
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
