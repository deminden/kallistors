use anyhow::{Result, anyhow};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

pub struct QuantArgs {
    pub index: PathBuf,
    pub out_dir: PathBuf,
    pub reads: Vec<PathBuf>,
    pub single: bool,
    pub fragment_length: Option<f64>,
    pub fragment_length_sd: Option<f64>,
    pub single_overhang: bool,
    pub fr_stranded: bool,
    pub rf_stranded: bool,
    pub bias: bool,
    pub transcripts: Option<PathBuf>,
    pub seed: u64,
    pub threads: usize,
    pub pseudobam: bool,
    pub genomebam: bool,
}

pub fn run(args: QuantArgs) -> Result<()> {
    let _ = args.seed;
    if args.fr_stranded && args.rf_stranded {
        return Err(anyhow!(
            "--fr-stranded and --rf-stranded are mutually exclusive"
        ));
    }
    if args.bias && args.transcripts.is_none() {
        return Err(anyhow!("--bias requires --transcripts"));
    }
    if args.pseudobam {
        return Err(anyhow!("--pseudobam is not supported yet"));
    }
    if args.genomebam {
        return Err(anyhow!("--genomebam is not supported yet"));
    }
    if args.single_overhang && !args.single {
        return Err(anyhow!("--single-overhang requires --single"));
    }

    let start_time = format_start_time();
    let strand_specific = if args.fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if args.rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else {
        None
    };

    let (single, reads1, reads2) = resolve_reads(&args)?;
    let fragment_length = if single {
        args.fragment_length
            .ok_or_else(|| anyhow!("--single requires -l/--fragment-length"))?
    } else if args.fragment_length.is_some() {
        return Err(anyhow!("-l/--fragment-length is only valid with --single"));
    } else {
        0.0
    };
    let fragment_length_sd = if single {
        args.fragment_length_sd
            .ok_or_else(|| anyhow!("--single requires -s/--fragment-length-sd"))?
    } else if args.fragment_length_sd.is_some() {
        return Err(anyhow!(
            "-s/--fragment-length-sd is only valid with --single"
        ));
    } else {
        0.0
    };

    let filter = if single && fragment_length > 0.0 {
        Some(kallistors::pseudoalign::FragmentFilter {
            fragment_length: fragment_length as u32,
            single_overhang: args.single_overhang,
        })
    } else {
        None
    };
    let load_positional_info = filter
        .map(|v| !v.single_overhang && v.fragment_length > 0)
        .unwrap_or(false);
    let index = if load_positional_info || !single {
        kallistors::pseudoalign::build_bifrost_index_with_positions(&args.index, true)
    } else {
        kallistors::pseudoalign::build_bifrost_index(&args.index)
    }
    .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;

    let mut reader1 = kallistors::io::open_fastq_reader(reads1)?;
    let mut reader2 = if let Some(reads2) = reads2 {
        Some(kallistors::io::open_fastq_reader(reads2)?)
    } else {
        None
    };

    let options = kallistors::pseudoalign::PseudoalignOptions {
        min_range: 1,
        do_union: false,
        dfk_onlist: false,
        strand_specific,
        no_jump: false,
        bias: args.bias,
        max_bias: 1_000_000,
    };

    let start = Instant::now();
    let threads = args.threads.max(1);
    let ec_counts = if let Some(reader2) = reader2.as_mut() {
        kallistors::pseudoalign::pseudoalign_paired_bifrost_with_options_threaded(
            &index,
            &mut reader1,
            reader2,
            kallistors::pseudoalign::Strand::Unstranded,
            options,
            threads,
        )
        .map_err(|err| anyhow!("pseudoalign failed: {err}"))?
    } else {
        kallistors::pseudoalign::pseudoalign_single_end_bifrost_with_options_threaded(
            &index,
            &mut reader1,
            kallistors::pseudoalign::Strand::Unstranded,
            filter,
            options,
            threads,
        )
        .map_err(|err| anyhow!("pseudoalign failed: {err}"))?
    };

    let index_meta = kallistors::index::Index::load(&args.index)?;
    let lengths = index_meta
        .transcripts
        .iter()
        .map(|t| t.length)
        .collect::<Vec<_>>();

    let transcript_seqs = if args.bias {
        Some(kallistors::quant::load_transcript_sequences(
            args.transcripts.as_ref().unwrap(),
            &index_meta
                .transcripts
                .iter()
                .map(|t| t.name.clone())
                .collect::<Vec<_>>(),
        )?)
    } else {
        None
    };

    let (mean_fragment_length, fragment_length_sd) = if single {
        (fragment_length, fragment_length_sd)
    } else {
        estimate_paired_fragment_lengths(&ec_counts).unwrap_or((200.0, 20.0))
    };

    let result = kallistors::quant::em_quantify(
        &kallistors::quant::EcCountsInput {
            ec_list: kallistors::ec::EcList {
                classes: ec_counts.ec_list.clone(),
            },
            counts: ec_counts.counts.clone(),
        },
        &lengths,
        transcript_seqs.as_deref(),
        ec_counts.bias.as_ref(),
        kallistors::quant::QuantOptions {
            mean_fragment_length,
            fragment_length_sd,
            bias: args.bias,
            strand_specific,
            ..kallistors::quant::QuantOptions::default()
        },
    )?;

    fs::create_dir_all(&args.out_dir)?;
    let abundance_path = args.out_dir.join("abundance.tsv");
    let run_info_path = args.out_dir.join("run_info.json");
    kallistors::quant::write_abundance_tsv(&abundance_path, &index_meta, &result)?;

    let unique = kallistors::pseudoalign::unique_pseudoaligned_reads(&ec_counts);
    let p_pseudoaligned = percent(ec_counts.reads_aligned, ec_counts.reads_processed);
    let p_unique = percent(unique, ec_counts.reads_processed);
    let run_info = kallistors::quant::RunInfo {
        n_targets: index_meta.transcripts.len(),
        n_bootstraps: 0,
        n_processed: ec_counts.reads_processed,
        n_pseudoaligned: ec_counts.reads_aligned,
        n_unique: unique,
        p_pseudoaligned,
        p_unique,
        kallisto_version: format!("kallistors {}", env!("CARGO_PKG_VERSION")),
        index_version: index_meta.index_version,
        kmer_length: index_meta.k,
        start_time,
        call: std::env::args().collect::<Vec<_>>().join(" "),
    };
    run_info.write_json(&run_info_path)?;
    let _ = ec_counts.fragment_length_hist;

    println!(
        "processed: {}, aligned: {}, ecs: {}, time: {:.2?}",
        ec_counts.reads_processed,
        ec_counts.reads_aligned,
        ec_counts.ec_list.len(),
        start.elapsed()
    );
    Ok(())
}

fn resolve_reads(args: &QuantArgs) -> Result<(bool, &Path, Option<&Path>)> {
    if args.single && args.reads.len() != 1 {
        return Err(anyhow!("--single expects exactly one reads file"));
    }
    if !args.single && args.reads.len() != 2 {
        return Err(anyhow!("paired-end quant expects two reads files"));
    }
    let reads1 = args.reads.first().ok_or_else(|| anyhow!("missing reads"))?;
    let reads2 = args.reads.get(1).map(|p| p.as_path());
    Ok((args.single, reads1.as_path(), reads2))
}

fn percent(numer: u64, denom: u64) -> f64 {
    if denom == 0 {
        0.0
    } else {
        numer as f64 * 100.0 / denom as f64
    }
}

fn estimate_paired_fragment_lengths(
    ec_counts: &kallistors::pseudoalign::EcCounts,
) -> Option<(f64, f64)> {
    let stats = ec_counts.fragment_length_stats.as_ref()?;
    let mean = stats.mean()?;
    let sd = stats.sd().unwrap_or(0.0);
    Some((mean, sd))
}

fn format_start_time() -> String {
    let format = time::format_description::parse(
        "[weekday repr:short] [month repr:short] [day padding:space] [hour]:[minute]:[second] [year]",
    )
    .unwrap_or_else(|_| Vec::new());
    match time::OffsetDateTime::now_local() {
        Ok(dt) => dt.format(&format).unwrap_or_else(|_| dt.to_string()),
        Err(_) => {
            let dt = time::OffsetDateTime::now_utc();
            dt.format(&format).unwrap_or_else(|_| dt.to_string())
        }
    }
}
