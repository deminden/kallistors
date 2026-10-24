use anyhow::{Result, anyhow};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

pub struct QuantArgs {
    pub index: PathBuf,
    pub out_dir: PathBuf,
    pub reads: Vec<PathBuf>,
    pub timings: bool,
    pub single: bool,
    pub fragment_length: Option<f64>,
    pub fragment_length_samples: usize,
    pub fragment_length_sd: Option<f64>,
    pub single_overhang: bool,
    pub fr_stranded: bool,
    pub rf_stranded: bool,
    pub bias: bool,
    pub transcripts: Option<PathBuf>,
    pub bootstrap_samples: usize,
    pub seed: u64,
    pub plaintext: bool,
    pub kallisto_enum: bool,
    pub kallisto_strict: bool,
    pub kallisto_local_fallback: bool,
    pub kallisto_fallback: bool,
    pub discard_special_only: bool,
    pub skip_overcrowded_minimizer: bool,
    pub kallisto_direct_kmer: bool,
    pub kallisto_bifrost_find: bool,
    pub kallisto_sparse_hits: bool,
    pub threads: usize,
    pub match_cache_size: kallistors::pseudoalign::MatchCacheSize,
    pub pseudobam: bool,
    pub genomebam: bool,
}

pub fn run(args: QuantArgs) -> Result<()> {
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

    let timings_enabled = args.timings || std::env::var_os("KALLISTORS_TIMINGS").is_some();
    kallistors::timing::set_enabled(timings_enabled);
    kallistors::timing::reset();
    // Whole-command boundaries complement the nested/overlapping diagnostic timers.
    // The benchmark harness timestamps these flushed stderr markers externally.
    let mark_phase = |phase: &str| {
        if timings_enabled {
            eprintln!("[quant-stage] {phase}");
        }
    };

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
    let options = kallistors::pseudoalign::PseudoalignOptions {
        min_range: 1,
        do_union: false,
        dfk_onlist: false,
        strand_specific,
        no_jump: false,
        kallisto_enum: args.kallisto_enum,
        kallisto_strict: args.kallisto_strict,
        kallisto_local_fallback: args.kallisto_local_fallback,
        kallisto_fallback: args.kallisto_fallback,
        discard_special_only: args.discard_special_only,
        skip_overcrowded_minimizer: args.skip_overcrowded_minimizer,
        kallisto_direct_kmer: args.kallisto_direct_kmer,
        kallisto_bifrost_find: args.kallisto_bifrost_find,
        kallisto_sparse_hits: args.kallisto_sparse_hits,
        bias: args.bias,
        max_bias: 1_000_000,
        fragment_length_samples: args.fragment_length_samples,
        match_cache_size: args.match_cache_size,
        investigation: super::investigation_options_from_env(),
    };

    let cache_allocation = args.match_cache_size.allocation(args.threads)?;
    eprintln!(
        "[quant] match cache: {} ({} bytes per worker; {} bytes total across {} workers)",
        args.match_cache_size,
        cache_allocation.bytes_per_worker(),
        cache_allocation.total_bytes(),
        cache_allocation.workers()
    );

    mark_phase("index_begin");
    let mut index =
        super::build_quant_index(&args.index, load_positional_info, args.threads, options)
            .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
    // Reuse the bounded exact candidate selectors for both read layouts.
    // Preparation is included in index loading; selected postings still receive
    // full verification with the same fragment and strand semantics.
    index.prepare_read_mapping(options, args.threads.max(1));
    mark_phase("index_end");

    let mut reader1 = kallistors::io::open_fastq_reader(reads1)?;
    let mut reader2 = if let Some(reads2) = reads2 {
        Some(kallistors::io::open_fastq_reader(reads2)?)
    } else {
        None
    };

    let start = Instant::now();
    let threads = args.threads.max(1);
    mark_phase("reads_begin");
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

    mark_phase("reads_end");
    let transcript_seqs = if args.bias {
        Some(kallistors::quant::load_transcript_sequences(
            args.transcripts.as_ref().unwrap(),
            &index.transcript_names,
        )?)
    } else {
        None
    };

    let (mean_fragment_length, fragment_length_sd) = if single {
        (fragment_length, fragment_length_sd)
    } else {
        estimate_paired_fragment_lengths(&ec_counts).unwrap_or((200.0, 20.0))
    };
    if !single && let Some(stats) = ec_counts.fragment_length_stats.as_ref() {
        eprintln!(
            "[quant] fragment length estimate: {} observations, mean {:.3}, sd {:.3}",
            stats.count(),
            mean_fragment_length,
            fragment_length_sd
        );
    }
    let unique = kallistors::pseudoalign::unique_pseudoaligned_reads(&ec_counts);
    let kallistors::pseudoalign::EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias,
        fragment_length_stats: _,
        fragment_length_hist,
    } = ec_counts;
    let ec_class_count = ec_list.len();
    let ec_input = kallistors::quant::EcCountsInput {
        ec_list: kallistors::ec::EcList { classes: ec_list },
        counts,
    };

    // Mapping is complete. Retain output/model metadata and release the graph,
    // minimizer tables and caches before allocating EM/output working arrays.
    let transcript_names = index.take_transcript_names();
    let transcript_lengths = std::mem::take(&mut index.transcript_lengths);
    let index_version = index.index_version;
    let k = index.k;
    index.retire();

    let quant_options = kallistors::quant::QuantOptions {
        mean_fragment_length,
        fragment_length_sd,
        bias: args.bias,
        strand_specific,
        ..kallistors::quant::QuantOptions::default()
    };
    let mut ec_input = Some(ec_input);
    mark_phase("em_begin");
    let result = {
        let _em_timing = kallistors::timing::scoped(kallistors::timing::Stage::Em);
        if args.bootstrap_samples == 0 {
            kallistors::quant::em_quantify_owned_with_fragment_histogram_threaded(
                ec_input.take().expect("owned EC input available"),
                &transcript_lengths,
                transcript_seqs.as_deref(),
                bias.as_ref(),
                if single {
                    None
                } else {
                    fragment_length_hist.as_deref()
                },
                quant_options,
                args.threads,
            )?
        } else {
            kallistors::quant::em_quantify_with_fragment_histogram_threaded(
                ec_input.as_ref().expect("bootstrap EC input retained"),
                &transcript_lengths,
                transcript_seqs.as_deref(),
                bias.as_ref(),
                if single {
                    None
                } else {
                    fragment_length_hist.as_deref()
                },
                quant_options,
                args.threads,
            )?
        }
    };
    mark_phase("em_end");

    let transcript_names = transcript_names.load()?;
    let output_timing = kallistors::timing::scoped(kallistors::timing::Stage::QuantOutput);
    fs::create_dir_all(&args.out_dir)?;
    let abundance_path = args.out_dir.join("abundance.tsv");
    let run_info_path = args.out_dir.join("run_info.json");
    let target_names = transcript_names
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let write_tsv = || {
        kallistors::quant::write_abundance_tsv_from_parts_threaded(
            &abundance_path,
            &target_names,
            &transcript_lengths,
            &result,
            if args.plaintext {
                args.threads
            } else {
                args.threads.saturating_sub(1).max(1)
            },
        )
    };

    let write_other_outputs = || -> Result<()> {
        let bootstrap_results = if args.bootstrap_samples == 0 {
            Vec::new()
        } else {
            kallistors::quant::bootstrap_quantify(
                ec_input.as_ref().expect("bootstrap EC input retained"),
                &result.eff_lengths,
                args.bootstrap_samples,
                args.seed,
                quant_options,
            )?
        };

        if args.bootstrap_samples > 0 && args.plaintext {
            for (idx, bootstrap) in bootstrap_results.iter().enumerate() {
                let path = args.out_dir.join(format!("bs_abundance_{idx}.tsv"));
                kallistors::quant::write_abundance_tsv_from_parts(
                    &path,
                    &target_names,
                    &transcript_lengths,
                    bootstrap,
                )?;
            }
        }

        if !args.plaintext {
            let h5_path = args.out_dir.join("abundance.h5");
            let bootstrap_counts = bootstrap_results
                .iter()
                .map(|res| res.est_counts.as_slice())
                .collect::<Vec<_>>();
            kallistors::quant::write_abundance_h5(
                &h5_path,
                kallistors::quant::H5QuantOutput {
                    target_names: &target_names,
                    target_lengths: &transcript_lengths,
                    result: &result,
                    bootstrap_est_counts: &bootstrap_counts,
                    reads_processed,
                    fragment_length_hist: fragment_length_hist.as_deref(),
                    bias_observed: bias.as_ref().map(|counts| counts.counts.as_slice()),
                    index_version,
                    kallisto_version: &format!("kallistors {}", env!("CARGO_PKG_VERSION")),
                    call: &std::env::args().collect::<Vec<_>>().join(" "),
                    start_time: &start_time,
                },
            )?;
        }

        Ok(())
    };
    if args.threads > 1 && !args.plaintext {
        std::thread::scope(|scope| -> Result<()> {
            let tsv = scope.spawn(write_tsv);
            let other = write_other_outputs();
            tsv.join()
                .map_err(|_| anyhow!("abundance TSV writer panicked"))??;
            other
        })?;
    } else {
        write_tsv()?;
        write_other_outputs()?;
    }

    let p_pseudoaligned = percent(reads_aligned, reads_processed);
    let p_unique = percent(unique, reads_processed);
    let run_info = kallistors::quant::RunInfo {
        n_targets: transcript_names.len(),
        n_bootstraps: args.bootstrap_samples,
        n_processed: reads_processed,
        n_pseudoaligned: reads_aligned,
        n_unique: unique,
        p_pseudoaligned,
        p_unique,
        kallisto_version: format!("kallistors {}", env!("CARGO_PKG_VERSION")),
        index_version,
        kmer_length: k as u32,
        start_time,
        call: std::env::args().collect::<Vec<_>>().join(" "),
    };
    run_info.write_json_with_match_cache(&run_info_path, cache_allocation)?;
    drop(output_timing);
    mark_phase("output_end");

    println!(
        "processed: {}, aligned: {}, ecs: {}, time: {:.2?}",
        reads_processed,
        reads_aligned,
        ec_class_count,
        start.elapsed()
    );
    if timings_enabled {
        eprintln!("[timings]");
        for timing in kallistors::timing::snapshot() {
            eprintln!(
                "{}\t{:.3}s",
                timing.stage.label(),
                timing.duration.as_secs_f64()
            );
        }
    }
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
    let format = time::format_description::parse_borrowed::<1>(
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
