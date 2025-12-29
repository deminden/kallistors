use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};

mod commands;

#[derive(Debug, Parser)]
#[command(
    name = "kallistors-cli",
    version,
    about = "kallistors command line interface"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Version,
    IndexInfo {
        #[arg(long)]
        index: std::path::PathBuf,
    },
    #[command(about = "DEBUG: dump index EC blocks in a kallisto-like format")]
    IndexEcDump {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
        #[arg(long, default_value_t = 1000)]
        limit: usize,
        #[arg(long)]
        unitig_list: Option<std::path::PathBuf>,
    },
    #[command(about = "DEBUG: dump unitig sequence windows for minimizer positions")]
    UnitigWindow {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        positions: std::path::PathBuf,
        #[arg(long)]
        minimizer: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
    },
    #[command(about = "DEBUG: inspect a kallisto EC text file (to be removed)")]
    EcStats {
        #[arg(long)]
        ec: std::path::PathBuf,
    },
    #[command(about = "DEBUG: rewrite EC text file in canonical order (to be removed)")]
    EcExport {
        #[arg(long)]
        ec: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
    },
    #[command(about = "DEBUG: compare two EC text files (to be removed)")]
    EcCompare {
        #[arg(long)]
        a: std::path::PathBuf,
        #[arg(long)]
        b: std::path::PathBuf,
    },
    #[command(about = "DEBUG: extract EC list from index (to be removed)")]
    EcFromIndex {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
    },
    Pseudoalign {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        reads: std::path::PathBuf,
        #[arg(long)]
        reads2: Option<std::path::PathBuf>,
        #[arg(long)]
        out: std::path::PathBuf,
        #[arg(long, default_value_t = 1)]
        threads: usize,
        #[arg(long)]
        debug_out: Option<std::path::PathBuf>,
        #[arg(long, default_value_t = 100)]
        debug_max: usize,
        #[arg(long, value_enum, default_value_t = Strand::Unstranded)]
        strand: Strand,
        #[arg(long)]
        fragment_length: Option<f64>,
        #[arg(long)]
        single_overhang: bool,
        #[arg(long, default_value_t = 1)]
        min_range: usize,
        #[arg(long)]
        union: bool,
        #[arg(long)]
        no_jump: bool,
        #[arg(long)]
        dfk_onlist: bool,
        #[arg(long)]
        fr_stranded: bool,
        #[arg(long)]
        rf_stranded: bool,
        #[arg(long)]
        bias: bool,
        #[arg(long)]
        bias_out: Option<std::path::PathBuf>,
        #[arg(long)]
        kallisto_enum: bool,
        #[arg(long)]
        kallisto_strict: bool,
        #[arg(long)]
        kallisto_local_fallback: bool,
        #[arg(long)]
        kallisto_fallback: bool,
        #[arg(long)]
        discard_special_only: bool,
        #[arg(long)]
        skip_overcrowded_minimizer: bool,
        #[arg(long)]
        kallisto_direct_kmer: bool,
        #[arg(long)]
        kallisto_bifrost_find: bool,
    },
    Quant {
        #[arg(short = 'i', long)]
        index: std::path::PathBuf,
        #[arg(short = 'o', long = "output-dir", alias = "out")]
        out_dir: std::path::PathBuf,
        #[arg(short = 't', long, default_value_t = 1)]
        threads: usize,
        #[arg(long)]
        single: bool,
        #[arg(short = 'l', long = "fragment-length")]
        fragment_length: Option<f64>,
        #[arg(short = 's', long = "fragment-length-sd")]
        fragment_length_sd: Option<f64>,
        #[arg(long)]
        single_overhang: bool,
        #[arg(long)]
        fr_stranded: bool,
        #[arg(long)]
        rf_stranded: bool,
        #[arg(long)]
        bias: bool,
        #[arg(long)]
        transcripts: Option<std::path::PathBuf>,
        #[arg(long, default_value_t = 42)]
        seed: u64,
        #[arg(long)]
        kallisto_enum: bool,
        #[arg(long)]
        kallisto_strict: bool,
        #[arg(long)]
        kallisto_local_fallback: bool,
        #[arg(long)]
        kallisto_fallback: bool,
        #[arg(long)]
        discard_special_only: bool,
        #[arg(long)]
        skip_overcrowded_minimizer: bool,
        #[arg(long)]
        kallisto_direct_kmer: bool,
        #[arg(long)]
        kallisto_bifrost_find: bool,
        #[arg(long)]
        pseudobam: bool,
        #[arg(long)]
        genomebam: bool,
        #[arg(required = true, num_args = 1..=2)]
        reads: Vec<std::path::PathBuf>,
    },
    #[command(about = "DEBUG: map transcript ids to names from index")]
    TranscriptLookup {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        ids: Option<std::path::PathBuf>,
        #[arg(long)]
        hits: Option<std::path::PathBuf>,
        #[arg(long)]
        out: Option<std::path::PathBuf>,
    },
    #[command(about = "DEBUG: trace per-read EC decisions (to be removed)")]
    TraceReads {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        reads: std::path::PathBuf,
        #[arg(long)]
        read_list: std::path::PathBuf,
        #[arg(long)]
        out: Option<std::path::PathBuf>,
        #[arg(long)]
        hits_out: Option<std::path::PathBuf>,
        #[arg(long)]
        hits_intersection_out: Option<std::path::PathBuf>,
        #[arg(long)]
        gene_map: Option<std::path::PathBuf>,
        #[arg(long)]
        minimizer_positions: Option<std::path::PathBuf>,
        #[arg(long)]
        minimizer_out: Option<std::path::PathBuf>,
        #[arg(long)]
        positions_visited_out: Option<std::path::PathBuf>,
        #[arg(long)]
        local_kmer_out: Option<std::path::PathBuf>,
        #[arg(long)]
        kallisto_enum: bool,
        #[arg(long)]
        kallisto_strict: bool,
        #[arg(long)]
        kallisto_local_fallback: bool,
        #[arg(long)]
        kallisto_fallback: bool,
        #[arg(long)]
        discard_special_only: bool,
        #[arg(long)]
        skip_overcrowded_minimizer: bool,
        #[arg(long)]
        kallisto_direct_kmer: bool,
        #[arg(long)]
        kallisto_bifrost_find: bool,
        #[arg(long, value_enum, default_value_t = Strand::Unstranded)]
        strand: Strand,
        #[arg(long)]
        fragment_length: Option<f64>,
        #[arg(long)]
        single_overhang: bool,
        #[arg(long, default_value_t = 1)]
        min_range: usize,
        #[arg(long)]
        union: bool,
        #[arg(long)]
        no_jump: bool,
        #[arg(long)]
        dfk_onlist: bool,
        #[arg(long)]
        fr_stranded: bool,
        #[arg(long)]
        rf_stranded: bool,
    },
    #[command(about = "DEBUG: lookup minimizer in MPH index")]
    MinimizerLookup {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        minimizer: String,
    },
    #[command(about = "DEBUG: scan minimizer bitmaps for a specific sequence")]
    MinimizerBitmapScan {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        minimizer: String,
        #[arg(long)]
        out: Option<std::path::PathBuf>,
        #[arg(long, default_value_t = 1000)]
        limit: usize,
    },
    #[command(about = "DEBUG: check MPH hits for minimizers in selected unitigs")]
    MinimizerUnitigMphfCheck {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        unitigs: std::path::PathBuf,
        #[arg(long)]
        out: Option<std::path::PathBuf>,
    },
    #[command(about = "DEBUG: dump MPHF final-hash key list")]
    MinimizerMphfKeys {
        #[arg(long)]
        index: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Version => commands::version::run(),
        Commands::IndexInfo { index } => commands::index_info::run(&index),
        Commands::IndexEcDump {
            index,
            out,
            limit,
            unitig_list,
        } => commands::index_ec_dump::run(&index, &out, limit, unitig_list.as_deref()),
        Commands::UnitigWindow {
            index,
            positions,
            minimizer,
            out,
        } => commands::unitig_window::run(&index, &positions, &minimizer, &out),
        Commands::EcStats { ec } => commands::ec_stats::run(&ec),
        Commands::EcExport { ec, out } => commands::ec_export::run(&ec, &out),
        Commands::EcCompare { a, b } => commands::ec_compare::run(&a, &b),
        Commands::EcFromIndex { index, out } => commands::ec_from_index::run(&index, &out),
        Commands::MinimizerLookup { index, minimizer } => {
            commands::minimizer_lookup::run(index, minimizer)
        }
        Commands::MinimizerBitmapScan {
            index,
            minimizer,
            out,
            limit,
        } => commands::minimizer_bitmap_scan::run(index, minimizer, out, limit),
        Commands::MinimizerUnitigMphfCheck {
            index,
            unitigs,
            out,
        } => commands::minimizer_unitig_mphf_check::run(index, unitigs, out),
        Commands::MinimizerMphfKeys { index, out } => {
            commands::minimizer_mphf_keys::run(index, out)
        }
        Commands::Pseudoalign {
            index,
            reads,
            reads2,
            out,
            threads,
            debug_out,
            debug_max,
            strand,
            fragment_length,
            single_overhang,
            min_range,
            union,
            no_jump,
            dfk_onlist,
            fr_stranded,
            rf_stranded,
            bias,
            bias_out,
            kallisto_enum,
            kallisto_strict,
            kallisto_local_fallback,
            kallisto_fallback,
            discard_special_only,
            skip_overcrowded_minimizer,
            kallisto_direct_kmer,
            kallisto_bifrost_find,
        } => commands::pseudoalign::run(
            &index,
            &reads,
            reads2.as_deref(),
            &out,
            threads,
            debug_out.as_deref(),
            debug_max,
            match strand {
                Strand::Unstranded => kallistors::pseudoalign::Strand::Unstranded,
                Strand::Forward => kallistors::pseudoalign::Strand::Forward,
                Strand::Reverse => kallistors::pseudoalign::Strand::Reverse,
            },
            fragment_length,
            single_overhang,
            min_range,
            union,
            no_jump,
            dfk_onlist,
            fr_stranded,
            rf_stranded,
            bias,
            bias_out.as_deref(),
            kallisto_enum,
            kallisto_strict,
            kallisto_local_fallback,
            kallisto_fallback,
            discard_special_only,
            skip_overcrowded_minimizer,
            kallisto_direct_kmer,
            kallisto_bifrost_find,
        ),
        Commands::Quant {
            index,
            out_dir,
            threads,
            single,
            bias,
            transcripts,
            fragment_length,
            fragment_length_sd,
            single_overhang,
            fr_stranded,
            rf_stranded,
            seed,
            kallisto_enum,
            kallisto_strict,
            kallisto_local_fallback,
            kallisto_fallback,
            discard_special_only,
            skip_overcrowded_minimizer,
            kallisto_direct_kmer,
            kallisto_bifrost_find,
            pseudobam,
            genomebam,
            reads,
        } => commands::quant::run(commands::quant::QuantArgs {
            index,
            out_dir,
            reads,
            single,
            fragment_length,
            fragment_length_sd,
            single_overhang,
            fr_stranded,
            rf_stranded,
            bias,
            transcripts,
            seed,
            kallisto_enum,
            kallisto_strict,
            kallisto_local_fallback,
            kallisto_fallback,
            discard_special_only,
            skip_overcrowded_minimizer,
            kallisto_direct_kmer,
            kallisto_bifrost_find,
            threads,
            pseudobam,
            genomebam,
        }),
        Commands::TranscriptLookup {
            index,
            ids,
            hits,
            out,
        } => commands::transcript_lookup::run(
            &index,
            ids.as_deref(),
            hits.as_deref(),
            out.as_deref(),
        ),
        Commands::TraceReads {
            index,
            reads,
            read_list,
            out,
            hits_out,
            hits_intersection_out,
            gene_map,
            minimizer_positions,
            minimizer_out,
            positions_visited_out,
            local_kmer_out,
            kallisto_enum,
            kallisto_strict,
            kallisto_local_fallback,
            kallisto_fallback,
            discard_special_only,
            skip_overcrowded_minimizer,
            kallisto_direct_kmer,
            kallisto_bifrost_find,
            strand,
            fragment_length,
            single_overhang,
            min_range,
            union,
            no_jump,
            dfk_onlist,
            fr_stranded,
            rf_stranded,
        } => commands::trace_reads::run(
            &index,
            &reads,
            &read_list,
            out.as_deref(),
            hits_out.as_deref(),
            hits_intersection_out.as_deref(),
            gene_map.as_deref(),
            minimizer_positions.as_deref(),
            minimizer_out.as_deref(),
            positions_visited_out.as_deref(),
            local_kmer_out.as_deref(),
            kallisto_enum,
            kallisto_strict,
            kallisto_local_fallback,
            kallisto_fallback,
            discard_special_only,
            skip_overcrowded_minimizer,
            kallisto_direct_kmer,
            kallisto_bifrost_find,
            match strand {
                Strand::Unstranded => kallistors::pseudoalign::Strand::Unstranded,
                Strand::Forward => kallistors::pseudoalign::Strand::Forward,
                Strand::Reverse => kallistors::pseudoalign::Strand::Reverse,
            },
            fragment_length,
            single_overhang,
            min_range,
            union,
            no_jump,
            dfk_onlist,
            fr_stranded,
            rf_stranded,
        ),
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Strand {
    Unstranded,
    Forward,
    Reverse,
}
