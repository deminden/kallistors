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
        pseudobam: bool,
        #[arg(long)]
        genomebam: bool,
        #[arg(required = true, num_args = 1..=2)]
        reads: Vec<std::path::PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Version => commands::version::run(),
        Commands::IndexInfo { index } => commands::index_info::run(&index),
        Commands::EcStats { ec } => commands::ec_stats::run(&ec),
        Commands::EcExport { ec, out } => commands::ec_export::run(&ec, &out),
        Commands::EcCompare { a, b } => commands::ec_compare::run(&a, &b),
        Commands::EcFromIndex { index, out } => commands::ec_from_index::run(&index, &out),
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
            threads,
            pseudobam,
            genomebam,
        }),
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Strand {
    Unstranded,
    Forward,
    Reverse,
}
