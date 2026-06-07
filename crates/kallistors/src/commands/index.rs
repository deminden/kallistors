use std::path::PathBuf;

use anyhow::{Context, Result};

pub struct Args {
    pub output: PathBuf,
    pub k: usize,
    pub minimizer_len: Option<usize>,
    pub threads: usize,
    pub timings: bool,
    pub make_unique: bool,
    pub aa: bool,
    pub ec_max_size: i32,
    pub fasta: Vec<PathBuf>,
}

pub fn run(args: Args) -> Result<()> {
    let options = kallistors::index::IndexBuildOptions {
        k: args.k,
        g: args.minimizer_len,
        threads: args.threads,
        make_unique: args.make_unique,
        ec_max_size: args.ec_max_size,
        aa: args.aa,
    };
    let timings_enabled = args.timings || std::env::var_os("KALLISTORS_TIMINGS").is_some();
    if timings_enabled {
        let report =
            kallistors::index::build_index_with_report(args.output.as_path(), &args.fasta, options)
                .with_context(|| format!("failed to build index {}", args.output.display()))?;
        print_timings(&report);
        Ok(())
    } else {
        kallistors::index::build_index(args.output.as_path(), &args.fasta, options)
            .with_context(|| format!("failed to build index {}", args.output.display()))
    }
}

fn print_timings(report: &kallistors::index::IndexBuildReport) {
    eprintln!("[index-build]");
    eprintln!("transcripts\t{}", report.transcripts);
    eprintln!("unitigs\t{}", report.unitigs);
    eprintln!("kmers\t{}", report.kmers);
    eprintln!("minimizers\t{}", report.minimizers);
    eprintln!("[timings]");
    eprintln!("fasta_parse\t{:.6}", report.fasta_parse.as_secs_f64());
    eprintln!("graph_build\t{:.6}", report.graph_build.as_secs_f64());
    eprintln!("ec_build\t{:.6}", report.ec_build.as_secs_f64());
    eprintln!("minimizer_mphf\t{:.6}", report.minimizer_mphf.as_secs_f64());
    eprintln!("write\t{:.6}", report.write.as_secs_f64());
    eprintln!("total\t{:.6}", report.total.as_secs_f64());
}
