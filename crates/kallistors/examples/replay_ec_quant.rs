//! Diagnostic EM replay with fixed EC counts and a saved fragment histogram.
//! Inputs: EC TSV (contiguous IDs), one length per line, one histogram bin per
//! line, output TSV, threads, max iterations or `default`, optional min rounds.
//! No reads or index needed.

use kallistors::quant::{
    QuantOptions, em_quantify_with_fragment_histogram_threaded, load_ec_counts,
};
use std::{env, fs, io::Write, path::Path};

fn integers(path: &str) -> anyhow::Result<Vec<u32>> {
    fs::read_to_string(path)?
        .lines()
        .map(|line| Ok(line.parse()?))
        .collect()
}

fn main() -> anyhow::Result<()> {
    let args: Vec<_> = env::args().collect();
    anyhow::ensure!(
        args.len() == 7 || args.len() == 8,
        "Usage: replay_ec_quant EC.tsv LENGTHS HISTOGRAM OUTPUT.tsv THREADS MAX_ITER_OR_DEFAULT [MIN_ROUNDS]"
    );
    let input = load_ec_counts(Path::new(&args[1]))?;
    let lengths = integers(&args[2])?;
    let histogram = integers(&args[3])?;
    let mut options = QuantOptions::default();
    if args[6] != "default" {
        options.max_iter = args[6].parse()?;
    }
    if args.len() == 8 {
        options.min_rounds = args[7].parse()?;
    }
    let result = em_quantify_with_fragment_histogram_threaded(
        &input,
        &lengths,
        None,
        None,
        Some(&histogram),
        options,
        args[5].parse()?,
    )?;
    let mut file = std::io::BufWriter::new(fs::File::create(&args[4])?);
    writeln!(file, "target_id\teff_length\test_counts\ttpm")?;
    for (i, ((eff, count), tpm)) in result
        .eff_lengths
        .iter()
        .zip(&result.est_counts)
        .zip(&result.tpm)
        .enumerate()
    {
        writeln!(file, "{i}\t{eff:.17}\t{count:.17}\t{tpm:.17}")?;
    }
    Ok(())
}
