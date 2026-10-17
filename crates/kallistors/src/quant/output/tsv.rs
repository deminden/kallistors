//! Buffered abundance table output.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::super::QuantResult;
use crate::Result;
use crate::index::Index;

/// Write kallisto-compatible `abundance.tsv`.
pub fn write_abundance_tsv(path: &Path, index: &Index, result: &QuantResult) -> Result<()> {
    let names = index
        .transcripts
        .iter()
        .map(|tr| tr.name.as_str())
        .collect::<Vec<_>>();
    let lengths = index
        .transcripts
        .iter()
        .map(|tr| tr.length)
        .collect::<Vec<_>>();
    write_abundance_tsv_from_parts(path, &names, &lengths, result)
}

/// Write kallisto-compatible `abundance.tsv` from already loaded transcript metadata.
pub fn write_abundance_tsv_from_parts(
    path: &Path,
    target_names: &[&str],
    target_lengths: &[u32],
    result: &QuantResult,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(1 << 20, file);
    writeln!(writer, "target_id\tlength\teff_length\test_counts\ttpm")?;
    for (i, name) in target_names.iter().enumerate() {
        let length = target_lengths.get(i).copied().unwrap_or(0);
        let eff_len = result.eff_lengths.get(i).copied().unwrap_or(0.0);
        let est = result.est_counts.get(i).copied().unwrap_or(0.0);
        let tpm = result.tpm.get(i).copied().unwrap_or(0.0);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            name, length, eff_len, est, tpm
        )?;
    }
    writer.flush()?;
    Ok(())
}
