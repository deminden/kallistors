//! Ordered abundance tables with bounded parallel formatting.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rayon::prelude::*;

use super::super::QuantResult;
use crate::index::Index;
use crate::{Error, Result};

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

/// Write an abundance table from loaded transcript metadata.
pub fn write_abundance_tsv_from_parts(
    path: &Path,
    target_names: &[&str],
    target_lengths: &[u32],
    result: &QuantResult,
) -> Result<()> {
    write_abundance_tsv_from_parts_threaded(path, target_names, target_lengths, result, 1)
}

/// Format independent row chunks concurrently, writing them in transcript order.
/// At most two chunks per worker are retained between formatting and writing.
pub fn write_abundance_tsv_from_parts_threaded(
    path: &Path,
    target_names: &[&str],
    target_lengths: &[u32],
    result: &QuantResult,
    threads: usize,
) -> Result<()> {
    const CHUNK: usize = 4096;
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(1 << 20, file);
    writeln!(writer, "target_id\tlength\teff_length\test_counts\ttpm")?;
    if threads > 1 && target_names.len() >= CHUNK * 4 {
        let threads = threads.min(target_names.len().div_ceil(CHUNK));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|error| {
                Error::InvalidFormat(format!("failed to build output thread pool: {error}"))
            })?;
        let wave_rows = CHUNK.saturating_mul(threads).saturating_mul(2);
        for (wave, names) in target_names.chunks(wave_rows).enumerate() {
            let chunks: std::io::Result<Vec<Vec<u8>>> = pool.install(|| {
                names
                    .par_chunks(CHUNK)
                    .enumerate()
                    .map(|(chunk, names)| {
                        let mut buffer = Vec::with_capacity(names.len() * 128);
                        write_rows(
                            &mut buffer,
                            wave * wave_rows + chunk * CHUNK,
                            names,
                            target_lengths,
                            result,
                        )?;
                        Ok(buffer)
                    })
                    .collect()
            });
            for chunk in chunks? {
                writer.write_all(&chunk)?;
            }
        }
    } else {
        write_rows(&mut writer, 0, target_names, target_lengths, result)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_rows<W: Write>(
    writer: &mut W,
    start: usize,
    target_names: &[&str],
    target_lengths: &[u32],
    result: &QuantResult,
) -> std::io::Result<()> {
    for (offset, name) in target_names.iter().enumerate() {
        let i = start + offset;
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
    Ok(())
}
