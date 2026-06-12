use anyhow::{Result, bail};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

pub(crate) struct BusRecord {
    pub(crate) barcode: u64,
    pub(crate) umi: u64,
    pub(crate) ec: i32,
    pub(crate) count: u32,
    pub(crate) flags: u32,
}

pub(crate) struct RunInfo {
    pub(crate) processed: u64,
    pub(crate) aligned: u64,
    pub(crate) unique: u64,
    pub(crate) targets: usize,
    pub(crate) k: usize,
    pub(crate) frame_clashes: Option<u64>,
}

pub(crate) fn write_bus_header(writer: &mut impl Write, bc_len: u32, umi_len: u32) -> Result<()> {
    writer.write_all(b"BUS\0")?;
    writer.write_all(&1u32.to_le_bytes())?;
    writer.write_all(&bc_len.to_le_bytes())?;
    writer.write_all(&umi_len.to_le_bytes())?;
    let text = b"BUS file produced by kallisto";
    writer.write_all(&(text.len() as u32).to_le_bytes())?;
    writer.write_all(text)?;
    Ok(())
}

pub(crate) fn write_bus_record(writer: &mut impl Write, record: &BusRecord) -> Result<()> {
    writer.write_all(&record.barcode.to_le_bytes())?;
    writer.write_all(&record.umi.to_le_bytes())?;
    writer.write_all(&record.ec.to_le_bytes())?;
    writer.write_all(&record.count.to_le_bytes())?;
    writer.write_all(&record.flags.to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;
    Ok(())
}

pub(crate) fn write_matrix_ec(path: &Path, ec_list: &[Vec<u32>]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (idx, ec) in ec_list.iter().enumerate() {
        write!(writer, "{idx}\t")?;
        for (tx_idx, transcript) in ec.iter().enumerate() {
            if tx_idx > 0 {
                writer.write_all(b",")?;
            }
            write!(writer, "{transcript}")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

pub(crate) fn write_transcripts(
    path: &Path,
    names: &[String],
    onlist: Option<&[bool]>,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (idx, name) in names.iter().enumerate() {
        if onlist.is_some_and(|onlist| !onlist.get(idx).copied().unwrap_or(false)) {
            continue;
        }
        writeln!(writer, "{name}")?;
    }
    Ok(())
}

pub(crate) fn write_cells<'a>(
    path: &Path,
    groups: impl IntoIterator<Item = (usize, Option<&'a str>)>,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (idx, id) in groups {
        let id = id.map_or_else(|| format!("batch{idx}"), ToString::to_string);
        writeln!(writer, "{id}")?;
    }
    Ok(())
}

pub(crate) fn write_sample_barcodes(
    path: &Path,
    batch_indices: impl IntoIterator<Item = u64>,
    len: u32,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    let len = usize::try_from(len).unwrap_or(usize::MAX);
    for batch_index in batch_indices {
        writeln!(
            writer,
            "{}",
            String::from_utf8_lossy(&fake_barcode(batch_index, len))
        )?;
    }
    Ok(())
}

pub(crate) fn write_unmapped_ratios(
    path: &Path,
    ratios: &[f64],
    trailing_comma: bool,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (idx, ratio) in ratios.iter().enumerate() {
        if idx > 0 {
            writer.write_all(b",")?;
        }
        write!(writer, "{ratio}")?;
    }
    if trailing_comma && !ratios.is_empty() {
        writer.write_all(b",")?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

pub(crate) fn write_novel_read(
    writer: &mut impl Write,
    label: &str,
    seq: &[u8],
    quality: &[u8],
) -> Result<()> {
    if quality.len() != seq.len() {
        bail!("novel read quality length does not match sequence length");
    }
    writeln!(writer, "@{label}")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(quality)?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub(crate) fn write_long_flens(
    path: &Path,
    transcript_lengths: &[u32],
    sums: &[u64],
    counts: &[u64],
    k: usize,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_long_flens_line(&mut writer, transcript_lengths, sums, counts, k)?;
    Ok(())
}

pub(crate) fn write_batch_long_flens(
    path: &Path,
    transcript_lengths: &[u32],
    batch_sums: &[Vec<u64>],
    batch_counts: &[Vec<u64>],
    k: usize,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (sums, counts) in batch_sums.iter().zip(batch_counts) {
        write_long_flens_line(&mut writer, transcript_lengths, sums, counts, k)?;
    }
    Ok(())
}

fn write_long_flens_line<W: Write>(
    writer: &mut W,
    transcript_lengths: &[u32],
    sums: &[u64],
    counts: &[u64],
    k: usize,
) -> Result<()> {
    for (idx, length) in transcript_lengths.iter().enumerate() {
        if idx > 0 {
            writer.write_all(b" ")?;
        }
        let value = if counts.get(idx).copied().unwrap_or(0) > 0 {
            sums[idx] as f64 / counts[idx] as f64 - k as f64
        } else {
            f64::from(*length) - k as f64
        }
        .abs();
        write!(writer, "{value}")?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

pub(crate) fn write_paired_flens(path: &Path, flens: &[u32]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_paired_flens_line(&mut writer, flens)?;
    Ok(())
}

pub(crate) fn write_batch_paired_flens(path: &Path, batch_flens: &[Vec<u32>]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for flens in batch_flens {
        write_paired_flens_line(&mut writer, flens)?;
    }
    Ok(())
}

fn write_paired_flens_line<W: Write>(writer: &mut W, flens: &[u32]) -> Result<()> {
    for (idx, value) in flens.iter().enumerate() {
        if idx > 0 {
            writer.write_all(b" ")?;
        }
        write!(writer, "{value}")?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

pub(crate) fn write_run_info(path: &Path, info: RunInfo) -> Result<()> {
    let p_pseudoaligned = if info.processed == 0 {
        0.0
    } else {
        info.aligned as f64 * 100.0 / info.processed as f64
    };
    let p_unique = if info.processed == 0 {
        0.0
    } else {
        info.unique as f64 * 100.0 / info.processed as f64
    };
    let start_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_secs().to_string())
        .unwrap_or_else(|_| "0".to_string());
    let call = std::env::args().collect::<Vec<_>>().join(" ");
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "{{")?;
    writeln!(writer, "  \"n_targets\": {},", info.targets)?;
    writeln!(writer, "  \"n_bootstraps\": 0,")?;
    writeln!(writer, "  \"n_processed\": {},", info.processed)?;
    writeln!(writer, "  \"n_pseudoaligned\": {},", info.aligned)?;
    writeln!(writer, "  \"n_unique\": {},", info.unique)?;
    writeln!(writer, "  \"p_pseudoaligned\": {p_pseudoaligned:.1},")?;
    writeln!(writer, "  \"p_unique\": {p_unique:.1},")?;
    writeln!(
        writer,
        "  \"kallisto_version\": \"kallistors {}\",",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(writer, "  \"index_version\": 13,")?;
    writeln!(writer, "  \"k-mer length\": {},", info.k)?;
    writeln!(writer, "  \"start_time\": \"{start_time}\",")?;
    if let Some(frame_clashes) = info.frame_clashes {
        writeln!(writer, "  \"call\": \"{}\",", json_escape(&call))?;
        writeln!(writer, "  \"n_frame_clashes\": {frame_clashes}")?;
    } else {
        writeln!(writer, "  \"call\": \"{}\"", json_escape(&call))?;
    }
    writeln!(writer, "}}")?;
    Ok(())
}

fn json_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            other => escaped.push(other),
        }
    }
    escaped
}

pub(crate) fn fake_barcode(value: u64, len: usize) -> Vec<u8> {
    binary_to_string(value, len.min(32)).into_bytes()
}

fn binary_to_string(value: u64, len: usize) -> String {
    let mut out = String::with_capacity(len);
    for idx in 0..len {
        let shift = 2 * (len - idx - 1);
        let base = match (value >> shift) & 0x03 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };
        out.push(base);
    }
    out
}
