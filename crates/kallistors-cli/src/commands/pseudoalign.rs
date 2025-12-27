use anyhow::{Result, anyhow};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

#[allow(clippy::too_many_arguments)]
pub fn run(
    index: &Path,
    reads: &Path,
    reads2: Option<&Path>,
    out: &Path,
    _threads: usize,
    debug_out: Option<&Path>,
    debug_max: usize,
    strand: kallistors::pseudoalign::Strand,
    fragment_length: Option<f64>,
    single_overhang: bool,
    min_range: usize,
    do_union: bool,
    no_jump: bool,
    dfk_onlist: bool,
    fr_stranded: bool,
    rf_stranded: bool,
    bias: bool,
    bias_out: Option<&Path>,
) -> Result<()> {
    if reads.extension().and_then(|s| s.to_str()) == Some("gz") {
        return Err(anyhow!("gzip FASTQ not supported yet"));
    }
    if reads2.is_some() && fragment_length.is_some() {
        return Err(anyhow!(
            "--fragment-length is only supported for single-end reads"
        ));
    }

    let file = File::open(reads)?;
    let mut reader = kallistors::io::FastqReader::new(BufReader::new(file));
    let mut reader2 = if let Some(path) = reads2 {
        let file = File::open(path)?;
        Some(kallistors::io::FastqReader::new(BufReader::new(file)))
    } else {
        None
    };
    if fr_stranded && rf_stranded {
        return Err(anyhow!(
            "--fr-stranded and --rf-stranded are mutually exclusive"
        ));
    }
    let strand_specific = if fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else {
        None
    };
    let options = kallistors::pseudoalign::PseudoalignOptions {
        min_range,
        do_union,
        dfk_onlist,
        strand_specific,
        no_jump,
        bias,
        max_bias: 1_000_000,
    };
    let filter =
        fragment_length
            .filter(|v| *v > 0.0)
            .map(|v| kallistors::pseudoalign::FragmentFilter {
                fragment_length: v as u32,
                single_overhang,
            });
    let load_positional_info = filter
        .map(|v| !v.single_overhang && v.fragment_length > 0)
        .unwrap_or(false);
    let index = if load_positional_info {
        kallistors::pseudoalign::build_bifrost_index_with_positions(index, true)
    } else {
        kallistors::pseudoalign::build_bifrost_index(index)
    }
    .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
    let (res, report) = if let Some(reader2) = reader2.as_mut() {
        if debug_out.is_some() {
            let (res, report) =
                kallistors::pseudoalign::pseudoalign_paired_bifrost_debug_with_options(
                    &index,
                    &mut reader,
                    reader2,
                    strand,
                    debug_max,
                    options,
                )
                .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
            (res, Some(report))
        } else {
            (
                kallistors::pseudoalign::pseudoalign_paired_bifrost_with_options(
                    &index,
                    &mut reader,
                    reader2,
                    strand,
                    options,
                )
                .map_err(|err| anyhow!("pseudoalign failed: {err}"))?,
                None,
            )
        }
    } else if debug_out.is_some() {
        let (res, report) =
            kallistors::pseudoalign::pseudoalign_single_end_bifrost_debug_with_options(
                &index,
                &mut reader,
                strand,
                debug_max,
                filter,
                options,
            )
            .map_err(|err| anyhow!("pseudoalign failed: {err}"))?;
        (res, Some(report))
    } else {
        (
            kallistors::pseudoalign::pseudoalign_single_end_bifrost_with_options(
                &index,
                &mut reader,
                strand,
                filter,
                options,
            )
            .map_err(|err| anyhow!("pseudoalign failed: {err}"))?,
            None,
        )
    };

    let out_file = File::create(out)?;
    let mut writer = BufWriter::new(out_file);
    writeln!(writer, "#ec_id\tcount\ttranscripts")?;
    for (i, ec) in res.ec_list.iter().enumerate() {
        let count = res.counts[i];
        write!(writer, "{}\t{}\t", i, count)?;
        for (j, t) in ec.iter().enumerate() {
            if j > 0 {
                write!(writer, ",")?;
            }
            write!(writer, "{}", t)?;
        }
        writeln!(writer)?;
    }

    println!(
        "processed: {}, aligned: {}, ecs: {}",
        res.reads_processed,
        res.reads_aligned,
        res.ec_list.len()
    );

    if let Some(path) = debug_out {
        if let Some(report) = report {
            write_debug_report(path, &res, &report)?;
        }
    }
    if bias {
        let path = bias_out
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| out.with_extension("bias.txt"));
        if let Some(bias_counts) = res.bias.as_ref() {
            bias_counts.write_text(&path)?;
        }
    }
    Ok(())
}

fn write_debug_report(
    path: &Path,
    res: &kallistors::pseudoalign::EcCounts,
    report: &kallistors::pseudoalign::DebugReport,
) -> Result<()> {
    let out_file = File::create(path)?;
    let mut writer = BufWriter::new(out_file);
    writeln!(writer, "#reads_processed\t{}", res.reads_processed)?;
    writeln!(writer, "#reads_aligned\t{}", res.reads_aligned)?;
    writeln!(writer, "#summary")?;
    writeln!(writer, "reason\tcount")?;
    let order = [
        kallistors::pseudoalign::DebugFailReason::NoValidKmer,
        kallistors::pseudoalign::DebugFailReason::NoMinimizerHit,
        kallistors::pseudoalign::DebugFailReason::NoPositions,
        kallistors::pseudoalign::DebugFailReason::NoKmerMatch,
        kallistors::pseudoalign::DebugFailReason::EmptyEc,
        kallistors::pseudoalign::DebugFailReason::IntersectionEmpty,
        kallistors::pseudoalign::DebugFailReason::Unknown,
    ];
    for reason in order {
        let count = report.counts.get(&reason).copied().unwrap_or(0);
        writeln!(writer, "{}\t{}", reason, count)?;
    }
    writeln!(writer, "#traces")?;
    writeln!(
        writer,
        "read\treason\tkmer_pos\tmin_pos\tkmer_seq\tmin_seq\tpositions\tsample_positions\tused_revcomp"
    )?;
    for trace in &report.traces {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            trace.header,
            trace.reason,
            trace
                .kmer_pos
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            trace
                .min_pos
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            trace.kmer_seq.as_deref().unwrap_or("-"),
            trace.min_seq.as_deref().unwrap_or("-"),
            trace
                .positions
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            trace.sample_positions.as_deref().unwrap_or("-"),
            trace.used_revcomp,
        )?;
    }
    Ok(())
}
