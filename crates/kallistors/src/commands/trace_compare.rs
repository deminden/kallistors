use anyhow::{Result, anyhow};
use kallistors::io::ReadSource;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn normalize_header(input: &[u8]) -> String {
    let text = String::from_utf8_lossy(input);
    text.trim()
        .strip_prefix('@')
        .unwrap_or(text.trim())
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn format_ec(ec: Option<&Vec<u32>>) -> String {
    match ec {
        Some(values) if !values.is_empty() => values
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(","),
        _ => "-".to_string(),
    }
}

fn format_reason(reason: Option<kallistors::pseudoalign::DebugFailReason>) -> String {
    reason
        .map(|value| value.to_string())
        .unwrap_or_else(|| "ok".to_string())
}

fn diff_vec_usize(left: Option<&Vec<usize>>, right: Option<&Vec<usize>>) -> String {
    match (left, right) {
        (Some(a), Some(b)) => {
            let max_len = a.len().max(b.len());
            for idx in 0..max_len {
                let va = a.get(idx);
                let vb = b.get(idx);
                if va != vb {
                    return format!(
                        "idx={idx};baseline={};fast={}",
                        va.map(|v| v.to_string()).unwrap_or_else(|| "-".to_string()),
                        vb.map(|v| v.to_string()).unwrap_or_else(|| "-".to_string())
                    );
                }
            }
            "-".to_string()
        }
        (None, None) => "-".to_string(),
        (Some(_), None) => "baseline_only".to_string(),
        (None, Some(_)) => "fast_only".to_string(),
    }
}

fn diff_candidates(
    left: Option<&Vec<kallistors::pseudoalign::MinimizerCandidateTrace>>,
    right: Option<&Vec<kallistors::pseudoalign::MinimizerCandidateTrace>>,
) -> String {
    match (left, right) {
        (Some(a), Some(b)) => {
            let max_len = a.len().max(b.len());
            for idx in 0..max_len {
                let va = a.get(idx);
                let vb = b.get(idx);
                if va.map(candidate_sig) != vb.map(candidate_sig) {
                    return format!(
                        "idx={idx};baseline={};fast={}",
                        va.map(candidate_sig).unwrap_or_else(|| "-".to_string()),
                        vb.map(candidate_sig).unwrap_or_else(|| "-".to_string())
                    );
                }
            }
            "-".to_string()
        }
        (None, None) => "-".to_string(),
        (Some(_), None) => "baseline_only".to_string(),
        (None, Some(_)) => "fast_only".to_string(),
    }
}

fn candidate_sig(value: &kallistors::pseudoalign::MinimizerCandidateTrace) -> String {
    format!(
        "{}:{}:{}:{}:{}:{}:{}",
        value.read_pos,
        value.min_pos,
        value.minimizer,
        value.mphf_hit as u8,
        value.positions_len,
        value.overcrowded as u8,
        value.matched as u8
    )
}

fn diff_jumps(
    left: Option<&Vec<kallistors::pseudoalign::JumpDecisionTrace>>,
    right: Option<&Vec<kallistors::pseudoalign::JumpDecisionTrace>>,
) -> String {
    match (left, right) {
        (Some(a), Some(b)) => {
            let max_len = a.len().max(b.len());
            for idx in 0..max_len {
                let va = a.get(idx);
                let vb = b.get(idx);
                if va.map(jump_sig) != vb.map(jump_sig) {
                    return format!(
                        "idx={idx};baseline={};fast={}",
                        va.map(jump_sig).unwrap_or_else(|| "-".to_string()),
                        vb.map(jump_sig).unwrap_or_else(|| "-".to_string())
                    );
                }
            }
            "-".to_string()
        }
        (None, None) => "-".to_string(),
        (Some(_), None) => "baseline_only".to_string(),
        (None, Some(_)) => "fast_only".to_string(),
    }
}

fn jump_sig(value: &kallistors::pseudoalign::JumpDecisionTrace) -> String {
    format!(
        "{}:{}:{}:{}:{}:{}:{}:{}",
        value.read_pos,
        value.matched_unitig_id,
        value.matched_block_idx,
        value.jump_distance,
        value.next_pos,
        value.next_hit_found as u8,
        value.jumped as u8,
        value.reason
    )
}

fn format_investigation(options: kallistors::pseudoalign::InvestigationOptions) -> String {
    let mut parts = Vec::new();
    if options.trace_fast_path {
        parts.push("trace_fast_path");
    }
    if options.bounded_incremental_scan {
        parts.push("bounded_incremental_scan");
    }
    if options.faster_probe_state {
        parts.push("faster_probe_state");
    }
    if options.altered_hit_bookkeeping {
        parts.push("altered_hit_bookkeeping");
    }
    if options.candidate_reuse {
        parts.push("candidate_reuse");
    }
    if parts.is_empty() {
        "-".to_string()
    } else {
        parts.join(",")
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    index: &Path,
    reads: &Path,
    reads2: Option<&Path>,
    read_list: &Path,
    out: &Path,
    strand: kallistors::pseudoalign::Strand,
    fragment_length: Option<f64>,
    single_overhang: bool,
    min_range: usize,
    do_union: bool,
    no_jump: bool,
    dfk_onlist: bool,
    fr_stranded: bool,
    rf_stranded: bool,
    kallisto_enum: bool,
    kallisto_strict: bool,
    kallisto_local_fallback: bool,
    kallisto_fallback: bool,
    discard_special_only: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_direct_kmer: bool,
    kallisto_bifrost_find: bool,
    kallisto_sparse_hits: bool,
) -> Result<()> {
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
        if kallisto_direct_kmer {
            kallistors::pseudoalign::build_bifrost_index_with_kmer_pos(index, true)
        } else if kallisto_fallback {
            kallistors::pseudoalign::build_bifrost_index_with_positions_and_kmer(index, true)
        } else {
            kallistors::pseudoalign::build_bifrost_index_with_positions(index, true)
        }
    } else if kallisto_direct_kmer {
        kallistors::pseudoalign::build_bifrost_index_with_kmer_pos(index, false)
    } else if kallisto_fallback {
        kallistors::pseudoalign::build_bifrost_index_with_kmer(index)
    } else {
        kallistors::pseudoalign::build_bifrost_index(index)
    }?;

    let baseline = kallistors::pseudoalign::PseudoalignOptions {
        min_range,
        do_union,
        dfk_onlist,
        strand_specific,
        no_jump,
        kallisto_enum,
        kallisto_strict,
        kallisto_local_fallback,
        kallisto_fallback,
        discard_special_only,
        skip_overcrowded_minimizer,
        kallisto_direct_kmer,
        kallisto_bifrost_find,
        kallisto_sparse_hits,
        bias: false,
        max_bias: 0,
        investigation: kallistors::pseudoalign::InvestigationOptions::default(),
    };
    let mut fast_investigation = super::investigation_options_from_env();
    fast_investigation.trace_fast_path = true;
    let fast = kallistors::pseudoalign::PseudoalignOptions {
        investigation: fast_investigation,
        ..baseline
    };

    let mut targets = HashSet::new();
    for line in BufReader::new(File::open(read_list)?).lines() {
        let key = line?;
        if !key.trim().is_empty() {
            targets.insert(key.trim().to_string());
        }
    }
    if targets.is_empty() {
        return Err(anyhow!("read list is empty"));
    }

    let mut reader1 = kallistors::io::open_fastq_reader(reads)?;
    let mut reader2 = if let Some(path) = reads2 {
        Some(kallistors::io::open_fastq_reader(path)?)
    } else {
        None
    };
    let mut writer = BufWriter::new(File::create(out)?);
    writeln!(
        writer,
        "read\tmode\taligned_baseline\taligned_fast\treason_baseline\treason_fast\tec_baseline\tec_fast\tprobe_diff\tcandidate_diff\tjump_diff\tleft_reason_baseline\tleft_reason_fast\tright_reason_baseline\tright_reason_fast\tleft_ec_baseline\tleft_ec_fast\tright_ec_baseline\tright_ec_fast\tinvestigation"
    )?;

    while let Some(record1) = reader1.next_record() {
        let record1 = record1?;
        let key = normalize_header(&record1.header);
        if !targets.contains(&key) {
            if let Some(reader2) = reader2.as_mut() {
                let _ = reader2.next_record().transpose()?;
            }
            continue;
        }
        if let Some(reader2) = reader2.as_mut() {
            let record2 = reader2
                .next_record()
                .transpose()?
                .ok_or_else(|| anyhow!("paired FASTQ length mismatch"))?;
            let baseline_trace = kallistors::pseudoalign::trace_read_pair_bifrost(
                &index,
                &record1.seq,
                &record2.seq,
                strand,
                baseline,
            );
            let fast_trace = kallistors::pseudoalign::trace_read_pair_bifrost(
                &index,
                &record1.seq,
                &record2.seq,
                strand,
                fast,
            );
            let probe_diff = format!(
                "left={};right={}",
                diff_vec_usize(
                    baseline_trace.left.positions_visited.as_ref(),
                    fast_trace.left.positions_visited.as_ref()
                ),
                diff_vec_usize(
                    baseline_trace.right.positions_visited.as_ref(),
                    fast_trace.right.positions_visited.as_ref()
                )
            );
            let candidate_diff = format!(
                "left={};right={}",
                diff_candidates(
                    baseline_trace.left.minimizer_candidates.as_ref(),
                    fast_trace.left.minimizer_candidates.as_ref()
                ),
                diff_candidates(
                    baseline_trace.right.minimizer_candidates.as_ref(),
                    fast_trace.right.minimizer_candidates.as_ref()
                )
            );
            let jump_diff = format!(
                "left={};right={}",
                diff_jumps(
                    baseline_trace.left.jump_decisions.as_ref(),
                    fast_trace.left.jump_decisions.as_ref()
                ),
                diff_jumps(
                    baseline_trace.right.jump_decisions.as_ref(),
                    fast_trace.right.jump_decisions.as_ref()
                )
            );
            writeln!(
                writer,
                "{}\tpaired\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                key,
                if baseline_trace.merged_reason.is_none() && baseline_trace.merged_ec.is_some() {
                    1
                } else {
                    0
                },
                if fast_trace.merged_reason.is_none() && fast_trace.merged_ec.is_some() {
                    1
                } else {
                    0
                },
                baseline_trace
                    .merged_reason
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "ok".to_string()),
                fast_trace
                    .merged_reason
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "ok".to_string()),
                format_ec(baseline_trace.merged_ec.as_ref()),
                format_ec(fast_trace.merged_ec.as_ref()),
                probe_diff,
                candidate_diff,
                jump_diff,
                format_reason(baseline_trace.left.reason),
                format_reason(fast_trace.left.reason),
                format_reason(baseline_trace.right.reason),
                format_reason(fast_trace.right.reason),
                format_ec(baseline_trace.left.ec_after_strand_filter.as_ref()),
                format_ec(fast_trace.left.ec_after_strand_filter.as_ref()),
                format_ec(baseline_trace.right.ec_after_strand_filter.as_ref()),
                format_ec(fast_trace.right.ec_after_strand_filter.as_ref()),
                format_investigation(fast_investigation)
            )?;
        } else {
            let (baseline_trace, _) = kallistors::pseudoalign::trace_read_bifrost_with_hits(
                &index,
                &record1.seq,
                strand,
                filter,
                baseline,
            );
            let (fast_trace, _) = kallistors::pseudoalign::trace_read_bifrost_with_hits(
                &index,
                &record1.seq,
                strand,
                filter,
                fast,
            );
            writeln!(
                writer,
                "{}\tsingle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t-\t-\t-\t-\t{}\t{}\t-\t-\t{}",
                key,
                if baseline_trace.reason.is_none() {
                    1
                } else {
                    0
                },
                if fast_trace.reason.is_none() { 1 } else { 0 },
                format_reason(baseline_trace.reason),
                format_reason(fast_trace.reason),
                format_ec(baseline_trace.ec_after_strand_filter.as_ref()),
                format_ec(fast_trace.ec_after_strand_filter.as_ref()),
                diff_vec_usize(
                    baseline_trace.positions_visited.as_ref(),
                    fast_trace.positions_visited.as_ref()
                ),
                diff_candidates(
                    baseline_trace.minimizer_candidates.as_ref(),
                    fast_trace.minimizer_candidates.as_ref()
                ),
                diff_jumps(
                    baseline_trace.jump_decisions.as_ref(),
                    fast_trace.jump_decisions.as_ref()
                ),
                format_ec(baseline_trace.ec_after_strand_filter.as_ref()),
                format_ec(fast_trace.ec_after_strand_filter.as_ref()),
                format_investigation(fast_investigation)
            )?;
        }
    }

    Ok(())
}
