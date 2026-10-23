//! Paired read traces and optional final-EC trace output.

use super::super::paired::estimate_fragment_length_for_pair;
use super::super::*;
use crate::io::FastqRecord;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn trace_read_pair_bifrost(
    index: &BifrostIndex,
    left_seq: &[u8],
    right_seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> PairedTraceResult {
    let mut dbg1 = ReadDebugState::default();
    let mut dbg2 = ReadDebugState::default();
    let mut ec1 =
        super::super::ec_for_read_bifrost(index, left_seq, strand, Some(&mut dbg1), options);
    let mut ec2 =
        super::super::ec_for_read_bifrost(index, right_seq, strand, Some(&mut dbg2), options);

    if let Some(mode) = options.strand_specific {
        let comprehensive = options.do_union || options.no_jump || index.use_shade;
        if let Some(read_ec) = ec1.as_mut()
            && let Some(filtered) = super::super::apply_strand_filter(
                index,
                read_ec,
                mode,
                true,
                comprehensive,
                &mut None,
                b"",
            )
        {
            read_ec.hard_reject_pair |= filtered.is_empty();
            read_ec.ec = filtered;
        }
        if let Some(read_ec) = ec2.as_mut()
            && let Some(filtered) = super::super::apply_strand_filter(
                index,
                read_ec,
                mode,
                false,
                comprehensive,
                &mut None,
                b"",
            )
        {
            read_ec.hard_reject_pair |= filtered.is_empty();
            read_ec.ec = filtered;
        }
    }

    let mut left_reason_override = None;
    let mut right_reason_override = None;
    if ec1.as_ref().is_some_and(|v| v.ec.is_empty())
        && ec1.as_ref().is_some_and(|v| !v.hard_reject_pair)
    {
        left_reason_override = Some(DebugFailReason::EmptyEc);
        ec1 = None;
    }
    if ec2.as_ref().is_some_and(|v| v.ec.is_empty())
        && ec2.as_ref().is_some_and(|v| !v.hard_reject_pair)
    {
        right_reason_override = Some(DebugFailReason::EmptyEc);
        ec2 = None;
    }

    let hard_reject_pair = ec1.as_ref().is_some_and(|v| v.hard_reject_pair)
        || ec2.as_ref().is_some_and(|v| v.hard_reject_pair);

    let mut merged = Vec::new();
    let mut had_offlist = false;
    let merged_reason = if hard_reject_pair {
        Some(DebugFailReason::IntersectionEmpty)
    } else {
        match (ec1.as_ref(), ec2.as_ref()) {
            (None, None) => {
                if left_reason_override.is_none() && right_reason_override.is_none() {
                    left_reason_override = Some(DebugFailReason::NoKmerMatch);
                    right_reason_override = Some(DebugFailReason::NoKmerMatch);
                }
                None
            }
            (Some(ec), None) | (None, Some(ec)) => {
                merged.extend_from_slice(&ec.ec);
                had_offlist = ec.had_offlist;
                None
            }
            (Some(left), Some(right)) => {
                super::super::intersect_sorted(&left.ec, &right.ec, &mut merged);
                had_offlist = left.had_offlist || right.had_offlist;
                if merged.is_empty() {
                    Some(DebugFailReason::IntersectionEmpty)
                } else {
                    None
                }
            }
        }
    };

    let left_before = ec1.as_ref().map(|v| v.ec.clone());
    let right_before = ec2.as_ref().map(|v| v.ec.clone());
    let estimated_fragment_length = if merged.len() == 1 {
        match (ec1.as_ref(), ec2.as_ref()) {
            (Some(left), Some(right)) => {
                estimate_fragment_length_for_pair(index, merged[0], left, right)
            }
            _ => None,
        }
    } else {
        None
    };
    let placement = if merged.len() == 1 {
        ec1.as_ref()
            .and_then(|left| {
                left.best_match.and_then(|best_match| {
                    super::super::placement_for_match(&merged, best_match, left_seq.len())
                })
            })
            .or_else(|| {
                ec2.as_ref().and_then(|right| {
                    right.best_match.and_then(|best_match| {
                        super::super::placement_for_match(&merged, best_match, right_seq.len())
                    })
                })
            })
    } else {
        None
    };
    let left_trace = super::super::trace_result_from_debug_state(
        index,
        left_seq,
        &dbg1,
        ec1.is_some(),
        left_before.clone(),
        left_before.clone(),
        left_before,
        left_reason_override,
    );
    let right_trace = super::super::trace_result_from_debug_state(
        index,
        right_seq,
        &dbg2,
        ec2.is_some(),
        right_before.clone(),
        right_before.clone(),
        right_before,
        right_reason_override,
    );

    PairedTraceResult {
        left: left_trace,
        right: right_trace,
        merged_ec: (!merged.is_empty()).then_some(merged),
        merged_reason,
        hard_reject_pair,
        had_offlist,
        estimated_fragment_length,
        placement,
    }
}

pub(in crate::pseudoalign) struct PairedFinalTraceWriter {
    writer: BufWriter<File>,
}

pub(in crate::pseudoalign) struct PairedFinalTraceRow<'a> {
    pub(in crate::pseudoalign) read_id: u64,
    pub(in crate::pseudoalign) left: &'a FastqRecord,
    pub(in crate::pseudoalign) right: &'a FastqRecord,
    pub(in crate::pseudoalign) merged_ec: Option<&'a [u32]>,
    pub(in crate::pseudoalign) reason: &'a str,
    pub(in crate::pseudoalign) left_ec: Option<&'a [u32]>,
    pub(in crate::pseudoalign) right_ec: Option<&'a [u32]>,
}

impl PairedFinalTraceWriter {
    pub(in crate::pseudoalign) fn from_env() -> Result<Option<Self>> {
        let Some(path) = std::env::var_os("KALLISTORS_PAIRED_FINAL_TRACE") else {
            return Ok(None);
        };
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "read_id\tread_name1\tread_name2\taligned\treason\tmerged_ec\tleft_ec\tright_ec"
        )?;
        Ok(Some(Self { writer }))
    }

    pub(in crate::pseudoalign) fn write(&mut self, row: PairedFinalTraceRow<'_>) -> Result<()> {
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.read_id,
            normalize_trace_header(&row.left.header),
            normalize_trace_header(&row.right.header),
            if row.merged_ec.is_some_and(|ec| !ec.is_empty()) {
                1
            } else {
                0
            },
            row.reason,
            format_trace_ec(row.merged_ec),
            format_trace_ec(row.left_ec),
            format_trace_ec(row.right_ec),
        )?;
        Ok(())
    }
}

fn normalize_trace_header(header: &[u8]) -> String {
    let text = String::from_utf8_lossy(header);
    text.trim()
        .strip_prefix('@')
        .unwrap_or(text.trim())
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn format_trace_ec(ec: Option<&[u32]>) -> String {
    let Some(ec) = ec else {
        return "-".to_string();
    };
    if ec.is_empty() {
        return "-".to_string();
    }
    ec.iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",")
}
