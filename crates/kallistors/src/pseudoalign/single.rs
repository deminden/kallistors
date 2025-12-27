use std::collections::HashMap;

use crate::Result;
use crate::bias::BiasCounts;
use crate::io::ReadSource;

use super::{
    BifrostIndex, DebugFailReason, DebugReport, EcCounts, FragmentFilter, KmerEcIndex,
    PseudoalignOptions, ReadDebugState, Strand,
};

/// Pseudoalign single-end reads using a naive k-mer map.
pub fn pseudoalign_single_end<R: ReadSource>(
    index: &KmerEcIndex,
    reader: &mut R,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;

    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();

    while let Some(record) = reader.next_record() {
        let record = record?;
        reads_processed += 1;
        if record.seq.len() < index.k {
            continue;
        }

        let mut has_hit = false;
        current.clear();

        for pos in 0..=record.seq.len() - index.k {
            if let Some((fwd, rev)) = super::encode_kmer_pair(&record.seq[pos..pos + index.k])
                && let Some(ec) =
                    super::lookup_ec(index, fwd).or_else(|| super::lookup_ec(index, rev))
            {
                if !has_hit {
                    current.extend_from_slice(ec);
                    has_hit = true;
                } else {
                    next.clear();
                    super::intersect_sorted(&current, ec, &mut next);
                    std::mem::swap(&mut current, &mut next);
                    if current.is_empty() {
                        break;
                    }
                }
            }
        }

        if !has_hit || current.is_empty() {
            continue;
        }

        reads_aligned += 1;
        let ec_id = match ec_map.get(&current) {
            Some(id) => *id,
            None => {
                let id = ec_list.len();
                ec_map.insert(current.clone(), id);
                ec_list.push(current.clone());
                counts.push(0);
                id
            }
        };
        counts[ec_id] += 1;
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias: None,
        fragment_length_stats: None,
        fragment_length_hist: None,
    })
}

pub fn pseudoalign_single_end_bifrost<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        Strand::Unstranded,
        None,
        None,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_single_end_bifrost_debug<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    max_traces: usize,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        Strand::Unstranded,
        Some(&mut report),
        None,
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_single_end_bifrost_with_strand<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        None,
        None,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_single_end_bifrost_debug_with_strand<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    max_traces: usize,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        Some(&mut report),
        None,
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_single_end_bifrost_with_strand_and_filter<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    filter: Option<FragmentFilter>,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        None,
        filter,
        PseudoalignOptions::default(),
    )
}

pub fn pseudoalign_single_end_bifrost_debug_with_strand_and_filter<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    max_traces: usize,
    filter: Option<FragmentFilter>,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        Some(&mut report),
        filter,
        PseudoalignOptions::default(),
    )?;
    Ok((counts, report))
}

pub fn pseudoalign_single_end_bifrost_with_options<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    pseudoalign_single_end_bifrost_inner(index, reader, strand, None, filter, options)
}

pub fn pseudoalign_single_end_bifrost_debug_with_options<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    max_traces: usize,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> Result<(EcCounts, DebugReport)> {
    let mut report = DebugReport::new(max_traces);
    let counts = pseudoalign_single_end_bifrost_inner(
        index,
        reader,
        strand,
        Some(&mut report),
        filter,
        options,
    )?;
    Ok((counts, report))
}

fn pseudoalign_single_end_bifrost_inner<R: ReadSource>(
    index: &BifrostIndex,
    reader: &mut R,
    strand: Strand,
    mut report: Option<&mut DebugReport>,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
) -> Result<EcCounts> {
    let mut ec_list: Vec<Vec<u32>> = Vec::new();
    let mut ec_map: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut counts: Vec<u32> = Vec::new();
    let mut reads_processed = 0u64;
    let mut reads_aligned = 0u64;
    let mut bias = if options.bias {
        Some(BiasCounts::new())
    } else {
        None
    };

    while let Some(record) = reader.next_record() {
        let record = record?;
        reads_processed += 1;

        let mut dbg = report.as_ref().map(|_| ReadDebugState::default());
        let ec = super::ec_for_read_bifrost(index, &record.seq, strand, dbg.as_mut(), options);
        let Some(mut read_ec) = ec else {
            if let Some(r) = report.as_deref_mut() {
                if let Some(state) = dbg {
                    let (reason, kmer_pos, min_pos) = super::debug_reason(&state);
                    let kmer_seq = kmer_pos.and_then(|pos| {
                        if pos + index.k <= record.seq.len() {
                            Some(
                                String::from_utf8_lossy(&record.seq[pos..pos + index.k])
                                    .into_owned(),
                            )
                        } else {
                            None
                        }
                    });
                    let min_seq = match (kmer_pos, min_pos) {
                        (Some(kp), Some(mp)) => {
                            let start = kp + mp;
                            let end = start + index.g;
                            if end <= record.seq.len() {
                                Some(String::from_utf8_lossy(&record.seq[start..end]).into_owned())
                            } else {
                                None
                            }
                        }
                        _ => None,
                    };
                    let (positions, sample_positions) = state
                        .first_no_match_positions
                        .as_ref()
                        .map(|v| (Some(v.2), Some(v.3.clone())))
                        .unwrap_or((None, None));
                    r.record(
                        &record.header,
                        reason,
                        kmer_pos,
                        min_pos,
                        kmer_seq,
                        min_seq,
                        positions,
                        sample_positions,
                        state.used_revcomp,
                    );
                } else {
                    r.record(
                        &record.header,
                        DebugFailReason::Unknown,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        false,
                    );
                }
            }
            continue;
        };
        if let Some(filter) = filter
            && !filter.single_overhang
            && filter.fragment_length > 0
            && let Some(best_match) = read_ec.best_match
        {
            read_ec.ec = super::filter_ec_by_fragment(
                index,
                &read_ec.ec,
                best_match,
                filter.fragment_length as i64,
            );
        }
        if let Some(mode) = options.strand_specific {
            let comprehensive = options.do_union || options.no_jump || index.use_shade;
            if let Some(filtered) = super::apply_strand_filter(
                index,
                &read_ec,
                mode,
                true,
                comprehensive,
                &mut report,
                &record.header,
            ) {
                read_ec.ec = filtered;
            }
        }
        if read_ec.ec.is_empty() {
            if let Some(r) = report.as_deref_mut() {
                if let Some(state) = dbg {
                    let (reason, kmer_pos, min_pos) = super::debug_reason(&state);
                    r.record(
                        &record.header,
                        reason,
                        kmer_pos,
                        min_pos,
                        None,
                        None,
                        None,
                        None,
                        state.used_revcomp,
                    );
                } else {
                    r.record(
                        &record.header,
                        DebugFailReason::Unknown,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        false,
                    );
                }
            }
            continue;
        }

        if let Some(bias_counts) = bias.as_mut()
            && bias_counts.total < options.max_bias as u64
            && let Some(best_match) = read_ec.best_match
            && let Some(hex) = super::bias_hexamer_for_match(index, best_match)
        {
            bias_counts.record(hex);
        }

        reads_aligned += 1;
        let ec_id = match ec_map.get(&read_ec.ec) {
            Some(id) => *id,
            None => {
                let id = ec_list.len();
                ec_map.insert(read_ec.ec.clone(), id);
                ec_list.push(read_ec.ec.clone());
                counts.push(0);
                id
            }
        };
        counts[ec_id] += 1;
    }

    Ok(EcCounts {
        ec_list,
        counts,
        reads_processed,
        reads_aligned,
        bias,
        fragment_length_stats: None,
        fragment_length_hist: None,
    })
}
