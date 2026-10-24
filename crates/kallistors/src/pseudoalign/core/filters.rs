use super::super::ec::merge_sorted_unique_vec;
use super::super::{
    BifrostIndex, DebugFailReason, DebugReport, Hit, MatchInfo, ReadEc, StrandSpecific,
};
use std::borrow::Cow;

// Fragment and strand filters interpret compact EC blocks back into transcript coordinates.
pub(crate) fn filter_ec_by_fragment(
    index: &BifrostIndex,
    ec: &[u32],
    best_match: MatchInfo,
    fragment_length: i64,
) -> Vec<u32> {
    let mut filtered = ec.to_vec();
    filter_ec_by_fragment_in_place(index, &mut filtered, best_match, fragment_length);
    filtered
}

// The scanner owns this EC already. Retain its storage and exact transcript order.
pub(crate) fn filter_ec_by_fragment_in_place(
    index: &BifrostIndex,
    ec: &mut Vec<u32>,
    best_match: MatchInfo,
    fragment_length: i64,
) {
    if let Some(keep) = fragment_predicate(index, best_match, fragment_length) {
        ec.retain(|&tr| keep(tr));
    }
}

// Keep an unchanged index EC borrowed. Allocate only after the first rejected
// transcript; owned intersections continue to filter in their existing buffer.
pub(crate) fn filter_ec_by_fragment_cow<'a>(
    index: &BifrostIndex,
    ec: Cow<'a, [u32]>,
    best_match: MatchInfo,
    fragment_length: i64,
) -> Cow<'a, [u32]> {
    let Some(keep) = fragment_predicate(index, best_match, fragment_length) else {
        return ec;
    };
    match ec {
        Cow::Owned(mut values) => {
            values.retain(|&tr| keep(tr));
            Cow::Owned(values)
        }
        Cow::Borrowed(values) => {
            let Some(first_reject) = values.iter().position(|&tr| !keep(tr)) else {
                return Cow::Borrowed(values);
            };
            let mut filtered = Vec::with_capacity(values.len().saturating_sub(1));
            filtered.extend_from_slice(&values[..first_reject]);
            filtered.extend(
                values[first_reject + 1..]
                    .iter()
                    .copied()
                    .filter(|&tr| keep(tr)),
            );
            Cow::Owned(filtered)
        }
    }
}

fn fragment_predicate(
    index: &BifrostIndex,
    best_match: MatchInfo,
    fragment_length: i64,
) -> Option<impl Fn(u32) -> bool + '_> {
    if fragment_length <= 0 {
        return None;
    }
    if index.transcript_lengths.is_empty() {
        return None;
    }
    if best_match.unitig_id >= index.flat_ec.owned_ec_blocks.len() {
        return None;
    }
    let blocks = &index.flat_ec.owned_ec_blocks[best_match.unitig_id];
    if blocks.is_empty() || blocks[0].positions.is_none() {
        return None;
    }
    let unitig_len = if best_match.unitig_id < index.unitig_count() {
        index
            .encoded_unitigs
            .unitig_bases
            .len(best_match.unitig_id)
            .unwrap_or(index.k)
    } else {
        index.k
    };
    // All candidates use the same matched unitig coordinate. Select its
    // annotation block once; individual transcripts still use exact lookup.
    let block = ec_blocks_leading_vals(blocks, best_match.unitig_pos as u32).last();
    Some(move |tr| {
        let Some(block) = block else {
            return false;
        };
        let Some(raw_positions) = block_positions(block, tr) else {
            return false;
        };
        let len = index
            .transcript_lengths
            .get(tr as usize)
            .copied()
            .unwrap_or(0) as i64;
        raw_positions.any(|raw| {
            map_transcript_position(
                blocks,
                tr,
                raw,
                best_match.unitig_pos,
                best_match.read_pos,
                best_match.used_revcomp,
                unitig_len,
                index.k,
            )
            .is_some_and(|(pos, forward)| {
                if forward {
                    pos + fragment_length <= len
                } else {
                    pos - fragment_length >= 0
                }
            })
        })
    })
}

pub(crate) fn apply_strand_filter(
    index: &BifrostIndex,
    read_ec: &ReadEc,
    mode: StrandSpecific,
    is_first_read: bool,
    comprehensive: bool,
    report: &mut Option<&mut DebugReport>,
    header: &[u8],
) -> Option<Vec<u32>> {
    let ec = &read_ec.ec;
    let target = match mode {
        StrandSpecific::FR => is_first_read,
        StrandSpecific::RF => !is_first_read,
    };

    if comprehensive {
        let mut union = Vec::new();
        for hit in &read_ec.hits {
            if let Some(filtered) = filter_ec_for_hit(index, ec, hit, target)
                && !filtered.is_empty()
            {
                merge_sorted_unique_vec(&mut union, &filtered);
            }
        }
        return Some(union);
    }

    let hit = read_ec.first_hit.as_ref()?;
    let filtered = filter_ec_for_hit(index, ec, hit, target)?;
    if filtered.len() < ec.len()
        && let Some(r) = report.as_deref_mut()
    {
        r.record(
            header,
            DebugFailReason::Unknown,
            None,
            None,
            None,
            None,
            None,
            None,
            read_ec
                .best_match
                .as_ref()
                .map(|m| m.used_revcomp)
                .unwrap_or(false),
        );
    }
    Some(filtered)
}

fn filter_ec_for_hit(
    index: &BifrostIndex,
    ec: &[u32],
    hit: &Hit,
    target: bool,
) -> Option<Vec<u32>> {
    let mut filtered = Vec::new();
    let block = index.flat_ec.owned_ec_blocks.get(hit.unitig_id)?;
    let block = block.get(hit.block_idx)?;
    let positions = block.positions.as_ref()?;
    let um_strand = !hit.used_revcomp;
    for &tr in ec {
        let idx = match block.ec.binary_search(&tr) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let sense = positions
            .get(idx)
            .and_then(|positions| positions.strand_summary())
            .unwrap_or(2);
        if sense == 2 || ((um_strand == (sense == 1)) == target) {
            filtered.push(tr);
        }
    }
    // An empty result rejects this strand; None means metadata was unavailable.
    Some(filtered)
}

// Ordinary strand filtering consumes one orientation hit, not a full hit list.
// Keep index-backed ECs borrowed until a transcript is actually removed.
pub(crate) fn filter_ec_by_strand_cow<'a>(
    index: &BifrostIndex,
    ec: Cow<'a, [u32]>,
    hit: &Hit,
    mode: StrandSpecific,
    is_first_read: bool,
) -> Cow<'a, [u32]> {
    let Some(block) = index
        .flat_ec
        .owned_ec_blocks
        .get(hit.unitig_id)
        .and_then(|blocks| blocks.get(hit.block_idx))
    else {
        return ec;
    };
    let Some(positions) = block.positions.as_ref() else {
        return ec;
    };
    let target = match mode {
        StrandSpecific::FR => is_first_read,
        StrandSpecific::RF => !is_first_read,
    };
    // The common one-block answer is the annotation EC itself. Its ordinals
    // already identify the matching positions; no per-transcript search is needed.
    let same_ec = std::ptr::eq(ec.as_ref(), block.ec.as_ref());
    if same_ec
        && positions.len() == ec.len()
        && let Some(sense) = positions.uniform_singleton_strand()
    {
        if ((!hit.used_revcomp) == (sense == 1)) == target {
            return ec;
        }
        return match ec {
            Cow::Borrowed(values) => Cow::Borrowed(&values[..0]),
            Cow::Owned(mut values) => {
                values.clear();
                Cow::Owned(values)
            }
        };
    }

    let mut ordinal = 0;
    let mut keep = |tr: u32| {
        let idx = if same_ec {
            let idx = ordinal;
            ordinal += 1;
            idx
        } else {
            let Ok(idx) = block.ec.binary_search(&tr) else {
                return false;
            };
            idx
        };
        let sense = positions.strand_summary_at(idx).unwrap_or(2);
        sense == 2 || (((!hit.used_revcomp) == (sense == 1)) == target)
    };
    match ec {
        Cow::Owned(mut values) => {
            values.retain(|&tr| keep(tr));
            Cow::Owned(values)
        }
        Cow::Borrowed(values) => {
            let Some(first) = values.iter().position(|&tr| !keep(tr)) else {
                return Cow::Borrowed(values);
            };
            let mut filtered = Vec::with_capacity(values.len().saturating_sub(1));
            filtered.extend_from_slice(&values[..first]);
            filtered.extend(values[first + 1..].iter().copied().filter(|&tr| keep(tr)));
            Cow::Owned(filtered)
        }
    }
}

// Single-end fragment and strand filters commonly use the same annotation
// block. Decode each target's coordinates once while keeping both predicates.
pub(crate) fn filter_ec_by_fragment_and_strand_cow<'a>(
    index: &BifrostIndex,
    ec: Cow<'a, [u32]>,
    best: MatchInfo,
    fragment_length: i64,
    hit: &Hit,
    mode: StrandSpecific,
) -> Cow<'a, [u32]> {
    let shared = index
        .flat_ec
        .owned_ec_blocks
        .get(best.unitig_id)
        .filter(|blocks| {
            fragment_length > 0
                && !index.transcript_lengths.is_empty()
                && blocks
                    .first()
                    .is_some_and(|block| block.positions.is_some())
        })
        .and_then(|blocks| {
            let block = ec_blocks_leading_vals(blocks, best.unitig_pos as u32).last()?;
            let strand_block = index
                .flat_ec
                .owned_ec_blocks
                .get(hit.unitig_id)?
                .get(hit.block_idx)?;
            std::ptr::eq(block, strand_block).then_some((blocks, block))
        });
    let Some((blocks, block)) = shared else {
        let ec = filter_ec_by_fragment_cow(index, ec, best, fragment_length);
        return filter_ec_by_strand_cow(index, ec, hit, mode, true);
    };
    let Some(positions) = block.positions.as_ref() else {
        let ec = filter_ec_by_fragment_cow(index, ec, best, fragment_length);
        return filter_ec_by_strand_cow(index, ec, hit, mode, true);
    };
    let unitig_len = if best.unitig_id < index.unitig_count() {
        index
            .encoded_unitigs
            .unitig_bases
            .len(best.unitig_id)
            .unwrap_or(index.k)
    } else {
        index.k
    };
    let target = matches!(mode, StrandSpecific::FR);
    let same_ec = std::ptr::eq(ec.as_ref(), block.ec.as_ref());
    let mut ordinal = 0;
    let mut keep = |tr: u32| {
        let idx = if same_ec {
            let idx = ordinal;
            ordinal += 1;
            idx
        } else {
            let Ok(idx) = block.ec.binary_search(&tr) else {
                return false;
            };
            idx
        };
        let Some(raw_positions) = positions.get(idx) else {
            return false;
        };
        let sense = raw_positions.strand_summary().unwrap_or(2);
        if sense != 2 && (((!hit.used_revcomp) == (sense == 1)) != target) {
            return false;
        }
        let len = index
            .transcript_lengths
            .get(tr as usize)
            .copied()
            .unwrap_or(0) as i64;
        raw_positions.any(|raw| {
            map_transcript_position(
                blocks,
                tr,
                raw,
                best.unitig_pos,
                best.read_pos,
                best.used_revcomp,
                unitig_len,
                index.k,
            )
            .is_some_and(|(pos, forward)| {
                if forward {
                    pos + fragment_length <= len
                } else {
                    pos - fragment_length >= 0
                }
            })
        })
    };
    match ec {
        Cow::Owned(mut values) => {
            values.retain(|&tr| keep(tr));
            Cow::Owned(values)
        }
        Cow::Borrowed(values) => {
            let Some(first) = values.iter().position(|&tr| !keep(tr)) else {
                return Cow::Borrowed(values);
            };
            let mut filtered = Vec::with_capacity(values.len().saturating_sub(1));
            filtered.extend_from_slice(&values[..first]);
            filtered.extend(values[first + 1..].iter().copied().filter(|&tr| keep(tr)));
            Cow::Owned(filtered)
        }
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn map_transcript_position(
    blocks: &[crate::index::EcBlock],
    tr: u32,
    rawpos: u32,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
    unitig_len: usize,
    k: usize,
) -> Option<(i64, bool)> {
    const POSITION_MASK: u32 = 0x7fff_ffff;
    let trpos_raw = rawpos & POSITION_MASK;
    if trpos_raw == POSITION_MASK {
        return None;
    }
    let trpos = trpos_raw as i64;
    let trsense = rawpos == trpos_raw;

    let csense = !used_revcomp;
    let um_dist = unitig_pos as i64;
    let um_size = unitig_len as i64;
    let p = read_pos as i64;
    let k = k as i64;

    if trsense {
        let mut padding = 0i64;
        if trpos == 0 && blocks.len() > 1 {
            let mut mc_cur = ec_block_at(blocks, unitig_pos as u32)?;
            let ecs = ec_blocks_leading_vals(blocks, unitig_pos as u32);
            for block in ecs.iter().rev().skip(1) {
                if !block_contains(block, tr) || !block_contains_rawpos(block, tr, trpos_raw) {
                    padding = mc_cur.0 as i64;
                    break;
                }
                mc_cur = prev_block_at(blocks, mc_cur.0)?;
            }
        }
        let mut pos = trpos + um_dist + 1 - padding;
        if csense {
            pos -= p;
        } else {
            pos += k - 1 + p;
        }
        Some((pos, csense))
    } else {
        let mut r_end = um_size - k;
        if trpos == 0 && blocks.len() > 1 {
            let mc = ec_block_at(blocks, unitig_pos as u32)?;
            let ecs = ec_blocks_trailing_vals(blocks, mc.1);
            let mut mc_cur = ec_block_at(blocks, mc.1)?;
            for block in ecs {
                if !block_contains(block, tr) || !block_contains_rawpos(block, tr, rawpos) {
                    r_end = mc_cur.0 as i64 - 1;
                    break;
                }
                mc_cur = ec_block_at(blocks, mc_cur.1)?;
            }
        }
        let mut pos = trpos + (r_end - um_dist) + 1;
        if csense {
            pos += k - 1 + p;
        } else {
            pos -= p;
        }
        Some((pos, !csense))
    }
}

pub(super) fn ec_blocks_leading_vals(
    blocks: &[crate::index::EcBlock],
    idx: u32,
) -> &[crate::index::EcBlock] {
    if blocks.is_empty() {
        return &[];
    }
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= idx {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    &blocks[..lo]
}

pub(super) fn ec_blocks_trailing_vals(
    blocks: &[crate::index::EcBlock],
    idx: u32,
) -> &[crate::index::EcBlock] {
    if blocks.is_empty() {
        return &[];
    }
    if blocks.len() == 1 {
        return if blocks[0].lb >= idx {
            &blocks[..1]
        } else {
            &[]
        };
    }
    &blocks[blocks.partition_point(|block| block.lb < idx)..]
}

pub(crate) fn block_index_for_position(
    blocks: &[crate::index::EcBlock],
    pos: usize,
) -> Option<usize> {
    if blocks.is_empty() {
        return None;
    }
    if blocks.len() == 1 {
        return Some(0);
    }
    let pos = pos as u32;
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= pos {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        return None;
    }
    Some(lo - 1)
}

pub(super) fn ec_block_at(blocks: &[crate::index::EcBlock], idx: u32) -> Option<(u32, u32)> {
    if blocks.is_empty() {
        return None;
    }
    if blocks.len() == 1 {
        return Some((blocks[0].lb, blocks[0].ub));
    }
    let mut lo = 0usize;
    let mut hi = blocks.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if blocks[mid].lb <= idx {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        return None;
    }
    let block = &blocks[lo - 1];
    Some((block.lb, block.ub))
}

pub(super) fn prev_block_at(blocks: &[crate::index::EcBlock], lb: u32) -> Option<(u32, u32)> {
    if lb == 0 {
        ec_block_at(blocks, u32::MAX)
    } else {
        ec_block_at(blocks, lb - 1)
    }
}

pub(super) fn block_contains(block: &crate::index::EcBlock, tr: u32) -> bool {
    block.ec.binary_search(&tr).is_ok()
}

fn block_positions(
    block: &crate::index::EcBlock,
    tr: u32,
) -> Option<crate::index::PositionSetRef<'_>> {
    let positions = block.positions.as_ref()?;
    let idx = block.ec.binary_search(&tr).ok()?;
    positions.get(idx)
}

fn block_contains_rawpos(block: &crate::index::EcBlock, tr: u32, rawpos: u32) -> bool {
    block_positions(block, tr).is_some_and(|positions| positions.contains(rawpos))
}
