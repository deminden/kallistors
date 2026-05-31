use super::super::ec::merge_sorted_unique_vec;
use super::super::{
    BifrostIndex, DebugFailReason, DebugReport, Hit, MatchInfo, ReadEc, StrandSpecific,
};

// Fragment and strand filters interpret compact EC blocks back into transcript coordinates.
pub(crate) fn filter_ec_by_fragment(
    index: &BifrostIndex,
    ec: &[u32],
    best_match: MatchInfo,
    fragment_length: i64,
) -> Vec<u32> {
    if fragment_length <= 0 {
        return ec.to_vec();
    }
    if index.transcript_lengths.is_empty() {
        return ec.to_vec();
    }
    if best_match.unitig_id >= index.ec_blocks.len() {
        return ec.to_vec();
    }
    let blocks = &index.ec_blocks[best_match.unitig_id];
    if blocks.is_empty() || blocks[0].positions.is_none() {
        return ec.to_vec();
    }
    let unitig_len = if best_match.unitig_id < index.unitigs.len() {
        index.unitigs[best_match.unitig_id].len()
    } else {
        index.k
    };
    let mut filtered = Vec::new();
    for &tr in ec {
        let positions = match find_positions_in_transcript(
            blocks,
            tr,
            best_match.unitig_pos,
            best_match.read_pos,
            best_match.used_revcomp,
            unitig_len,
            index.k,
        ) {
            Some(positions) => positions,
            None => continue,
        };
        let len = index
            .transcript_lengths
            .get(tr as usize)
            .copied()
            .unwrap_or(0) as i64;
        for (pos, forward) in positions {
            if forward && pos + fragment_length <= len {
                filtered.push(tr);
                break;
            }
            if !forward && pos - fragment_length >= 0 {
                filtered.push(tr);
                break;
            }
        }
    }
    filtered
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
    let block = index.ec_blocks.get(hit.unitig_id)?;
    let block = block.get(hit.block_idx)?;
    let strands = block.strands.as_ref()?;
    let um_strand = !hit.used_revcomp;
    for &tr in ec {
        let idx = match block.ec.binary_search(&tr) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let sense = strands.get(idx).copied().unwrap_or(2);
        if sense == 2 || ((um_strand == (sense == 1)) == target) {
            filtered.push(tr);
        }
    }
    if filtered.is_empty() {
        None
    } else {
        Some(filtered)
    }
}

pub(super) fn find_positions_in_transcript(
    blocks: &[crate::index::EcBlock],
    tr: u32,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
    unitig_len: usize,
    k: usize,
) -> Option<Vec<(i64, bool)>> {
    let idx = unitig_pos as u32;
    let ecs = ec_blocks_leading_vals(blocks, idx);
    if ecs.is_empty() {
        return None;
    }
    let v_ec = ecs.last().copied()?;
    let raw_positions = block_positions(v_ec, tr)?;
    let mut positions = Vec::new();
    raw_positions.for_each(|rawpos| {
        if let Some(mapped) = map_transcript_position(
            blocks,
            tr,
            rawpos,
            unitig_pos,
            read_pos,
            used_revcomp,
            unitig_len,
            k,
        ) {
            positions.push(mapped);
        }
    });
    Some(positions)
}

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
) -> Vec<&crate::index::EcBlock> {
    if blocks.is_empty() {
        return Vec::new();
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
    blocks[..lo].iter().collect()
}

pub(super) fn ec_blocks_trailing_vals(
    blocks: &[crate::index::EcBlock],
    idx: u32,
) -> Vec<&crate::index::EcBlock> {
    if blocks.is_empty() {
        return Vec::new();
    }
    if blocks.len() == 1 {
        return if blocks[0].lb >= idx {
            vec![&blocks[0]]
        } else {
            Vec::new()
        };
    }
    blocks.iter().filter(|block| block.lb >= idx).collect()
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

fn block_positions(block: &crate::index::EcBlock, tr: u32) -> Option<&crate::index::PositionSet> {
    let positions = block.positions.as_ref()?;
    let idx = block.ec.binary_search(&tr).ok()?;
    positions.get(idx)
}

fn block_contains_rawpos(block: &crate::index::EcBlock, tr: u32, rawpos: u32) -> bool {
    block_positions(block, tr).is_some_and(|positions| positions.contains(rawpos))
}
