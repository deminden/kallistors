use super::debug::{ReadDebugState, format_positions_sample};
use super::minimizers::{fill_revcomp, minimizers_for_kmer};
use super::{BifrostIndex, PseudoalignOptions};
use super::{DebugFailReason, DebugReport, Hit, MatchInfo, ReadEc, Strand, StrandSpecific};
use super::{ec_has_offlist, filter_shades, intersect_sorted, merge_sorted_unique_vec};

pub(super) fn ec_for_read_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    if seq.len() < index.k {
        if let Some(state) = dbg {
            state.saw_valid_kmer = false;
        }
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut rev_buf: Vec<u8> = Vec::new();
    let allow_forward = strand != Strand::Reverse;
    let allow_rev = strand != Strand::Forward;
    let mut has_hit = false;
    let diff = index.k.saturating_sub(index.g);
    let mut best_match: Option<MatchInfo> = None;
    let mut first_hit: Option<Hit> = None;
    let mut hits: Vec<Hit> = Vec::new();

    let mut pos = 0usize;
    let last_pos = seq.len() - index.k;
    while pos <= last_pos {
        let kmer = &seq[pos..pos + index.k];
        let mut matched = None;
        let min_candidates = match minimizers_for_kmer(kmer, index.g) {
            Some(v) => {
                if let Some(state) = dbg.as_deref_mut() {
                    state.saw_valid_kmer = true;
                }
                v
            }
            None => {
                pos += 1;
                continue;
            }
        };
        let mut any_mphf = false;
        let mut any_positions = false;
        for (min_bytes, min_pos) in min_candidates {
            let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
                if let Some(state) = dbg.as_deref_mut()
                    && state.first_mphf_miss.is_none()
                {
                    state.first_mphf_miss = Some((pos, min_pos));
                }
                continue;
            };
            any_mphf = true;
            let positions = &index.minz_positions[min_idx as usize];
            if positions.is_empty() {
                if let Some(state) = dbg.as_deref_mut()
                    && state.first_no_positions.is_none()
                {
                    state.first_no_positions = Some((pos, min_pos));
                }
                continue;
            }
            any_positions = true;
            if let Some(state) = dbg.as_deref_mut()
                && state.first_no_match_positions.is_none()
            {
                state.first_no_match_positions = Some((
                    pos,
                    min_pos,
                    positions.len(),
                    format_positions_sample(positions),
                ));
            }
            let mut local_match = None;
            for &pos_id in positions {
                let (unitig_id, rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
                if is_km {
                    let uid = unitig_id as usize;
                    if uid >= index.km_unitigs.len() {
                        continue;
                    }
                    let km_pos = rel_pos as usize;
                    let rev_match = km_pos <= diff && min_pos == diff - km_pos;
                    let fwd_match = min_pos == km_pos;
                    let has_rev = rev_match && allow_rev;
                    let has_fwd = fwd_match && allow_forward;
                    if !has_rev && !has_fwd {
                        continue;
                    }
                    if !fill_revcomp(kmer, &mut rev_buf) {
                        continue;
                    }
                    let canonical = if rev_buf.as_slice() < kmer {
                        rev_buf.as_slice()
                    } else {
                        kmer
                    };
                    if index.km_unitigs[uid] == canonical {
                        local_match =
                            Some((index.unitigs.len() + uid, 0usize, pos, min_pos, rev_match));
                        if rev_match && let Some(state) = dbg.as_deref_mut() {
                            state.used_revcomp = true;
                        }
                        break;
                    }
                } else {
                    let uid = unitig_id as usize;
                    if uid >= index.unitigs.len() {
                        continue;
                    }
                    let rel_pos = rel_pos as isize;
                    if allow_forward {
                        let start = rel_pos - min_pos as isize;
                        if start >= 0 {
                            let start = start as usize;
                            if start + index.k <= index.unitigs[uid].len()
                                && &index.unitigs[uid][start..start + index.k] == kmer
                            {
                                local_match = Some((uid, start, pos, min_pos, false));
                                break;
                            }
                        }
                    }
                    if allow_rev && fill_revcomp(kmer, &mut rev_buf) {
                        let start = rel_pos - diff as isize + min_pos as isize;
                        if start >= 0 {
                            let start = start as usize;
                            if start + index.k <= index.unitigs[uid].len()
                                && index.unitigs[uid][start..start + index.k] == rev_buf
                            {
                                local_match = Some((uid, start, pos, min_pos, true));
                                if let Some(state) = dbg.as_deref_mut() {
                                    state.used_revcomp = true;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            if let Some(hit) = local_match {
                matched = Some(hit);
                break;
            }
            if let Some(state) = dbg.as_deref_mut()
                && state.first_no_match.is_none()
            {
                state.first_no_match = Some((pos, min_pos));
            }
        }
        if let Some(state) = dbg.as_deref_mut() {
            if any_mphf {
                state.saw_mphf_hit = true;
            }
            if any_positions {
                state.saw_positions = true;
            }
        }

        let Some((uid, start, kmer_pos, min_pos, used_revcomp)) = matched else {
            pos += 1;
            continue;
        };
        if let Some(state) = dbg.as_deref_mut() {
            state.saw_match = true;
        }
        if best_match.map(|m| kmer_pos < m.read_pos).unwrap_or(true) {
            best_match = Some(MatchInfo {
                unitig_id: uid,
                unitig_pos: start,
                read_pos: kmer_pos,
                used_revcomp,
            });
        }
        let Some(block_idx) = block_index_for_position(&index.ec_blocks[uid], start) else {
            pos += 1;
            continue;
        };
        let ec = &index.ec_blocks[uid][block_idx].ec;
        if ec.is_empty() {
            if let Some(state) = dbg.as_deref_mut()
                && state.first_empty_ec.is_none()
            {
                state.first_empty_ec = Some((kmer_pos, min_pos));
            }
        } else if let Some(state) = dbg.as_deref_mut() {
            state.saw_ec = true;
        }
        hits.push(Hit {
            unitig_id: uid,
            read_pos: kmer_pos,
            block_idx,
            used_revcomp,
        });
        if first_hit
            .as_ref()
            .map(|h| kmer_pos < h.read_pos)
            .unwrap_or(true)
        {
            first_hit = Some(Hit {
                unitig_id: uid,
                read_pos: kmer_pos,
                block_idx,
                used_revcomp,
            });
        }
        has_hit = true;
        if !options.no_jump
            && let Some(dist) = jump_distance_for_match(index, uid, block_idx, start, used_revcomp)
        {
            let mut next_pos = pos + dist;
            if next_pos > last_pos {
                next_pos = last_pos;
            }
            if next_pos > pos {
                let kmer_next = &seq[next_pos..next_pos + index.k];
                if let Some((uid2, _start2, _rev2, block_idx2)) = match_kmer_at_pos(
                    index,
                    kmer_next,
                    allow_forward,
                    allow_rev,
                    diff,
                    &mut rev_buf,
                ) && uid2 == uid
                    && index.ec_blocks[uid][block_idx].ec == index.ec_blocks[uid2][block_idx2].ec
                {
                    pos = next_pos;
                    continue;
                }
            }
        }
        pos += 1;
    }

    if !has_hit || hits.is_empty() {
        return None;
    }

    hits.sort_by(|a, b| {
        a.unitig_id
            .cmp(&b.unitig_id)
            .then_with(|| a.read_pos.cmp(&b.read_pos))
    });

    let mut minpos = usize::MAX;
    let mut maxpos = 0usize;
    let mut last_ec: Vec<u32> = Vec::new();
    let mut last_unitig = None;
    let mut found_nonempty = false;
    let mut current_offlist = false;
    let use_shade = index.use_shade && !options.do_union;
    let mut shade_scratch: Vec<u32> = Vec::new();
    let mut ec_scratch: Vec<u32> = Vec::new();
    for hit in &hits {
        minpos = minpos.min(hit.read_pos);
        maxpos = maxpos.max(hit.read_pos);
        let ec = &index.ec_blocks[hit.unitig_id][hit.block_idx].ec;
        let ec = if use_shade {
            filter_shades(ec, &index.shade_sequences, &mut shade_scratch);
            shade_scratch.as_slice()
        } else {
            ec
        };
        let ec_offlist = if let Some(onlist) = index.onlist.as_deref() {
            ec_has_offlist(ec, onlist)
        } else {
            false
        };
        if !found_nonempty {
            if !ec.is_empty() {
                if options.do_union {
                    current.extend_from_slice(ec);
                    current.sort_unstable();
                    current.dedup();
                } else {
                    current.extend_from_slice(ec);
                    last_ec = ec.to_vec();
                    last_unitig = Some(hit.unitig_id);
                }
                current_offlist = ec_offlist;
                found_nonempty = true;
            }
            continue;
        }
        if last_unitig == Some(hit.unitig_id) && ec == last_ec.as_slice() {
            continue;
        }
        if ec.is_empty() {
            continue;
        }
        if options.do_union {
            if options.dfk_onlist
                && (current_offlist || ec_offlist)
                && let Some(dummy) = index.onlist.as_ref().map(|v| v.len() as u32)
                && current.last().copied() != Some(dummy)
            {
                current.push(dummy);
            }
            merge_sorted_unique_vec(&mut current, ec);
            current_offlist = current_offlist || ec_offlist;
        } else {
            if current.is_empty() {
                current.extend_from_slice(ec);
                last_ec = ec.to_vec();
                last_unitig = Some(hit.unitig_id);
                current_offlist = ec_offlist;
                continue;
            }
            let ec_slice = if options.dfk_onlist && (current_offlist || ec_offlist) {
                if let Some(onlist) = index.onlist.as_deref() {
                    let dummy = onlist.len() as u32;
                    ec_scratch.clear();
                    ec_scratch.extend_from_slice(ec);
                    if ec_scratch.last().copied() != Some(dummy) {
                        ec_scratch.push(dummy);
                    }
                    if current.last().copied() != Some(dummy) {
                        current.push(dummy);
                    }
                    ec_scratch.sort_unstable();
                    ec_scratch.dedup();
                    ec_scratch.as_slice()
                } else {
                    ec
                }
            } else {
                ec
            };
            next.clear();
            intersect_sorted(&current, ec_slice, &mut next);
            std::mem::swap(&mut current, &mut next);
            if current.is_empty() {
                if let Some(state) = dbg.as_deref_mut() {
                    state.intersection_empty = true;
                }
                return None;
            }
            last_ec = ec.to_vec();
            last_unitig = Some(hit.unitig_id);
            current_offlist = current_offlist || ec_offlist;
        }
    }

    if current.is_empty() {
        return None;
    }
    let min_range = options.min_range.max(1);
    if maxpos >= minpos && (maxpos - minpos + index.k) < min_range {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
        if options.dfk_onlist && current_offlist && !current.is_empty() {
            let dummy = onlist.len() as u32;
            if current.last().copied() != Some(dummy) {
                current.push(dummy);
            }
        }
    }
    if current.is_empty() {
        None
    } else {
        let mut shade_union = Vec::new();
        if index.use_shade {
            for hit in &hits {
                let ec = &index.ec_blocks[hit.unitig_id][hit.block_idx].ec;
                for &tr in ec {
                    let idx = tr as usize;
                    if idx < index.shade_sequences.len() && index.shade_sequences[idx] {
                        shade_union.push(tr);
                    }
                }
            }
            shade_union.sort_unstable();
            shade_union.dedup();
        }
        Some(ReadEc {
            ec: current,
            best_match,
            first_hit,
            hits,
            had_offlist: current_offlist,
            shade_union,
        })
    }
}

fn match_kmer_at_pos(
    index: &BifrostIndex,
    kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
    diff: usize,
    rev_buf: &mut Vec<u8>,
) -> Option<(usize, usize, bool, usize)> {
    let min_candidates = minimizers_for_kmer(kmer, index.g)?;
    for (min_bytes, min_pos) in min_candidates {
        let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
            continue;
        };
        let positions = &index.minz_positions[min_idx as usize];
        if positions.is_empty() {
            continue;
        }
        for &pos_id in positions {
            let (unitig_id, rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
            if is_km {
                let uid = unitig_id as usize;
                if uid >= index.km_unitigs.len() {
                    continue;
                }
                let km_pos = rel_pos as usize;
                let rev_match = km_pos <= diff && min_pos == diff - km_pos;
                let fwd_match = min_pos == km_pos;
                let has_rev = rev_match && allow_rev;
                let has_fwd = fwd_match && allow_forward;
                if !has_rev && !has_fwd {
                    continue;
                }
                if !fill_revcomp(kmer, rev_buf) {
                    continue;
                }
                let canonical = if rev_buf.as_slice() < kmer {
                    rev_buf.as_slice()
                } else {
                    kmer
                };
                if index.km_unitigs[uid] == canonical {
                    let unitig_id = index.unitigs.len() + uid;
                    let block_idx = block_index_for_position(&index.ec_blocks[unitig_id], 0)?;
                    return Some((unitig_id, 0usize, rev_match, block_idx));
                }
            } else {
                let uid = unitig_id as usize;
                if uid >= index.unitigs.len() {
                    continue;
                }
                let rel_pos = rel_pos as isize;
                if allow_forward {
                    let start = rel_pos - min_pos as isize;
                    if start >= 0 {
                        let start = start as usize;
                        if start + index.k <= index.unitigs[uid].len()
                            && &index.unitigs[uid][start..start + index.k] == kmer
                        {
                            let block_idx = block_index_for_position(&index.ec_blocks[uid], start)?;
                            return Some((uid, start, false, block_idx));
                        }
                    }
                }
                if allow_rev && fill_revcomp(kmer, rev_buf) {
                    let start = rel_pos - diff as isize + min_pos as isize;
                    if start >= 0 {
                        let start = start as usize;
                        if start + index.k <= index.unitigs[uid].len()
                            && &index.unitigs[uid][start..start + index.k] == rev_buf.as_slice()
                        {
                            let block_idx = block_index_for_position(&index.ec_blocks[uid], start)?;
                            return Some((uid, start, true, block_idx));
                        }
                    }
                }
            }
        }
    }
    None
}

fn jump_distance_for_match(
    index: &BifrostIndex,
    unitig_id: usize,
    block_idx: usize,
    start: usize,
    used_revcomp: bool,
) -> Option<usize> {
    let blocks = index.ec_blocks.get(unitig_id)?;
    let block = blocks.get(block_idx)?;
    if block.ub <= block.lb {
        return None;
    }
    let contig_start = block.lb as isize;
    let contig_len = (block.ub - block.lb) as isize;
    let start = start as isize;
    let forward = !used_revcomp;
    let dist = if forward {
        contig_len - 1 - (start - contig_start)
    } else {
        start - contig_start
    };
    if dist >= 2 { Some(dist as usize) } else { None }
}

pub(super) fn filter_ec_by_fragment(
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
        let (pos, forward) = match find_position_in_transcript(
            blocks,
            tr,
            best_match.unitig_pos,
            best_match.read_pos,
            best_match.used_revcomp,
            unitig_len,
            index.k,
        ) {
            Some((pos, forward)) => (pos, forward),
            None => {
                filtered.push(tr);
                continue;
            }
        };
        let len = index
            .transcript_lengths
            .get(tr as usize)
            .copied()
            .unwrap_or(0) as i64;
        if forward {
            if pos + fragment_length - 1 <= len {
                filtered.push(tr);
            }
        } else if pos - fragment_length >= 0 {
            filtered.push(tr);
        }
    }
    filtered
}

pub(super) fn apply_strand_filter(
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
                super::merge_sorted_unique_vec(&mut union, &filtered);
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

pub(super) fn find_position_in_transcript(
    blocks: &[crate::index::EcBlock],
    tr: u32,
    unitig_pos: usize,
    read_pos: usize,
    used_revcomp: bool,
    unitig_len: usize,
    k: usize,
) -> Option<(i64, bool)> {
    let idx = unitig_pos as u32;
    let mc = ec_block_at(blocks, idx)?;
    let ecs = ec_blocks_leading_vals(blocks, idx);
    if ecs.is_empty() {
        return None;
    }
    let v_ec = ecs.last().copied()?;
    let rawpos = block_min_pos(v_ec, tr)?;
    let trpos = (rawpos & 0x7fff_ffff) as i64;
    let trsense = rawpos == (rawpos & 0x7fff_ffff);

    let csense = !used_revcomp;
    let um_dist = unitig_pos as i64;
    let um_size = unitig_len as i64;
    let p = read_pos as i64;
    let k = k as i64;
    let mc_first = mc.0 as i64;
    let mc_second = mc.1 as i64;
    let _um_dist_block = um_dist - mc_first;

    if trsense {
        if csense {
            let mut padding = 0i64;
            if trpos == 0 {
                let mut mc_cur = mc;
                for block in ecs.iter().rev().skip(1) {
                    if !block_contains(block, tr) {
                        padding = mc_cur.0 as i64;
                        break;
                    }
                    mc_cur = prev_block_at(blocks, mc_cur.0)?;
                }
            }
            let pos = trpos - p + um_dist + 1 - padding;
            Some((pos, csense))
        } else {
            let mut mc_cur = mc;
            let mut right_one = 0i64;
            let mut left_one = 0i64;
            let initial = mc_second;
            for (i, block) in ecs.iter().enumerate().rev() {
                if i == ecs.len() - 1 {
                    right_one = mc_cur.1 as i64;
                }
                if !block_contains(block, tr) {
                    left_one = mc_cur.1 as i64;
                    break;
                } else if i == 0 {
                    left_one = 0;
                }
                mc_cur = prev_block_at(blocks, mc_cur.0)?;
            }
            let padding = -(left_one + right_one - um_size + k - 1);
            let pos = trpos + p + k - (um_size - k - um_dist) + initial - 1 + padding;
            Some((pos, csense))
        }
    } else if csense {
        let mut left_one = 0i64;
        let mut right_one = 0i64;
        let mut unmapped_len = 0i64;
        let mut found_first_mapped = false;
        let ecs_all = ec_blocks_leading_vals(blocks, u32::MAX);
        let mut curr_mc = 0u32;
        for block in ecs_all {
            let mc_ = ec_block_at(blocks, curr_mc)?;
            if !block_contains(block, tr) && found_first_mapped {
                if unmapped_len == 0 {
                    left_one = mc_.0 as i64;
                }
                right_one = mc_.1 as i64;
                unmapped_len += mc_.1 as i64 - mc_.0 as i64;
            }
            if block_contains(block, tr) {
                found_first_mapped = true;
            }
            curr_mc = mc_.1;
        }
        let mut start = 0i64;
        start -= right_one - left_one;
        start += um_size - k;
        let pos = trpos + (-(um_dist - start)) + k + p;
        Some((pos, !csense))
    } else {
        let mut left_one = 0i64;
        let mut right_one = 0i64;
        let mut unmapped_len = 0i64;
        let mut found_first_mapped = false;
        let ecs_all = ec_blocks_leading_vals(blocks, u32::MAX);
        let mut curr_mc = 0u32;
        for block in ecs_all {
            let mc_ = ec_block_at(blocks, curr_mc)?;
            if !block_contains(block, tr) && found_first_mapped {
                if unmapped_len == 0 {
                    left_one = mc_.0 as i64;
                }
                right_one = mc_.1 as i64;
                unmapped_len += mc_.1 as i64 - mc_.0 as i64;
            }
            if block_contains(block, tr) {
                found_first_mapped = true;
            }
            curr_mc = mc_.1;
        }
        let unmapped_len = right_one - left_one;
        let padding = um_size - um_dist - unmapped_len - k + 1;
        let pos = trpos + padding - p;
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

pub(super) fn block_index_for_position(
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

pub(super) fn block_min_pos(block: &crate::index::EcBlock, tr: u32) -> Option<u32> {
    let positions = block.positions.as_ref()?;
    let idx = block.ec.binary_search(&tr).ok()?;
    positions.get(idx).copied()
}
