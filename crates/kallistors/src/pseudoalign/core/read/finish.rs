//! Final hit folding, shade/on-list filtering and general fallback results.

use crate::pseudoalign::core::local::{
    collect_unitigs_for_read, ec_for_read_kmer_index, ec_for_read_local_kmer,
};
use crate::pseudoalign::core::runtime::{ec_slice, onlist_cardinality};
use crate::pseudoalign::debug::ReadDebugState;
use crate::pseudoalign::ec::{
    ec_has_offlist, filter_shades, intersect_sorted, merge_sorted_unique_vec,
};
use crate::pseudoalign::{BifrostIndex, Hit, MatchInfo, PseudoalignOptions, ReadEc};

pub(super) struct ScannedRead {
    pub(super) shade_scratch: Vec<u32>,
    pub(super) best_match: Option<MatchInfo>,
    pub(super) first_hit: Option<Hit>,
    pub(super) hits: Vec<Hit>,
    pub(super) has_hit: bool,
    pub(super) saw_special_hit: bool,
    pub(super) saw_non_special_hit: bool,
    pub(super) append_dummy_hit: bool,
}

pub(super) fn finish_read(
    index: &BifrostIndex,
    seq: &[u8],
    mut dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
    scanned: ScannedRead,
) -> Option<ReadEc> {
    let ScannedRead {
        mut shade_scratch,
        best_match,
        first_hit,
        hits,
        has_hit,
        saw_special_hit,
        saw_non_special_hit,
        append_dummy_hit,
    } = scanned;
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let use_shade = index.use_shade && !options.do_union;
    let dlist_onlist = index.onlist.as_deref();
    if options.discard_special_only && saw_special_hit && !saw_non_special_hit {
        if let Some(state) = dbg.as_deref_mut() {
            state.special_only = true;
        }
        return None;
    }

    if !has_hit || hits.is_empty() {
        if (options.kallisto_local_fallback || options.kallisto_bifrost_find)
            && !options.kallisto_strict
        {
            let unitigs = collect_unitigs_for_read(
                index,
                seq,
                options.kallisto_enum,
                options.kallisto_strict,
                options.skip_overcrowded_minimizer,
                options.kallisto_bifrost_find,
                256,
            );
            if !unitigs.is_empty()
                && let Some(ec) =
                    ec_for_read_local_kmer(index, seq, &unitigs, true, options.discard_special_only)
            {
                if let Some(state) = dbg {
                    state.saw_valid_kmer = true;
                    state.saw_ec = true;
                }
                return Some(ReadEc {
                    ec,
                    best_match: None,
                    first_hit: None,
                    hits: Vec::new(),
                    had_offlist: false,
                    shade_union: Vec::new(),
                    hard_reject_pair: false,
                });
            }
        }
        if options.kallisto_fallback
            && let Some(kmer_index) = index.kmer_index.as_ref()
            && let Some(ec) = ec_for_read_kmer_index(kmer_index, seq)
        {
            if let Some(state) = dbg.as_deref_mut() {
                state.saw_valid_kmer = true;
                state.saw_ec = true;
            }
            return Some(ReadEc {
                ec,
                best_match: None,
                first_hit: None,
                hits: Vec::new(),
                had_offlist: false,
                shade_union: Vec::new(),
                hard_reject_pair: false,
            });
        }
        return None;
    }

    let mut minpos = usize::MAX;
    let mut maxpos = 0usize;
    let mut last_unitig = None;
    let mut last_block_idx = None;
    let mut found_nonempty = false;
    let mut current_offlist = false;
    let mut ec_scratch: Vec<u32> = Vec::new();
    for hit in &hits {
        minpos = minpos.min(hit.read_pos);
        maxpos = maxpos.max(hit.read_pos);
        let ec = ec_slice(index, hit.unitig_id, hit.block_idx);
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
                    last_unitig = Some(hit.unitig_id);
                    last_block_idx = Some(hit.block_idx);
                }
                current_offlist = ec_offlist;
                found_nonempty = true;
            }
            continue;
        }
        if last_unitig == Some(hit.unitig_id) && last_block_idx == Some(hit.block_idx) {
            continue;
        }
        if ec.is_empty() {
            continue;
        }
        if options.do_union {
            if options.dfk_onlist
                && (current_offlist || ec_offlist)
                && let Some(dummy) = index.onlist.as_ref().map(|v| onlist_cardinality(v))
                && current.last().copied() != Some(dummy)
            {
                current.push(dummy);
            }
            merge_sorted_unique_vec(&mut current, ec);
            current_offlist = current_offlist || ec_offlist;
        } else {
            if current.is_empty() {
                current.extend_from_slice(ec);
                last_unitig = Some(hit.unitig_id);
                last_block_idx = Some(hit.block_idx);
                current_offlist = ec_offlist;
                continue;
            }
            let ec_slice = if options.dfk_onlist && (current_offlist || ec_offlist) {
                if let Some(onlist) = index.onlist.as_deref() {
                    let dummy = onlist_cardinality(onlist);
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
                if let Some(state) = dbg {
                    state.intersection_empty = true;
                }
                return Some(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: current_offlist,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                });
            }
            last_unitig = Some(hit.unitig_id);
            last_block_idx = Some(hit.block_idx);
            current_offlist = current_offlist || ec_offlist;
        }
    }

    // Apply synthetic D-list dummy hits (kallisto um_dummy behavior). They must
    // participate in the fold before on-list filtering because they can reject it.
    if append_dummy_hit && let Some(onlist) = dlist_onlist {
        let dummy = onlist_cardinality(onlist);
        let ec = [dummy];
        if options.do_union {
            merge_sorted_unique_vec(&mut current, &ec);
            current_offlist = true;
        } else if current.is_empty() {
            current.extend_from_slice(&ec);
            current_offlist = true;
        } else {
            next.clear();
            intersect_sorted(&current, &ec, &mut next);
            std::mem::swap(&mut current, &mut next);
            if current.is_empty() {
                if let Some(state) = dbg {
                    state.intersection_empty = true;
                }
                return Some(ReadEc {
                    ec: Vec::new(),
                    best_match,
                    first_hit,
                    hits,
                    had_offlist: current_offlist,
                    shade_union: Vec::new(),
                    hard_reject_pair: true,
                });
            }
            current_offlist = true;
        }
    }

    if current.is_empty() {
        return None;
    }
    let min_range = options.min_range.max(1);
    // Coverage spans the leftmost hit's start through the rightmost hit's end,
    // not the number of probes; successful jumps deliberately skip probes.
    if maxpos >= minpos && (maxpos - minpos + index.k) < min_range {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
        if options.dfk_onlist && current_offlist && !current.is_empty() {
            let dummy = onlist_cardinality(onlist);
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
                let ec = ec_slice(index, hit.unitig_id, hit.block_idx);
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
            hard_reject_pair: false,
        })
    }
}
