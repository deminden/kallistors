use super::super::BifrostIndex;
use super::super::ec::encode_kmer_pair;
use super::super::minimizers::{
    fill_revcomp, minhash_next_after_hash, minhash_primary_for_kmer, minimizer_for_kmer_strict,
    minimizer_tail_for_kmer, minimizers_for_kmer, minimizers_ranked_for_kmer,
};
use super::runtime::block_index_for_position_fast;

// Shared k-mer matching primitives used by both the fast and parity-oriented scanners.
#[allow(clippy::too_many_arguments)]
pub(crate) fn match_unitig_candidate(
    unitig: &[u8],
    k: usize,
    rel_pos: usize,
    diff: usize,
    min_pos_fwd: usize,
    min_pos_rev: usize,
    kmer: &[u8],
    rev_kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
    allow_relaxed: bool,
) -> Option<(usize, bool, bool, bool)> {
    if unitig.len() < k {
        return None;
    }

    let rel = rel_pos as isize;
    let start_fwd = rel - min_pos_fwd as isize;
    if allow_forward && start_fwd >= 0 {
        let start = start_fwd as usize;
        if start + k <= unitig.len() && &unitig[start..start + k] == kmer {
            return Some((start, false, false, true));
        }
    }

    let start_rev = rel - diff as isize + min_pos_rev as isize;
    if allow_rev && start_rev >= 0 {
        let start = start_rev as usize;
        if start + k <= unitig.len() && &unitig[start..start + k] == rev_kmer {
            return Some((start, true, false, false));
        }
    }

    if !allow_relaxed {
        return None;
    }

    // Relaxed fallback: allow any k-mer start whose window can contain
    // the selected minimizer occurrence (rel_pos - (k-g) .. rel_pos).
    let max_start = unitig.len() - k;
    let start_lo = rel_pos.saturating_sub(diff).min(max_start);
    let start_hi = rel_pos.min(max_start);
    if start_lo > start_hi {
        return None;
    }
    for start in start_lo..=start_hi {
        if allow_forward && &unitig[start..start + k] == kmer {
            return Some((start, false, true, true));
        }
        if allow_rev && &unitig[start..start + k] == rev_kmer {
            return Some((start, true, true, false));
        }
    }

    None
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn match_unitig_candidate_encoded(
    index: &BifrostIndex,
    unitig_id: usize,
    rel_pos: usize,
    diff: usize,
    min_pos_fwd: usize,
    min_pos_rev: usize,
    read_fwd: u64,
    read_rev: u64,
    allow_forward: bool,
    allow_rev: bool,
    allow_relaxed: bool,
) -> Option<(usize, bool, bool, bool)> {
    let unitig_len = index.encoded_unitigs.unitig_bases.len(unitig_id)?;
    if unitig_len < index.k {
        return None;
    }

    let rel = rel_pos as isize;
    let start_fwd = rel - min_pos_fwd as isize;
    if allow_forward && start_fwd >= 0 {
        let start = start_fwd as usize;
        if start + index.k <= unitig_len
            && index
                .encoded_unitigs
                .unitig_bases
                .encode_kmer_at(unitig_id, start, index.k)
                == Some(read_fwd)
        {
            return Some((start, false, false, true));
        }
    }

    let start_rev = rel - diff as isize + min_pos_rev as isize;
    if allow_rev && start_rev >= 0 {
        let start = start_rev as usize;
        if start + index.k <= unitig_len
            && index
                .encoded_unitigs
                .unitig_bases
                .encode_kmer_at(unitig_id, start, index.k)
                == Some(read_rev)
        {
            return Some((start, true, false, false));
        }
    }
    if !allow_relaxed {
        return None;
    }
    let max_start = unitig_len - index.k;
    let start_lo = rel_pos.saturating_sub(diff).min(max_start);
    let start_hi = rel_pos.min(max_start);
    if start_lo > start_hi {
        return None;
    }
    for start in start_lo..=start_hi {
        let Some(code) = index
            .encoded_unitigs
            .unitig_bases
            .encode_kmer_at(unitig_id, start, index.k)
        else {
            continue;
        };
        if allow_forward && code == read_fwd {
            return Some((start, false, true, true));
        }
        if allow_rev && code == read_rev {
            return Some((start, true, true, false));
        }
    }

    None
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn match_kmer_at_pos(
    index: &BifrostIndex,
    kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
    _diff: usize,
    rev_buf: &mut Vec<u8>,
    kallisto_direct_kmer: bool,
    kallisto_enum: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_bifrost_find: bool,
    allow_relaxed: bool,
) -> Option<(usize, usize, bool, usize, bool)> {
    if !fill_revcomp(kmer, rev_buf) {
        return None;
    }
    let use_revcomp = rev_buf.as_slice() < kmer;
    let kmer_canon = if use_revcomp {
        rev_buf.as_slice()
    } else {
        kmer
    };
    let min_input = if kallisto_bifrost_find || kallisto_strict {
        kmer_canon
    } else {
        kmer
    };
    if kallisto_direct_kmer
        && let Some((uid, start, used_revcomp)) =
            match_kmer_direct(index, kmer, allow_forward, allow_rev)
    {
        let block_idx = block_index_for_position_fast(index, uid, start)?;
        return Some((uid, start, used_revcomp, block_idx, false));
    }
    let (mut min_candidates, mut min_hash_current) = if kallisto_bifrost_find {
        let (min_hash, candidates) = minhash_primary_for_kmer(min_input, index.g)?;
        (candidates, Some(min_hash))
    } else {
        let candidates = if kallisto_strict {
            vec![minimizer_for_kmer_strict(min_input, index.g)?]
        } else if kallisto_enum {
            minimizers_ranked_for_kmer(min_input, index.g, 2)?
        } else {
            minimizers_for_kmer(min_input, index.g)?
        };
        (candidates, None)
    };
    if !kallisto_bifrost_find
        && !kallisto_strict
        && min_candidates.iter().any(|(min_bytes, _)| {
            index
                .mphf
                .lookup(min_bytes)
                .map(|min_idx| {
                    index.minz_positions[min_idx as usize]
                        .iter()
                        .any(|&pos_id| {
                            (pos_id >> 32) as u32 == u32::MAX && (pos_id & 0x8000_0000) != 0
                        })
                })
                .unwrap_or(false)
        })
    {
        return match_kmer_at_pos(
            index,
            kmer,
            allow_forward,
            allow_rev,
            _diff,
            rev_buf,
            kallisto_direct_kmer,
            kallisto_enum,
            kallisto_strict,
            skip_overcrowded_minimizer,
            true,
            allow_relaxed,
        );
    }
    let mut tried_tail_minimizer = false;
    'min_outer: loop {
        let mut request_next_min = false;
        for (min_bytes, min_pos) in min_candidates.iter().copied() {
            let min_pos_fwd = if use_revcomp {
                _diff.saturating_sub(min_pos)
            } else {
                min_pos
            };
            let min_pos_rev = if use_revcomp {
                min_pos
            } else {
                _diff.saturating_sub(min_pos)
            };
            let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
                continue;
            };
            let positions = &index.minz_positions[min_idx as usize];
            if positions.is_empty() {
                continue;
            }
            for &pos_id in positions {
                let unitig_id_raw = (pos_id >> 32) as u32;
                if unitig_id_raw == u32::MAX {
                    if kallisto_strict {
                        continue;
                    }
                    let overcrowded = (pos_id & 0x8000_0000) != 0;
                    if !kallisto_bifrost_find && overcrowded && skip_overcrowded_minimizer {
                        continue;
                    }
                    let special_uid = if kallisto_bifrost_find {
                        special_unitig_for_kmer_raw(index, kmer)
                    } else {
                        special_unitig_for_kmer(index, kmer)
                    };
                    if let Some(uid) = special_uid {
                        let used_revcomp = kmer_canon != kmer;
                        if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                            continue;
                        }
                        let block_idx = block_index_for_position_fast(index, uid, 0)?;
                        return Some((uid, 0usize, used_revcomp, block_idx, false));
                    }
                    if kallisto_bifrost_find && overcrowded {
                        request_next_min = true;
                    }
                    continue;
                }
                let (unitig_id, rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
                if is_km {
                    let uid = unitig_id as usize;
                    if uid >= index.km_unitigs.len() {
                        continue;
                    }
                    let km_pos = rel_pos as usize;
                    let km_seq = index.km_unitigs[uid].as_slice();
                    if km_seq != kmer_canon {
                        continue;
                    }
                    if min_pos != km_pos && min_pos + km_pos != _diff {
                        continue;
                    }
                    let used_revcomp = kmer_canon != kmer;
                    if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                        continue;
                    }
                    let unitig_id = index.unitigs.len() + uid;
                    let block_idx = block_index_for_position_fast(index, unitig_id, 0)?;
                    return Some((unitig_id, 0usize, used_revcomp, block_idx, false));
                } else {
                    let uid = unitig_id as usize;
                    if uid >= index.unitigs.len() {
                        continue;
                    }
                    if let Some((start, used_revcomp, matched_relaxed, _forward_strand)) =
                        match_unitig_candidate(
                            index.unitigs[uid].as_slice(),
                            index.k,
                            rel_pos as usize,
                            _diff,
                            min_pos_fwd,
                            min_pos_rev,
                            kmer,
                            rev_buf.as_slice(),
                            allow_forward,
                            allow_rev,
                            allow_relaxed,
                        )
                    {
                        let block_idx = block_index_for_position_fast(index, uid, start)?;
                        return Some((uid, start, used_revcomp, block_idx, matched_relaxed));
                    }
                }
            }
        }
        if !kallisto_bifrost_find
            && !kallisto_strict
            && !kallisto_enum
            && !tried_tail_minimizer
            && let Some((tail_bytes, tail_pos)) = minimizer_tail_for_kmer(min_input, index.g)
            && !min_candidates
                .iter()
                .any(|&(min_bytes, min_pos)| min_bytes == tail_bytes && min_pos == tail_pos)
        {
            tried_tail_minimizer = true;
            min_candidates = vec![(tail_bytes, tail_pos)];
            continue 'min_outer;
        }
        if !kallisto_bifrost_find {
            break;
        }
        if request_next_min
            && let Some(curr_hash) = min_hash_current
            && let Some((next_hash, next_min)) =
                minhash_next_after_hash(min_input, index.g, curr_hash)
        {
            min_hash_current = Some(next_hash);
            min_candidates = vec![next_min];
            continue 'min_outer;
        }
        break;
    }
    None
}

pub(crate) fn match_kmer_direct(
    index: &BifrostIndex,
    kmer: &[u8],
    allow_forward: bool,
    allow_rev: bool,
) -> Option<(usize, usize, bool)> {
    let (fwd, rev) = encode_kmer_pair(kmer)?;
    if let Some(map) = index.kmer_pos_index.as_ref() {
        if let Some(entries) = map.get(&fwd) {
            for &(unitig_id, pos, used_revcomp) in entries {
                if used_revcomp {
                    continue;
                }
                if allow_forward {
                    return Some((unitig_id, pos, false));
                }
            }
        }
        if let Some(entries) = map.get(&rev) {
            for &(unitig_id, pos, used_revcomp) in entries {
                if !used_revcomp {
                    continue;
                }
                if allow_rev {
                    return Some((unitig_id, pos, true));
                }
            }
        }
    }

    let code = fwd.min(rev);
    for &(unitig_id, seq_fwd) in index.selective_fallback_kmers.lookup(code)? {
        if allow_forward && seq_fwd == fwd {
            return Some((unitig_id, 0, false));
        }
        if allow_rev && seq_fwd == rev {
            return Some((unitig_id, 0, true));
        }
    }
    None
}

pub(crate) fn special_unitig_for_kmer(index: &BifrostIndex, kmer: &[u8]) -> Option<usize> {
    let map = index.h_kmer_map.as_ref()?;
    let (fwd, rev) = encode_kmer_pair(kmer)?;
    let code = fwd.min(rev);
    map.get(&code).copied()
}

pub(crate) fn special_unitig_for_kmer_raw(index: &BifrostIndex, kmer: &[u8]) -> Option<usize> {
    let map = index.h_kmer_map.as_ref()?;
    let (fwd, rev) = encode_kmer_pair(kmer)?;
    let code = fwd.min(rev);
    map.get(&code).copied()
}

pub(crate) fn jump_distance_for_match(
    index: &BifrostIndex,
    unitig_id: usize,
    block_idx: usize,
    start: usize,
    forward_strand: bool,
) -> Option<usize> {
    let blocks = index.ec_blocks.get(unitig_id)?;
    let block = blocks.get(block_idx)?;
    if block.ub <= block.lb {
        return None;
    }
    let contig_start = block.lb as isize;
    let contig_len = (block.ub - block.lb) as isize;
    let start = start as isize;
    let dist = if forward_strand {
        contig_len - 1 - (start - contig_start)
    } else {
        start - contig_start
    };
    if dist >= 2 { Some(dist as usize) } else { None }
}
