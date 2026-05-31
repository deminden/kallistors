use std::collections::{HashMap, HashSet};

use super::super::ec::{encode_kmer_pair, intersect_sorted, lookup_ec, merge_sorted_unique};
use super::super::minimizers::{
    fill_revcomp, minhash_next_after_hash, minhash_primary_for_kmer, minimizer_for_kmer_strict,
    minimizers_for_kmer, minimizers_ranked_for_kmer,
};
use super::super::{BifrostIndex, KmerEcIndex};
use super::matcher::{special_unitig_for_kmer, special_unitig_for_kmer_raw};
use super::runtime::block_index_for_position_fast;

// Local k-mer fallback narrows expensive rescues to the unitigs already suggested by minimizers.
pub(crate) fn collect_unitigs_for_read(
    index: &BifrostIndex,
    seq: &[u8],
    kallisto_enum: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_bifrost_find: bool,
    max_unitigs: usize,
) -> Vec<usize> {
    if seq.len() < index.k {
        return Vec::new();
    }
    let mut out: HashSet<usize> = HashSet::new();
    let mut rev_buf: Vec<u8> = Vec::new();
    for pos in 0..=seq.len() - index.k {
        let kmer = &seq[pos..pos + index.k];
        if !fill_revcomp(kmer, &mut rev_buf) {
            continue;
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
        let (mut min_candidates, mut min_hash_current) = if kallisto_bifrost_find {
            let (min_hash, candidates) = match minhash_primary_for_kmer(min_input, index.g) {
                Some(v) => v,
                None => continue,
            };
            (candidates, Some(min_hash))
        } else {
            let candidates = if kallisto_strict {
                match minimizer_for_kmer_strict(min_input, index.g) {
                    Some(v) => vec![v],
                    None => continue,
                }
            } else if kallisto_enum {
                match minimizers_ranked_for_kmer(min_input, index.g, 2) {
                    Some(v) => v,
                    None => continue,
                }
            } else {
                match minimizers_for_kmer(min_input, index.g) {
                    Some(v) => v,
                    None => continue,
                }
            };
            (candidates, None)
        };
        'min_loop: loop {
            let mut request_next_min = false;
            for (min_bytes, _min_pos) in min_candidates.iter() {
                let Some(min_idx) = index.mphf.lookup(min_bytes) else {
                    continue;
                };
                let positions = &index.minz_positions[min_idx as usize];
                if positions.is_empty() {
                    continue;
                }
                for &pos_id in positions {
                    let unitig_id_raw = (pos_id >> 32) as u32;
                    if unitig_id_raw == u32::MAX {
                        if !kallisto_strict {
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
                                out.insert(uid);
                                if out.len() >= max_unitigs {
                                    return out.into_iter().collect();
                                }
                            }
                            if kallisto_bifrost_find && overcrowded {
                                request_next_min = true;
                            }
                        }
                        continue;
                    }
                    let (unitig_id, _rel_pos, is_km) = crate::index::bifrost::decode_pos_id(pos_id);
                    let uid = unitig_id as usize;
                    let real_id = if is_km {
                        index.unitigs.len() + uid
                    } else {
                        uid
                    };
                    out.insert(real_id);
                    if out.len() >= max_unitigs {
                        return out.into_iter().collect();
                    }
                }
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
                continue 'min_loop;
            }
            break;
        }
    }
    out.into_iter().collect()
}

fn build_local_kmer_map(index: &BifrostIndex, unitig_ids: &[usize]) -> HashMap<u64, Vec<u32>> {
    let mut map: HashMap<u64, Vec<u32>> = HashMap::new();
    for &unitig_id in unitig_ids {
        let seq_ref: Option<&[u8]> = if unitig_id < index.unitigs.len() {
            Some(index.unitigs[unitig_id].as_slice())
        } else {
            let km_idx = unitig_id - index.unitigs.len();
            index.km_unitigs.get(km_idx).map(|v| v.as_slice())
        };
        let Some(seq_ref) = seq_ref else { continue };
        let blocks = match index.ec_blocks.get(unitig_id) {
            Some(v) => v,
            None => continue,
        };
        for block in blocks {
            let start = block.lb as usize;
            let end = block.ub as usize;
            for pos in start..end {
                if pos + index.k > seq_ref.len() {
                    break;
                }
                if let Some((fwd, rev)) = encode_kmer_pair(&seq_ref[pos..pos + index.k]) {
                    map.entry(fwd)
                        .and_modify(|existing| merge_sorted_unique(existing, &block.ec))
                        .or_insert_with(|| block.ec.clone());
                    if rev != fwd {
                        map.entry(rev)
                            .and_modify(|existing| merge_sorted_unique(existing, &block.ec))
                            .or_insert_with(|| block.ec.clone());
                    }
                }
            }
        }
    }
    map
}

pub fn local_kmer_hits(
    index: &BifrostIndex,
    seq: &[u8],
    kallisto_enum: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_bifrost_find: bool,
    max_unitigs: usize,
) -> Vec<(usize, String, bool, usize)> {
    if seq.len() < index.k {
        return Vec::new();
    }
    let unitigs = collect_unitigs_for_read(
        index,
        seq,
        kallisto_enum,
        kallisto_strict,
        skip_overcrowded_minimizer,
        kallisto_bifrost_find,
        max_unitigs,
    );
    let map = build_local_kmer_map(index, &unitigs);
    let mut out = Vec::new();
    for pos in 0..=seq.len() - index.k {
        let kmer = &seq[pos..pos + index.k];
        let mut hit = false;
        let mut ec_size = 0usize;
        if let Some((fwd, rev)) = encode_kmer_pair(kmer)
            && let Some(ec) = map.get(&fwd).or_else(|| map.get(&rev))
        {
            hit = true;
            ec_size = ec.len();
        } else if !kallisto_strict && let Some(ec) = special_ec_for_kmer(index, kmer) {
            hit = true;
            ec_size = ec.len();
        }
        out.push((
            pos,
            String::from_utf8_lossy(kmer).into_owned(),
            hit,
            ec_size,
        ));
    }
    out
}

fn special_ec_for_kmer<'a>(index: &'a BifrostIndex, kmer: &'a [u8]) -> Option<&'a [u32]> {
    let uid = special_unitig_for_kmer(index, kmer)?;
    let _ = index.ec_blocks.get(uid)?;
    let block_idx = block_index_for_position_fast(index, uid, 0)?;
    Some(index.flat_ec.ec(uid, block_idx))
}

pub(crate) fn ec_for_read_local_kmer(
    index: &BifrostIndex,
    seq: &[u8],
    unitig_ids: &[usize],
    allow_special: bool,
    discard_special_only: bool,
) -> Option<Vec<u32>> {
    if seq.len() < index.k {
        return None;
    }
    let map = build_local_kmer_map(index, unitig_ids);
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut has_hit = false;
    let mut saw_non_special = false;
    let mut saw_special = false;
    for pos in 0..=seq.len() - index.k {
        if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + index.k])
            && let Some(ec) = map.get(&fwd).or_else(|| map.get(&rev))
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                intersect_sorted(&current, ec, &mut next);
                std::mem::swap(&mut current, &mut next);
                if current.is_empty() {
                    break;
                }
            }
            saw_non_special = true;
        } else if allow_special
            && let Some(ec) = special_ec_for_kmer(index, &seq[pos..pos + index.k])
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                intersect_sorted(&current, ec, &mut next);
                std::mem::swap(&mut current, &mut next);
                if current.is_empty() {
                    break;
                }
            }
            saw_special = true;
        }
    }
    if discard_special_only && saw_special && !saw_non_special {
        return None;
    }
    if !has_hit || current.is_empty() {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
    }
    if current.is_empty() {
        None
    } else {
        Some(current)
    }
}
pub(crate) fn ec_for_read_kmer_index(index: &KmerEcIndex, seq: &[u8]) -> Option<Vec<u32>> {
    if seq.len() < index.k {
        return None;
    }
    let mut current: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = Vec::new();
    let mut has_hit = false;
    for pos in 0..=seq.len() - index.k {
        if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + index.k])
            && let Some(ec) = lookup_ec(index, fwd).or_else(|| lookup_ec(index, rev))
        {
            if !has_hit {
                current.extend_from_slice(ec);
                has_hit = true;
            } else {
                next.clear();
                intersect_sorted(&current, ec, &mut next);
                std::mem::swap(&mut current, &mut next);
                if current.is_empty() {
                    break;
                }
            }
        }
    }
    if !has_hit || current.is_empty() {
        return None;
    }
    if let Some(onlist) = index.onlist.as_deref() {
        current.retain(|&t| (t as usize) < onlist.len() && onlist[t as usize]);
    }
    if current.is_empty() {
        None
    } else {
        Some(current)
    }
}
