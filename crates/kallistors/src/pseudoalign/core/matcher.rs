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

#[inline(always)]
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
    match_unitig_candidate_encoded_fixed::<0>(
        index,
        unitig_id,
        rel_pos,
        diff,
        min_pos_fwd,
        min_pos_rev,
        read_fwd,
        read_rev,
        allow_forward,
        allow_rev,
        allow_relaxed,
    )
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn match_unitig_candidate_encoded_typed<const BUS: bool>(
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
    if BUS && index.k == 31 {
        match_unitig_candidate_encoded_fixed::<31>(
            index,
            unitig_id,
            rel_pos,
            diff,
            min_pos_fwd,
            min_pos_rev,
            read_fwd,
            read_rev,
            allow_forward,
            allow_rev,
            allow_relaxed,
        )
    } else {
        match_unitig_candidate_encoded_fixed::<0>(
            index,
            unitig_id,
            rel_pos,
            diff,
            min_pos_fwd,
            min_pos_rev,
            read_fwd,
            read_rev,
            allow_forward,
            allow_rev,
            allow_relaxed,
        )
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn match_unitig_candidate_encoded_fixed<const K: usize>(
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
    let k = if K == 0 { index.k } else { K };
    let sequence = index.encoded_unitigs.unitig_bases.sequence(unitig_id)?;
    let unitig_len = sequence.len;
    if unitig_len < k {
        return None;
    }

    // A length-k node has exactly one possible reference word. A mismatch in
    // both allowed orientations proves that neither strict nor relaxed search
    // can succeed. Keep hit handling below so its priority and flags are exact.
    if unitig_len == k {
        let code = sequence.kmer_at(0, k)?;
        if (!allow_forward || code != read_fwd) && (!allow_rev || code != read_rev) {
            return None;
        }
    }

    let rel = rel_pos as isize;
    let start_fwd = rel - min_pos_fwd as isize;
    let start_rev = rel - diff as isize + min_pos_rev as isize;
    if start_fwd == start_rev {
        // Complementary minimizer offsets address one reference word for both
        // strand checks. Preserve forward priority while decoding it only once.
        if (allow_forward || allow_rev) && start_fwd >= 0 {
            let start = start_fwd as usize;
            if start + k <= unitig_len {
                let code = sequence.kmer_at(start, k);
                if allow_forward && code == Some(read_fwd) {
                    return Some((start, false, false, true));
                }
                if allow_rev && code == Some(read_rev) {
                    return Some((start, true, false, false));
                }
            }
        }
    } else {
        if allow_forward && start_fwd >= 0 {
            let start = start_fwd as usize;
            if start + k <= unitig_len && sequence.kmer_at(start, k) == Some(read_fwd) {
                return Some((start, false, false, true));
            }
        }

        if allow_rev && start_rev >= 0 {
            let start = start_rev as usize;
            if start + k <= unitig_len && sequence.kmer_at(start, k) == Some(read_rev) {
                return Some((start, true, false, false));
            }
        }
    }
    if !allow_relaxed {
        return None;
    }
    let max_start = unitig_len - k;
    let start_lo = rel_pos.saturating_sub(diff).min(max_start);
    let start_hi = rel_pos.min(max_start);
    if start_lo > start_hi {
        return None;
    }
    if let Some(windows) = sequence.range(start_lo, start_hi, k) {
        return windows
            .first_match(read_fwd, read_rev, allow_forward, allow_rev)
            .map(|(start, rev)| (start, rev, true, !rev));
    }
    for start in start_lo..=start_hi {
        let Some(code) = sequence.kmer_at(start, k) else {
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
    let encoded_match_codes = (index.k <= 32 && kmer.iter().all(u8::is_ascii_uppercase))
        .then(|| encode_kmer_pair(kmer))
        .flatten();
    let min_codes = encoded_match_codes.map(|(f, r)| {
        if use_revcomp && (kallisto_bifrost_find || kallisto_strict) {
            (r, f)
        } else {
            (f, r)
        }
    });
    let (mut min_candidates, mut min_hash_current) = if kallisto_bifrost_find {
        let (min_hash, candidates) = if let Some((f, r)) = min_codes {
            super::super::minimizers::minhash_primary_for_encoded(min_input.len(), index.g, f, r)?
        } else {
            minhash_primary_for_kmer(min_input, index.g)?
        };
        (candidates, Some(min_hash))
    } else {
        let candidates = if kallisto_strict {
            vec![minimizer_for_kmer_strict(min_input, index.g)?]
        } else if kallisto_enum {
            minimizers_ranked_for_kmer(min_input, index.g, 2)?
        } else if let Some((f, r)) = min_codes {
            let mut candidates = Vec::new();
            if !super::super::minimizers::minimizers_for_encoded_kmer_into(
                min_input,
                index.g,
                f,
                r,
                &mut candidates,
            ) {
                return None;
            }
            candidates
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
    let fingerprint_codes = (index.k == 31 && index.g == 23 && !index.minz_fingerprints.is_empty())
        .then(|| encoded_match_codes.or_else(|| encode_kmer_pair(kmer)))
        .flatten();
    let mut query_fingerprints = super::super::fingerprints::QueryFingerprintCache::default();

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
            if let Some((forward, reverse)) = fingerprint_codes {
                let fingerprint = index
                    .minz_fingerprints
                    .get(min_idx as usize)
                    .copied()
                    .unwrap_or(0);
                if fingerprint != 0
                    && fingerprint != super::super::fingerprints::DISABLED
                    && fingerprint
                        != super::super::fingerprints::word_fingerprint(u64::from_le_bytes(
                            min_bytes,
                        ))
                    && !query_fingerprints.contains(forward, reverse, fingerprint)
                {
                    continue;
                }
            }
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
                    if uid >= index.km_unitig_count() {
                        continue;
                    }
                    let km_pos = rel_pos as usize;
                    if let Some((forward, reverse)) = encoded_match_codes {
                        // Packed graphs already retain this complete k-mer. Do
                        // not rebuild its ASCII sequence for each rescue probe.
                        if index
                            .encoded_unitigs
                            .km_bases
                            .encode_kmer_at(uid, 0, index.k)
                            != Some(forward.min(reverse))
                        {
                            continue;
                        }
                    } else {
                        let km_seq = index.km_unitig_sequence(uid)?;
                        if km_seq.as_ref() != kmer_canon {
                            continue;
                        }
                    }
                    if min_pos != km_pos && min_pos + km_pos != _diff {
                        continue;
                    }
                    let used_revcomp = kmer_canon != kmer;
                    if (used_revcomp && !allow_rev) || (!used_revcomp && !allow_forward) {
                        continue;
                    }
                    let unitig_id = index.unitig_count() + uid;
                    let block_idx = block_index_for_position_fast(index, unitig_id, 0)?;
                    return Some((unitig_id, 0usize, used_revcomp, block_idx, false));
                } else {
                    let uid = unitig_id as usize;
                    if uid >= index.unitig_count() {
                        continue;
                    }
                    let matched = if let Some((forward, reverse)) = encoded_match_codes {
                        match_unitig_candidate_encoded(
                            index,
                            uid,
                            rel_pos as usize,
                            _diff,
                            min_pos_fwd,
                            min_pos_rev,
                            forward,
                            reverse,
                            allow_forward,
                            allow_rev,
                            allow_relaxed,
                        )
                    } else {
                        match_unitig_candidate(
                            index.unitig_sequence(uid)?.as_ref(),
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
                    };
                    if let Some((start, used_revcomp, matched_relaxed, _forward_strand)) = matched {
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
            min_candidates.clear();
            min_candidates.push((tail_bytes, tail_pos));
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
            min_candidates.clear();
            min_candidates.push(next_min);
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
        // Both orientations are indexed under the query's encoded bases;
        // the entry flag records orientation relative to the stored unitig.
        if let Some(entries) = map.get(&fwd) {
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
    // The flat table has identical bounds and is also available in compact quant indexes.
    let (lb, ub) = index.flat_ec.block_bounds(unitig_id, block_idx)?;
    if ub <= lb {
        return None;
    }
    let contig_start = lb as isize;
    let contig_len = (ub - lb) as isize;
    let start = start as isize;
    let dist = if forward_strand {
        contig_len - 1 - (start - contig_start)
    } else {
        start - contig_start
    };
    if dist >= 2 { Some(dist as usize) } else { None }
}

#[cfg(test)]
mod encoded_probe_tests {
    use super::*;

    #[test]
    fn shared_and_distinct_strict_words_match_byte_oracle() {
        let temp = tempfile::TempDir::new().unwrap();
        let fasta = temp.path().join("transcripts.fa");
        let path = temp.path().join("index.idx");
        let mut state = 17u64;
        let sequence: String = (0..500)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[((state >> 32) & 3) as usize] as char
            })
            .collect();
        std::fs::write(&fasta, format!(">tx\n{sequence}\n>short\nCGGGTGTGGTGGCGGGCACCTGTAGTCCCAG\n>two\nCGAGGTTCGGCGGCTAGATCTAGTTGTTACGT\n")).unwrap();
        crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default())
            .unwrap();
        let index = crate::pseudoalign::build_bifrost_index(&path).unwrap();
        assert!(index.unitigs.iter().any(|unitig| unitig.len() == 31));
        assert!(index.unitigs.iter().any(|unitig| unitig.len() == 32));
        for (uid, unitig) in index.unitigs.iter().enumerate() {
            if unitig.len() < index.k {
                continue;
            }
            for query in unitig
                .windows(index.k)
                .step_by(if unitig.len() <= 33 { 1 } else { 41 })
                .flat_map(|word| [word.to_vec(), crate::util::reverse_complement(word)])
            {
                for alter in [None, Some(13), Some(30)] {
                    let mut read = query.to_vec();
                    if let Some(position) = alter {
                        read[position] = if read[position] == b'A' { b'C' } else { b'A' };
                    }
                    let (fwd, rev) = encode_kmer_pair(&read).unwrap();
                    let mut reverse = Vec::new();
                    assert!(fill_revcomp(&read, &mut reverse));
                    for rel in [0, 3, 8, unitig.len() - index.k, unitig.len() - 1] {
                        for (min_fwd, min_rev) in [(0, 8), (2, 6), (8, 0), (1, 1)] {
                            for flags in 0..8 {
                                let forward = flags & 1 != 0;
                                let reverse_allowed = flags & 2 != 0;
                                let relaxed = flags & 4 != 0;
                                let expected = match_unitig_candidate(
                                    unitig,
                                    index.k,
                                    rel,
                                    8,
                                    min_fwd,
                                    min_rev,
                                    &read,
                                    &reverse,
                                    forward,
                                    reverse_allowed,
                                    relaxed,
                                );
                                let actual = match_unitig_candidate_encoded(
                                    &index,
                                    uid,
                                    rel,
                                    8,
                                    min_fwd,
                                    min_rev,
                                    fwd,
                                    rev,
                                    forward,
                                    reverse_allowed,
                                    relaxed,
                                );
                                let bus_actual = match_unitig_candidate_encoded_typed::<true>(
                                    &index,
                                    uid,
                                    rel,
                                    8,
                                    min_fwd,
                                    min_rev,
                                    fwd,
                                    rev,
                                    forward,
                                    reverse_allowed,
                                    relaxed,
                                );
                                assert_eq!(
                                    bus_actual, expected,
                                    "BUS uid={uid} rel={rel} min={min_fwd}/{min_rev} flags={flags}"
                                );
                                assert_eq!(
                                    actual, expected,
                                    "uid={uid} rel={rel} min={min_fwd}/{min_rev} flags={flags}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
