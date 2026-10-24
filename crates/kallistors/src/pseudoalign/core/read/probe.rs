//! Ordered general minimizer probes and their diagnostic observations.

use crate::pseudoalign::core::matcher::{
    match_kmer_at_pos, match_kmer_direct, match_unitig_candidate, special_unitig_for_kmer,
    special_unitig_for_kmer_raw,
};
use crate::pseudoalign::core::runtime::{ALLOW_TAIL_MINIMIZER_IN_BACKOFF, env_flag, minimizer_hex};
use crate::pseudoalign::debug::{MinimizerCandidateTrace, ReadDebugState, format_positions_sample};
use crate::pseudoalign::minimizers::{
    minhash_next_after_hash, minhash_primary_for_kmer, minimizer_for_kmer_strict,
    minimizer_tail_for_kmer, minimizers_for_kmer, minimizers_ranked_for_kmer,
};
use crate::pseudoalign::{BifrostIndex, PseudoalignOptions, Strand};

pub(super) type CandidateMatch = (usize, usize, usize, usize, bool, bool, bool, bool, bool);

pub(super) struct ProbeCandidates {
    pub(super) matched: Option<CandidateMatch>,
    pub(super) min_candidates: Vec<([u8; 8], usize)>,
    pub(super) any_mphf: bool,
    pub(super) any_positions: bool,
}

pub(super) enum ProbeOutcome {
    Skip,
    Candidates(ProbeCandidates),
}

#[allow(clippy::too_many_arguments)]
pub(super) fn probe_kmer(
    index: &BifrostIndex,
    kmer: &[u8],
    rev_buf: &[u8],
    pos: usize,
    strand: Strand,
    options: PseudoalignOptions,
    in_backoff: bool,
    mut dbg: Option<&mut ReadDebugState>,
) -> Option<ProbeOutcome> {
    let allow_forward_base = strand != Strand::Reverse;
    let allow_rev_base = strand != Strand::Forward;
    let diff = index.k.saturating_sub(index.g);
    let mut matched = None;
    let mut any_mphf = false;
    let mut any_positions = false;
    let use_revcomp = rev_buf < kmer;
    let kmer_canon = if use_revcomp { rev_buf } else { kmer };
    let min_input = if options.kallisto_bifrost_find || options.kallisto_strict {
        kmer_canon
    } else {
        kmer
    };

    if options.kallisto_direct_kmer {
        if let Some((uid, start, used_revcomp)) =
            match_kmer_direct(index, kmer, allow_forward_base, allow_rev_base)
        {
            matched = Some((
                uid,
                start,
                pos,
                0,
                used_revcomp,
                false,
                false,
                false,
                !used_revcomp,
            ));
        }
        if matched.is_none() {
            return Some(ProbeOutcome::Skip);
        }
    }

    let (mut min_candidates, mut min_hash_current) = if options.kallisto_bifrost_find {
        match minhash_primary_for_kmer(min_input, index.g) {
            Some((min_hash, candidates)) => (candidates, Some(min_hash)),
            None => {
                return Some(ProbeOutcome::Skip);
            }
        }
    } else {
        let candidates = if options.kallisto_strict {
            minimizer_for_kmer_strict(min_input, index.g).map(|v| vec![v])
        } else if options.kallisto_enum {
            minimizers_ranked_for_kmer(min_input, index.g, 2)
        } else {
            minimizers_for_kmer(min_input, index.g)
        };
        let Some(candidates) = candidates else {
            return Some(ProbeOutcome::Skip);
        };
        (candidates, None)
    };
    if let Some(state) = dbg.as_deref_mut() {
        state.saw_valid_kmer = true;
    }
    let mut tried_tail_minimizer = false;
    let mut saw_special_overcrowded_unresolved = false;
    'min_outer: loop {
        let mut request_next_min = false;
        for (min_bytes, min_pos) in min_candidates.iter().copied() {
            if matched.is_some() {
                break;
            }
            let mut candidate = MinimizerCandidateTrace {
                read_pos: pos,
                min_pos,
                minimizer: minimizer_hex(min_bytes),
                mphf_hit: false,
                positions_len: 0,
                sample_positions: "-".to_string(),
                has_special: false,
                overcrowded: false,
                matched: false,
            };
            let min_pos_fwd = if use_revcomp {
                diff.saturating_sub(min_pos)
            } else {
                min_pos
            };
            let min_pos_rev = if use_revcomp {
                min_pos
            } else {
                diff.saturating_sub(min_pos)
            };
            let Some(min_idx) = index.mphf.lookup(&min_bytes) else {
                if let Some(state) = dbg.as_deref_mut()
                    && state.first_mphf_miss.is_none()
                {
                    state.first_mphf_miss = Some((pos, min_pos));
                }
                if let Some(state) = dbg.as_deref_mut()
                    && state.minimizer_candidates.len() < 1024
                {
                    state.minimizer_candidates.push(candidate);
                }
                continue;
            };
            candidate.mphf_hit = true;
            any_mphf = true;
            let positions = &index.minz_positions[min_idx as usize];
            if positions.is_empty() {
                if let Some(state) = dbg.as_deref_mut()
                    && state.first_no_positions.is_none()
                {
                    state.first_no_positions = Some((pos, min_pos));
                }
                if let Some(state) = dbg.as_deref_mut()
                    && state.minimizer_candidates.len() < 1024
                {
                    state.minimizer_candidates.push(candidate);
                }
                continue;
            }
            let mut minimizer_overcrowded = false;
            if options.skip_overcrowded_minimizer {
                for &pos_id in positions {
                    let unitig_id_raw = (pos_id >> 32) as u32;
                    if unitig_id_raw == u32::MAX && (pos_id & 0x8000_0000) != 0 {
                        minimizer_overcrowded = true;
                        break;
                    }
                }
            }
            candidate.positions_len = positions.len();
            candidate.sample_positions = format_positions_sample(positions);
            candidate.overcrowded = minimizer_overcrowded;
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
                let unitig_id_raw = (pos_id >> 32) as u32;
                if unitig_id_raw == u32::MAX {
                    candidate.has_special = true;
                    if options.kallisto_strict {
                        continue;
                    }
                    let overcrowded = (pos_id & 0x8000_0000) != 0;
                    candidate.overcrowded = overcrowded;

                    let special_uid = if options.kallisto_bifrost_find {
                        special_unitig_for_kmer_raw(index, kmer)
                    } else {
                        special_unitig_for_kmer(index, kmer)
                    };
                    if let Some(uid) = special_uid {
                        let used_revcomp = kmer_canon != kmer;
                        if (used_revcomp && !allow_rev_base)
                            || (!used_revcomp && !allow_forward_base)
                        {
                            continue;
                        }
                        candidate.matched = true;
                        local_match = Some((
                            uid,
                            0usize,
                            pos,
                            min_pos,
                            used_revcomp,
                            true,
                            minimizer_overcrowded,
                            false,
                            !used_revcomp,
                        ));
                        break;
                    }
                    if overcrowded {
                        saw_special_overcrowded_unresolved = true;
                    }
                    if options.kallisto_bifrost_find && overcrowded {
                        request_next_min = true;
                    } else if options.skip_overcrowded_minimizer && overcrowded {
                        if let Some(state) = dbg.as_deref_mut() {
                            state.overcrowded_jump_skip = true;
                        }
                        continue;
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
                    let km_seq = index.km_unitig_sequence(uid)?;
                    if km_seq.as_ref() != kmer_canon {
                        continue;
                    }
                    if min_pos != km_pos && min_pos + km_pos != diff {
                        continue;
                    }
                    let used_revcomp = kmer_canon != kmer;
                    if (used_revcomp && !allow_rev_base) || (!used_revcomp && !allow_forward_base) {
                        continue;
                    }
                    candidate.matched = true;
                    local_match = Some((
                        index.unitig_count() + uid,
                        0usize,
                        pos,
                        min_pos,
                        used_revcomp,
                        false,
                        minimizer_overcrowded,
                        false,
                        !used_revcomp,
                    ));
                    if used_revcomp && let Some(state) = dbg.as_deref_mut() {
                        state.used_revcomp = true;
                    }
                    break;
                } else {
                    let uid = unitig_id as usize;
                    if uid >= index.unitig_count() {
                        continue;
                    }
                    if let Some((start, used_revcomp, matched_relaxed, forward_strand)) =
                        match_unitig_candidate(
                            index.unitig_sequence(uid)?.as_ref(),
                            index.k,
                            rel_pos as usize,
                            diff,
                            min_pos_fwd,
                            min_pos_rev,
                            kmer,
                            rev_buf,
                            allow_forward_base,
                            allow_rev_base,
                            true,
                        )
                    {
                        candidate.matched = true;
                        local_match = Some((
                            uid,
                            start,
                            pos,
                            min_pos,
                            used_revcomp,
                            false,
                            minimizer_overcrowded,
                            matched_relaxed,
                            forward_strand,
                        ));
                        if used_revcomp && let Some(state) = dbg.as_deref_mut() {
                            state.used_revcomp = true;
                        }
                        break;
                    }
                }
            }
            if let Some(hit) = local_match {
                matched = Some(hit);
                if let Some(state) = dbg.as_deref_mut()
                    && state.minimizer_candidates.len() < 1024
                {
                    state.minimizer_candidates.push(candidate);
                }
                break;
            }
            if let Some(state) = dbg.as_deref_mut()
                && state.first_no_match.is_none()
            {
                state.first_no_match = Some((pos, min_pos));
            }
            if let Some(state) = dbg.as_deref_mut()
                && state.minimizer_candidates.len() < 1024
            {
                state.minimizer_candidates.push(candidate);
            }
        }
        if matched.is_some() {
            break;
        }
        if !options.kallisto_bifrost_find
            && !options.kallisto_strict
            && saw_special_overcrowded_unresolved
        {
            let mut probe_rev_buf: Vec<u8> = Vec::new();
            if let Some((uid, start, used_revcomp, _block_idx, matched_relaxed)) = match_kmer_at_pos(
                index,
                kmer,
                allow_forward_base,
                allow_rev_base,
                diff,
                &mut probe_rev_buf,
                options.kallisto_direct_kmer,
                options.kallisto_enum,
                options.kallisto_strict,
                options.skip_overcrowded_minimizer,
                true,
                false,
            ) {
                matched = Some((
                    uid,
                    start,
                    pos,
                    0usize,
                    used_revcomp,
                    false,
                    false,
                    matched_relaxed,
                    !used_revcomp,
                ));
                break;
            }
        }
        if !options.kallisto_bifrost_find
            && !options.kallisto_strict
            && !options.kallisto_enum
            && !tried_tail_minimizer
            && (!in_backoff
                || env_flag(
                    &ALLOW_TAIL_MINIMIZER_IN_BACKOFF,
                    "KALLISTORS_ALLOW_TAIL_MINIMIZER_IN_BACKOFF",
                ))
            && let Some((tail_bytes, tail_pos)) = minimizer_tail_for_kmer(min_input, index.g)
            && !min_candidates
                .iter()
                .any(|&(min_bytes, min_pos)| min_bytes == tail_bytes && min_pos == tail_pos)
        {
            tried_tail_minimizer = true;
            min_candidates = vec![(tail_bytes, tail_pos)];
            continue 'min_outer;
        }
        if !options.kallisto_bifrost_find {
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
    if let Some(state) = dbg {
        if any_mphf {
            state.saw_mphf_hit = true;
        }
        if any_positions {
            state.saw_positions = true;
        }
    }

    Some(ProbeOutcome::Candidates(ProbeCandidates {
        matched,
        min_candidates,
        any_mphf,
        any_positions,
    }))
}
