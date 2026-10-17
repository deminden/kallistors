//! Diagnostic minimizer probes and candidate matching.

use super::super::core::{
    block_index_for_position, special_unitig_for_kmer, special_unitig_for_kmer_raw,
};
use super::super::debug::format_positions_sample;
use super::super::minimizers::{
    fill_revcomp, minhash_candidates_for_kmer, minimizer_for_kmer_strict, minimizers_for_kmer,
};
use super::super::{BifrostIndex, MinimizerDebugRow, Strand, encode_kmer_pair};

/// Debug minimizer lookups and unitig matches for specific read positions.
pub fn debug_minimizer_lookup(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    positions: &[usize],
    kallisto_bifrost_find: bool,
    kallisto_strict: bool,
    skip_overcrowded_minimizer: bool,
) -> Vec<MinimizerDebugRow> {
    let allow_forward_base = strand != Strand::Reverse;
    let allow_rev_base = strand != Strand::Forward;
    let diff = index.k.saturating_sub(index.g);
    let mut rev_buf: Vec<u8> = Vec::new();
    let mut out = Vec::new();
    let format_mphf_key = |key: &[u8; 8]| -> String {
        let mut out = String::with_capacity(16);
        for byte in key {
            use std::fmt::Write;
            let _ = write!(out, "{byte:02x}");
        }
        out
    };

    for &pos in positions {
        if pos + index.k > seq.len() {
            out.push(MinimizerDebugRow {
                read_pos: pos,
                kmer: "-".to_string(),
                kmer_canon: "-".to_string(),
                used_revcomp: false,
                min_pos: None,
                min_seq: None,
                mphf_key: None,
                mphf_level: None,
                mphf_hash_raw: None,
                mphf_hash_domain: None,
                mphf_bucket: None,
                mphf_rank: None,
                mphf_final_hash_hit: None,
                mphf_hit: false,
                positions: None,
                sample_positions: None,
                matched: false,
                match_unitig_id: None,
                match_unitig_pos: None,
                match_block_idx: None,
                is_km: None,
                special_code: None,
                special_in_dlist: None,
                special_in_hmap: None,
                special_uid: None,
                note: "out_of_range".to_string(),
            });
            continue;
        }

        let kmer = &seq[pos..pos + index.k];
        let kmer_text = String::from_utf8_lossy(kmer).into_owned();
        if !fill_revcomp(kmer, &mut rev_buf) {
            out.push(MinimizerDebugRow {
                read_pos: pos,
                kmer: kmer_text.clone(),
                kmer_canon: kmer_text,
                used_revcomp: false,
                min_pos: None,
                min_seq: None,
                mphf_key: None,
                mphf_level: None,
                mphf_hash_raw: None,
                mphf_hash_domain: None,
                mphf_bucket: None,
                mphf_rank: None,
                mphf_final_hash_hit: None,
                mphf_hit: false,
                positions: None,
                sample_positions: None,
                matched: false,
                match_unitig_id: None,
                match_unitig_pos: None,
                match_block_idx: None,
                is_km: None,
                special_code: None,
                special_in_dlist: None,
                special_in_hmap: None,
                special_uid: None,
                note: "invalid_kmer".to_string(),
            });
            continue;
        }

        let use_revcomp = rev_buf.as_slice() < kmer;

        let kmer_canon = if use_revcomp {
            rev_buf.as_slice()
        } else {
            kmer
        };
        let kmer_canon_text = String::from_utf8_lossy(kmer_canon).into_owned();
        let min_input = if kallisto_bifrost_find || kallisto_strict {
            kmer_canon
        } else {
            kmer
        };
        let mut special_code: Option<u64> = None;
        let mut special_in_dlist: Option<bool> = None;
        let mut special_in_hmap: Option<bool> = None;
        let mut special_uid: Option<usize> = None;
        if let Some((fwd, rev)) = encode_kmer_pair(kmer) {
            let code = fwd.min(rev);
            special_code = Some(code);
            if let Some(dlist) = index.dlist.as_ref() {
                special_in_dlist = Some(dlist.contains(&code));
            }
            if let Some(map) = index.h_kmer_map.as_ref() {
                if let Some(&uid) = map.get(&code) {
                    special_in_hmap = Some(true);
                    special_uid = Some(uid);
                } else {
                    special_in_hmap = Some(false);
                }
            }
        }
        let (mut min_candidates, mut fallback_min) = if kallisto_bifrost_find {
            match minhash_candidates_for_kmer(min_input, index.g) {
                Some(v) => v,
                None => {
                    out.push(MinimizerDebugRow {
                        read_pos: pos,
                        kmer: kmer_text.clone(),
                        kmer_canon: kmer_canon_text,
                        used_revcomp: use_revcomp,
                        min_pos: None,
                        min_seq: None,
                        mphf_key: None,
                        mphf_level: None,
                        mphf_hash_raw: None,
                        mphf_hash_domain: None,
                        mphf_bucket: None,
                        mphf_rank: None,
                        mphf_final_hash_hit: None,
                        mphf_hit: false,
                        positions: None,
                        sample_positions: None,
                        matched: false,
                        match_unitig_id: None,
                        match_unitig_pos: None,
                        match_block_idx: None,
                        is_km: None,
                        special_code,
                        special_in_dlist,
                        special_in_hmap,
                        special_uid,
                        note: "no_minimizer".to_string(),
                    });
                    continue;
                }
            }
        } else {
            let candidates = if kallisto_strict {
                match minimizer_for_kmer_strict(min_input, index.g) {
                    Some(v) => vec![v],
                    None => {
                        out.push(MinimizerDebugRow {
                            read_pos: pos,
                            kmer: kmer_text.clone(),
                            kmer_canon: kmer_canon_text,
                            used_revcomp: use_revcomp,
                            min_pos: None,
                            min_seq: None,
                            mphf_key: None,
                            mphf_level: None,
                            mphf_hash_raw: None,
                            mphf_hash_domain: None,
                            mphf_bucket: None,
                            mphf_rank: None,
                            mphf_final_hash_hit: None,
                            mphf_hit: false,
                            positions: None,
                            sample_positions: None,
                            matched: false,
                            match_unitig_id: None,
                            match_unitig_pos: None,
                            match_block_idx: None,
                            is_km: None,
                            special_code,
                            special_in_dlist,
                            special_in_hmap,
                            special_uid,
                            note: "no_minimizer".to_string(),
                        });
                        continue;
                    }
                }
            } else {
                match minimizers_for_kmer(min_input, index.g) {
                    Some(v) => v,
                    None => {
                        out.push(MinimizerDebugRow {
                            read_pos: pos,
                            kmer: kmer_text.clone(),
                            kmer_canon: kmer_canon_text,
                            used_revcomp: use_revcomp,
                            min_pos: None,
                            min_seq: None,
                            mphf_key: None,
                            mphf_level: None,
                            mphf_hash_raw: None,
                            mphf_hash_domain: None,
                            mphf_bucket: None,
                            mphf_rank: None,
                            mphf_final_hash_hit: None,
                            mphf_hit: false,
                            positions: None,
                            sample_positions: None,
                            matched: false,
                            match_unitig_id: None,
                            match_unitig_pos: None,
                            match_block_idx: None,
                            is_km: None,
                            special_code,
                            special_in_dlist,
                            special_in_hmap,
                            special_uid,
                            note: "no_minimizer".to_string(),
                        });
                        continue;
                    }
                }
            };
            (candidates, None)
        };

        let mut fallback_active = false;
        'min_loop: loop {
            let mut skip_minimizer = false;
            for (min_bytes, min_pos) in min_candidates.iter().copied() {
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
                let min_seq = kmer_canon
                    .get(min_pos..min_pos + index.g)
                    .map(|v| String::from_utf8_lossy(v).into_owned());
                let (min_idx, mphf_debug) = index.mphf.debug_lookup(&min_bytes);
                let mphf_key = Some(format_mphf_key(&min_bytes));
                let Some(min_idx) = min_idx else {
                    let mut note = "mphf_miss".to_string();
                    if fallback_active {
                        note.push_str("|minhash_fallback");
                    }
                    out.push(MinimizerDebugRow {
                        read_pos: pos,
                        kmer: kmer_text.clone(),
                        kmer_canon: kmer_canon_text.clone(),
                        used_revcomp: use_revcomp,
                        min_pos: Some(min_pos),
                        min_seq,
                        mphf_key,
                        mphf_level: Some(mphf_debug.level),
                        mphf_hash_raw: Some(mphf_debug.hash_raw),
                        mphf_hash_domain: Some(mphf_debug.hash_domain),
                        mphf_bucket: mphf_debug.bucket,
                        mphf_rank: mphf_debug.rank,
                        mphf_final_hash_hit: Some(mphf_debug.final_hash_hit),
                        mphf_hit: false,
                        positions: None,
                        sample_positions: None,
                        matched: false,
                        match_unitig_id: None,
                        match_unitig_pos: None,
                        match_block_idx: None,
                        is_km: None,
                        special_code,
                        special_in_dlist,
                        special_in_hmap,
                        special_uid,
                        note,
                    });
                    continue;
                };

                let positions = &index.minz_positions[min_idx as usize];
                if positions.is_empty() {
                    let mut note = "no_positions".to_string();
                    if fallback_active {
                        note.push_str("|minhash_fallback");
                    }
                    out.push(MinimizerDebugRow {
                        read_pos: pos,
                        kmer: kmer_text.clone(),
                        kmer_canon: kmer_canon_text.clone(),
                        used_revcomp: use_revcomp,
                        min_pos: Some(min_pos),
                        min_seq,
                        mphf_key,
                        mphf_level: Some(mphf_debug.level),
                        mphf_hash_raw: Some(mphf_debug.hash_raw),
                        mphf_hash_domain: Some(mphf_debug.hash_domain),
                        mphf_bucket: mphf_debug.bucket,
                        mphf_rank: mphf_debug.rank,
                        mphf_final_hash_hit: Some(mphf_debug.final_hash_hit),
                        mphf_hit: true,
                        positions: Some(0),
                        sample_positions: Some("-".to_string()),
                        matched: false,
                        match_unitig_id: None,
                        match_unitig_pos: None,
                        match_block_idx: None,
                        is_km: None,
                        special_code,
                        special_in_dlist,
                        special_in_hmap,
                        special_uid,
                        note,
                    });
                    continue;
                }

                let mut matched = None;
                let mut special_note: Option<&str> = None;
                for &pos_id in positions {
                    let unitig_id_raw = (pos_id >> 32) as u32;
                    if unitig_id_raw == u32::MAX {
                        if kallisto_strict {
                            special_note = Some("special_strict_skip");
                            continue;
                        }
                        let overcrowded = (pos_id & 0x8000_0000) != 0;
                        if overcrowded {
                            special_note = Some("special_overcrowded");
                            if skip_overcrowded_minimizer {
                                if kallisto_bifrost_find {
                                    skip_minimizer = true;
                                    break;
                                }
                                continue;
                            }
                        }
                        let special_uid = if kallisto_bifrost_find {
                            special_unitig_for_kmer_raw(index, kmer)
                        } else {
                            special_unitig_for_kmer(index, kmer)
                        };
                        if let Some(uid) = special_uid {
                            let used_revcomp = kmer_canon != kmer;
                            if (used_revcomp && !allow_rev_base)
                                || (!used_revcomp && !allow_forward_base)
                            {
                                special_note = Some("special_abundant_strand_skip");
                                continue;
                            }
                            let block_idx =
                                block_index_for_position(index.ec_blocks[uid].as_slice(), 0);
                            matched = Some((uid, 0usize, block_idx, None, used_revcomp));
                            special_note = Some("special_abundant_match");
                            break;
                        }
                        special_note = Some("special_abundant_miss");
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
                        if !kallisto_bifrost_find && min_pos != km_pos && min_pos + km_pos != diff {
                            continue;
                        }
                        let used_revcomp = kmer_canon != kmer;
                        if (used_revcomp && !allow_rev_base)
                            || (!used_revcomp && !allow_forward_base)
                        {
                            continue;
                        }
                        let unitig_id = index.unitigs.len() + uid;
                        let block_idx =
                            block_index_for_position(index.ec_blocks[unitig_id].as_slice(), 0);
                        matched = Some((unitig_id, 0usize, block_idx, Some(true), used_revcomp));
                        break;
                    } else {
                        let uid = unitig_id as usize;
                        if uid >= index.unitigs.len() {
                            continue;
                        }
                        let unitig = &index.unitigs[uid];
                        if unitig.len() < index.k {
                            continue;
                        }
                        let rel_pos = rel_pos as usize;
                        let mut found = false;

                        let start_fwd = rel_pos as isize - min_pos_fwd as isize;
                        if allow_forward_base && start_fwd >= 0 {
                            let start = start_fwd as usize;
                            if start + index.k <= unitig.len()
                                && &unitig[start..start + index.k] == kmer
                            {
                                let block_idx = block_index_for_position(
                                    index.ec_blocks[uid].as_slice(),
                                    start,
                                );
                                matched = Some((uid, start, block_idx, Some(false), false));
                                found = true;
                            }
                        }
                        if !found {
                            let start_rev = rel_pos as isize - diff as isize + min_pos_rev as isize;
                            if allow_rev_base && start_rev >= 0 {
                                let start = start_rev as usize;
                                if start + index.k <= unitig.len()
                                    && &unitig[start..start + index.k] == rev_buf.as_slice()
                                {
                                    let block_idx = block_index_for_position(
                                        index.ec_blocks[uid].as_slice(),
                                        start,
                                    );
                                    matched = Some((uid, start, block_idx, Some(false), true));
                                    found = true;
                                }
                            }
                        }
                        if !found {
                            let start_low = rel_pos.saturating_sub(diff);
                            let start_high = rel_pos.min(unitig.len().saturating_sub(index.k));
                            for start in start_low..=start_high {
                                if allow_forward_base && &unitig[start..start + index.k] == kmer {
                                    let block_idx = block_index_for_position(
                                        index.ec_blocks[uid].as_slice(),
                                        start,
                                    );
                                    matched = Some((uid, start, block_idx, Some(false), false));
                                    special_note = Some("match_relaxed");
                                    found = true;
                                    break;
                                }
                                if allow_rev_base
                                    && &unitig[start..start + index.k] == rev_buf.as_slice()
                                {
                                    let block_idx = block_index_for_position(
                                        index.ec_blocks[uid].as_slice(),
                                        start,
                                    );
                                    matched = Some((uid, start, block_idx, Some(false), true));
                                    special_note = Some("match_relaxed");
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if found {
                            break;
                        }
                    }
                }

                let (match_unitig_id, match_unitig_pos, match_block_idx, is_km, note, match_rev) =
                    if let Some((uid, start, block_idx, is_km, match_rev)) = matched {
                        (
                            Some(uid),
                            Some(start),
                            block_idx,
                            is_km,
                            "match",
                            Some(match_rev),
                        )
                    } else if skip_minimizer {
                        (None, None, None, None, "special_overcrowded_skip", None)
                    } else if let Some(note) = special_note {
                        (None, None, None, None, note, None)
                    } else {
                        (None, None, None, None, "no_match", None)
                    };
                let used_revcomp = match_rev.unwrap_or(use_revcomp);
                let mut note_text = note.to_string();
                if fallback_active {
                    note_text.push_str("|minhash_fallback");
                }

                out.push(MinimizerDebugRow {
                    read_pos: pos,
                    kmer: kmer_text.clone(),
                    kmer_canon: kmer_canon_text.clone(),
                    used_revcomp,
                    min_pos: Some(min_pos),
                    min_seq,
                    mphf_key,
                    mphf_level: Some(mphf_debug.level),
                    mphf_hash_raw: Some(mphf_debug.hash_raw),
                    mphf_hash_domain: Some(mphf_debug.hash_domain),
                    mphf_bucket: mphf_debug.bucket,
                    mphf_rank: mphf_debug.rank,
                    mphf_final_hash_hit: Some(mphf_debug.final_hash_hit),
                    mphf_hit: true,
                    positions: Some(positions.len()),
                    sample_positions: Some(format_positions_sample(positions)),
                    matched: match_unitig_id.is_some(),
                    match_unitig_id,
                    match_unitig_pos,
                    match_block_idx,
                    is_km,
                    special_code,
                    special_in_dlist,
                    special_in_hmap,
                    special_uid,
                    note: note_text,
                });
                if skip_minimizer {
                    break;
                }
            }
            if !kallisto_bifrost_find {
                break;
            }
            if skip_minimizer && let Some((bytes, pos)) = fallback_min.take() {
                min_candidates = vec![(bytes, pos)];
                fallback_active = true;
                continue 'min_loop;
            }
            break;
        }
    }
    out
}
