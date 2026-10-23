//! Minimizer enumeration and grouped hash generation.

use super::*;
pub(super) fn build_unitig_minimizer_index(
    unitigs: &[Vec<u8>],
    k: usize,
    g: usize,
    threads: usize,
) -> Result<MinimizerIndexParts> {
    let mut keys = Vec::new();
    let total_positions = unitigs.iter().try_fold(0usize, |acc, seq| {
        acc.checked_add(seq.len().saturating_sub(g) + 1)
            .ok_or_else(|| Error::InvalidFormat("too many minimizer positions".into()))
    })?;
    let bitmap_count = (total_positions as u64 >> 32) as usize + 1;
    let mut bitmap_blocks = vec![Vec::new(); bitmap_count.max(1)];
    if threads > 1 && total_positions >= PARALLEL_KMER_COLLECTION_MIN {
        let mut bases = Vec::with_capacity(unitigs.len());
        let mut global = 0u64;
        for seq in unitigs {
            bases.push(global);
            global += minimizer_position_count(seq.len(), g) as u64;
        }
        let pool = IndexWorkers::new(threads)?;
        // Bounded groups amortize allocation without retaining one output per unitig.
        for first in (0..unitigs.len()).step_by(4096 * threads) {
            let end = (first + 4096 * threads).min(unitigs.len());
            let shards = pool.install(|| {
                unitigs[first..end]
                    .par_chunks(512)
                    .enumerate()
                    .map(|(chunk, sequences)| {
                        let mut local_keys = Vec::new();
                        let mut local_blocks = vec![Vec::new(); bitmap_count.max(1)];
                        let mut scratch = MinimizerScratch::default();
                        for (relative, seq) in sequences.iter().enumerate() {
                            insert_selected_minimizers(
                                seq,
                                k,
                                g,
                                bases[first + chunk * 512 + relative],
                                &mut local_keys,
                                &mut local_blocks,
                                &mut scratch,
                            )?;
                        }
                        Ok::<_, Error>((local_keys, local_blocks))
                    })
                    .collect::<Result<Vec<_>>>()
            })?;
            for (local_keys, local_blocks) in shards {
                keys.extend(local_keys);
                for (dest, source) in bitmap_blocks.iter_mut().zip(local_blocks) {
                    dest.extend(source);
                }
            }
        }
    } else {
        let mut global = 0u64;
        let mut scratch = MinimizerScratch::default();
        for seq in unitigs {
            insert_selected_minimizers(
                seq,
                k,
                g,
                global,
                &mut keys,
                &mut bitmap_blocks,
                &mut scratch,
            )?;
            global += minimizer_position_count(seq.len(), g) as u64;
        }
    }

    // Keep exact numeric keys throughout collection and sorting, avoiding a
    // second full-size allocation solely to change element alignment.
    if threads > 1 && keys.len() >= PARALLEL_KMER_COLLECTION_MIN {
        let pool = IndexWorkers::new(threads)?;
        pool.install(|| keys.radix_sort_unstable());
    } else {
        keys.sort_unstable();
    }
    let written = keys.len();
    keys.dedup();
    super::super::scratch_pages::discard_unused(&mut keys, written);
    Ok((keys, bitmap_blocks))
}

#[derive(Default)]
pub(super) struct MinimizerScratch {
    pub(super) window: [(usize, u64); 32],
}

pub(super) fn insert_selected_minimizers(
    seq: &[u8],
    k: usize,
    g: usize,
    global_base: u64,
    keys: &mut Vec<u64>,
    bitmap_blocks: &mut [Vec<u32>],
    scratch: &mut MinimizerScratch,
) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    let Some(bounds) = bifrost_minimizer_bounds(k, g) else {
        return Ok(());
    };
    if seq
        .chunks(32)
        .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none())
    {
        return Err(Error::InvalidFormat("invalid minimizer".into()));
    }
    let window = &mut scratch.window;
    let mut head = 0usize;
    let mut length = 0usize;
    debug_assert!(bounds.last - bounds.start < 32);
    let (mut forward, mut reverse) =
        minimizer_rep_hash_words(&seq[bounds.start..bounds.start + g])?;
    let remove_forward = REP_HASH_VALS.map(|value| value.rotate_left(g as u32));
    let remove_reverse = REP_HASH_VALS.map(|value| value.rotate_right(1));
    let add_reverse = REP_HASH_VALS.map(|value| value.rotate_left((g - 1) as u32));
    let mut next_pos = bounds.start;
    let mut started = false;
    let mut last_selected = None;
    let mut previous_minimum = None;
    for kmer_start in 0..=seq.len() - k {
        let window_start = kmer_start + bounds.start;
        let window_end = kmer_start + bounds.last;
        while length != 0 && window[head].0 < window_start {
            head = (head + 1) & 31;
            length -= 1;
        }
        while next_pos <= window_end {
            if started {
                let old = seq[next_pos - 1] as usize;
                let new = seq[next_pos + g - 1] as usize;
                forward = forward.rotate_left(1)
                    ^ remove_forward[HASH_CODES[old] as usize]
                    ^ REP_HASH_VALS[HASH_CODES[new] as usize];
                reverse = reverse.rotate_right(1)
                    ^ remove_reverse[HASH_COMP_CODES[old] as usize]
                    ^ add_reverse[HASH_COMP_CODES[new] as usize];
            } else {
                started = true;
            }
            let hash = finish_minimizer_rep_hash(forward, reverse);
            while length != 0 && window[(head + length - 1) & 31].1 > hash {
                length -= 1;
            }
            window[(head + length) & 31] = (next_pos, hash);
            length += 1;
            next_pos += 1;
        }
        if length == 0 {
            continue;
        }
        let minimum = window[head].1;
        // With an unchanged minimum, all old ties were emitted already. Only
        // the one newly admitted position can add another tied minimum.
        let first = if previous_minimum == Some(minimum) {
            length - 1
        } else {
            0
        };
        previous_minimum = Some(minimum);
        for slot in first..length {
            let (pos, hash) = window[(head + slot) & 31];
            if hash != minimum {
                break;
            }
            if last_selected.is_none_or(|last| pos > last) {
                let key = encode_minimizer_rep(&seq[pos..pos + g])
                    .ok_or_else(|| Error::InvalidFormat("invalid minimizer".into()))?;
                add_minimizer_position(global_base + pos as u64, key, keys, bitmap_blocks);
                last_selected = Some(pos);
            }
        }
    }
    Ok(())
}

pub(super) fn add_minimizer_position(
    pos: u64,
    bytes: [u8; 8],
    keys: &mut Vec<u64>,
    bitmap_blocks: &mut [Vec<u32>],
) {
    keys.push(u64::from_be_bytes(bytes));
    bitmap_blocks[(pos >> 32) as usize].push(pos as u32);
}

pub(super) fn minimizer_position_count(seq_len: usize, g: usize) -> usize {
    if g == 0 || seq_len < g {
        0
    } else {
        seq_len - g + 1
    }
}

pub(super) fn bifrost_minimizer_bounds(k: usize, g: usize) -> Option<RangeInclusive<usize>> {
    let shift = 1usize;
    if k < g + shift + 1 {
        return None;
    }
    let end = k.checked_sub(g + shift)?;
    if end < shift {
        None
    } else {
        Some(RangeInclusive {
            start: shift,
            last: end,
        })
    }
}

pub(super) fn minimizer_rep_hash_words(seq: &[u8]) -> Result<(u64, u64)> {
    if (1..=32).contains(&seq.len()) {
        let (forward, reverse) = crate::dna::encode_kmer_pair(seq)
            .ok_or_else(|| Error::InvalidFormat("invalid minimizer".into()))?;
        return Ok((
            grouped_state(forward, seq.len()),
            grouped_state(reverse, seq.len()),
        ));
    }
    let mut forward = 0u64;
    let mut reverse = 0u64;
    for i in 0..seq.len() {
        let f = HASH_CODES[seq[i] as usize];
        if f == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid minimizer".into()));
        }
        forward = forward.rotate_left(1) ^ REP_HASH_VALS[f as usize];

        let r = HASH_COMP_CODES[seq[seq.len() - 1 - i] as usize];
        if r == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid minimizer".into()));
        }
        reverse = reverse.rotate_left(1) ^ REP_HASH_VALS[r as usize];
    }
    Ok((forward, reverse))
}

pub(super) fn finish_minimizer_rep_hash(forward: u64, reverse: u64) -> u64 {
    let mut buf = [0u8; 16];
    buf[..8].copy_from_slice(&forward.min(reverse).to_le_bytes());
    buf[8..].copy_from_slice(&forward.max(reverse).to_le_bytes());
    wyhash(&buf, 0)
}

#[cfg(test)]
pub(super) fn minimizer_rep_hashes_into(seq: &[u8], g: usize, hashes: &mut Vec<u64>) -> Result<()> {
    let (mut forward, mut reverse) = minimizer_rep_hash_words(&seq[..g])?;
    let remove_forward = REP_HASH_VALS.map(|value| value.rotate_left(g as u32));
    let remove_reverse = REP_HASH_VALS.map(|value| value.rotate_right(1));
    let add_reverse = REP_HASH_VALS.map(|value| value.rotate_left((g - 1) as u32));
    hashes.clear();
    hashes.reserve(seq.len() - g + 1);
    hashes.push(finish_minimizer_rep_hash(forward, reverse));
    for pos in g..seq.len() {
        let old = seq[pos - g] as usize;
        let new = seq[pos] as usize;
        let new_forward = HASH_CODES[new];
        if new_forward == INVALID_BASE_CODE {
            return Err(Error::InvalidFormat("invalid minimizer".into()));
        }
        // Remove the outgoing base and add the incoming base in both orientations.
        // This is the same REP hash as re-encoding the complete overlapping window.
        forward = forward.rotate_left(1)
            ^ remove_forward[HASH_CODES[old] as usize]
            ^ REP_HASH_VALS[new_forward as usize];
        reverse = reverse.rotate_right(1)
            ^ remove_reverse[HASH_COMP_CODES[old] as usize]
            ^ add_reverse[HASH_COMP_CODES[new] as usize];
        hashes.push(finish_minimizer_rep_hash(forward, reverse));
    }
    Ok(())
}

#[cfg(test)]
pub(super) fn minimizer_rep_hashes(seq: &[u8], g: usize) -> Result<Vec<u64>> {
    let mut hashes = Vec::new();
    minimizer_rep_hashes_into(seq, g, &mut hashes)?;
    Ok(hashes)
}

#[cfg(test)]
pub(super) fn minimizer_rep_hash(seq: &[u8]) -> Result<u64> {
    let (forward, reverse) = minimizer_rep_hash_words(seq)?;
    Ok(finish_minimizer_rep_hash(forward, reverse))
}
