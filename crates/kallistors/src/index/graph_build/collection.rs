//! Canonical k-mer collection and bounded parallel merging.

use super::*;

pub(super) struct KmerStore {
    pub(super) sorted: Vec<u64>,
}

// Keep the private word oracle independent of destructive scratch retirement.
#[cfg(test)]
pub(super) fn collect_kmers(
    transcripts: &[Transcript],
    k: usize,
    threads: usize,
) -> Result<KmerStore> {
    collect_kmers_and_pack(&mut transcripts.to_vec(), k, threads).map(|(words, _)| words)
}

pub(super) fn collect_kmers_and_pack(
    transcripts: &mut [Transcript],
    k: usize,
    threads: usize,
) -> Result<(KmerStore, Vec<PackedSequence>)> {
    let total_positions = transcripts
        .iter()
        .map(|tr| kmer_position_count(tr.seq.len(), k))
        .sum::<usize>();
    if threads > 1 && transcripts.len() > 1 && total_positions >= PARALLEL_KMER_COLLECTION_MIN {
        let pool = IndexWorkers::new(threads)?;
        let mut sorted = Vec::with_capacity(total_positions);
        struct Scratch {
            words: Vec<u64>,
            written: usize,
            packed: Vec<(usize, Vec<PackedSequence>)>,
        }
        let buffers = {
            let collected = std::sync::Mutex::new(&mut sorted);
            // Moving disjoint chunks out of this iterator is safe; workers do
            // all packing only after their group's borrowed anchors are gone.
            let groups = std::sync::Mutex::new(
                bounded_transcript_groups(transcripts, k, 1 << 24).enumerate(),
            );
            pool.install(|| {
                rayon::broadcast(|_| -> Result<Scratch> {
                    let mut words = Vec::new();
                    let mut written = 0usize;
                    let mut packed = Vec::new();
                    let mut prefixes = AnchorMap::default();
                    let mut suffixes = AnchorMap::default();
                    loop {
                        let Some((group, chunk)) = groups.lock().unwrap().next() else {
                            break;
                        };
                        let count = chunk
                            .iter()
                            .map(|tr| kmer_position_count(tr.seq.len(), k))
                            .sum();
                        if words.capacity() < count {
                            super::super::scratch_pages::discard_unused(&mut words, written);
                            drop(std::mem::take(&mut words));
                            words = Vec::with_capacity(count);
                            written = 0;
                        }
                        prefixes.clear();
                        suffixes.clear();
                        // These values borrow this group. A fresh map bounds their
                        // lifetime before its ASCII sequences are packed and freed.
                        let mut interior = AnchorMap::default();
                        for i in 0..chunk.len() {
                            let sequence = anchored_transcript_region(
                                chunk,
                                i,
                                k,
                                &mut prefixes,
                                &mut suffixes,
                            );
                            insert_interior_unique_kmers::<4>(
                                &mut words,
                                sequence,
                                k,
                                i,
                                &mut interior,
                            )?;
                        }
                        drop(interior);
                        let group_packed = pack_transcript_sequences(chunk, 1)?;
                        if words.len()
                            >= total_positions
                                .div_ceil(threads)
                                .max(PARALLEL_KMER_COLLECTION_MIN)
                        {
                            words.radix_sort_unstable();
                        } else {
                            words
                                .radix_sort_builder()
                                .with_parallel(false)
                                .with_single_threaded_tuner()
                                .sort();
                        }
                        written = written.max(words.len());
                        words.dedup();
                        collected.lock().unwrap().extend_from_slice(&words);
                        packed.push((group, group_packed));
                        words.clear();
                    }
                    Ok(Scratch {
                        words,
                        written,
                        packed,
                    })
                })
            })
        };
        let mut groups = Vec::new();
        for buffer in buffers {
            let Scratch {
                mut words,
                written,
                packed,
            } = buffer?;
            super::super::scratch_pages::discard_unused(&mut words, written);
            groups.extend(packed);
        }
        groups.sort_unstable_by_key(|(group, _)| *group);
        let packed = concatenate_shards(groups.into_iter().map(|(_, packed)| packed).collect());
        // All original ASCII groups and their borrowed maps are gone now.
        // Reclaim on their owners before the global sorting peak, once only.
        pool.collect_unused(&sorted);
        pool.install(|| sorted.radix_sort_unstable());
        let written = sorted.len();
        sorted.dedup();
        super::super::scratch_pages::discard_unused(&mut sorted, written);
        return Ok((KmerStore { sorted }, packed));
    }
    let mut sorted = Vec::with_capacity(total_positions);
    for transcript in transcripts.iter() {
        insert_transcript_kmers(&mut sorted, &transcript.seq, k)?;
    }
    sorted.sort_unstable();
    sorted.dedup();
    sorted.shrink_to_fit();
    let packed = pack_transcript_sequences(transcripts, threads)?;
    if sorted.len() >= PARALLEL_KMER_COLLECTION_MIN {
        IndexWorkers::new(threads)?.collect_unused(&sorted);
    }
    Ok((KmerStore { sorted }, packed))
}

#[cfg(test)]
pub(super) fn novel_transcript_region<'a>(seq: &'a [u8], previous: &[u8], k: usize) -> &'a [u8] {
    &seq[novel_transcript_bounds(seq, previous, k)]
}

#[cfg(test)]
pub(super) fn novel_transcript_history<'a>(
    seq: &'a [u8],
    previous: &[Transcript],
    k: usize,
) -> &'a [u8] {
    if seq.len() < k {
        return seq;
    }
    let mut first = 0;
    let mut end = seq.len();
    for prior in previous.iter().rev().take(2) {
        let bounds = novel_transcript_bounds(seq, &prior.seq, k);
        first = first.max(bounds.start);
        end = end.min(bounds.end);
        if first + k > end {
            return &seq[..0];
        }
    }
    &seq[first..end]
}

pub(super) fn anchored_transcript_region<'a>(
    chunk: &'a [Transcript],
    i: usize,
    k: usize,
    prefixes: &mut AnchorMap<usize>,
    suffixes: &mut AnchorMap<usize>,
) -> &'a [u8] {
    let seq = &chunk[i].seq;
    if seq.len() < k {
        return seq;
    }
    let prefix = if seq.len() >= 32 {
        crate::dna::encode_kmer_pair(&seq[..32]).map(|pair| pair.0)
    } else {
        None
    };
    let suffix = if seq.len() >= 32 {
        crate::dna::encode_kmer_pair(&seq[seq.len() - 32..]).map(|pair| pair.0)
    } else {
        None
    };
    let candidates = [
        i.checked_sub(1),
        prefix.and_then(|key| prefixes.get(&key).copied()),
        suffix.and_then(|key| suffixes.get(&key).copied()),
    ];
    let mut first = 0;
    let mut end = seq.len();
    for (slot, id) in candidates.iter().copied().enumerate() {
        let Some(id) = id else {
            continue;
        };
        if candidates[..slot].contains(&Some(id)) {
            continue;
        }
        let bounds = novel_transcript_bounds(seq, &chunk[id].seq, k);
        first = first.max(bounds.start);
        end = end.min(bounds.end);
    }
    if let Some(key) = prefix {
        prefixes.insert(key, i);
    }
    if let Some(key) = suffix {
        suffixes.insert(key, i);
    }
    if first + k > end {
        &seq[..0]
    } else {
        &seq[first..end]
    }
}

pub(super) fn shared_prefix_bytes(a: &[u8], b: &[u8]) -> usize {
    let length = a.len().min(b.len());
    let mut matched = 0;
    while matched + 32 <= length && a[matched..matched + 32] == b[matched..matched + 32] {
        matched += 32;
    }
    while matched + 8 <= length {
        let difference = u64::from_le_bytes(a[matched..matched + 8].try_into().unwrap())
            ^ u64::from_le_bytes(b[matched..matched + 8].try_into().unwrap());
        if difference != 0 {
            return matched + difference.trailing_zeros() as usize / 8;
        }
        matched += 8;
    }
    matched
        + a[matched..]
            .iter()
            .zip(&b[matched..])
            .take_while(|(a, b)| a == b)
            .count()
}

pub(super) fn shared_suffix_bytes(a: &[u8], b: &[u8]) -> usize {
    let length = a.len().min(b.len());
    let mut matched = 0;
    while matched + 32 <= length
        && a[a.len() - matched - 32..a.len() - matched]
            == b[b.len() - matched - 32..b.len() - matched]
    {
        matched += 32;
    }
    while matched + 8 <= length {
        let difference = u64::from_le_bytes(
            a[a.len() - matched - 8..a.len() - matched]
                .try_into()
                .unwrap(),
        ) ^ u64::from_le_bytes(
            b[b.len() - matched - 8..b.len() - matched]
                .try_into()
                .unwrap(),
        );
        if difference != 0 {
            return matched + difference.leading_zeros() as usize / 8;
        }
        matched += 8;
    }
    matched
        + a[..a.len() - matched]
            .iter()
            .rev()
            .zip(b[..b.len() - matched].iter().rev())
            .take_while(|(a, b)| a == b)
            .count()
}

pub(super) fn novel_transcript_bounds(
    seq: &[u8],
    previous: &[u8],
    k: usize,
) -> std::ops::Range<usize> {
    if seq.len() < k || previous.len() < k {
        return 0..seq.len();
    }
    let prefix = shared_prefix_bytes(seq, previous);
    if prefix == seq.len() {
        return 0..0;
    }
    let suffix = shared_suffix_bytes(seq, previous);
    let first = prefix.saturating_sub(k - 1);
    let end = kmer_position_count(seq.len(), k).saturating_sub(suffix.saturating_sub(k - 1));
    if first >= end {
        0..0
    } else {
        first..end + k - 1
    }
}

// Anchors choose only candidates. Complete byte equality proves each skipped
// interval; every seam window remains and every transcript retains its metadata.
pub(super) fn insert_interior_unique_kmers<'a, const ANCHOR_BITS: u32>(
    entries: &mut Vec<u64>,
    seq: &'a [u8],
    k: usize,
    transcript: usize,
    anchors: &mut AnchorMap<(usize, &'a [u8])>,
) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    if seq
        .chunks(32)
        .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none())
    {
        return Err(Error::InvalidFormat("invalid DNA base".into()));
    }
    let mask = kmer_mask(k);
    let shift = 64 - 2 * k;
    let (mut forward, mut reverse) = encode_kmer_pair(&seq[..k])?;
    let mut pos = 0;
    let anchor_mask = (1u64 << ANCHOR_BITS) - 1;
    while pos <= seq.len() - k {
        let hash = (forward ^ (forward >> 27)).wrapping_mul(0x9e3779b97f4a7c15) >> 32;
        if hash & anchor_mask == 0 {
            let (prior_id, prior) = *anchors.entry(forward).or_insert((transcript, &seq[pos..]));
            if prior_id < transcript {
                let current = &seq[pos..];
                let matched = shared_prefix_bytes(current, prior);
                if matched >= k {
                    pos += matched - k + 1;
                    if pos > seq.len() - k {
                        break;
                    }
                    (forward, reverse) = encode_kmer_pair(&seq[pos..pos + k])?;
                    continue;
                }
            }
        }
        entries.push(forward.min(reverse));
        pos += 1;
        if pos > seq.len() - k {
            break;
        }
        let base = seq[pos + k - 1];
        let code = u64::from(((base >> 1) ^ (base >> 2)) & 3);
        forward = ((forward << 2) & mask) | (code << shift);
        reverse = ((reverse >> 2) | ((3 - code) << 62)) & mask;
    }
    Ok(())
}

pub(super) fn insert_transcript_kmers(entries: &mut Vec<u64>, seq: &[u8], k: usize) -> Result<()> {
    if seq.len() < k {
        return Ok(());
    }
    let mask = kmer_mask(k);
    let shift = 64 - 2 * k;
    let (mut fwd, mut rev) = encode_kmer_pair(&seq[..k])?;
    entries.push(fwd.min(rev));
    // Transcript sequences have already been normalized by the FASTA loader.
    // Validate the tail once, allowing the hot rolling loop to use bit arithmetic.
    if seq[k..]
        .chunks(32)
        .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none())
    {
        return Err(Error::InvalidFormat("invalid DNA base".into()));
    }
    entries.extend(seq[k..].iter().map(|&base| {
        let code = u64::from(((base >> 1) ^ (base >> 2)) & 3);
        fwd = ((fwd << 2) & mask) | (code << shift);
        rev = ((rev >> 2) | ((3 - code) << 62)) & mask;
        fwd.min(rev)
    }));

    Ok(())
}

// Keep ordinary groups at their established target count while bounding raw
// scratch for long transcripts. A single oversized transcript remains intact,
// preserving every seam window; the bound is budget plus one transcript.
pub(super) fn bounded_transcript_groups(
    transcripts: &mut [Transcript],
    k: usize,
    position_budget: usize,
) -> impl Iterator<Item = &mut [Transcript]> {
    let mut targets = 0usize;
    let mut positions = 0usize;
    transcripts.split_inclusive_mut(move |transcript| {
        targets += 1;
        positions = positions.saturating_add(kmer_position_count(transcript.seq.len(), k));
        if targets >= 2048 || positions >= position_budget {
            targets = 0;
            positions = 0;
            true
        } else {
            false
        }
    })
}

pub(super) fn kmer_position_count(seq_len: usize, k: usize) -> usize {
    if seq_len >= k { seq_len - k + 1 } else { 0 }
}
