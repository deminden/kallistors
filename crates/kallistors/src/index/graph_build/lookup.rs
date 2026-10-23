//! Exact rank lookup, visited claims and cache hints.

use super::*;

#[derive(Clone, Copy, Default)]
pub(super) struct KmerLoc(pub(super) u64);

pub(super) struct KmerLookup {
    pub(super) directory: KmerDirectory,
    pub(super) shift: usize,
    pub(super) padding: usize,
    pub(super) suffix_bits: usize,
    pub(super) suffix_mask: u64,
    pub(super) rank_bits: usize,
    pub(super) rank_mask: u64,
    pub(super) padding_mask: u64,
    pub(super) pair_bits: usize,
    pub(super) pair_mask: u64,
    pub(super) triple_bits: usize,
    pub(super) triple_mask: u64,
    pub(super) quad_bits: usize,
    pub(super) quad_mask: u64,
}

pub(super) enum KmerDirectory {
    Ranks(RankDirectory),
    Inline(Vec<u64>),
}

pub(super) struct RankDirectory {
    pub(super) delta_bits: usize,
    pub(super) bases: Vec<u32>,
    pub(super) deltas: Vec<u8>,
    pub(super) wide: Vec<u32>,
}

impl RankDirectory {
    pub(super) fn new(offsets: Vec<u32>, parallel: bool) -> Self {
        if offsets.last().copied().unwrap_or(0) >= (1 << 31) {
            return Self {
                delta_bits: 12,
                bases: Vec::new(),
                deltas: Vec::new(),
                wide: offsets,
            };
        }
        if offsets.last().copied().unwrap_or(0) > 1 << 28 {
            Self::packed::<16>(offsets, parallel)
        } else {
            Self::packed::<12>(offsets, parallel)
        }
    }
    pub(super) fn packed<const BITS: usize>(offsets: Vec<u32>, parallel: bool) -> Self {
        let mut bases = Vec::with_capacity(offsets.len().div_ceil(256));
        let mut deltas = vec![0u8; offsets.len().div_ceil(2) * (BITS / 4) + 1];
        let mut wide = Vec::new();
        for chunk in offsets.chunks(256) {
            let base = chunk[0];
            if chunk[chunk.len() - 1] - base < (1 << BITS) {
                bases.push(base);
            } else {
                bases.push((1 << 31) | wide.len() as u32);
                wide.extend_from_slice(chunk);
            }
        }
        let fill = |(dest, source): (&mut [u8], &[u32])| {
            let base = source[0];
            if source[source.len() - 1] - base < (1 << BITS) {
                for (dest, pair) in dest.chunks_mut(BITS / 4).zip(source.chunks(2)) {
                    let first = pair[0] - base;
                    let second = pair.get(1).map_or(0, |value| value - base);
                    let packed = (first | (second << BITS)).to_le_bytes();
                    dest.copy_from_slice(&packed[..BITS / 4]);
                }
            }
        };
        if parallel {
            deltas
                .par_chunks_mut(256 * BITS / 8)
                .zip(offsets.par_chunks(256))
                .for_each(fill);
        } else {
            deltas
                .chunks_mut(256 * BITS / 8)
                .zip(offsets.chunks(256))
                .for_each(fill);
        }
        // The raw prefix table has been copied exactly into compact blocks.
        // Retire its initialized pages before allocating the edge cache.
        super::super::scratch_pages::retire(offsets);
        Self {
            delta_bits: BITS,
            bases,
            deltas,
            wide,
        }
    }
    #[inline]
    pub(super) fn prefetch(&self, bucket: usize) {
        if !self.deltas.is_empty() {
            let first = if self.delta_bits == 16 {
                bucket * 2
            } else {
                (bucket >> 1) * 3
            };
            prefetch_rank_delta(&self.deltas[first]);
            prefetch_rank_delta(&self.bases[bucket >> 8]);
        }
    }
    #[inline]
    pub(super) fn get(&self, bucket: usize) -> usize {
        if self.bases.is_empty() {
            return self.wide[bucket] as usize;
        }
        let base = self.bases[bucket >> 8];
        if base & (1 << 31) != 0 {
            self.wide[(base & !(1 << 31)) as usize + (bucket & 255)] as usize
        } else {
            let first = if self.delta_bits == 16 {
                bucket * 2
            } else {
                (bucket >> 1) * 3 + (bucket & 1)
            };
            let value = u16::from_le_bytes(
                self.deltas[first..first + 2]
                    .try_into()
                    .expect("packed rank pair"),
            );
            let value = if self.delta_bits == 16 {
                value
            } else {
                (value >> ((bucket & 1) * 4)) & 4095
            };
            (base + u32::from(value)) as usize
        }
    }
    #[inline(always)]
    pub(super) fn range(&self, bucket: usize) -> (usize, usize) {
        if self.bases.is_empty() {
            let pair = &self.wide[bucket..bucket + 2];
            return (pair[0] as usize, pair[1] as usize);
        }
        if bucket & 255 == 255 {
            return (self.get(bucket), self.get(bucket + 1));
        }
        let base = self.bases[bucket >> 8];
        if base & (1 << 31) != 0 {
            let index = (base & !(1 << 31)) as usize + (bucket & 255);
            let pair = &self.wide[index..index + 2];
            (pair[0] as usize, pair[1] as usize)
        } else {
            let first = if self.delta_bits == 16 {
                bucket * 2
            } else {
                (bucket >> 1) * 3 + (bucket & 1)
            };
            let packed = u32::from_le_bytes(
                self.deltas[first..first + 4]
                    .try_into()
                    .expect("packed rank range"),
            );
            let shift = if self.delta_bits == 16 {
                0
            } else {
                (bucket & 1) * 4
            };
            let mask = (1 << self.delta_bits) - 1;
            (
                (base + ((packed >> shift) & mask)) as usize,
                (base + ((packed >> (self.delta_bits + shift)) & mask)) as usize,
            )
        }
    }
}

pub(super) struct VisitedKmers {
    pub(super) words: Vec<u64>,
}

impl KmerLoc {
    pub(super) fn new(unitig_id: usize, offset: usize, unitig_is_canonical: bool) -> Result<Self> {
        let unitig_id = u32::try_from(unitig_id)
            .map_err(|_| Error::InvalidFormat("too many unitigs for v13 index".into()))?;
        if offset >= (1usize << 31) {
            return Err(Error::InvalidFormat(
                "unitig is too long for packed v13 index metadata".into(),
            ));
        }
        Ok(Self(
            ((unitig_id as u64) << 32) | ((offset as u64) << 1) | u64::from(unitig_is_canonical),
        ))
    }

    pub(super) fn unitig_id(self) -> usize {
        (self.0 >> 32) as usize
    }

    pub(super) fn offset(self) -> usize {
        ((self.0 >> 1) & 0x7fff_ffff) as usize
    }

    pub(super) fn unitig_is_canonical(self) -> bool {
        self.0 & 1 != 0
    }
}

impl VisitedKmers {
    pub(super) fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
        }
    }

    pub(super) fn insert_index(&mut self, idx: usize) -> bool {
        let word = &mut self.words[idx >> 6];
        let mask = 1u64 << (idx & 63);
        let previous = *word;
        // Re-setting an existing bit is inert. Publish the same exact claim
        // without a conditional store in the ordered hot loop.
        *word = previous | mask;
        previous & mask == 0
    }

    pub(super) fn contains_index(&self, idx: usize) -> bool {
        self.words[idx >> 6] & (1u64 << (idx & 63)) != 0
    }
}

#[inline(always)]
pub(super) fn prefetch_rank_delta<T>(delta: &T) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: the live in-bounds reference supplies a cache hint only.
    unsafe {
        core::arch::asm!("prfm pldl1keep, [{address}]",address=in(reg) delta as *const T,options(nostack,readonly,preserves_flags));
    }
    #[cfg(not(target_arch = "aarch64"))]
    let _ = delta;
}

#[inline(always)]
pub(super) fn prefetch_kmer(word: &u64) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: The reference identifies an in-bounds live key. PRFM is a
    // nonblocking cache hint; it neither reads a Rust value nor changes memory.
    unsafe {
        core::arch::asm!("prfm pldl1keep, [{address}]",address=in(reg) word as *const u64,options(nostack,readonly,preserves_flags));
    }
    #[cfg(not(target_arch = "aarch64"))]
    let _ = word;
}

impl KmerLookup {
    #[cfg(test)]
    pub(super) fn new(entries: &[u64]) -> Result<Self> {
        Self::with_threads(entries, 1, 32)
    }

    pub(super) fn with_threads(entries: &[u64], threads: usize, k: usize) -> Result<Self> {
        if entries.len() > u32::MAX as usize {
            return Err(Error::InvalidFormat(
                "too many distinct k-mers for compact lookup".into(),
            ));
        }
        let target = entries.len().div_ceil(2).max(1);
        let bits = (usize::BITS as usize - target.saturating_sub(1).leading_zeros() as usize)
            .clamp(MIN_KMER_LOOKUP_PREFIX_BITS, MAX_KMER_LOOKUP_PREFIX_BITS);
        let suffix_bits = (2 * k).saturating_sub(bits);
        let padding = 64 - 2 * k;
        let rank_bits =
            usize::BITS as usize - entries.len().saturating_sub(1).leading_zeros() as usize;
        // Keep the directory within512MiB and leave the descriptor's high bit
        // free to distinguish a complete singleton key from a rank range.
        let inline = rank_bits + suffix_bits <= 63
            && entries.len() < (1 << 27)
            && (1usize << bits) <= 1 << 26;
        let prefix = if inline {
            bits
        } else {
            kmer_lookup_prefix_bits(entries.len())
        };
        let shift = 64 - prefix;
        let count = 1usize << prefix;
        let suffix_mask = if suffix_bits == 0 {
            0
        } else {
            (1u64 << suffix_bits) - 1
        };
        let pair_bits = (61 - rank_bits) / 2;
        let pair_mask = (1u64 << pair_bits) - 1;
        let triple_bits = (60 - rank_bits) / 3;
        let triple_mask = (1u64 << triple_bits) - 1;
        let quad_bits = (59 - rank_bits) / 4;
        let quad_mask = (1u64 << quad_bits) - 1;
        let pool = if threads > 1 && entries.len() >= PARALLEL_KMER_COLLECTION_MIN {
            Some(IndexWorkers::new(threads)?)
        } else {
            None
        };
        let directory = if inline {
            let mut slots = vec![0u64; count];
            let fill = |(chunk, dest): (usize, &mut [u64])| {
                let first = chunk << 12;
                let mut index = entries.partition_point(|&word| (word >> shift) < first as u64);
                for (local, slot) in dest.iter_mut().enumerate() {
                    let start = index;
                    let end_bucket = (first + local + 1) as u64;
                    while index < entries.len() && (entries[index] >> shift) < end_bucket {
                        index += 1;
                    }
                    *slot = if index == start {
                        0
                    } else if index == start + 1 {
                        (1u64 << 63)
                            | ((start as u64) << suffix_bits)
                            | ((entries[start] >> padding) & suffix_mask)
                    } else if index == start + 2 {
                        let first = (entries[start] >> padding) & pair_mask;
                        let second = (entries[start + 1] >> padding) & pair_mask;
                        (1u64 << 61)
                            | start as u64
                            | (first << rank_bits)
                            | (second << (rank_bits + pair_bits))
                    } else if index == start + 3 {
                        let mut descriptor = (1u64 << 60) | start as u64;
                        for slot in 0..3 {
                            descriptor |= ((entries[start + slot] >> padding) & triple_mask)
                                << (rank_bits + slot * triple_bits);
                        }
                        descriptor
                    } else if index == start + 4 {
                        let mut descriptor = (1u64 << 59) | start as u64;
                        for slot in 0..4 {
                            descriptor |= ((entries[start + slot] >> padding) & quad_mask)
                                << (rank_bits + slot * quad_bits);
                        }
                        descriptor
                    } else {
                        ((index as u64) << 32) | start as u64
                    };
                }
            };
            if let Some(pool) = &pool {
                pool.install(|| slots.par_chunks_mut(1 << 12).enumerate().for_each(fill));
            } else {
                slots.chunks_mut(1 << 12).enumerate().for_each(fill);
            }
            KmerDirectory::Inline(slots)
        } else {
            let mut offsets = vec![0u32; count + 1];
            let fill = |(chunk, dest): (usize, &mut [u32])| {
                let first = chunk << 18;
                let mut index = entries.partition_point(|&word| (word >> shift) < first as u64);
                for (local, slot) in dest.iter_mut().enumerate() {
                    let bucket = (first + local) as u64;
                    while index < entries.len() && (entries[index] >> shift) < bucket {
                        index += 1;
                    }
                    *slot = index as u32;
                }
            };
            if let Some(pool) = &pool {
                pool.install(|| offsets.par_chunks_mut(1 << 18).enumerate().for_each(fill));
            } else {
                offsets.chunks_mut(1 << 18).enumerate().for_each(fill);
            }
            KmerDirectory::Ranks(if let Some(pool) = &pool {
                pool.install(|| RankDirectory::new(offsets, true))
            } else {
                RankDirectory::new(offsets, false)
            })
        };
        Ok(Self {
            directory,
            shift,
            padding,
            suffix_bits,
            suffix_mask,
            rank_bits,
            rank_mask: (1u64 << rank_bits) - 1,
            padding_mask: (1u64 << padding) - 1,
            pair_bits,
            pair_mask,
            triple_bits,
            triple_mask,
            quad_bits,
            quad_mask,
        })
    }

    #[inline]
    pub(super) fn search(entries: &[u64], word: u64, start: usize, end: usize) -> Option<usize> {
        let slice = &entries[start..end];
        if slice.len() == 1 {
            return (slice[0] == word).then_some(start);
        }
        if slice.len() == 2 {
            return if slice[0] == word {
                Some(start)
            } else {
                (slice[1] == word).then_some(start + 1)
            };
        }
        if slice.len() <= LINEAR_KMER_LOOKUP_MAX {
            slice
                .iter()
                .position(|&entry| entry == word)
                .map(|idx| start + idx)
        } else {
            slice.binary_search(&word).ok().map(|idx| start + idx)
        }
    }

    #[inline]
    pub(super) fn guarded_result<const N: usize>(
        &self,
        entries: &[u64],
        word: u64,
        descriptor: u64,
        bits: usize,
        mask: u64,
    ) -> Option<usize> {
        let start = (descriptor & self.rank_mask) as usize;
        let guard = (word >> self.padding) & mask;
        for slot in 0..N {
            if guard == (descriptor >> (self.rank_bits + slot * bits)) & mask
                && entries[start + slot] == word
            {
                return Some(start + slot);
            }
        }
        None
    }

    #[inline(always)]
    pub(super) fn inline_result(
        &self,
        entries: &[u64],
        word: u64,
        descriptor: u64,
    ) -> Option<usize> {
        if descriptor == 0 {
            return None;
        }
        match descriptor.leading_zeros() {
            0 => {
                if word & self.padding_mask != 0
                    || (word >> self.padding) & self.suffix_mask != descriptor & self.suffix_mask
                {
                    return None;
                }
                Some(((descriptor & !(1 << 63)) >> self.suffix_bits) as usize)
            }
            2 => {
                let start = (descriptor & self.rank_mask) as usize;
                let fingerprint = (word >> self.padding) & self.pair_mask;
                if fingerprint == (descriptor >> self.rank_bits) & self.pair_mask
                    && entries[start] == word
                {
                    return Some(start);
                }
                if fingerprint == (descriptor >> (self.rank_bits + self.pair_bits)) & self.pair_mask
                    && entries[start + 1] == word
                {
                    return Some(start + 1);
                }
                None
            }
            3 => self.guarded_result::<3>(
                entries,
                word,
                descriptor,
                self.triple_bits,
                self.triple_mask,
            ),
            4 => {
                self.guarded_result::<4>(entries, word, descriptor, self.quad_bits, self.quad_mask)
            }
            _ => Self::search(
                entries,
                word,
                (descriptor as u32) as usize,
                (descriptor >> 32) as usize,
            ),
        }
    }

    #[inline]
    pub(super) fn find_many<const N: usize>(
        &self,
        entries: &[u64],
        words: [u64; N],
    ) -> [Option<usize>; N] {
        self.find_many_count(entries, words, N)
    }

    #[inline(always)]
    pub(super) fn find_many_count<const N: usize>(
        &self,
        entries: &[u64],
        words: [u64; N],
        count: usize,
    ) -> [Option<usize>; N] {
        debug_assert!(count <= N);
        let mut found = [None; N];
        match &self.directory {
            KmerDirectory::Inline(slots) => {
                let mut descriptors = [0u64; N];
                for i in 0..count {
                    descriptors[i] = slots[(words[i] >> self.shift) as usize];
                }
                for i in 0..count {
                    found[i] = self.inline_result(entries, words[i], descriptors[i]);
                }
            }
            KmerDirectory::Ranks(offsets) => {
                for &word in &words[..count] {
                    offsets.prefetch((word >> self.shift) as usize);
                }
                let mut bounds = [(0usize, 0usize); N];
                for i in 0..count {
                    let bucket = (words[i] >> self.shift) as usize;
                    bounds[i] = offsets.range(bucket);
                }
                // Hint independent exact key loads before searching each range.
                // Empty buckets issue no hint; large buckets start at midpoint.
                for &(first, end) in &bounds[..count] {
                    if first < end {
                        let index = if end - first > LINEAR_KMER_LOOKUP_MAX {
                            first + (end - first) / 2
                        } else {
                            first
                        };
                        prefetch_kmer(&entries[index]);
                    }
                }
                for i in 0..count {
                    found[i] = Self::search(entries, words[i], bounds[i].0, bounds[i].1);
                }
            }
        }
        found
    }

    pub(super) fn find(&self, entries: &[u64], word: u64) -> Option<usize> {
        let bucket = (word >> self.shift) as usize;
        match &self.directory {
            KmerDirectory::Inline(slots) => self.inline_result(entries, word, slots[bucket]),
            KmerDirectory::Ranks(offsets) => {
                let (start, end) = offsets.range(bucket);
                Self::search(entries, word, start, end)
            }
        }
    }
}

pub(super) fn kmer_lookup_prefix_bits(entry_count: usize) -> usize {
    if entry_count == 0 {
        return MIN_KMER_LOOKUP_PREFIX_BITS;
    }
    let target_buckets = entry_count.div_ceil(TARGET_KMERS_PER_LOOKUP_BUCKET);
    let bits = usize::BITS as usize - target_buckets.saturating_sub(1).leading_zeros() as usize;
    bits.clamp(MIN_KMER_LOOKUP_PREFIX_BITS, MAX_KMER_LOOKUP_PREFIX_BITS)
}
