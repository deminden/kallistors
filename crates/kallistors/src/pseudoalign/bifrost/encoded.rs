//! Packed graph sequences and deferred transcript metadata.

use super::*;

pub(super) enum SequenceMetadata {
    Compact(Vec<u64>),
    // Transient BUS layout: low bit marks a complete left-aligned 31-mer.
    // Other entries contain a 31-bit word offset and the original 32-bit length.
    InlineShort(Vec<u64>),
    Wide(Vec<(usize, usize)>),
}

impl SequenceMetadata {
    pub(super) fn retire(self) {
        match self {
            Self::Compact(v) | Self::InlineShort(v) => crate::index::scratch_pages::retire(v),
            Self::Wide(v) => crate::index::scratch_pages::retire(v),
        }
    }

    pub(super) fn new(offsets: Vec<usize>, lengths: Vec<usize>) -> Self {
        if offsets.last().copied().unwrap_or(0) <= u32::MAX as usize
            && lengths.iter().all(|&len| len <= u32::MAX as usize)
        {
            // A u64 collection can otherwise reuse the geometrically grown
            // usize input allocation. Keep only exact metadata capacity and
            // return the consumed scratch pages before mapping starts.
            let pairs = offsets.iter().zip(&lengths);
            let mut records = Vec::with_capacity(pairs.len());
            records.extend(pairs.map(|(&start, &len)| start as u64 | ((len as u64) << 32)));
            crate::index::scratch_pages::retire(offsets);
            crate::index::scratch_pages::retire(lengths);
            Self::Compact(records)
        } else {
            Self::Wide(offsets.into_iter().zip(lengths).collect())
        }
    }

    #[inline]
    pub(super) fn get(&self, index: usize) -> Option<(usize, usize)> {
        match self {
            Self::Compact(values) => values
                .get(index)
                .map(|&entry| (entry as u32 as usize, (entry >> 32) as usize)),
            // Inline words have no offset in the packed sequence arena.
            // EncodedBaseStore handles their sequence view and length directly.
            Self::InlineShort(values) => {
                let &entry = values.get(index)?;
                (entry & 1 == 0).then_some((
                    ((entry >> 1) & 0x7fff_ffff) as usize,
                    (entry >> 32) as usize,
                ))
            }
            Self::Wide(values) => values.get(index).copied(),
        }
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::Compact(v) | Self::InlineShort(v) => v.len(),
            Self::Wide(v) => v.len(),
        }
    }
    pub(super) fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub(crate) struct EncodedBaseStore {
    pub(super) metadata: SequenceMetadata,
    pub(super) values: Vec<u64>,
}

pub(crate) struct EncodedSequence<'a> {
    pub(super) words: &'a [u64],
    pub(crate) len: usize,
}

#[inline]
fn canonical_minimizer_bytes(word: u64, g: usize) -> [u8; 8] {
    let shift = (32 - g) * 2;
    let forward = word << shift;
    let reversed = word.reverse_bits();
    let reverse = !(((reversed & 0x5555_5555_5555_5555) << 1)
        | ((reversed & 0xaaaa_aaaa_aaaa_aaaa) >> 1))
        & (u64::MAX << shift);
    forward.min(reverse).to_le_bytes()
}

impl EncodedSequence<'_> {
    #[inline]
    pub(super) fn minimizer_at(&self, start: usize, g: usize) -> Option<[u8; 8]> {
        if g == 0 || g > 32 {
            return None;
        }
        Some(canonical_minimizer_bytes(self.kmer_at(start, g)?, g))
    }

    #[inline]
    pub(crate) fn kmer_at(&self, start: usize, k: usize) -> Option<u64> {
        if k > 32 || start.checked_add(k)? > self.len {
            return None;
        }
        if k == 0 {
            return Some(0);
        }
        let offset = start / 32;
        let within = start % 32;
        let shift = within * 2;
        // SAFETY: k is nonzero and start+k <= len. The sequence owns
        // ceil(len/32) words; the crossing guard also bounds the second word.
        let mut word = unsafe { *self.words.get_unchecked(offset) } << shift;
        if within + k > 32 {
            word |= unsafe { *self.words.get_unchecked(offset + 1) } >> (64 - shift);
        }
        Some(word >> ((32 - k) * 2))
    }

    #[inline]
    pub(crate) fn range(&self, first: usize, last: usize, k: usize) -> Option<PackedKmerRange> {
        if !(1..=32).contains(&k) || first > last || last.checked_add(k)? > self.len {
            return None;
        }
        let span = last.checked_sub(first)?.checked_add(k)?;
        if span > 64 {
            return None;
        }
        let word_idx = first / 32;
        let within = first % 32;
        // SAFETY: first <= last, k >= 1, and last+k <= len bound every
        // touched word; the span guards below select only required words.
        let first_word = unsafe { *self.words.get_unchecked(word_idx) };
        let second = if within + span > 32 {
            unsafe { *self.words.get_unchecked(word_idx + 1) }
        } else {
            0
        };
        let third = if within + span > 64 {
            unsafe { *self.words.get_unchecked(word_idx + 2) }
        } else {
            0
        };
        let shift = within * 2;
        let (high, low) = if shift == 0 {
            (first_word, second)
        } else {
            (
                (first_word << shift) | (second >> (64 - shift)),
                (second << shift) | (third >> (64 - shift)),
            )
        };
        Some(PackedKmerRange {
            high,
            low,
            right_shift: (64 - 2 * k) as u32,
            next: first,
            end: last,
        })
    }
}

pub(crate) struct PackedKmerRange {
    pub(super) high: u64,
    pub(super) low: u64,
    pub(super) right_shift: u32,
    pub(super) next: usize,
    pub(super) end: usize,
}

impl PackedKmerRange {
    #[inline]
    pub(crate) fn first_match(
        mut self,
        forward: u64,
        reverse: u64,
        allow_forward: bool,
        allow_reverse: bool,
    ) -> Option<(usize, bool)> {
        #[cfg(target_arch = "aarch64")]
        if self.right_shift == 2 {
            use std::arch::aarch64::*;
            // SAFETY: NEON is baseline on AArch64. All vector reads/writes use
            // complete local arrays; reference sequence words were bounded by
            // range construction. No vector access reads outside that range.
            unsafe {
                let shifts01 = vld1q_s64([0, 2].as_ptr());
                let shifts23 = vld1q_s64([4, 6].as_ptr());
                let right01 = vld1q_s64([-64, -62].as_ptr());
                let right23 = vld1q_s64([-60, -58].as_ptr());
                let indices = vld1q_u32([0, 1, 2, 3].as_ptr());
                let wanted_forward = vdupq_n_u64(forward);
                let wanted_reverse = vdupq_n_u64(reverse);
                while self.next <= self.end && self.end - self.next >= 3 {
                    let high = vdupq_n_u64(self.high);
                    let low = vdupq_n_u64(self.low);
                    let first = vshrq_n_u64::<2>(vorrq_u64(
                        vshlq_u64(high, shifts01),
                        vshlq_u64(low, right01),
                    ));
                    let second = vshrq_n_u64::<2>(vorrq_u64(
                        vshlq_u64(high, shifts23),
                        vshlq_u64(low, right23),
                    ));
                    let f0 = if allow_forward {
                        vceqq_u64(first, wanted_forward)
                    } else {
                        vdupq_n_u64(0)
                    };
                    let f1 = if allow_forward {
                        vceqq_u64(second, wanted_forward)
                    } else {
                        vdupq_n_u64(0)
                    };
                    let r0 = if allow_reverse {
                        vceqq_u64(first, wanted_reverse)
                    } else {
                        vdupq_n_u64(0)
                    };
                    let r1 = if allow_reverse {
                        vceqq_u64(second, wanted_reverse)
                    } else {
                        vdupq_n_u64(0)
                    };
                    let mask =
                        vcombine_u32(vmovn_u64(vorrq_u64(f0, r0)), vmovn_u64(vorrq_u64(f1, r1)));
                    // Most candidate windows miss. Reduce to an earliest lane
                    // only when at least one of the four exact comparisons hit.
                    if vmaxvq_u32(mask) != 0 {
                        let offset = vminvq_u32(vbslq_u32(mask, indices, vdupq_n_u32(u32::MAX)));
                        if offset < 4 {
                            let mut codes = [0u64; 4];
                            vst1q_u64(codes.as_mut_ptr(), first);
                            vst1q_u64(codes.as_mut_ptr().add(2), second);
                            return Some((
                                self.next + offset as usize,
                                !allow_forward || codes[offset as usize] != forward,
                            ));
                        }
                    }
                    self.high = (self.high << 8) | (self.low >> 56);
                    self.low <<= 8;
                    self.next += 4;
                }
            }
        }
        for (start, code) in self {
            if allow_forward && code == forward {
                return Some((start, false));
            }
            if allow_reverse && code == reverse {
                return Some((start, true));
            }
        }
        None
    }
}

impl Iterator for PackedKmerRange {
    type Item = (usize, u64);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.next > self.end {
            return None;
        }
        let result = (self.next, self.high >> self.right_shift);
        self.next += 1;
        self.high = (self.high << 2) | (self.low >> 62);
        self.low <<= 2;
        Some(result)
    }
}

impl EncodedBaseStore {
    pub(super) fn from_fixed_words(values: Vec<u64>, k: usize) -> Self {
        let metadata = if values.len() <= u32::MAX as usize {
            SequenceMetadata::Compact(
                (0..values.len())
                    .map(|i| i as u64 | ((k as u64) << 32))
                    .collect(),
            )
        } else {
            SequenceMetadata::Wide((0..values.len()).map(|i| (i, k)).collect())
        };
        Self { metadata, values }
    }

    pub(super) fn retire(self) {
        self.metadata.retire();
        crate::index::scratch_pages::retire(self.values);
    }

    /// Rewrite the existing metadata allocation only when every backed offset
    /// fits losslessly. Short words use one otherwise unused low padding bit;
    /// no valid substring can observe the marker. The original arena stays
    /// intact, so long-sequence offsets and all serialized index data are stable.
    pub(super) fn inline_short_words(&mut self) {
        let SequenceMetadata::Compact(records) = &mut self.metadata else {
            return;
        };
        if !records.iter().any(|&entry| entry >> 32 == 31)
            || records
                .iter()
                .any(|&entry| entry >> 32 != 31 && entry as u32 > 0x7fff_ffff)
        {
            return;
        }
        for entry in records.iter_mut() {
            let offset = *entry as u32 as usize;
            let len = (*entry >> 32) as usize;
            *entry = if len == 31 {
                self.values[offset] | 1
            } else {
                ((offset as u64) << 1) | ((len as u64) << 32)
            };
        }
        let records = std::mem::take(records);
        self.metadata = SequenceMetadata::InlineShort(records);
    }

    #[inline]
    pub(crate) fn sequence(&self, id: usize) -> Option<EncodedSequence<'_>> {
        if let SequenceMetadata::InlineShort(records) = &self.metadata
            && *records.get(id)? & 1 != 0
        {
            return Some(EncodedSequence {
                words: &records[id..id + 1],
                len: 31,
            });
        }
        let (start, len) = self.metadata.get(id)?;
        Some(EncodedSequence {
            words: &self.values[start..start + len.div_ceil(32)],
            len,
        })
    }
    #[cfg(test)]
    pub(crate) fn kmer_range(
        &self,
        seq_idx: usize,
        first: usize,
        last: usize,
        k: usize,
    ) -> Option<PackedKmerRange> {
        self.sequence(seq_idx)?.range(first, last, k)
    }

    pub(super) fn from_sequences(seqs: &[Vec<u8>], pool: Option<&ThreadPool>) -> Self {
        let mut offsets = Vec::with_capacity(seqs.len() + 1);
        offsets.push(0);
        let mut total = 0usize;
        for seq in seqs {
            total += seq.len().div_ceil(32);
            offsets.push(total);
        }
        let lengths = seqs.iter().map(Vec::len).collect();

        let build_sequential = || {
            let mut values = vec![0u64; total];
            for (seq_idx, seq) in seqs.iter().enumerate() {
                let start = offsets[seq_idx];
                let end = offsets[seq_idx + 1];
                if start == end {
                    continue;
                }
                pack_bases_into(seq, &mut values[start..end]);
            }
            values
        };

        let values = if let Some(pool) = pool {
            let mut values = vec![0u64; total];
            let base_ptr = values.as_mut_ptr() as usize;
            pool.install(|| {
                (0..seqs.len())
                    .into_par_iter()
                    .with_min_len(2048)
                    .for_each(|seq_idx| {
                        let start = offsets[seq_idx];
                        let end = offsets[seq_idx + 1];
                        if start == end {
                            return;
                        }
                        // SAFETY: each sequence owns a disjoint prefix-summed slice.
                        let out = unsafe {
                            std::slice::from_raw_parts_mut(
                                (base_ptr as *mut u64).add(start),
                                end - start,
                            )
                        };
                        pack_bases_into(&seqs[seq_idx], out);
                    });
            });
            values
        } else {
            build_sequential()
        };

        Self {
            metadata: SequenceMetadata::new(offsets, lengths),
            values,
        }
    }

    #[inline]
    // Loader sequences are decoded from the two-bit index and contain only A/C/G/T.
    // Derive only bitmap-selected minimizers; no per-position table is retained.
    pub(super) fn minimizer_at(&self, seq_idx: usize, start: usize, g: usize) -> Option<[u8; 8]> {
        if g == 0 || g > 32 {
            return None;
        }
        let word = self.encode_kmer_at(seq_idx, start, g)?;
        Some(canonical_minimizer_bytes(word, g))
    }

    #[inline]
    pub(crate) fn len(&self, seq_idx: usize) -> Option<usize> {
        if let SequenceMetadata::InlineShort(records) = &self.metadata
            && *records.get(seq_idx)? & 1 != 0
        {
            return Some(31);
        }
        self.metadata.get(seq_idx).map(|(_, len)| len)
    }

    #[inline]
    pub(crate) fn encode_kmer_at(&self, seq_idx: usize, start: usize, k: usize) -> Option<u64> {
        if matches!(&self.metadata, SequenceMetadata::InlineShort(_)) {
            return self.sequence(seq_idx)?.kmer_at(start, k);
        }
        if k > 32 {
            return None;
        }
        let (seq_start, seq_len) = self.metadata.get(seq_idx)?;
        if start.checked_add(k)? > seq_len {
            return None;
        }
        if k == 0 {
            return Some(0);
        }
        let word_idx = seq_start + start / 32;
        let within_word = start % 32;
        let shift = within_word * 2;
        let mut word = self.values[word_idx] << shift;
        if within_word + k > 32 {
            word |= self.values[word_idx + 1] >> (64 - shift);
        }
        Some(word >> ((32 - k) * 2))
    }
}

pub(crate) struct EncodedUnitigs {
    pub(super) g: usize,
    pub(crate) unitig_bases: EncodedBaseStore,
    pub(crate) km_bases: EncodedBaseStore,
}

impl EncodedUnitigs {
    #[inline]
    pub(super) fn unitig_minimizer_count(&self, seq_idx: usize) -> usize {
        minimizer_count(
            self.unitig_bases.len(seq_idx).expect("valid unitig"),
            self.g,
        )
    }

    #[inline]
    pub(super) fn unitig_minimizer_at(&self, seq_idx: usize, pos: usize) -> Option<[u8; 8]> {
        self.unitig_bases.minimizer_at(seq_idx, pos, self.g)
    }
}
