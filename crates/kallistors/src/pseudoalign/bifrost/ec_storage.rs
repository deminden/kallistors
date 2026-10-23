//! Compact EC blocks, fallback k-mers and minimizer offsets.

use super::*;

#[derive(Clone, Copy)]
pub(super) struct CompactEcBlock {
    pub(super) lb: u32,
    pub(super) ub: u32,
    pub(super) ec_start: u32,
    pub(super) ec_end: u32,
}

#[derive(Clone, Copy)]
pub(super) struct WideEcBlock {
    pub(super) lb: u32,
    pub(super) ub: u32,
    pub(super) ec_start: usize,
    pub(super) ec_end: usize,
}

pub(super) enum FlatEcBlocks {
    Compact {
        values: Vec<u64>,
        overflow: Vec<CompactEcBlock>,
    },
    Wide(Vec<WideEcBlock>),
}

impl FlatEcBlocks {
    pub(super) fn retire(self) {
        match self {
            Self::Compact { values, overflow } => {
                crate::index::scratch_pages::retire(values);
                crate::index::scratch_pages::retire(overflow);
            }
            Self::Wide(v) => crate::index::scratch_pages::retire(v),
        }
    }

    pub(super) fn new() -> Self {
        Self::Compact {
            values: Vec::new(),
            overflow: Vec::new(),
        }
    }
    pub(super) fn len(&self) -> usize {
        match self {
            Self::Compact { values, .. } => values.len(),
            Self::Wide(v) => v.len(),
        }
    }
    pub(super) fn push(&mut self, lb: u32, ub: u32, ec_start: usize, ec_end: usize) {
        if ec_end > u32::MAX as usize && matches!(self, Self::Compact { .. }) {
            let wide = (0..self.len())
                .map(|i| self.get(i).expect("valid EC block"))
                .collect();
            *self = Self::Wide(wide);
        }
        match self {
            Self::Compact { values, overflow } => {
                let len = ec_end - ec_start;
                let packed = if lb < 1 << 12 && ub < 1 << 12 && ec_start < 1 << 26 && len < 1 << 13
                {
                    u64::from(lb)
                        | (u64::from(ub) << 12)
                        | ((ec_start as u64) << 24)
                        | ((len as u64) << 50)
                } else {
                    let id = overflow.len() as u64;
                    overflow.push(CompactEcBlock {
                        lb,
                        ub,
                        ec_start: ec_start as u32,
                        ec_end: ec_end as u32,
                    });
                    (1u64 << 63) | id
                };
                values.push(packed);
            }
            Self::Wide(values) => values.push(WideEcBlock {
                lb,
                ub,
                ec_start,
                ec_end,
            }),
        }
    }
    #[inline]
    pub(super) fn get(&self, index: usize) -> Option<WideEcBlock> {
        match self {
            Self::Compact { values, overflow } => {
                let &packed = values.get(index)?;
                if packed >> 63 != 0 {
                    let b = overflow.get((packed & !(1u64 << 63)) as usize)?;
                    Some(WideEcBlock {
                        lb: b.lb,
                        ub: b.ub,
                        ec_start: b.ec_start as usize,
                        ec_end: b.ec_end as usize,
                    })
                } else {
                    let start = ((packed >> 24) & ((1 << 26) - 1)) as usize;
                    Some(WideEcBlock {
                        lb: (packed & 4095) as u32,
                        ub: ((packed >> 12) & 4095) as u32,
                        ec_start: start,
                        ec_end: start + ((packed >> 50) & 8191) as usize,
                    })
                }
            }
            Self::Wide(values) => values.get(index).copied(),
        }
    }
}

pub(crate) struct FlatEcIndex {
    pub(crate) owned_ec_blocks: Vec<Box<[crate::index::EcBlock]>>,
    pub(super) unitig_block_offsets: MinzOffsets,
    pub(super) blocks: FlatEcBlocks,
    pub(super) ec_values: Vec<u32>,
    pub(super) last_ec: std::ops::Range<usize>,
}

impl FlatEcIndex {
    pub(super) fn new(unitig_capacity: usize) -> Self {
        let mut offsets = Vec::with_capacity(unitig_capacity + 1);
        offsets.push(0);
        Self {
            unitig_block_offsets: MinzOffsets::Compact(offsets),
            owned_ec_blocks: Vec::new(),
            blocks: FlatEcBlocks::new(),
            ec_values: Vec::new(),
            last_ec: 0..0,
        }
    }

    pub(super) fn finish_ec_block(&mut self, lb: u32, ub: u32, start: usize) {
        let end = self.ec_values.len();
        let previous = self.last_ec.clone();
        if end - start == previous.len()
            && self.ec_values[previous.clone()] == self.ec_values[start..end]
        {
            self.ec_values.truncate(start);
            self.blocks.push(lb, ub, previous.start, previous.end);
        } else {
            self.last_ec = start..end;
            self.blocks.push(lb, ub, start, end);
        }
    }

    pub(super) fn push_buffered_node<R: BufRead>(
        &mut self,
        reader: &mut R,
        size: usize,
        scratch: &mut Vec<u8>,
    ) -> Result<()> {
        let available = reader.fill_buf()?;
        if available.len() >= size {
            self.push_serialized_node(&available[..size])?;
            reader.consume(size);
        } else {
            scratch.resize(size, 0);
            reader.read_exact(scratch)?;
            self.push_serialized_node(scratch)?;
        }
        Ok(())
    }

    pub(super) fn push_serialized_node(&mut self, node: &[u8]) -> Result<()> {
        let blocks = node.get(4..).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "failed to fill whole buffer",
            )
        })?;
        crate::index::read_block_array_slice(blocks, |lb, ub, bitmap| {
            let start = self.ec_values.len();
            unsafe { crate::index::append_roaring_to_vec(bitmap, &mut self.ec_values) }
                .ok_or_else(|| {
                    Error::InvalidFormat("failed to deserialize roaring bitmap".into())
                })?;
            self.finish_ec_block(lb, ub, start);
            Ok(())
        })?;
        self.unitig_block_offsets.push(self.blocks.len());
        Ok(())
    }

    pub(super) fn push_unitig_blocks(&mut self, blocks: &[crate::index::EcBlock]) {
        for block in blocks {
            let start = self.ec_values.len();
            self.ec_values.extend_from_slice(&block.ec);
            self.finish_ec_block(block.lb, block.ub, start);
        }
        self.unitig_block_offsets.push(self.blocks.len());
    }

    pub(crate) fn block_index_for_position(&self, unitig_id: usize, pos: usize) -> Option<usize> {
        if !self.owned_ec_blocks.is_empty() {
            let blocks = self.owned_ec_blocks.get(unitig_id)?;
            let mut lo = 0;
            let mut hi = blocks.len();
            let pos = pos as u32;
            while lo < hi {
                let mid = (lo + hi) / 2;
                if pos < blocks[mid].lb {
                    hi = mid;
                } else if pos >= blocks[mid].ub {
                    lo = mid + 1;
                } else {
                    return Some(mid);
                }
            }
            return None;
        }
        let start = self.unitig_block_offsets.value_at(unitig_id)?;
        let end = self
            .unitig_block_offsets
            .value_at(unitig_id.checked_add(1)?)?;
        let pos = pos as u32;
        let mut lo = start;
        let mut hi = end;
        while lo < hi {
            let mid = (lo + hi) / 2;
            let block = self.blocks.get(mid)?;
            if pos < block.lb {
                hi = mid;
            } else if pos >= block.ub {
                lo = mid + 1;
            } else {
                return Some(mid - start);
            }
        }
        None
    }

    pub(crate) fn ec(&self, unitig_id: usize, block_idx: usize) -> &[u32] {
        if !self.owned_ec_blocks.is_empty() {
            return &self.owned_ec_blocks[unitig_id][block_idx].ec;
        }
        let global_idx = self
            .unitig_block_offsets
            .value_at(unitig_id)
            .expect("valid unitig")
            + block_idx;
        let block = self.blocks.get(global_idx).expect("valid EC block");
        &self.ec_values[block.ec_start..block.ec_end]
    }

    pub(crate) fn block_bounds(&self, unitig_id: usize, block_idx: usize) -> Option<(u32, u32)> {
        if !self.owned_ec_blocks.is_empty() {
            let block = self.owned_ec_blocks.get(unitig_id)?.get(block_idx)?;
            return Some((block.lb, block.ub));
        }
        let start = self.unitig_block_offsets.value_at(unitig_id)?;
        let end = self
            .unitig_block_offsets
            .value_at(unitig_id.checked_add(1)?)?;
        let global_idx = start.checked_add(block_idx)?;
        if global_idx >= end {
            return None;
        }
        let block = self.blocks.get(global_idx)?;
        Some((block.lb, block.ub))
    }
}

pub(crate) struct SelectiveFallbackKmerIndex {
    pub(super) entries: HashMap<u64, (usize, u64)>,
    pub(super) duplicates: HashMap<u64, Vec<(usize, u64)>>,
}

impl SelectiveFallbackKmerIndex {
    pub(super) fn from_codes(
        codes: impl Iterator<Item = (usize, u64, u64)>,
        capacity: usize,
    ) -> Self {
        let mut result = Self {
            entries: HashMap::with_capacity(capacity),
            duplicates: HashMap::new(),
        };
        for (id, fwd, rev) in codes {
            let key = fwd.min(rev);
            let value = (id, fwd);
            match result.entries.entry(key) {
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(value);
                }
                std::collections::hash_map::Entry::Occupied(slot) => {
                    result
                        .duplicates
                        .entry(key)
                        .or_insert_with(|| vec![*slot.get()])
                        .push(value);
                }
            }
        }
        result
    }
    pub(super) fn build(unitigs_len: usize, km_unitigs: &[Vec<u8>], h_kmers: &[Vec<u8>]) -> Self {
        let codes = km_unitigs
            .iter()
            .chain(h_kmers)
            .enumerate()
            .filter_map(|(id, seq)| {
                encode_kmer_pair(seq).map(|(fwd, rev)| (unitigs_len + id, fwd, rev))
            });
        Self::from_codes(codes, km_unitigs.len() + h_kmers.len())
    }
    pub(super) fn build_packed(
        unitigs_len: usize,
        km: &EncodedBaseStore,
        h_words: &[u64],
        k: usize,
    ) -> Self {
        let count = km.metadata.len();
        let km_codes = (0..count).map(|id| {
            let fwd = km.encode_kmer_at(id, 0, k).expect("validated fixed k-mer");
            (
                unitigs_len + id,
                fwd,
                crate::dna::reverse_complement_word(fwd, k),
            )
        });
        let h_codes = h_words.iter().enumerate().map(|(id, &word)| {
            let fwd = word.checked_shr((64 - 2 * k) as u32).unwrap_or(0);
            (
                unitigs_len + count + id,
                fwd,
                crate::dna::reverse_complement_word(fwd, k),
            )
        });
        Self::from_codes(km_codes.chain(h_codes), count + h_words.len())
    }

    pub(crate) fn lookup(&self, code: u64) -> Option<&[(usize, u64)]> {
        let entry = self.entries.get(&code)?;
        if let Some(values) = self.duplicates.get(&code) {
            Some(values)
        } else {
            Some(std::slice::from_ref(entry))
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub(super) enum MinzOffsets {
    Compact(Vec<u32>),
    Wide(Vec<usize>),
}

impl MinzOffsets {
    pub(super) fn retire(self) {
        match self {
            Self::Compact(v) => crate::index::scratch_pages::retire(v),
            Self::Wide(v) => crate::index::scratch_pages::retire(v),
        }
    }

    #[inline]
    pub(super) fn value_at(&self, index: usize) -> Option<usize> {
        match self {
            Self::Compact(v) => v.get(index).map(|&x| x as usize),
            Self::Wide(v) => v.get(index).copied(),
        }
    }

    pub(super) fn push(&mut self, value: usize) {
        if let Self::Compact(values) = self {
            if value <= u32::MAX as usize {
                values.push(value as u32);
                return;
            }
            *self = Self::Wide(values.iter().map(|&x| x as usize).collect());
        }
        if let Self::Wide(values) = self {
            values.push(value);
        }
    }

    pub(super) fn from_counts(mut counts: Vec<usize>) -> (Self, usize) {
        // Compute prefixes in the original buffer before choosing a narrower layout.
        let mut total = 0usize;
        for count in &mut counts {
            let next = total + *count;
            *count = total;
            total = next;
        }
        counts.push(total);
        let offsets = if total <= u32::MAX as usize {
            let mut compact = vec![0u32; counts.len()];
            for (slot, offset) in compact.iter_mut().zip(counts) {
                *slot = offset as u32;
            }
            Self::Compact(compact)
        } else {
            Self::Wide(counts)
        };
        (offsets, total)
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::Compact(v) => v.len(),
            Self::Wide(v) => v.len(),
        }
    }

    #[inline]
    pub(super) fn bounds(&self, idx: usize) -> (usize, usize) {
        match self {
            Self::Compact(v) => (v[idx] as usize, v[idx + 1] as usize),
            Self::Wide(v) => (v[idx], v[idx + 1]),
        }
    }

    pub(super) fn cursors(&self, buckets: usize) -> Vec<usize> {
        match self {
            Self::Compact(v) => v[..buckets].iter().map(|&p| p as usize).collect(),
            Self::Wide(v) => v[..buckets].to_vec(),
        }
    }
}

pub struct MinzPositionIndex {
    pub(super) offsets: MinzOffsets,
    pub(super) values: Vec<u64>,
    pub(super) inline: Option<Vec<u64>>,
    pub(super) inline_building: bool,
}

impl MinzPositionIndex {
    pub(super) fn retire(self) {
        self.offsets.retire();
        crate::index::scratch_pages::retire(self.values);
        if let Some(entries) = self.inline {
            crate::index::scratch_pages::retire(entries);
        }
    }

    pub(super) fn from_counts(counts: Vec<usize>) -> Self {
        let (offsets, total) = MinzOffsets::from_counts(counts);
        Self {
            offsets,
            values: vec![0; total],
            inline: None,
            inline_building: false,
        }
    }

    pub(super) fn from_counts_inline(
        mut counts: Vec<usize>,
        escapes: &[usize],
        enabled: bool,
    ) -> Self {
        if !enabled || usize::BITS != 64 {
            return Self::from_counts(counts);
        }
        let total: usize = counts.iter().sum();
        // Bit61 stores the complete special-marker predicate. Keep the
        // existing total/offset domain, falling back only for a per-bucket
        // count that cannot fit its29-bit field.
        if total > 0x3fff_ffff || counts.iter().any(|&count| count > 0x1fff_ffff) {
            return Self::from_counts(counts);
        }
        let mut overflow = 0usize;
        let mut next_escape = 0usize;
        for (bucket, count) in counts.iter_mut().enumerate() {
            while next_escape < escapes.len() && escapes[next_escape] < bucket {
                next_escape += 1;
            }
            let escaped = escapes.get(next_escape) == Some(&bucket);
            let n = *count;
            if n == 1 && !escaped {
                *count = 0;
            } else {
                *count = ((1u64 << 63) | ((n as u64) << 32) | overflow as u64) as usize;
                overflow += n;
            }
        }
        let entries = counts.into_iter().map(|x| x as u64).collect();
        Self {
            offsets: MinzOffsets::Compact(Vec::new()),
            values: vec![0; overflow],
            inline: Some(entries),
            inline_building: true,
        }
    }

    pub(super) fn len(&self) -> usize {
        self.inline
            .as_ref()
            .map_or_else(|| self.offsets.len().saturating_sub(1), Vec::len)
    }

    pub(super) fn for_each_large_bucket(
        &self,
        minimum: usize,
        mut visit: impl FnMut(usize, &[u64]),
    ) {
        if let Some(entries) = &self.inline {
            let scanned;
            #[cfg(target_arch = "aarch64")]
            {
                let mut offset_base = 0;
                // Most descriptors are small. Reject eight at once using a
                // conservative unsigned bound, then recheck every selected
                // descriptor scalarly so the overcrowded flag cannot affect
                // inclusion, ordering or the supported minimum==0 case.
                let (blocks, _) = entries.as_chunks::<8>();
                for block in blocks {
                    if inline_block_may_contain_large(block, minimum) {
                        for (offset, &entry) in block.iter().enumerate() {
                            if entry >> 63 != 0 && ((entry >> 32) & 0x1fff_ffff) as usize >= minimum
                            {
                                visit(offset_base + offset, self.get(offset_base + offset));
                            }
                        }
                    }
                    offset_base += 8;
                }
                scanned = offset_base;
            }
            #[cfg(not(target_arch = "aarch64"))]
            {
                scanned = 0;
            }
            for (bucket, &entry) in entries.iter().enumerate().skip(scanned) {
                if entry >> 63 != 0 && ((entry >> 32) & 0x1fff_ffff) as usize >= minimum {
                    visit(bucket, self.get(bucket));
                }
            }
        } else {
            for bucket in 0..self.len() {
                let positions = self.get(bucket);
                if positions.len() >= minimum {
                    visit(bucket, positions);
                }
            }
        }
    }

    pub(super) fn cursors(&self) -> Vec<usize> {
        if self.inline.is_some() {
            Vec::new()
        } else {
            self.offsets.cursors(self.len())
        }
    }

    pub fn get(&self, idx: usize) -> &[u64] {
        if let Some(entries) = &self.inline {
            debug_assert!(!self.inline_building);
            let entry = entries[idx];
            if entry >> 63 == 0 {
                &entries[idx..idx + 1]
            } else {
                let start = entry as u32 as usize;
                let len = ((entry >> 32) & 0x1fff_ffff) as usize;
                &self.values[start..start + len]
            }
        } else {
            let (start, end) = self.offsets.bounds(idx);
            &self.values[start..end]
        }
    }

    /// Return the existing positions and their precomputed special-marker flag.
    /// The flag occupies an unused descriptor bit; CSR indexes retain the scan.
    #[inline]
    pub(crate) fn get_with_overcrowded(&self, idx: usize) -> (&[u64], bool) {
        let positions = self.get(idx);
        let overcrowded = if let Some(entries) = &self.inline {
            let entry = entries[idx];
            entry >> 63 != 0 && entry & (1u64 << 62) != 0
        } else {
            positions
                .iter()
                .any(|&pos_id| (pos_id >> 32) as u32 == u32::MAX && (pos_id & 0x8000_0000) != 0)
        };
        (positions, overcrowded)
    }

    /// Complete special-marker predicate, including zero-count markers.
    /// Both summary bits are transient; serialized index bytes stay unchanged.
    #[inline]
    pub(crate) fn get_with_special(&self, idx: usize) -> (&[u64], bool) {
        let positions = self.get(idx);
        let special = if let Some(entries) = &self.inline {
            let entry = entries[idx];
            entry >> 63 != 0 && entry & (1u64 << 61) != 0
        } else {
            positions
                .iter()
                .any(|&position| (position >> 32) as u32 == u32::MAX)
        };
        (positions, special)
    }

    #[cfg(target_arch = "aarch64")]
    #[inline]
    pub(super) fn prefetch_inline_descriptor(&self, idx: usize) {
        if let Some(entries) = &self.inline
            && let Some(entry) = entries.get(idx)
        {
            // The descriptor stays allocated during the ordered fill. This
            // write-intent hint does not access or change any Rust value.
            unsafe {
                core::arch::asm!("prfm pstl1keep, [{address}]",
                address = in(reg) entry as *const u64,
                options(nostack, readonly, preserves_flags));
            }
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[inline]
    pub(super) fn prefetch_inline_overflow(&self, idx: usize) {
        if let Some(entries) = &self.inline
            && let Some(&entry) = entries.get(idx)
            && entry >> 63 != 0
            && let Some(value) = self.values.get(entry as u32 as usize)
        {
            // A future descriptor may advance before use; prefetching an older
            // in-bounds slot is harmless. Ordered writes remain in push_at.
            unsafe {
                core::arch::asm!("prfm pstl1keep, [{address}]",
                address = in(reg) value as *const u64,
                options(nostack, readonly, preserves_flags));
            }
        }
    }

    pub(super) fn push_at(&mut self, cursors: &mut [usize], idx: usize, pos_id: u64) {
        if let Some(entries) = &mut self.inline {
            debug_assert!(self.inline_building);
            let entry = &mut entries[idx];
            if *entry >> 63 == 0 {
                debug_assert_eq!(
                    pos_id >> 63,
                    0,
                    "top-bit singleton requires overflow storage"
                );
                *entry = pos_id;
            } else {
                self.values[*entry as u32 as usize] = pos_id;
                *entry += 1;
                if (pos_id >> 32) as u32 == u32::MAX {
                    *entry |= 1u64 << 61;
                    if (pos_id & 0x8000_0000) != 0 {
                        *entry |= 1u64 << 62;
                    }
                }
            }
        } else {
            let cursor = &mut cursors[idx];
            self.values[*cursor] = pos_id;
            *cursor += 1;
        }
    }

    pub(super) fn finish_inline(&mut self) {
        if let Some(entries) = &mut self.inline {
            for entry in entries {
                if *entry >> 63 != 0 {
                    *entry -= (*entry >> 32) & 0x1fff_ffff;
                }
            }
            self.inline_building = false;
        }
    }
}

impl Index<usize> for MinzPositionIndex {
    type Output = [u64];

    fn index(&self, index: usize) -> &Self::Output {
        self.get(index)
    }
}

#[cfg(target_arch = "aarch64")]
#[inline]
fn inline_block_may_contain_large(block: &[u64], minimum: usize) -> bool {
    use std::arch::aarch64::*;
    if minimum > 0x1fff_ffff {
        return false;
    }
    debug_assert_eq!(block.len(), 8);
    // SAFETY: the caller passes a complete eight-u64 slice. All four unaligned
    // loads are in bounds; NEON is baseline on AArch64. This only skips blocks
    // with no possible match; exact descriptor checks remain in the caller.
    unsafe {
        let threshold = vdupq_n_u64((1u64 << 63) | ((minimum as u64) << 32));
        let a = vcgeq_u64(vld1q_u64(block.as_ptr()), threshold);
        let b = vcgeq_u64(vld1q_u64(block.as_ptr().add(2)), threshold);
        let c = vcgeq_u64(vld1q_u64(block.as_ptr().add(4)), threshold);
        let d = vcgeq_u64(vld1q_u64(block.as_ptr().add(6)), threshold);
        let any = vorrq_u64(vorrq_u64(a, b), vorrq_u64(c, d));
        vgetq_lane_u64::<0>(any) != 0 || vgetq_lane_u64::<1>(any) != 0
    }
}

#[cfg(test)]
mod large_bucket_tests {
    use super::*;

    #[test]
    #[cfg(target_pointer_width = "64")]
    fn large_bucket_scan_preserves_flags_boundaries_tail_and_csr_order() {
        let counts = [1usize, 2, 63, 64, 65, 1, 1, 64, 2, 1, 65];
        let mut index = MinzPositionIndex::from_counts_inline(counts.to_vec(), &[], true);
        let mut cursors = index.cursors();
        for (bucket, &count) in counts.iter().enumerate() {
            for ordinal in 0..count {
                let position = if count > 1 && ordinal == 0 && bucket % 2 == 0 {
                    (u64::from(u32::MAX) << 32) | 0x8000_0000
                } else {
                    ((bucket as u64) << 32) | ordinal as u64
                };
                index.push_at(&mut cursors, bucket, position);
            }
        }
        index.finish_inline();
        for minimum in [0, 1, 2, 63, 64, 65, 66, 0x4000_0000, usize::MAX] {
            let expected: Vec<_> = (0..index.len())
                .filter(|&bucket| {
                    let entry = index.inline.as_ref().unwrap()[bucket];
                    entry >> 63 != 0 && ((entry >> 32) & 0x1fff_ffff) as usize >= minimum
                })
                .map(|bucket| (bucket, index.get(bucket).to_vec()))
                .collect();
            let mut actual = Vec::new();
            index.for_each_large_bucket(minimum, |bucket, values| {
                actual.push((bucket, values.to_vec()))
            });
            assert_eq!(actual, expected, "minimum={minimum}");
        }
        let csr = MinzPositionIndex::from_counts(counts.to_vec());
        for minimum in [0, 1, 2, 64, 66] {
            let expected: Vec<_> = (0..csr.len())
                .filter(|&bucket| csr.get(bucket).len() >= minimum)
                .collect();
            let mut actual = Vec::new();
            csr.for_each_large_bucket(minimum, |bucket, _| actual.push(bucket));
            assert_eq!(actual, expected);
        }
    }
}
