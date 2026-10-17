use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::ops::Index;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::index::bifrost::{BooPhf, read_bitcontainer_values};
use crate::index::read_block_array_blocks;
use crate::timing::{self, Stage};
use crate::{Error, Result};

use super::ec::encode_kmer_pair;
use super::io::{read_i32_le, read_u32_le, read_u64_le};
use super::{KmerEcIndex, read_onlist_with_lengths};

pub type KmerPosIndex = HashMap<u64, Vec<(usize, usize, bool)>>;

const LOAD_BATCH_POSITIONS: usize = 1 << 20;
const LOAD_MIN_PAR_CHUNK: usize = 1 << 14;
const INVALID_BASE_CODE: u8 = 4;
const BASE_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'g' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes[b't' as usize] = 3;
    codes
};

const DECODED_BYTE_BASES: [[u8; 4]; 256] = {
    let mut table = [[0; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut pos = 0;
        while pos < 4 {
            table[byte][pos] = b"ACGT"[(byte >> (pos * 2)) & 3];
            pos += 1;
        }
        byte += 1;
    }
    table
};

enum SequenceMetadata {
    Compact(Vec<(u32, u32)>),
    Wide(Vec<(usize, usize)>),
}

impl SequenceMetadata {
    fn new(offsets: Vec<usize>, lengths: Vec<usize>) -> Self {
        if offsets.last().copied().unwrap_or(0) <= u32::MAX as usize
            && lengths.iter().all(|&len| len <= u32::MAX as usize)
        {
            Self::Compact(
                offsets
                    .into_iter()
                    .zip(lengths)
                    .map(|(start, len)| (start as u32, len as u32))
                    .collect(),
            )
        } else {
            Self::Wide(offsets.into_iter().zip(lengths).collect())
        }
    }

    #[inline]
    fn get(&self, index: usize) -> Option<(usize, usize)> {
        match self {
            Self::Compact(values) => values
                .get(index)
                .map(|&(start, len)| (start as usize, len as usize)),
            Self::Wide(values) => values.get(index).copied(),
        }
    }

    fn is_empty(&self) -> bool {
        match self {
            Self::Compact(v) => v.is_empty(),
            Self::Wide(v) => v.is_empty(),
        }
    }
}

pub(crate) struct EncodedBaseStore {
    metadata: SequenceMetadata,
    values: Vec<u64>,
}

pub(crate) struct PackedKmerRange {
    high: u64,
    low: u64,
    right_shift: u32,
    next: usize,
    end: usize,
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
    pub(crate) fn kmer_range(
        &self,
        seq_idx: usize,
        first: usize,
        last: usize,
        k: usize,
    ) -> Option<PackedKmerRange> {
        if !(1..=32).contains(&k) || first > last {
            return None;
        }
        let span = last.checked_sub(first)?.checked_add(k)?;
        if span > 64 {
            return None;
        }
        let (seq_start, seq_len) = self.metadata.get(seq_idx)?;
        if last.checked_add(k)? > seq_len {
            return None;
        }
        let word_idx = seq_start + first / 32;
        let within = first % 32;
        let first_word = *self.values.get(word_idx)?;
        let second = if within + span > 32 {
            *self.values.get(word_idx + 1)?
        } else {
            0
        };
        let third = if within + span > 64 {
            *self.values.get(word_idx + 2)?
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

    fn from_sequences(seqs: &[Vec<u8>], pool: Option<&ThreadPool>) -> Self {
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
    fn minimizer_at(&self, seq_idx: usize, start: usize, g: usize) -> Option<[u8; 8]> {
        if g == 0 || g > 32 {
            return None;
        }
        let word = self.encode_kmer_at(seq_idx, start, g)?;
        let shift = (32 - g) * 2;
        let forward = word << shift;
        let reversed = word.reverse_bits();
        let reverse = !(((reversed & 0x5555_5555_5555_5555) << 1)
            | ((reversed & 0xaaaa_aaaa_aaaa_aaaa) >> 1))
            & (u64::MAX << shift);
        Some(forward.min(reverse).to_le_bytes())
    }

    #[inline]
    pub(crate) fn len(&self, seq_idx: usize) -> Option<usize> {
        self.metadata.get(seq_idx).map(|(_, len)| len)
    }

    #[inline]
    pub(crate) fn encode_kmer_at(&self, seq_idx: usize, start: usize, k: usize) -> Option<u64> {
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
    g: usize,
    pub(crate) unitig_bases: EncodedBaseStore,
    pub(crate) km_bases: EncodedBaseStore,
}

impl EncodedUnitigs {
    #[inline]
    fn unitig_minimizer_count(&self, seq_idx: usize) -> usize {
        minimizer_count(
            self.unitig_bases.len(seq_idx).expect("valid unitig"),
            self.g,
        )
    }

    #[inline]
    fn unitig_minimizer_at(&self, seq_idx: usize, pos: usize) -> Option<[u8; 8]> {
        self.unitig_bases.minimizer_at(seq_idx, pos, self.g)
    }
}

#[derive(Clone, Copy)]
struct CompactEcBlock {
    lb: u32,
    ub: u32,
    ec_start: u32,
    ec_end: u32,
}

#[derive(Clone, Copy)]
struct WideEcBlock {
    lb: u32,
    ub: u32,
    ec_start: usize,
    ec_end: usize,
}

enum FlatEcBlocks {
    Compact(Vec<CompactEcBlock>),
    Wide(Vec<WideEcBlock>),
}

impl FlatEcBlocks {
    fn new() -> Self {
        Self::Compact(Vec::new())
    }

    fn len(&self) -> usize {
        match self {
            Self::Compact(v) => v.len(),
            Self::Wide(v) => v.len(),
        }
    }

    fn push(&mut self, lb: u32, ub: u32, ec_start: usize, ec_end: usize) {
        if let Self::Compact(values) = self {
            if ec_end <= u32::MAX as usize {
                values.push(CompactEcBlock {
                    lb,
                    ub,
                    ec_start: ec_start as u32,
                    ec_end: ec_end as u32,
                });
                return;
            }
            *self = Self::Wide(
                values
                    .iter()
                    .map(|b| WideEcBlock {
                        lb: b.lb,
                        ub: b.ub,
                        ec_start: b.ec_start as usize,
                        ec_end: b.ec_end as usize,
                    })
                    .collect(),
            );
        }
        if let Self::Wide(values) = self {
            values.push(WideEcBlock {
                lb,
                ub,
                ec_start,
                ec_end,
            });
        }
    }

    #[inline]
    fn get(&self, index: usize) -> Option<WideEcBlock> {
        match self {
            Self::Compact(values) => values.get(index).map(|b| WideEcBlock {
                lb: b.lb,
                ub: b.ub,
                ec_start: b.ec_start as usize,
                ec_end: b.ec_end as usize,
            }),
            Self::Wide(values) => values.get(index).copied(),
        }
    }
}

pub(crate) struct FlatEcIndex {
    unitig_block_offsets: MinzOffsets,
    blocks: FlatEcBlocks,
    ec_values: Vec<u32>,
}

impl FlatEcIndex {
    fn new(unitig_capacity: usize) -> Self {
        let mut offsets = Vec::with_capacity(unitig_capacity + 1);
        offsets.push(0);
        Self {
            unitig_block_offsets: MinzOffsets::Compact(offsets),
            blocks: FlatEcBlocks::new(),
            ec_values: Vec::new(),
        }
    }

    fn push_serialized_node(&mut self, node: &[u8]) -> Result<()> {
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
            self.blocks.push(lb, ub, start, self.ec_values.len());
            Ok(())
        })?;
        self.unitig_block_offsets.push(self.blocks.len());
        Ok(())
    }

    fn push_unitig_blocks(&mut self, blocks: &[crate::index::EcBlock]) {
        for block in blocks {
            let start = self.ec_values.len();
            self.ec_values.extend_from_slice(&block.ec);
            self.blocks
                .push(block.lb, block.ub, start, self.ec_values.len());
        }
        self.unitig_block_offsets.push(self.blocks.len());
    }

    pub(crate) fn block_index_for_position(&self, unitig_id: usize, pos: usize) -> Option<usize> {
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
        let global_idx = self
            .unitig_block_offsets
            .value_at(unitig_id)
            .expect("valid unitig")
            + block_idx;
        let block = self.blocks.get(global_idx).expect("valid EC block");
        &self.ec_values[block.ec_start..block.ec_end]
    }

    pub(crate) fn block_bounds(&self, unitig_id: usize, block_idx: usize) -> Option<(u32, u32)> {
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
    entries: HashMap<u64, Vec<(usize, u64)>>,
}

impl SelectiveFallbackKmerIndex {
    fn build(unitigs_len: usize, km_unitigs: &[Vec<u8>], h_kmers: &[Vec<u8>]) -> Self {
        let mut entries = HashMap::<u64, Vec<(usize, u64)>>::new();
        for (km_id, seq) in km_unitigs.iter().enumerate() {
            if let Some((fwd, rev)) = encode_kmer_pair(seq) {
                entries
                    .entry(fwd.min(rev))
                    .or_default()
                    .push((unitigs_len + km_id, fwd));
            }
        }
        let base_h = unitigs_len + km_unitigs.len();
        for (hk_id, seq) in h_kmers.iter().enumerate() {
            if let Some((fwd, rev)) = encode_kmer_pair(seq) {
                entries
                    .entry(fwd.min(rev))
                    .or_default()
                    .push((base_h + hk_id, fwd));
            }
        }
        Self { entries }
    }

    pub(crate) fn lookup(&self, code: u64) -> Option<&[(usize, u64)]> {
        self.entries.get(&code).map(Vec::as_slice)
    }
}

#[derive(Debug, PartialEq, Eq)]
enum MinzOffsets {
    Compact(Vec<u32>),
    Wide(Vec<usize>),
}

impl MinzOffsets {
    #[inline]
    fn value_at(&self, index: usize) -> Option<usize> {
        match self {
            Self::Compact(v) => v.get(index).map(|&x| x as usize),
            Self::Wide(v) => v.get(index).copied(),
        }
    }

    fn push(&mut self, value: usize) {
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

    fn from_counts(mut counts: Vec<usize>) -> (Self, usize) {
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

    fn len(&self) -> usize {
        match self {
            Self::Compact(v) => v.len(),
            Self::Wide(v) => v.len(),
        }
    }

    #[inline]
    fn bounds(&self, idx: usize) -> (usize, usize) {
        match self {
            Self::Compact(v) => (v[idx] as usize, v[idx + 1] as usize),
            Self::Wide(v) => (v[idx], v[idx + 1]),
        }
    }

    fn cursors(&self, buckets: usize) -> Vec<usize> {
        match self {
            Self::Compact(v) => v[..buckets].iter().map(|&p| p as usize).collect(),
            Self::Wide(v) => v[..buckets].to_vec(),
        }
    }
}

pub struct MinzPositionIndex {
    offsets: MinzOffsets,
    values: Vec<u64>,
    inline: Option<Vec<u64>>,
    inline_building: bool,
}

impl MinzPositionIndex {
    fn from_counts(counts: Vec<usize>) -> Self {
        let (offsets, total) = MinzOffsets::from_counts(counts);
        Self {
            offsets,
            values: vec![0; total],
            inline: None,
            inline_building: false,
        }
    }

    fn from_counts_inline(mut counts: Vec<usize>, escapes: &[usize], enabled: bool) -> Self {
        if !enabled || usize::BITS != 64 {
            return Self::from_counts(counts);
        }
        let total: usize = counts.iter().sum();
        if total > 0x3fff_ffff {
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

    fn len(&self) -> usize {
        self.inline
            .as_ref()
            .map_or_else(|| self.offsets.len().saturating_sub(1), Vec::len)
    }

    fn cursors(&self) -> Vec<usize> {
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
                let len = ((entry >> 32) & 0x3fff_ffff) as usize;
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

    fn push_at(&mut self, cursors: &mut [usize], idx: usize, pos_id: u64) {
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
                if (pos_id >> 32) as u32 == u32::MAX && (pos_id & 0x8000_0000) != 0 {
                    *entry |= 1u64 << 62;
                }
            }
        } else {
            let cursor = &mut cursors[idx];
            self.values[*cursor] = pos_id;
            *cursor += 1;
        }
    }

    fn finish_inline(&mut self) {
        if let Some(entries) = &mut self.inline {
            for entry in entries {
                if *entry >> 63 != 0 {
                    *entry -= (*entry >> 32) & 0x3fff_ffff;
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

static NEXT_CACHE_ID: AtomicU64 = AtomicU64::new(1);

pub struct BifrostIndex {
    pub(crate) cache_id: u64,
    pub index_version: u64,
    pub k: usize,
    pub g: usize,
    pub unitigs: Vec<Vec<u8>>,
    pub km_unitigs: Vec<Vec<u8>>,
    pub h_kmers: Vec<Vec<u8>>,
    pub h_kmer_map: Option<HashMap<u64, usize>>,
    pub dlist: Option<HashSet<u64>>,
    #[allow(dead_code)]
    pub(crate) encoded_unitigs: EncodedUnitigs,
    pub kmer_pos_index: Option<KmerPosIndex>,
    pub kmer_index: Option<KmerEcIndex>,
    pub(crate) ec_blocks: Vec<Vec<crate::index::EcBlock>>,
    pub(crate) flat_ec: FlatEcIndex,
    pub minz_positions: MinzPositionIndex,
    pub mphf: BooPhf,
    pub(crate) selective_fallback_kmers: SelectiveFallbackKmerIndex,
    pub onlist: Option<Vec<bool>>,
    pub transcript_lengths: Vec<u32>,
    pub transcript_names: Vec<String>,
    pub shade_to_color: Vec<Option<u32>>,
    pub shade_sequences: Vec<bool>,
    pub use_shade: bool,
}

impl BifrostIndex {
    pub fn has_strand_annotations(&self) -> bool {
        self.ec_blocks.iter().flatten().any(|block| {
            block
                .strands
                .as_ref()
                .is_some_and(|strands| !strands.is_empty())
        })
    }
}

pub fn build_bifrost_index(path: &Path) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, false, 1)
}

fn associate_abundant_node(
    kmers: &mut [Vec<u8>],
    map: &mut HashMap<u64, usize>,
    base: usize,
    node_row: usize,
    node_kmer: &[u8],
) -> Result<()> {
    let (forward, reverse) = encode_kmer_pair(node_kmer)
        .ok_or_else(|| Error::InvalidFormat("invalid abundant node k-mer".into()))?;
    let id = map
        .get_mut(&forward.min(reverse))
        .ok_or_else(|| Error::InvalidFormat("abundant node k-mer absent from graph".into()))?;
    let slot = node_row
        .checked_sub(base)
        .and_then(|offset| kmers.get_mut(offset))
        .ok_or_else(|| Error::InvalidFormat("abundant node out of range".into()))?;
    slot.copy_from_slice(node_kmer);
    *id = node_row;
    Ok(())
}

/// Load a flat counting index without duplicate EC blocks. Positional loading retains
/// the nested metadata needed for fragment/strand filters. Alternate k-mer scanners
/// and detailed hit tracing should use the full index loader instead.
pub fn build_bifrost_index_quant_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(
        path,
        load_positional_info,
        threads,
        load_positional_info,
    )
}

pub fn build_bifrost_index_with_kmer_pos(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index_with_positions_threaded(path, load_positional_info, 1)?;
    index.kmer_pos_index = Some(build_kmer_pos_index(&index));
    Ok(index)
}

pub fn build_bifrost_index_with_positions(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, load_positional_info, 1)
}

pub fn build_bifrost_index_threaded(path: &Path, threads: usize) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, false, threads)
}

pub fn build_bifrost_index_with_kmer_pos_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    let mut index =
        build_bifrost_index_with_positions_threaded(path, load_positional_info, threads)?;
    index.kmer_pos_index = Some(build_kmer_pos_index(&index));
    Ok(index)
}

pub fn build_bifrost_index_with_positions_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(path, load_positional_info, threads, true)
}

fn build_bifrost_index_with_positions_threaded_impl(
    path: &Path,
    load_positional_info: bool,
    _threads: usize,
    store_ec_blocks: bool,
) -> Result<BifrostIndex> {
    let _header_timing = timing::scoped(Stage::IndexHeaderParse);
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);

    let index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;
    drop(_header_timing);

    let _graph_timing = timing::scoped(Stage::GraphDecode);
    let file_format_version = read_u64_le(&mut reader)?;
    let k = read_i32_le(&mut reader)?;
    let g = read_i32_le(&mut reader)?;
    let _ = file_format_version;

    if k <= 0 || g <= 0 {
        return Err(Error::InvalidFormat("invalid k/g".into()));
    }
    let k = k as usize;
    let g = g as usize;
    if k > 32 || g > 32 {
        return Err(Error::UnsupportedFeature(
            "k/g > 32 not supported in Bifrost path".into(),
        ));
    }

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut unitigs = Vec::with_capacity(v_unitigs_sz);
    let mut compressed_scratch = Vec::new();
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        compressed_scratch.resize(data_sz, 0);
        reader.read_exact(&mut compressed_scratch)?;
        unitigs.push(decode_compressed_sequence(&compressed_scratch, len));
    }

    drop(compressed_scratch);
    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut km_unitigs = Vec::with_capacity(km_unitigs_sz);
    for _ in 0..km_unitigs_sz {
        let mut buf = [0u8; 8];
        reader.read_exact(&mut buf)?;
        km_unitigs.push(decode_kmer(&buf, k));
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    let mut h_kmers = Vec::with_capacity(h_kmers_ccov_sz);
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = [0u8; 8];
        reader.read_exact(&mut buf)?;
        h_kmers.push(decode_kmer(&buf, k));
    }
    let mut h_kmer_map = if h_kmers.is_empty() {
        None
    } else {
        let mut map = HashMap::with_capacity(h_kmers.len());
        let base = unitigs.len() + km_unitigs.len();
        for (idx, seq) in h_kmers.iter().enumerate() {
            if let Some((fwd, rev)) = encode_kmer_pair(seq) {
                let code = fwd.min(rev);
                map.insert(code, base + idx);
            }
        }
        Some(map)
    };

    let _bfg_header = read_u64_le(&mut reader)?;
    let _checksum = read_u64_le(&mut reader)?;
    let _v_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _km_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _h_kmers_sz2 = read_u64_le(&mut reader)?;
    let hmap_min_unitigs_sz = read_u64_le(&mut reader)? as usize;

    let nb_bmp_unitigs = read_u64_le(&mut reader)? as usize;
    let mut bmp_unitigs = Vec::with_capacity(nb_bmp_unitigs);
    for _ in 0..nb_bmp_unitigs {
        bmp_unitigs.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_bmp_km = read_u64_le(&mut reader)? as usize;
    let mut bmp_km = Vec::with_capacity(nb_bmp_km);
    for _ in 0..nb_bmp_km {
        bmp_km.push(read_bitcontainer_values(&mut reader)?);
    }

    let nb_special_minz = read_u64_le(&mut reader)? as usize;
    let nb_bmp_special = (nb_special_minz >> 32) + 1;
    let mut abundant = HashSet::new();
    let mut overcrowded = HashSet::new();
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            abundant.insert(id as usize);
        }
    }
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            overcrowded.insert(id as usize);
        }
    }

    let mut special_minz = Vec::with_capacity(nb_special_minz);
    for i in 0..nb_special_minz {
        let minz = read_u64_le(&mut reader)?;
        let is_abundant = abundant.contains(&i);
        let is_overcrowded = overcrowded.contains(&i);
        let pos_id = if is_abundant || is_overcrowded {
            let mut id = 0xffff_ffff_0000_0000 | if is_overcrowded { 0x8000_0000 } else { 0 };
            if is_abundant {
                let count = read_u32_le(&mut reader)? as u64;
                id |= count;
            }
            id
        } else {
            read_u64_le(&mut reader)?
        };
        special_minz.push((minz, pos_id));
    }

    drop(abundant);
    drop(overcrowded);
    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;

    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mphf = BooPhf::load(&mut mphf_cur)?;
    drop(mphf_buf);

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let mut dlist = if dlist_size > 0 {
        HashSet::with_capacity(dlist_size as usize)
    } else {
        HashSet::new()
    };
    for _ in 0..dlist_size {
        let mut buf = [0u8; 8];
        reader.read_exact(&mut buf)?;
        let seq = decode_kmer(&buf, k);
        if let Some((fwd, rev)) = encode_kmer_pair(&seq) {
            dlist.insert(fwd.min(rev));
        }
    }
    let dlist = if dlist.is_empty() { None } else { Some(dlist) };

    let node_count = read_u64_le(&mut reader)? as usize;
    let mut ec_blocks = if store_ec_blocks {
        Vec::with_capacity(node_count)
    } else {
        Vec::new()
    };
    let mut flat_ec = FlatEcIndex::new(node_count);
    let mut node_buf = Vec::new();
    let abundant_base = unitigs.len() + km_unitigs.len();
    if node_count != abundant_base + h_kmers.len() {
        return Err(Error::InvalidFormat(
            "unitig and node counts mismatch".into(),
        ));
    }
    let mut node_kmer = vec![0; k];
    for node_row in 0..node_count {
        if node_row < abundant_base {
            // Vector-backed graph entries retain their serialized order.
            reader.seek_relative(k as i64)?;
        } else {
            // Bifrost reinserts abundant k-mers into a hash table when loading
            // the graph. Its iteration order can differ from the graph file's
            // order. Node records carry the authoritative k-mer/EC association.
            reader.read_exact(&mut node_kmer)?;
            associate_abundant_node(
                &mut h_kmers,
                h_kmer_map.as_mut().expect("nonempty abundant k-mers"),
                abundant_base,
                node_row,
                &node_kmer,
            )?;
        }
        let node_size = read_u32_le(&mut reader)? as usize;
        node_buf.resize(node_size, 0);
        reader.read_exact(&mut node_buf)?;
        if !store_ec_blocks && !load_positional_info {
            flat_ec.push_serialized_node(&node_buf)?;
            continue;
        }
        let mut cur = std::io::Cursor::new(node_buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let blocks = if load_positional_info {
            crate::index::read_block_array_blocks_with_positions(&mut cur)?
        } else {
            read_block_array_blocks(&mut cur)?
        };
        flat_ec.push_unitig_blocks(&blocks);
        if store_ec_blocks {
            ec_blocks.push(blocks);
        }
    }
    drop(node_buf);
    let selective_fallback_kmers =
        SelectiveFallbackKmerIndex::build(unitigs.len(), &km_unitigs, &h_kmers);
    let load_threads = _threads.max(1);
    let rayon_pool = if load_threads > 1 {
        Some(
            ThreadPoolBuilder::new()
                .num_threads(load_threads)
                .build()
                .map_err(|e| Error::InvalidFormat(format!("failed to build rayon pool: {e}")))?,
        )
    } else {
        None
    };
    let encoded_unitigs = precompute_encoded_unitigs(&unitigs, &km_unitigs, g, rayon_pool.as_ref());

    let mut minz_counts = Vec::with_capacity(hmap_min_unitigs_sz + 1);
    minz_counts.resize(hmap_min_unitigs_sz, 0);
    drop(_graph_timing);
    let _count_timing = timing::scoped(Stage::MinimizerCountPass);
    let mut unitig_bucket_cache = (bmp_unitigs.len() == 1
        && hmap_min_unitigs_sz < u32::MAX as usize
        && !store_ec_blocks
        && !load_positional_info)
        .then(|| Vec::with_capacity(bmp_unitigs[0].len()));
    populate_unitig_minimizer_counts(
        &mut minz_counts,
        &bmp_unitigs,
        &encoded_unitigs,
        &mphf,
        rayon_pool.as_ref(),
        unitig_bucket_cache.as_mut(),
    );

    let km_glen = k.saturating_sub(g) + 1;
    let kmer_load_ctx = KmerLoadContext {
        km_bases: &encoded_unitigs.km_bases,
        g,
        km_glen,
        mphf: &mphf,
    };
    populate_kmer_minimizer_counts(
        &mut minz_counts,
        &bmp_km,
        &kmer_load_ctx,
        rayon_pool.as_ref(),
    );

    let mut inline_escapes = Vec::new();
    for (minz, pos_id) in special_minz.iter().copied() {
        let bytes = minz.to_le_bytes();
        if let Some(idx) = mphf.lookup(&bytes) {
            minz_counts[idx as usize] += 1;
            if pos_id >> 63 != 0 {
                inline_escapes.push(idx as usize);
            }
        }
    }
    inline_escapes.sort_unstable();
    inline_escapes.dedup();
    let inline_enabled = !store_ec_blocks
        && !load_positional_info
        && bmp_unitigs.len() <= 1
        && bmp_km.len() <= 1
        && unitigs.len() < (1usize << 31)
        && km_unitigs.len() < (1usize << 31)
        && hmap_min_unitigs_sz < u32::MAX as usize
        && (bmp_unitigs.is_empty() || unitig_bucket_cache.is_some() || rayon_pool.is_none());
    let mut minz_positions =
        MinzPositionIndex::from_counts_inline(minz_counts, &inline_escapes, inline_enabled);
    let mut minz_cursors = minz_positions.cursors();
    drop(_count_timing);
    let _fill_timing = timing::scoped(Stage::MinimizerFillPass);
    populate_unitig_minimizer_storage(
        &mut minz_positions,
        &mut minz_cursors,
        &bmp_unitigs,
        &encoded_unitigs,
        &mphf,
        rayon_pool.as_ref(),
        unitig_bucket_cache.as_deref(),
    );
    drop(unitig_bucket_cache);
    drop(bmp_unitigs);
    populate_kmer_minimizer_storage(
        &mut minz_positions,
        &mut minz_cursors,
        &bmp_km,
        &kmer_load_ctx,
        rayon_pool.as_ref(),
    );

    drop(bmp_km);
    for (minz, pos_id) in special_minz {
        let bytes = minz.to_le_bytes();
        if let Some(idx) = mphf.lookup(&bytes) {
            minz_positions.push_at(&mut minz_cursors, idx as usize, pos_id);
        }
    }

    drop(minz_cursors);
    minz_positions.finish_inline();
    let onlist = read_onlist_with_lengths(&mut reader)?;

    let cache_id = loop {
        let id = NEXT_CACHE_ID.load(Ordering::Relaxed);
        let next = id
            .checked_add(1)
            .expect("index cache identity space exhausted");
        if NEXT_CACHE_ID
            .compare_exchange_weak(id, next, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
        {
            break id;
        }
    };
    Ok(BifrostIndex {
        cache_id,
        index_version,
        k,
        g,
        unitigs,
        km_unitigs,
        h_kmers,
        h_kmer_map,
        dlist,
        encoded_unitigs,
        kmer_pos_index: None,
        kmer_index: None,
        ec_blocks,
        flat_ec,
        minz_positions,
        mphf,
        selective_fallback_kmers,
        onlist: onlist.onlist,
        transcript_lengths: onlist.lengths,
        transcript_names: onlist.names,
        shade_to_color: onlist.shade_to_color,
        shade_sequences: onlist.shade_sequences,
        use_shade: onlist.use_shade,
    })
}

struct KmerLoadContext<'a> {
    km_bases: &'a EncodedBaseStore,
    g: usize,
    km_glen: usize,
    mphf: &'a BooPhf,
}

#[cfg(test)]
fn base_code(base: u8) -> Option<u64> {
    let code = BASE_CODES[base as usize];
    if code == INVALID_BASE_CODE {
        None
    } else {
        Some(u64::from(code))
    }
}

#[cfg(test)]
fn encode_word_pair(seq: &[u8]) -> Option<(u64, u64)> {
    if seq.is_empty() {
        return Some((0, 0));
    }
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let last = seq.len() - 1;
    for (i, &base) in seq.iter().enumerate() {
        let code = base_code(base)?;
        fwd |= code << (62 - ((i & 0x1f) << 1));
        rev |= (3 - code) << (62 - (((last - i) & 0x1f) << 1));
    }
    Some((fwd, rev))
}

#[inline]
fn minimizer_count(seq_len: usize, g: usize) -> usize {
    if g == 0 || seq_len < g {
        0
    } else {
        seq_len - g + 1
    }
}

#[cfg(test)]
fn compute_minimizer_reps_into(seq: &[u8], g: usize, out: &mut [[u8; 8]]) {
    if out.is_empty() {
        return;
    }

    debug_assert_eq!(out.len(), minimizer_count(seq.len(), g));

    let Some((mut fwd, mut rev)) = encode_word_pair(&seq[..g]) else {
        return;
    };
    let tail_shift = 64usize.saturating_sub(g * 2);
    let high_mask = if g == 32 {
        u64::MAX
    } else {
        (!0u64) << tail_shift
    };

    out[0] = fwd.min(rev).to_le_bytes();

    for (dst, &base) in out.iter_mut().skip(1).zip(&seq[g..]) {
        let code = BASE_CODES[base as usize];
        if code == INVALID_BASE_CODE {
            return;
        }
        let code = u64::from(code);
        let comp = 3 - code;
        fwd = ((fwd << 2) & high_mask) | (code << tail_shift);
        rev = ((rev >> 2) & high_mask) | (comp << 62);
        *dst = fwd.min(rev).to_le_bytes();
    }
}

fn pack_bases_into(seq: &[u8], out: &mut [u64]) {
    debug_assert_eq!(seq.len().div_ceil(32), out.len());
    for (chunk, word) in seq.chunks(32).zip(out) {
        let mut encoded = 0u64;
        for &base in chunk {
            encoded = (encoded << 2) | u64::from(BASE_CODES[base as usize] & 3);
        }
        *word = encoded << ((32 - chunk.len()) * 2);
    }
}

#[cfg(test)]
mod packed_base_tests {
    use super::SequenceMetadata;
    use super::{BASE_CODES, EncodedBaseStore, compute_minimizer_reps_into, minimizer_count};

    #[test]
    fn rolling_packed_windows_match_independent_windows() {
        let seqs = vec![
            (0..193)
                .map(|i| b"ACGT"[(i * i + i / 7) % 4])
                .collect::<Vec<_>>(),
        ];
        let store = EncodedBaseStore::from_sequences(&seqs, None);
        for k in 1..=32 {
            for first in 0..=193 - k {
                let last = (first + 64 - k).min(193 - k);
                let got: Vec<_> = store.kmer_range(0, first, last, k).unwrap().collect();
                let expected: Vec<_> = (first..=last)
                    .map(|pos| (pos, store.encode_kmer_at(0, pos, k).unwrap()))
                    .collect();
                assert_eq!(got, expected, "k={k} first={first}");
            }
        }
        assert!(store.kmer_range(0, 0, 64, 1).is_none());
        assert!(store.kmer_range(0, 4, 3, 31).is_none());
        assert!(store.kmer_range(0, 180, 190, 31).is_none());
        assert!(store.kmer_range(0, usize::MAX, usize::MAX, 31).is_none());
    }

    #[test]
    fn packed_metadata_preserves_wide_boundaries() {
        let compact = SequenceMetadata::new(vec![0, 2], vec![63]);
        assert!(matches!(compact, SequenceMetadata::Compact(_)));
        assert_eq!(compact.get(0), Some((0, 63)));
        assert_eq!(compact.get(1), None);
        let Some(big) = (u32::MAX as usize).checked_add(1) else {
            return;
        };
        let wide = SequenceMetadata::new(vec![0, big, big + 1], vec![big, 1]);
        assert!(matches!(wide, SequenceMetadata::Wide(_)));
        assert_eq!(wide.get(1), Some((big, 1)));
    }

    #[test]
    fn packed_minimizers_equal_byte_reference_for_all_g_and_word_boundaries() {
        let mut state = 43u64;
        let sequences: Vec<Vec<u8>> = (0..=129)
            .map(|len| {
                (0..len)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        b"ACGTacgt"[(state & 7) as usize]
                    })
                    .collect()
            })
            .chain([vec![b'A'; 97], vec![b'T'; 97], b"ACGT".repeat(32)])
            .collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        for store in [
            EncodedBaseStore::from_sequences(&sequences, None),
            EncodedBaseStore::from_sequences(&sequences, Some(&pool)),
        ] {
            for (seq_id, sequence) in sequences.iter().enumerate() {
                for g in 1..=32 {
                    let mut expected = vec![[0; 8]; minimizer_count(sequence.len(), g)];
                    compute_minimizer_reps_into(sequence, g, &mut expected);
                    for (pos, bytes) in expected.iter().enumerate() {
                        assert_eq!(
                            store.minimizer_at(seq_id, pos, g).as_ref(),
                            Some(bytes),
                            "sequence {seq_id}, position {pos}, g {g}"
                        );
                    }
                    assert_eq!(store.minimizer_at(seq_id, expected.len(), g), None);
                    assert_eq!(store.minimizer_at(seq_id, usize::MAX, g), None);
                }
                assert_eq!(store.minimizer_at(seq_id, 0, 0), None);
                assert_eq!(store.minimizer_at(seq_id, 0, 33), None);
            }
            assert_eq!(store.minimizer_at(sequences.len(), 0, 1), None);
        }
    }

    #[test]
    fn packed_windows_equal_original_byte_encoding_at_every_boundary() {
        let sequences = vec![
            Vec::new(),
            b"A".to_vec(),
            b"ACGTNacgtnX".to_vec(),
            (0..129).map(|i| b"ACGTacgtN"[i % 9]).collect(),
        ];
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .expect("thread pool");
        for store in [
            EncodedBaseStore::from_sequences(&sequences, None),
            EncodedBaseStore::from_sequences(&sequences, Some(&pool)),
        ] {
            for (seq_id, sequence) in sequences.iter().enumerate() {
                assert_eq!(store.len(seq_id), Some(sequence.len()));
                for start in 0..=sequence.len() + 1 {
                    for k in 0..=33 {
                        let expected = if k <= 32 && start + k <= sequence.len() {
                            Some(sequence[start..start + k].iter().fold(0u64, |word, &base| {
                                (word << 2) | u64::from(BASE_CODES[base as usize] & 3)
                            }))
                        } else {
                            None
                        };
                        assert_eq!(store.encode_kmer_at(seq_id, start, k), expected);
                    }
                }
                assert_eq!(store.encode_kmer_at(seq_id, usize::MAX, 1), None);
            }
            assert_eq!(store.len(sequences.len()), None);
            assert_eq!(store.encode_kmer_at(sequences.len(), 0, 0), None);
        }
    }
}

fn precompute_encoded_unitigs(
    unitigs: &[Vec<u8>],
    km_unitigs: &[Vec<u8>],
    g: usize,
    pool: Option<&ThreadPool>,
) -> EncodedUnitigs {
    EncodedUnitigs {
        g,
        unitig_bases: EncodedBaseStore::from_sequences(unitigs, pool),
        km_bases: EncodedBaseStore::from_sequences(km_unitigs, pool),
    }
}

#[derive(Clone, Copy)]
struct UnitigBlockState {
    unitig_id: usize,
    total_unitig_len: u64,
    current_unitig_len: u64,
}

fn advance_unitig_state(state: &mut UnitigBlockState, encoded_unitigs: &EncodedUnitigs, pos: u64) {
    while pos >= state.total_unitig_len + state.current_unitig_len {
        state.unitig_id += 1;
        state.total_unitig_len += state.current_unitig_len;
        state.current_unitig_len = encoded_unitigs.unitig_minimizer_count(state.unitig_id) as u64;
    }
}

fn precompute_unitig_block_states(
    bmp_unitigs: &[Vec<u32>],
    encoded_unitigs: &EncodedUnitigs,
) -> Vec<UnitigBlockState> {
    let mut out = Vec::with_capacity(bmp_unitigs.len());
    if encoded_unitigs.unitig_bases.metadata.is_empty() {
        return out;
    }

    let mut state = UnitigBlockState {
        unitig_id: 0,
        total_unitig_len: 0,
        current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
    };
    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let base = (i as u64) << 32;
        advance_unitig_state(&mut state, encoded_unitigs, base);
        out.push(state);
        if let Some(&last_pos) = bitmap.last() {
            advance_unitig_state(&mut state, encoded_unitigs, base + last_pos as u64);
        }
    }
    out
}

fn populate_unitig_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_unitigs: &[Vec<u32>],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
    mut bucket_cache: Option<&mut Vec<u32>>,
) {
    if bmp_unitigs.is_empty() || encoded_unitigs.unitig_bases.metadata.is_empty() {
        return;
    }

    let Some(pool) = pool else {
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        for (i, bitmap) in bmp_unitigs.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                let bucket = encoded_unitigs
                    .unitig_minimizer_at(state.unitig_id, rel)
                    .and_then(|bytes| mphf.lookup(&bytes));
                if let Some(cache) = bucket_cache.as_mut() {
                    cache.push(bucket.map_or(u32::MAX, |bucket| bucket as u32));
                }
                if let Some(bucket) = bucket {
                    minz_counts[bucket as usize] += 1;
                }
            }
        }
        return;
    };

    let states = precompute_unitig_block_states(bmp_unitigs, encoded_unitigs);
    let mut batch_pos_ids = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut flush_batch = |batch_pos_ids: &mut Vec<u64>| {
        if batch_pos_ids.is_empty() {
            return;
        }
        if let Some(cache) = bucket_cache.as_mut() {
            let chunk_len = (batch_pos_ids.len() / (pool.current_num_threads() * 4).max(1))
                .max(LOAD_MIN_PAR_CHUNK);
            let old_len = cache.len();
            cache.reserve(batch_pos_ids.len());
            let destinations = &mut cache.spare_capacity_mut()[..batch_pos_ids.len()];
            let chunks: Vec<Vec<u32>> = pool.install(|| {
                destinations
                    .par_chunks_mut(chunk_len)
                    .zip(batch_pos_ids.par_chunks(chunk_len))
                    .map(|(destinations, positions)| {
                        let mut indices = Vec::with_capacity(positions.len());
                        for (destination, &pos_id) in destinations.iter_mut().zip(positions) {
                            let bucket = encoded_unitigs
                                .unitig_minimizer_at(
                                    (pos_id >> 32) as usize,
                                    pos_id as u32 as usize,
                                )
                                .and_then(|bytes| mphf.lookup(&bytes))
                                .map_or(u32::MAX, |bucket| bucket as u32);
                            destination.write(bucket);
                            if bucket != u32::MAX {
                                indices.push(bucket);
                            }
                        }
                        indices.sort_unstable();
                        indices
                    })
                    .collect()
            });
            // SAFETY: the indexed parallel zip writes every destination exactly once.
            // A panic unwinds before this publication and leaves the old length intact.
            unsafe {
                cache.set_len(old_len + batch_pos_ids.len());
            }
            for chunk in chunks {
                let mut at = 0;
                while at < chunk.len() && chunk[at] != u32::MAX {
                    let bucket = chunk[at];
                    let mut end = at + 1;
                    while end < chunk.len() && chunk[end] == bucket {
                        end += 1;
                    }
                    minz_counts[bucket as usize] += end - at;
                    at = end;
                }
            }
        } else {
            let grouped = collect_unitig_batch_counts(batch_pos_ids, encoded_unitigs, mphf, pool);
            merge_count_batches(minz_counts, grouped);
        }
        batch_pos_ids.clear();
    };

    for (i, bitmap) in bmp_unitigs.iter().enumerate() {
        let mut state = states[i];
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            let rel = (pos - state.total_unitig_len) as usize;
            batch_pos_ids.push(((state.unitig_id as u64) << 32) | rel as u64);
            if batch_pos_ids.len() >= LOAD_BATCH_POSITIONS {
                flush_batch(&mut batch_pos_ids);
            }
        }
    }
    flush_batch(&mut batch_pos_ids);
}

fn populate_unitig_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_unitigs: &[Vec<u32>],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: Option<&ThreadPool>,
    bucket_cache: Option<&[u32]>,
) {
    if bmp_unitigs.is_empty() || encoded_unitigs.unitig_bases.metadata.is_empty() {
        return;
    }

    if let Some(cache) = bucket_cache {
        debug_assert_eq!(bmp_unitigs.len(), 1);
        debug_assert_eq!(cache.len(), bmp_unitigs[0].len());
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        for (&pos, &bucket) in bmp_unitigs[0].iter().zip(cache) {
            let pos = u64::from(pos);
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            if bucket != u32::MAX {
                let rel = pos - state.total_unitig_len;
                minz_positions.push_at(
                    minz_cursors,
                    bucket as usize,
                    ((state.unitig_id as u64) << 32) | rel,
                );
            }
        }
        return;
    }
    let Some(pool) = pool else {
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        for (i, bitmap) in bmp_unitigs.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizer_at(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(&bytes)
                {
                    let pos_id = ((state.unitig_id as u64) << 32) | rel as u64;
                    minz_positions.push_at(minz_cursors, idx as usize, pos_id);
                }
            }
        }
        return;
    };

    if bmp_unitigs.len() == 1 {
        // A single bitmap otherwise gives the outer parallel loop only one task.
        // Resolve MPHF buckets concurrently, then append in the original order.
        // Bucket order affects which exact matching candidate is visited first.
        let mut state = UnitigBlockState {
            unitig_id: 0,
            total_unitig_len: 0,
            current_unitig_len: encoded_unitigs.unitig_minimizer_count(0) as u64,
        };
        let mut batch = Vec::with_capacity(LOAD_BATCH_POSITIONS);
        let mut flush = |batch: &mut Vec<u64>| {
            let chunk_len =
                (batch.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
            let events: Vec<Vec<(usize, u64)>> = pool.install(|| {
                batch
                    .par_chunks(chunk_len)
                    .map(|chunk| {
                        chunk
                            .iter()
                            .filter_map(|&pos_id| {
                                let uid = (pos_id >> 32) as usize;
                                let rel = (pos_id & 0xffff_ffff) as usize;
                                let bytes = encoded_unitigs.unitig_minimizer_at(uid, rel)?;
                                Some((mphf.lookup(&bytes)? as usize, pos_id))
                            })
                            .collect()
                    })
                    .collect()
            });
            for chunk in events {
                for (bucket, pos_id) in chunk {
                    minz_positions.push_at(minz_cursors, bucket, pos_id);
                }
            }
            batch.clear();
        };
        for &pos in &bmp_unitigs[0] {
            let pos = u64::from(pos);
            advance_unitig_state(&mut state, encoded_unitigs, pos);
            let rel = pos - state.total_unitig_len;
            batch.push(((state.unitig_id as u64) << 32) | rel);
            if batch.len() >= LOAD_BATCH_POSITIONS {
                flush(&mut batch);
            }
        }
        flush(&mut batch);
        return;
    }

    let states = precompute_unitig_block_states(bmp_unitigs, encoded_unitigs);
    let cursors = minz_cursors
        .iter()
        .copied()
        .map(AtomicUsize::new)
        .collect::<Vec<_>>();
    let values_ptr = minz_positions.values.as_mut_ptr() as usize;
    pool.install(|| {
        bmp_unitigs.par_iter().enumerate().for_each(|(i, bitmap)| {
            let mut state = states[i];
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                advance_unitig_state(&mut state, encoded_unitigs, pos);
                let rel = (pos - state.total_unitig_len) as usize;
                if let Some(bytes) = encoded_unitigs.unitig_minimizer_at(state.unitig_id, rel)
                    && let Some(idx) = mphf.lookup(&bytes)
                {
                    let pos_id = ((state.unitig_id as u64) << 32) | rel as u64;
                    let slot = cursors[idx as usize].fetch_add(1, Ordering::Relaxed);
                    // SAFETY: the precomputed bucket counts make each fetched slot unique,
                    // so parallel workers never write the same position.
                    unsafe {
                        (values_ptr as *mut u64).add(slot).write(pos_id);
                    }
                }
            }
        });
    });
    for (dst, src) in minz_cursors.iter_mut().zip(cursors) {
        *dst = src.into_inner();
    }
}

fn populate_kmer_minimizer_counts(
    minz_counts: &mut [usize],
    bmp_km: &[Vec<u32>],
    ctx: &KmerLoadContext<'_>,
    pool: Option<&ThreadPool>,
) {
    let Some(pool) = pool else {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                    && let Some(idx) = ctx.mphf.lookup(&bytes)
                {
                    minz_counts[idx as usize] += 1;
                }
            }
        }
        return;
    };

    let mut batch_pos_ids = Vec::with_capacity(LOAD_BATCH_POSITIONS);
    let mut flush_batch = |batch_pos_ids: &mut Vec<u64>| {
        if batch_pos_ids.is_empty() {
            return;
        }
        let grouped = collect_kmer_batch_counts(batch_pos_ids, ctx, pool);
        merge_count_batches(minz_counts, grouped);
        batch_pos_ids.clear();
    };

    for (i, bitmap) in bmp_km.iter().enumerate() {
        let base = (i as u64) << 32;
        for &pos_bmp in bitmap {
            let pos = base + pos_bmp as u64;
            let km_id = (pos / ctx.km_glen as u64) as usize;
            let km_pos = (pos % ctx.km_glen as u64) as usize;
            batch_pos_ids.push(((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64);
            if batch_pos_ids.len() >= LOAD_BATCH_POSITIONS {
                flush_batch(&mut batch_pos_ids);
            }
        }
    }
    flush_batch(&mut batch_pos_ids);
}

fn populate_kmer_minimizer_storage(
    minz_positions: &mut MinzPositionIndex,
    minz_cursors: &mut [usize],
    bmp_km: &[Vec<u32>],
    ctx: &KmerLoadContext<'_>,
    pool: Option<&ThreadPool>,
) {
    // A single bitmap supplies only one outer Rayon task, so atomics do not
    // provide parallelism and duplicate the entire reference-sized cursor table.
    let Some(pool) = pool.filter(|_| bmp_km.len() > 1) else {
        for (i, bitmap) in bmp_km.iter().enumerate() {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                    && let Some(idx) = ctx.mphf.lookup(&bytes)
                {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64;
                    minz_positions.push_at(minz_cursors, idx as usize, pos_id);
                }
            }
        }
        return;
    };

    let cursors = minz_cursors
        .iter()
        .copied()
        .map(AtomicUsize::new)
        .collect::<Vec<_>>();
    let values_ptr = minz_positions.values.as_mut_ptr() as usize;
    pool.install(|| {
        bmp_km.par_iter().enumerate().for_each(|(i, bitmap)| {
            let base = (i as u64) << 32;
            for &pos_bmp in bitmap {
                let pos = base + pos_bmp as u64;
                let km_id = (pos / ctx.km_glen as u64) as usize;
                let km_pos = (pos % ctx.km_glen as u64) as usize;
                if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                    && let Some(idx) = ctx.mphf.lookup(&bytes)
                {
                    let pos_id = ((km_id as u64) << 32) | 0x8000_0000 | km_pos as u64;
                    let slot = cursors[idx as usize].fetch_add(1, Ordering::Relaxed);
                    // SAFETY: the precomputed bucket counts make each fetched slot unique,
                    // so parallel workers never write the same position.
                    unsafe {
                        (values_ptr as *mut u64).add(slot).write(pos_id);
                    }
                }
            }
        });
    });
    for (dst, src) in minz_cursors.iter_mut().zip(cursors) {
        *dst = src.into_inner();
    }
}

fn collect_unitig_batch_counts(
    pos_ids: &[u64],
    encoded_unitigs: &EncodedUnitigs,
    mphf: &BooPhf,
    pool: &ThreadPool,
) -> Vec<Vec<usize>> {
    let chunk_len =
        (pos_ids.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    pool.install(|| {
        pos_ids
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut indices = Vec::with_capacity(chunk.len());
                for &pos_id in chunk {
                    let unitig_id = (pos_id >> 32) as usize;
                    let rel = (pos_id & 0xffff_ffff) as usize;
                    if let Some(bytes) = encoded_unitigs.unitig_minimizer_at(unitig_id, rel)
                        && let Some(idx) = mphf.lookup(&bytes)
                    {
                        indices.push(idx as usize);
                    }
                }
                indices.sort_unstable();
                indices
            })
            .collect()
    })
}

fn collect_kmer_batch_counts(
    pos_ids: &[u64],
    ctx: &KmerLoadContext<'_>,
    pool: &ThreadPool,
) -> Vec<Vec<usize>> {
    let chunk_len =
        (pos_ids.len() / (pool.current_num_threads() * 4).max(1)).max(LOAD_MIN_PAR_CHUNK);
    pool.install(|| {
        pos_ids
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut indices = Vec::with_capacity(chunk.len());
                for &pos_id in chunk {
                    let km_id = (pos_id >> 32) as usize;
                    let km_pos = (pos_id & 0x7fff_ffff) as usize;
                    if let Some(bytes) = ctx.km_bases.minimizer_at(km_id, km_pos, ctx.g)
                        && let Some(idx) = ctx.mphf.lookup(&bytes)
                    {
                        indices.push(idx as usize);
                    }
                }
                indices.sort_unstable();
                indices
            })
            .collect()
    })
}

fn merge_count_batches(minz_counts: &mut [usize], grouped: Vec<Vec<usize>>) {
    for indices in grouped {
        let mut iter = indices.into_iter();
        let Some(mut idx) = iter.next() else {
            continue;
        };
        let mut count = 1usize;
        for next_idx in iter {
            if next_idx == idx {
                count += 1;
            } else {
                minz_counts[idx] += count;
                idx = next_idx;
                count = 1;
            }
        }
        minz_counts[idx] += count;
    }
}

fn build_kmer_pos_index(index: &BifrostIndex) -> KmerPosIndex {
    let mut map: KmerPosIndex = HashMap::new();
    let k = index.k;
    for (unitig_id, seq) in index.unitigs.iter().enumerate() {
        if seq.len() < k {
            continue;
        }
        for pos in 0..=seq.len() - k {
            if let Some((fwd, rev)) = encode_kmer_pair(&seq[pos..pos + k]) {
                map.entry(fwd).or_default().push((unitig_id, pos, false));
                if rev != fwd {
                    map.entry(rev).or_default().push((unitig_id, pos, true));
                }
            }
        }
    }
    let base_km = index.unitigs.len();
    for (km_id, seq) in index.km_unitigs.iter().enumerate() {
        if seq.len() != k {
            continue;
        }
        if let Some((fwd, rev)) = encode_kmer_pair(seq) {
            let unitig_id = base_km + km_id;
            map.entry(fwd).or_default().push((unitig_id, 0, false));
            if rev != fwd {
                map.entry(rev).or_default().push((unitig_id, 0, true));
            }
        }
    }
    let base_h = index.unitigs.len() + index.km_unitigs.len();
    for (hk_id, seq) in index.h_kmers.iter().enumerate() {
        if seq.len() != k {
            continue;
        }
        if let Some((fwd, rev)) = encode_kmer_pair(seq) {
            let unitig_id = base_h + hk_id;
            map.entry(fwd).or_default().push((unitig_id, 0, false));
            if rev != fwd {
                map.entry(rev).or_default().push((unitig_id, 0, true));
            }
        }
    }
    map
}

pub(super) fn read_unitig_sequences<R: Read + Seek>(
    reader: &mut R,
    dbg_start: u64,
    dbg_size: u64,
    k: u32,
    kmer_bytes: usize,
) -> Result<(Vec<Vec<u8>>, usize)> {
    reader.seek(SeekFrom::Start(dbg_start))?;
    let _file_format_version = read_u64_le(reader)?;
    let _k = read_i32_le(reader)?;
    let _g = read_i32_le(reader)?;

    let v_unitigs_sz = read_u64_le(reader)?;
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(v_unitigs_sz as usize);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        seqs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(reader)?;
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    let abundant_base = seqs.len();
    let h_kmers_ccov_sz = read_u64_le(reader)?;
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    Ok((seqs, abundant_base))
}

fn decode_compressed_sequence(data: &[u8], len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    let full = len / 4;
    for &byte in &data[..full] {
        out.extend_from_slice(&DECODED_BYTE_BASES[byte as usize]);
    }
    let tail = len % 4;
    if tail != 0 {
        out.extend_from_slice(&DECODED_BYTE_BASES[data[full] as usize][..tail]);
    }
    out
}

fn decode_kmer(bytes: &[u8], k: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(k);
    if k <= 32
        && let Some(chunk) = bytes.first_chunk::<8>()
    {
        let word = u64::from_le_bytes(*chunk);
        for i in 0..k {
            out.push(b"ACGT"[((word >> (62 - 2 * i)) & 3) as usize]);
        }
        return out;
    }
    let mut longs: Vec<u64> = Vec::with_capacity(bytes.len() / 8);
    for chunk in bytes.as_chunks::<8>().0 {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(chunk);
        longs.push(u64::from_le_bytes(arr));
    }
    for i in 0..k {
        let idx = i >> 5;
        let shift = 62 - ((i & 0x1F) << 1);
        let v = (longs[idx] >> shift) & 0x3;
        out.push(match v {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}

#[cfg(test)]
mod loader_order_tests {
    use super::*;

    #[test]
    fn abundant_nodes_use_node_keys_after_hash_table_reordering() {
        let first = b"ACCTCAGCCTCCCAAAATGCTAGGATTACAG";
        let second = b"CGGGTGTGGTGGCGGGCACCTGTAGTCCCAG";
        let canonical = |seq: &[u8]| {
            let (forward, reverse) = encode_kmer_pair(seq).unwrap();
            forward.min(reverse)
        };
        let base = 19;
        // The graph's hash-table serialization has the opposite order to
        // the node records, each of which supplies its own EC and head k-mer.
        let mut kmers = vec![second.to_vec(), first.to_vec()];
        let mut map = HashMap::from([(canonical(second), base), (canonical(first), base + 1)]);
        associate_abundant_node(&mut kmers, &mut map, base, base, first).unwrap();
        associate_abundant_node(&mut kmers, &mut map, base, base + 1, second).unwrap();
        assert_eq!(kmers, [first.to_vec(), second.to_vec()]);
        assert_eq!(map[&canonical(first)], base);
        assert_eq!(map[&canonical(second)], base + 1);
        let fallback = SelectiveFallbackKmerIndex::build(base, &[], &kmers);
        assert_eq!(fallback.lookup(canonical(first)).unwrap()[0].0, base);
        assert_eq!(fallback.lookup(canonical(second)).unwrap()[0].0, base + 1);
        assert!(
            associate_abundant_node(
                &mut kmers,
                &mut map,
                base,
                base,
                b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            )
            .is_err()
        );
    }

    #[test]
    fn singleton_bucket_storage_preserves_zero_specials_and_tie_order() {
        if usize::BITS != 64 {
            return;
        }
        let mut index = MinzPositionIndex::from_counts_inline(vec![0, 1, 2, 1, 0], &[3], true);
        assert!(index.inline.is_some());
        let mut cursors = index.cursors();
        assert!(cursors.is_empty());
        index.push_at(&mut cursors, 1, 0);
        index.push_at(&mut cursors, 2, 10);
        index.push_at(&mut cursors, 3, u64::MAX);
        index.push_at(&mut cursors, 2, 20);
        index.finish_inline();
        assert_eq!(index.get(0), &[]);
        assert_eq!(index.get(1), &[0]);
        assert_eq!(index.get(2), &[10, 20]);
        assert_eq!(index.get(3), &[u64::MAX]);
        assert_eq!(index.get(4), &[]);
        for bucket in 0..5 {
            assert_eq!(index.get_with_overcrowded(bucket).1, bucket == 3);
        }
    }

    #[test]
    fn descriptor_marker_bit_preserves_specials_and_high_singleton_bits() {
        if usize::BITS != 64 {
            return;
        }
        let counts = vec![1, 3, 1, 1];
        let mut inline = MinzPositionIndex::from_counts_inline(counts.clone(), &[2, 3], true);
        let mut csr = MinzPositionIndex::from_counts(counts);
        let mut inline_cursors = inline.cursors();
        let mut csr_cursors = csr.cursors();
        for (bucket, value) in [
            (0, (1u64 << 62) | 123),
            (1, 11),
            (1, u64::MAX),
            (1, 22),
            (2, 0xffff_ffff_0000_0001),
            (3, 0xffff_ffff_8000_0000),
        ] {
            inline.push_at(&mut inline_cursors, bucket, value);
            csr.push_at(&mut csr_cursors, bucket, value);
        }
        inline.finish_inline();
        for bucket in 0..4 {
            assert_eq!(
                inline.get_with_overcrowded(bucket),
                csr.get_with_overcrowded(bucket)
            );
            assert_eq!(
                inline.get_with_overcrowded(bucket).1,
                bucket == 1 || bucket == 3
            );
        }
    }

    #[test]
    fn ec_metadata_promotes_without_narrowing_large_ranges() {
        let mut blocks = FlatEcBlocks::new();
        blocks.push(1, 3, 0, 2);
        let Some(high) = (u32::MAX as usize).checked_add(1) else {
            return;
        };
        blocks.push(4, 9, high, high + 3);
        assert!(matches!(blocks, FlatEcBlocks::Wide(_)));
        assert_eq!(blocks.get(0).unwrap().ec_end, 2);
        assert_eq!(blocks.get(1).unwrap().ec_start, high);
        let mut offsets = MinzOffsets::Compact(vec![0, 2]);
        offsets.push(high);
        assert!(matches!(offsets, MinzOffsets::Wide(_)));
        assert_eq!(offsets.value_at(1), Some(2));
        assert_eq!(offsets.value_at(2), Some(high));
    }

    #[test]
    fn grouped_sequence_decode_equals_byte_oracle_for_all_tail_lengths() {
        let bytes: Vec<u8> = (0..=255).collect();
        for len in 0..=bytes.len() * 4 {
            let expected: Vec<u8> = (0..len)
                .map(|i| b"ACGT"[((bytes[i / 4] >> (2 * (i % 4))) & 3) as usize])
                .collect();
            assert_eq!(decode_compressed_sequence(&bytes, len), expected);
        }
        let mut word = 91u64;
        for _ in 0..1024 {
            word ^= word << 13;
            word ^= word >> 7;
            word ^= word << 17;
            for k in 0..=32 {
                let expected: Vec<u8> = (0..k)
                    .map(|i| b"ACGT"[((word >> (62 - 2 * i)) & 3) as usize])
                    .collect();
                assert_eq!(decode_kmer(&word.to_le_bytes(), k), expected);
            }
        }
    }

    #[test]
    fn compact_prefixes_and_wide_count_buffer_reuse_preserve_offsets() {
        for counts in [vec![], vec![0], vec![2, 0, 3, 1]] {
            let mut reserved = Vec::with_capacity(counts.len() + 1);
            reserved.extend_from_slice(&counts);
            let index = MinzPositionIndex::from_counts(reserved);
            assert!(matches!(index.offsets, MinzOffsets::Compact(_)));
            assert_eq!(index.len(), counts.len());
            assert_eq!(index.values.len(), counts.iter().sum());
            for (bucket, &count) in counts.iter().enumerate() {
                assert_eq!(index.get(bucket).len(), count);
            }
        }
        if let Some(total) = (u32::MAX as usize).checked_add(1) {
            let mut counts = Vec::with_capacity(3);
            counts.extend_from_slice(&[u32::MAX as usize, 1]);
            let pointer = counts.as_ptr();
            let (offsets, sum) = MinzOffsets::from_counts(counts);
            assert_eq!(sum, total);
            assert_eq!(offsets.bounds(0), (0, u32::MAX as usize));
            assert_eq!(offsets.bounds(1), (u32::MAX as usize, total));
            match offsets {
                MinzOffsets::Wide(v) => assert_eq!(v.as_ptr(), pointer),
                MinzOffsets::Compact(_) => panic!("large offsets must retain usize width"),
            }
        }
    }

    #[test]
    fn parallel_fill_preserves_every_bucket_in_serial_order() {
        let temp = tempfile::tempdir().unwrap();
        let fasta = temp.path().join("transcripts.fa");
        let mut state = 19u64;
        let mut text = String::new();
        for transcript in 0..5 {
            text.push_str(&format!(">tx{transcript}\n"));
            for _ in 0..2000 {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                text.push(b"ACGT"[(state & 3) as usize] as char);
            }
            text.push('\n');
        }
        std::fs::write(&fasta, text).unwrap();
        let path = temp.path().join("index.idx");
        crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default())
            .unwrap();
        let serial = build_bifrost_index_threaded(&path, 1).unwrap();
        let parallel = build_bifrost_index_threaded(&path, 4).unwrap();
        assert_eq!(
            serial.minz_positions.offsets,
            parallel.minz_positions.offsets
        );
        assert_eq!(serial.minz_positions.values, parallel.minz_positions.values);
        let compact = build_bifrost_index_quant_threaded(&path, false, 4).unwrap();
        assert!(compact.ec_blocks.is_empty());
        assert_eq!(serial.minz_positions.len(), compact.minz_positions.len());
        for bucket in 0..serial.minz_positions.len() {
            assert_eq!(
                serial.minz_positions.get(bucket),
                compact.minz_positions.get(bucket)
            );
            assert_eq!(
                serial.minz_positions.get_with_overcrowded(bucket),
                compact.minz_positions.get_with_overcrowded(bucket)
            );
        }
        for (uid, blocks) in serial.ec_blocks.iter().enumerate() {
            for (bid, block) in blocks.iter().enumerate() {
                assert_eq!(compact.flat_ec.ec(uid, bid), block.ec);
                assert_eq!(
                    compact.flat_ec.block_bounds(uid, bid),
                    Some((block.lb, block.ub))
                );
            }
            assert_eq!(compact.flat_ec.block_bounds(uid, blocks.len()), None);
        }
        assert_eq!(compact.flat_ec.block_bounds(usize::MAX, 0), None);
        let positional = build_bifrost_index_quant_threaded(&path, true, 2).unwrap();
        assert!(!positional.ec_blocks.is_empty());
        assert!(positional.has_strand_annotations());
        for seq in &serial.unitigs {
            for read in seq.windows(100).step_by(37) {
                for strand in [
                    super::super::Strand::Unstranded,
                    super::super::Strand::Reverse,
                ] {
                    let mut options = super::super::PseudoalignOptions::default();
                    options.investigation.trace_fast_path = true;
                    for mutate in [false, true] {
                        let mut query = read.to_vec();
                        if mutate {
                            query[40] = b'N';
                        }
                        let live = super::super::ec_for_read_bifrost(
                            &compact, &query, strand, None, options,
                        );
                        let mut debug = super::super::ReadDebugState::default();
                        let traced = super::super::ec_for_read_bifrost(
                            &compact,
                            &query,
                            strand,
                            Some(&mut debug),
                            options,
                        );
                        assert_eq!(format!("{live:?}"), format!("{traced:?}"));
                    }
                    for general in [false, true] {
                        let mut debug_full = super::super::ReadDebugState::default();
                        let mut debug_compact = super::super::ReadDebugState::default();
                        let full = super::super::ec_for_read_bifrost(
                            &serial,
                            read,
                            strand,
                            general.then_some(&mut debug_full),
                            super::super::PseudoalignOptions::default(),
                        );
                        let lean = super::super::ec_for_read_bifrost(
                            &compact,
                            read,
                            strand,
                            general.then_some(&mut debug_compact),
                            super::super::PseudoalignOptions::default(),
                        );
                        assert_eq!(format!("{full:?}"), format!("{lean:?}"));
                    }
                }
            }
        }
    }
}
