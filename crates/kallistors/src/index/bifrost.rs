use super::binary::{read_f64_le, read_i32_le, read_u16_le, read_u64_array, read_u64_le};
pub(crate) use super::bitmap::deserialize_compatible_to_vec as deserialize_roaring_to_vec;
use super::hashing::xorshift_next;
pub use super::hashing::{encode_minimizer_rep, minimizer_hash};
pub(crate) use super::hashing::{wyhash, wyhash16_words};

use std::collections::HashMap;
use std::io::Read;
#[cfg(test)]
use strength_reduce::StrengthReducedU64;

use crate::{Error, Result};

const MASK_CONTIG_TYPE: u64 = 0x8000_0000;
const MASK_CONTIG_POS: u64 = 0x7fff_ffff;
const FLAG_MASK: u64 = 0x7;
const POINTER_MASK: u64 = 0xffff_ffff_ffff_fff8;
const SHIFT_MASK_BITS: u64 = 3;
const FLAG_LOCAL_TINY: u64 = 0x0;
const FLAG_LOCAL_BITVECTOR: u64 = 0x1;
const FLAG_LOCAL_SINGLE: u64 = 0x2;
const FLAG_PTR_BITMAP: u64 = 0x3;

const TINY_BMP_MODE: u16 = 0x0000;
const TINY_LIST_MODE: u16 = 0x0002;
const TINY_RLE_MODE: u16 = 0x0004;
const TINY_MODE_MASK: u16 = 0x0006;
const TINY_SIZE_MASK: u16 = 0xfff8;

#[derive(Clone)]
pub struct BooPhf {
    _gamma: f64,
    nb_levels: i32,
    lastbitsetrank: u64,
    nelem: u64,
    levels: Vec<BooLevel>,
    final_hash: HashMap<[u8; 8], u64>,
    compact_levels: bool,
    absolute_levels: bool,
}

#[derive(Clone, Debug)]
pub struct MphfDebugInfo {
    pub level: u32,
    pub hash_raw: u64,
    pub hash_domain: u64,
    pub bucket: Option<u64>,
    pub rank: Option<u64>,
    pub final_hash_hit: bool,
}

// Word-aligned domains keep the low six bits unchanged. Their remaining
// numerator has only 58 bits, allowing a rounded-up reciprocal with a proven
// error below 1/divisor and no remainder correction. Other domains retain the
// full-u64 reciprocal and its single exact correction.
#[derive(Clone, Copy)]
pub(super) struct ReducedDomain {
    divisor: u64,
    reciprocal: u64,
    word_shift: Option<u32>,
}
impl ReducedDomain {
    pub(super) fn new(divisor: u64) -> Self {
        assert_ne!(divisor, 0);
        if divisor >= 4096 && divisor.is_multiple_of(64) {
            let words = divisor >> 6;
            let bits = u64::BITS - (words - 1).leading_zeros();
            let exponent = 58 + bits;
            let reciprocal = (1u128 << exponent).div_ceil(u128::from(words)) as u64;
            Self {
                divisor: words,
                reciprocal,
                word_shift: Some(exponent - 64),
            }
        } else {
            let reciprocal = if divisor == 1 {
                0
            } else {
                ((1u128 << 64) / u128::from(divisor)) as u64
            };
            Self {
                divisor,
                reciprocal,
                word_shift: None,
            }
        }
    }
    #[inline]
    pub(super) fn remainder(self, value: u64) -> u64 {
        if let Some(shift) = self.word_shift {
            let numerator = value >> 6;
            let quotient =
                (((u128::from(numerator) * u128::from(self.reciprocal)) >> 64) as u64) >> shift;
            return ((numerator - quotient * self.divisor) << 6) | (value & 63);
        }
        if self.divisor == 1 {
            return 0;
        }
        let quotient = ((u128::from(value) * u128::from(self.reciprocal)) >> 64) as u64;
        let remainder = value - quotient * self.divisor;
        if remainder >= self.divisor {
            remainder - self.divisor
        } else {
            remainder
        }
    }
}

#[derive(Clone)]
struct BooLevel {
    bitset: BitVector,
    hash_domain: u64,
    reduced_domain: Option<ReducedDomain>,
}

impl BooLevel {
    #[inline(always)]
    fn rank_for_hash<const COMPACT: bool, const ABSOLUTE: bool>(&self, hash: u64) -> Option<u64> {
        let bucket = self.bucket_for_hash(hash);
        if ABSOLUTE {
            self.bitset.get_rank_absolute(bucket)
        } else if COMPACT {
            self.bitset.get_rank_compact(bucket)
        } else {
            self.bitset.get_rank(bucket)
        }
    }

    #[inline(always)]
    fn bucket_for_hash(&self, hash: u64) -> u64 {
        match self.reduced_domain {
            Some(domain) => domain.remainder(hash),
            None => hash % self.hash_domain,
        }
    }
}

#[derive(Clone)]
#[repr(C, packed(2))]
struct RankedBitWord {
    bits: u64,
    rank: u16,
}

// Sixteen-byte, naturally aligned words avoid a second dependent base load.
// Smaller references keep the ten-byte compact representation.
#[derive(Clone)]
#[repr(C, align(16))]
struct AbsoluteRankedBitWord {
    bits: u64,
    rank: u64,
}

#[derive(Clone)]
struct BitVector {
    _size: u64,
    _nchar: u64,
    bit_array: Vec<u64>,
    ranks: Vec<u64>,
    ranked_words: Vec<RankedBitWord>,
    rank_bases: Vec<u32>,
    absolute_words: Vec<AbsoluteRankedBitWord>,
}

impl BitVector {
    #[cfg(test)]
    fn load<R: Read>(reader: &mut R) -> Result<Self> {
        Self::load_with_layout(reader, false)
    }

    fn load_with_layout<R: Read>(reader: &mut R, absolute: bool) -> Result<Self> {
        let size = read_u64_le(reader)?;
        let nchar = read_u64_le(reader)?;
        let bit_array = read_u64_array(reader, nchar as usize)?;
        let ranks_len = read_u64_le(reader)? as usize;
        let ranks = read_u64_array(reader, ranks_len)?;
        let mut ranked_words = Vec::with_capacity(if absolute { 0 } else { bit_array.len() });
        let mut absolute_words = Vec::with_capacity(if absolute { bit_array.len() } else { 0 });
        let mut completed = 0usize;
        let mut rank_bases = Vec::new();
        'blocks: for (words, &base) in bit_array.chunks(8).zip(&ranks) {
            let mut rank = base;
            for &bits in words {
                let Ok(compact_rank) = u32::try_from(rank) else {
                    // Very large rank domains retain the independent full-width
                    // bitmap/rank implementation, including its original order.
                    ranked_words = Vec::new();
                    absolute_words = Vec::new();
                    rank_bases = Vec::new();
                    break 'blocks;
                };
                if completed.is_multiple_of(1024) {
                    rank_bases.push(compact_rank);
                }
                // A valid block spans at most1023 preceding64-bit words:
                // its relative rank cannot exceed65472. Malformed/nonmonotone
                // rank tables keep the established full-width fallback.
                let Some(relative) = compact_rank
                    .checked_sub(*rank_bases.last().expect("rank block"))
                    .and_then(|rank| u16::try_from(rank).ok())
                else {
                    ranked_words = Vec::new();
                    absolute_words = Vec::new();
                    rank_bases = Vec::new();
                    break 'blocks;
                };
                if absolute {
                    absolute_words.push(AbsoluteRankedBitWord {
                        bits,
                        rank: u64::from(compact_rank),
                    });
                } else {
                    ranked_words.push(RankedBitWord {
                        bits,
                        rank: relative,
                    });
                }
                completed += 1;
                rank += u64::from(bits.count_ones());
            }
        }
        if absolute {
            rank_bases = Vec::new();
        }
        Ok(Self {
            _size: size,
            _nchar: nchar,
            bit_array,
            ranks,
            ranked_words,
            rank_bases,
            absolute_words,
        })
    }

    fn promote_absolute_words(&mut self) {
        if !self.absolute_words.is_empty() {
            assert!(self.absolute_words.len() as u64 >= self._size.div_ceil(64));
            return;
        }
        assert!(self.ranked_words.len() as u64 >= self._size.div_ceil(64));
        self.absolute_words = self
            .ranked_words
            .iter()
            .enumerate()
            .map(|(index, word)| AbsoluteRankedBitWord {
                bits: word.bits,
                rank: u64::from(self.compact_word_rank(index, word)),
            })
            .collect();
        self.ranked_words = Vec::new();
        self.rank_bases = Vec::new();
    }

    #[inline(always)]
    fn get_rank_absolute(&self, pos: u64) -> Option<u64> {
        let word = &self.absolute_words[(pos >> 6) as usize];
        let offset = (pos & 63) as u32;
        if word.bits & (1u64 << offset) == 0 {
            return None;
        }
        let before = if offset == 0 {
            0
        } else {
            (word.bits << (64 - offset)).count_ones()
        };
        Some(word.rank + u64::from(before))
    }

    #[inline(always)]
    fn rank_word_at(&self, pos: u64, word: &RankedBitWord) -> Option<u64> {
        let offset = (pos & 63) as u32;
        if word.bits & (1u64 << offset) == 0 {
            return None;
        }
        let before = if offset == 0 {
            0
        } else {
            (word.bits << (64 - offset)).count_ones()
        };
        Some(u64::from(self.compact_word_rank((pos >> 6) as usize, word)) + u64::from(before))
    }
    #[inline(always)]
    fn get_rank_compact(&self, pos: u64) -> Option<u64> {
        // The loader proved complete compact coverage before enabling this
        // representation. Ordinary Rust bounds checks remain authoritative.
        self.rank_word_at(pos, &self.ranked_words[(pos >> 6) as usize])
    }
    #[inline]
    fn compact_word_rank(&self, index: usize, word: &RankedBitWord) -> u32 {
        self.rank_bases[index >> 10] + u32::from(word.rank)
    }

    #[inline]
    fn get_rank(&self, pos: u64) -> Option<u64> {
        if !self.absolute_words.is_empty() {
            return self.get_rank_absolute(pos);
        }
        if let Some(word) = self.ranked_words.get((pos >> 6) as usize) {
            self.rank_word_at(pos, word)
        } else {
            self.get_rank_fallback(pos)
        }
    }

    // Compact rank words cover ordinary indexes. Keep the full-width path
    // available for large/incomplete tables without expanding each hot probe.
    #[cold]
    #[inline(never)]
    fn get_rank_fallback(&self, pos: u64) -> Option<u64> {
        self.get(pos).then(|| self.rank(pos))
    }

    fn get(&self, pos: u64) -> bool {
        if !self.absolute_words.is_empty() {
            return (self.absolute_words[(pos >> 6) as usize].bits >> (pos & 63)) & 1 != 0;
        }
        if self.bit_array.is_empty() && !self.ranked_words.is_empty() {
            return (self.ranked_words[(pos >> 6) as usize].bits >> (pos & 63)) & 1 != 0;
        }
        let idx = (pos >> 6) as usize;
        let bit = pos & 63;
        ((self.bit_array[idx] >> bit) & 1) != 0
    }

    fn domain(&self) -> u64 {
        self._size
    }

    fn rank(&self, pos: u64) -> u64 {
        if !self.absolute_words.is_empty() {
            let word = &self.absolute_words[(pos >> 6) as usize];
            let offset = (pos & 63) as u32;
            return word.rank
                + u64::from(if offset == 0 {
                    0
                } else {
                    (word.bits << (64 - offset)).count_ones()
                });
        }
        if self.bit_array.is_empty() && !self.ranked_words.is_empty() {
            let word = &self.ranked_words[(pos >> 6) as usize];
            let offset = (pos & 63) as u32;
            return u64::from(self.compact_word_rank((pos >> 6) as usize, word))
                + u64::from(if offset == 0 {
                    0
                } else {
                    (word.bits << (64 - offset)).count_ones()
                });
        }
        let word_idx = pos / 64;
        let word_offset = pos % 64;
        let block = pos / 512;
        let mut r = self.ranks[block as usize];
        let start = (block * 512) / 64;
        for w in start..word_idx {
            r += self.bit_array[w as usize].count_ones() as u64;
        }
        let final_word = self.bit_array[word_idx as usize];
        if word_offset > 0 {
            r += (final_word << (64 - word_offset)).count_ones() as u64;
        }
        r
    }
}

// Absence is independent of every rank value, including the full-u64 fallback.
// Avoid materializing a two-word Option for every intermediate loader lane.
pub(crate) struct BatchLookupResult<const LANES: usize> {
    ranks: [u64; LANES],
    found: u64,
}
pub(crate) struct BatchLookupIter<const LANES: usize> {
    ranks: std::array::IntoIter<u64, LANES>,
    found: u64,
}
impl<const LANES: usize> Iterator for BatchLookupIter<LANES> {
    type Item = Option<u64>;
    fn next(&mut self) -> Option<Self::Item> {
        self.ranks.next().map(|rank| {
            let present = self.found & 1 != 0;
            self.found >>= 1;
            present.then_some(rank)
        })
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.ranks.size_hint()
    }
}
impl<const LANES: usize> ExactSizeIterator for BatchLookupIter<LANES> {}
impl<const LANES: usize> IntoIterator for BatchLookupResult<LANES> {
    type Item = Option<u64>;
    type IntoIter = BatchLookupIter<LANES>;
    fn into_iter(self) -> Self::IntoIter {
        BatchLookupIter {
            ranks: self.ranks.into_iter(),
            found: self.found,
        }
    }
}
impl<const LANES: usize> BatchLookupResult<LANES> {
    // Borrowed found lanes retain original order and every full-u64 rank.
    // No sentinel rank value represents absence; only the found mask does.
    #[inline]
    pub(crate) fn found_ranks(&self) -> impl Iterator<Item = u64> + '_ {
        let mut found = self.found;
        std::iter::from_fn(move || {
            if found == 0 {
                return None;
            }
            let lane = found.trailing_zeros() as usize;
            found &= found - 1;
            Some(self.ranks[lane])
        })
    }
}

#[cfg(test)]
impl<const LANES: usize> BatchLookupResult<LANES> {
    fn into_options(self) -> [Option<u64>; LANES] {
        std::array::from_fn(|lane| (self.found & (1u64 << lane) != 0).then_some(self.ranks[lane]))
    }
}

// Keep absent lanes explicit: every possible eight-byte key remains valid.
pub(crate) struct MphfBatchKeys<const LANES: usize> {
    words: [u64; LANES],
    present: u64,
}
impl<const LANES: usize> Default for MphfBatchKeys<LANES> {
    fn default() -> Self {
        assert!(LANES <= 64);
        Self {
            words: [0; LANES],
            present: 0,
        }
    }
}
impl<const LANES: usize> MphfBatchKeys<LANES> {
    pub(crate) fn new() -> Self {
        Self::default()
    }
    pub(crate) fn set(&mut self, lane: usize, key: Option<[u8; 8]>) {
        assert!(lane < LANES);
        if let Some(key) = key {
            self.words[lane] = u64::from_le_bytes(key);
            self.present |= 1u64 << lane;
        } else {
            self.present &= !(1u64 << lane);
        }
    }
    pub(crate) fn clear(&mut self) {
        self.present = 0;
    }
    pub(crate) fn iter(&self) -> impl ExactSizeIterator<Item = Option<[u8; 8]>> + '_ {
        self.words
            .iter()
            .enumerate()
            .map(|(lane, word)| (self.present & (1u64 << lane) != 0).then_some(word.to_le_bytes()))
    }
}
#[cfg(test)]
impl<const LANES: usize> From<&[Option<[u8; 8]>; LANES]> for MphfBatchKeys<LANES> {
    fn from(keys: &[Option<[u8; 8]>; LANES]) -> Self {
        let mut batch = Self::new();
        for (lane, key) in keys.iter().copied().enumerate() {
            batch.set(lane, key);
        }
        batch
    }
}

impl BooPhf {
    // Compact words and their small block bases retain the complete bitmap/rank.
    // Wider or incomplete tables retain their original full-width storage.
    pub(crate) fn retire_storage(self) {
        for level in self.levels {
            super::scratch_pages::retire(level.bitset.bit_array);
            super::scratch_pages::retire(level.bitset.ranks);
            super::scratch_pages::retire(level.bitset.ranked_words);
            super::scratch_pages::retire(level.bitset.rank_bases);
            super::scratch_pages::retire(level.bitset.absolute_words);
        }
    }

    pub(crate) fn release_redundant_rank_storage(&mut self) {
        for level in &mut self.levels {
            let bits = &mut level.bitset;
            if (bits.ranked_words.len() == bits._nchar as usize && !bits.ranked_words.is_empty())
                || (bits.absolute_words.len() as u64 >= bits._size.div_ceil(64)
                    && !bits.absolute_words.is_empty())
            {
                bits.bit_array = Vec::new();
                bits.ranks = Vec::new();
            }
        }
    }
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let gamma = read_f64_le(reader)?;
        let nb_levels = read_i32_le(reader)?;
        let lastbitsetrank = read_u64_le(reader)?;
        let nelem = read_u64_le(reader)?;
        let mut levels = Vec::with_capacity(nb_levels as usize);
        for level in 0..nb_levels {
            let prefer_absolute = nelem >= (1 << 28) && level < nb_levels.saturating_sub(1);
            let bitset = BitVector::load_with_layout(reader, prefer_absolute)?;
            let hash_domain = bitset.domain();
            levels.push(BooLevel {
                bitset,
                hash_domain,
                reduced_domain: (hash_domain != 0).then(|| ReducedDomain::new(hash_domain)),
            });
        }
        let final_hash_size = read_u64_le(reader)? as usize;
        let mut final_hash = HashMap::with_capacity(final_hash_size);
        for _ in 0..final_hash_size {
            let key = read_u64_le(reader)?;
            let val = read_u64_le(reader)?;
            final_hash.insert(key.to_le_bytes(), val);
        }
        let compact_levels = levels
            .iter()
            .take(nb_levels.saturating_sub(1).max(0) as usize)
            .all(|level| {
                level.hash_domain != 0
                    && (level.bitset.ranked_words.len() as u64 >= level.hash_domain.div_ceil(64)
                        || level.bitset.absolute_words.len() as u64
                            >= level.hash_domain.div_ceil(64))
            });
        let absolute_levels = nelem >= (1 << 28) && compact_levels;
        if absolute_levels {
            for level in levels
                .iter_mut()
                .take(nb_levels.saturating_sub(1).max(0) as usize)
            {
                level.bitset.promote_absolute_words();
            }
        }
        let phf = Self {
            _gamma: gamma,
            nb_levels,
            lastbitsetrank,
            nelem,
            levels,
            final_hash,
            compact_levels,
            absolute_levels,
        };
        Ok(phf)
    }

    #[inline]
    pub fn lookup(&self, key: &[u8; 8]) -> Option<u64> {
        if self.absolute_levels {
            self.lookup_impl::<true, true>(key)
        } else if self.compact_levels {
            self.lookup_impl::<true, false>(key)
        } else {
            self.lookup_impl::<false, false>(key)
        }
    }
    #[inline(never)]
    fn lookup_impl<const COMPACT: bool, const ABSOLUTE: bool>(&self, key: &[u8; 8]) -> Option<u64> {
        let base_hash = minimizer_hash(key, 0);
        let mut s0 = base_hash ^ 0xAAAAAAAA55555555;
        let mut s1 = base_hash ^ 0x33333333CCCCCCCC;
        let max_level = self.nb_levels.saturating_sub(1).max(0) as usize;
        let mut levels = self.levels.iter().take(max_level);
        // The first two seeds do not depend on the xorshift recurrence.
        if let Some(level) = levels.next()
            && let Some(rank) = level.rank_for_hash::<COMPACT, ABSOLUTE>(s0)
        {
            return Some(rank);
        }
        if let Some(level) = levels.next()
            && let Some(rank) = level.rank_for_hash::<COMPACT, ABSOLUTE>(s1)
        {
            return Some(rank);
        }
        for level in levels {
            if let Some(rank) =
                level.rank_for_hash::<COMPACT, ABSOLUTE>(xorshift_next(&mut s0, &mut s1))
            {
                return Some(rank);
            }
        }
        self.lookup_final(key)
    }
    #[cold]
    #[inline(never)]
    fn lookup_final(&self, key: &[u8; 8]) -> Option<u64> {
        self.final_hash.get(key).map(|v| v + self.lastbitsetrank)
    }

    // Independent loader probes can advance through one MPHF level together,
    // exposing unrelated rank-table loads without changing any hash or rank.
    // Missing/tail lanes never inspect the MPHF, and output order is unchanged.
    pub(crate) fn lookup_batch<const LANES: usize>(
        &self,
        keys: &MphfBatchKeys<LANES>,
    ) -> BatchLookupResult<LANES> {
        if self.absolute_levels {
            self.lookup_batch_impl::<LANES, true>(keys)
        } else {
            self.lookup_batch_impl::<LANES, false>(keys)
        }
    }

    fn lookup_batch_impl<const LANES: usize, const ABSOLUTE: bool>(
        &self,
        keys: &MphfBatchKeys<LANES>,
    ) -> BatchLookupResult<LANES> {
        assert!(LANES <= 64);
        let mut result = [0u64; LANES];
        let mut found = 0u64;
        let mut first = [0; LANES];
        let mut second = [0; LANES];
        let mut pending = keys.present;
        if ABSOLUTE && LANES == 64 && pending == u64::MAX {
            // Complete large-reference batches have no holes. A contiguous
            // traversal exposes independent exact hashes to the compiler;
            // sparse, missing and tail lanes retain mask traversal below.
            for (&word, (first, second)) in keys
                .words
                .iter()
                .zip(first.iter_mut().zip(second.iter_mut()))
            {
                let base = minimizer_hash(&word.to_le_bytes(), 0);
                *first = base ^ 0xAAAAAAAA55555555;
                *second = base ^ 0x33333333CCCCCCCC;
            }
        } else {
            let mut to_hash = pending;
            while to_hash != 0 {
                let lane = to_hash.trailing_zeros() as usize;
                to_hash &= to_hash - 1;
                let base = minimizer_hash(&keys.words[lane].to_le_bytes(), 0);
                first[lane] = base ^ 0xAAAAAAAA55555555;
                second[lane] = base ^ 0x33333333CCCCCCCC;
            }
        }
        let max_level = self.nb_levels.saturating_sub(1).max(0) as usize;
        let mut buckets = [0; LANES];
        let mut words = [0u64; LANES];
        let mut ranks = [0u32; LANES];
        for (number, level) in self.levels.iter().take(max_level).enumerate() {
            if pending == 0 {
                break;
            }
            // Resolve arithmetic for every independent lane before the
            // random rank-word loads. The subsequent gather keeps lane order.
            if ABSOLUTE
                && let Some(ReducedDomain {
                    divisor,
                    reciprocal,
                    word_shift: Some(shift),
                }) = level.reduced_domain
            {
                let mut address_lanes = pending;
                while address_lanes != 0 {
                    let lane = address_lanes.trailing_zeros() as usize;
                    address_lanes &= address_lanes - 1;
                    let hash = match number {
                        0 => first[lane],
                        1 => second[lane],
                        _ => xorshift_next(&mut first[lane], &mut second[lane]),
                    };
                    buckets[lane] = {
                        let numerator = hash >> 6;
                        let quotient = (((u128::from(numerator) * u128::from(reciprocal)) >> 64)
                            as u64)
                            >> shift;
                        ((numerator - quotient * divisor) << 6) | (hash & 63)
                    };
                }
            } else {
                let mut address_lanes = pending;
                while address_lanes != 0 {
                    let lane = address_lanes.trailing_zeros() as usize;
                    address_lanes &= address_lanes - 1;
                    let hash = match number {
                        0 => first[lane],
                        1 => second[lane],
                        _ => xorshift_next(&mut first[lane], &mut second[lane]),
                    };
                    buckets[lane] = level.bucket_for_hash(hash);
                }
            }
            let mut loaded = 0u64;
            // Gather every independent rank word before branching on hits.
            // Testing each lane immediately would serialize the random loads.
            let mut to_load = pending;
            while to_load != 0 {
                let lane = to_load.trailing_zeros() as usize;
                to_load &= to_load - 1;
                let bucket = buckets[lane];
                if ABSOLUTE {
                    let word = &level.bitset.absolute_words[(bucket >> 6) as usize];
                    words[lane] = word.bits;
                    ranks[lane] = word.rank as u32;
                } else if let Some(word) = level.bitset.ranked_words.get((bucket >> 6) as usize) {
                    words[lane] = word.bits;
                    ranks[lane] = level.bitset.compact_word_rank((bucket >> 6) as usize, word);
                    loaded |= 1 << lane;
                }
            }
            let mut to_resolve = pending;
            while to_resolve != 0 {
                let lane = to_resolve.trailing_zeros() as usize;
                to_resolve &= to_resolve - 1;
                let bucket = buckets[lane];
                let bits = words[lane];
                let base = ranks[lane];
                // Absolute promotion proves full hash-domain coverage, so
                // every active lane was loaded above. Only compact layouts
                // need the explicit mask and legacy rank fallback.
                let rank = if ABSOLUTE || loaded & (1 << lane) != 0 {
                    let offset = (bucket & 63) as u32;
                    if bits & (1u64 << offset) == 0 {
                        None
                    } else {
                        let before = if offset == 0 {
                            0
                        } else {
                            (bits << (64 - offset)).count_ones()
                        };
                        Some(u64::from(base) + u64::from(before))
                    }
                } else {
                    level.bitset.get_rank_fallback(bucket)
                };
                if let Some(rank) = rank {
                    result[lane] = rank;
                    found |= 1u64 << lane;
                    pending &= !(1 << lane);
                }
            }
        }
        // Valid loader probes usually resolve before the final hash table.
        // Absent lanes need neither key materialization nor a table probe.
        while pending != 0 {
            let lane = pending.trailing_zeros() as usize;
            pending &= pending - 1;
            if let Some(value) = self.final_hash.get(&keys.words[lane].to_le_bytes()) {
                result[lane] = value + self.lastbitsetrank;
                found |= 1u64 << lane;
            }
        }
        BatchLookupResult {
            ranks: result,
            found,
        }
    }

    pub fn debug_lookup(&self, key: &[u8; 8]) -> (Option<u64>, MphfDebugInfo) {
        let mut s0 = 0u64;
        let mut s1 = 0u64;
        let mut level = 0;
        let mut hash_raw = 0u64;
        let max_level = self.nb_levels.saturating_sub(1);
        for i in 0..max_level {
            if i == 0 {
                s0 = minimizer_hash(key, 0xAAAAAAAA55555555);
                hash_raw = s0;
            } else if i == 1 {
                s1 = minimizer_hash(key, 0x33333333CCCCCCCC);
                hash_raw = s1;
            } else {
                let next = xorshift_next(&mut s0, &mut s1);
                hash_raw = next;
            }
            if self.levels[i as usize]
                .bitset
                .get(hash_raw % self.levels[i as usize].hash_domain)
            {
                break;
            }
            level += 1;
        }

        if level == (self.nb_levels - 1) {
            let final_hash_hit = self.final_hash.get(key).copied();
            let rank = final_hash_hit.map(|v| v + self.lastbitsetrank);
            let info = MphfDebugInfo {
                level: level as u32,
                hash_raw,
                hash_domain: 0,
                bucket: None,
                rank,
                final_hash_hit: final_hash_hit.is_some(),
            };
            return (rank, info);
        }

        let level_usize = level as usize;
        let hash_domain = self.levels[level_usize].hash_domain;
        let bucket = hash_raw % hash_domain;
        let rank = self.levels[level_usize].bitset.rank(bucket);
        let info = MphfDebugInfo {
            level: level as u32,
            hash_raw,
            hash_domain,
            bucket: Some(bucket),
            rank: Some(rank),
            final_hash_hit: false,
        };
        (Some(rank), info)
    }

    pub fn size(&self) -> u64 {
        self.nelem
    }

    pub fn final_hash_entries(&self) -> Vec<([u8; 8], u64)> {
        self.final_hash.iter().map(|(k, v)| (*k, *v)).collect()
    }
}

pub fn read_bitcontainer_values<R: Read>(reader: &mut R) -> Result<Vec<u32>> {
    let set_bits = read_u64_le(reader)?;
    let flag = set_bits & FLAG_MASK;
    match flag {
        FLAG_PTR_BITMAP => {
            let expected_sz = (set_bits >> SHIFT_MASK_BITS) as usize;
            let mut buf = vec![0u8; expected_sz];
            reader.read_exact(&mut buf)?;
            let vals = unsafe { deserialize_roaring_to_vec(&buf) }
                .ok_or_else(|| Error::InvalidFormat("invalid roaring bitmap".into()))?;
            Ok(vals)
        }
        FLAG_LOCAL_TINY => read_tinybitmap_values(reader),
        FLAG_LOCAL_SINGLE => Ok(vec![(set_bits >> SHIFT_MASK_BITS) as u32]),
        FLAG_LOCAL_BITVECTOR => {
            let mut out = Vec::new();
            let bits = set_bits & POINTER_MASK;
            for i in 0..61u64 {
                if ((bits >> (i + SHIFT_MASK_BITS)) & 1) != 0 {
                    out.push(i as u32);
                }
            }
            Ok(out)
        }
        _ => Err(Error::InvalidFormat("invalid bitcontainer flag".into())),
    }
}

fn read_tinybitmap_values<R: Read>(reader: &mut R) -> Result<Vec<u32>> {
    let header = read_u16_le(reader)?;
    let sz = (header & TINY_SIZE_MASK) >> 3;
    if sz == 0 {
        return Ok(Vec::new());
    }
    let mut buf = vec![0u16; sz as usize];
    buf[0] = header;
    for slot in buf.iter_mut().skip(1) {
        *slot = read_u16_le(reader)?;
    }
    Ok(tinybitmap_iter_values(&buf))
}

fn tinybitmap_iter_values(bmp: &[u16]) -> Vec<u32> {
    let sz = (bmp[0] & TINY_SIZE_MASK) >> 3;
    let mode = bmp[0] & TINY_MODE_MASK;
    let card = bmp[1];
    let offset = (bmp[2] as u32) << 16;
    if card == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(card as usize);
    if mode == TINY_BMP_MODE {
        let mut i = 2usize;
        while (i as u16) < sz && out.len() < card as usize {
            i += 1;
            let mut j = 0u32;
            let mut e = bmp[i];
            while e != 0 && out.len() < card as usize {
                if (e & 1) != 0 {
                    let val = offset | (((i - 3) as u32) << 4) | j;
                    out.push(val);
                }
                e >>= 1;
                j += 1;
            }
        }
    } else if mode == TINY_LIST_MODE {
        for i in 0..card as usize {
            out.push(offset | bmp[i + 3] as u32);
        }
    } else if mode == TINY_RLE_MODE {
        let mut i = 3usize;
        let mut j = 4usize;
        let mut val = offset | bmp[i] as u32;
        while i < (card as usize + 3) && out.len() < card as usize {
            out.push(val);
            val += 1;
            if (val & 0xffff) > bmp[j] as u32 {
                i += 2;
                j += 2;
                if i >= card as usize + 3 {
                    break;
                }
                val = offset | bmp[i] as u32;
            }
        }
    }
    out
}

pub fn decode_pos_id(pos_id: u64) -> (u32, u32, bool) {
    let unitig = (pos_id >> 32) as u32;
    let low = pos_id & 0xffff_ffff;
    let is_km = (low & MASK_CONTIG_TYPE) != 0;
    let pos = (low & MASK_CONTIG_POS) as u32;
    (unitig, pos, is_km)
}

#[cfg(test)]
mod mphf_lookup_tests;

#[cfg(test)]
mod exact_rep_hash_tests;
