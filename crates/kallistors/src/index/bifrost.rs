use std::collections::HashMap;
use std::io::Read;
#[cfg(test)]
use strength_reduce::StrengthReducedU64;

use croaring_sys::{
    roaring_bitmap_deserialize, roaring_bitmap_free, roaring_bitmap_get_cardinality,
    roaring_bitmap_portable_deserialize_safe, roaring_bitmap_to_uint32_array,
};

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

const WYHASH_SECRET: [u64; 4] = [
    0xa0761d6478bd642f,
    0xe7037ed1a0b428db,
    0x8ebc6af09c88c6e3,
    0x589965cc75374cc3,
];
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

#[derive(Clone)]
pub struct BooPhf {
    _gamma: f64,
    nb_levels: i32,
    lastbitsetrank: u64,
    nelem: u64,
    levels: Vec<BooLevel>,
    final_hash: HashMap<[u8; 8], u64>,
    compact_levels: bool,
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
    fn rank_for_hash<const COMPACT: bool>(&self, hash: u64) -> Option<u64> {
        let bucket = self.bucket_for_hash(hash);
        if COMPACT {
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

#[derive(Clone)]
struct BitVector {
    _size: u64,
    _nchar: u64,
    bit_array: Vec<u64>,
    ranks: Vec<u64>,
    ranked_words: Vec<RankedBitWord>,
    rank_bases: Vec<u32>,
}

impl BitVector {
    fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let size = read_u64_le(reader)?;
        let nchar = read_u64_le(reader)?;
        let bit_array = read_u64_array(reader, nchar as usize)?;
        let ranks_len = read_u64_le(reader)? as usize;
        let ranks = read_u64_array(reader, ranks_len)?;
        let mut ranked_words = Vec::with_capacity(bit_array.len());
        let mut rank_bases = Vec::new();
        'blocks: for (words, &base) in bit_array.chunks(8).zip(&ranks) {
            let mut rank = base;
            for &bits in words {
                let Ok(compact_rank) = u32::try_from(rank) else {
                    // Very large rank domains retain the independent full-width
                    // bitmap/rank implementation, including its original order.
                    ranked_words = Vec::new();
                    rank_bases = Vec::new();
                    break 'blocks;
                };
                if ranked_words.len().is_multiple_of(1024) {
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
                    rank_bases = Vec::new();
                    break 'blocks;
                };
                ranked_words.push(RankedBitWord {
                    bits,
                    rank: relative,
                });
                rank += u64::from(bits.count_ones());
            }
        }
        Ok(Self {
            _size: size,
            _nchar: nchar,
            bit_array,
            ranks,
            ranked_words,
            rank_bases,
        })
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

impl BooPhf {
    // Compact words and their small block bases retain the complete bitmap/rank.
    // Wider or incomplete tables retain their original full-width storage.
    pub(crate) fn retire_storage(self) {
        for level in self.levels {
            super::scratch_pages::retire(level.bitset.bit_array);
            super::scratch_pages::retire(level.bitset.ranks);
            super::scratch_pages::retire(level.bitset.ranked_words);
            super::scratch_pages::retire(level.bitset.rank_bases);
        }
    }

    pub(crate) fn release_redundant_rank_storage(&mut self) {
        for level in &mut self.levels {
            let bits = &mut level.bitset;
            if bits.ranked_words.len() == bits._nchar as usize && !bits.ranked_words.is_empty() {
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
        for _ in 0..nb_levels {
            let bitset = BitVector::load(reader)?;
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
                    && level.bitset.ranked_words.len() as u64 >= level.hash_domain.div_ceil(64)
            });
        let phf = Self {
            _gamma: gamma,
            nb_levels,
            lastbitsetrank,
            nelem,
            levels,
            final_hash,
            compact_levels,
        };
        Ok(phf)
    }

    #[inline]
    pub fn lookup(&self, key: &[u8; 8]) -> Option<u64> {
        if self.compact_levels {
            self.lookup_impl::<true>(key)
        } else {
            self.lookup_impl::<false>(key)
        }
    }
    #[inline(never)]
    fn lookup_impl<const COMPACT: bool>(&self, key: &[u8; 8]) -> Option<u64> {
        let base_hash = minimizer_hash(key, 0);
        let mut s0 = base_hash ^ 0xAAAAAAAA55555555;
        let mut s1 = base_hash ^ 0x33333333CCCCCCCC;
        let max_level = self.nb_levels.saturating_sub(1).max(0) as usize;
        let mut levels = self.levels.iter().take(max_level);
        // The first two seeds do not depend on the xorshift recurrence.
        if let Some(level) = levels.next()
            && let Some(rank) = level.rank_for_hash::<COMPACT>(s0)
        {
            return Some(rank);
        }
        if let Some(level) = levels.next()
            && let Some(rank) = level.rank_for_hash::<COMPACT>(s1)
        {
            return Some(rank);
        }
        for level in levels {
            if let Some(rank) = level.rank_for_hash::<COMPACT>(xorshift_next(&mut s0, &mut s1)) {
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
    pub(crate) fn lookup_batch8(&self, keys: &[Option<[u8; 8]>; 8]) -> [Option<u64>; 8] {
        let mut result = [None; 8];
        let mut first = [0; 8];
        let mut second = [0; 8];
        let mut pending = 0u8;
        for (lane, key) in keys.iter().enumerate() {
            if let Some(key) = key {
                let base = minimizer_hash(key, 0);
                first[lane] = base ^ 0xAAAAAAAA55555555;
                second[lane] = base ^ 0x33333333CCCCCCCC;
                pending |= 1 << lane;
            }
        }
        let max_level = self.nb_levels.saturating_sub(1).max(0) as usize;
        let mut buckets = [0; 8];
        let mut words = [0u64; 8];
        let mut ranks = [0u32; 8];
        for (number, level) in self.levels.iter().take(max_level).enumerate() {
            if pending == 0 {
                break;
            }
            let mut loaded = 0u8;
            // Gather every independent rank word before branching on hits.
            // Testing each lane immediately would serialize the random loads.
            for (lane, (bucket, (bits, rank))) in buckets
                .iter_mut()
                .zip(words.iter_mut().zip(ranks.iter_mut()))
                .enumerate()
            {
                if pending & (1 << lane) != 0 {
                    let hash = match number {
                        0 => first[lane],
                        1 => second[lane],
                        _ => xorshift_next(&mut first[lane], &mut second[lane]),
                    };
                    *bucket = level.bucket_for_hash(hash);
                    if let Some(word) = level.bitset.ranked_words.get((*bucket >> 6) as usize) {
                        *bits = word.bits;
                        *rank = level
                            .bitset
                            .compact_word_rank((*bucket >> 6) as usize, word);
                        loaded |= 1 << lane;
                    }
                }
            }
            for (lane, (&bucket, (&bits, &base))) in
                buckets.iter().zip(words.iter().zip(&ranks)).enumerate()
            {
                if pending & (1 << lane) == 0 {
                    continue;
                }
                let rank = if loaded & (1 << lane) != 0 {
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
                    result[lane] = Some(rank);
                    pending &= !(1 << lane);
                }
            }
        }
        for (lane, (result, key)) in result.iter_mut().zip(keys).enumerate() {
            if pending & (1 << lane) != 0 {
                *result = self
                    .final_hash
                    .get(key.as_ref().expect("pending query has a key"))
                    .map(|value| value + self.lastbitsetrank);
            }
        }
        result
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

pub fn encode_minimizer_rep(seq: &[u8]) -> Option<[u8; 8]> {
    let (fwd, rev) = encode_word_pair(seq)?;
    let rep = fwd.min(rev);
    Some(rep.to_le_bytes())
}

/// The serialized minimizer key is exactly eight little-endian bytes.
/// This is the fixed-length form of wyhash(key,0), followed by the original
/// external seed XOR. Keep generic byte hashing as the independent oracle.
#[inline]
pub fn minimizer_hash(key: &[u8; 8], seed: u64) -> u64 {
    let word = u64::from_le_bytes(*key);
    wymix(
        WYHASH_SECRET[1] ^ 8,
        wymix(
            word.rotate_left(32) ^ WYHASH_SECRET[1],
            word ^ WYHASH_SECRET[0],
        ),
    ) ^ seed
}

fn encode_word_pair(seq: &[u8]) -> Option<(u64, u64)> {
    if (1..=32).contains(&seq.len()) {
        let (forward, reverse) = crate::dna::encode_kmer_pair(seq)?;
        let shift = 64 - 2 * seq.len();
        return Some((forward << shift, reverse << shift));
    }

    if seq.is_empty() {
        return Some((0, 0));
    }
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let last = seq.len() - 1;
    for (i, &b) in seq.iter().enumerate() {
        let v = BASE_CODES[b as usize];
        if v == INVALID_BASE_CODE {
            return None;
        }
        let v = u64::from(v);
        let shift = 62 - ((i & 0x1f) << 1);
        fwd |= v << shift;
        let rev_shift = 62 - (((last - i) & 0x1f) << 1);
        rev |= (3 - v) << rev_shift;
    }
    Some((fwd, rev))
}

fn xorshift_next(s0: &mut u64, s1: &mut u64) -> u64 {
    let mut x = *s0;
    let y = *s1;
    *s0 = y;
    x ^= x << 23;
    *s1 = x ^ y ^ (x >> 17) ^ (y >> 26);
    s1.wrapping_add(y)
}

/// Exact seed-zero wyhash for two little-endian words (the16-byte REP key).
/// Keep the generic byte hash as the independent format/compatibility oracle.
#[inline]
pub(crate) fn wyhash16_words(first: u64, second: u64) -> u64 {
    let a = (first << 32) | u64::from(second as u32);
    let b = (second & 0xffff_ffff_0000_0000) | (first >> 32);
    wymix(
        WYHASH_SECRET[1] ^ 16,
        wymix(a ^ WYHASH_SECRET[1], b ^ WYHASH_SECRET[0]),
    )
}

pub(crate) fn wyhash(key: &[u8], seed: u64) -> u64 {
    let mut seed = seed ^ WYHASH_SECRET[0];
    let len = key.len();
    let a: u64;
    let b: u64;
    if len <= 16 {
        if len >= 4 {
            a = (wyr4(&key[0..]) << 32) | wyr4(&key[((len >> 3) << 2)..]);
            b = (wyr4(&key[len - 4..]) << 32) | wyr4(&key[len - 4 - ((len >> 3) << 2)..]);
        } else if len > 0 {
            a = wyr3(key, len);
            b = 0;
        } else {
            a = 0;
            b = 0;
        }
    } else {
        let mut i = len;
        let mut p = 0usize;
        if i > 48 {
            let mut see1 = seed;
            let mut see2 = seed;
            while i > 48 {
                seed = wymix(
                    wyr8(&key[p..]) ^ WYHASH_SECRET[1],
                    wyr8(&key[p + 8..]) ^ seed,
                );
                see1 = wymix(
                    wyr8(&key[p + 16..]) ^ WYHASH_SECRET[2],
                    wyr8(&key[p + 24..]) ^ see1,
                );
                see2 = wymix(
                    wyr8(&key[p + 32..]) ^ WYHASH_SECRET[3],
                    wyr8(&key[p + 40..]) ^ see2,
                );
                p += 48;
                i -= 48;
            }
            seed ^= see1 ^ see2;
        }
        while i > 16 {
            seed = wymix(
                wyr8(&key[p..]) ^ WYHASH_SECRET[1],
                wyr8(&key[p + 8..]) ^ seed,
            );
            p += 16;
            i -= 16;
        }
        a = wyr8(&key[p + i - 16..]);
        b = wyr8(&key[p + i - 8..]);
    }
    wymix(
        WYHASH_SECRET[1] ^ len as u64,
        wymix(a ^ WYHASH_SECRET[1], b ^ seed),
    )
}

fn wymix(a: u64, b: u64) -> u64 {
    let r = (a as u128).wrapping_mul(b as u128);
    let lo = r as u64;
    let hi = (r >> 64) as u64;
    lo ^ hi
}

fn wyr8(p: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    buf.copy_from_slice(&p[0..8]);
    u64::from_le_bytes(buf)
}

fn wyr4(p: &[u8]) -> u64 {
    let mut buf = [0u8; 4];
    buf.copy_from_slice(&p[0..4]);
    u64::from(u32::from_le_bytes(buf))
}

fn wyr3(p: &[u8], k: usize) -> u64 {
    ((p[0] as u64) << 16) | ((p[k >> 1] as u64) << 8) | (p[k - 1] as u64)
}

pub(crate) unsafe fn deserialize_roaring_to_vec(buf: &[u8]) -> Option<Vec<u32>> {
    let ptr = unsafe { roaring_bitmap_portable_deserialize_safe(buf.as_ptr().cast(), buf.len()) };
    let ptr = if ptr.is_null() {
        unsafe { roaring_bitmap_deserialize(buf.as_ptr().cast()) }
    } else {
        ptr
    };
    if ptr.is_null() {
        return None;
    }
    let card = unsafe { roaring_bitmap_get_cardinality(ptr) };
    let mut out = vec![0u32; card as usize];
    unsafe {
        roaring_bitmap_to_uint32_array(ptr, out.as_mut_ptr());
        roaring_bitmap_free(ptr);
    }
    Some(out)
}

fn read_u64_array<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u64>> {
    let mut words = vec![0u64; len];
    // SAFETY: every u64 bit pattern is valid. The byte view spans exactly
    // initialized exclusive Vec storage, and cannot outlive that borrow.
    let bytes = unsafe {
        std::slice::from_raw_parts_mut(
            words.as_mut_ptr().cast::<u8>(),
            std::mem::size_of_val(words.as_slice()),
        )
    };
    reader.read_exact(bytes)?;
    #[cfg(target_endian = "big")]
    for word in &mut words {
        *word = u64::from_le(*word);
    }
    Ok(words)
}

fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u16_le<R: Read>(reader: &mut R) -> Result<u16> {
    let mut buf = [0u8; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_i32_le<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_f64_le<R: Read>(reader: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

pub fn decode_pos_id(pos_id: u64) -> (u32, u32, bool) {
    let unitig = (pos_id >> 32) as u32;
    let low = pos_id & 0xffff_ffff;
    let is_km = (low & MASK_CONTIG_TYPE) != 0;
    let pos = (low & MASK_CONTIG_POS) as u32;
    (unitig, pos, is_km)
}

#[cfg(test)]
mod mphf_lookup_tests {
    use super::*;

    #[test]
    fn bulk_bitvector_words_preserve_little_endian_and_truncated_input_errors() {
        let words = [0u64, 1, u64::MAX, 1 << 63, 0x0123456789abcdef];
        let bytes: Vec<_> = words.iter().flat_map(|word| word.to_le_bytes()).collect();
        assert_eq!(
            read_u64_array(&mut std::io::Cursor::new(&bytes), words.len()).unwrap(),
            words
        );
        assert!(
            read_u64_array(&mut std::io::Cursor::new([]), 0)
                .unwrap()
                .is_empty()
        );
        for end in 0..bytes.len() {
            assert!(read_u64_array(&mut std::io::Cursor::new(&bytes[..end]), words.len()).is_err());
        }
    }

    #[test]
    fn reduced_domains_preserve_exact_u64_remainders_at_boundaries() {
        let mut state = 17u64;
        for domain in [
            1,
            2,
            3,
            7,
            64,
            65,
            4096,
            4160,
            8128,
            8192,
            8256,
            1 << 20,
            50_102_144,
            50_102_171,
            u32::MAX as u64,
            1u64 << 63,
            u64::MAX - 63,
            u64::MAX,
        ] {
            let reduced = StrengthReducedU64::new(domain);
            let reciprocal = ReducedDomain::new(domain);
            for value in [0, 1, domain - 1, domain, domain.saturating_add(1), u64::MAX] {
                assert_eq!(value % reduced, value % domain);
                assert_eq!(reciprocal.remainder(value), value % domain);
            }
            for _ in 0..20_000 {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                assert_eq!(state % reduced, state % domain);
                assert_eq!(reciprocal.remainder(state), state % domain);
            }
        }
    }

    #[test]
    fn lookup_matches_diagnostic_oracle_across_hash_levels() {
        assert_eq!(std::mem::size_of::<RankedBitWord>(), 10);
        // Independent full-width rank oracle at both1024-word base transitions,
        // including the largest possible16-bit relative rank.
        let mut words = vec![u64::MAX; 2051];
        words[1023] = 0;
        words[1024] = 5;
        words[2048] = 1;
        let mut rank = 17u64;
        let ranks: Vec<_> = words
            .chunks(8)
            .map(|chunk| {
                let base = rank;
                rank += chunk.iter().map(|v| u64::from(v.count_ones())).sum::<u64>();
                base
            })
            .collect();
        let mut bytes = Vec::new();
        for value in [words.len() as u64 * 64, words.len() as u64] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        for value in &words {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes.extend_from_slice(&(ranks.len() as u64).to_le_bytes());
        for value in &ranks {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        let mut compact_words = BitVector::load(&mut std::io::Cursor::new(bytes)).unwrap();
        let probes: Vec<_> = [0usize, 1, 1022, 1023, 1024, 1025, 2047, 2048, 2050]
            .into_iter()
            .flat_map(|word| [0usize, 1, 17, 63].map(move |offset| (word * 64 + offset) as u64))
            .collect();
        let expected: Vec<_> = probes
            .iter()
            .map(|&pos| compact_words.get(pos).then(|| compact_words.rank(pos)))
            .collect();
        assert_eq!(compact_words.rank_bases.len(), 3);
        compact_words.bit_array.clear();
        compact_words.ranks.clear();
        for (&pos, expected) in probes.iter().zip(expected) {
            assert_eq!(compact_words.get_rank(pos), expected);
        }
        for base in [u64::from(u32::MAX) - 1, u64::from(u32::MAX) + 1] {
            let bytes: Vec<u8> = [64u64, 1, 0b101, 1, base]
                .into_iter()
                .flat_map(u64::to_le_bytes)
                .collect();
            let bitset = BitVector::load(&mut std::io::Cursor::new(bytes)).unwrap();
            for pos in 0..64 {
                assert_eq!(
                    bitset.get_rank(pos),
                    bitset.get(pos).then(|| bitset.rank(pos))
                );
            }
        }

        let temp = tempfile::tempdir().unwrap();
        let fasta = temp.path().join("transcripts.fa");
        let mut state = 137u64;
        let mut text = String::from(">tx\n");
        for _ in 0..5000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            text.push(b"ACGT"[(state & 3) as usize] as char);
        }
        text.push('\n');
        std::fs::write(&fasta, text).unwrap();
        let path = temp.path().join("index.idx");
        crate::index::build_index(&path, &[fasta], crate::index::IndexBuildOptions::default())
            .unwrap();
        let data = std::fs::read(path).unwrap();
        let mut reader = std::io::Cursor::new(data.as_slice());
        let _version = read_u64_le(&mut reader).unwrap();
        let graph_bytes = read_u64_le(&mut reader).unwrap();
        reader.set_position(16 + graph_bytes);
        let _mphf_bytes = read_u64_le(&mut reader).unwrap();
        let mphf = BooPhf::load(&mut reader).unwrap();
        let mut compact = mphf.clone();
        compact.release_redundant_rank_storage();
        for (old, new) in mphf.levels.iter().zip(&compact.levels) {
            for pos in 0..old.hash_domain {
                assert_eq!(new.bitset.get_rank(pos), old.bitset.get_rank(pos));
                assert_eq!(new.bitset.get(pos), old.bitset.get(pos));
                assert_eq!(new.bitset.rank(pos), old.bitset.rank(pos));
            }
        }
        for level in &mphf.levels {
            for pos in 0..level.hash_domain {
                assert_eq!(
                    level.bitset.get_rank(pos),
                    level.bitset.get(pos).then(|| level.bitset.rank(pos))
                );
            }
        }
        assert!(mphf.compact_levels);
        let mut fallback = mphf.clone();
        fallback.compact_levels = false;
        for level in &mut fallback.levels {
            level.bitset.ranked_words.clear();
        }
        for batch in 0..2500 {
            let keys = std::array::from_fn(|lane| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (lane < batch % 9).then_some(state.to_le_bytes())
            });
            let expected = keys.map(|key| key.and_then(|key| mphf.debug_lookup(&key).0));
            assert_eq!(mphf.lookup_batch8(&keys), expected);
            assert_eq!(compact.lookup_batch8(&keys), expected);
            assert_eq!(fallback.lookup_batch8(&keys), expected);
            for key in keys.into_iter().flatten() {
                assert_eq!(mphf.lookup(&key), mphf.debug_lookup(&key).0);
                assert_eq!(compact.lookup(&key), mphf.debug_lookup(&key).0);
                assert_eq!(fallback.lookup(&key), mphf.debug_lookup(&key).0);
            }
        }
        let key = 1729u64.to_le_bytes();
        let mut tail = mphf;
        tail.nb_levels = 1;
        tail.levels.clear();
        tail.final_hash = HashMap::from([(key, 12)]);
        tail.lastbitsetrank = 7;
        let keys = [
            Some(key),
            None,
            Some([0; 8]),
            None,
            Some(key),
            None,
            None,
            None,
        ];
        assert_eq!(
            tail.lookup_batch8(&keys),
            [Some(19), None, None, None, Some(19), None, None, None]
        );
    }
}

#[cfg(test)]
mod exact_rep_hash_tests {
    #[test]
    fn fixed_minimizer_keys_equal_generic_bytes_for_all_bits_and_seeds() {
        let check = |word: u64| {
            let bytes = word.to_le_bytes();
            for seed in [0, 1, u64::MAX] {
                assert_eq!(
                    super::minimizer_hash(&bytes, seed),
                    super::wyhash(&bytes, 0) ^ seed
                );
            }
        };
        for bit in 0..64 {
            check(1u64 << bit);
            check(!(1u64 << bit));
        }
        for word in [0, u64::MAX, u32::MAX as u64, 1u64 << 32] {
            check(word);
        }
        let mut state = 59u64;
        for _ in 0..4096 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            check(state);
        }
    }

    #[test]
    fn fixed_words_equal_generic_byte_hash_in_both_orders() {
        let mut state = 17u64;
        for first in [0, 1, u32::MAX as u64, 1u64 << 32, u64::MAX] {
            for second in [0, 1, u64::MAX] {
                check(first, second);
            }
        }
        for _ in 0..4096 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let first = state;
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            check(first, state);
            check(state, first);
        }
    }
    fn check(first: u64, second: u64) {
        let mut bytes = [0u8; 16];
        bytes[..8].copy_from_slice(&first.to_le_bytes());
        bytes[8..].copy_from_slice(&second.to_le_bytes());
        assert_eq!(
            super::wyhash16_words(first, second),
            super::wyhash(&bytes, 0)
        );
    }
}
