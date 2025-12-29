use std::collections::HashMap;
use std::io::Read;

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

#[derive(Clone)]
pub struct BooPhf {
    gamma: f64,
    nb_levels: i32,
    lastbitsetrank: u64,
    nelem: u64,
    levels: Vec<BooLevel>,
    final_hash: HashMap<[u8; 8], u64>,
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

#[derive(Clone)]
struct BooLevel {
    bitset: BitVector,
    hash_domain: u64,
}

#[derive(Clone)]
struct BitVector {
    _size: u64,
    _nchar: u64,
    bit_array: Vec<u64>,
    ranks: Vec<u64>,
}

impl BitVector {
    fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let size = read_u64_le(reader)?;
        let nchar = read_u64_le(reader)?;
        let mut bit_array = vec![0u64; nchar as usize];
        for slot in &mut bit_array {
            *slot = read_u64_le(reader)?;
        }
        let ranks_len = read_u64_le(reader)? as usize;
        let mut ranks = vec![0u64; ranks_len];
        for slot in &mut ranks {
            *slot = read_u64_le(reader)?;
        }
        Ok(Self {
            _size: size,
            _nchar: nchar,
            bit_array,
            ranks,
        })
    }

    fn get(&self, pos: u64) -> bool {
        let idx = (pos >> 6) as usize;
        let bit = pos & 63;
        ((self.bit_array[idx] >> bit) & 1) != 0
    }

    fn rank(&self, pos: u64) -> u64 {
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
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let gamma = read_f64_le(reader)?;
        let nb_levels = read_i32_le(reader)?;
        let lastbitsetrank = read_u64_le(reader)?;
        let nelem = read_u64_le(reader)?;
        let mut levels = Vec::with_capacity(nb_levels as usize);
        for _ in 0..nb_levels {
            let bitset = BitVector::load(reader)?;
            levels.push(BooLevel {
                bitset,
                hash_domain: 0,
            });
        }
        let final_hash_size = read_u64_le(reader)? as usize;
        let mut final_hash = HashMap::with_capacity(final_hash_size);
        for _ in 0..final_hash_size {
            let key = read_u64_le(reader)?;
            let val = read_u64_le(reader)?;
            final_hash.insert(key.to_le_bytes(), val);
        }
        let mut phf = Self {
            gamma,
            nb_levels,
            lastbitsetrank,
            nelem,
            levels,
            final_hash,
        };
        phf.recompute_hash_domains();
        Ok(phf)
    }

    fn recompute_hash_domains(&mut self) {
        let nelem = self.nelem as f64;
        let gamma = self.gamma;
        let proba_collision = 1.0 - (((gamma * nelem - 1.0) / (gamma * nelem)).powf(nelem - 1.0));
        let hash_domain = (gamma * nelem).ceil() as u64;
        for (i, level) in self.levels.iter_mut().enumerate() {
            let raw = (hash_domain as f64 * proba_collision.powi(i as i32)).ceil() as u64;
            let mut domain = raw.div_ceil(64) * 64;
            if domain == 0 {
                domain = 64;
            }
            level.hash_domain = domain;
        }
    }

    pub fn lookup(&self, key: &[u8; 8]) -> Option<u64> {
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
            let v = self.final_hash.get(key)?;
            return Some(v + self.lastbitsetrank);
        }

        let level = level as usize;
        let non_min = hash_raw % self.levels[level].hash_domain;
        Some(self.levels[level].bitset.rank(non_min))
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

pub fn minimizer_hash(key: &[u8; 8], seed: u64) -> u64 {
    wyhash(key, 0) ^ seed
}

fn encode_word_pair(seq: &[u8]) -> Option<(u64, u64)> {
    let mut fwd = 0u64;
    for (i, &b) in seq.iter().enumerate() {
        let v = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        let shift = 62 - ((i & 0x1f) << 1);
        fwd |= v << shift;
    }
    let mut rev = 0u64;
    for (i, &b) in seq.iter().rev().enumerate() {
        let v = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        let comp = 3 - v;
        let shift = 62 - ((i & 0x1f) << 1);
        rev |= comp << shift;
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
