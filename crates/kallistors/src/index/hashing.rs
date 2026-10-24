//! Exact native-index minimizer hashing and MPHF level state updates.

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

pub(super) fn xorshift_next(s0: &mut u64, s1: &mut u64) -> u64 {
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
