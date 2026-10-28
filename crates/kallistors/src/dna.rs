//! Shared exact DNA word encoding for index construction and read matching.

pub(crate) const INVALID_BASE_CODE: u8 = 4;
pub(crate) const BASE_CODES: [u8; 256] = {
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

#[inline]
fn zero_byte_high_bits(word: u64) -> u64 {
    !((word & 0x7f7f_7f7f_7f7f_7f7f).wrapping_add(0x7f7f_7f7f_7f7f_7f7f) | word)
        & 0x8080_8080_8080_8080
}

#[inline]
fn gather_base_pairs(mut word: u64) -> u64 {
    word = (word | (word >> 6)) & 0x000f_000f_000f_000f;
    word = (word | (word >> 12)) & 0x0000_00ff_0000_00ff;
    (word | (word >> 24)) & 0xffff
}

#[cfg(all(
    target_arch = "aarch64",
    target_feature = "neon",
    target_endian = "little"
))]
#[inline]
pub(crate) fn encode_sixteen_bases(bytes: &[u8; 16]) -> Option<u64> {
    use std::arch::aarch64::*;
    // SAFETY: the target enables NEON; each load covers a complete array.
    unsafe {
        let input = vld1q_u8(bytes.as_ptr());
        let upper = vandq_u8(input, vdupq_n_u8(0xdf));
        let valid = vorrq_u8(
            vorrq_u8(
                vceqq_u8(upper, vdupq_n_u8(b'A')),
                vceqq_u8(upper, vdupq_n_u8(b'C')),
            ),
            vorrq_u8(
                vceqq_u8(upper, vdupq_n_u8(b'G')),
                vceqq_u8(upper, vdupq_n_u8(b'T')),
            ),
        );
        if vminvq_u8(valid) != 255 {
            return None;
        }
        let codes = vandq_u8(
            veorq_u8(vshrq_n_u8::<1>(input), vshrq_n_u8::<2>(input)),
            vdupq_n_u8(3),
        );
        const WEIGHTS: [u8; 16] = [64, 16, 4, 1, 64, 16, 4, 1, 64, 16, 4, 1, 64, 16, 4, 1];
        let weighted = vmulq_u8(codes, vld1q_u8(WEIGHTS.as_ptr()));
        let groups = vpaddlq_u16(vpaddlq_u8(weighted));
        let narrow = vmovn_u16(vcombine_u16(vmovn_u32(groups), vdup_n_u16(0)));
        Some(u64::from(
            vget_lane_u32::<0>(vreinterpret_u32_u8(narrow)).swap_bytes(),
        ))
    }
}

#[inline]
pub(crate) fn reverse_complement_word(forward: u64, len: usize) -> u64 {
    if len == 0 {
        return 0;
    }
    let reversed = forward.reverse_bits();
    // Reversing bits also swaps the two bits within each base. Undo that swap
    // before complementing (A/C/G/T codes complement by XOR with 3).
    let complement =
        !(((reversed & 0x5555_5555_5555_5555) << 1) | ((reversed & 0xaaaa_aaaa_aaaa_aaaa) >> 1));
    complement >> ((32 - len) * 2)
}

// Both words are right-aligned, with the first base in the highest used pair
// of bits. Validate every byte before using the ASCII-to-two-bit shortcut.
pub(crate) fn encode_kmer_pair(seq: &[u8]) -> Option<(u64, u64)> {
    if seq.len() > 32 {
        return None;
    }
    #[cfg(all(
        target_arch = "aarch64",
        target_feature = "neon",
        target_endian = "little"
    ))]
    if seq.len() >= 16 {
        let first = encode_sixteen_bases(seq[..16].try_into().expect("sixteen bases"))?;
        let remaining = seq.len() - 16;
        let forward = if remaining == 0 {
            first
        } else {
            let last =
                encode_sixteen_bases(seq[seq.len() - 16..].try_into().expect("sixteen bases"))?;
            (first << (2 * remaining)) | (last & ((1u64 << (2 * remaining)) - 1))
        };
        return Some((forward, reverse_complement_word(forward, seq.len())));
    }
    let mut forward = 0u64;
    let (chunks, tail) = seq.as_chunks::<8>();
    for &bytes in chunks {
        let word = u64::from_be_bytes(bytes);
        let upper = word & 0xdfdf_dfdf_dfdf_dfdf;
        let valid = zero_byte_high_bits(upper ^ 0x4141_4141_4141_4141)
            | zero_byte_high_bits(upper ^ 0x4343_4343_4343_4343)
            | zero_byte_high_bits(upper ^ 0x4747_4747_4747_4747)
            | zero_byte_high_bits(upper ^ 0x5454_5454_5454_5454);
        if valid != 0x8080_8080_8080_8080 {
            return None;
        }
        let codes = ((word >> 1) ^ (word >> 2)) & 0x0303_0303_0303_0303;
        forward = (forward << 16) | gather_base_pairs(codes);
    }
    for &base in tail {
        let code = BASE_CODES[base as usize];
        if code == INVALID_BASE_CODE {
            return None;
        }
        forward = (forward << 2) | u64::from(code);
    }
    Some((forward, reverse_complement_word(forward, seq.len())))
}
