//! Byte-exact read-cache keys with a lossless uppercase DNA representation.
// The common packed cDNA key fits in one cache entry. Raw/non-DNA and
// longer keys retain their exact byte representation and allocation reuse.
pub(super) enum StoredKey {
    Inline { header: [u8; 2], bytes: [u8; 24] },
    Heap(Box<[u8]>),
}

pub(super) struct QueryKey<'a> {
    bytes: &'a [u8],
    header: [u8; 2],
}

impl<'a> QueryKey<'a> {
    pub(super) fn new(sequence: &'a [u8], scratch: &'a mut [u8; 64]) -> Self {
        let packed = (sequence.len() <= 256)
            .then(|| pack(sequence, scratch))
            .flatten();
        let header = (sequence.len() as u16) | if packed.is_some() { 1 << 15 } else { 0 };
        Self {
            bytes: packed.map_or(sequence, |len| &scratch[..len]),
            header: header.to_le_bytes(),
        }
    }

    pub(super) fn hash(&self) -> usize {
        let mut hash = u64::from(u16::from_le_bytes(self.header));
        let (words, tail) = self.bytes.as_chunks::<8>();
        for word in words {
            hash = (hash ^ u64::from_le_bytes(*word))
                .wrapping_mul(0x9e3779b97f4a7c15)
                .rotate_left(27);
        }
        for &base in tail {
            hash = (hash ^ u64::from(base)).wrapping_mul(0x100000001b3);
        }
        (hash ^ (hash >> 32)) as usize
    }

    pub(super) fn matches(&self, stored: &StoredKey) -> bool {
        match stored {
            StoredKey::Inline { header, bytes } => {
                self.header == *header && bytes.get(..self.bytes.len()) == Some(self.bytes)
            }
            StoredKey::Heap(bytes) => self.matches_bytes(bytes),
        }
    }

    fn matches_bytes(&self, stored: &[u8]) -> bool {
        stored.get(..2) == Some(self.header.as_slice()) && stored.get(2..) == Some(self.bytes)
    }

    pub(super) fn stored(&self) -> StoredKey {
        if self.bytes.len() <= 24 {
            let mut bytes = [0; 24];
            bytes[..self.bytes.len()].copy_from_slice(self.bytes);
            StoredKey::Inline {
                header: self.header,
                bytes,
            }
        } else {
            StoredKey::Heap(self.boxed())
        }
    }

    pub(super) fn boxed(&self) -> Box<[u8]> {
        let mut value = Vec::with_capacity(2 + self.bytes.len());
        value.extend_from_slice(&self.header);
        value.extend_from_slice(self.bytes);
        value.into_boxed_slice()
    }

    pub(super) fn replace(&self, stored: &mut StoredKey) {
        match stored {
            StoredKey::Inline { header, bytes } if self.bytes.len() <= 24 => {
                *header = self.header;
                bytes[..self.bytes.len()].copy_from_slice(self.bytes);
            }
            StoredKey::Heap(bytes)
                if self.bytes.len() > 24 && bytes.len() == 2 + self.bytes.len() =>
            {
                bytes[..2].copy_from_slice(&self.header);
                bytes[2..].copy_from_slice(self.bytes);
            }
            _ => *stored = self.stored(),
        }
    }
}

fn pack(sequence: &[u8], output: &mut [u8; 64]) -> Option<usize> {
    let (blocks, tail) = sequence.as_chunks::<16>();
    for (number, block) in blocks.iter().enumerate() {
        output[number * 4..number * 4 + 4].copy_from_slice(&pack_sixteen(block)?.to_be_bytes());
    }
    let mut value = 0u8;
    let start = blocks.len() * 4;
    for (number, &byte) in tail.iter().enumerate() {
        let code = match byte {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        value |= code << (6 - 2 * (number & 3));
        if number & 3 == 3 {
            output[start + number / 4] = value;
            value = 0;
        }
    }
    if tail.len() & 3 != 0 {
        output[start + tail.len() / 4] = value;
    }
    Some(sequence.len().div_ceil(4))
}

#[inline]
fn pack_sixteen(bytes: &[u8; 16]) -> Option<u32> {
    #[cfg(all(target_arch = "aarch64", target_endian = "little"))]
    {
        use std::arch::aarch64::*;
        // SAFETY: both loads cover complete sixteen-byte arrays. NEON is
        // baseline on AArch64; no alignment or optional ISA is required.
        unsafe {
            let input = vld1q_u8(bytes.as_ptr());
            let valid = vorrq_u8(
                vorrq_u8(
                    vceqq_u8(input, vdupq_n_u8(b'A')),
                    vceqq_u8(input, vdupq_n_u8(b'C')),
                ),
                vorrq_u8(
                    vceqq_u8(input, vdupq_n_u8(b'G')),
                    vceqq_u8(input, vdupq_n_u8(b'T')),
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
            Some(u32::from_be_bytes(
                vget_lane_u32::<0>(vreinterpret_u32_u8(narrow)).to_le_bytes(),
            ))
        }
    }
    #[cfg(not(all(target_arch = "aarch64", target_endian = "little")))]
    {
        let mut value = 0u32;
        for &byte in bytes {
            let code = match byte {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            value = (value << 2) | code;
        }
        Some(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_keys_decode_to_exact_original_bytes_at_every_boundary() {
        for len in [
            0, 1, 2, 3, 4, 7, 8, 15, 16, 17, 23, 24, 25, 31, 32, 33, 91, 95, 96, 97, 99, 100, 101,
            255, 256,
        ] {
            let sequence: Vec<_> = (0..len).map(|i| b"ACGT"[(i * 7 + i / 5) & 3]).collect();
            for offset in [0, len / 2, len.saturating_sub(1)] {
                for byte in 0..=255u8 {
                    let mut original = sequence.clone();
                    if len != 0 {
                        original[offset] = byte;
                    }
                    let mut scratch = [0; 64];
                    let key = QueryKey::new(&original, &mut scratch);
                    let stored = key.boxed();
                    assert!(key.matches_bytes(&stored));
                    let mut inline_or_heap = key.stored();
                    assert!(key.matches(&inline_or_heap));
                    // Every byte and boundary must survive both representation
                    // transitions, including reuse of a longer inline payload.
                    let mut large_scratch = [0; 64];
                    let large = [b'N'; 256];
                    QueryKey::new(&large, &mut large_scratch).replace(&mut inline_or_heap);
                    assert!(QueryKey::new(&large, &mut large_scratch).matches(&inline_or_heap));
                    key.replace(&mut inline_or_heap);
                    assert!(key.matches(&inline_or_heap));
                    let header = u16::from_le_bytes(stored[..2].try_into().unwrap());
                    assert_eq!(usize::from(header & 0x7fff), len);
                    let decoded: Vec<_> = if header & 0x8000 != 0 {
                        (0..len)
                            .map(|i| {
                                b"ACGT"[((stored[2 + i / 4] >> (6 - 2 * (i & 3))) & 3) as usize]
                            })
                            .collect()
                    } else {
                        stored[2..].to_vec()
                    };
                    assert_eq!(decoded, original);
                    let other = if len == 0 {
                        vec![b'A']
                    } else {
                        original[..len - 1].to_vec()
                    };
                    let mut other_scratch = [0; 64];
                    let other_key = QueryKey::new(&other, &mut other_scratch);
                    assert!(!other_key.matches_bytes(&stored));
                    assert!(!other_key.matches(&inline_or_heap));
                    other_key.replace(&mut inline_or_heap);
                    assert!(!key.matches(&inline_or_heap));
                    key.replace(&mut inline_or_heap);
                    assert!(key.matches(&inline_or_heap));
                }
            }
        }
    }
}
