//! Lossless normalized transcript scratch, read directly during graph mapping.
use crate::{Error, Result};
use std::sync::Arc;

pub(super) struct PackedSequence {
    len: usize,
    bytes: Option<Arc<Vec<u8>>>,
    offset: usize,
}
impl PackedSequence {
    #[cfg(test)]
    pub(super) fn new(seq: &[u8]) -> Result<Self> {
        Self::group(std::iter::once(seq)).map(|mut group| group.remove(0))
    }
    pub(super) fn group<'a>(
        sequences: impl ExactSizeIterator<Item = &'a [u8]> + Clone,
    ) -> Result<Vec<Self>> {
        // Large shared groups separate packed scratch from the small ASCII
        // allocations being retired, and avoid one allocator object per sequence.
        let mut bytes =
            Vec::with_capacity(sequences.clone().map(|seq| seq.len().div_ceil(4)).sum());
        let mut ranges = Vec::with_capacity(sequences.len());
        for seq in sequences {
            ranges.push((bytes.len(), seq.len()));
            Self::append(&mut bytes, seq)?;
        }
        let bytes = Arc::new(bytes);
        Ok(ranges
            .into_iter()
            .map(|(offset, len)| Self {
                len,
                offset,
                bytes: Some(Arc::clone(&bytes)),
            })
            .collect())
    }
    fn append(bytes: &mut Vec<u8>, seq: &[u8]) -> Result<()> {
        #[cfg(all(
            target_arch = "aarch64",
            target_feature = "neon",
            target_endian = "little"
        ))]
        let tail = {
            let (chunks, tail) = seq.as_chunks::<16>();
            for chunk in chunks {
                let word = crate::dna::encode_sixteen_bases(chunk)
                    .ok_or_else(|| Error::InvalidFormat("invalid DNA base".into()))?;
                bytes.extend_from_slice(&(word as u32).to_be_bytes());
            }
            tail
        };
        #[cfg(not(all(
            target_arch = "aarch64",
            target_feature = "neon",
            target_endian = "little"
        )))]
        let tail = {
            let (chunks, tail) = seq.as_chunks::<32>();
            for chunk in chunks {
                let (word, _) = crate::dna::encode_kmer_pair(chunk)
                    .ok_or_else(|| Error::InvalidFormat("invalid DNA base".into()))?;
                bytes.extend_from_slice(&word.to_be_bytes());
            }
            tail
        };
        for chunk in tail.chunks(4) {
            let (word, _) = crate::dna::encode_kmer_pair(chunk)
                .ok_or_else(|| Error::InvalidFormat("invalid DNA base".into()))?;
            bytes.push((word << (2 * (4 - chunk.len()))) as u8);
        }
        Ok(())
    }
    pub(super) fn len(&self) -> usize {
        self.len
    }
    #[inline]
    pub(super) fn code(&self, index: usize) -> u64 {
        u64::from(
            (self.bytes.as_ref().expect("live packed sequence")[self.offset + (index >> 2)]
                >> (6 - 2 * (index & 3)))
                & 3,
        )
    }
    #[cfg(test)]
    pub(super) fn base(&self, index: usize) -> u8 {
        b"ACGT"[self.code(index) as usize]
    }
    #[inline]
    pub(super) fn word(&self, position: usize, k: usize) -> u64 {
        debug_assert!(k > 0 && k <= 32 && position + k <= self.len);
        let bytes = &self.bytes.as_ref().expect("live packed sequence")
            [self.offset + (position >> 2)..self.offset + self.len.div_ceil(4)];
        let mut word = if bytes.len() >= 8 {
            u64::from_be_bytes(bytes[..8].try_into().expect("eight packed bytes"))
        } else {
            let mut padded = [0; 8];
            padded[..bytes.len()].copy_from_slice(bytes);
            u64::from_be_bytes(padded)
        };
        let shift = (position & 3) * 2;
        word <<= shift;
        if shift != 0 && bytes.len() > 8 {
            word |= u64::from(bytes[8]) >> (8 - shift);
        }
        word & (!0u64 << (64 - 2 * k))
    }
    pub(super) fn clear(&mut self) {
        self.bytes = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn packed_words_match_byte_encoder_at_every_offset_and_tail() {
        let mut state = 42u64;
        for len in 0..129 {
            let seq = (0..len)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect::<Vec<_>>();
            let packed = PackedSequence::new(&seq).unwrap();
            for (i, &base) in seq.iter().enumerate() {
                assert_eq!(packed.base(i), base);
            }
            for k in (3..=31).step_by(2).filter(|&k| k <= len) {
                for start in 0..=len - k {
                    let (forward, _) =
                        crate::dna::encode_kmer_pair(&seq[start..start + k]).unwrap();
                    assert_eq!(
                        packed.word(start, k),
                        forward << (64 - 2 * k),
                        "len={len},start={start},k={k}"
                    );
                }
            }
        }
        assert!(PackedSequence::new(b"ACNT").is_err());
        let input = [
            b"ACG".as_slice(),
            b"T",
            b"",
            b"TGCAACGTTGCATGCATGCATGCATGCATGCATGCA",
        ];
        let mut grouped = PackedSequence::group(input.into_iter()).unwrap();
        grouped[1].clear();
        for i in [0, 2, 3] {
            for (j, &base) in input[i].iter().enumerate() {
                assert_eq!(grouped[i].base(j), base);
            }
            for k in (3..=31).step_by(2).filter(|&k| k <= input[i].len()) {
                for start in 0..=input[i].len() - k {
                    let (word, _) =
                        crate::dna::encode_kmer_pair(&input[i][start..start + k]).unwrap();
                    assert_eq!(grouped[i].word(start, k), word << (64 - 2 * k));
                }
            }
        }
    }
}
