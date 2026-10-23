//! Streaming compressed sequence decoding and packed-word conversion.

use super::*;

// Decode a complete buffered graph record in place. Records crossing the
// input buffer retain Read::read_exact truncation and tail behavior.
pub(super) fn append_buffered_compressed<R: BufRead>(
    reader: &mut R,
    len: usize,
    out: &mut Vec<u64>,
    scratch: &mut Vec<u8>,
) -> Result<()> {
    let bytes = len.div_ceil(4);
    let available = reader.fill_buf()?;
    if available.len() >= bytes {
        append_compressed_words(&available[..bytes], len, out);
        reader.consume(bytes);
    } else {
        scratch.resize(bytes, 0);
        reader.read_exact(scratch)?;
        append_compressed_words(scratch, len, out);
    }
    Ok(())
}

// Bifrost stores the first base in the low pair of each byte; packed
// matchers place it in the high pair. Convert directly without ASCII buffers.
pub(super) fn append_compressed_words(data: &[u8], len: usize, out: &mut Vec<u64>) {
    const REVERSE_PAIRS: [u8; 256] = {
        let mut values = [0; 256];
        let mut i = 0;
        while i < 256 {
            values[i] =
                (((i & 3) << 6) | ((i & 12) << 2) | ((i & 48) >> 2) | ((i & 192) >> 6)) as u8;
            i += 1;
        }
        values
    };
    let (words, tail) = data.as_chunks::<8>();
    out.reserve(words.len() + usize::from(!tail.is_empty()));
    for bytes in words {
        // Reverse base-pair order with one word reversal, retaining bit order
        // within each two-bit base. Bifrost's byte order is little endian.
        let reversed = u64::from_le_bytes(*bytes).reverse_bits();
        out.push(
            ((reversed & 0x5555_5555_5555_5555) << 1) | ((reversed >> 1) & 0x5555_5555_5555_5555),
        );
    }
    if !tail.is_empty() {
        let mut bytes = [0; 8];
        for (slot, &byte) in bytes.iter_mut().zip(tail) {
            *slot = REVERSE_PAIRS[byte as usize];
        }
        out.push(u64::from_be_bytes(bytes));
    }
    // Padding belongs to the file format and is outside the sequence length.
    if !len.is_multiple_of(32)
        && let Some(last) = out.last_mut()
    {
        *last &= u64::MAX << ((32 - len % 32) * 2);
    }
}

pub(super) fn decode_compressed_sequence(data: &[u8], len: usize) -> Vec<u8> {
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

pub(super) fn decode_kmer(bytes: &[u8], k: usize) -> Vec<u8> {
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
