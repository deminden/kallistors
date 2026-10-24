//! Shared native-index binary primitives with explicit little-endian encoding.

use crate::Result;
use std::io::{Read, Write};

pub(super) fn read_u64_array<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u64>> {
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

pub(super) fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

pub(super) fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

pub(super) fn read_u16_le<R: Read>(reader: &mut R) -> Result<u16> {
    let mut buf = [0u8; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

pub(super) fn read_u32_le<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(super) fn read_i32_le<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

pub(super) fn read_f64_le<R: Read>(reader: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

pub(super) fn write_u8<W: Write>(writer: &mut W, v: u8) -> Result<()> {
    writer.write_all(&[v])?;
    Ok(())
}

pub(super) fn write_u32<W: Write>(writer: &mut W, v: u32) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

pub(super) fn write_i32<W: Write>(writer: &mut W, v: i32) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

pub(super) fn write_u64<W: Write>(writer: &mut W, v: u64) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}

pub(super) fn write_f64<W: Write>(writer: &mut W, v: f64) -> Result<()> {
    writer.write_all(&v.to_le_bytes())?;
    Ok(())
}
