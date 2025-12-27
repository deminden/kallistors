//! Sequence-specific bias utilities.

use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Observed 5' hexamer bias counts.
#[derive(Debug, Clone)]
pub struct BiasCounts {
    pub counts: Vec<u32>,
    pub total: u64,
}

impl BiasCounts {
    pub fn new() -> Self {
        Self {
            counts: vec![0; 4096],
            total: 0,
        }
    }

    pub fn record(&mut self, hex: usize) {
        if hex >= self.counts.len() {
            return;
        }
        self.counts[hex] = self.counts[hex].saturating_add(1);
        self.total = self.total.saturating_add(1);
    }

    pub fn load_text(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
        let reader = BufReader::new(file);
        let mut counts = vec![0u32; 4096];
        let mut total = 0u64;
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut parts = line.splitn(2, '\t');
            let idx: usize = parts
                .next()
                .ok_or_else(|| {
                    Error::InvalidFormat(format!("missing hexamer at line {}", line_no + 1))
                })?
                .parse()
                .map_err(|_| {
                    Error::InvalidFormat(format!("invalid hexamer at line {}", line_no + 1))
                })?;
            let count: u32 = parts
                .next()
                .ok_or_else(|| {
                    Error::InvalidFormat(format!("missing count at line {}", line_no + 1))
                })?
                .parse()
                .map_err(|_| {
                    Error::InvalidFormat(format!("invalid count at line {}", line_no + 1))
                })?;
            if idx >= counts.len() {
                return Err(Error::InvalidFormat(format!(
                    "hexamer out of range at line {}",
                    line_no + 1
                )));
            }
            counts[idx] = count;
            total += count as u64;
        }
        Ok(Self { counts, total })
    }

    pub fn write_text(&self, path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "#hexamer\tcount")?;
        for (i, count) in self.counts.iter().enumerate() {
            writeln!(writer, "{}\t{}", i, count)?;
        }
        Ok(())
    }
}

impl Default for BiasCounts {
    fn default() -> Self {
        Self::new()
    }
}

/// Convert a 6-mer into its integer representation.
pub fn hexamer_to_int(seq: &[u8], revcomp: bool) -> Option<usize> {
    if seq.len() < 6 {
        return None;
    }
    let mut hex: usize = 0;
    if !revcomp {
        for &b in &seq[..6] {
            hex <<= 2;
            match b & 0xDF {
                b'A' => {}
                b'C' => hex |= 1,
                b'G' => hex |= 2,
                b'T' => hex |= 3,
                _ => return None,
            }
        }
    } else {
        for (i, &b) in seq[..6].iter().enumerate() {
            let val = match b & 0xDF {
                b'A' => 3,
                b'C' => 2,
                b'G' => 1,
                b'T' => 0,
                _ => return None,
            };
            hex |= val << (2 * i);
        }
    }
    Some(hex)
}

/// Update a hexamer rolling hash with the next base.
pub fn update_hexamer(mut hex: usize, base: u8, revcomp: bool) -> Option<usize> {
    if !revcomp {
        hex = (hex & 0x3ff) << 2;
        let val = match base & 0xDF {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        Some(hex | val)
    } else {
        hex >>= 2;
        let val = match base & 0xDF {
            b'A' => 3,
            b'C' => 2,
            b'G' => 1,
            b'T' => 0,
            _ => return None,
        };
        Some(hex | (val << 10))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hexamer_to_int_forward() {
        let seq = b"ACGTAC";
        let hex = hexamer_to_int(seq, false).unwrap();
        assert_eq!(hex, 0b00_01_10_11_00_01);
    }

    #[test]
    fn hexamer_to_int_revcomp() {
        let seq = b"ACGTAC";
        let hex = hexamer_to_int(seq, true).unwrap();
        assert_eq!(hex, 2843);
    }
}
