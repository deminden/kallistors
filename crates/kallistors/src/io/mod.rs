//! I/O utilities and FASTQ parsing.

use std::io::BufRead;

use crate::Result;
use crate::error::Error;

/// A FASTQ record (raw bytes, no validation).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub header: Vec<u8>,
    pub seq: Vec<u8>,
    pub plus: Vec<u8>,
    pub qual: Vec<u8>,
}

/// Streaming source for FASTQ records.
pub trait ReadSource {
    fn next_record(&mut self) -> Option<Result<FastqRecord>>;
}

/// A buffered FASTQ reader for plain text input.
pub struct FastqReader<R: BufRead> {
    reader: R,
    buf: Vec<u8>,
}

impl<R: BufRead> FastqReader<R> {
    /// Create a new FASTQ reader.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buf: Vec::with_capacity(256),
        }
    }

    fn read_line(&mut self) -> Result<Option<Vec<u8>>> {
        self.buf.clear();
        let bytes = self.reader.read_until(b'\n', &mut self.buf)?;
        if bytes == 0 {
            return Ok(None);
        }
        if self.buf.ends_with(b"\n") {
            self.buf.pop();
            if self.buf.ends_with(b"\r") {
                self.buf.pop();
            }
        }
        Ok(Some(self.buf.clone()))
    }
}

impl<R: BufRead> ReadSource for FastqReader<R> {
    fn next_record(&mut self) -> Option<Result<FastqRecord>> {
        let header = match self.read_line() {
            Ok(Some(line)) => line,
            Ok(None) => return None,
            Err(err) => return Some(Err(err)),
        };
        let seq = match self.read_line() {
            Ok(Some(line)) => line,
            Ok(None) => return Some(Err(Error::InvalidFormat("truncated FASTQ".into()))),
            Err(err) => return Some(Err(err)),
        };
        let plus = match self.read_line() {
            Ok(Some(line)) => line,
            Ok(None) => return Some(Err(Error::InvalidFormat("truncated FASTQ".into()))),
            Err(err) => return Some(Err(err)),
        };
        let qual = match self.read_line() {
            Ok(Some(line)) => line,
            Ok(None) => return Some(Err(Error::InvalidFormat("truncated FASTQ".into()))),
            Err(err) => return Some(Err(err)),
        };

        Some(Ok(FastqRecord {
            header,
            seq,
            plus,
            qual,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_single_record() {
        let data = b"@r1\nACGT\n+\n!!!!\n";
        let mut reader = FastqReader::new(Cursor::new(data));
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.header, b"@r1");
        assert_eq!(record.seq, b"ACGT");
        assert_eq!(record.plus, b"+");
        assert_eq!(record.qual, b"!!!!");
        assert!(reader.next_record().is_none());
    }
}
