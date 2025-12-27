//! I/O utilities and FASTQ parsing.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

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

/// Open a FASTQ reader from a plain or gzip-compressed file.
pub fn open_fastq_reader(path: &Path) -> Result<FastqReader<Box<dyn BufRead>>> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let reader: Box<dyn BufRead> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        let decoder = flate2::read::MultiGzDecoder::new(file);
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(BufReader::new(file))
    };
    Ok(FastqReader::new(reader))
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

/// ReadSource over an owned in-memory batch of records.
pub struct VecReadSource {
    records: std::vec::IntoIter<FastqRecord>,
}

impl VecReadSource {
    pub fn new(records: Vec<FastqRecord>) -> Self {
        Self {
            records: records.into_iter(),
        }
    }
}

impl ReadSource for VecReadSource {
    fn next_record(&mut self) -> Option<Result<FastqRecord>> {
        self.records.next().map(Ok)
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

    #[test]
    fn open_fastq_reader_gz() {
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("reads.fq.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            File::create(&path).expect("create gz"),
            flate2::Compression::default(),
        );
        encoder
            .write_all(b"@r1\nACGT\n+\n!!!!\n")
            .expect("write gz");
        encoder.finish().expect("finish gz");

        let mut reader = open_fastq_reader(&path).expect("open gz reader");
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.header, b"@r1");
        assert_eq!(record.seq, b"ACGT");
        assert!(reader.next_record().is_none());
    }
}
