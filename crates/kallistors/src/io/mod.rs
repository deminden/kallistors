//! I/O utilities and FASTQ parsing.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::Result;
use crate::error::Error;

const FASTQ_IO_BUFFER_BYTES: usize = 1 << 20;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PackedSpan {
    start: u32,
    len: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PackedFastqRecord {
    header: PackedSpan,
    seq: PackedSpan,
    plus: PackedSpan,
    qual: PackedSpan,
}

#[derive(Debug, Clone)]
pub struct PackedFastqBatch {
    storage: Vec<u8>,
    records: Vec<PackedFastqRecord>,
}

#[derive(Debug, Clone, Copy)]
pub struct PackedFastqRecordRef<'a> {
    pub header: &'a [u8],
    pub seq: &'a [u8],
    pub plus: &'a [u8],
    pub qual: &'a [u8],
}

pub trait PackedBatchSource {
    fn next_packed_batch(&mut self, batch_size: usize) -> Result<Option<PackedFastqBatch>>;
}

impl PackedFastqBatch {
    fn with_capacity(records: usize, bytes: usize) -> Self {
        Self {
            storage: Vec::with_capacity(bytes),
            records: Vec::with_capacity(records),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn records(&self) -> impl Iterator<Item = PackedFastqRecordRef<'_>> + '_ {
        self.records.iter().map(|record| PackedFastqRecordRef {
            header: self.slice(record.header),
            seq: self.slice(record.seq),
            plus: self.slice(record.plus),
            qual: self.slice(record.qual),
        })
    }

    fn slice(&self, span: PackedSpan) -> &[u8] {
        let start = span.start as usize;
        let end = start + span.len as usize;
        &self.storage[start..end]
    }
}

pub enum FastqInput {
    Gzip(Box<BufReader<flate2::read::MultiGzDecoder<BufReader<File>>>>),
    Plain(BufReader<File>),
}

impl Read for FastqInput {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Gzip(reader) => reader.read(buf),
            Self::Plain(reader) => reader.read(buf),
        }
    }
}

impl BufRead for FastqInput {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            Self::Gzip(reader) => reader.fill_buf(),
            Self::Plain(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::Gzip(reader) => reader.consume(amt),
            Self::Plain(reader) => reader.consume(amt),
        }
    }
}

/// Open a FASTQ reader from a plain or gzip-compressed file.
pub fn open_fastq_reader(path: &Path) -> Result<FastqReader<FastqInput>> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let reader = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        let file = BufReader::with_capacity(FASTQ_IO_BUFFER_BYTES, file);
        let decoder = flate2::read::MultiGzDecoder::new(file);
        FastqInput::Gzip(Box::new(BufReader::with_capacity(
            FASTQ_IO_BUFFER_BYTES,
            decoder,
        )))
    } else {
        FastqInput::Plain(BufReader::with_capacity(FASTQ_IO_BUFFER_BYTES, file))
    };
    Ok(FastqReader::new(reader))
}

/// A buffered FASTQ reader for plain text input.
pub struct FastqReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    /// Create a new FASTQ reader.
    pub fn new(reader: R) -> Self {
        Self { reader }
    }

    fn read_line_into(&mut self, out: &mut Vec<u8>) -> Result<bool> {
        out.clear();
        let bytes = self.reader.read_until(b'\n', out)?;
        if bytes == 0 {
            return Ok(false);
        }
        if out.ends_with(b"\n") {
            out.pop();
            if out.ends_with(b"\r") {
                out.pop();
            }
        }
        Ok(true)
    }

    fn read_line_span_into(&mut self, storage: &mut Vec<u8>) -> Result<Option<PackedSpan>> {
        let start = storage.len();
        let bytes = self.reader.read_until(b'\n', storage)?;
        if bytes == 0 {
            return Ok(None);
        }
        if storage.ends_with(b"\n") {
            storage.pop();
            if storage.ends_with(b"\r") {
                storage.pop();
            }
        }
        let len = storage.len() - start;
        Ok(Some(PackedSpan {
            start: start as u32,
            len: len as u32,
        }))
    }
}

impl<R: BufRead> ReadSource for FastqReader<R> {
    fn next_record(&mut self) -> Option<Result<FastqRecord>> {
        let mut header = Vec::with_capacity(128);
        let mut seq = Vec::with_capacity(256);
        let mut plus = Vec::with_capacity(8);
        let mut qual = Vec::with_capacity(256);

        match self.read_line_into(&mut header) {
            Ok(true) => {}
            Ok(false) => return None,
            Err(err) => return Some(Err(err)),
        };
        match self.read_line_into(&mut seq) {
            Ok(true) => {}
            Ok(false) => return Some(Err(Error::InvalidFormat("truncated FASTQ".into()))),
            Err(err) => return Some(Err(err)),
        };
        match self.read_line_into(&mut plus) {
            Ok(true) => {}
            Ok(false) => return Some(Err(Error::InvalidFormat("truncated FASTQ".into()))),
            Err(err) => return Some(Err(err)),
        };
        match self.read_line_into(&mut qual) {
            Ok(true) => {}
            Ok(false) => return Some(Err(Error::InvalidFormat("truncated FASTQ".into()))),
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

impl<R: BufRead> PackedBatchSource for FastqReader<R> {
    fn next_packed_batch(&mut self, batch_size: usize) -> Result<Option<PackedFastqBatch>> {
        let mut batch = PackedFastqBatch::with_capacity(batch_size, batch_size * 700);
        while batch.records.len() < batch_size {
            let Some(header) = self.read_line_span_into(&mut batch.storage)? else {
                break;
            };
            let Some(seq) = self.read_line_span_into(&mut batch.storage)? else {
                return Err(Error::InvalidFormat("truncated FASTQ".into()));
            };
            let Some(plus) = self.read_line_span_into(&mut batch.storage)? else {
                return Err(Error::InvalidFormat("truncated FASTQ".into()));
            };
            let Some(qual) = self.read_line_span_into(&mut batch.storage)? else {
                return Err(Error::InvalidFormat("truncated FASTQ".into()));
            };
            batch.records.push(PackedFastqRecord {
                header,
                seq,
                plus,
                qual,
            });
        }
        if batch.is_empty() {
            Ok(None)
        } else {
            Ok(Some(batch))
        }
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
