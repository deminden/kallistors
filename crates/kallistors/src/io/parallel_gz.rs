//! Parallel decoding of gzip files whose members record their own size.
//!
//! Block-compressed FASTQ writers (fasterp's `FP` subfield, BGZF's `BC`) emit
//! many small gzip members and store each member's compressed length in its
//! header. A splitter thread cuts the file into whole members without inflating
//! them, a few workers inflate members concurrently, and the reader hands the
//! output back in file order. Every stage is connected by a bounded channel,
//! so a slow consumer stops the splitter instead of growing memory.
//!
//! The output is byte-identical to sequential multi-member decoding: each
//! member is inflated and CRC-checked from its exact start, and members are
//! concatenated in order. A member whose stated size disagrees with what the
//! decoder consumes is an error. If a later member has no size field, the
//! splitter decodes the rest of the file sequentially itself.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, Sender, bounded};

/// Inflate workers per file. Two keep a paired-end run's producers ahead of
/// eight mapping workers on the measured hardware.
const WORKERS: usize = 2;
/// Members in flight per worker, bounding memory to a few MiB per file.
const QUEUE_PER_WORKER: usize = 4;
/// Output chunk size for the sequential fallback.
const FALLBACK_CHUNK: usize = 1 << 20;

/// Fixed gzip header bytes up to and including XLEN.
const HEADER_WITH_XLEN: usize = 12;

enum Job {
    /// One complete gzip member, in file order.
    Member(u64, Vec<u8>),
}

type Output = (u64, std::io::Result<Vec<u8>>);

/// Compressed length of a gzip member, read from its FEXTRA subfields:
/// fasterp's `FP` (four-byte total member size) or BGZF's `BC` (BSIZE, the
/// total size minus one). `None` when the header carries neither.
pub(crate) fn member_size(header: &[u8], extra: &[u8]) -> Option<usize> {
    if header.len() < HEADER_WITH_XLEN
        || header[0] != 0x1f
        || header[1] != 0x8b
        || header[2] != 8
        || header[3] & 4 == 0
    {
        return None;
    }
    let mut fields = extra;
    while fields.len() >= 4 {
        let len = u16::from_le_bytes([fields[2], fields[3]]) as usize;
        let data = fields.get(4..4 + len)?;
        match (&fields[..2], len) {
            (b"FP", 4) => {
                return Some(u32::from_le_bytes(data.try_into().ok()?) as usize);
            }
            (b"BC", 2) => return Some(u16::from_le_bytes([data[0], data[1]]) as usize + 1),
            _ => fields = &fields[4 + len..],
        }
    }
    None
}

/// Whether the first member of `file` records its size, without consuming it.
pub(crate) fn has_sized_members(file: &mut BufReader<File>) -> std::io::Result<bool> {
    let buffer = file.fill_buf()?;
    let Some(header) = buffer.get(..HEADER_WITH_XLEN) else {
        return Ok(false);
    };
    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
    let Some(extra) = buffer.get(HEADER_WITH_XLEN..HEADER_WITH_XLEN + xlen) else {
        return Ok(false);
    };
    Ok(member_size(header, extra).is_some_and(|size| size >= HEADER_WITH_XLEN + xlen + 8))
}

pub struct ParallelGzReader {
    outputs: Receiver<Output>,
    pending: BTreeMap<u64, Vec<u8>>,
    next: u64,
    current: Vec<u8>,
    position: usize,
    threads: Vec<JoinHandle<()>>,
}

impl ParallelGzReader {
    pub(crate) fn new(file: BufReader<File>) -> Self {
        let (jobs_tx, jobs_rx) = bounded::<Job>(WORKERS * QUEUE_PER_WORKER);
        let (outputs_tx, outputs_rx) = bounded::<Output>(WORKERS * QUEUE_PER_WORKER);
        let mut threads = Vec::with_capacity(WORKERS + 1);
        for _ in 0..WORKERS {
            let jobs = jobs_rx.clone();
            let outputs = outputs_tx.clone();
            threads.push(std::thread::spawn(move || inflate_members(&jobs, &outputs)));
        }
        threads.push(std::thread::spawn(move || {
            split_members(file, &jobs_tx, &outputs_tx)
        }));
        Self {
            outputs: outputs_rx,
            pending: BTreeMap::new(),
            next: 0,
            current: Vec::new(),
            position: 0,
            threads,
        }
    }

    /// Make the next member's output current. `false` at end of input.
    fn advance(&mut self) -> std::io::Result<bool> {
        loop {
            if let Some(output) = self.pending.remove(&self.next) {
                self.next += 1;
                self.current = output;
                self.position = 0;
                return Ok(true);
            }
            match self.outputs.recv() {
                Ok((sequence, Ok(output))) => {
                    self.pending.insert(sequence, output);
                }
                Ok((_, Err(error))) => return Err(error),
                // Every sender is gone and nothing is pending: end of input.
                Err(_) => return Ok(false),
            }
        }
    }
}

impl Drop for ParallelGzReader {
    fn drop(&mut self) {
        // Closing the output channel makes the splitter and workers stop.
        let (_, closed) = bounded(0);
        self.outputs = closed;
        for thread in self.threads.drain(..) {
            let _ = thread.join();
        }
    }
}

impl Read for ParallelGzReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let available = self.fill_buf()?;
        let count = available.len().min(buf.len());
        buf[..count].copy_from_slice(&available[..count]);
        self.consume(count);
        Ok(count)
    }
}

impl BufRead for ParallelGzReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        while self.position == self.current.len() {
            if !self.advance()? {
                return Ok(&[]);
            }
        }
        Ok(&self.current[self.position..])
    }

    fn consume(&mut self, amount: usize) {
        self.position = (self.position + amount).min(self.current.len());
    }
}

fn invalid(message: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, message.to_string())
}

/// Cut the file into members by their recorded sizes. A member without a
/// size switches the rest of the file to sequential decoding here.
fn split_members(mut file: BufReader<File>, jobs: &Sender<Job>, outputs: &Sender<Output>) {
    let mut sequence = 0u64;
    loop {
        let mut header = [0u8; HEADER_WITH_XLEN];
        match read_full(&mut file, &mut header) {
            Ok(0) => return,
            Ok(HEADER_WITH_XLEN) => {}
            Ok(_) => {
                let _ = outputs.send((sequence, Err(invalid("truncated gzip member header"))));
                return;
            }
            Err(error) => {
                let _ = outputs.send((sequence, Err(error)));
                return;
            }
        }
        let xlen = if header[3] & 4 != 0 {
            u16::from_le_bytes([header[10], header[11]]) as usize
        } else {
            0
        };
        let mut extra = vec![0u8; xlen];
        if header[3] & 4 != 0
            && let Err(error) = file.read_exact(&mut extra)
        {
            let _ = outputs.send((sequence, Err(error)));
            return;
        }
        let Some(size) =
            member_size(&header, &extra).filter(|&size| size >= HEADER_WITH_XLEN + xlen + 8)
        else {
            // No usable size: decode this member and the rest in order here.
            let mut prefix = header.to_vec();
            if header[3] & 4 != 0 {
                prefix.extend_from_slice(&extra);
            } else {
                // Without FEXTRA the last two header bytes start the next field.
                prefix.truncate(10);
                prefix.extend_from_slice(&header[10..]);
            }
            let rest = std::io::Cursor::new(prefix).chain(file);
            decode_sequentially(rest, sequence, outputs);
            return;
        };
        let mut member = Vec::with_capacity(size);
        member.extend_from_slice(&header);
        member.extend_from_slice(&extra);
        member.resize(size, 0);
        if let Err(error) = file.read_exact(&mut member[HEADER_WITH_XLEN + xlen..]) {
            let _ = outputs.send((sequence, Err(error)));
            return;
        }
        if jobs.send(Job::Member(sequence, member)).is_err() {
            return;
        }
        sequence += 1;
    }
}

fn decode_sequentially(input: impl Read, mut sequence: u64, outputs: &Sender<Output>) {
    let mut decoder = flate2::read::MultiGzDecoder::new(input);
    loop {
        let mut chunk = vec![0u8; FALLBACK_CHUNK];
        let mut filled = 0;
        while filled < chunk.len() {
            match decoder.read(&mut chunk[filled..]) {
                Ok(0) => break,
                Ok(count) => filled += count,
                Err(error) => {
                    let _ = outputs.send((sequence, Err(error)));
                    return;
                }
            }
        }
        if filled == 0 {
            return;
        }
        chunk.truncate(filled);
        if outputs.send((sequence, Ok(chunk))).is_err() {
            return;
        }
        sequence += 1;
    }
}

fn inflate_members(jobs: &Receiver<Job>, outputs: &Sender<Output>) {
    while let Ok(Job::Member(sequence, member)) = jobs.recv() {
        if outputs.send((sequence, inflate_member(&member))).is_err() {
            return;
        }
    }
}

/// Inflate exactly one member; it must span the whole slice.
fn inflate_member(member: &[u8]) -> std::io::Result<Vec<u8>> {
    let size_hint = member
        .len()
        .checked_sub(4)
        .and_then(|start| member.get(start..))
        .map_or(0, |isize| {
            u32::from_le_bytes(isize.try_into().expect("four bytes")) as usize
        });
    let mut input = member;
    let mut decoder = flate2::bufread::GzDecoder::new(&mut input);
    let mut output = Vec::with_capacity(size_hint);
    decoder.read_to_end(&mut output)?;
    drop(decoder);
    if !input.is_empty() {
        return Err(invalid("gzip member is shorter than its recorded size"));
    }
    Ok(output)
}

fn read_full(file: &mut impl Read, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match file.read(&mut buf[filled..])? {
            0 => break,
            count => filled += count,
        }
    }
    Ok(filled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn member(payload: &[u8], tag: &[u8; 2]) -> Vec<u8> {
        let mut body = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
        body.write_all(payload).unwrap();
        let body = body.finish().unwrap();
        let data_len = if tag == b"FP" { 4 } else { 2 };
        let total = HEADER_WITH_XLEN + 4 + data_len + body.len() + 8;
        let mut out = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff];
        out.extend_from_slice(&((4 + data_len) as u16).to_le_bytes());
        out.extend_from_slice(tag);
        out.extend_from_slice(&(data_len as u16).to_le_bytes());
        if tag == b"FP" {
            out.extend_from_slice(&(total as u32).to_le_bytes());
        } else {
            out.extend_from_slice(&((total - 1) as u16).to_le_bytes());
        }
        out.extend_from_slice(&body);
        let mut crc = flate2::Crc::new();
        crc.update(payload);
        out.extend_from_slice(&crc.sum().to_le_bytes());
        out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        assert_eq!(out.len(), total);
        out
    }

    fn decode(bytes: &[u8]) -> std::io::Result<Vec<u8>> {
        let path = std::env::temp_dir().join(format!(
            "kallistors-pgz-{}-{}.gz",
            std::process::id(),
            bytes.len()
        ));
        std::fs::write(&path, bytes).unwrap();
        let mut file = BufReader::new(File::open(&path).unwrap());
        assert!(has_sized_members(&mut file).unwrap());
        let mut out = Vec::new();
        let result = ParallelGzReader::new(file)
            .read_to_end(&mut out)
            .map(|_| out);
        std::fs::remove_file(&path).unwrap();
        result
    }

    #[test]
    fn sized_members_decode_in_order_with_sequential_tail() {
        let payloads: Vec<Vec<u8>> = (0..40u8)
            .map(|i| {
                (0..(i as usize * 997 + 1))
                    .map(|j| b"ACGT\n"[(j + i as usize) % 5])
                    .collect()
            })
            .collect();
        let mut file = Vec::new();
        let mut expected = Vec::new();
        for (i, payload) in payloads.iter().enumerate() {
            file.extend(member(payload, if i % 2 == 0 { b"FP" } else { b"BC" }));
            expected.extend_from_slice(payload);
        }
        // An ordinary member with no size field: the rest decodes in order.
        let mut plain = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        plain.write_all(b"@tail\nACGT\n+\nIIII\n").unwrap();
        file.extend(plain.finish().unwrap());
        expected.extend_from_slice(b"@tail\nACGT\n+\nIIII\n");
        file.extend(member(b"after", b"FP"));
        expected.extend_from_slice(b"after");
        assert_eq!(decode(&file).unwrap(), expected);
    }

    #[test]
    fn wrong_recorded_size_is_an_error() {
        let mut file = member(b"hello world", b"FP");
        let mut second = member(b"next", b"FP");
        // Claim four extra bytes, which belong to the following member.
        let claimed = file.len() as u32 + 4;
        file[16..20].copy_from_slice(&claimed.to_le_bytes());
        file.append(&mut second);
        assert!(decode(&file).is_err());
    }
}
