//! Bounded framing with independent positional decoders and ordered publication.
use super::positional_arena::EcBlocksBuilder;
use crate::Result;
use rayon::prelude::*;
use std::io::{BufRead, Seek};
use std::ops::Range;
// The pinned C library's x86 ISA cache is atomic on C11 platforms, but
// its MSVC/nonatomic fallback also needs serial initialization before these
// decoder tasks. The C entry point takes no pointers or ownership.
#[cfg(target_arch = "x86_64")]
fn prepare_decoder_isa() {
    unsafe extern "C" {
        fn croaring_hardware_support() -> std::ffi::c_int;
    }
    static INIT: std::sync::Once = std::sync::Once::new();
    INIT.call_once(|| {
        // SAFETY: linked croaring-sys5.2.2 declares this no-argument C function
        // on x86-64; it returns CPU feature flags and has no caller preconditions.
        // Finish ISA-cache initialization before any new parallel decode task.
        let _ = unsafe { croaring_hardware_support() };
    });
}
#[cfg(not(target_arch = "x86_64"))]
fn prepare_decoder_isa() {}
const RAW_BYTES: usize = 64 * 1024 * 1024;
const FRAMES: usize = 1 << 16;
const TASK: usize = 1 << 12;
struct Frames {
    bytes: Vec<u8>,
    ranges: Vec<Range<usize>>,
    raw: usize,
    frames: usize,
    task: usize,
}
impl Frames {
    fn new(raw: usize, frames: usize, task: usize) -> Self {
        assert!(raw > 0 && frames > 0 && task > 0);
        Self {
            bytes: Vec::with_capacity(raw),
            ranges: Vec::with_capacity(frames),
            raw,
            frames,
            task,
        }
    }
    fn flush(&mut self, pool: &rayon::ThreadPool, target: &mut EcBlocksBuilder) -> Result<()> {
        if self.ranges.is_empty() {
            self.bytes.clear();
            return Ok(());
        }
        // Indexed collection preserves the first serialized error, regardless
        // of which independent decoder finishes first. Never short-circuit
        // a parallel Result reduction whose error order would be unspecified.
        let results: Vec<_> = pool.install(|| {
            self.ranges
                .par_chunks(self.task)
                .map(|rows| {
                    let mut builder = EcBlocksBuilder::new(rows.len());
                    for range in rows {
                        builder.append_payload(&self.bytes[range.clone()], true)?;
                    }
                    Ok(builder.finish_for_transfer())
                })
                .collect::<Vec<Result<_>>>()
        });
        self.ranges.clear();
        self.bytes.clear();
        for result in results {
            target.append_finished(result?)?;
        }
        Ok(())
    }
}
pub(super) fn read_prefix<R: BufRead + Seek>(
    reader: &mut R,
    nodes: usize,
    k: usize,
    pool: &rayon::ThreadPool,
    target: &mut EcBlocksBuilder,
) -> Result<()> {
    prepare_decoder_isa();
    read_prefix_limited(
        reader,
        nodes,
        k,
        pool,
        target,
        Frames::new(RAW_BYTES, FRAMES, TASK),
    )
}
fn read_prefix_limited<R: BufRead + Seek>(
    reader: &mut R,
    nodes: usize,
    k: usize,
    pool: &rayon::ThreadPool,
    target: &mut EcBlocksBuilder,
    mut batch: Frames,
) -> Result<()> {
    let mut oversized = Vec::new();
    let framing = (|| -> Result<()> {
        for _ in 0..nodes {
            reader.seek_relative(k as i64)?;
            let size = super::read_u32_le(reader)? as usize;
            if size > batch.raw {
                batch.flush(pool, target)?;
                super::loader::with_buffered_node(reader, size, &mut oversized, |bytes| {
                    target.append_payload(bytes, true)
                })?;
                continue;
            }
            if batch.ranges.len() == batch.frames || size > batch.raw - batch.bytes.len() {
                batch.flush(pool, target)?;
            }
            let begin = batch.bytes.len();
            let end = begin + size;
            batch.bytes.resize(end, 0);
            reader.read_exact(&mut batch.bytes[begin..end])?;
            batch.ranges.push(begin..end);
        }
        Ok(())
    })();
    // A later framing/I/O error must not hide an earlier invalid body that
    // was already fully buffered. Private partial index state is dropped on
    // either error; abundant-node association stays on the original path.
    batch.flush(pool, target)?;
    framing
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Error;
    use std::io::{BufReader, Cursor, Read};
    fn bitmap(word: u32) -> Vec<u8> {
        let mut out = vec![1];
        out.extend_from_slice(&1u32.to_le_bytes());
        out.extend_from_slice(&word.to_le_bytes());
        out
    }
    fn payload(id: u32, blocks: usize) -> Vec<u8> {
        let mut out = id.to_le_bytes().to_vec();
        if blocks < 2 {
            out.push(blocks as u8);
        } else {
            out.push(2);
            out.extend_from_slice(&(blocks as u64).to_le_bytes());
        }
        for b in 0..blocks {
            out.extend_from_slice(&(b as u32 * 31).to_le_bytes());
            out.extend_from_slice(&(b as u32 * 31 + 29).to_le_bytes());
            let ids = bitmap(id + b as u32);
            out.extend_from_slice(&(ids.len() as u64).to_le_bytes());
            out.extend_from_slice(&ids);
            out.extend_from_slice(&1u64.to_le_bytes());
            let pos = bitmap(0x80000000 | b as u32);
            out.extend_from_slice(&(pos.len() as u64).to_le_bytes());
            out.extend_from_slice(&pos);
        }
        out
    }
    fn frame(out: &mut Vec<u8>, payload: &[u8]) {
        out.extend_from_slice(&[b'A'; 31]);
        out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        out.extend_from_slice(payload);
    }
    #[test]
    fn parallel_positional_frames_preserve_complete_fields_order_boundaries_and_errors() {
        let nodes = TASK * 2 + 17;
        let mut bytes = Vec::new();
        for id in 0..nodes {
            frame(&mut bytes, &payload(id as u32, id % 3));
        }
        bytes.push(0x5a);
        for workers in [1, 2, 3, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            for (capacity, limits) in [(7, None), (1 << 20, None), (7, Some((64, 3, 2)))] {
                let mut reader = BufReader::with_capacity(capacity, Cursor::new(bytes.as_slice()));
                let mut builder = EcBlocksBuilder::new(nodes);
                if let Some((raw, frames, task)) = limits {
                    read_prefix_limited(
                        &mut reader,
                        nodes,
                        31,
                        &pool,
                        &mut builder,
                        Frames::new(raw, frames, task),
                    )
                    .unwrap();
                } else {
                    read_prefix(&mut reader, nodes, 31, &pool, &mut builder).unwrap();
                }
                let owned = builder.finish();
                assert_eq!(owned.len(), nodes);
                for id in 0..nodes {
                    let blocks = &owned[id];
                    assert_eq!(blocks.len(), id % 3);
                    for (b, block) in blocks.iter().enumerate() {
                        assert_eq!((block.lb, block.ub), (b as u32 * 31, b as u32 * 31 + 29));
                        assert_eq!(&*block.ec, &[id as u32 + b as u32]);
                        let positions = block.positions.as_ref().unwrap();
                        assert_eq!(positions.len(), 1);
                        assert!(
                            matches!(positions.iter().next().unwrap(),crate::index::PositionSetRef::Single(word) if word==(0x80000000|b as u32))
                        );
                    }
                }
                let mut tail = [0];
                reader.read_exact(&mut tail).unwrap();
                assert_eq!(tail, [0x5a]);
            }
            // A fully read invalid body precedes a truncated next header.
            let mut bad = Vec::new();
            frame(&mut bad, &[0, 0, 0, 0, 3]);
            let mut reader = BufReader::new(Cursor::new(bad));
            let mut builder = EcBlocksBuilder::new(2);
            assert!(
                matches!(read_prefix(&mut reader,2,31,&pool,&mut builder),Err(Error::InvalidFormat(message)) if message=="invalid block array flag")
            );
            // Frame-count bounds apply even when payload size is zero.
            let mut bad = Vec::new();
            for _ in 0..=FRAMES {
                frame(&mut bad, &[]);
            }
            let mut reader = BufReader::new(Cursor::new(bad));
            let mut builder = EcBlocksBuilder::new(FRAMES + 1);
            assert!(matches!(
                read_prefix(&mut reader, FRAMES + 1, 31, &pool, &mut builder),
                Err(Error::Io(_))
            ));
        }
    }
}
