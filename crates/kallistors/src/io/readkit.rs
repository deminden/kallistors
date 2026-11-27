//! Sequence batches from readkit's columnar read batches.
//!
//! readkit keeps each mate's bases in one contiguous buffer with an offset
//! per read, which is the layout of [`PackedSeqBatch`]. [`take_sequences`]
//! moves those buffers out of a batch instead of copying them; the batch's
//! pool gets empty vectors back. [`copy_sequences`] leaves the batch intact
//! for callers that still need it, for example to forward it to a QC stage.

use std::ops::Range;

use readkit_batch::{Batch, Layout};

use super::PackedSeqBatch;
use crate::{Error, Result};

/// One sequence batch per segment (mate) of `batch`, moving the bases out.
pub fn take_sequences<L: Layout>(batch: Batch<L>) -> Result<Vec<PackedSeqBatch>> {
    batch
        .into_parts()
        .segments
        .into_iter()
        .map(|mut segment| {
            let spans = spans(&segment.offsets, 0)?;
            PackedSeqBatch::from_parts(std::mem::take(segment.seqs.vec_mut()), spans)
        })
        .collect()
}

/// One sequence batch per segment (mate) of `batch`, copying the bases.
pub fn copy_sequences<L: Layout>(batch: &Batch<L>) -> Result<Vec<PackedSeqBatch>> {
    batch
        .columns()
        .map(|column| {
            let seqs = column.seqs();
            let offsets = seqs.offsets();
            let base = offsets.first().copied().unwrap_or(0);
            let spans = spans(offsets, base)?;
            let end = offsets.last().copied().unwrap_or(base);
            let bytes = usize::try_from(base)
                .ok()
                .zip(usize::try_from(end).ok())
                .and_then(|(start, end)| seqs.values().get(start..end))
                .ok_or_else(|| Error::InvalidFormat("readkit offsets outside the batch".into()))?;
            PackedSeqBatch::from_parts(bytes.to_vec(), spans)
        })
        .collect()
}

fn spans(offsets: &[i32], base: i32) -> Result<Vec<Range<u32>>> {
    offsets
        .windows(2)
        .map(|w| {
            match (
                u32::try_from(w[0] - base).ok(),
                u32::try_from(w[1] - base).ok(),
            ) {
                (Some(start), Some(end)) if start <= end => Ok(start..end),
                _ => Err(Error::InvalidFormat("readkit offsets out of order".into())),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use readkit_batch::{BatchBuilder, Paired, Pool, Read};

    fn paired(pool: &Pool) -> Batch<Paired> {
        let mut builder = BatchBuilder::<Paired>::new(pool);
        builder
            .push_pair(Read::new(b"ACGT", b"IIII"), Read::new(b"GG", b"II"))
            .unwrap();
        builder
            .push_pair(Read::new(b"", b""), Read::new(b"TTTAA", b"IIIII"))
            .unwrap();
        builder.finish()
    }

    fn seqs(batch: &PackedSeqBatch) -> Vec<Vec<u8>> {
        batch.records().map(|r| r.seq.to_vec()).collect()
    }

    #[test]
    fn taken_and_copied_sequences_match_the_batch() {
        let pool = Pool::new();
        let copied = copy_sequences(&paired(&pool)).unwrap();
        let taken = take_sequences(paired(&pool)).unwrap();
        for mates in [&copied, &taken] {
            assert_eq!(mates.len(), 2);
            assert_eq!(seqs(&mates[0]), [b"ACGT".to_vec(), Vec::new()]);
            assert_eq!(seqs(&mates[1]), [b"GG".to_vec(), b"TTTAA".to_vec()]);
        }
    }
}
