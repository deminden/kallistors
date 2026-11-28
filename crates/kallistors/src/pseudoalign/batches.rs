//! Mapping reads from readkit batch sources (`readkit` feature).
//!
//! These take the place of the FASTQ readers when reads arrive already
//! decoded, for example from a trimming stage in the same process. Each
//! batch's sequence buffers move into the mapping pipeline without a copy.
//! The fragment-length sample depends only on the order of the pairs, not on
//! how the source cuts them into batches.

use readkit_batch::{Paired, Single};
use readkit_io::BatchSource;

use crate::io::PackedSeqBatch;
use crate::io::readkit::take_sequences;
use crate::{Error, Result};

use super::threaded::{paired_pipeline, single_end_pipeline, with_producer};
use super::{BifrostIndex, EcCounts, FragmentFilter, PseudoalignOptions, Strand};

fn source_error<E: std::error::Error + Send + Sync + 'static>(err: E) -> Error {
    Error::Source(Box::new(err))
}

fn next_pair<S: BatchSource<Paired>>(
    source: &mut S,
) -> Result<Option<(PackedSeqBatch, PackedSeqBatch)>> {
    loop {
        let Some(batch) = source.next_batch().map_err(source_error)? else {
            return Ok(None);
        };
        if batch.is_empty() {
            continue;
        }
        let mut mates = take_sequences(batch)?.into_iter();
        let (Some(left), Some(right)) = (mates.next(), mates.next()) else {
            return Err(Error::InvalidFormat(
                "paired batch without two mates".into(),
            ));
        };
        return Ok(Some((left, right)));
    }
}

fn next_single<S: BatchSource<Single>>(source: &mut S) -> Result<Option<PackedSeqBatch>> {
    loop {
        let Some(batch) = source.next_batch().map_err(source_error)? else {
            return Ok(None);
        };
        if batch.is_empty() {
            continue;
        }
        return take_sequences(batch)?
            .pop()
            .map(Some)
            .ok_or_else(|| Error::InvalidFormat("single-end batch without reads".into()));
    }
}

/// Maps paired reads from a readkit source.
///
/// A source error ends the run and comes back as [`Error::Source`], with the
/// source's own error inside.
pub fn pseudoalign_paired_batches_threaded<S: BatchSource<Paired> + Send>(
    index: &BifrostIndex,
    source: &mut S,
    strand: Strand,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    if threads <= 1 {
        return paired_pipeline(index, &mut || next_pair(source), strand, options, threads);
    }
    with_producer(
        |tx| {
            while let Some(pair) = next_pair(source)? {
                if tx.send(pair).is_err() {
                    break;
                }
            }
            Ok(())
        },
        |next| paired_pipeline(index, next, strand, options, threads),
    )
}

/// Maps single-end reads from a readkit source.
///
/// A source error ends the run and comes back as [`Error::Source`], with the
/// source's own error inside.
pub fn pseudoalign_single_end_batches_threaded<S: BatchSource<Single> + Send>(
    index: &BifrostIndex,
    source: &mut S,
    strand: Strand,
    filter: Option<FragmentFilter>,
    options: PseudoalignOptions,
    threads: usize,
) -> Result<EcCounts> {
    if threads <= 1 {
        return single_end_pipeline(
            index,
            &mut || next_single(source),
            strand,
            filter,
            options,
            threads,
        );
    }
    with_producer(
        |tx| {
            while let Some(batch) = next_single(source)? {
                if tx.send(batch).is_err() {
                    break;
                }
            }
            Ok(())
        },
        |next| single_end_pipeline(index, next, strand, filter, options, threads),
    )
}
