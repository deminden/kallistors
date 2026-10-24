//! Read-layout entry points and fast/general mapper dispatch.

use super::fast::ec_for_read_bifrost_fast;
use super::runtime::{self, fast_path_enabled};
use crate::pseudoalign::debug::ReadDebugState;
use crate::pseudoalign::{BifrostIndex, BorrowedReadEc, PseudoalignOptions, ReadEc, Strand};

mod finish;
mod general;
mod probe;

// Counting callers have finished placement sampling and do not consume hit metadata.
#[inline(always)]
pub(crate) fn ec_for_read_bifrost_counting_borrowed<'index>(
    index: &'index BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<BorrowedReadEc<'index>> {
    if index.k > 0 && seq.len() >= index.k && fast_path_enabled(index, false, options) {
        runtime::ensure_cache_context(index);
        if options.strand_specific.is_some() {
            super::fast::ec_for_read_bifrost_fast_counting_borrowed::<true>(
                index, seq, strand, options,
            )
        } else {
            super::fast::ec_for_read_bifrost_fast_counting_borrowed::<false>(
                index, seq, strand, options,
            )
        }
    } else {
        ec_for_read_bifrost(index, seq, strand, None, options).map(BorrowedReadEc::owned)
    }
}

pub(crate) fn ec_for_read_bifrost_counting(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    ec_for_read_bifrost_counting_borrowed(index, seq, strand, options)
        .map(BorrowedReadEc::into_owned)
}

// Single-end fragment filtering needs the best placement, but ordinary
// unstranded counting does not consume the complete hit vector.
pub(crate) fn ec_for_read_bifrost_single(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    ec_for_read_bifrost_single_borrowed(index, seq, strand, options).map(BorrowedReadEc::into_owned)
}

pub(crate) fn ec_for_read_bifrost_single_borrowed<'index>(
    index: &'index BifrostIndex,
    seq: &[u8],
    strand: Strand,
    options: PseudoalignOptions,
) -> Option<BorrowedReadEc<'index>> {
    if index.k > 0
        && seq.len() >= index.k
        && !options.investigation.any_enabled()
        && fast_path_enabled(index, false, options)
    {
        runtime::ensure_cache_context(index);
        super::fast::ec_for_read_bifrost_fast_counting_borrowed::<true>(index, seq, strand, options)
    } else {
        ec_for_read_bifrost(index, seq, strand, None, options).map(BorrowedReadEc::owned)
    }
}

// Main parity scanner: stream hits, keep an online EC intersection, and jump when safe.
pub(crate) fn ec_for_read_bifrost(
    index: &BifrostIndex,
    seq: &[u8],
    strand: Strand,
    dbg: Option<&mut ReadDebugState>,
    options: PseudoalignOptions,
) -> Option<ReadEc> {
    if seq.len() < index.k {
        if let Some(state) = dbg {
            state.saw_valid_kmer = false;
        }
        return None;
    }
    runtime::ensure_cache_context(index);
    if fast_path_enabled(index, dbg.is_some(), options) {
        return ec_for_read_bifrost_fast(index, seq, strand, dbg, options);
    }
    general::scan_read(index, seq, strand, dbg, options)
}
