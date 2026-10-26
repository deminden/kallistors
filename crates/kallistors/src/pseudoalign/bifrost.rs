mod buckets;
mod candidate_cache;
pub(crate) use candidate_cache::CachedBucket;
mod decoding;
mod ec_storage;
mod positional_arena;
mod positional_parallel;
use positional_arena::{EcBlocksBuilder, OwnedEcBlocks};
mod encoded;
mod large_one_pass;
mod loader;
mod loader_probe;

use buckets::{
    UnitigBlockState, advance_unitig_state, build_kmer_pos_index, populate_kmer_minimizer_counts,
    populate_kmer_minimizer_storage, populate_unitig_minimizer_counts,
    populate_unitig_minimizer_storage, precompute_encoded_unitigs,
};
#[cfg(test)]
use decoding::append_compressed_words;
use decoding::{append_buffered_compressed, decode_compressed_sequence, decode_kmer};
#[cfg(test)]
use ec_storage::FlatEcBlocks;
use ec_storage::{FlatEcIndex, MinzOffsets, MinzPositionIndex, SelectiveFallbackKmerIndex};
#[cfg(test)]
use encoded::EncodedSequence;
use encoded::{EncodedBaseStore, EncodedUnitigs, SequenceMetadata};
use loader::build_bifrost_index_with_positions_threaded_impl;

#[path = "bifrost/positions.rs"]
mod load_positions;
#[path = "bifrost/one_pass.rs"]
mod one_pass;
use load_positions::BitmapPositions;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::ops::Index;
use std::path::Path;
use std::sync::atomic::{AtomicU16, AtomicU64, AtomicUsize, Ordering};

use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use rdst::RadixSort;

use crate::index::bifrost::{BooPhf, read_bitcontainer_values};
use crate::index::read_block_array_blocks;
use crate::timing::{self, Stage};
use crate::{Error, Result};

use super::ec::encode_kmer_pair;
use super::fingerprints::{BucketCountKey, TaggedBucket};
use super::io::{read_i32_le, read_u32_le, read_u64_le};
use super::{KmerEcIndex, read_onlist_with_lengths};

pub type KmerPosIndex = HashMap<u64, Vec<(usize, usize, bool)>>;

const LOAD_BATCH_POSITIONS: usize = 1 << 20;
const LOAD_MIN_PAR_CHUNK: usize = 1 << 14;
const INVALID_BASE_CODE: u8 = 4;
const BASE_CODES: [u8; 256] = {
    let mut codes = [INVALID_BASE_CODE; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'g' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes[b't' as usize] = 3;
    codes
};

const DECODED_BYTE_BASES: [[u8; 4]; 256] = {
    let mut table = [[0; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut pos = 0;
        while pos < 4 {
            table[byte][pos] = b"ACGT"[(byte >> (pos * 2)) & 3];
            pos += 1;
        }
        byte += 1;
    }
    table
};

static NEXT_CACHE_ID: AtomicU64 = AtomicU64::new(1);

pub struct BifrostIndex {
    pub(crate) cache_id: u64,
    pub index_version: u64,
    pub k: usize,
    pub g: usize,
    pub unitigs: Vec<Vec<u8>>,
    pub km_unitigs: Vec<Vec<u8>>,
    pub h_kmers: Vec<Vec<u8>>,
    pub h_kmer_map: Option<HashMap<u64, usize>>,
    pub dlist: Option<HashSet<u64>>,
    pub(crate) encoded_unitigs: EncodedUnitigs,
    pub kmer_pos_index: Option<KmerPosIndex>,
    pub kmer_index: Option<KmerEcIndex>,
    pub(crate) flat_ec: FlatEcIndex,
    pub minz_positions: MinzPositionIndex,
    large_candidate_cache: Vec<candidate_cache::CachedBucket>,
    candidate_cache_slots: Vec<u32>,
    pub(crate) minz_fingerprints: Vec<u16>,
    pub mphf: BooPhf,
    pub(crate) selective_fallback_kmers: SelectiveFallbackKmerIndex,
    pub onlist: Option<Vec<bool>>,
    pub transcript_lengths: Vec<u32>,
    pub transcript_names: Vec<String>,
    deferred_names: Option<super::onlist::TranscriptNames>,
    pub shade_to_color: Vec<Option<u32>>,
    pub shade_sequences: Vec<bool>,
    pub use_shade: bool,
}

impl BifrostIndex {
    #[doc(hidden)]
    pub fn take_transcript_names(&mut self) -> super::onlist::TranscriptNames {
        self.deferred_names.take().unwrap_or_else(|| {
            super::onlist::TranscriptNames::eager(std::mem::take(&mut self.transcript_names))
        })
    }
    /// Release graph storage after all mapping workers have joined. All bytes
    /// remain owned through normal deallocation; only unused pages are reset.
    #[doc(hidden)]
    pub fn retire(self) {
        self.encoded_unitigs.unitig_bases.retire();
        self.encoded_unitigs.km_bases.retire();
        self.flat_ec.unitig_block_offsets.retire();
        self.flat_ec.blocks.retire();
        crate::index::scratch_pages::retire(self.flat_ec.ec_values);
        self.minz_positions.retire();
        crate::index::scratch_pages::retire(self.minz_fingerprints);
        self.mphf.retire_storage();
        crate::index::scratch_pages::retire(self.unitigs);
        crate::index::scratch_pages::retire(self.km_unitigs);
        crate::index::scratch_pages::retire(self.h_kmers);
        self.flat_ec.owned_ec_blocks.retire();
    }

    /// Number of ordinary graph nodes, independent of decoded sequence storage.
    #[inline]
    pub fn unitig_count(&self) -> usize {
        self.encoded_unitigs.unitig_bases.metadata.len()
    }

    /// Return ordinary graph bases, decoding packed storage only when the
    /// counting loader has released its duplicate ASCII sequence.
    pub fn unitig_sequence(&self, id: usize) -> Option<Cow<'_, [u8]>> {
        decoded_graph_sequence(&self.unitigs, &self.encoded_unitigs.unitig_bases, id)
    }
    /// Number of fixed-length graph nodes, independent of ASCII storage.
    #[inline]
    pub fn km_unitig_count(&self) -> usize {
        self.encoded_unitigs.km_bases.metadata.len()
    }
    /// Return fixed-length graph bases, decoding packed storage on demand.
    pub fn km_unitig_sequence(&self, id: usize) -> Option<Cow<'_, [u8]>> {
        decoded_graph_sequence(&self.km_unitigs, &self.encoded_unitigs.km_bases, id)
    }

    pub fn has_strand_annotations(&self) -> bool {
        self.flat_ec
            .owned_ec_blocks
            .iter()
            .flat_map(|blocks| blocks.iter())
            .any(|block| {
                block
                    .positions
                    .as_ref()
                    .is_some_and(|positions| positions.len() > 0)
            })
    }
}

fn decoded_graph_sequence<'a>(
    raw: &'a [Vec<u8>],
    packed: &'a EncodedBaseStore,
    id: usize,
) -> Option<Cow<'a, [u8]>> {
    if let Some(raw) = raw.get(id)
        && !raw.is_empty()
    {
        return Some(Cow::Borrowed(raw));
    }
    let seq = packed.sequence(id)?;
    Some(Cow::Owned(
        (0..seq.len)
            .map(|i| b"ACGT"[seq.kmer_at(i, 1).expect("validated packed base") as usize])
            .collect(),
    ))
}

pub fn build_bifrost_index(path: &Path) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, false, 1)
}

fn associate_abundant_node(
    kmers: &mut [Vec<u8>],
    map: &mut HashMap<u64, usize>,
    base: usize,
    node_row: usize,
    node_kmer: &[u8],
) -> Result<()> {
    let (forward, reverse) = encode_kmer_pair(node_kmer)
        .ok_or_else(|| Error::InvalidFormat("invalid abundant node k-mer".into()))?;
    let id = map
        .get_mut(&forward.min(reverse))
        .ok_or_else(|| Error::InvalidFormat("abundant node k-mer absent from graph".into()))?;
    let slot = node_row
        .checked_sub(base)
        .and_then(|offset| kmers.get_mut(offset))
        .ok_or_else(|| Error::InvalidFormat("abundant node out of range".into()))?;
    slot.copy_from_slice(node_kmer);
    *id = node_row;
    Ok(())
}

fn associate_abundant_word(
    words: &mut [u64],
    map: &mut HashMap<u64, usize>,
    base: usize,
    node_row: usize,
    node_kmer: &[u8],
) -> Result<()> {
    let (forward, reverse) = encode_kmer_pair(node_kmer)
        .ok_or_else(|| Error::InvalidFormat("invalid abundant node k-mer".into()))?;
    let id = map
        .get_mut(&forward.min(reverse))
        .ok_or_else(|| Error::InvalidFormat("abundant node k-mer absent from graph".into()))?;
    let slot = node_row
        .checked_sub(base)
        .and_then(|offset| words.get_mut(offset))
        .ok_or_else(|| Error::InvalidFormat("abundant node out of range".into()))?;
    *slot = forward
        .checked_shl((64 - 2 * node_kmer.len()) as u32)
        .unwrap_or(0);
    *id = node_row;
    Ok(())
}

/// Load a flat counting index without duplicate EC blocks. Positional loading retains
/// the nested metadata needed for fragment/strand filters. Alternate k-mer scanners
/// and detailed hit tracing should use the full index loader instead.
pub fn build_bifrost_index_quant_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(
        path,
        load_positional_info,
        threads,
        load_positional_info,
        false,
        false,
    )
}

/// Load a counting index using packed graph sequences without ASCII copies.
/// Use `unitig_count` / `km_unitig_count` for identifiers and the matching
/// sequence helpers for cold decoding. Full loaders retain ASCII graph fields.
#[doc(hidden)]
pub fn build_bifrost_index_quant_packed_threaded(
    path: &Path,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(path, false, threads, false, true, false)
}

#[doc(hidden)]
pub fn build_bifrost_index_quant_deferred_threaded(
    path: &Path,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(path, false, threads, false, true, true)
}

/// Load positional and strand annotations while keeping graph bases packed.
/// Diagnostic and alternate scanners retain their full graph loader.
#[doc(hidden)]
pub fn build_bifrost_index_positional_packed_threaded(
    path: &Path,
    threads: usize,
    defer_names: bool,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(path, true, threads, true, true, defer_names)
}

pub fn build_bifrost_index_with_kmer_pos(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    let mut index = build_bifrost_index_with_positions_threaded(path, load_positional_info, 1)?;
    index.kmer_pos_index = Some(build_kmer_pos_index(&index));
    Ok(index)
}

pub fn build_bifrost_index_with_positions(
    path: &Path,
    load_positional_info: bool,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, load_positional_info, 1)
}

pub fn build_bifrost_index_threaded(path: &Path, threads: usize) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded(path, false, threads)
}

pub fn build_bifrost_index_with_kmer_pos_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    let mut index =
        build_bifrost_index_with_positions_threaded(path, load_positional_info, threads)?;
    index.kmer_pos_index = Some(build_kmer_pos_index(&index));
    Ok(index)
}

pub fn build_bifrost_index_with_positions_threaded(
    path: &Path,
    load_positional_info: bool,
    threads: usize,
) -> Result<BifrostIndex> {
    build_bifrost_index_with_positions_threaded_impl(
        path,
        load_positional_info,
        threads,
        true,
        false,
        false,
    )
}

struct KmerLoadContext<'a> {
    km_bases: &'a EncodedBaseStore,
    g: usize,
    km_glen: usize,
    mphf: &'a BooPhf,
}

#[cfg(test)]
fn base_code(base: u8) -> Option<u64> {
    let code = BASE_CODES[base as usize];
    if code == INVALID_BASE_CODE {
        None
    } else {
        Some(u64::from(code))
    }
}

#[cfg(test)]
fn encode_word_pair(seq: &[u8]) -> Option<(u64, u64)> {
    if seq.is_empty() {
        return Some((0, 0));
    }
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let last = seq.len() - 1;
    for (i, &base) in seq.iter().enumerate() {
        let code = base_code(base)?;
        fwd |= code << (62 - ((i & 0x1f) << 1));
        rev |= (3 - code) << (62 - (((last - i) & 0x1f) << 1));
    }
    Some((fwd, rev))
}

#[inline]
fn minimizer_count(seq_len: usize, g: usize) -> usize {
    if g == 0 || seq_len < g {
        0
    } else {
        seq_len - g + 1
    }
}

#[cfg(test)]
fn compute_minimizer_reps_into(seq: &[u8], g: usize, out: &mut [[u8; 8]]) {
    if out.is_empty() {
        return;
    }

    debug_assert_eq!(out.len(), minimizer_count(seq.len(), g));

    let Some((mut fwd, mut rev)) = encode_word_pair(&seq[..g]) else {
        return;
    };
    let tail_shift = 64usize.saturating_sub(g * 2);
    let high_mask = if g == 32 {
        u64::MAX
    } else {
        (!0u64) << tail_shift
    };

    out[0] = fwd.min(rev).to_le_bytes();

    for (dst, &base) in out.iter_mut().skip(1).zip(&seq[g..]) {
        let code = BASE_CODES[base as usize];
        if code == INVALID_BASE_CODE {
            return;
        }
        let code = u64::from(code);
        let comp = 3 - code;
        fwd = ((fwd << 2) & high_mask) | (code << tail_shift);
        rev = ((rev >> 2) & high_mask) | (comp << 62);
        *dst = fwd.min(rev).to_le_bytes();
    }
}

fn pack_bases_into(seq: &[u8], out: &mut [u64]) {
    debug_assert_eq!(seq.len().div_ceil(32), out.len());
    for (chunk, word) in seq.chunks(32).zip(out) {
        let mut encoded = 0u64;
        for &base in chunk {
            encoded = (encoded << 2) | u64::from(BASE_CODES[base as usize] & 3);
        }
        *word = encoded << ((32 - chunk.len()) * 2);
    }
}

#[cfg(test)]
#[path = "bifrost/packed_base_tests.rs"]
mod packed_base_tests;

pub(super) fn read_unitig_sequences<R: Read + Seek>(
    reader: &mut R,
    dbg_start: u64,
    dbg_size: u64,
    k: u32,
    kmer_bytes: usize,
) -> Result<(Vec<Vec<u8>>, usize)> {
    reader.seek(SeekFrom::Start(dbg_start))?;
    let _file_format_version = read_u64_le(reader)?;
    let _k = read_i32_le(reader)?;
    let _g = read_i32_le(reader)?;

    let v_unitigs_sz = read_u64_le(reader)?;
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(v_unitigs_sz as usize);
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(reader)? as usize;
        let data_sz = len.div_ceil(4);
        let mut data = vec![0u8; data_sz];
        reader.read_exact(&mut data)?;
        seqs.push(decode_compressed_sequence(&data, len));
    }

    let km_unitigs_sz = read_u64_le(reader)?;
    for _ in 0..km_unitigs_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    let abundant_base = seqs.len();
    let h_kmers_ccov_sz = read_u64_le(reader)?;
    for _ in 0..h_kmers_ccov_sz {
        let mut buf = vec![0u8; kmer_bytes];
        reader.read_exact(&mut buf)?;
        seqs.push(decode_kmer(&buf, k as usize));
    }

    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;
    Ok((seqs, abundant_base))
}

#[cfg(test)]
#[path = "bifrost/loader_order_tests.rs"]
mod loader_order_tests;
