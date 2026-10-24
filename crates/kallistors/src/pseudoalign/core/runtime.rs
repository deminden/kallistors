//! Per-thread lookup context and mapper policy dispatch.

mod block_cache;
mod budget;
mod match_cache;
mod minimizer_cache;

use block_cache::{BlockLookupCache, BlockLookupCacheEntry};
use budget::MATCH_CACHE_DEFAULT_ENTRIES;
pub(crate) use budget::MatchCacheBudgetScope;
pub(crate) use match_cache::FastKmerMatch;
use match_cache::{FastKmerMatchCache, FastKmerMatchCacheEntry};
use minimizer_cache::{MinimizerCandidateCache, MinimizerCandidateCacheEntry};

#[cfg(test)]
mod cache_context_tests;
#[cfg(test)]
mod paired_match_prefetch_tests;

use std::cell::{Cell, RefCell};
use std::fmt::Write as _;
use std::sync::OnceLock;

use super::super::minimizers::minimizers_for_encoded_kmer_into;
use super::super::{BifrostIndex, PseudoalignOptions};

// Runtime switches and tiny per-thread caches keep hot read loops allocation-light.
pub(crate) static DISABLE_FAST_PATH: OnceLock<bool> = OnceLock::new();
pub(crate) static DISABLE_MINIMIZER_CACHE: OnceLock<bool> = OnceLock::new();
pub(crate) static DISABLE_MATCH_CACHE: OnceLock<bool> = OnceLock::new();
pub(crate) static ALLOW_TAIL_MINIMIZER_IN_BACKOFF: OnceLock<bool> = OnceLock::new();
pub(crate) static BACKOFF_DIRECT_RESCUE: OnceLock<bool> = OnceLock::new();
pub(crate) static ACCEPT_RELAXED_STREAM: OnceLock<bool> = OnceLock::new();

#[inline]
pub(crate) fn env_flag(cache: &'static OnceLock<bool>, name: &'static str) -> bool {
    *cache.get_or_init(|| std::env::var_os(name).is_some())
}

#[derive(Clone, Copy, Default)]
struct MphfCacheEntry {
    key: [u8; 8],
    // Original u32 index, full u16 fingerprint and two-bit presence state.
    value: u64,
}

pub(crate) struct MphfLookupCache {
    entries: Vec<MphfCacheEntry>,
    matches: FastKmerMatchCache,
    minimizers: MinimizerCandidateCache,
    primary_cursor: super::super::minimizers::PrimaryMinimizerCursor,
    query_fingerprints: super::super::fingerprints::QueryFingerprintCache,
    pub(super) scan_scratch: Option<Box<super::fast::FastScannerScratch>>,
    pub(super) slow_rev_buf: Vec<u8>,
}

impl Default for MphfLookupCache {
    fn default() -> Self {
        Self {
            entries: vec![MphfCacheEntry::default(); 2048],
            matches: FastKmerMatchCache::new(),
            minimizers: MinimizerCandidateCache::new(),
            primary_cursor: Default::default(),
            query_fingerprints: super::super::fingerprints::QueryFingerprintCache::default(),
            scan_scratch: Some(Box::default()),
            slow_rev_buf: Vec::new(),
        }
    }
}

impl MphfLookupCache {
    #[inline]
    pub(crate) fn query_fingerprint_matches(
        &mut self,
        forward: u64,
        reverse: u64,
        fingerprint: u16,
    ) -> bool {
        self.query_fingerprints
            .contains(forward, reverse, fingerprint)
    }
}

thread_local! {
    static MATCH_CACHE_ENTRIES: Cell<usize> = const { Cell::new(MATCH_CACHE_DEFAULT_ENTRIES) };
    static MPHF_LOOKUP_CACHE: RefCell<MphfLookupCache> = RefCell::new(MphfLookupCache::default());
    static CACHE_CONTEXT: Cell<Option<(u64, usize, usize)>> = const { Cell::new(None) };
    static BLOCK_LOOKUP_CACHE: RefCell<BlockLookupCache> =
        RefCell::new(BlockLookupCache::new());
}

#[inline]
pub(crate) fn ensure_cache_context(index: &BifrostIndex) {
    ensure_cache_context_key((index.cache_id, index.k, index.g));
}

#[inline]
fn ensure_cache_context_key(context: (u64, usize, usize)) {
    CACHE_CONTEXT.with(|owner| {
        let previous = owner.get();
        if previous != Some(context) {
            if previous.is_some() {
                MPHF_LOOKUP_CACHE
                    .with(|cache| cache.borrow_mut().entries.fill(MphfCacheEntry::default()));
                MPHF_LOOKUP_CACHE.with(|cache| {
                    cache
                        .borrow_mut()
                        .matches
                        .entries
                        .fill(FastKmerMatchCacheEntry::default());
                });
                BLOCK_LOOKUP_CACHE.with(|cache| {
                    cache
                        .borrow_mut()
                        .entries
                        .fill(BlockLookupCacheEntry::default())
                });
                MPHF_LOOKUP_CACHE.with(|cache| {
                    cache
                        .borrow_mut()
                        .minimizers
                        .entries
                        .fill(MinimizerCandidateCacheEntry::default());
                });
            }
            owner.set(Some(context));
        }
    });
}

pub(crate) fn reset_thread_local_caches() {
    MPHF_LOOKUP_CACHE.with(|cache| *cache.borrow_mut() = MphfLookupCache::default());
    BLOCK_LOOKUP_CACHE.with(|tl| *tl.borrow_mut() = BlockLookupCache::new());
}

pub(crate) fn with_mphf_lookup_cache<T>(f: impl FnOnce(&mut MphfLookupCache) -> T) -> T {
    MPHF_LOOKUP_CACHE.with(|cache| f(&mut cache.borrow_mut()))
}

#[cfg(test)]
pub(crate) fn read_fast_match_cache(key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
    MPHF_LOOKUP_CACHE.with(|tl| tl.borrow().matches.get(key, flags))
}

#[cfg(test)]
pub(crate) fn write_fast_match_cache(key: u64, flags: u8, value: Option<FastKmerMatch>) {
    MPHF_LOOKUP_CACHE.with(|tl| tl.borrow_mut().matches.put(key, flags, value));
}

pub(crate) fn minimizer_hex(bytes: [u8; 8]) -> String {
    let mut out = String::with_capacity(16);
    for b in bytes {
        let _ = write!(out, "{b:02x}");
    }
    out
}

#[inline]
pub(crate) fn block_index_for_position_fast(
    index: &BifrostIndex,
    unitig_id: usize,
    pos: usize,
) -> Option<usize> {
    if let Some(cached) = BLOCK_LOOKUP_CACHE.with(|tl| tl.borrow().get(unitig_id, pos)) {
        return Some(cached);
    }
    index
        .flat_ec
        .block_index_for_position(unitig_id, pos)
        .inspect(|&block_idx| {
            if let Some((lb, ub)) = index.flat_ec.block_bounds(unitig_id, block_idx) {
                BLOCK_LOOKUP_CACHE
                    .with(|tl| tl.borrow_mut().put(unitig_id, pos, block_idx, lb, ub));
            }
        })
}

#[inline]
pub(crate) fn ec_slice(index: &BifrostIndex, unitig_id: usize, block_idx: usize) -> &[u32] {
    index.flat_ec.ec(unitig_id, block_idx)
}

#[inline]
pub(crate) fn onlist_cardinality(onlist: &[bool]) -> u32 {
    onlist.iter().filter(|&&v| v).count() as u32
}

#[inline]
pub(crate) fn fast_path_enabled(
    index: &BifrostIndex,
    dbg: bool,
    options: PseudoalignOptions,
) -> bool {
    if env_flag(&DISABLE_FAST_PATH, "KALLISTORS_DISABLE_FAST_PATH") {
        return false;
    }
    (!dbg || options.investigation.trace_fast_path)
        && !index.use_shade
        && !options.do_union
        && !options.dfk_onlist
        && !options.no_jump
        && !options.kallisto_enum
        && !options.kallisto_strict
        && !options.kallisto_local_fallback
        && !options.kallisto_fallback
        && !options.discard_special_only
        && !options.kallisto_direct_kmer
        && !options.kallisto_bifrost_find
}

#[inline]
fn mphf_lookup_cached(
    index: &BifrostIndex,
    cache: &mut MphfLookupCache,
    min_bytes: &[u8; 8],
) -> Option<(usize, u16)> {
    // Canonical minimizers are left-aligned, so their low padding bits carry no entropy.
    let word = u64::from_le_bytes(*min_bytes);
    let slot = (word.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 53) as usize;
    let entry = &cache.entries[slot];
    let state = entry.value >> 48;
    if state != 0 && entry.key == *min_bytes {
        return (state == 1).then_some((entry.value as u32 as usize, (entry.value >> 32) as u16));
    }
    let idx = index
        .mphf
        .lookup(min_bytes)
        .and_then(|v| usize::try_from(v).ok());
    let fingerprint = idx
        .and_then(|i| index.minz_fingerprints.get(i))
        .copied()
        .unwrap_or(0);
    cache.entries[slot] = MphfCacheEntry {
        key: *min_bytes,
        value: u64::from(idx.unwrap_or(0) as u32)
            | (u64::from(fingerprint) << 32)
            | ((if idx.is_some() { 1u64 } else { 2u64 }) << 48),
    };
    idx.map(|idx| (idx, fingerprint))
}

#[inline]
pub(crate) fn minimizer_bucket_idx(
    index: &BifrostIndex,
    cache: &mut MphfLookupCache,
    min_bytes: &[u8; 8],
) -> Option<usize> {
    mphf_lookup_cached(index, cache, min_bytes).map(|(idx, _)| idx)
}

#[inline]
pub(crate) fn minimizer_bucket_with_fingerprint(
    index: &BifrostIndex,
    cache: &mut MphfLookupCache,
    key: &[u8; 8],
) -> Option<(usize, u16)> {
    mphf_lookup_cached(index, cache, key)
}

#[inline]
pub(crate) fn minimizer_candidates_cached_into_with_code(
    cache: &mut MphfLookupCache,
    kmer: &[u8],
    g: usize,
    fwd_code: u64,
    rev_code: u64,
    out: &mut Vec<([u8; 8], usize)>,
    cache_disabled: bool,
) -> bool {
    if !cache_disabled {
        let hit = cache.minimizers.get_into(fwd_code, out);
        if hit {
            return true;
        }
    }
    let valid = if kmer.len() == 31 && g == 23 {
        cache.primary_cursor.fill(fwd_code, rev_code, out)
    } else {
        minimizers_for_encoded_kmer_into(kmer, g, fwd_code, rev_code, out)
    };
    if !valid {
        return false;
    }
    if !cache_disabled {
        cache.minimizers.put(fwd_code, out);
    }
    true
}

#[inline]
pub(crate) fn minimizer_candidates_single_cell(
    cache: &mut MphfLookupCache,
    kmer: &[u8],
    g: usize,
    forward: u64,
    reverse: u64,
    out: &mut Vec<([u8; 8], usize)>,
) -> bool {
    if kmer.len() == 31 && g == 23 {
        cache.primary_cursor.fill_bus(forward, reverse, out)
    } else {
        minimizers_for_encoded_kmer_into(kmer, g, forward, reverse, out)
    }
}

#[inline]
pub(crate) fn minimizer_bucket_has_overcrowded_marker(
    index: &BifrostIndex,
    cache: &mut MphfLookupCache,
    min_bytes: &[u8; 8],
) -> bool {
    let Some(bucket_idx) = minimizer_bucket_idx(index, cache, min_bytes) else {
        return false;
    };
    index.minz_positions.get_with_overcrowded(bucket_idx).1
}
