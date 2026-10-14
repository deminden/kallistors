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
    idx: u32,
    state: u8,
}

pub(crate) struct MphfLookupCache {
    entries: [MphfCacheEntry; 128],
}

impl Default for MphfLookupCache {
    fn default() -> Self {
        Self {
            entries: [MphfCacheEntry::default(); 128],
        }
    }
}

#[derive(Clone, Copy, Default)]
struct BlockLookupCacheEntry {
    unitig_id: usize,
    block_idx: u32,
    lb: u32,
    ub: u32,
    state: u8,
}

struct BlockLookupCache {
    entries: [BlockLookupCacheEntry; 128],
}

impl BlockLookupCache {
    fn new() -> Self {
        Self {
            entries: [BlockLookupCacheEntry::default(); 128],
        }
    }

    #[inline]
    fn slot(unitig_id: usize, pos: usize) -> usize {
        (unitig_id ^ pos) & 127
    }

    #[inline]
    fn get(&self, unitig_id: usize, pos: usize) -> Option<usize> {
        let entry = &self.entries[Self::slot(unitig_id, pos)];
        let pos = pos as u32;
        if entry.state != 0 && entry.unitig_id == unitig_id && pos >= entry.lb && pos < entry.ub {
            return Some(entry.block_idx as usize);
        }
        None
    }

    #[inline]
    fn put(&mut self, unitig_id: usize, pos: usize, block_idx: usize, lb: u32, ub: u32) {
        self.entries[Self::slot(unitig_id, pos)] = BlockLookupCacheEntry {
            unitig_id,
            block_idx: block_idx as u32,
            lb,
            ub,
            state: 1,
        };
    }
}

#[derive(Clone, Copy, Default)]
struct MinimizerCandidateCacheEntry {
    key: u64,
    count: u8,
    state: u8,
    values: [([u8; 8], usize); 8],
}

struct MinimizerCandidateCache {
    entries: Vec<MinimizerCandidateCacheEntry>,
}

impl MinimizerCandidateCache {
    fn new() -> Self {
        Self {
            entries: vec![MinimizerCandidateCacheEntry::default(); 1 << 16],
        }
    }

    #[inline]
    fn slot(&self, key: u64) -> usize {
        (key ^ (key >> 32)) as usize & (self.entries.len() - 1)
    }

    #[inline]
    fn get(&self, key: u64) -> Option<&[([u8; 8], usize)]> {
        let entry = &self.entries[self.slot(key)];
        if entry.state != 0 && entry.key == key {
            return Some(&entry.values[..entry.count as usize]);
        }
        None
    }

    #[inline]
    fn put(&mut self, key: u64, values: &[([u8; 8], usize)]) {
        // More ties than the inline capacity must be recomputed in full.
        if values.len() > 8 {
            let slot = self.slot(key);
            self.entries[slot].state = 0;
            return;
        }
        let mut entry = MinimizerCandidateCacheEntry {
            key,
            count: values.len().min(8) as u8,
            state: 1,
            values: [([0u8; 8], 0usize); 8],
        };
        for (dst, src) in entry.values.iter_mut().zip(values.iter().copied()) {
            *dst = src;
        }
        let slot = self.slot(key);
        self.entries[slot] = entry;
    }
}

#[derive(Clone, Copy, Default)]
pub(crate) struct FastKmerMatch {
    pub(crate) unitig_id: usize,
    pub(crate) start: usize,
    pub(crate) block_idx: usize,
    pub(crate) used_revcomp: bool,
    pub(crate) forward_strand: bool,
    pub(crate) matched_relaxed: bool,
}

#[derive(Clone, Copy, Default)]
struct FastKmerMatchCacheEntry {
    key: u64,
    value: FastKmerMatch,
    flags: u8,
    state: u8,
}

struct FastKmerMatchCache {
    entries: Vec<FastKmerMatchCacheEntry>,
}

impl FastKmerMatchCache {
    fn new() -> Self {
        Self {
            entries: vec![FastKmerMatchCacheEntry::default(); 1 << 14],
        }
    }

    #[inline]
    fn slot(&self, key: u64) -> usize {
        (key ^ (key >> 32)) as usize & (self.entries.len() - 1)
    }

    #[inline]
    fn get(&self, key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
        let entry = &self.entries[self.slot(key)];
        if entry.state != 0 && entry.key == key && entry.flags == flags {
            return Some((entry.state == 1).then_some(entry.value));
        }
        None
    }

    #[inline]
    fn put(&mut self, key: u64, flags: u8, value: Option<FastKmerMatch>) {
        let state = if value.is_some() { 1 } else { 2 };
        let slot = self.slot(key);
        let entry = &mut self.entries[slot];
        entry.key = key;
        entry.flags = flags;
        entry.state = state;
        entry.value = value.unwrap_or_default();
    }
}

thread_local! {
    static MPHF_LOOKUP_CACHE: RefCell<MphfLookupCache> = RefCell::new(MphfLookupCache::default());
    static CACHE_CONTEXT: Cell<Option<(u64, usize, usize)>> = const { Cell::new(None) };
    static FAST_KMER_MATCH_CACHE: RefCell<FastKmerMatchCache> =
        RefCell::new(FastKmerMatchCache::new());
    static BLOCK_LOOKUP_CACHE: RefCell<BlockLookupCache> =
        RefCell::new(BlockLookupCache::new());
    static MINIMIZER_CANDIDATE_CACHE: RefCell<MinimizerCandidateCache> =
        RefCell::new(MinimizerCandidateCache::new());
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
                FAST_KMER_MATCH_CACHE.with(|cache| {
                    cache
                        .borrow_mut()
                        .entries
                        .fill(FastKmerMatchCacheEntry::default())
                });
                BLOCK_LOOKUP_CACHE.with(|cache| {
                    cache
                        .borrow_mut()
                        .entries
                        .fill(BlockLookupCacheEntry::default())
                });
                MINIMIZER_CANDIDATE_CACHE.with(|cache| {
                    cache
                        .borrow_mut()
                        .entries
                        .fill(MinimizerCandidateCacheEntry::default())
                });
            }
            owner.set(Some(context));
        }
    });
}

pub(crate) fn reset_thread_local_caches() {
    MPHF_LOOKUP_CACHE.with(|cache| *cache.borrow_mut() = MphfLookupCache::default());
    FAST_KMER_MATCH_CACHE.with(|tl| *tl.borrow_mut() = FastKmerMatchCache::new());
    BLOCK_LOOKUP_CACHE.with(|tl| *tl.borrow_mut() = BlockLookupCache::new());
    MINIMIZER_CANDIDATE_CACHE.with(|tl| *tl.borrow_mut() = MinimizerCandidateCache::new());
}

pub(crate) fn with_mphf_lookup_cache<T>(f: impl FnOnce(&mut MphfLookupCache) -> T) -> T {
    MPHF_LOOKUP_CACHE.with(|cache| f(&mut cache.borrow_mut()))
}

pub(crate) fn read_fast_match_cache(key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
    FAST_KMER_MATCH_CACHE.with(|tl| tl.borrow().get(key, flags))
}

pub(crate) fn write_fast_match_cache(key: u64, flags: u8, value: Option<FastKmerMatch>) {
    FAST_KMER_MATCH_CACHE.with(|tl| tl.borrow_mut().put(key, flags, value));
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
) -> Option<usize> {
    let slot = usize::from(min_bytes[0] ^ min_bytes[7]) & (cache.entries.len() - 1);
    let entry = &cache.entries[slot];
    if entry.state != 0 && entry.key == *min_bytes {
        return (entry.state == 1).then_some(entry.idx as usize);
    }
    let idx = index
        .mphf
        .lookup(min_bytes)
        .and_then(|v| usize::try_from(v).ok());
    cache.entries[slot] = MphfCacheEntry {
        key: *min_bytes,
        idx: idx.unwrap_or(0) as u32,
        state: if idx.is_some() { 1 } else { 2 },
    };
    idx
}

#[inline]
pub(crate) fn minimizer_bucket_idx(
    index: &BifrostIndex,
    cache: &mut MphfLookupCache,
    min_bytes: &[u8; 8],
) -> Option<usize> {
    mphf_lookup_cached(index, cache, min_bytes)
}

#[inline]
pub(crate) fn minimizer_candidates_cached_into_with_code(
    kmer: &[u8],
    g: usize,
    fwd_code: u64,
    rev_code: u64,
    out: &mut Vec<([u8; 8], usize)>,
    cache_disabled: bool,
) -> bool {
    if !cache_disabled {
        let hit = MINIMIZER_CANDIDATE_CACHE.with(|tl| {
            let cache = tl.borrow();
            if let Some(cached) = cache.get(fwd_code) {
                out.clear();
                out.extend_from_slice(cached);
                true
            } else {
                false
            }
        });
        if hit {
            return true;
        }
    }
    if !minimizers_for_encoded_kmer_into(kmer, g, fwd_code, rev_code, out) {
        return false;
    }
    if !cache_disabled {
        MINIMIZER_CANDIDATE_CACHE.with(|tl| tl.borrow_mut().put(fwd_code, out));
    }
    true
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
    index
        .minz_positions
        .get(bucket_idx)
        .iter()
        .any(|&pos_id| (pos_id >> 32) as u32 == u32::MAX && (pos_id & 0x8000_0000) != 0)
}

#[cfg(test)]
mod cache_context_tests {
    use super::*;

    #[test]
    fn cache_entries_never_cross_index_or_encoding_context() {
        reset_thread_local_caches();
        ensure_cache_context_key((100, 31, 23));
        write_fast_match_cache(
            42,
            15,
            Some(FastKmerMatch {
                unitig_id: 7,
                ..FastKmerMatch::default()
            }),
        );
        MINIMIZER_CANDIDATE_CACHE.with(|cache| cache.borrow_mut().put(42, &[([1; 8], 3)]));
        BLOCK_LOOKUP_CACHE.with(|cache| cache.borrow_mut().put(7, 3, 0, 0, 10));
        ensure_cache_context_key((100, 31, 23));
        assert_eq!(
            read_fast_match_cache(42, 15).flatten().unwrap().unitig_id,
            7
        );
        ensure_cache_context_key((101, 31, 23));
        assert!(read_fast_match_cache(42, 15).is_none());
        assert!(MINIMIZER_CANDIDATE_CACHE.with(|cache| cache.borrow().get(42).is_none()));
        assert!(BLOCK_LOOKUP_CACHE.with(|cache| cache.borrow().get(7, 3).is_none()));
        write_fast_match_cache(42, 15, None);
        ensure_cache_context_key((101, 31, 21));
        assert!(read_fast_match_cache(42, 15).is_none());
    }

    #[test]
    fn candidate_cache_never_truncates_a_tie_list() {
        let mut cache = MinimizerCandidateCache::new();
        cache.put(42, &[([1; 8], 3); 8]);
        cache.put(42, &[([1; 8], 3); 9]);
        assert!(cache.get(42).is_none());
        cache.put(43, &[([1; 8], 3); 8]);
        assert_eq!(cache.get(43).unwrap().len(), 8);
    }
}
