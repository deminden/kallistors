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
    fingerprint: u16,
}

pub(crate) struct MphfLookupCache {
    entries: Vec<MphfCacheEntry>,
    matches: FastKmerMatchCache,
    minimizers: MinimizerCandidateCache,
    primary_cursor: super::super::minimizers::PrimaryMinimizerCursor,
    query_fingerprints: super::super::fingerprints::QueryFingerprintCache,
    pub(super) scan_scratch: super::fast::FastScannerScratch,
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
            scan_scratch: super::fast::FastScannerScratch::default(),
            slow_rev_buf: Vec::new(),
        }
    }
}

impl MphfLookupCache {
    #[inline]
    pub(crate) fn match_get(&self, key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
        self.matches.get(key, flags)
    }
    #[inline]
    pub(crate) fn match_put(&mut self, key: u64, flags: u8, value: Option<FastKmerMatch>) {
        self.matches.put(key, flags, value);
    }

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
    // Singletons use the canonical word's eight zero padding bits for a tag
    // and a seven-bit position. Other words/positions are recomputed exactly.
    // Zero = invalid; one = empty; 2..=8 = complete sparse tie list.
    value: u64,
}

struct MinimizerCandidateCache {
    entries: Vec<MinimizerCandidateCacheEntry>,
    ties: std::collections::HashMap<usize, Vec<([u8; 8], usize)>>,
}

impl MinimizerCandidateCache {
    fn new() -> Self {
        Self {
            entries: vec![MinimizerCandidateCacheEntry::default(); 1 << 16],
            ties: std::collections::HashMap::new(),
        }
    }

    #[inline]
    fn slot(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as usize & (self.entries.len() - 1)
    }

    #[inline]
    fn get_into(&self, key: u64, out: &mut Vec<([u8; 8], usize)>) -> bool {
        let slot = self.slot(key);
        let entry = &self.entries[slot];
        if entry.value == 0 || entry.key != key {
            return false;
        }
        if entry.value & 128 != 0 {
            out.clear();
            out.push((
                (entry.value & !255).to_le_bytes(),
                (entry.value & 127) as usize,
            ));
            true
        } else if entry.value == 1 {
            out.clear();
            true
        } else if let Some(ties) = self.ties.get(&slot) {
            out.clear();
            out.extend_from_slice(ties);
            true
        } else {
            false
        }
    }

    #[cfg(test)]
    fn get(&self, key: u64) -> Option<Vec<([u8; 8], usize)>> {
        let mut out = Vec::new();
        self.get_into(key, &mut out).then_some(out)
    }

    #[inline]
    fn put(&mut self, key: u64, values: &[([u8; 8], usize)]) {
        let slot = self.slot(key);
        let value = match values {
            [] => 1,
            [(bytes, position)] => {
                let word = u64::from_le_bytes(*bytes);
                if word & 255 != 0 || *position > 127 {
                    self.entries[slot].value = 0;
                    return;
                }
                word | 128 | *position as u64
            }
            ties if ties.len() <= 8 => {
                let cached = self.ties.entry(slot).or_default();
                cached.clear();
                cached.extend_from_slice(ties);
                ties.len() as u64
            }
            _ => {
                self.entries[slot].value = 0;
                return;
            }
        };
        self.entries[slot] = MinimizerCandidateCacheEntry { key, value };
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
    unitig_id: u32,
    start: u32,
    block_idx: u32,
    flags: u8,
    state: u8,
    match_bits: u8,
    negative_flags: u8,
    negative_key: u64,
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
        (key.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as usize & (self.entries.len() - 1)
    }

    #[inline]
    fn get(&self, key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
        let entry = &self.entries[self.slot(key)];
        if entry.state == 1 && entry.key == key && entry.flags == flags {
            return Some(Some(FastKmerMatch {
                unitig_id: entry.unitig_id as usize,
                start: entry.start as usize,
                block_idx: entry.block_idx as usize,
                used_revcomp: entry.match_bits & 1 != 0,
                forward_strand: entry.match_bits & 2 != 0,
                matched_relaxed: entry.match_bits & 4 != 0,
            }));
        }
        if entry.match_bits & 8 != 0 && entry.negative_key == key && entry.negative_flags == flags {
            return Some(None);
        }
        None
    }

    #[inline]
    fn put(&mut self, key: u64, flags: u8, value: Option<FastKmerMatch>) {
        let slot = self.slot(key);
        let entry = &mut self.entries[slot];
        if let Some(value) = value {
            if entry.negative_key == key && entry.negative_flags == flags {
                entry.match_bits &= !8;
            }
            let (Ok(unitig_id), Ok(start), Ok(block_idx)) = (
                u32::try_from(value.unitig_id),
                u32::try_from(value.start),
                u32::try_from(value.block_idx),
            ) else {
                entry.state = 0;
                return;
            };
            entry.key = key;
            entry.flags = flags;
            entry.state = 1;
            entry.unitig_id = unitig_id;
            entry.start = start;
            entry.block_idx = block_idx;
            entry.match_bits = (entry.match_bits & 8)
                | u8::from(value.used_revcomp)
                | (u8::from(value.forward_strand) << 1)
                | (u8::from(value.matched_relaxed) << 2);
        } else {
            if entry.key == key && entry.flags == flags {
                entry.state = 0;
            }
            entry.negative_key = key;
            entry.negative_flags = flags;
            entry.match_bits |= 8;
        }
    }
}

thread_local! {
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
    if entry.state != 0 && entry.key == *min_bytes {
        return (entry.state == 1).then_some((entry.idx as usize, entry.fingerprint));
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
        idx: idx.unwrap_or(0) as u32,
        state: if idx.is_some() { 1 } else { 2 },
        fingerprint,
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

#[cfg(test)]
mod cache_context_tests {
    use super::*;

    #[test]
    fn compact_match_entries_keep_all_bits_and_reject_wide_coordinates() {
        assert_eq!(std::mem::size_of::<FastKmerMatchCacheEntry>(), 32);
        let mut cache = FastKmerMatchCache::new();
        for bits in 0..8u8 {
            let value = FastKmerMatch {
                unitig_id: u32::MAX as usize,
                start: 71,
                block_idx: 19,
                used_revcomp: bits & 1 != 0,
                forward_strand: bits & 2 != 0,
                matched_relaxed: bits & 4 != 0,
            };
            cache.put(123, 7, Some(value));
            let got = cache.get(123, 7).unwrap().unwrap();
            assert_eq!(
                (got.unitig_id, got.start, got.block_idx),
                (value.unitig_id, value.start, value.block_idx)
            );
            assert_eq!(
                (got.used_revcomp, got.forward_strand, got.matched_relaxed),
                (
                    value.used_revcomp,
                    value.forward_strand,
                    value.matched_relaxed
                )
            );
            assert!(cache.get(123, 6).is_none());
        }
        let collision = (124..u64::MAX)
            .find(|&key| cache.slot(key) == cache.slot(123))
            .unwrap();
        cache.put(collision, 255, None);
        assert!(cache.get(123, 7).unwrap().is_some());
        assert!(matches!(cache.get(collision, 255), Some(None)));
        assert!(cache.get(collision, 15).is_none());
        cache.put(123, 7, None);
        assert!(matches!(cache.get(123, 7), Some(None)));
        if let Some(wide) = (u32::MAX as usize).checked_add(1) {
            for which in 0..3 {
                let mut value = FastKmerMatch::default();
                match which {
                    0 => value.unitig_id = wide,
                    1 => value.start = wide,
                    _ => value.block_idx = wide,
                }
                cache.put(123, 7, Some(value));
                assert!(cache.get(123, 7).is_none());
            }
        }
    }

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
        MPHF_LOOKUP_CACHE.with(|cache| cache.borrow_mut().minimizers.put(42, &[([0; 8], 3)]));
        BLOCK_LOOKUP_CACHE.with(|cache| cache.borrow_mut().put(7, 3, 0, 0, 10));
        ensure_cache_context_key((100, 31, 23));
        assert_eq!(
            MPHF_LOOKUP_CACHE.with(|cache| cache.borrow().minimizers.get(42)),
            Some(vec![([0; 8], 3)])
        );
        assert_eq!(
            read_fast_match_cache(42, 15).flatten().unwrap().unitig_id,
            7
        );
        ensure_cache_context_key((101, 31, 23));
        assert!(read_fast_match_cache(42, 15).is_none());
        assert!(MPHF_LOOKUP_CACHE.with(|cache| cache.borrow().minimizers.get(42).is_none()));
        assert!(BLOCK_LOOKUP_CACHE.with(|cache| cache.borrow().get(7, 3).is_none()));
        write_fast_match_cache(42, 15, None);
        ensure_cache_context_key((101, 31, 21));
        assert!(read_fast_match_cache(42, 15).is_none());
    }

    #[test]
    fn candidate_cache_never_truncates_a_tie_list() {
        assert_eq!(std::mem::size_of::<MinimizerCandidateCacheEntry>(), 16);
        let mut cache = MinimizerCandidateCache::new();
        for word in [0, 0x8000_0000_0000_0000u64, 0xffff_ffff_ffff_ff00] {
            for position in [0, 1, 31, 127] {
                cache.put(40, &[(word.to_le_bytes(), position)]);
                assert_eq!(cache.get(40), Some(vec![(word.to_le_bytes(), position)]));
            }
        }
        for (word, position) in [(0u64, 128), (1, 31), (255, 127), (0, usize::MAX)] {
            cache.put(40, &[(word.to_le_bytes(), position)]);
            assert!(cache.get(40).is_none());
        }
        cache.put(41, &[]);
        assert_eq!(cache.get(41), Some(vec![]));
        cache.put(42, &[([1; 8], usize::MAX); 8]);
        cache.put(42, &[([1; 8], 3); 9]);
        assert!(cache.get(42).is_none());
        cache.put(43, &[([1; 8], usize::MAX); 8]);
        assert_eq!(cache.get(43), Some(vec![([1; 8], usize::MAX); 8]));
        let single = (0x1234_5600u64.to_le_bytes(), 5);
        cache.put(43, &[single]);
        assert_eq!(cache.get(43), Some(vec![single]));
        let collision = (44..u64::MAX)
            .find(|&key| cache.slot(key) == cache.slot(43))
            .unwrap();
        cache.put(collision, &[(0x6789_ab00u64.to_le_bytes(), 7)]);
        assert!(cache.get(43).is_none());
        assert_eq!(
            cache.get(collision),
            Some(vec![(0x6789_ab00u64.to_le_bytes(), 7)])
        );
    }
}
