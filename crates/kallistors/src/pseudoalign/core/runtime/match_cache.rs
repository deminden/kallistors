//! Exact positive/negative k-mer records and read-only hardware hints.

use super::{MATCH_CACHE_ENTRIES, MphfLookupCache};
use std::cell::Cell;

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
pub(super) struct FastKmerMatchCacheEntry {
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

pub(super) struct FastKmerMatchCache {
    pub(super) entries: Vec<FastKmerMatchCacheEntry>,
}

impl FastKmerMatchCache {
    pub(super) fn new() -> Self {
        Self {
            // Full PE working-set evidence supports this bounded table.
            // Complete query/flag keys still decide every positive/negative hit.
            entries: vec![FastKmerMatchCacheEntry::default(); MATCH_CACHE_ENTRIES.with(Cell::get)],
        }
    }

    #[inline]
    pub(super) fn slot(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as usize & (self.entries.len() - 1)
    }

    #[inline]
    pub(super) fn get(&self, key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
        // Outer None means uncached; Some(None) is a confirmed miss. Positive
        // and negative records share a slot but retain independent exact keys.
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
    pub(super) fn put(&mut self, key: u64, flags: u8, value: Option<FastKmerMatch>) {
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
    pub(crate) fn prefetch_match_capacity(&self) -> bool {
        cfg!(any(target_arch = "aarch64", target_arch = "x86_64"))
            && self.matches.entries.len() >= (1 << 21)
    }

    #[inline]
    pub(crate) fn prefetch_initial_match_words(&self, left: &[u8], right: &[u8]) {
        for sequence in [left, right] {
            if let Some(prefix) = sequence.get(..31)
                && let Some((forward, _)) = crate::dna::encode_kmer_pair(prefix)
            {
                let slot = self.matches.slot(forward);
                // SAFETY: cache construction/scoped resizing ensures a nonzero
                // power-of-two table; slot masks into this live allocation.
                // The hint retains no pointer and performs no software lookup.
                let pointer = unsafe { self.matches.entries.as_ptr().add(slot) }.cast::<u8>();
                prefetch_match_cache_line(pointer);
            }
        }
    }
}

#[inline]
fn prefetch_match_cache_line(pointer: *const u8) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: the caller supplies a pointer into a live immutable table. PRFM
    // is a hardware read hint; it does not store data or retain the pointer.
    unsafe {
        std::arch::asm!("prfm pldl1keep, [{pointer}]",pointer=in(reg)pointer,options(readonly,nostack,preserves_flags));
    }
    #[cfg(target_arch = "x86_64")]
    // SAFETY: SSE is part of the x86-64 baseline and the address is valid.
    unsafe {
        std::arch::x86_64::_mm_prefetch(pointer.cast::<i8>(), std::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    let _ = pointer;
}
