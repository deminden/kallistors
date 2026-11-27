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
    /// When set, every worker reads and writes this one table instead of
    /// `entries`, so a k-mer resolved by any worker is a hit for all of them.
    pub(super) shared: Option<std::sync::Arc<SharedMatchTable>>,
}

/// One 32-byte slot holding the same positive and negative records as
/// [`FastKmerMatchCacheEntry`], packed into four words guarded by a seqlock.
///
/// - `w[0]`: positive key; `w[1]`: negative key;
/// - `w[2]`: unitig id (low 32 bits) and start (high 32 bits);
/// - `w[3]`: block index (bits 0-31), flags (32-35), negative flags (36-39),
///   match bits (40-43), positive state (44) and the sequence (45-63).
///
/// Writers take the slot by making the sequence odd with a CAS and skip the
/// insert when another writer holds it. Readers accept a record only if the
/// sequence was even and unchanged across their loads, so a hit is always a
/// complete record that some worker stored for exactly this key and flags.
#[repr(C, align(32))]
pub(crate) struct SharedMatchEntry {
    w: [std::sync::atomic::AtomicU64; 4],
}

const SEQ_ONE: u64 = 1 << 45;
const LOW_MASK: u64 = SEQ_ONE - 1;

pub(crate) struct SharedMatchTable {
    pointer: std::ptr::NonNull<SharedMatchEntry>,
    len: usize,
}

// SAFETY: the table only exposes atomics and frees its allocation on drop.
unsafe impl Send for SharedMatchTable {}
// SAFETY: every access goes through the per-word atomics.
unsafe impl Sync for SharedMatchTable {}

impl SharedMatchTable {
    pub(crate) fn new(entries: usize) -> crate::Result<std::sync::Arc<Self>> {
        let len = entries.clamp(1, 1 << 32);
        let layout = std::alloc::Layout::array::<SharedMatchEntry>(len)
            .map_err(|error| crate::Error::Io(std::io::Error::other(error.to_string())))?;
        // SAFETY: the layout is nonzero; all-zero words are a valid empty slot.
        let raw = unsafe { std::alloc::alloc_zeroed(layout) }.cast::<SharedMatchEntry>();
        let pointer = std::ptr::NonNull::new(raw).ok_or_else(|| {
            crate::Error::Io(std::io::Error::other(format!(
                "cannot allocate shared exact-match cache ({len} entries, {} bytes)",
                layout.size()
            )))
        })?;
        Ok(std::sync::Arc::new(Self { pointer, len }))
    }

    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    #[inline]
    fn slot(&self, key: u64) -> usize {
        let hash = key.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32;
        ((hash * self.len as u64) >> 32) as usize
    }

    #[inline]
    fn entry(&self, key: u64) -> &SharedMatchEntry {
        // SAFETY: slot() < len and the allocation lives as long as self.
        unsafe { &*self.pointer.as_ptr().add(self.slot(key)) }
    }

    #[inline]
    pub(super) fn slot_pointer(&self, key: u64) -> *const u8 {
        std::ptr::from_ref(self.entry(key)).cast()
    }

    #[inline]
    fn get(&self, key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
        use std::sync::atomic::{Ordering, fence};
        let entry = self.entry(key);
        let before = entry.w[3].load(Ordering::Acquire);
        if before & SEQ_ONE != 0 {
            return None;
        }
        let positive_key = entry.w[0].load(Ordering::Relaxed);
        let negative_key = entry.w[1].load(Ordering::Relaxed);
        let location = entry.w[2].load(Ordering::Relaxed);
        fence(Ordering::Acquire);
        if entry.w[3].load(Ordering::Relaxed) != before {
            return None;
        }
        let entry_flags = ((before >> 32) & 0xf) as u8;
        let negative_flags = ((before >> 36) & 0xf) as u8;
        let match_bits = ((before >> 40) & 0xf) as u8;
        if before & (1 << 44) != 0 && positive_key == key && entry_flags == flags {
            return Some(Some(FastKmerMatch {
                unitig_id: (location & 0xffff_ffff) as usize,
                start: (location >> 32) as usize,
                block_idx: (before & 0xffff_ffff) as usize,
                used_revcomp: match_bits & 1 != 0,
                forward_strand: match_bits & 2 != 0,
                matched_relaxed: match_bits & 4 != 0,
            }));
        }
        if match_bits & 8 != 0 && negative_key == key && negative_flags == flags {
            return Some(None);
        }
        None
    }

    #[inline]
    fn put(&self, key: u64, flags: u8, value: Option<FastKmerMatch>) {
        use std::sync::atomic::{Ordering, fence};
        debug_assert!(flags < 16);
        let entry = self.entry(key);
        let current = entry.w[3].load(Ordering::Relaxed);
        if current & SEQ_ONE != 0
            || entry.w[3]
                .compare_exchange(
                    current,
                    current.wrapping_add(SEQ_ONE),
                    Ordering::Acquire,
                    Ordering::Relaxed,
                )
                .is_err()
        {
            return;
        }
        fence(Ordering::Release);
        let mut positive_key = entry.w[0].load(Ordering::Relaxed);
        let mut negative_key = entry.w[1].load(Ordering::Relaxed);
        let mut location = entry.w[2].load(Ordering::Relaxed);
        let mut block_idx = current & 0xffff_ffff;
        let mut entry_flags = ((current >> 32) & 0xf) as u8;
        let mut negative_flags = ((current >> 36) & 0xf) as u8;
        let mut match_bits = ((current >> 40) & 0xf) as u8;
        let mut state = current & (1 << 44) != 0;
        if let Some(value) = value {
            if negative_key == key && negative_flags == flags {
                match_bits &= !8;
            }
            if let (Ok(unitig_id), Ok(start), Ok(block)) = (
                u32::try_from(value.unitig_id),
                u32::try_from(value.start),
                u32::try_from(value.block_idx),
            ) {
                positive_key = key;
                entry_flags = flags;
                state = true;
                location = u64::from(unitig_id) | (u64::from(start) << 32);
                block_idx = u64::from(block);
                match_bits = (match_bits & 8)
                    | u8::from(value.used_revcomp)
                    | (u8::from(value.forward_strand) << 1)
                    | (u8::from(value.matched_relaxed) << 2);
            } else {
                state = false;
            }
        } else {
            if positive_key == key && entry_flags == flags {
                state = false;
            }
            negative_key = key;
            negative_flags = flags;
            match_bits |= 8;
        }
        entry.w[0].store(positive_key, Ordering::Relaxed);
        entry.w[1].store(negative_key, Ordering::Relaxed);
        entry.w[2].store(location, Ordering::Relaxed);
        let low = block_idx
            | (u64::from(entry_flags & 0xf) << 32)
            | (u64::from(negative_flags & 0xf) << 36)
            | (u64::from(match_bits & 0xf) << 40)
            | (u64::from(state) << 44);
        let sequence = current.wrapping_add(2 * SEQ_ONE) & !LOW_MASK;
        entry.w[3].store(sequence | low, Ordering::Release);
    }
}

impl Drop for SharedMatchTable {
    fn drop(&mut self) {
        let layout = std::alloc::Layout::array::<SharedMatchEntry>(self.len).expect("layout");
        // SAFETY: allocated in new() with this exact layout.
        unsafe { std::alloc::dealloc(self.pointer.as_ptr().cast(), layout) };
    }
}

impl FastKmerMatchCache {
    pub(super) fn new() -> Self {
        Self {
            // Full PE working-set evidence supports this bounded table.
            // Complete query/flag keys still decide every positive/negative hit.
            entries: vec![FastKmerMatchCacheEntry::default(); MATCH_CACHE_ENTRIES.with(Cell::get)],
            shared: None,
        }
    }

    #[inline]
    pub(super) fn slot(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as usize & (self.entries.len() - 1)
    }

    #[inline]
    pub(super) fn get(&self, key: u64, flags: u8) -> Option<Option<FastKmerMatch>> {
        if let Some(shared) = &self.shared {
            return shared.get(key, flags);
        }
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
        if let Some(shared) = &self.shared {
            return shared.put(key, flags, value);
        }
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
            && self
                .matches
                .shared
                .as_ref()
                .map_or(self.matches.entries.len(), |shared| shared.len())
                >= (1 << 21)
    }

    /// Hint the slot of one forward k-mer word. Read-only.
    #[inline]
    pub(crate) fn prefetch_match_key(&self, forward: u64, canonical: u64) {
        if let Some(filter) = &self.kmer_filter {
            prefetch_match_cache_line(filter.block_pointer(canonical));
        }
        let pointer = if let Some(shared) = &self.matches.shared {
            shared.slot_pointer(forward)
        } else {
            let slot = self.matches.slot(forward);
            // SAFETY: slot masks into the live power-of-two table.
            unsafe { self.matches.entries.as_ptr().add(slot) }.cast::<u8>()
        };
        prefetch_match_cache_line(pointer);
    }

    /// Hint the slots of a read's first and last forward k-mer words.
    /// Read-only: no lookup, insertion or retained pointer.
    #[inline]
    pub(crate) fn prefetch_match_word_ends(&self, sequence: &[u8]) {
        if sequence.len() < 31 {
            return;
        }
        for pos in [0, sequence.len() - 31] {
            if let Some((forward, _)) = crate::dna::encode_kmer_pair(&sequence[pos..pos + 31]) {
                let pointer = if let Some(shared) = &self.matches.shared {
                    shared.slot_pointer(forward)
                } else {
                    let slot = self.matches.slot(forward);
                    // SAFETY: cache construction/scoped resizing ensures a nonzero
                    // power-of-two table; slot masks into this live allocation.
                    unsafe { self.matches.entries.as_ptr().add(slot) }.cast::<u8>()
                };
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
