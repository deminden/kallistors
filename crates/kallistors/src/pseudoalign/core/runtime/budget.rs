//! Scoped exact-match allocation and restoration.

use super::{FastKmerMatchCacheEntry, MATCH_CACHE_ENTRIES, MPHF_LOOKUP_CACHE};

// Each worker owns a private table. The coordinator resolves the aggregate
// allowance once; retention alone changes, and full keys/flags decide hits.
pub(super) const MATCH_CACHE_DEFAULT_ENTRIES: usize = 1 << 18;
const _: () = assert!(
    std::mem::size_of::<FastKmerMatchCacheEntry>()
        == crate::pseudoalign::cache_budget::MATCH_CACHE_ENTRY_BYTES
);

pub(crate) struct MatchCacheBudgetScope {
    previous: usize,
}

impl MatchCacheBudgetScope {
    pub(crate) fn enter(entries: usize) -> crate::Result<Self> {
        assert!(entries.is_power_of_two());
        // Bootstrap an uninitialized TLS cache with a small table. Larger
        // requested allocations use the fallible reserve below; smaller ones
        // do not briefly allocate the default table first.
        let previous = MATCH_CACHE_ENTRIES
            .with(|current| current.replace(entries.min(MATCH_CACHE_DEFAULT_ENTRIES)));
        if let Err(error) = resize_match_cache(entries) {
            // An unsuccessful reserve leaves an empty table. Restore the
            // previous usable capacity before returning to library callers.
            MATCH_CACHE_ENTRIES.with(|current| current.set(previous));
            resize_match_cache(previous)
                .expect("restore exact-match cache after allocation failure");
            return Err(error);
        }
        MATCH_CACHE_ENTRIES.with(|current| current.set(entries));
        Ok(Self { previous })
    }
}

impl Drop for MatchCacheBudgetScope {
    fn drop(&mut self) {
        MATCH_CACHE_ENTRIES.with(|current| current.set(self.previous));
        resize_match_cache(self.previous).expect("restore exact-match cache scope");
    }
}

fn resize_match_cache(entries: usize) -> crate::Result<()> {
    MPHF_LOOKUP_CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        if cache.matches.entries.len() != entries {
            // Slot placement changes with capacity. Release old pages before
            // reserving the replacement, then start with an empty exact cache.
            crate::index::scratch_pages::retire(std::mem::take(&mut cache.matches.entries));
            cache
                .matches
                .entries
                .try_reserve_exact(entries)
                .map_err(|error| {
                    crate::Error::Io(std::io::Error::other(format!(
                        "cannot allocate exact-match cache ({entries} entries, {} bytes): {error}",
                        entries * std::mem::size_of::<FastKmerMatchCacheEntry>()
                    )))
                })?;
            cache
                .matches
                .entries
                .resize(entries, FastKmerMatchCacheEntry::default());
        }
        Ok(())
    })
}
