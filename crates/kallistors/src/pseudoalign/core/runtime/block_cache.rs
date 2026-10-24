//! Exact annotation-block bounds for repeated unitig positions.

#[derive(Clone, Copy, Default)]
pub(super) struct BlockLookupCacheEntry {
    unitig_id: usize,
    block_idx: u32,
    lb: u32,
    ub: u32,
    state: u8,
}

pub(super) struct BlockLookupCache {
    pub(super) entries: [BlockLookupCacheEntry; 128],
}

impl BlockLookupCache {
    pub(super) fn new() -> Self {
        Self {
            entries: [BlockLookupCacheEntry::default(); 128],
        }
    }

    #[inline]
    fn slot(unitig_id: usize, pos: usize) -> usize {
        (unitig_id ^ pos) & 127
    }

    #[inline]
    pub(super) fn get(&self, unitig_id: usize, pos: usize) -> Option<usize> {
        let entry = &self.entries[Self::slot(unitig_id, pos)];
        let pos = pos as u32;
        if entry.state != 0 && entry.unitig_id == unitig_id && pos >= entry.lb && pos < entry.ub {
            return Some(entry.block_idx as usize);
        }
        None
    }

    #[inline]
    pub(super) fn put(&mut self, unitig_id: usize, pos: usize, block_idx: usize, lb: u32, ub: u32) {
        self.entries[Self::slot(unitig_id, pos)] = BlockLookupCacheEntry {
            unitig_id,
            block_idx: block_idx as u32,
            lb,
            ub,
            state: 1,
        };
    }
}
