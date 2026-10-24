//! Lossless cached minimizer candidate lists, including complete ties.

#[derive(Clone, Copy, Default)]
pub(super) struct MinimizerCandidateCacheEntry {
    key: u64,
    // Singletons use the canonical word's eight zero padding bits for a tag
    // and a seven-bit position. Other words/positions are recomputed exactly.
    // Zero = invalid; one = empty; 2..=8 = complete sparse tie list.
    value: u64,
}

pub(super) struct MinimizerCandidateCache {
    pub(super) entries: Vec<MinimizerCandidateCacheEntry>,
    ties: std::collections::HashMap<usize, Vec<([u8; 8], usize)>>,
}

impl MinimizerCandidateCache {
    pub(super) fn new() -> Self {
        Self {
            entries: Vec::new(),
            ties: std::collections::HashMap::new(),
        }
    }

    #[inline]
    pub(super) fn slot(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as usize & (self.entries.len() - 1)
    }

    #[inline]
    pub(super) fn get_into(&self, key: u64, out: &mut Vec<([u8; 8], usize)>) -> bool {
        if self.entries.is_empty() {
            return false;
        }
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
    pub(super) fn get(&self, key: u64) -> Option<Vec<([u8; 8], usize)>> {
        let mut out = Vec::new();
        self.get_into(key, &mut out).then_some(out)
    }

    #[inline]
    pub(super) fn put(&mut self, key: u64, values: &[([u8; 8], usize)]) {
        // BUS can bypass this cache entirely. Allocate its table only on an
        // actual insertion; bulk retains the same capacity and cache answers.
        if self.entries.is_empty() {
            self.entries
                .resize(1 << 16, MinimizerCandidateCacheEntry::default());
        }
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
