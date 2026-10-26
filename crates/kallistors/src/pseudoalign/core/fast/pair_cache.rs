//! Exact immutable annotation-pair results owned by one mapping invocation.
use super::BifrostIndex;

const MINIMUM: usize = 128;
const MAX_OWNED: usize = 4096;
const METADATA_BUDGET: usize = 6 * 1024 * 1024;
const OWNED_BUDGET: usize = 64 * 1024 * 1024;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Key {
    a: usize,
    b: usize,
    a_len: u32,
    b_len: u32,
}
impl Key {
    pub(super) fn eligible(a: &[u32], b: &[u32]) -> Option<Self> {
        if a.len().min(b.len()) < MINIMUM
            || std::ptr::eq(a, b)
            || a[0] > *b.last()?
            || b[0] > *a.last()?
        {
            return None;
        }
        let a = (a.as_ptr() as usize, u32::try_from(a.len()).ok()?);
        let b = (b.as_ptr() as usize, u32::try_from(b.len()).ok()?);
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        Some(Self {
            a: a.0,
            b: b.0,
            a_len: a.1,
            b_len: b.1,
        })
    }
    fn hash(self) -> usize {
        let mixed = ((self.a as u64) >> 3).wrapping_mul(0x9e37_79b9_7f4a_7c15)
            ^ ((self.b as u64) >> 3).rotate_left(27)
            ^ u64::from(self.a_len).wrapping_mul(0x85eb_ca6b)
            ^ u64::from(self.b_len).rotate_left(13);
        (mixed ^ (mixed >> 29)) as usize
    }
}

// Each64bytes of the hot plane hold two complete key/kind records. Only owned
// hits touch the separate answer-pointer plane; subset/empty hits need none.
#[derive(Default)]
#[repr(C)]
struct Entry {
    key: Key,
    kind: u8,
}
fn slot_metadata_bytes() -> usize {
    std::mem::size_of::<Entry>() + std::mem::size_of::<Option<Box<[u32]>>>()
}
#[derive(Clone, Copy)]
struct Plan {
    slots: usize,
    owned_bytes: usize,
}
impl Plan {
    fn for_workers(workers: usize) -> Option<Self> {
        let workers = workers.max(1);
        let limit = (METADATA_BUDGET / slot_metadata_bytes() / workers).min(65536);
        if limit < 2048 {
            return None;
        }
        let slots = 1usize << (usize::BITS - 1 - limit.leading_zeros());
        Some(Self {
            slots,
            owned_bytes: (OWNED_BUDGET / workers).min(64 * 1024 * 1024),
        })
    }
}
pub(super) struct PairData {
    plan: Plan,
    entries: Vec<Entry>,
    values: Vec<Option<Box<[u32]>>>,
    owned_bytes: usize,
}
impl PairData {
    fn new(plan: Plan) -> Self {
        Self {
            plan,
            entries: Vec::new(),
            values: Vec::new(),
            owned_bytes: 0,
        }
    }
    #[inline(always)]
    pub(super) fn apply<'index>(
        &self,
        key: Key,
        a: &'index [u32],
        b: &'index [u32],
        initial: &mut Option<&'index [u32]>,
        owned: &mut Vec<u32>,
    ) -> bool {
        let slot = key.hash() & (self.plan.slots - 1);
        let Some(entry) = self.entries.get(slot) else {
            return false;
        };
        if entry.key != key {
            return false;
        }
        debug_assert!(owned.is_empty());
        match entry.kind {
            1 => {
                *initial = Some(if a.len() <= b.len() { a } else { b });
            }
            2 => {
                *initial = None;
            }
            3 => {
                let Some(values) = self.values.get(slot).and_then(Option::as_deref) else {
                    return false;
                };
                owned.extend_from_slice(values);
                *initial = None;
            }
            _ => return false,
        }
        true
    }
    pub(super) fn remember(
        &mut self,
        key: Key,
        a: &[u32],
        b: &[u32],
        initial: Option<&[u32]>,
        owned: &[u32],
    ) {
        let subset = if let Some(value) = initial {
            if !std::ptr::eq(value, a) && !std::ptr::eq(value, b) {
                return;
            }
            true
        } else {
            false
        };
        let bytes = if subset {
            0
        } else {
            if owned.len() > MAX_OWNED {
                return;
            }
            owned.len() * 4
        };
        let slot = key.hash() & (self.plan.slots - 1);
        let old_bytes = self
            .values
            .get(slot)
            .and_then(Option::as_deref)
            .map_or(0, |v| v.len() * 4);
        let retained = self.owned_bytes - old_bytes;
        if bytes > self.plan.owned_bytes.saturating_sub(retained) {
            return;
        }
        if self.entries.is_empty() {
            self.entries.resize_with(self.plan.slots, Entry::default);
        }
        let kind = if subset {
            1
        } else if owned.is_empty() {
            2
        } else {
            3
        };
        if bytes != 0 && self.values.is_empty() {
            self.values.resize_with(self.plan.slots, || None);
        }
        if !self.values.is_empty() {
            self.values[slot] = if bytes != 0 { Some(owned.into()) } else { None };
        }
        self.entries[slot] = Entry { key, kind };
        self.owned_bytes = retained + bytes;
    }
}

/// A fresh immutable-index borrow is required for each invocation/worker.
/// Data is never shared across mutable-between-call index changes.
pub(crate) struct ScopedEcPairs<'index> {
    index: &'index BifrostIndex,
    data: Option<Box<PairData>>,
    in_batch: bool,
}
impl<'index> ScopedEcPairs<'index> {
    pub(crate) fn new(index: &'index BifrostIndex, workers: usize) -> Self {
        Self {
            index,
            data: Plan::for_workers(workers).map(|p| Box::new(PairData::new(p))),
            in_batch: false,
        }
    }
    pub(super) fn enter(&mut self, index: &BifrostIndex, slot: &mut Option<Box<PairData>>) -> bool {
        // An unexpected owner always retains uncached computation. Never
        // reinterpret/dereference stored integer addresses as references.
        *slot = None;
        if self.in_batch || !std::ptr::eq(self.index, index) {
            return false;
        }
        std::mem::swap(&mut self.data, slot);
        self.in_batch = true;
        true
    }
    pub(super) fn leave(&mut self, slot: &mut Option<Box<PairData>>) {
        debug_assert!(self.in_batch && self.data.is_none());
        std::mem::swap(&mut self.data, slot);
        self.in_batch = false;
    }
}

#[cfg(test)]
mod tests {
    use super::super::OnlineIntersection;
    use super::*;

    fn check(cache: &mut PairData, a: &[u32], b: &[u32], dirty: &[u32]) {
        let mut expected_owned = dirty.to_vec();
        let mut expected = OnlineIntersection {
            initial: Some(a),
            owned: &mut expected_owned,
            pairs: None,
        };
        expected.intersect_original(b);
        let wanted = expected.to_vec();
        let borrowed = expected.initial.is_some();
        let mut actual_owned = dirty.to_vec();
        {
            let mut actual = OnlineIntersection {
                initial: Some(a),
                owned: &mut actual_owned,
                pairs: Some(cache),
            };
            actual.intersect(b);
            assert_eq!(&*actual, wanted.as_slice());
            assert_eq!(actual.initial.is_some(), borrowed);
            if borrowed {
                assert!(std::ptr::eq(
                    actual.initial.unwrap(),
                    expected.initial.unwrap()
                ));
            }
        }
        // The caller can mutate a materialized answer after the cache stores
        // it. Subsequent hits must not alias this writable read scratch.
        actual_owned.push(u32::MAX);
    }
    #[test]
    fn exact_pair_results_match_original_owned_borrowed_duplicate_and_empty_answers() {
        let a: Vec<_> = (0..256).collect();
        let b: Vec<_> = (64..320).collect();
        let equal = a.clone();
        let even: Vec<_> = (0..512).step_by(2).collect();
        let odd: Vec<_> = (1..512).step_by(2).collect();
        let duplicates: Vec<_> = (0..384).map(|i| i / 3).collect();
        let other_duplicates: Vec<_> = (64..448).map(|i| i / 2).collect();
        let large: Vec<_> = (0..(MAX_OWNED * 4) as u32).collect();
        let other_large: Vec<_> = ((MAX_OWNED * 2) as u32..(MAX_OWNED * 6) as u32).collect();
        let fixtures: [(&[u32], &[u32]); 9] = [
            (&a, &b),
            (&a, &equal),
            (&a, &a[..128]),
            (&even, &odd),
            (&duplicates, &other_duplicates),
            (&a[..128], &b[..128]),
            (&a, &[1, 4]),
            (&large, &other_large),
            (&a[64..], &a),
        ];
        let mut cache = PairData::new(Plan {
            slots: 16,
            owned_bytes: 4096,
        });
        // This input exceeds the actual128-ID admission threshold. Verify
        // the integrated wrapper remembered a usable owned answer; a test
        // that silently stayed uncached would miss a broken accelerator.
        check(&mut cache, &a, &b, &[]);
        let key = Key::eligible(&a, &b).unwrap();
        let mut initial = Some(a.as_slice());
        let mut owned = Vec::new();
        assert!(cache.apply(key, &a, &b, &mut initial, &mut owned));
        assert!(initial.is_none());
        assert_eq!(owned, (64..256).collect::<Vec<_>>());
        for (left, right) in fixtures {
            for _ in 0..3 {
                check(&mut cache, left, right, &[]);
                check(&mut cache, right, left, &[]);
            }
            check(&mut cache, left, right, &[19, 21]);
            assert!(cache.owned_bytes <= cache.plan.owned_bytes);
        }
    }
    #[test]
    fn budgets_keys_eviction_and_oversize_fallback_are_bounded() {
        for workers in [0usize, 1, 2, 3, 8, 16, 64, 100, usize::MAX] {
            if let Some(plan) = Plan::for_workers(workers) {
                let workers = workers.max(1);
                assert!(plan.slots.is_power_of_two() && (2048..=65536).contains(&plan.slots));
                assert!(plan.slots * slot_metadata_bytes() * workers <= METADATA_BUDGET);
                assert!(plan.owned_bytes * workers <= OWNED_BUDGET);
            }
        }
        let a: Vec<_> = (0..256).collect();
        let b: Vec<_> = (64..320).collect();
        let key = Key::eligible(&a, &b).unwrap();
        assert_eq!(key, Key::eligible(&b, &a).unwrap());
        assert_ne!(key, Key::eligible(&a[..255], &b).unwrap());
        let mut cache = PairData::new(Plan {
            slots: 1,
            owned_bytes: 1024,
        });
        let output: Vec<_> = (64..256).collect();
        cache.remember(key, &a, &b, None, &output);
        assert_eq!(cache.owned_bytes, 768);
        let other = Key::eligible(&a[..255], &b).unwrap();
        cache.remember(other, &a[..255], &b, None, &[7; 384]);
        assert_eq!(cache.entries[0].key, key);
        assert_eq!(cache.owned_bytes, 768);
        cache.remember(other, &a[..255], &b, None, &[7; MAX_OWNED + 1]);
        assert_eq!(cache.entries[0].key, key);
        cache.remember(other, &a[..255], &b, None, &[7; 16]);
        assert_eq!(cache.owned_bytes, 64);
        assert_eq!(cache.entries[0].key, other);
        cache.remember(key, &a, &b, Some(&a), &[]);
        assert_eq!(cache.owned_bytes, 0);
        assert_eq!(cache.entries[0].key, key);
        // Give this valid oversized intersection enough payload allowance:
        // refusal must exercise MAX_OWNED, not the separate byte-budget gate.
        let big_a: Vec<_> = (0..(MAX_OWNED + 2) as u32).collect();
        let big_b: Vec<_> = (1..(MAX_OWNED + 3) as u32).collect();
        let oversized: Vec<_> = (1..(MAX_OWNED + 2) as u32).collect();
        let big_key = Key::eligible(&big_a, &big_b).unwrap();
        let mut roomy = PairData::new(Plan {
            slots: 1,
            owned_bytes: oversized.len() * 4 + 4,
        });
        roomy.remember(big_key, &big_a, &big_b, None, &oversized);
        assert!(roomy.entries.is_empty() && roomy.values.is_empty());
        assert_eq!(roomy.owned_bytes, 0);
    }
    #[test]
    fn invocation_owner_swap_and_new_scope_never_reuse_old_answers() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("reference.fa");
        let path = dir.path().join("reference.idx");
        std::fs::write(&fasta, ">tx\nCGGGTGTGGTGGCGGGCACCTGTAGTCCCAG\n").unwrap();
        crate::index::build_index(&path, &[fasta], Default::default()).unwrap();
        let index = crate::pseudoalign::build_bifrost_index(&path).unwrap();
        let other = crate::pseudoalign::build_bifrost_index(&path).unwrap();
        let a: Vec<_> = (0..256).collect();
        let b: Vec<_> = (64..320).collect();
        let key = Key::eligible(&a, &b).unwrap();
        let mut slot = None;
        {
            let mut scope = ScopedEcPairs::new(&index, 1);
            assert!(!scope.enter(&other, &mut slot));
            assert!(slot.is_none());
            assert!(scope.enter(&index, &mut slot));
            assert!(scope.data.is_none());
            slot.as_deref_mut()
                .unwrap()
                .remember(key, &a, &b, None, &[17, 18]);
            scope.leave(&mut slot);
            assert!(slot.is_none());
            assert!(!scope.data.as_ref().unwrap().entries.is_empty());
        }
        let scope = ScopedEcPairs::new(&index, 1);
        assert!(scope.data.as_ref().unwrap().entries.is_empty());
    }
}
