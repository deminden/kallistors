//! Exact adaptive prefixes for narrow temporary graph-edge ranks.
use super::{KmerDirectory, KmerLookup};

pub(super) const BRANCH: u32 = 1 << 31;

pub(super) struct PrefixPlan {
    pub(super) roots: Box<[u32]>,
    pub(super) branches: Box<[[u32; 4]]>,
}

impl PrefixPlan {
    pub(super) fn new(
        entries: &[u64],
        lookup: &KmerLookup,
        k: usize,
        root_bases: usize,
        limit: u32,
    ) -> Self {
        assert!(entries.len() < BRANCH as usize && limit > 0 && root_bases > 0);
        let bits = 2 * k.min(root_bases);
        let mut builder = Builder {
            entries,
            lookup,
            k,
            limit: limit as usize,
            branches: Vec::new(),
        };
        let mut roots = Vec::with_capacity((1usize << bits) + 1);
        let mut first = 0;
        for prefix in 0..1u64 << bits {
            let end = builder.boundary(prefix + 1, bits, first, entries.len());
            roots.push(builder.split(first, end, prefix, bits));
            first = end;
        }
        roots.push(entries.len() as u32);
        Self {
            roots: roots.into_boxed_slice(),
            branches: builder.branches.into_boxed_slice(),
        }
    }
}

struct Builder<'a> {
    entries: &'a [u64],
    lookup: &'a KmerLookup,
    k: usize,
    limit: usize,
    branches: Vec<[u32; 4]>,
}

impl Builder<'_> {
    fn boundary(&self, prefix: u64, bits: usize, first: usize, end: usize) -> usize {
        if let KmerDirectory::Ranks(offsets) = &self.lookup.directory
            && self.lookup.shift <= 64 - bits
        {
            offsets.get((prefix as usize) << (64 - bits - self.lookup.shift))
        } else if prefix == 1u64 << bits {
            self.entries.len()
        } else {
            let lower = prefix << (64 - bits);
            first + self.entries[first..end].partition_point(|&word| word < lower)
        }
    }

    fn split(&mut self, first: usize, end: usize, prefix: u64, bits: usize) -> u32 {
        if end - first <= self.limit {
            return first as u32;
        }
        // Distinct complete k-mers always separate before the last base.
        assert!(bits < 2 * self.k);
        let branch = self.branches.len();
        assert!(branch < BRANCH as usize);
        self.branches.push([0; 4]);
        let mut child_first = first;
        for code in 0..4 {
            let child_prefix = (prefix << 2) | code as u64;
            let child_end = self.boundary(child_prefix + 1, bits + 2, child_first, end);
            let child = self.split(child_first, child_end, child_prefix, bits + 2);
            self.branches[branch][code] = child;
            child_first = child_end;
        }
        BRANCH | branch as u32
    }
}

#[inline]
pub(super) fn start(roots: &[u32], branches: &[[u32; 4]], shift: usize, word: u64) -> u32 {
    let mut node = roots[(word >> shift) as usize];
    let mut shift = shift;
    while node & BRANCH != 0 {
        shift -= 2;
        node = branches[(node & !BRANCH) as usize][((word >> shift) & 3) as usize];
    }
    node
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adaptive_leaves_preserve_exact_dense_sparse_and_empty_ranges() {
        for k in [3, 5, 15, 31] {
            let mut random = 71293u64;
            let mut entries = (0..2048)
                .map(|i| {
                    random ^= random << 13;
                    random ^= random >> 7;
                    random ^= random << 17;
                    let value = if i % 3 == 0 { random & 0xffff } else { random };
                    value & super::super::kmer_mask(k)
                })
                .collect::<Vec<_>>();
            entries.extend([0, super::super::kmer_mask(k)]);
            entries.sort_unstable();
            entries.dedup();
            let lookup = KmerLookup::with_threads(&entries, 1, k).unwrap();
            for root_bases in [1, 3] {
                for limit in [1, 3, 15, 127] {
                    let plan = PrefixPlan::new(&entries, &lookup, k, root_bases, limit);
                    let bits = 2 * k.min(root_bases);
                    let mut pending = plan.roots[..plan.roots.len() - 1]
                        .iter()
                        .enumerate()
                        .map(|(prefix, &node)| (node, prefix as u64, bits))
                        .collect::<Vec<_>>();
                    while let Some((node, prefix, bits)) = pending.pop() {
                        if node & BRANCH != 0 {
                            for (code, &child) in
                                plan.branches[(node & !BRANCH) as usize].iter().enumerate()
                            {
                                pending.push((child, (prefix << 2) | code as u64, bits + 2));
                            }
                        } else {
                            let first =
                                entries.partition_point(|&word| word >> (64 - bits) < prefix);
                            let end =
                                entries.partition_point(|&word| word >> (64 - bits) <= prefix);
                            assert_eq!(node as usize, first);
                            assert!(end - first <= limit as usize);
                            for &word in &entries[first..end] {
                                assert_eq!(
                                    start(
                                        &plan.roots,
                                        &plan.branches,
                                        64 - 2 * k.min(root_bases),
                                        word
                                    ),
                                    node
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
