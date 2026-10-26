//! Stable layout of independent EM transcript/EC components.

use super::{EmMultiEc, EmWork};

mod parallel;
mod scalar;

pub(super) use parallel::{step, tile as parallel_tile};
pub(super) use scalar::{bootstrap_tile, scalar_step, scalar_tile};

#[cfg(test)]
mod tests;

struct Components {
    parents: Vec<usize>,
    ranks: Vec<u8>,
}

impl Components {
    fn new(len: usize) -> Self {
        Self {
            parents: (0..len).collect(),
            ranks: vec![0; len],
        }
    }

    fn root(&mut self, mut tr: usize) -> usize {
        while self.parents[tr] != tr {
            self.parents[tr] = self.parents[self.parents[tr]];
            tr = self.parents[tr];
        }
        tr
    }

    fn join(&mut self, a: usize, b: usize) {
        let (mut a, mut b) = (self.root(a), self.root(b));
        if a == b {
            return;
        }
        if self.ranks[a] < self.ranks[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parents[b] = a;
        if self.ranks[a] == self.ranks[b] {
            self.ranks[a] += 1;
        }
    }
}

pub(super) struct EmComponent {
    fixed: bool,
    parallel: bool,
    last_changed: bool,
    // Owned by the certified scalar scheduler; other engines ignore it.
    scalar_finite: bool,
    // Each stable recursion interval stores its bit on its first component.
    // Ancestors share that component but use a different depth bit.
    fixed_subtrees: u64,
    transcripts: core::range::Range<usize>,
    ecs: core::range::Range<usize>,
}

#[inline]
fn witness(groups: &[EmComponent]) -> Option<usize> {
    groups
        .iter()
        .enumerate()
        .filter(|(_, group)| !group.fixed)
        .max_by_key(|(_, group)| (group.last_changed, group.ecs.end - group.ecs.start))
        .map(|(index, _)| index)
}

impl EmWork {
    pub(super) fn group_components(&mut self) -> Vec<EmComponent> {
        let mut groups = Components::new(self.singleton_alpha.len());
        for ec in &self.multi_ecs {
            let ids = &self.multi_transcripts[ec.transcript_range];
            for &tr in &ids[1..] {
                groups.join(ids[0], tr);
            }
        }
        for &tr in &self.active_transcripts {
            let root = groups.root(tr);
            groups.parents[tr] = root;
        }
        let roots = groups.parents;
        self.active_transcripts
            .sort_unstable_by_key(|&tr| (roots[tr], tr));
        let mut order: Vec<_> = (0..self.multi_ecs.len()).collect();
        order.sort_unstable_by_key(|&id| {
            (
                roots[self.multi_transcripts[self.multi_ecs[id].transcript_range.start]],
                id,
            )
        });
        let mut components = Vec::new();
        let mut tr_start = 0;
        let mut ec_start = 0;
        while tr_start < self.active_transcripts.len() {
            let root = roots[self.active_transcripts[tr_start]];
            let mut tr_end = tr_start + 1;
            while tr_end < self.active_transcripts.len()
                && roots[self.active_transcripts[tr_end]] == root
            {
                tr_end += 1;
            }
            let mut ec_end = ec_start;
            while ec_end < order.len()
                && roots
                    [self.multi_transcripts[self.multi_ecs[order[ec_end]].transcript_range.start]]
                    == root
            {
                ec_end += 1;
            }
            components.push(EmComponent {
                fixed: false,
                parallel: false,
                last_changed: true,
                scalar_finite: false,
                fixed_subtrees: 0,
                transcripts: core::range::Range {
                    start: tr_start,
                    end: tr_end,
                },
                ecs: core::range::Range {
                    start: ec_start,
                    end: ec_end,
                },
            });
            tr_start = tr_end;
            ec_start = ec_end;
        }
        debug_assert_eq!(ec_start, order.len());
        let mut ecs = Vec::with_capacity(self.multi_ecs.len());
        let mut transcripts = Vec::with_capacity(self.multi_transcripts.len());
        let mut weights = Vec::with_capacity(self.multi_weights.len());
        for id in order {
            let ec = &self.multi_ecs[id];
            let start = transcripts.len();
            transcripts.extend_from_slice(&self.multi_transcripts[ec.transcript_range]);
            weights.extend_from_slice(&self.multi_weights[ec.transcript_range]);
            ecs.push(EmMultiEc {
                transcript_range: core::range::Range {
                    start,
                    end: transcripts.len(),
                },
                count: ec.count,
            });
        }
        self.multi_ecs = ecs;
        self.multi_transcripts = transcripts;
        self.multi_weights = weights;
        components
    }
}

/// Small components retain local iterations; a component larger than a worker's
/// fair share exposes its independent EC rows and transcript sums to the pool.
/// The absolute floor avoids barriers for tiny graphs, regardless of thread count.
pub(super) fn configure_parallel(groups: &mut [EmComponent], work: &EmWork, threads: usize) {
    let threshold = work
        .multi_transcripts
        .len()
        .div_ceil(threads.max(1))
        .max(262_144);
    for group in groups {
        let ecs = &work.multi_ecs[group.ecs];
        let edges = ecs.first().zip(ecs.last()).map_or(0, |(first, last)| {
            last.transcript_range.end - first.transcript_range.start
        });
        group.parallel = threads > 1 && edges > threshold;
    }
}

// Thresholding changes alpha after the global stop decision. Every component
// must execute that final update even if it was previously an exact fixed point.
pub(super) fn reset_fixed(groups: &mut [EmComponent]) {
    for group in groups {
        group.fixed = false;
        group.last_changed = true;
        group.scalar_finite = false;
        group.fixed_subtrees = 0;
    }
}
