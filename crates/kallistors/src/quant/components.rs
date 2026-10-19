//! Stable layout of independent EM transcript/EC components.

use super::{EmMultiEc, EmWork};
use rayon::prelude::*;

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
    // Each stable recursion interval stores its bit on its first component.
    // Ancestors share that component but use a different depth bit.
    fixed_subtrees: u64,
    transcripts: core::range::Range<usize>,
    ecs: core::range::Range<usize>,
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

#[allow(clippy::too_many_arguments)]
pub(super) fn step<I: super::IncomingEcId + Sync, const TRACK: bool>(
    pool: &rayon::ThreadPool,
    groups: &mut [EmComponent],
    work: &EmWork,
    incoming: &super::EmIncoming,
    ids: &[I],
    alpha: &mut [f64],
    scales: &mut super::EcScales,
    rounds: usize,
) -> usize {
    debug_assert!(!TRACK || rounds == 1);
    pool.install(|| {
        step_range::<I, TRACK>(
            groups,
            work,
            incoming,
            ids,
            alpha,
            &mut scales.values,
            &mut scales.present,
            0,
            0,
            rounds,
            0,
        )
    })
}

#[allow(clippy::too_many_arguments)]
fn step_range<I: super::IncomingEcId + Sync, const TRACK: bool>(
    groups: &mut [EmComponent],
    work: &EmWork,
    incoming: &super::EmIncoming,
    ids: &[I],
    alpha: &mut [f64],
    values: &mut [f64],
    present: &mut [u8],
    tr_base: usize,
    ec_base: usize,
    rounds: usize,
    depth: u32,
) -> usize {
    let Some(first) = groups.first() else {
        return 0;
    };
    let subtree = 1u64 << depth;
    if first.fixed_subtrees & subtree != 0 {
        return 0;
    }
    if groups.len() > 32 {
        let middle = groups.len() / 2;
        let tr_split = groups[middle].transcripts.start - tr_base;
        let ec_split = groups[middle].ecs.start - ec_base;
        let (left_alpha, right_alpha) = alpha.split_at_mut(tr_split);
        let (left_values, right_values) = values.split_at_mut(ec_split);
        let (left_present, right_present) = present.split_at_mut(ec_split);
        let (left_groups, right_groups) = groups.split_at_mut(middle);
        let (left, right) = rayon::join(
            || {
                step_range::<I, TRACK>(
                    left_groups,
                    work,
                    incoming,
                    ids,
                    left_alpha,
                    left_values,
                    left_present,
                    tr_base,
                    ec_base,
                    rounds,
                    depth + 1,
                )
            },
            || {
                step_range::<I, TRACK>(
                    right_groups,
                    work,
                    incoming,
                    ids,
                    right_alpha,
                    right_values,
                    right_present,
                    tr_base + tr_split,
                    ec_base + ec_split,
                    rounds,
                    depth + 1,
                )
            },
        );
        let children = 1u64 << (depth + 1);
        if groups[0].fixed_subtrees & children != 0 && groups[middle].fixed_subtrees & children != 0
        {
            groups[0].fixed_subtrees |= subtree;
        }
        return left + right;
    }
    let mut changes = 0;
    let mut all_fixed = true;
    for group in groups.iter_mut() {
        if group.fixed {
            continue;
        }
        let alpha = &mut alpha[group.transcripts.start - tr_base..group.transcripts.end - tr_base];
        let values = &mut values[group.ecs.start - ec_base..group.ecs.end - ec_base];
        let present = &mut present[group.ecs.start - ec_base..group.ecs.end - ec_base];
        if group.parallel {
            changes += step_parallel::<I, TRACK>(
                group, work, incoming, ids, alpha, values, present, rounds,
            );
            all_fixed &= group.fixed;
            continue;
        }
        for _ in 0..rounds {
            let mut unchanged = true;
            for (ec, (scale, present)) in work.multi_ecs[group.ecs]
                .iter()
                .zip(values.iter_mut().zip(present.iter_mut()))
            {
                let mut denominator = 0.0;
                for (&tr, &weight) in ids[ec.transcript_range]
                    .iter()
                    .zip(&work.multi_weights[ec.transcript_range])
                {
                    denominator += alpha[tr.as_usize() - group.transcripts.start] * weight;
                }
                if denominator < super::TOLERANCE {
                    *present = 0;
                    *scale = 0.0;
                } else {
                    *present = 1;
                    *scale = ec.count / denominator;
                }
            }
            for (local, value) in alpha.iter_mut().enumerate() {
                let tr = group.transcripts.start + local;
                let old = *value;
                let range = incoming.offsets[tr]..incoming.offsets[tr + 1];
                let weights = &incoming.weights[range.clone()];
                let next = match &incoming.ec_indices {
                    super::IncomingEcIds::Compact(ids) => sum_incoming(
                        &ids[range],
                        weights,
                        values,
                        present,
                        group.ecs.start,
                        work.singleton_alpha[tr],
                        old,
                    ),
                    super::IncomingEcIds::Wide(ids) => sum_incoming(
                        &ids[range],
                        weights,
                        values,
                        present,
                        group.ecs.start,
                        work.singleton_alpha[tr],
                        old,
                    ),
                };
                *value = next;
                unchanged &= next.to_bits() == old.to_bits();
                if TRACK {
                    changes += usize::from(
                        next > super::ALPHA_CHANGE_LIMIT
                            && (next - old).abs() / next > super::ALPHA_CHANGE,
                    );
                }
            }
            if unchanged {
                group.fixed = true;
                break;
            }
        }
        all_fixed &= group.fixed;
    }
    if all_fixed {
        groups[0].fixed_subtrees |= subtree;
    }
    changes
}

#[allow(clippy::too_many_arguments)]
fn step_parallel<I: super::IncomingEcId + Sync, const TRACK: bool>(
    group: &mut EmComponent,
    work: &EmWork,
    incoming: &super::EmIncoming,
    ids: &[I],
    alpha: &mut [f64],
    values: &mut [f64],
    present: &mut [u8],
    rounds: usize,
) -> usize {
    let tr_base = group.transcripts.start;
    let ec_base = group.ecs.start;
    let mut changes = 0;
    for _ in 0..rounds {
        // Complete every denominator before any alpha update. Each scalar sum
        // retains the serial edge order; only independent results run together.
        values
            .par_iter_mut()
            .zip(present.par_iter_mut())
            .zip(&work.multi_ecs[group.ecs])
            .for_each(|((scale, present), ec)| {
                let mut denominator = 0.0;
                for (&tr, &weight) in ids[ec.transcript_range]
                    .iter()
                    .zip(&work.multi_weights[ec.transcript_range])
                {
                    denominator += alpha[tr.as_usize() - tr_base] * weight;
                }
                if denominator < super::TOLERANCE {
                    *present = 0;
                    *scale = 0.0;
                } else {
                    *present = 1;
                    *scale = ec.count / denominator;
                }
            });
        let (changed, unchanged) = alpha
            .par_iter_mut()
            .enumerate()
            .map(|(local, value)| {
                let tr = tr_base + local;
                let old = *value;
                let range = incoming.offsets[tr]..incoming.offsets[tr + 1];
                let weights = &incoming.weights[range.clone()];
                let next = match &incoming.ec_indices {
                    super::IncomingEcIds::Compact(ids) => sum_incoming(
                        &ids[range],
                        weights,
                        values,
                        present,
                        ec_base,
                        work.singleton_alpha[tr],
                        old,
                    ),
                    super::IncomingEcIds::Wide(ids) => sum_incoming(
                        &ids[range],
                        weights,
                        values,
                        present,
                        ec_base,
                        work.singleton_alpha[tr],
                        old,
                    ),
                };
                *value = next;
                let changed = usize::from(
                    TRACK
                        && next > super::ALPHA_CHANGE_LIMIT
                        && (next - old).abs() / next > super::ALPHA_CHANGE,
                );
                (changed, next.to_bits() == old.to_bits())
            })
            .reduce(
                || (0, true),
                |(a, fixed_a), (b, fixed_b)| (a + b, fixed_a && fixed_b),
            );
        changes += changed;
        if unchanged {
            group.fixed = true;
            break;
        }
    }
    changes
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn sum_incoming<I: super::IncomingEcId>(
    ids: &[I],
    weights: &[f64],
    values: &[f64],
    present: &[u8],
    base: usize,
    singleton: f64,
    alpha: f64,
) -> f64 {
    let mut sum = singleton;
    for (&ec, &weight) in ids.iter().zip(weights) {
        let ec = ec.as_usize() - base;
        if present[ec] != 0 {
            sum += weight * alpha * values[ec];
        }
    }
    sum
}

// Thresholding changes alpha after the global stop decision. Every component
// must execute that final update even if it was previously an exact fixed point.
pub(super) fn reset_fixed(groups: &mut [EmComponent]) {
    for group in groups {
        group.fixed = false;
        group.fixed_subtrees = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_fixed_components_resume_after_alpha_reset() {
        for (parallel, wide) in [(false, false), (true, false), (true, true)] {
            let n = 70;
            let classes: Vec<_> = (0..n)
                .map(|i| vec![(i * 2) as u32, (i * 2 + 1) as u32])
                .collect();
            let lengths: Vec<_> = (0..n).flat_map(|_| [1.0, 2.0]).collect();
            let mut work = super::super::prepare_em_work(&classes, &vec![1; n], &lengths);
            let mut groups = work.group_components();
            work.compact_transcript_ids();
            for group in &mut groups {
                group.parallel = parallel;
            }
            let mut incoming = super::super::EmIncoming::new(&work, false);
            if wide && let super::super::IncomingEcIds::Compact(ids) = &incoming.ec_indices {
                incoming.ec_indices =
                    super::super::IncomingEcIds::Wide(ids.iter().map(|&id| id as usize).collect());
            }
            let mut scales = super::super::EcScales::new(n);
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(2)
                .build()
                .unwrap();
            // Finish one complete branch while the other keeps changing. The
            // cached branch must skip work without suppressing live neighbors.
            let mut mixed: Vec<f64> = (0..n)
                .flat_map(|i| if i < n / 2 { [1.0, 0.0] } else { [0.5, 0.5] })
                .collect();
            for _ in 0..3 {
                let mut expected = mixed.clone();
                for pair in expected.as_chunks_mut::<2>().0 {
                    let denominator = 0.0 + pair[0] * 1.0 + pair[1] * 0.5;
                    let scale = 1.0 / denominator;
                    let first = 0.0 + 1.0 * pair[0] * scale;
                    let second = 0.0 + 0.5 * pair[1] * scale;
                    pair.copy_from_slice(&[first, second]);
                }
                step::<_, true>(
                    &pool,
                    &mut groups,
                    &work,
                    &incoming,
                    &work.multi_transcripts,
                    &mut mixed,
                    &mut scales,
                    1,
                );
                assert_eq!(
                    mixed.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert_ne!(groups[0].fixed_subtrees, 0);
                assert_eq!(groups[0].fixed_subtrees & 1, 0);
            }
            reset_fixed(&mut groups);
            let mut alpha: Vec<_> = (0..n).flat_map(|_| [1.0, 0.0]).collect();
            step::<_, true>(
                &pool,
                &mut groups,
                &work,
                &incoming,
                &work.multi_transcripts,
                &mut alpha,
                &mut scales,
                1,
            );
            assert!(groups.iter().all(|g| g.fixed));
            assert_ne!(groups[0].fixed_subtrees & 1, 0);
            alpha.fill(0.5);
            reset_fixed(&mut groups);
            assert!(groups.iter().all(|g| !g.fixed && g.fixed_subtrees == 0));
            step::<_, true>(
                &pool,
                &mut groups,
                &work,
                &incoming,
                &work.multi_transcripts,
                &mut alpha,
                &mut scales,
                1,
            );
            assert_eq!(
                alpha,
                (0..n)
                    .flat_map(|_| [2.0 / 3.0, 1.0 / 3.0])
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn balanced_union_trees_preserve_every_transcripts_ec_order() {
        let mut classes = Vec::new();
        for width in [1, 2, 4, 8, 16] {
            for first in (0..32).step_by(width * 2) {
                for component in 0..2 {
                    classes.push(vec![
                        (first * 2 + component) as u32,
                        ((first + width) * 2 + component) as u32,
                    ]);
                }
            }
        }
        let counts: Vec<_> = (1..=classes.len() as u32).collect();
        let mut work = super::super::prepare_em_work(&classes, &counts, &[1.0; 64]);
        let visits = |work: &EmWork| {
            let mut visits = vec![Vec::new(); 64];
            for ec in &work.multi_ecs {
                for &tr in &work.multi_transcripts[ec.transcript_range] {
                    visits[tr].push(ec.count.to_bits());
                }
            }
            visits
        };
        let before = visits(&work);
        let _ = work.group_components();
        assert_eq!(visits(&work), before);
        assert!(work.active_transcripts[..32].iter().all(|tr| tr % 2 == 0));
        assert!(work.active_transcripts[32..].iter().all(|tr| tr % 2 == 1));
    }
}
