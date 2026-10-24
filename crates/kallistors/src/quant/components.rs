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
    last_changed: bool,
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
                last_changed: true,
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
        group.last_changed = true;
        group.fixed_subtrees = 0;
    }
}

/// Skip only components whose entire alpha vector is bitwise unchanged.
/// Their weights and counts are immutable and no other component contributes
/// to their transcripts. Global final pruning resets these certificates.
pub(super) fn scalar_step(
    groups: &mut [EmComponent],
    work: &EmWork,
    active: &[usize],
    alpha: &mut [f64],
    next_alpha: &mut [f64],
    products: &mut [f64],
) -> usize {
    // Scalar component callers group and compact before entering this path.
    debug_assert!(active.iter().copied().eq(0..alpha.len()));
    let mut changes = 0;
    for group in groups {
        if group.fixed {
            continue;
        }
        let (count, identical) = work.scalar_dense_step_mode::<2>(
            group.ecs,
            group.transcripts,
            alpha,
            next_alpha,
            products,
        );
        changes += count;
        group.last_changed = count != 0;
        group.fixed = identical;
    }
    changes
}

// A nonzero component count proves the global count is nonzero. Advance one
// witness to its first eligible zero (or the bounded cap), then advance every
// other independent component to that same boundary. No possible global stop
// is crossed, so neither an alpha snapshot nor speculative replay is needed.
#[allow(clippy::too_many_arguments)]
pub(super) fn scalar_tile(
    groups: &mut [EmComponent],
    work: &EmWork,
    active: &[usize],
    alpha: &mut [f64],
    next_alpha: &mut [f64],
    products: &mut [f64],
    first_iter: usize,
    min_rounds: usize,
    remaining: usize,
) -> (usize, usize) {
    debug_assert!(active.iter().copied().eq(0..alpha.len()));
    let cap = remaining.min(64);
    debug_assert!(cap > 0);
    if cap == 1 {
        return (
            1,
            scalar_step(groups, work, active, alpha, next_alpha, products),
        );
    }
    let witness = groups
        .iter()
        .enumerate()
        .filter(|(_, group)| !group.fixed)
        .max_by_key(|(_, group)| (group.last_changed, group.ecs.end - group.ecs.start))
        .map(|(index, _)| index);
    let Some(witness) = witness else {
        let taken = (0..cap)
            .position(|offset| first_iter + offset > min_rounds)
            .map_or(cap, |offset| offset + 1);
        return (taken, 0);
    };
    let group = &mut groups[witness];
    let mut taken = 0;
    let mut changes = 0;
    for offset in 0..cap {
        let eligible = first_iter + offset > min_rounds;
        let (count, identical) = if group.fixed {
            (0, true)
        } else if offset + 1 == cap {
            work.scalar_dense_step_mode::<2>(
                group.ecs,
                group.transcripts,
                alpha,
                next_alpha,
                products,
            )
        } else if eligible {
            work.scalar_dense_step_mode::<1>(
                group.ecs,
                group.transcripts,
                alpha,
                next_alpha,
                products,
            )
        } else {
            work.scalar_dense_step_mode::<0>(
                group.ecs,
                group.transcripts,
                alpha,
                next_alpha,
                products,
            )
        };
        group.fixed = identical;
        taken = offset + 1;
        changes = count;
        if eligible && count == 0 {
            break;
        }
    }
    // At a cap the full count was computed; at an earlier exit any-change
    // returned zero, which is also the exact full count. Both are observable.
    group.last_changed = changes != 0;
    for (index, group) in groups.iter_mut().enumerate() {
        if index == witness || group.fixed {
            continue;
        }
        let mut last_count = 0;
        for offset in 0..taken {
            let (count, identical) = if offset + 1 == taken {
                work.scalar_dense_step_mode::<2>(
                    group.ecs,
                    group.transcripts,
                    alpha,
                    next_alpha,
                    products,
                )
            } else {
                work.scalar_dense_step_mode::<0>(
                    group.ecs,
                    group.transcripts,
                    alpha,
                    next_alpha,
                    products,
                )
            };
            last_count = count;
            group.fixed = identical;
            if identical {
                break;
            }
        }
        group.last_changed = last_count != 0;
        changes += last_count;
    }
    (taken, changes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_tiles_preserve_first_stop_and_reset() {
        // Includes independent components that reach relative convergence at
        // different times. Every returned boundary is checked against the
        // original one-round full-count equations, including all intervening
        // possible stops, rather than against the witness implementation.
        let fixtures = [
            (
                vec![vec![0, 1]],
                vec![5],
                vec![1.0, 2.0],
                vec![4.99975, 0.00025],
            ),
            (
                vec![
                    vec![0, 2],
                    vec![1, 3],
                    vec![4],
                    vec![0, 2, 5],
                    vec![1, 3, 6],
                ],
                vec![5, 11, 3, 7, 2],
                vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
                vec![0.25; 7],
            ),
            (
                vec![vec![0, 1], vec![0, 1], vec![0, 1], vec![2, 3]],
                vec![1, 2, 2, 19],
                vec![1.0, 2.0, 1.0, 2.0],
                vec![4.99975, 0.00025, 0.5, 0.5],
            ),
            (
                vec![vec![0], vec![1]],
                vec![3, 7],
                vec![1.0, 1.0],
                vec![0.5, 0.5],
            ),
            (
                vec![vec![2, 8], vec![8, 11], vec![4, 9]],
                vec![5, 7, 11],
                (0..12).map(|id| 1.0 + id as f64 * 0.25).collect(),
                vec![0.25; 12],
            ),
            (vec![], vec![], vec![], vec![]),
        ];
        let mut early_boundary_checked = false;
        let mut witness_switched = false;
        for (classes, counts, lengths, initial) in fixtures {
            for fill in [
                None,
                Some(0.0),
                Some(-0.0),
                Some(f64::from_bits(1)),
                Some(f64::INFINITY),
                Some(f64::from_bits(0x7ff8_0000_0000_0023)),
            ] {
                for (first_iter, min_rounds) in [(0, 0), (0, 2), (49, 50), (51, 50), (0, 200)] {
                    for remaining in [1, 2, 7, 8, 9, 63, 64, 65, 130] {
                        let mut work = super::super::prepare_em_work(&classes, &counts, &lengths);
                        let mut groups = work.group_components();
                        let mut actual: Vec<_> = work
                            .active_transcripts
                            .iter()
                            .map(|&tr| fill.unwrap_or(initial[tr]))
                            .collect();
                        work.compact_transcript_ids();
                        let active: Vec<_> = (0..actual.len()).collect();
                        let mut expected = actual.clone();
                        let mut next = vec![0.0; actual.len()];
                        let max_ec = work
                            .multi_ecs
                            .iter()
                            .map(|ec| ec.transcript_range.end - ec.transcript_range.start)
                            .max()
                            .unwrap_or(0);
                        let mut products = vec![0.0; max_ec.saturating_mul(2)];
                        let mut completed = 0;
                        while completed < remaining {
                            let wanted_witness = groups
                                .iter()
                                .enumerate()
                                .filter(|(_, group)| !group.fixed)
                                .max_by_key(|(_, group)| {
                                    (group.last_changed, group.ecs.end - group.ecs.start)
                                })
                                .map(|(index, _)| index);
                            let (taken, count) = scalar_tile(
                                &mut groups,
                                &work,
                                &active,
                                &mut actual,
                                &mut next,
                                &mut products,
                                first_iter + completed,
                                min_rounds,
                                remaining - completed,
                            );
                            assert!(taken > 0 && taken <= 64 && taken <= remaining - completed);
                            let mut expected_count = 0;
                            for offset in 0..taken {
                                expected_count = work
                                    .scalar_step(
                                        core::range::Range {
                                            start: 0,
                                            end: work.multi_ecs.len(),
                                        },
                                        &active,
                                        &mut expected,
                                        &mut next,
                                        &mut products,
                                    )
                                    .0;
                                if offset + 1 < taken
                                    && first_iter + completed + offset > min_rounds
                                {
                                    assert_ne!(
                                        expected_count, 0,
                                        "tile crossed the first valid global stop"
                                    );
                                }
                            }
                            assert_eq!(count, expected_count);
                            assert!(
                                actual
                                    .iter()
                                    .zip(&expected)
                                    .all(|(a, b)| a.to_bits() == b.to_bits())
                            );
                            early_boundary_checked |=
                                taken < (remaining - completed).min(64) && fill.is_none();
                            let next_witness = groups
                                .iter()
                                .enumerate()
                                .filter(|(_, group)| !group.fixed)
                                .max_by_key(|(_, group)| {
                                    (group.last_changed, group.ecs.end - group.ecs.start)
                                })
                                .map(|(index, _)| index);
                            witness_switched |= wanted_witness.is_some()
                                && next_witness.is_some()
                                && wanted_witness != next_witness;
                            completed += taken;
                            if count == 0 && first_iter + completed - 1 > min_rounds {
                                break;
                            }
                        }
                        // Independent iteration-by-iteration stop oracle.
                        let mut reference: Vec<_> = work
                            .active_transcripts
                            .iter()
                            .map(|&tr| fill.unwrap_or(initial[tr]))
                            .collect();
                        let mut stop = 0;
                        for offset in 0..remaining {
                            let count = work
                                .scalar_step(
                                    core::range::Range {
                                        start: 0,
                                        end: work.multi_ecs.len(),
                                    },
                                    &active,
                                    &mut reference,
                                    &mut next,
                                    &mut products,
                                )
                                .0;
                            stop += 1;
                            if count == 0 && first_iter + offset > min_rounds {
                                break;
                            }
                        }
                        assert_eq!(completed, stop);
                        assert!(
                            actual
                                .iter()
                                .zip(&reference)
                                .all(|(a, b)| a.to_bits() == b.to_bits())
                        );
                        // Final pruning must reset exact fixed certificates.
                        for (a, b) in actual.iter_mut().zip(&mut expected) {
                            if *a < super::super::ALPHA_LIMIT / 10.0 {
                                *a = 0.0;
                            }
                            if *b < super::super::ALPHA_LIMIT / 10.0 {
                                *b = 0.0;
                            }
                        }
                        reset_fixed(&mut groups);
                        assert!(
                            groups
                                .iter()
                                .all(|group| !group.fixed && group.last_changed)
                        );
                        let expected_count = work
                            .scalar_step(
                                core::range::Range {
                                    start: 0,
                                    end: work.multi_ecs.len(),
                                },
                                &active,
                                &mut expected,
                                &mut next,
                                &mut products,
                            )
                            .0;
                        assert_eq!(
                            scalar_step(
                                &mut groups,
                                &work,
                                &active,
                                &mut actual,
                                &mut next,
                                &mut products
                            ),
                            expected_count
                        );
                        assert!(
                            actual
                                .iter()
                                .zip(&expected)
                                .all(|(a, b)| a.to_bits() == b.to_bits())
                        );
                    }
                }
            }
        }
        assert!(early_boundary_checked);
        assert!(witness_switched);
        // The additional scheduling bool fits the existing layout padding.
        assert_eq!(std::mem::size_of::<EmComponent>(), 48);
    }

    #[test]
    fn exact_fixed_components_resume_after_alpha_reset() {
        for (scalar, parallel, wide) in [
            (true, false, false),
            (false, false, false),
            (false, true, false),
            (false, true, true),
        ] {
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
            let active: Vec<_> = (0..work.active_transcripts.len()).collect();
            let mut next = vec![0.0; active.len()];
            let mut products = vec![0.0; 2];
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
                if scalar {
                    scalar_step(
                        &mut groups,
                        &work,
                        &active,
                        &mut mixed,
                        &mut next,
                        &mut products,
                    );
                } else {
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
                }
                assert_eq!(
                    mixed.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                if scalar {
                    assert!(groups[..n / 2].iter().all(|g| g.fixed));
                } else {
                    assert_ne!(groups[0].fixed_subtrees, 0);
                    assert_eq!(groups[0].fixed_subtrees & 1, 0);
                }
            }
            reset_fixed(&mut groups);
            let mut alpha: Vec<_> = (0..n).flat_map(|_| [1.0, 0.0]).collect();
            if scalar {
                scalar_step(
                    &mut groups,
                    &work,
                    &active,
                    &mut alpha,
                    &mut next,
                    &mut products,
                );
            } else {
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
            }
            assert!(groups.iter().all(|g| g.fixed));
            if !scalar {
                assert_ne!(groups[0].fixed_subtrees & 1, 0);
            }
            alpha.fill(0.5);
            reset_fixed(&mut groups);
            assert!(groups.iter().all(|g| !g.fixed && g.fixed_subtrees == 0));
            if scalar {
                scalar_step(
                    &mut groups,
                    &work,
                    &active,
                    &mut alpha,
                    &mut next,
                    &mut products,
                );
            } else {
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
            }
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

        // Unequal row lengths, overlapping targets, odd final rows and
        // near-zero denominators check the scalar two-row schedule against
        // the original row-at-a-time arithmetic, bit for bit.
        for rows in 1..=5 {
            let classes = [
                vec![0, 1, 2, 3],
                vec![1, 3],
                vec![0, 2, 3],
                vec![0, 1],
                vec![1, 2, 3],
            ];
            let counts = [17, 3, 101, 2, 7];
            let work = super::super::prepare_em_work(
                &classes[..rows],
                &counts[..rows],
                &[1.125, 2.75, 11.25, 17.0],
            );
            let active = work.active_transcripts.clone();
            for initial in [
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [1e-16, 1e-16, 1e-16, 1e-16],
                [0.125, 0.25, 0.375, 0.25],
                [-0.0; 4],
                [f64::from_bits(1); 4],
                [f64::INFINITY; 4],
                [f64::from_bits(0x7ff8_0000_0000_0023); 4],
            ] {
                let mut expected = initial.to_vec();
                let mut actual = initial.to_vec();
                let mut next = vec![0.0; 4];
                let mut products = vec![0.0; 8];
                for iteration in 0..20 {
                    if iteration == 10 {
                        expected[0] = 0.0;
                        actual[0] = 0.0;
                    }
                    let mut wanted = work.singleton_alpha.clone();
                    for ec in &work.multi_ecs {
                        let ids = &work.multi_transcripts[ec.transcript_range];
                        let weights = &work.multi_weights[ec.transcript_range];
                        let mut denominator = 0.0;
                        let products: Vec<_> = ids
                            .iter()
                            .zip(weights)
                            .map(|(&id, &weight)| {
                                let product = expected[id] * weight;
                                denominator += product;
                                product
                            })
                            .collect();
                        if denominator < super::super::TOLERANCE {
                            continue;
                        }
                        let scale = ec.count / denominator;
                        for (&id, product) in ids.iter().zip(products) {
                            wanted[id] += product * scale;
                        }
                    }
                    for &id in &active {
                        expected[id] = wanted[id];
                    }
                    work.scalar_step(
                        core::range::Range {
                            start: 0,
                            end: rows,
                        },
                        &active,
                        &mut actual,
                        &mut next,
                        &mut products,
                    );
                    assert_eq!(
                        actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                    );
                }
            }
        }
    }
}
