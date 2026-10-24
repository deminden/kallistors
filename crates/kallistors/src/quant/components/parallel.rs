//! Ordered parallel component updates and incoming-edge sums.

use super::{EmComponent, EmWork};
use rayon::prelude::*;

#[allow(clippy::too_many_arguments)]
pub(in crate::quant) fn step<I: crate::quant::IncomingEcId + Sync, const TRACK: bool>(
    pool: &rayon::ThreadPool,
    groups: &mut [EmComponent],
    work: &EmWork,
    incoming: &crate::quant::EmIncoming,
    ids: &[I],
    alpha: &mut [f64],
    scales: &mut crate::quant::EcScales,
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
fn step_range<I: crate::quant::IncomingEcId + Sync, const TRACK: bool>(
    groups: &mut [EmComponent],
    work: &EmWork,
    incoming: &crate::quant::EmIncoming,
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
                if denominator < crate::quant::TOLERANCE {
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
                    crate::quant::IncomingEcIds::Compact(ids) => sum_incoming(
                        &ids[range],
                        weights,
                        values,
                        present,
                        group.ecs.start,
                        work.singleton_alpha[tr],
                        old,
                    ),
                    crate::quant::IncomingEcIds::Wide(ids) => sum_incoming(
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
                        next > crate::quant::ALPHA_CHANGE_LIMIT
                            && (next - old).abs() / next > crate::quant::ALPHA_CHANGE,
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
fn step_parallel<I: crate::quant::IncomingEcId + Sync, const TRACK: bool>(
    group: &mut EmComponent,
    work: &EmWork,
    incoming: &crate::quant::EmIncoming,
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
                if denominator < crate::quant::TOLERANCE {
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
                    crate::quant::IncomingEcIds::Compact(ids) => sum_incoming(
                        &ids[range],
                        weights,
                        values,
                        present,
                        ec_base,
                        work.singleton_alpha[tr],
                        old,
                    ),
                    crate::quant::IncomingEcIds::Wide(ids) => sum_incoming(
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
                        && next > crate::quant::ALPHA_CHANGE_LIMIT
                        && (next - old).abs() / next > crate::quant::ALPHA_CHANGE,
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
fn sum_incoming<I: crate::quant::IncomingEcId>(
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
