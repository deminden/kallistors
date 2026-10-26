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
            None,
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
    skip_tr: Option<usize>,
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
                    skip_tr,
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
                    skip_tr,
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
        if skip_tr == Some(group.transcripts.start) {
            all_fixed = false;
            continue;
        }
        let trs = group.transcripts;
        let ecs = group.ecs;
        changes += step_group::<I, TRACK>(
            group,
            work,
            incoming,
            ids,
            &mut alpha[trs.start - tr_base..trs.end - tr_base],
            &mut values[ecs.start - ec_base..ecs.end - ec_base],
            &mut present[ecs.start - ec_base..ecs.end - ec_base],
            rounds,
        );
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
        changes = 0;
        // Complete every denominator before any alpha update. Each scalar sum
        // retains the serial edge order; only independent results run together.
        let finite_scales = values
            .par_chunks_mut(2)
            .zip(present.par_chunks_mut(2))
            .zip(work.multi_ecs[group.ecs].par_chunks(2))
            .map(|((scale_rows, present_rows), rows)| {
                let denominators =
                    paired_denominators(rows, ids, &work.multi_weights, alpha, tr_base);
                for (((scale, present), ec), denominator) in scale_rows
                    .iter_mut()
                    .zip(present_rows)
                    .zip(rows)
                    .zip(denominators)
                {
                    if denominator < crate::quant::TOLERANCE {
                        *present = 0;
                        *scale = 0.0;
                    } else {
                        *present = 1;
                        *scale = ec.count / denominator;
                    }
                }
                scale_rows.iter().all(|value| value.is_finite())
            })
            // This reduction still evaluates all denominator chunks.
            .reduce(|| true, |left, right| left & right);
        let zero_safe = incoming.finite_weights && finite_scales;
        let (changed, unchanged) = alpha
            .par_chunks_mut(2)
            .enumerate()
            .map(|(local, pair)| {
                let tr = tr_base + local * 2;
                let next = match &incoming.ec_indices {
                    crate::quant::IncomingEcIds::Compact(ids) => adjacent_incoming(
                        incoming,
                        ids,
                        work,
                        tr,
                        pair,
                        IncomingScales {
                            values,
                            present,
                            base: ec_base,
                            zero_safe,
                        },
                    ),
                    crate::quant::IncomingEcIds::Wide(ids) => adjacent_incoming(
                        incoming,
                        ids,
                        work,
                        tr,
                        pair,
                        IncomingScales {
                            values,
                            present,
                            base: ec_base,
                            zero_safe,
                        },
                    ),
                };
                let mut changes = 0;
                let mut identical = true;
                for (value, next) in pair.iter_mut().zip(next) {
                    let old = *value;
                    *value = next;
                    changes += usize::from(
                        TRACK
                            && next > crate::quant::ALPHA_CHANGE_LIMIT
                            && (next - old).abs() / next > crate::quant::ALPHA_CHANGE,
                    );
                    identical &= next.to_bits() == old.to_bits();
                }
                (changes, identical)
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

#[allow(clippy::too_many_arguments)]
fn step_group<I: crate::quant::IncomingEcId + Sync, const TRACK: bool>(
    group: &mut EmComponent,
    work: &EmWork,
    incoming: &crate::quant::EmIncoming,
    ids: &[I],
    alpha: &mut [f64],
    values: &mut [f64],
    present: &mut [u8],
    rounds: usize,
) -> usize {
    // Select static arithmetic kernels once for the complete catch-up tile.
    // No intermediate count escapes this component. Fixed-point exits have
    // zero changes; otherwise only the final tracked update is observable.
    if TRACK && rounds > 1 {
        step_group::<I, false>(
            group,
            work,
            incoming,
            ids,
            alpha,
            values,
            present,
            rounds - 1,
        );
        if group.fixed {
            group.last_changed = false;
            return 0;
        }
        return step_group::<I, true>(group, work, incoming, ids, alpha, values, present, 1);
    }
    if group.parallel {
        let changes =
            step_parallel::<I, TRACK>(group, work, incoming, ids, alpha, values, present, rounds);
        if TRACK {
            group.last_changed = changes != 0;
        }
        return changes;
    }
    let mut changes = 0;
    for _ in 0..rounds {
        changes = 0;
        let mut unchanged = true;
        let mut finite_scales = true;
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
            finite_scales &= scale.is_finite();
        }
        let zero_safe = incoming.finite_weights && finite_scales;
        for (local, value) in alpha.iter_mut().enumerate() {
            let tr = group.transcripts.start + local;
            let old = *value;
            if zero_safe && old.to_bits() == 0 && work.singleton_alpha[tr].to_bits() == 0 {
                continue;
            }
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
    if TRACK {
        group.last_changed = changes != 0;
    }
    changes
}

/// An independent nonzero witness excludes every earlier global stop. Other
/// components catch up only to that proven boundary, using the same kernels.
#[allow(clippy::too_many_arguments)]
pub(in crate::quant) fn tile<I: crate::quant::IncomingEcId + Sync>(
    pool: &rayon::ThreadPool,
    groups: &mut [EmComponent],
    work: &EmWork,
    incoming: &crate::quant::EmIncoming,
    ids: &[I],
    alpha: &mut [f64],
    scales: &mut crate::quant::EcScales,
    first_iter: usize,
    minimum: usize,
    remaining: usize,
) -> (usize, usize) {
    let cap = remaining.min(64);
    debug_assert!(cap > 0);
    let Some(witness) = super::witness(groups) else {
        let taken = (0..cap)
            .position(|offset| first_iter + offset > minimum)
            .map_or(cap, |offset| offset + 1);
        return (taken, 0);
    };
    pool.install(|| {
        let group = &mut groups[witness];
        let start = group.transcripts.start;
        let mut taken = 0;
        let mut changes = 0;
        for offset in 0..cap {
            let trs = group.transcripts;
            let ecs = group.ecs;
            changes = if group.fixed {
                0
            } else {
                step_group::<I, true>(
                    group,
                    work,
                    incoming,
                    ids,
                    &mut alpha[trs.start..trs.end],
                    &mut scales.values[ecs.start..ecs.end],
                    &mut scales.present[ecs.start..ecs.end],
                    1,
                )
            };
            taken = offset + 1;
            if first_iter + offset > minimum && changes == 0 {
                break;
            }
        }
        // Preserve original tree/depth ownership. Computing the witness through
        // a singleton step_range would attach subtree bits to the wrong tree.
        let others = step_range::<I, true>(
            groups,
            work,
            incoming,
            ids,
            alpha,
            &mut scales.values,
            &mut scales.present,
            0,
            0,
            taken,
            0,
            Some(start),
        );
        (taken, changes + others)
    })
}

#[inline]
fn paired_denominators<I: crate::quant::IncomingEcId>(
    rows: &[crate::quant::EmMultiEc],
    ids: &[I],
    weights: &[f64],
    alpha: &[f64],
    base: usize,
) -> [f64; 2] {
    debug_assert!(!rows.is_empty() && rows.len() <= 2);
    let left = &rows[0];
    let left_ids = &ids[left.transcript_range];
    let left_weights = &weights[left.transcript_range];
    let mut sums = [0.0; 2];
    if let Some(right) = rows.get(1) {
        let right_ids = &ids[right.transcript_range];
        let right_weights = &weights[right.transcript_range];
        let common = left_ids.len().min(right_ids.len());
        for ((&left, &lw), (&right, &rw)) in left_ids[..common]
            .iter()
            .zip(&left_weights[..common])
            .zip(right_ids[..common].iter().zip(&right_weights[..common]))
        {
            sums[0] += alpha[left.as_usize() - base] * lw;
            sums[1] += alpha[right.as_usize() - base] * rw;
        }
        for (&id, &weight) in left_ids[common..].iter().zip(&left_weights[common..]) {
            sums[0] += alpha[id.as_usize() - base] * weight;
        }
        for (&id, &weight) in right_ids[common..].iter().zip(&right_weights[common..]) {
            sums[1] += alpha[id.as_usize() - base] * weight;
        }
    } else {
        for (&id, &weight) in left_ids.iter().zip(left_weights) {
            sums[0] += alpha[id.as_usize() - base] * weight;
        }
    }
    sums
}

#[derive(Clone, Copy)]
struct IncomingScales<'a> {
    values: &'a [f64],
    present: &'a [u8],
    base: usize,
    zero_safe: bool,
}

struct IncomingRow<'a, I> {
    ids: &'a [I],
    weights: &'a [f64],
    singleton: f64,
    alpha: f64,
}

#[inline]
fn paired_incoming<I: crate::quant::IncomingEcId>(
    left: IncomingRow<'_, I>,
    right: IncomingRow<'_, I>,
    scales: IncomingScales<'_>,
) -> [f64; 2] {
    let IncomingScales {
        values,
        present,
        base,
        ..
    } = scales;
    let common = left.ids.len().min(right.ids.len());
    let mut sums = [left.singleton, right.singleton];
    for ((&left_id, &lw), (&right_id, &rw)) in left.ids[..common]
        .iter()
        .zip(&left.weights[..common])
        .zip(right.ids[..common].iter().zip(&right.weights[..common]))
    {
        let le = left_id.as_usize() - base;
        if present[le] != 0 {
            sums[0] += lw * left.alpha * values[le];
        }
        let re = right_id.as_usize() - base;
        if present[re] != 0 {
            sums[1] += rw * right.alpha * values[re];
        }
    }
    sums[0] = sum_incoming(
        &left.ids[common..],
        &left.weights[common..],
        values,
        present,
        base,
        sums[0],
        left.alpha,
    );
    sums[1] = sum_incoming(
        &right.ids[common..],
        &right.weights[common..],
        values,
        present,
        base,
        sums[1],
        right.alpha,
    );
    sums
}

#[inline]
fn adjacent_incoming<I: crate::quant::IncomingEcId>(
    incoming: &crate::quant::EmIncoming,
    ids: &[I],
    work: &EmWork,
    tr: usize,
    old: &[f64],
    scales: IncomingScales<'_>,
) -> [f64; 2] {
    let IncomingScales {
        values,
        present,
        base,
        ..
    } = scales;
    let row = |slot: usize, alpha: f64| {
        let singleton = work.singleton_alpha[slot];
        if scales.zero_safe && alpha.to_bits() == 0 && singleton.to_bits() == 0 {
            IncomingRow {
                ids: &ids[..0],
                weights: &incoming.weights[..0],
                singleton,
                alpha,
            }
        } else {
            let range = incoming.offsets[slot]..incoming.offsets[slot + 1];
            IncomingRow {
                ids: &ids[range.clone()],
                weights: &incoming.weights[range],
                singleton,
                alpha,
            }
        }
    };
    let left = row(tr, old[0]);
    if old.len() == 1 {
        return [
            sum_incoming(
                left.ids,
                left.weights,
                values,
                present,
                base,
                left.singleton,
                left.alpha,
            ),
            0.0,
        ];
    }
    paired_incoming(left, row(tr + 1, old[1]), scales)
}

#[cfg(test)]
mod incoming_tests {
    use super::*;
    #[test]
    fn incoming_pairs_preserve_order_missing_rows_and_exceptional_bits() {
        let ids = [2usize, 3, 4, 5, 6];
        let compact = ids.map(|x| x as u32);
        for left_len in 0..=5 {
            for right_len in 0..=5 {
                for alpha in [[0.0, -0.0], [f64::MAX, 1.0], [f64::NAN, f64::INFINITY]] {
                    for present in [[1, 1, 1, 1, 1], [0, 1, 0, 1, 0], [0, 0, 0, 0, 0]] {
                        let weights = [-0.0, f64::from_bits(1), 2.0, f64::INFINITY, f64::NAN];
                        let values = [1.0, 2.0, f64::MAX, f64::INFINITY, 0.0];
                        let mut expected = [0.0, 3.0];
                        for (side, len) in [left_len, right_len].into_iter().enumerate() {
                            for edge in 0..len {
                                if present[edge] != 0 {
                                    expected[side] += weights[edge] * alpha[side] * values[edge];
                                }
                            }
                        }
                        let wide = paired_incoming(
                            IncomingRow {
                                ids: &ids[..left_len],
                                weights: &weights[..left_len],
                                singleton: 0.0,
                                alpha: alpha[0],
                            },
                            IncomingRow {
                                ids: &ids[..right_len],
                                weights: &weights[..right_len],
                                singleton: 3.0,
                                alpha: alpha[1],
                            },
                            IncomingScales {
                                values: &values,
                                present: &present,
                                base: 2,
                                zero_safe: false,
                            },
                        );
                        let narrow = paired_incoming(
                            IncomingRow {
                                ids: &compact[..left_len],
                                weights: &weights[..left_len],
                                singleton: 0.0,
                                alpha: alpha[0],
                            },
                            IncomingRow {
                                ids: &compact[..right_len],
                                weights: &weights[..right_len],
                                singleton: 3.0,
                                alpha: alpha[1],
                            },
                            IncomingScales {
                                values: &values,
                                present: &present,
                                base: 2,
                                zero_safe: false,
                            },
                        );
                        assert_eq!(wide.map(f64::to_bits), expected.map(f64::to_bits));
                        assert_eq!(narrow.map(f64::to_bits), expected.map(f64::to_bits));
                    }
                }
            }
        }
    }
}
#[cfg(test)]
mod denominator_tests {
    use super::*;
    #[test]
    fn paired_denominators_preserve_literal_results_and_odd_tail() {
        let ids = vec![1usize, 3, 3, 2, 4, 1, 2];
        let compact: Vec<_> = ids.iter().map(|&id| id as u32).collect();
        let rows = [
            crate::quant::EmMultiEc {
                transcript_range: core::range::Range { start: 0, end: 4 },
                count: 3.0,
            },
            crate::quant::EmMultiEc {
                transcript_range: core::range::Range { start: 4, end: 7 },
                count: 5.0,
            },
        ];
        for pattern in [
            [1.0, 2.0, 3.0, 4.0],
            [-0.0, f64::from_bits(1), f64::MAX, 1.0],
            [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0],
        ] {
            let pattern = std::hint::black_box(pattern);
            let weights: Vec<_> = (0..ids.len()).map(|i| pattern[i % 4]).collect();
            for alpha in [
                [1.0, 2.0, 3.0, 4.0],
                [0.0, -0.0, f64::MAX, 1.0],
                [
                    f64::INFINITY,
                    0.0,
                    f64::from_bits(0x7ff8_0000_0000_0023),
                    -1.0,
                ],
            ] {
                let alpha = std::hint::black_box(alpha);
                for len in [1, 2] {
                    let mut expected = [0.0; 2];
                    for (i, row) in rows[..len].iter().enumerate() {
                        for edge in row.transcript_range {
                            expected[i] += alpha[ids[edge] - 1] * weights[edge];
                        }
                    }
                    for actual in [
                        paired_denominators(&rows[..len], &ids, &weights, &alpha, 1),
                        paired_denominators(&rows[..len], &compact, &weights, &alpha, 1),
                    ] {
                        for (actual, expected) in actual.into_iter().zip(expected) {
                            // Arithmetic NaN payload/sign are unspecified by Rust.
                            // Every non-NaN result, including signed zero, stays exact.
                            if expected.is_nan() {
                                assert!(actual.is_nan());
                            } else {
                                assert_eq!(actual.to_bits(), expected.to_bits());
                            }
                        }
                    }
                }
            }
        }
    }
}
