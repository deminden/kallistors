//! Exact component-local scalar scheduling and convergence witnesses.

use super::super::IncomingEcId;
use super::{EmComponent, EmWork};

/// Skip only components whose entire alpha vector is bitwise unchanged.
/// Their weights and counts are immutable and no other component contributes
/// to their transcripts. Global final pruning resets these certificates.
pub(in crate::quant) fn scalar_step(
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

// Static policies share the controller and update loop. Bootstrap must keep
// literal operand order whenever the current component has exceptional values.
#[allow(clippy::too_many_arguments)]
fn component_step<I: IncomingEcId, const STRICT: bool, const MODE: u8>(
    group: &mut EmComponent,
    work: &EmWork,
    ids: &[I],
    alpha: &mut [f64],
    next: &mut [f64],
    products: &mut [f64],
    finite_weights: bool,
) -> (usize, bool) {
    if STRICT {
        let (changes, identical, finite) = if finite_weights && group.scalar_finite {
            work.scalar_dense_step_certified::<_, MODE, false>(
                ids,
                group.ecs,
                group.transcripts,
                alpha,
                next,
                products,
            )
        } else {
            work.scalar_dense_step_certified::<_, MODE, true>(
                ids,
                group.ecs,
                group.transcripts,
                alpha,
                next,
                products,
            )
        };
        group.scalar_finite = finite;
        (changes, identical)
    } else {
        work.scalar_dense_step_mode::<MODE>(group.ecs, group.transcripts, alpha, next, products)
    }
}

#[allow(clippy::too_many_arguments)]
pub(in crate::quant) fn bootstrap_tile<I: IncomingEcId>(
    groups: &mut [EmComponent],
    work: &EmWork,
    ids: &[I],
    active: &[usize],
    alpha: &mut [f64],
    next: &mut [f64],
    products: &mut [f64],
    first_iter: usize,
    min_rounds: usize,
    remaining: usize,
    finite_weights: bool,
) -> (usize, usize) {
    scalar_tile_policy::<_, true>(
        groups,
        work,
        ids,
        active,
        alpha,
        next,
        products,
        first_iter,
        min_rounds,
        remaining,
        finite_weights,
    )
}

// A nonzero component count proves the global count is nonzero. Advance one
// witness to its first eligible zero (or the bounded cap), then advance every
// other independent component to that same boundary. No possible global stop
// is crossed, so neither an alpha snapshot nor speculative replay is needed.
#[allow(clippy::too_many_arguments)]
pub(in crate::quant) fn scalar_tile(
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
    scalar_tile_policy::<_, false>(
        groups,
        work,
        &work.multi_transcripts,
        active,
        alpha,
        next_alpha,
        products,
        first_iter,
        min_rounds,
        remaining,
        true,
    )
}

#[allow(clippy::too_many_arguments)]
fn scalar_tile_policy<I: IncomingEcId, const STRICT: bool>(
    groups: &mut [EmComponent],
    work: &EmWork,
    ids: &[I],
    active: &[usize],
    alpha: &mut [f64],
    next_alpha: &mut [f64],
    products: &mut [f64],
    first_iter: usize,
    min_rounds: usize,
    remaining: usize,
    finite_weights: bool,
) -> (usize, usize) {
    debug_assert!(active.iter().copied().eq(0..alpha.len()));
    let cap = remaining.min(64);
    debug_assert!(cap > 0);
    if cap == 1 {
        let mut changes = 0;
        for group in groups {
            if !group.fixed {
                let (count, identical) = component_step::<_, STRICT, 2>(
                    group,
                    work,
                    ids,
                    alpha,
                    next_alpha,
                    products,
                    finite_weights,
                );
                changes += count;
                group.last_changed = count != 0;
                group.fixed = identical;
            }
        }
        return (1, changes);
    }
    let witness = super::witness(groups);
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
            component_step::<_, STRICT, 2>(
                group,
                work,
                ids,
                alpha,
                next_alpha,
                products,
                finite_weights,
            )
        } else if eligible {
            component_step::<_, STRICT, 1>(
                group,
                work,
                ids,
                alpha,
                next_alpha,
                products,
                finite_weights,
            )
        } else {
            component_step::<_, STRICT, 0>(
                group,
                work,
                ids,
                alpha,
                next_alpha,
                products,
                finite_weights,
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
                component_step::<_, STRICT, 2>(
                    group,
                    work,
                    ids,
                    alpha,
                    next_alpha,
                    products,
                    finite_weights,
                )
            } else {
                component_step::<_, STRICT, 0>(
                    group,
                    work,
                    ids,
                    alpha,
                    next_alpha,
                    products,
                    finite_weights,
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
