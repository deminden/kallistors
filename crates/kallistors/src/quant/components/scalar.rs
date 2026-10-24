//! Exact component-local scalar scheduling and convergence witnesses.

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
