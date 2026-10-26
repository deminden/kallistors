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
                            if offset + 1 < taken && first_iter + completed + offset > min_rounds {
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

#[test]
fn certified_scalar_components_track_nonfinite_fixed_groups_and_resets() {
    let classes = vec![vec![0], vec![1], vec![0, 1], vec![2, 3], vec![4]];
    let mut work = super::super::prepare_em_work(
        &classes,
        &[5, 7, 11, 13, 3],
        &[1.0, 2.0, f64::NAN, 1.0, 1.0],
    );
    let mut groups = work.group_components();
    work.compact_transcript_ids();
    let mut alpha = vec![0.2; work.singleton_alpha.len()];
    let mut next = vec![0.0; alpha.len()];
    let active: Vec<_> = (0..alpha.len()).collect();
    let mut products = vec![0.0; 8];
    for _ in 0..4 {
        let _ = bootstrap_tile(
            &mut groups,
            &work,
            &work.multi_transcripts,
            &active,
            &mut alpha,
            &mut next,
            &mut products,
            0,
            0,
            1,
            false,
        );
        assert_eq!(
            groups.iter().all(|g| g.scalar_finite),
            alpha.iter().all(|value| value.is_finite())
        );
    }
    assert!(
        groups
            .iter()
            .any(|group| group.fixed && group.scalar_finite)
    );
    assert!(
        groups
            .iter()
            .any(|group| group.fixed && !group.scalar_finite)
    );
    alpha.fill(0.0);
    reset_fixed(&mut groups);
    let _ = bootstrap_tile(
        &mut groups,
        &work,
        &work.multi_transcripts,
        &active,
        &mut alpha,
        &mut next,
        &mut products,
        0,
        0,
        1,
        false,
    );
    assert_eq!(
        groups.iter().all(|g| g.scalar_finite),
        alpha.iter().all(|value| value.is_finite())
    );
    assert!(alpha.iter().any(|value| value.is_nan()));
}

#[test]
fn parallel_witness_preserves_first_stop_tree_ownership_and_reset() {
    use super::super::{EcScales, EmIncoming, IncomingEcIds, prepare_em_work};
    for threads in [1, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for parallel in [false, true] {
            for wide in [false, true] {
                for minimum in [0, 50] {
                    let mut classes: Vec<_> = (0..70)
                        .map(|i| vec![(2 * i) as u32, (2 * i + 1) as u32])
                        .collect();
                    // Leading largest component exposes incorrect depth-bit
                    // ownership if a singleton subtree is used as witness.
                    classes.extend((0..8).map(|_| vec![0, 1]));
                    let counts = vec![1; classes.len()];
                    let lengths: Vec<_> = (0..70).flat_map(|_| [1.0, 2.0]).collect();
                    let mut work = prepare_em_work(&classes, &counts, &lengths);
                    let mut actual_groups = work.group_components();
                    let mut reference_groups = work.group_components();
                    work.compact_transcript_ids();
                    for group in actual_groups.iter_mut().chain(&mut reference_groups) {
                        group.parallel = parallel;
                    }
                    let mut incoming = EmIncoming::new(&work, false);
                    if wide && let IncomingEcIds::Compact(ids) = &incoming.ec_indices {
                        incoming.ec_indices =
                            IncomingEcIds::Wide(ids.iter().map(|&x| x as usize).collect());
                    }
                    let mut actual_scales = EcScales::new(work.multi_ecs.len());
                    let mut reference_scales = EcScales::new(work.multi_ecs.len());
                    for epoch in 0..2 {
                        let mut actual: Vec<_> = work
                            .active_transcripts
                            .iter()
                            .map(|&id| {
                                if id / 2 != 0 && id / 2 < 35 {
                                    if id % 2 == 0 { 1.0 } else { 0.0 }
                                } else {
                                    0.5
                                }
                            })
                            .collect();
                        let mut expected = actual.clone();
                        reset_fixed(&mut actual_groups);
                        reset_fixed(&mut reference_groups);
                        let mut completed = 0;
                        while completed < 200 {
                            let (taken, changed) = parallel::tile(
                                &pool,
                                &mut actual_groups,
                                &work,
                                &incoming,
                                &work.multi_transcripts,
                                &mut actual,
                                &mut actual_scales,
                                completed,
                                minimum,
                                200 - completed,
                            );
                            assert!(taken > 0 && taken <= 64 && taken <= 200 - completed);
                            let mut expected_changed = 0;
                            for offset in 0..taken {
                                expected_changed = step::<_, true>(
                                    &pool,
                                    &mut reference_groups,
                                    &work,
                                    &incoming,
                                    &work.multi_transcripts,
                                    &mut expected,
                                    &mut reference_scales,
                                    1,
                                );
                                if offset + 1 < taken && completed + offset > minimum {
                                    assert_ne!(
                                        expected_changed, 0,
                                        "crossed first global stop, epoch{epoch}"
                                    );
                                }
                            }
                            assert_eq!(changed, expected_changed);
                            assert_eq!(
                                actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                                expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                            );
                            completed += taken;
                            if changed == 0 && completed - 1 > minimum {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn zero_incoming_certificate_preserves_signed_zero_exceptional_and_refresh() {
    use super::super::{EcScales, EmIncoming, IncomingEcIds, prepare_em_work};
    // Force both scheduling forms and both incoming ID widths. The independent
    // original incoming equations explicitly retain every multiplication.
    for threads in [1, 2, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for parallel in [false, true] {
            for wide in [false, true] {
                for length in [1.0, -2.0, 0.0, f64::INFINITY, f64::NAN] {
                    let classes = vec![vec![0, 1], vec![0, 1], vec![2, 3], vec![3]];
                    let mut work =
                        prepare_em_work(&classes, &[5, 7, 11, 3], &[length, 2.0, 1.0, 3.0]);
                    let mut groups = work.group_components();
                    work.compact_transcript_ids();
                    for group in &mut groups {
                        group.parallel = parallel;
                    }
                    let mut incoming = EmIncoming::new(&work, true);
                    if wide && let IncomingEcIds::Compact(ids) = &incoming.ec_indices {
                        incoming.ec_indices =
                            IncomingEcIds::Wide(ids.iter().map(|&id| id as usize).collect());
                    }
                    let mut alpha: Vec<f64> = work
                        .active_transcripts
                        .iter()
                        .map(|&id| match id {
                            0 => 0.0,
                            1 => 1.0,
                            2 => -0.0,
                            _ => 0.5,
                        })
                        .collect();
                    let mut scales = EcScales::new(work.multi_ecs.len());
                    for iteration in 0..5 {
                        if iteration == 2 {
                            work.multi_weights.fill(1.0);
                            incoming.refresh_weights(&work);
                            reset_fixed(&mut groups);
                            alpha.fill(0.0);
                        }
                        if iteration == 3 {
                            work.multi_weights[0] = f64::INFINITY;
                            incoming.refresh_weights(&work);
                            reset_fixed(&mut groups);
                        }
                        let mut values = vec![0.0; work.multi_ecs.len()];
                        let mut present = vec![false; values.len()];
                        for (i, ec) in work.multi_ecs.iter().enumerate() {
                            let mut denominator = 0.0;
                            for edge in ec.transcript_range {
                                denominator +=
                                    alpha[work.multi_transcripts[edge]] * work.multi_weights[edge];
                            }
                            if denominator >= super::super::TOLERANCE || denominator.is_nan() {
                                present[i] = true;
                                values[i] = ec.count / denominator;
                            }
                        }
                        let mut expected = work.singleton_alpha.clone();
                        for tr in 0..alpha.len() {
                            for (i, ec) in work.multi_ecs.iter().enumerate() {
                                for edge in ec.transcript_range {
                                    if work.multi_transcripts[edge] == tr && present[i] {
                                        expected[tr] +=
                                            work.multi_weights[edge] * alpha[tr] * values[i];
                                    }
                                }
                            }
                        }
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
                            alpha.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                            expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                        );
                    }
                }
            }
        }
    }
}
