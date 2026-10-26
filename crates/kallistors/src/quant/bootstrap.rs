//! Bootstrap resampling and quantification with fixed effective lengths.

use super::*;

mod sampling;
use sampling::SamplingPlan;

#[inline(never)]
pub fn bootstrap_quantify(
    input: &EcCountsInput,
    eff_lens: &[f64],
    bootstrap_samples: usize,
    seed: u64,
    options: QuantOptions,
) -> Result<Vec<QuantResult>> {
    if input.ec_list.classes.len() != input.counts.len() {
        return Err(Error::InvalidFormat(
            "EC list and counts length mismatch".into(),
        ));
    }

    if bootstrap_samples == 0 {
        return Ok(Vec::new());
    }
    let sampling = SamplingPlan::new(&input.counts).with_exact_inverse();
    let mut seed_rng = SplitMix64::new(seed);
    let mut out = Vec::with_capacity(bootstrap_samples);
    for _ in 0..bootstrap_samples {
        let mut sample_rng = SplitMix64::new(seed_rng.next_u64());
        let counts = sampling.sample(&mut sample_rng);
        out.push(em_quantify_fixed_eff_lens(
            &input.ec_list.classes,
            &counts,
            eff_lens,
            options.max_iter,
            options.min_rounds,
        ));
    }
    Ok(out)
}

const BOOTSTRAP_SCRATCH_BUDGET: usize = 2 * 1024 * 1024 * 1024;

/// Quantify independent bootstrap samples with bounded requested workers.
/// Seeds and collected results retain the serial sample order.
pub fn bootstrap_quantify_threaded(
    input: &EcCountsInput,
    eff_lens: &[f64],
    bootstrap_samples: usize,
    seed: u64,
    options: QuantOptions,
    threads: usize,
) -> Result<Vec<QuantResult>> {
    if threads <= 1 || bootstrap_samples == 0 {
        return bootstrap_quantify(input, eff_lens, bootstrap_samples, seed, options);
    }
    if input.ec_list.classes.len() != input.counts.len() {
        return Err(Error::InvalidFormat(
            "EC list and counts length mismatch".into(),
        ));
    }
    let edges = input
        .ec_list
        .classes
        .iter()
        .fold(0usize, |sum, ec| sum.saturating_add(ec.len()));
    let scratch = bootstrap_scratch_bytes(edges, input.counts.len(), eff_lens.len());
    // Independent-sample scheduling cannot use spare workers for one sample.
    // Large models reuse the ordinary ordered component engine in one pool;
    // keep tiny and oversized cases on the existing serial route.
    if bootstrap_samples == 1 && edges >= 262_144 && scratch <= BOOTSTRAP_SCRATCH_BUDGET {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|error| {
                Error::InvalidFormat(format!("failed to build bootstrap thread pool: {error}"))
            })?;
        let sampling = SamplingPlan::new(&input.counts);
        let mut seed_rng = SplitMix64::new(seed);
        let mut sample_rng = SplitMix64::new(seed_rng.next_u64());
        let counts = sampling.sample_with_pool(&mut sample_rng, &pool);
        return Ok(vec![em_quantify_fixed_eff_lens_with_pool::<true>(
            &input.ec_list.classes,
            &counts,
            eff_lens,
            options.max_iter,
            options.min_rounds,
            Some(&pool),
        )]);
    }
    let workers = bootstrap_workers(scratch, threads, bootstrap_samples);
    if workers <= 1 {
        return bootstrap_quantify(input, eff_lens, bootstrap_samples, seed, options);
    }
    let parallel_within = bootstrap_inner_parallel(
        edges,
        input.counts.len(),
        eff_lens.len(),
        bootstrap_samples,
        threads,
        workers,
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(if parallel_within { threads } else { workers })
        .build()
        .map_err(|error| {
            Error::InvalidFormat(format!("failed to build bootstrap thread pool: {error}"))
        })?;
    let sampling = SamplingPlan::new(&input.counts);
    let mut seed_rng = SplitMix64::new(seed);
    let seeds: Vec<_> = (0..bootstrap_samples)
        .map(|_| seed_rng.next_u64())
        .collect();
    // Indexed collection keeps sample numbers stable; each sample owns exactly
    // the same RNG and scalar equations as the unchanged serial API.
    Ok(pool.install(|| {
        seeds
            .into_par_iter()
            .map(|sample_seed| {
                let mut sample_rng = SplitMix64::new(sample_seed);
                let counts = if parallel_within {
                    sampling.sample_with_pool(&mut sample_rng, &pool)
                } else {
                    sampling.sample(&mut sample_rng)
                };
                if parallel_within {
                    // All samples share this pool: nested component work steals
                    // spare workers without creating a pool per sample.
                    em_quantify_fixed_eff_lens_with_pool::<true>(
                        &input.ec_list.classes,
                        &counts,
                        eff_lens,
                        options.max_iter,
                        options.min_rounds,
                        Some(&pool),
                    )
                } else {
                    em_quantify_fixed_eff_lens(
                        &input.ec_list.classes,
                        &counts,
                        eff_lens,
                        options.max_iter,
                        options.min_rounds,
                    )
                }
            })
            .collect()
    }))
}

fn bootstrap_scratch_bytes(edges: usize, classes: usize, transcripts: usize) -> usize {
    // Conservative capacity/reallocation allowances for per-sample prepared
    // edges, EC records, sampling arrays and alpha vectors. Shared inputs and
    // returned samples are not duplicated into this worker scratch.
    edges
        .saturating_mul(48)
        .saturating_add(classes.saturating_mul(128))
        .saturating_add(transcripts.saturating_mul(64))
        .saturating_add(1 << 20)
}

fn bootstrap_inner_parallel(
    edges: usize,
    classes: usize,
    transcripts: usize,
    samples: usize,
    threads: usize,
    workers: usize,
) -> bool {
    let scratch = bootstrap_scratch_bytes(edges, classes, transcripts)
        .saturating_add(edges.saturating_mul(12))
        .saturating_add(transcripts.saturating_mul(16))
        .saturating_add(classes.saturating_mul(9));
    workers > 1
        && workers == samples
        && workers < threads
        && edges >= 262_144
        && scratch.saturating_mul(workers) <= BOOTSTRAP_SCRATCH_BUDGET
}

fn bootstrap_workers(scratch: usize, requested: usize, samples: usize) -> usize {
    let memory_workers = (BOOTSTRAP_SCRATCH_BUDGET / scratch.max(1)).max(1);
    requested.max(1).min(samples).min(memory_workers)
}

pub(super) fn em_quantify_fixed_eff_lens(
    ec_classes: &[Vec<u32>],
    counts: &[u32],
    eff_lens: &[f64],
    max_iter: usize,
    min_rounds: usize,
) -> QuantResult {
    em_quantify_fixed_eff_lens_with_pool::<false>(
        ec_classes, counts, eff_lens, max_iter, min_rounds, None,
    )
}

// Scheduling specializations share the same bootstrap prior, stop and final
// update. Arithmetic stays in the existing component kernels.
#[inline(never)]
fn em_quantify_fixed_eff_lens_with_pool<const PARALLEL: bool>(
    ec_classes: &[Vec<u32>],
    counts: &[u32],
    eff_lens: &[f64],
    max_iter: usize,
    min_rounds: usize,
    pool: Option<&rayon::ThreadPool>,
) -> QuantResult {
    let num_trans = eff_lens.len();
    let mut em_work = prepare_em_work(ec_classes, counts, eff_lens);
    // The zero-iteration result includes the uniform prior for inactive IDs.
    if max_iter == 0 {
        let alpha = vec![1.0 / num_trans as f64; num_trans];
        let tpm = counts_to_tpm(&alpha, eff_lens);
        return QuantResult {
            est_counts: alpha,
            tpm,
            eff_lengths: eff_lens.to_vec(),
            post_bias: None,
        };
    }
    let mut groups = em_work.group_components();
    if PARALLEL {
        components::configure_parallel(
            &mut groups,
            &em_work,
            pool.expect("parallel bootstrap owns its requested pool")
                .current_num_threads(),
        );
    }
    em_work.compact_transcript_ids();
    let scalar_ids = if !PARALLEL
        && em_work.multi_transcripts.len() >= 262_144
        && em_work.active_transcripts.len() <= u32::MAX as usize
    {
        Some(CompactTranscriptIds::new(std::mem::take(
            &mut em_work.multi_transcripts,
        )))
    } else {
        None
    };
    let mut alpha = vec![1.0 / num_trans as f64; em_work.singleton_alpha.len()];
    let mut next_alpha = if PARALLEL {
        Vec::new()
    } else {
        vec![0.0; alpha.len()]
    };
    let dense_active: Vec<_> = if PARALLEL {
        Vec::new()
    } else {
        (0..alpha.len()).collect()
    };
    let incoming = if PARALLEL {
        Some(EmIncoming::new(&em_work, false))
    } else {
        None
    };
    let mut scales = EcScales::new(if PARALLEL { em_work.multi_ecs.len() } else { 0 });

    let finite_weights = em_work.multi_weights.iter().all(|value| value.is_finite());
    let mut products = if !PARALLEL && finite_weights && max_iter != 0 {
        vec![
            0.0;
            em_work
                .multi_ecs
                .iter()
                .map(|ec| ec.transcript_range.end - ec.transcript_range.start)
                .max()
                .unwrap_or(0)
                .saturating_mul(2)
        ]
    } else {
        Vec::new()
    };
    // No global stop can occur at or before min_rounds. Reuse the main
    // engine's exact component-local mandatory rounds without change counts.
    let mut iter = if PARALLEL {
        let rounds = max_iter.min(min_rounds.saturating_add(1));
        components::step::<_, false>(
            pool.expect("parallel bootstrap owns its requested pool"),
            &mut groups,
            &em_work,
            incoming
                .as_ref()
                .expect("parallel bootstrap owns incoming edges"),
            &em_work.multi_transcripts,
            &mut alpha,
            &mut scales,
            rounds,
        );
        rounds
    } else {
        0
    };
    let mut final_round = false;
    while iter < max_iter {
        let remaining = if final_round { 1 } else { max_iter - iter };
        let (taken, chcount) = if PARALLEL {
            components::parallel_tile(
                pool.expect("parallel bootstrap owns its requested pool"),
                &mut groups,
                &em_work,
                incoming
                    .as_ref()
                    .expect("parallel bootstrap owns incoming edges"),
                &em_work.multi_transcripts,
                &mut alpha,
                &mut scales,
                iter,
                min_rounds,
                remaining,
            )
        } else {
            match &scalar_ids {
                Some(ids) => components::bootstrap_tile(
                    &mut groups,
                    &em_work,
                    ids,
                    &dense_active,
                    &mut alpha,
                    &mut next_alpha,
                    &mut products,
                    iter,
                    min_rounds,
                    remaining,
                    finite_weights,
                ),
                None => components::bootstrap_tile(
                    &mut groups,
                    &em_work,
                    &em_work.multi_transcripts,
                    &dense_active,
                    &mut alpha,
                    &mut next_alpha,
                    &mut products,
                    iter,
                    min_rounds,
                    remaining,
                    finite_weights,
                ),
            }
        };
        iter += taken;

        let stop_em = chcount == 0 && iter - 1 > min_rounds;
        if final_round {
            break;
        }
        if stop_em {
            final_round = true;
            for val in &mut alpha {
                if *val < ALPHA_LIMIT / 10.0 {
                    *val = 0.0;
                }
            }
            components::reset_fixed(&mut groups);
        }
    }

    let mut scattered = vec![0.0; num_trans];
    for (&original, value) in em_work.active_transcripts.iter().zip(alpha) {
        scattered[original] = value;
    }
    let alpha = scattered;
    let tpm = counts_to_tpm(&alpha, eff_lens);
    QuantResult {
        est_counts: alpha,
        tpm,
        eff_lengths: eff_lens.to_vec(),
        post_bias: None,
    }
}

pub(super) struct SplitMix64 {
    pub(super) state: u64,
}

impl SplitMix64 {
    pub(super) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn advance(&mut self, draws: u64) {
        self.state = self
            .state
            .wrapping_add(0x9e37_79b9_7f4a_7c15u64.wrapping_mul(draws));
    }

    pub(super) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bootstrap_workers_bound_requested_samples_and_saturating_scratch() {
        for requested in [0usize, 1, 2, 3, 8, 64, usize::MAX] {
            for samples in [0usize, 1, 2, 7, 100, usize::MAX] {
                for bytes in [0usize, 1, 1 << 20, 1 << 28, 1 << 30, usize::MAX] {
                    let workers = bootstrap_workers(bytes, requested, samples);
                    assert!(workers <= requested.max(1) && workers <= samples);
                    if samples != 0 {
                        assert!(workers >= 1);
                    }
                    if workers > 1 {
                        assert!(workers.saturating_mul(bytes) <= BOOTSTRAP_SCRATCH_BUDGET);
                    }
                }
            }
        }
        assert_eq!(bootstrap_workers(usize::MAX, 8, 100), 1);
        assert_eq!(
            bootstrap_scratch_bytes(usize::MAX, usize::MAX, usize::MAX),
            usize::MAX
        );
        assert_eq!(
            bootstrap_workers(bootstrap_scratch_bytes(16, 4, 8), 8, 7),
            7
        );
        assert!(bootstrap_inner_parallel(300_000, 100_000, 670_670, 4, 8, 4));
        assert!(!bootstrap_inner_parallel(
            300_000, 100_000, 670_670, 4, 4, 4
        ));
        assert!(!bootstrap_inner_parallel(
            300_000, 100_000, 670_670, 4, 8, 2
        ));
        assert!(!bootstrap_inner_parallel(16, 4, 8, 4, 8, 4));
        assert!(!bootstrap_inner_parallel(
            usize::MAX,
            usize::MAX,
            usize::MAX,
            4,
            8,
            4
        ));
        for threads in [1, 2, 8, 16, 64] {
            for samples in [0, 1, 2, 4, 8, 32, 100] {
                let workers = bootstrap_workers(1 << 20, threads, samples);
                if bootstrap_inner_parallel(300_000, 100_000, 670_670, samples, threads, workers) {
                    assert!(workers > 1 && workers == samples && samples < threads);
                }
            }
        }
    }

    #[test]
    fn bootstrap_shared_kernel_and_exceptional_step_match_ordered_equations() {
        fn step<I: IncomingEcId>(
            work: &EmWork,
            ids: &[I],
            prior: &[f64],
            reuse: bool,
        ) -> (Vec<f64>, bool) {
            let mut alpha = prior.to_vec();
            let mut next = work.singleton_alpha.clone();
            let mut products = vec![0.0; 16];
            let ecs = Range {
                start: 0,
                end: work.multi_ecs.len(),
            };
            let active = Range {
                start: 0,
                end: alpha.len(),
            };
            let (_, _, finite) = if reuse {
                work.scalar_dense_step_certified::<_, 2, false>(
                    ids,
                    ecs,
                    active,
                    &mut alpha,
                    &mut next,
                    &mut products,
                )
            } else {
                work.scalar_dense_step_certified::<_, 2, true>(
                    ids,
                    ecs,
                    active,
                    &mut alpha,
                    &mut next,
                    &mut products,
                )
            };
            (alpha, finite)
        }
        let classes = vec![vec![0, 1], vec![1, 2], vec![0, 2, 3, 3], vec![1]];
        let counts = [3, 7, 11, 5];
        let weights = [
            vec![1.0],
            vec![-0.0, 1.0],
            vec![f64::from_bits(1), 1.0],
            vec![f64::MAX],
            vec![f64::INFINITY, 1.0],
            vec![f64::NEG_INFINITY, -0.0],
            vec![f64::from_bits(0x7ff8_0000_0000_0023), 1.0],
        ];
        let initial = [
            [0.0, 1.0, 2.0, 3.0],
            [-0.0, 1.0, 2.0, f64::MAX],
            [f64::from_bits(1), f64::MAX, f64::MAX, -1.0],
            [f64::MAX; 4],
            [f64::from_bits(0x7ff8_0000_0000_0049), 1.0, 0.0, 2.0],
            [f64::INFINITY, 0.0, f64::NEG_INFINITY, -0.0],
        ];
        let mut exercised = [false; 2];
        for pattern in weights {
            let mut work = prepare_em_work(&classes, &counts, &[1.0; 4]);
            for (i, weight) in work.multi_weights.iter_mut().enumerate() {
                *weight = pattern[i % pattern.len()];
            }
            for alpha in initial {
                let mut expected = work.singleton_alpha.clone();
                // Independent literal ordered equations; no product helper or
                // finite certificate is shared with the implementation.
                for ec in &work.multi_ecs {
                    let mut denominator = 0.0;
                    for edge in ec.transcript_range {
                        denominator +=
                            alpha[work.multi_transcripts[edge]] * work.multi_weights[edge];
                    }
                    if denominator < TOLERANCE {
                        continue;
                    }
                    let scale = ec.count / denominator;
                    for edge in ec.transcript_range {
                        let tr = work.multi_transcripts[edge];
                        expected[tr] += work.multi_weights[edge] * alpha[tr] * scale;
                    }
                }
                let reuse = work.multi_weights.iter().all(|x| x.is_finite())
                    && alpha.iter().all(|x| x.is_finite());
                exercised[usize::from(reuse)] = true;
                let narrow: Vec<u32> = work.multi_transcripts.iter().map(|&id| id as u32).collect();
                for (actual, finite) in [
                    step(&work, &work.multi_transcripts, &alpha, reuse),
                    step(&work, &narrow, &alpha, reuse),
                ] {
                    assert_eq!(finite, expected.iter().all(|value| value.is_finite()));
                    for (actual, expected) in actual.iter().zip(&expected) {
                        // Rust leaves arithmetic NaN payload/sign unspecified.
                        // All other results, including infinities and -0, stay bit-exact.
                        if expected.is_nan() {
                            assert!(actual.is_nan());
                        } else {
                            assert_eq!(actual.to_bits(), expected.to_bits());
                        }
                    }
                }
            }
        }
        assert_eq!(exercised, [true, true]);
    }

    #[test]
    fn compact_bootstrap_preserves_original_prior_and_scattered_zero_targets() {
        // IDs are intentionally sparse. Equal effective lengths make the
        // symmetric two-target solution independent of the convergence route.
        let classes = vec![vec![1, 5], vec![3]];
        let lengths = [10.0; 8];
        for limit in [0, 1, 2, 500] {
            let result = em_quantify_fixed_eff_lens(&classes, &[16, 4], &lengths, limit, 0);
            let expected: Vec<f64> = if limit == 0 {
                vec![0.125; 8]
            } else {
                vec![0.0, 8.0, 0.0, 4.0, 0.0, 8.0, 0.0, 0.0]
            };
            assert_eq!(
                result
                    .est_counts
                    .iter()
                    .map(|x| x.to_bits())
                    .collect::<Vec<_>>(),
                expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
            assert_eq!(result.eff_lengths, lengths);
        }
        let empty = em_quantify_fixed_eff_lens(&classes, &[0, 0], &lengths, 3, 0);
        assert!(empty.est_counts.iter().all(|x| x.to_bits() == 0));
    }
    #[test]
    fn bootstrap_witness_preserves_first_stop_pruning_and_exceptional_order() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        // A complete equation oracle for the scheduler contract, deliberately
        // without component grouping, cached products or fixed-point skipping.
        let classes = vec![
            vec![1],
            vec![1],
            vec![5],
            vec![1, 3],
            vec![3, 7],
            vec![2, 4],
            vec![2, 4],
            vec![6, 8],
            vec![0, 9],
        ];
        let counts = [7, 3, 11, 17, 19, 23, 29, 0, 2];
        for exceptional in [false, true] {
            let mut lengths = [11.0, 31.0, 51.0, 17.0, 8.0, 19.0, 13.0, 91.0, 5.0, 7.0];
            if exceptional {
                lengths[2] = f64::NAN;
                lengths[4] = 0.0;
            }
            for max_iter in [0, 1, 2, 3, 63, 64, 65, 127, 129] {
                for minimum in [0, 1, 50, 128, usize::MAX] {
                    let mut alpha = vec![0.1; lengths.len()];
                    let mut final_round = false;
                    for iteration in 0..max_iter {
                        let mut next = vec![0.0; lengths.len()];
                        for (class, &count) in classes.iter().zip(&counts) {
                            if class.len() == 1 {
                                next[class[0] as usize] = count as f64;
                            }
                        }
                        for (class, &count) in classes.iter().zip(&counts) {
                            if class.len() <= 1 || count == 0 {
                                continue;
                            }
                            let mut denominator = 0.0;
                            for &id in class {
                                let index = id as usize;
                                denominator += alpha[index] * (count as f64 / lengths[index]);
                            }
                            if denominator < TOLERANCE {
                                continue;
                            }
                            let scale = count as f64 / denominator;
                            for &id in class {
                                let index = id as usize;
                                next[index] +=
                                    (count as f64 / lengths[index]) * alpha[index] * scale;
                            }
                        }
                        let changed = alpha.iter().zip(&next).any(|(&old, &new)| {
                            new > ALPHA_CHANGE_LIMIT && (new - old).abs() / new > ALPHA_CHANGE
                        });
                        alpha = next;
                        if final_round {
                            break;
                        }
                        if !changed && iteration > minimum {
                            final_round = true;
                            for value in &mut alpha {
                                if *value < ALPHA_LIMIT / 10.0 {
                                    *value = 0.0;
                                }
                            }
                        }
                    }
                    for result in [
                        em_quantify_fixed_eff_lens(&classes, &counts, &lengths, max_iter, minimum),
                        em_quantify_fixed_eff_lens_with_pool::<true>(
                            &classes,
                            &counts,
                            &lengths,
                            max_iter,
                            minimum,
                            Some(&pool),
                        ),
                    ] {
                        for (actual, expected) in result.est_counts.iter().zip(&alpha) {
                            // Rust arithmetic does not specify NaN payload/sign.
                            // Finite values, infinities and both zero signs stay exact.
                            assert!(
                                (actual.is_nan() && expected.is_nan())
                                    || actual.to_bits() == expected.to_bits(),
                                "exceptional={exceptional}, limit={max_iter}, minimum={minimum}"
                            );
                        }
                    }
                }
            }
        }
    }
}
