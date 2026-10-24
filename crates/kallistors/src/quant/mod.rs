//! Quantification routines (EM over ECs with optional sequence bias).

mod bootstrap;
mod input;
mod models;
mod work;

pub use bootstrap::bootstrap_quantify;
pub use input::{load_ec_counts, load_transcript_sequences};
use models::{
    calc_eff_lens, empirical_frag_lens_by_transcript, mean_frag_lens_by_transcript, update_eff_lens,
};
#[cfg(test)]
use work::sum_scaled_edges;
use work::{
    EcScales, EmIncoming, EmMultiEc, EmWork, IncomingEcId, IncomingEcIds, prepare_em_work,
    refresh_multi_weights, weighted_alpha_sum,
};

use core::range::Range;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::bias::{BiasCounts, hexamer_to_int, update_hexamer};
use crate::ec::EcList;
use crate::pseudoalign::StrandSpecific;
use crate::{Error, Result};
use rayon::prelude::*;

mod compact_ids;
mod components;
use compact_ids::CompactTranscriptIds;
mod output;
pub use output::{
    H5QuantOutput, RunInfo, write_abundance_h5, write_abundance_tsv,
    write_abundance_tsv_from_parts, write_abundance_tsv_from_parts_threaded,
};

const MIN_ALPHA: f64 = 1e-8;
const ALPHA_LIMIT: f64 = 1e-7;
const ALPHA_CHANGE_LIMIT: f64 = 1e-2;
const ALPHA_CHANGE: f64 = 1e-2;
const TOLERANCE: f64 = f64::from_bits(1);

#[derive(Debug, Clone, Copy)]
pub struct QuantOptions {
    pub mean_fragment_length: f64,
    pub fragment_length_sd: f64,
    pub max_iter: usize,
    pub min_rounds: usize,
    pub bias: bool,
    pub strand_specific: Option<StrandSpecific>,
}

impl Default for QuantOptions {
    fn default() -> Self {
        Self {
            mean_fragment_length: 200.0,
            fragment_length_sd: 20.0,
            max_iter: 10_000,
            min_rounds: 50,
            bias: false,
            strand_specific: None,
        }
    }
}

pub struct QuantResult {
    pub est_counts: Vec<f64>,
    pub tpm: Vec<f64>,
    pub eff_lengths: Vec<f64>,
    pub post_bias: Option<Vec<f64>>,
}

pub struct EcCountsInput {
    pub ec_list: EcList,
    pub counts: Vec<u32>,
}

pub fn em_quantify(
    input: &EcCountsInput,
    transcript_lengths: &[u32],
    transcript_sequences: Option<&[Vec<u8>]>,
    bias_counts: Option<&BiasCounts>,
    options: QuantOptions,
) -> Result<QuantResult> {
    em_quantify_with_fragment_histogram(
        input,
        transcript_lengths,
        transcript_sequences,
        bias_counts,
        None,
        options,
    )
}

/// Quantify using observed fragment lengths for paired-end input, or the supplied
/// Gaussian fragment model when no empirical histogram is available.
pub fn em_quantify_with_fragment_histogram(
    input: &EcCountsInput,
    transcript_lengths: &[u32],
    transcript_sequences: Option<&[Vec<u8>]>,
    bias_counts: Option<&BiasCounts>,
    fragment_length_hist: Option<&[u32]>,
    options: QuantOptions,
) -> Result<QuantResult> {
    em_quantify_with_fragment_histogram_threaded(
        input,
        transcript_lengths,
        transcript_sequences,
        bias_counts,
        fragment_length_hist,
        options,
        1,
    )
}

/// Quantify with parallel EC denominators and independent transcript updates.
/// Each transcript retains its original contribution order for exact f64 sums.
pub fn em_quantify_with_fragment_histogram_threaded(
    input: &EcCountsInput,
    transcript_lengths: &[u32],
    transcript_sequences: Option<&[Vec<u8>]>,
    bias_counts: Option<&BiasCounts>,
    fragment_length_hist: Option<&[u32]>,
    options: QuantOptions,
    threads: usize,
) -> Result<QuantResult> {
    em_quantify_input_threaded(
        QuantInput::Borrowed(input),
        transcript_lengths,
        transcript_sequences,
        bias_counts,
        fragment_length_hist,
        options,
        threads,
    )
}

/// Quantify owned EC counts, releasing their duplicate transcript lists after
/// preparation. Use the borrowed API when counts are needed for bootstrapping.
pub fn em_quantify_owned_with_fragment_histogram_threaded(
    input: EcCountsInput,
    transcript_lengths: &[u32],
    transcript_sequences: Option<&[Vec<u8>]>,
    bias_counts: Option<&BiasCounts>,
    fragment_length_hist: Option<&[u32]>,
    options: QuantOptions,
    threads: usize,
) -> Result<QuantResult> {
    em_quantify_input_threaded(
        QuantInput::Owned(input),
        transcript_lengths,
        transcript_sequences,
        bias_counts,
        fragment_length_hist,
        options,
        threads,
    )
}

enum QuantInput<'a> {
    Borrowed(&'a EcCountsInput),
    Owned(EcCountsInput),
}
impl std::ops::Deref for QuantInput<'_> {
    type Target = EcCountsInput;
    fn deref(&self) -> &EcCountsInput {
        match self {
            Self::Borrowed(input) => input,
            Self::Owned(input) => input,
        }
    }
}

fn em_quantify_input_threaded(
    input: QuantInput<'_>,
    transcript_lengths: &[u32],
    transcript_sequences: Option<&[Vec<u8>]>,
    bias_counts: Option<&BiasCounts>,
    fragment_length_hist: Option<&[u32]>,
    options: QuantOptions,
    threads: usize,
) -> Result<QuantResult> {
    let num_trans = transcript_lengths.len();
    if num_trans == 0 {
        return Err(Error::InvalidFormat("no transcripts provided".into()));
    }
    if input.ec_list.classes.len() != input.counts.len() {
        return Err(Error::InvalidFormat(
            "EC list and counts length mismatch".into(),
        ));
    }

    let means = if let Some(histogram) = fragment_length_hist {
        empirical_frag_lens_by_transcript(transcript_lengths, histogram)?
    } else {
        mean_frag_lens_by_transcript(
            transcript_lengths,
            options.mean_fragment_length,
            options.fragment_length_sd,
        )
    };
    let mut eff_lens = calc_eff_lens(transcript_lengths, &means);
    let mut post_bias = None;

    let mut alpha = vec![1.0 / num_trans as f64; num_trans];
    let mut next_alpha = vec![0.0f64; num_trans];

    let mut em_work = prepare_em_work(&input.ec_list.classes, &input.counts, &eff_lens);
    drop(input);

    let em_pool = if threads > 1 && !em_work.multi_ecs.is_empty() {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|err| {
                    Error::InvalidFormat(format!("failed to build EM thread pool: {err}"))
                })?,
        )
    } else {
        None
    };
    let compact = em_pool.is_some() && !options.bias && options.max_iter != 0;
    // Preserve the ordered scalar equations while keeping active alpha values
    // contiguous. Original target IDs are retained for the final scatter.
    let serial_dense = em_pool.is_none()
        && !em_work.multi_ecs.is_empty()
        && !options.bias
        && options.max_iter != 0;
    let mut scalar_components = None;
    let mut components = if compact {
        let mut components = em_work.group_components();
        components::configure_parallel(&mut components, &em_work, threads);
        em_work.compact_transcript_ids();
        alpha.truncate(em_work.active_transcripts.len());
        // Compact updates use only their own previous alpha after the denominator barrier.
        next_alpha = Vec::new();
        Some(components)
    } else if serial_dense {
        // Stable component layout groups independent ECs and active targets.
        // Each target keeps exactly its original contribution order.
        scalar_components = Some(em_work.group_components());
        em_work.compact_transcript_ids();
        alpha.truncate(em_work.active_transcripts.len());
        next_alpha.truncate(em_work.active_transcripts.len());
        None
    } else {
        None
    };
    let dense_active: Vec<usize> = if serial_dense {
        (0..em_work.active_transcripts.len()).collect()
    } else {
        Vec::new()
    };
    let mut incoming = em_pool
        .as_ref()
        .map(|_| EmIncoming::new(&em_work, options.bias));
    let compact_transcripts = if compact && em_work.active_transcripts.len() <= u32::MAX as usize {
        Some(CompactTranscriptIds::new(std::mem::take(
            &mut em_work.multi_transcripts,
        )))
    } else {
        None
    };
    // The scalar update visits each EC twice. Retain its weighted alpha
    // products in one reusable buffer so the second visit avoids gathering
    // alpha and multiplying it by the same weight again.
    let mut scalar_products = if em_pool.is_none() {
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
    let mut scales = EcScales::new(if em_pool.is_some() {
        em_work.multi_ecs.len()
    } else {
        0
    });
    let mut next_active = if em_pool.is_some() && !compact {
        vec![0.0; em_work.active_transcripts.len()]
    } else {
        Vec::new()
    };
    // Convergence cannot terminate at iter <= min_rounds. Independent
    // components can execute these mandatory rounds while their data is hot.
    let burn_in = if let (Some(pool), Some(incoming), Some(components)) =
        (&em_pool, &incoming, &mut components)
    {
        let rounds = options.max_iter.min(options.min_rounds.saturating_add(1));
        match &compact_transcripts {
            Some(ids) => {
                components::step::<_, false>(
                    pool,
                    components,
                    &em_work,
                    incoming,
                    ids,
                    &mut alpha,
                    &mut scales,
                    rounds,
                );
            }
            None => {
                components::step::<_, false>(
                    pool,
                    components,
                    &em_work,
                    incoming,
                    &em_work.multi_transcripts,
                    &mut alpha,
                    &mut scales,
                    rounds,
                );
            }
        }
        rounds
    } else {
        0
    };
    // Keep the existing model-size eligibility. Witness scheduling needs
    // no speculative alpha snapshot; parallel and bias routes are unchanged.
    let scalar_tiles = serial_dense && !alpha.is_empty() && alpha.len() <= 4 * 1024 * 1024 / 8;
    let mut final_round = false;
    let mut iterations = burn_in;
    let mut iter = burn_in;
    while iter < options.max_iter {
        iterations = iter + 1;
        if options.bias
            && (iter == options.min_rounds || iter == options.min_rounds + 500)
            && let (Some(bias_counts), Some(transcript_sequences)) =
                (bias_counts, transcript_sequences)
        {
            let (updated_eff, bias_model) = update_eff_lens(
                &means,
                bias_counts,
                transcript_sequences,
                &alpha,
                &eff_lens,
                options.strand_specific,
            );
            eff_lens = updated_eff;
            post_bias = Some(bias_model);
            refresh_multi_weights(&mut em_work, &eff_lens);
            if let Some(incoming) = &mut incoming {
                incoming.refresh_weights(&em_work);
            }
        }

        let mut step_changes = None;
        if let (Some(pool), Some(incoming), Some(components)) =
            (&em_pool, &incoming, &mut components)
        {
            step_changes = Some(match &compact_transcripts {
                Some(ids) => components::step::<_, true>(
                    pool,
                    components,
                    &em_work,
                    incoming,
                    ids,
                    &mut alpha,
                    &mut scales,
                    1,
                ),
                None => components::step::<_, true>(
                    pool,
                    components,
                    &em_work,
                    incoming,
                    &em_work.multi_transcripts,
                    &mut alpha,
                    &mut scales,
                    1,
                ),
            });
        } else if let (Some(pool), Some(incoming)) = (&em_pool, &incoming) {
            pool.install(|| {
                scales
                    .values
                    .par_iter_mut()
                    .zip(&mut scales.present)
                    .zip(&em_work.multi_ecs)
                    .for_each(|((scale, present), ec)| {
                        let weights = &em_work.multi_weights[ec.transcript_range];
                        let denom = match &compact_transcripts {
                            Some(ids) => {
                                weighted_alpha_sum(&ids[ec.transcript_range], weights, &alpha)
                            }
                            None => weighted_alpha_sum(
                                &em_work.multi_transcripts[ec.transcript_range],
                                weights,
                                &alpha,
                            ),
                        };
                        if denom < TOLERANCE {
                            *present = 0;
                            *scale = 0.0;
                        } else {
                            *present = 1;
                            *scale = ec.count / denom;
                        }
                    });
            });
            if compact {
                step_changes = Some(pool.install(|| {
                    alpha
                        .par_iter_mut()
                        .enumerate()
                        .map(|(tr, value)| {
                            let old = *value;
                            let next = incoming.transcript_sum(
                                tr,
                                em_work.singleton_alpha[tr],
                                old,
                                &scales,
                            );
                            *value = next;
                            usize::from(
                                next > ALPHA_CHANGE_LIMIT
                                    && (next - old).abs() / next > ALPHA_CHANGE,
                            )
                        })
                        .sum::<usize>()
                }));
            } else {
                pool.install(|| {
                    next_active
                        .par_iter_mut()
                        .zip(&em_work.active_transcripts)
                        .for_each(|(next, &tr)| {
                            *next = incoming.transcript_sum(
                                tr,
                                em_work.singleton_alpha[tr],
                                alpha[tr],
                                &scales,
                            );
                        });
                });
                for (&tr, &next) in em_work.active_transcripts.iter().zip(&next_active) {
                    next_alpha[tr] = next;
                }
            }
        } else {
            let active = if serial_dense {
                &dense_active
            } else {
                &em_work.active_transcripts
            };
            step_changes = Some(if let Some(groups) = &mut scalar_components {
                if !final_round && scalar_tiles {
                    let (taken, changes) = components::scalar_tile(
                        groups,
                        &em_work,
                        active,
                        &mut alpha,
                        &mut next_alpha,
                        &mut scalar_products,
                        iter,
                        options.min_rounds,
                        options.max_iter - iter,
                    );
                    iter += taken - 1;
                    iterations = iter + 1;
                    changes
                } else {
                    components::scalar_step(
                        groups,
                        &em_work,
                        active,
                        &mut alpha,
                        &mut next_alpha,
                        &mut scalar_products,
                    )
                }
            } else {
                em_work
                    .scalar_step(
                        Range {
                            start: 0,
                            end: em_work.multi_ecs.len(),
                        },
                        active,
                        &mut alpha,
                        &mut next_alpha,
                        &mut scalar_products,
                    )
                    .0
            });
            if iter == 0 && !serial_dense {
                // Inactive targets lose their initial uniform prior after iteration zero.
                alpha.copy_from_slice(&next_alpha);
            }
        }

        let chcount = if let Some(changes) = step_changes {
            changes
        } else {
            let mut changes = 0;
            let active = if serial_dense {
                &dense_active
            } else {
                &em_work.active_transcripts
            };
            for &i in active {
                if next_alpha[i] > ALPHA_CHANGE_LIMIT {
                    let rel = (next_alpha[i] - alpha[i]).abs() / next_alpha[i];
                    if rel > ALPHA_CHANGE {
                        changes += 1;
                    }
                }
                alpha[i] = next_alpha[i];
            }
            if iter == 0 {
                // Unobserved targets become zero after their uniform initial prior.
                alpha.copy_from_slice(&next_alpha);
            }
            changes
        };

        let stop_em = chcount == 0 && iter > options.min_rounds;
        if final_round {
            break;
        }
        if stop_em {
            final_round = true;
            if let Some(groups) = &mut components {
                components::reset_fixed(groups);
            }
            if let Some(groups) = &mut scalar_components {
                components::reset_fixed(groups);
            }
            for val in &mut alpha {
                if *val < ALPHA_LIMIT / 10.0 {
                    *val = 0.0;
                }
            }
        }
        iter += 1;
    }

    if crate::timing::is_enabled() {
        eprintln!("[quant-em] iterations={iterations}");
    }

    if compact || serial_dense {
        let mut full_alpha = vec![0.0; num_trans];
        for (&original, &value) in em_work.active_transcripts.iter().zip(&alpha) {
            full_alpha[original] = value;
        }
        alpha = full_alpha;
    }

    let mut rates = Vec::with_capacity(num_trans);
    let mut rate_sum = 0.0;
    for i in 0..num_trans {
        let rate = if eff_lens[i] > 0.0 {
            alpha[i] / eff_lens[i]
        } else {
            0.0
        };
        rates.push(rate);
        rate_sum += rate;
    }

    let mut tpm = vec![0.0; num_trans];
    if rate_sum > 0.0 {
        for i in 0..num_trans {
            tpm[i] = rates[i] / rate_sum * 1_000_000.0;
        }
    }

    Ok(QuantResult {
        est_counts: alpha,
        tpm,
        eff_lengths: eff_lens,
        post_bias,
    })
}

fn counts_to_tpm(est_counts: &[f64], eff_lens: &[f64]) -> Vec<f64> {
    let mut tpm = vec![0.0; est_counts.len()];
    let mut total_mass = 0.0;
    for i in 0..est_counts.len() {
        let eff_len = eff_lens.get(i).copied().unwrap_or(0.0);
        if eff_len > 0.0 {
            tpm[i] = est_counts[i] / eff_len;
            total_mass += tpm[i];
        }
    }
    if total_mass > 0.0 {
        for val in &mut tpm {
            *val = *val / total_mass * 1_000_000.0;
        }
    }
    tpm
}

#[cfg(test)]
mod incoming_id_tests {
    use super::*;

    #[test]
    fn compact_and_wide_incoming_ids_preserve_skip_and_sum_order() {
        let ids32 = [0u32, 1, 2];
        let ids64 = [0usize, 1, 2];
        let mut weights = [2.0, 4.0, 6.0];
        let mut scales = EcScales {
            values: vec![0.5, 0.0, 2.0],
            present: vec![1, 0, 1],
        };
        assert_eq!(sum_scaled_edges(&ids32, &weights, &scales, 1.0, 2.0), 27.0);
        assert_eq!(
            sum_scaled_edges(&ids32, &weights, &scales, 1.0, 2.0),
            sum_scaled_edges(&ids64, &weights, &scales, 1.0, 2.0)
        );
        weights[2] = f64::INFINITY;
        scales.present[2] = 0;
        assert_eq!(sum_scaled_edges(&ids32, &weights, &scales, 1.0, 0.0), 1.0);
        assert_eq!(sum_scaled_edges(&ids64, &weights, &scales, 1.0, 0.0), 1.0);
        scales.present[2] = 1;
        scales.values[2] = 1.0;
        assert!(sum_scaled_edges(&ids32, &weights, &scales, 1.0, 0.0).is_nan());
        assert!(sum_scaled_edges(&ids64, &weights, &scales, 1.0, 0.0).is_nan());
        assert!(matches!(
            IncomingEcIds::new(3, 3),
            IncomingEcIds::Compact(_)
        ));
        if let Some(wide_ecs) = (u32::MAX as usize).checked_add(1) {
            assert!(matches!(
                IncomingEcIds::new(3, wide_ecs),
                IncomingEcIds::Wide(_)
            ));
        }
    }
}

#[cfg(test)]
mod scale_presence_tests {
    use super::*;

    #[test]
    fn split_scales_preserve_skip_and_exceptional_values() {
        let options = [
            None,
            Some(0.0),
            Some(2.0),
            Some(f64::INFINITY),
            Some(f64::NAN),
        ];
        for alpha in [0.0, 1.0, f64::INFINITY, f64::NAN] {
            for selected in 0..options.len() {
                let mut scales = EcScales::new(options.len());
                for (i, value) in options.iter().enumerate() {
                    scales.present[i] = u8::from(value.is_some());
                    scales.values[i] = value.unwrap_or(f64::NAN);
                }
                let mut expected = 7.0;
                if let Some(scale) = options[selected] {
                    expected += 0.5 * alpha * scale;
                }
                let result = sum_scaled_edges(&[selected as u32], &[0.5], &scales, 7.0, alpha);
                assert!(
                    result.to_bits() == expected.to_bits() || result.is_nan() && expected.is_nan()
                );
            }
        }
    }
}
