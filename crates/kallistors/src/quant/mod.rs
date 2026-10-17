//! Quantification routines (EM over ECs with optional sequence bias).

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

struct EmMultiEc {
    transcript_range: Range<usize>,
    count: f64,
}

struct EmWork {
    singleton_alpha: Vec<f64>,
    active_transcripts: Vec<usize>,
    multi_ecs: Vec<EmMultiEc>,
    multi_transcripts: Vec<usize>,
    multi_weights: Vec<f64>,
}

impl EmWork {
    fn compact_transcript_ids(&mut self) {
        let mut dense_ids = vec![0usize; self.singleton_alpha.len()];
        let mut singleton_alpha = Vec::with_capacity(self.active_transcripts.len());
        for (dense, &original) in self.active_transcripts.iter().enumerate() {
            dense_ids[original] = dense;
            singleton_alpha.push(self.singleton_alpha[original]);
        }
        for tr in &mut self.multi_transcripts {
            *tr = dense_ids[*tr];
        }
        self.singleton_alpha = singleton_alpha;
    }
}

enum IncomingEcIds {
    Compact(Vec<u32>),
    Wide(Vec<usize>),
}

impl IncomingEcIds {
    fn new(edges: usize, ecs: usize) -> Self {
        if ecs <= u32::MAX as usize {
            Self::Compact(vec![0; edges])
        } else {
            Self::Wide(vec![0; edges])
        }
    }

    fn set(&mut self, slot: usize, ec: usize) {
        match self {
            Self::Compact(ids) => ids[slot] = ec as u32,
            Self::Wide(ids) => ids[slot] = ec,
        }
    }
}

trait IncomingEcId: Copy {
    fn as_usize(self) -> usize;
}

impl IncomingEcId for u32 {
    #[inline]
    fn as_usize(self) -> usize {
        self as usize
    }
}

impl IncomingEcId for usize {
    #[inline]
    fn as_usize(self) -> usize {
        self
    }
}

#[inline]
fn weighted_alpha_sum<I: IncomingEcId>(ids: &[I], weights: &[f64], alpha: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (&tr, &weight) in ids.iter().zip(weights) {
        sum += alpha[tr.as_usize()] * weight;
    }
    sum
}

struct EcScales {
    values: Vec<f64>,
    present: Vec<u8>,
}

impl EcScales {
    fn new(len: usize) -> Self {
        Self {
            values: vec![0.0; len],
            present: vec![0; len],
        }
    }
}

#[inline]
fn sum_scaled_edges<I: IncomingEcId>(
    ids: &[I],
    weights: &[f64],
    scales: &EcScales,
    singleton: f64,
    alpha: f64,
) -> f64 {
    let mut sum = singleton;
    for (&ec, &weight) in ids.iter().zip(weights) {
        let ec = ec.as_usize();
        if scales.present[ec] != 0 {
            sum += weight * alpha * scales.values[ec];
        }
    }
    sum
}

// Each transcript's incoming edges retain their original EC/edge order.
// Independent transcript sums can then run concurrently without reordering f64 additions.
struct EmIncoming {
    offsets: Vec<usize>,
    edge_indices: Vec<usize>,
    ec_indices: IncomingEcIds,
    weights: Vec<f64>,
}

impl EmIncoming {
    fn new(work: &EmWork, track_edges: bool) -> Self {
        let mut offsets = vec![0; work.singleton_alpha.len() + 1];
        for &tr in &work.multi_transcripts {
            offsets[tr + 1] += 1;
        }
        for tr in 0..work.singleton_alpha.len() {
            offsets[tr + 1] += offsets[tr];
        }
        let mut cursors = offsets[..work.singleton_alpha.len()].to_vec();
        let edges = work.multi_transcripts.len();
        let mut incoming = Self {
            offsets,
            edge_indices: if track_edges {
                vec![0; edges]
            } else {
                Vec::new()
            },
            ec_indices: IncomingEcIds::new(edges, work.multi_ecs.len()),
            weights: vec![0.0; edges],
        };
        for (ec_id, ec) in work.multi_ecs.iter().enumerate() {
            for edge in ec.transcript_range {
                let tr = work.multi_transcripts[edge];
                let slot = cursors[tr];
                cursors[tr] += 1;
                if track_edges {
                    incoming.edge_indices[slot] = edge;
                }
                incoming.ec_indices.set(slot, ec_id);
                incoming.weights[slot] = work.multi_weights[edge];
            }
        }
        incoming
    }

    #[inline]
    fn transcript_sum(&self, tr: usize, singleton: f64, alpha: f64, scales: &EcScales) -> f64 {
        let range = self.offsets[tr]..self.offsets[tr + 1];
        let weights = &self.weights[range.clone()];
        match &self.ec_indices {
            IncomingEcIds::Compact(ids) => {
                sum_scaled_edges(&ids[range], weights, scales, singleton, alpha)
            }
            IncomingEcIds::Wide(ids) => {
                sum_scaled_edges(&ids[range], weights, scales, singleton, alpha)
            }
        }
    }

    fn refresh_weights(&mut self, work: &EmWork) {
        for (weight, &edge) in self.weights.iter_mut().zip(&self.edge_indices) {
            *weight = work.multi_weights[edge];
        }
    }
}

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

pub fn load_ec_counts(path: &Path) -> Result<EcCountsInput> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let reader = BufReader::new(file);
    let mut classes: Vec<Vec<u32>> = Vec::new();
    let mut counts: Vec<u32> = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(3, '\t');
        let ec_id: usize = parts
            .next()
            .ok_or_else(|| Error::InvalidFormat(format!("missing EC id at line {}", line_no + 1)))?
            .parse()
            .map_err(|_| Error::InvalidFormat(format!("invalid EC id at line {}", line_no + 1)))?;
        let count: u32 = parts
            .next()
            .ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC count at line {}", line_no + 1))
            })?
            .parse()
            .map_err(|_| {
                Error::InvalidFormat(format!("invalid EC count at line {}", line_no + 1))
            })?;
        let list_str = parts
            .next()
            .ok_or_else(|| {
                Error::InvalidFormat(format!("missing EC list at line {}", line_no + 1))
            })?
            .trim();
        let mut ids = Vec::new();
        if !list_str.is_empty() {
            for token in list_str.split(',') {
                if token.is_empty() {
                    continue;
                }
                let val: u32 = token.parse().map_err(|_| {
                    Error::InvalidFormat(format!("invalid transcript id at line {}", line_no + 1))
                })?;
                ids.push(val);
            }
        }
        if ec_id != classes.len() {
            return Err(Error::InvalidFormat(format!(
                "EC id out of order at line {}",
                line_no + 1
            )));
        }
        classes.push(ids);
        counts.push(count);
    }

    Ok(EcCountsInput {
        ec_list: EcList { classes },
        counts,
    })
}

pub fn load_transcript_sequences(fasta: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
    let file = File::open(fasta).map_err(|_| Error::MissingFile(fasta.to_path_buf()))?;
    let reader = BufReader::new(file);
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                seqs.insert(name, std::mem::take(&mut current_seq));
            }
            current_name = Some(rest.trim().to_string());
        } else {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = current_name.take() {
        seqs.insert(name, current_seq);
    }

    let mut out = Vec::with_capacity(names.len());
    for name in names {
        let seq = seqs
            .get(name)
            .cloned()
            .ok_or_else(|| Error::InvalidFormat(format!("missing transcript {}", name)))?;
        out.push(seq);
    }
    Ok(out)
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
    let mut components = if compact {
        let components = em_work.group_components();
        em_work.compact_transcript_ids();
        alpha.truncate(em_work.active_transcripts.len());
        // Compact updates use only their own previous alpha after the denominator barrier.
        next_alpha = Vec::new();
        Some(components)
    } else {
        None
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
    let mut final_round = false;
    let mut iterations = burn_in;
    for iter in burn_in..options.max_iter {
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

        let mut parallel_changes = None;
        if let (Some(pool), Some(incoming), Some(components)) =
            (&em_pool, &incoming, &mut components)
        {
            parallel_changes = Some(match &compact_transcripts {
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
                parallel_changes = Some(pool.install(|| {
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
            for &tr in &em_work.active_transcripts {
                next_alpha[tr] = em_work.singleton_alpha[tr];
            }
            for ec in &em_work.multi_ecs {
                let transcripts = &em_work.multi_transcripts[ec.transcript_range];
                let weights = &em_work.multi_weights[ec.transcript_range];
                let mut denom = 0.0;
                for (&tr, &weight) in transcripts.iter().zip(weights) {
                    denom += alpha[tr] * weight;
                }
                if denom < TOLERANCE {
                    continue;
                }
                let scale = ec.count / denom;
                for (&tr, &weight) in transcripts.iter().zip(weights) {
                    next_alpha[tr] += weight * alpha[tr] * scale;
                }
            }
        }

        let chcount = if let Some(changes) = parallel_changes {
            changes
        } else {
            let mut changes = 0;
            for &i in &em_work.active_transcripts {
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
            for val in &mut alpha {
                if *val < ALPHA_LIMIT / 10.0 {
                    *val = 0.0;
                }
            }
        }
    }

    if crate::timing::is_enabled() {
        eprintln!("[quant-em] iterations={iterations}");
    }

    if compact {
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

    let mut seed_rng = SplitMix64::new(seed);
    let mut out = Vec::with_capacity(bootstrap_samples);
    for _ in 0..bootstrap_samples {
        let mut sample_rng = SplitMix64::new(seed_rng.next_u64());
        let counts = sample_ec_counts(&input.counts, &mut sample_rng);
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

fn em_quantify_fixed_eff_lens(
    ec_classes: &[Vec<u32>],
    counts: &[u32],
    eff_lens: &[f64],
    max_iter: usize,
    min_rounds: usize,
) -> QuantResult {
    let num_trans = eff_lens.len();
    let mut alpha = vec![1.0 / num_trans as f64; num_trans];
    let mut next_alpha = vec![0.0f64; num_trans];
    let em_work = prepare_em_work(ec_classes, counts, eff_lens);

    let mut final_round = false;
    for iter in 0..max_iter {
        for &tr in &em_work.active_transcripts {
            next_alpha[tr] = em_work.singleton_alpha[tr];
        }

        for ec in &em_work.multi_ecs {
            let transcripts = &em_work.multi_transcripts[ec.transcript_range];
            let weights = &em_work.multi_weights[ec.transcript_range];
            let mut denom = 0.0;
            for (&tr, &weight) in transcripts.iter().zip(weights) {
                denom += alpha[tr] * weight;
            }
            if denom < TOLERANCE {
                continue;
            }
            let scale = ec.count / denom;
            for (&tr, &weight) in transcripts.iter().zip(weights) {
                next_alpha[tr] += weight * alpha[tr] * scale;
            }
        }

        let mut chcount = 0;
        for &i in &em_work.active_transcripts {
            if next_alpha[i] > ALPHA_CHANGE_LIMIT {
                let rel = (next_alpha[i] - alpha[i]).abs() / next_alpha[i];
                if rel > ALPHA_CHANGE {
                    chcount += 1;
                }
            }
            alpha[i] = next_alpha[i];
        }
        if iter == 0 {
            alpha.copy_from_slice(&next_alpha);
        }

        let stop_em = chcount == 0 && iter > min_rounds;
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
        }
    }

    let tpm = counts_to_tpm(&alpha, eff_lens);
    QuantResult {
        est_counts: alpha,
        tpm,
        eff_lengths: eff_lens.to_vec(),
        post_bias: None,
    }
}

fn prepare_em_work(ec_classes: &[Vec<u32>], counts: &[u32], eff_lens: &[f64]) -> EmWork {
    let mut singleton_alpha = vec![0.0f64; eff_lens.len()];
    let mut multi_ecs = Vec::new();
    let mut multi_transcripts = Vec::new();
    let mut multi_weights = Vec::new();
    for (ec_id, ec) in ec_classes.iter().enumerate() {
        let count = counts.get(ec_id).copied().unwrap_or(0) as f64;
        if ec.len() == 1 {
            let tr = ec[0] as usize;
            if tr < singleton_alpha.len() {
                singleton_alpha[tr] = count;
            }
        } else if ec.len() > 1 && count != 0.0 {
            let start = multi_transcripts.len();
            for &tr in ec {
                let tr = tr as usize;
                let len = eff_lens.get(tr).copied().unwrap_or(1.0);
                multi_transcripts.push(tr);
                multi_weights.push(count / len);
            }
            multi_ecs.push(EmMultiEc {
                transcript_range: Range {
                    start,
                    end: multi_transcripts.len(),
                },
                count,
            });
        }
    }
    let mut active = vec![false; eff_lens.len()];
    for (tr, &count) in singleton_alpha.iter().enumerate() {
        active[tr] = count != 0.0;
    }
    for &tr in &multi_transcripts {
        active[tr] = true;
    }
    let active_transcripts = active
        .into_iter()
        .enumerate()
        .filter_map(|(tr, observed)| observed.then_some(tr))
        .collect();
    EmWork {
        singleton_alpha,
        active_transcripts,
        multi_ecs,
        multi_transcripts,
        multi_weights,
    }
}

fn sample_ec_counts(counts: &[u32], rng: &mut SplitMix64) -> Vec<u32> {
    let total = counts.iter().map(|&v| u64::from(v)).sum::<u64>();
    let mut sampled = vec![0u32; counts.len()];
    if total == 0 {
        return sampled;
    }

    let mut cumulative = Vec::with_capacity(counts.len());
    let mut running = 0u64;
    for &count in counts {
        running = running.saturating_add(u64::from(count));
        cumulative.push(running);
    }

    for _ in 0..total {
        let draw = rng.gen_range(total);
        let idx = cumulative.partition_point(|&bound| bound <= draw);
        if let Some(slot) = sampled.get_mut(idx) {
            *slot = slot.saturating_add(1);
        }
    }
    sampled
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

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn gen_range(&mut self, upper: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % upper;
        loop {
            let value = self.next_u64();
            if value < zone {
                return value % upper;
            }
        }
    }
}

fn refresh_multi_weights(work: &mut EmWork, eff_lens: &[f64]) {
    for ec in &work.multi_ecs {
        for i in ec.transcript_range {
            let tr = work.multi_transcripts[i];
            let len = eff_lens.get(tr).copied().unwrap_or(1.0);
            work.multi_weights[i] = ec.count / len;
        }
    }
}

fn calc_eff_lens(lengths: &[u32], means: &[f64]) -> Vec<f64> {
    let mut eff = Vec::with_capacity(lengths.len());
    for (i, &len) in lengths.iter().enumerate() {
        let mean = means.get(i).copied().unwrap_or(0.0);
        let mut cur = len as f64 - mean + 1.0;
        if cur < 1.0 {
            cur = len as f64;
        }
        eff.push(cur);
    }
    eff
}

fn mean_frag_lens_by_transcript(lengths: &[u32], mean: f64, sd: f64) -> Vec<f64> {
    const MAX_FRAG_LEN: usize = 1000;
    if sd <= 0.0 {
        return vec![mean; lengths.len()];
    }
    let mean_trunc = trunc_gaussian_fld(0, MAX_FRAG_LEN, mean, sd);
    let marginal = *mean_trunc.last().unwrap_or(&mean);
    lengths
        .iter()
        .map(|&len| {
            let len = len as usize;
            if len >= MAX_FRAG_LEN {
                marginal
            } else {
                mean_trunc[len]
            }
        })
        .collect()
}

fn empirical_frag_lens_by_transcript(lengths: &[u32], histogram: &[u32]) -> Result<Vec<f64>> {
    let mut count = 0u64;
    let mut mass = 0.0;
    let mut truncated_means = Vec::with_capacity(histogram.len());
    for (length, &observations) in histogram.iter().enumerate() {
        count += u64::from(observations);
        mass += length as f64 * f64::from(observations);
        truncated_means.push(if count > 0 { mass / count as f64 } else { 0.0 });
    }
    if count == 0 {
        return Err(Error::InvalidFormat(
            "fragment length histogram has no observations".into(),
        ));
    }
    let marginal_mean = *truncated_means.last().expect("nonempty histogram");
    Ok(lengths
        .iter()
        .map(|&length| {
            truncated_means
                .get(length as usize)
                .copied()
                .unwrap_or(marginal_mean)
        })
        .collect())
}

fn trunc_gaussian_fld(start: usize, stop: usize, mean: f64, sd: f64) -> Vec<f64> {
    let n = stop.saturating_sub(start);
    let mut mean_fl = vec![0.0; n];
    let mut total_mass = 0.0;
    let mut total_density = 0.0;

    for (i, slot) in mean_fl.iter_mut().enumerate() {
        let mut x = (start + i) as f64;
        x = (x - mean) / sd;
        let cur_density = (-0.5 * x * x).exp() / sd;
        total_mass += cur_density * i as f64;
        total_density += cur_density;
        if total_mass > 0.0 {
            *slot = total_mass / total_density;
        }
    }
    mean_fl
}

fn update_eff_lens(
    means: &[f64],
    bias_counts: &BiasCounts,
    transcript_seqs: &[Vec<u8>],
    alpha: &[f64],
    eff_lens: &[f64],
    strand_specific: Option<StrandSpecific>,
) -> (Vec<f64>, Vec<f64>) {
    let num6 = 4096usize;
    let bias_data_norm: f64 = bias_counts.counts.iter().map(|&v| v as f64).sum();
    if bias_data_norm == 0.0 {
        return (eff_lens.to_vec(), vec![1.0; num6]);
    }

    let mut dbias5 = vec![0.0; num6];
    let mut bias_alpha_norm = 0.0;

    for i in 0..transcript_seqs.len() {
        if i >= means.len() || i >= alpha.len() || i >= eff_lens.len() {
            break;
        }
        if alpha[i] < MIN_ALPHA {
            continue;
        }
        let mean = means[i];
        let seq = &transcript_seqs[i];
        if seq.len() < 6 || (seq.len() as f64) < mean {
            continue;
        }
        let mut contrib = 0.5 * alpha[i] / eff_lens[i];
        if strand_specific.is_some() {
            contrib = alpha[i] / eff_lens[i];
        }
        let seqlen = seq.len() as isize;

        if (strand_specific.is_none() || strand_specific == Some(StrandSpecific::FR))
            && let Some(mut hex) = hexamer_to_int(seq, false)
        {
            let fwlimit = (seqlen as f64 - mean - 6.0).max(0.0) as isize;
            for j in 0..fwlimit {
                dbias5[hex] += contrib;
                if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], false) {
                    hex = next;
                } else {
                    break;
                }
            }
        }

        if strand_specific.is_none() || strand_specific == Some(StrandSpecific::RF) {
            let bwlimit = (mean - 6.0).max(0.0) as isize;
            let start = bwlimit as usize;
            if start + 6 <= seq.len()
                && let Some(mut hex) = hexamer_to_int(&seq[start..], true)
            {
                for j in bwlimit..(seqlen - 6) {
                    dbias5[hex] += contrib;
                    if (j as usize) < seq.len() - 6 {
                        if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], true) {
                            hex = next;
                        } else {
                            break;
                        }
                    }
                }
            }
        }
    }

    for val in &dbias5 {
        bias_alpha_norm += *val;
    }
    if bias_alpha_norm == 0.0 {
        return (eff_lens.to_vec(), dbias5);
    }

    let mut bias_lens = vec![0.0; transcript_seqs.len()];
    for i in 0..transcript_seqs.len() {
        if i >= means.len() || i >= alpha.len() || i >= eff_lens.len() {
            break;
        }
        let mean = means[i];
        let seq = &transcript_seqs[i];
        if seq.len() < 6 || (seq.len() as f64) < mean || alpha[i] < MIN_ALPHA {
            bias_lens[i] = eff_lens[i];
            continue;
        }
        let seqlen = seq.len() as isize;
        let mut efflen = 0.0;

        if (strand_specific.is_none() || strand_specific == Some(StrandSpecific::FR))
            && let Some(mut hex) = hexamer_to_int(seq, false)
        {
            let fwlimit = (seqlen as f64 - mean - 6.0).max(0.0) as isize;
            for j in 0..fwlimit {
                let denom = dbias5[hex];
                if denom > 0.0 {
                    efflen += bias_counts.counts[hex] as f64 / denom;
                }
                if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], false) {
                    hex = next;
                } else {
                    break;
                }
            }
        }

        if strand_specific.is_none() || strand_specific == Some(StrandSpecific::RF) {
            let bwlimit = (mean - 6.0).max(0.0) as isize;
            let start = bwlimit as usize;
            if start + 6 <= seq.len()
                && let Some(mut hex) = hexamer_to_int(&seq[start..], true)
            {
                for j in bwlimit..(seqlen - 6) {
                    let denom = dbias5[hex];
                    if denom > 0.0 {
                        efflen += bias_counts.counts[hex] as f64 / denom;
                    }
                    if (j as usize) < seq.len() - 6 {
                        if let Some(next) = update_hexamer(hex, seq[(j + 6) as usize], true) {
                            hex = next;
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        if strand_specific.is_none() {
            efflen *= 0.5 * bias_alpha_norm / bias_data_norm;
        } else {
            efflen *= bias_alpha_norm / bias_data_norm;
        }

        if efflen > mean {
            bias_lens[i] = efflen;
        } else {
            bias_lens[i] = eff_lens[i];
        }
    }

    (bias_lens, dbias5)
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
