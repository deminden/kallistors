//! Bootstrap resampling and quantification with fixed effective lengths.

use super::*;

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

pub(super) fn em_quantify_fixed_eff_lens(
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

pub(super) fn sample_ec_counts(counts: &[u32], rng: &mut SplitMix64) -> Vec<u32> {
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

pub(super) struct SplitMix64 {
    pub(super) state: u64,
}

impl SplitMix64 {
    pub(super) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(super) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    pub(super) fn gen_range(&mut self, upper: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % upper;
        loop {
            let value = self.next_u64();
            if value < zone {
                return value % upper;
            }
        }
    }
}
