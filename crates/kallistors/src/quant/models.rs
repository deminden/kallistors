//! Fragment-length models, effective lengths and sequence bias.

use super::*;

pub(super) fn calc_eff_lens(lengths: &[u32], means: &[f64]) -> Vec<f64> {
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

pub(super) fn mean_frag_lens_by_transcript(lengths: &[u32], mean: f64, sd: f64) -> Vec<f64> {
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

pub(super) fn empirical_frag_lens_by_transcript(
    lengths: &[u32],
    histogram: &[u32],
) -> Result<Vec<f64>> {
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

pub(super) fn trunc_gaussian_fld(start: usize, stop: usize, mean: f64, sd: f64) -> Vec<f64> {
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

pub(super) fn update_eff_lens(
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
