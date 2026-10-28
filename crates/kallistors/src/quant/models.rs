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
    // A short transcript can only contain fragments up to its own length.
    // Store each prefix's conditional mean, rather than applying one global mean.
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
    update_eff_lens_impl::<false>(
        means,
        bias_counts,
        transcript_seqs,
        alpha,
        eff_lens,
        strand_specific,
        None,
    )
}

pub(super) fn update_eff_lens_with_pool(
    means: &[f64],
    bias_counts: &BiasCounts,
    transcript_seqs: &[Vec<u8>],
    alpha: &[f64],
    eff_lens: &[f64],
    strand_specific: Option<StrandSpecific>,
    pool: Option<&rayon::ThreadPool>,
) -> (Vec<f64>, Vec<f64>) {
    match pool {
        Some(pool) => update_eff_lens_impl::<true>(
            means,
            bias_counts,
            transcript_seqs,
            alpha,
            eff_lens,
            strand_specific,
            Some(pool),
        ),
        None => update_eff_lens(
            means,
            bias_counts,
            transcript_seqs,
            alpha,
            eff_lens,
            strand_specific,
        ),
    }
}

// Two scheduling specializations share the same ordered floating-point code.
// Outlining prevents parallel scheduler machinery from expanding scalar callers.
#[inline(never)]
fn update_eff_lens_impl<const PARALLEL: bool>(
    means: &[f64],
    bias_counts: &BiasCounts,
    transcript_seqs: &[Vec<u8>],
    alpha: &[f64],
    eff_lens: &[f64],
    strand_specific: Option<StrandSpecific>,
    pool: Option<&rayon::ThreadPool>,
) -> (Vec<f64>, Vec<f64>) {
    let num6 = 4096usize;
    let bias_data_norm: f64 = bias_counts.counts.iter().map(|&v| v as f64).sum();
    if bias_data_norm == 0.0 {
        return (eff_lens.to_vec(), vec![1.0; num6]);
    }

    let mut dbias5 = vec![0.0; num6];
    let mut bias_alpha_norm = 0.0;

    // These shared histogram additions stay in transcript order. Parallelism is
    // reserved for the independent per-transcript rescaling after normalization.
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
    let count = transcript_seqs
        .len()
        .min(means.len())
        .min(alpha.len())
        .min(eff_lens.len());
    let rescale = |i: usize| {
        let mean = means[i];
        let seq = &transcript_seqs[i];
        if seq.len() < 6 || (seq.len() as f64) < mean || alpha[i] < MIN_ALPHA {
            return eff_lens[i];
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

        if efflen > mean { efflen } else { eff_lens[i] }
    };
    if PARALLEL {
        let pool = pool.expect("parallel rescale requires its existing EM pool");
        pool.install(|| {
            bias_lens[..count]
                .par_iter_mut()
                .enumerate()
                .for_each(|(i, output)| *output = rescale(i))
        });
    } else {
        for (i, output) in bias_lens[..count].iter_mut().enumerate() {
            *output = rescale(i);
        }
    }

    (bias_lens, dbias5)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Direct six-base windows independently specify the ordered model. No
    // rolling encoder, cached ratio, parallel helper or production walker.
    fn windows(seq: &[u8], start: usize, end: usize, reverse: bool) -> Vec<usize> {
        let mut out = Vec::new();
        for position in start..end {
            let mut hex = 0;
            for (i, &base) in seq[position..position + 6].iter().enumerate() {
                let code = match base {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => return out,
                };
                if reverse {
                    hex |= (3 - code) << (2 * i);
                } else {
                    hex = (hex << 2) | code;
                }
            }
            out.push(hex);
        }
        out
    }

    #[test]
    fn bias_model_matches_ordered_direct_windows_and_partial_inputs() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap();
        let seqs = vec![
            b"acGTACGTACGTACGTacGTACGTA".to_vec(),
            b"ACGTNACGTACGTACGTACGTACG".to_vec(),
            b"TTTCCCAAAGGGTTTCCCAAAGGG".to_vec(),
            b"ACG".to_vec(),
            b"ACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        for strand in [None, Some(StrandSpecific::FR), Some(StrandSpecific::RF)] {
            for exceptional in [false, true] {
                for available in [0, 3, 5] {
                    let means = [8.25, 7.0, 8.0, 8.0, f64::NAN];
                    let alpha = [
                        1.5,
                        0.0,
                        if exceptional { f64::INFINITY } else { 4.0 },
                        7.0,
                        0.5,
                    ];
                    let eff = [18.5, 20.0, 15.0, 3.0, 16.0];
                    for width in [0, 1, 4096, 4103] {
                        let mut counts = BiasCounts::new();
                        counts.counts = (0..width)
                            .map(|i| if i % 13 == 0 { 0 } else { (i % 31) as u32 })
                            .collect();
                        // Short arrays are permitted if no missing positive
                        // denominator is visited; the empty histogram fast
                        // return is a distinct partial-input contract.
                        if width < 4096 {
                            counts.counts.fill(0);
                        }
                        let n = available;
                        let (actual_lens, actual_hist) = update_eff_lens_with_pool(
                            &means[..n],
                            &counts,
                            &seqs,
                            &alpha,
                            &eff,
                            strand,
                            None,
                        );
                        let (parallel_lens, parallel_hist) = update_eff_lens_with_pool(
                            &means[..n],
                            &counts,
                            &seqs,
                            &alpha,
                            &eff,
                            strand,
                            Some(&pool),
                        );
                        assert_eq!(
                            actual_lens.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                            parallel_lens
                                .iter()
                                .map(|x| x.to_bits())
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            actual_hist.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                            parallel_hist
                                .iter()
                                .map(|x| x.to_bits())
                                .collect::<Vec<_>>()
                        );
                        let norm: f64 = counts.counts.iter().map(|&x| x as f64).sum();
                        let mut histogram = vec![0.0; 4096];
                        let mut visits = vec![Vec::new(); seqs.len()];
                        for i in 0..n {
                            let seq = &seqs[i];
                            if alpha[i] < MIN_ALPHA
                                || seq.len() < 6
                                || (seq.len() as f64) < means[i]
                            {
                                continue;
                            }
                            if strand != Some(StrandSpecific::RF) {
                                visits[i].extend(windows(
                                    seq,
                                    0,
                                    (seq.len() as f64 - means[i] - 6.0).max(0.0) as usize,
                                    false,
                                ));
                            }
                            if strand != Some(StrandSpecific::FR) {
                                let start = (means[i] - 6.0).max(0.0) as usize;
                                if start + 6 <= seq.len() {
                                    visits[i].extend(windows(seq, start, seq.len() - 6, true));
                                }
                            }
                            let contrib = if strand.is_some() {
                                alpha[i] / eff[i]
                            } else {
                                0.5 * alpha[i] / eff[i]
                            };
                            for &hex in &visits[i] {
                                histogram[hex] += contrib;
                            }
                        }
                        let mut expected_lens = eff.to_vec();
                        let mut expected_hist = vec![1.0; 4096];
                        if norm != 0.0 {
                            let mut hnorm = 0.0;
                            for &v in &histogram {
                                hnorm += v;
                            }
                            expected_hist = histogram.clone();
                            if hnorm != 0.0 {
                                expected_lens = vec![0.0; seqs.len()];
                                for i in 0..n {
                                    if seqs[i].len() < 6
                                        || (seqs[i].len() as f64) < means[i]
                                        || alpha[i] < MIN_ALPHA
                                    {
                                        expected_lens[i] = eff[i];
                                        continue;
                                    }
                                    let mut length = 0.0;
                                    for &hex in &visits[i] {
                                        if histogram[hex] > 0.0 {
                                            length += counts.counts[hex] as f64 / histogram[hex];
                                        }
                                    }
                                    length *= if strand.is_some() {
                                        hnorm / norm
                                    } else {
                                        0.5 * hnorm / norm
                                    };
                                    expected_lens[i] =
                                        if length > means[i] { length } else { eff[i] };
                                }
                            }
                        }
                        assert_eq!(
                            actual_lens.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                            expected_lens
                                .iter()
                                .map(|x| x.to_bits())
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            actual_hist.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                            expected_hist
                                .iter()
                                .map(|x| x.to_bits())
                                .collect::<Vec<_>>()
                        );
                    }
                }
            }
        }
    }
}
