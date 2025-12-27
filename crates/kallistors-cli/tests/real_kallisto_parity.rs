use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

const INDEX_PATH: &str = "data/Human_kallisto_index";
const READS_PATH: &str = "data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz";

// Tight tolerances: allow minor float/EM differences while enforcing strong parity.
const P_PSEUDOALIGNED_TOL: f64 = 5e-5;
const P_UNIQUE_TOL: f64 = 5e-5;
const MAX_EFF_LENGTH_DIFF: f64 = 5e-3;

const REL_P50_MAX: f64 = 5e-5;
const REL_P90_MAX: f64 = 2e-4;
const REL_P99_MAX: f64 = 2e-3;
const REL_MAX_MAX: f64 = 2e-2;

#[test]
fn real_data_parity_with_kallisto() {
    if !kallisto_available() {
        return;
    }
    let index = Path::new(INDEX_PATH);
    let reads = Path::new(READS_PATH);
    if !index.exists() || !reads.exists() {
        return;
    }

    let tempdir = TempDir::new().expect("tempdir");
    let kallisto_out = tempdir.path().join("kallisto_out");
    let ours_out = tempdir.path().join("kallistors_out");

    let status = Command::new("kallisto")
        .arg("quant")
        .arg("--single")
        .arg("-l")
        .arg("200")
        .arg("-s")
        .arg("20")
        .arg("-i")
        .arg(index)
        .arg("-o")
        .arg(&kallisto_out)
        .arg("-t")
        .arg("1")
        .arg(reads)
        .status()
        .expect("kallisto quant");
    assert!(status.success(), "kallisto quant failed");

    let status = run_kallistors_quant(index, &ours_out, reads).expect("kallistors quant");
    assert!(status.success(), "kallistors quant failed");

    let kallisto_info =
        fs::read_to_string(kallisto_out.join("run_info.json")).expect("read kallisto run_info");
    let ours_info =
        fs::read_to_string(ours_out.join("run_info.json")).expect("read kallistors run_info");

    assert_eq!(
        json_u64(&kallisto_info, "n_targets"),
        json_u64(&ours_info, "n_targets")
    );
    assert_eq!(
        json_u64(&kallisto_info, "n_processed"),
        json_u64(&ours_info, "n_processed")
    );
    assert_eq!(
        json_u64(&kallisto_info, "n_pseudoaligned"),
        json_u64(&ours_info, "n_pseudoaligned")
    );
    assert_eq!(
        json_u64(&kallisto_info, "n_unique"),
        json_u64(&ours_info, "n_unique")
    );

    let p_kallisto = json_f64(&kallisto_info, "p_pseudoaligned").unwrap_or(0.0);
    let p_ours = json_f64(&ours_info, "p_pseudoaligned").unwrap_or(0.0);
    assert!(
        (p_kallisto - p_ours).abs() <= P_PSEUDOALIGNED_TOL,
        "p_pseudoaligned diff too large: {p_kallisto} vs {p_ours}"
    );

    let pu_kallisto = json_f64(&kallisto_info, "p_unique").unwrap_or(0.0);
    let pu_ours = json_f64(&ours_info, "p_unique").unwrap_or(0.0);
    assert!(
        (pu_kallisto - pu_ours).abs() <= P_UNIQUE_TOL,
        "p_unique diff too large: {pu_kallisto} vs {pu_ours}"
    );

    let kallisto_abund = read_abundance(kallisto_out.join("abundance.tsv"));
    let ours_abund = read_abundance(ours_out.join("abundance.tsv"));
    assert_eq!(kallisto_abund.len(), ours_abund.len());

    let mut est_errors = Vec::with_capacity(kallisto_abund.len());
    let mut tpm_errors = Vec::with_capacity(kallisto_abund.len());
    let mut max_eff_length_diff: f64 = 0.0;

    for (target, (len_k, eff_k, est_k, tpm_k)) in kallisto_abund {
        let (len_o, eff_o, est_o, tpm_o) = ours_abund
            .get(&target)
            .copied()
            .unwrap_or((0, 0.0, 0.0, 0.0));
        assert_eq!(len_k, len_o, "length mismatch for {target}");
        max_eff_length_diff = max_eff_length_diff.max((eff_k - eff_o).abs());
        est_errors.push(rel_error(est_k, est_o));
        tpm_errors.push(rel_error(tpm_k, tpm_o));
    }

    est_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    tpm_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Use percentiles to avoid a single outlier from dominating parity checks.
    let est_p50 = percentile(&est_errors, 50.0);
    let est_p90 = percentile(&est_errors, 90.0);
    let est_p99 = percentile(&est_errors, 99.0);
    let est_max = est_errors.last().copied().unwrap_or(0.0);

    let tpm_p50 = percentile(&tpm_errors, 50.0);
    let tpm_p90 = percentile(&tpm_errors, 90.0);
    let tpm_p99 = percentile(&tpm_errors, 99.0);
    let tpm_max = tpm_errors.last().copied().unwrap_or(0.0);

    assert!(
        max_eff_length_diff <= MAX_EFF_LENGTH_DIFF,
        "max eff_length diff too large: {max_eff_length_diff}"
    );

    assert!(
        est_p50 <= REL_P50_MAX
            && est_p90 <= REL_P90_MAX
            && est_p99 <= REL_P99_MAX
            && est_max <= REL_MAX_MAX,
        "est_counts errors too large: p50={est_p50} p90={est_p90} p99={est_p99} max={est_max}"
    );
    assert!(
        tpm_p50 <= REL_P50_MAX
            && tpm_p90 <= REL_P90_MAX
            && tpm_p99 <= REL_P99_MAX
            && tpm_max <= REL_MAX_MAX,
        "tpm errors too large: p50={tpm_p50} p90={tpm_p90} p99={tpm_p99} max={tpm_max}"
    );
}

fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn run_kallistors_quant(
    index: &Path,
    out_dir: &Path,
    reads: &Path,
) -> std::io::Result<std::process::ExitStatus> {
    if let Ok(bin) = std::env::var("CARGO_BIN_EXE_kallistors-cli") {
        Command::new(bin)
            .arg("quant")
            .arg("--single")
            .arg("-l")
            .arg("200")
            .arg("-s")
            .arg("20")
            .arg("-i")
            .arg(index)
            .arg("-o")
            .arg(out_dir)
            .arg("-t")
            .arg("1")
            .arg(reads)
            .status()
    } else {
        Command::new("cargo")
            .arg("run")
            .arg("-p")
            .arg("kallistors-cli")
            .arg("--quiet")
            .arg("--")
            .arg("quant")
            .arg("--single")
            .arg("-l")
            .arg("200")
            .arg("-s")
            .arg("20")
            .arg("-i")
            .arg(index)
            .arg("-o")
            .arg(out_dir)
            .arg("-t")
            .arg("1")
            .arg(reads)
            .status()
    }
}

fn json_u64(input: &str, key: &str) -> Option<u64> {
    json_value(input, key)?.parse().ok()
}

fn json_f64(input: &str, key: &str) -> Option<f64> {
    json_value(input, key)?.parse().ok()
}

fn json_value(input: &str, key: &str) -> Option<String> {
    let needle = format!("\"{}\"", key);
    let idx = input.find(&needle)?;
    let after = &input[idx + needle.len()..];
    let colon = after.find(':')?;
    let mut rest = after[colon + 1..].trim_start();
    if rest.starts_with('"') {
        rest = &rest[1..];
        let end = rest.find('"')?;
        return Some(rest[..end].to_string());
    }
    let mut out = String::new();
    for ch in rest.chars() {
        if ch.is_ascii_digit() || ch == '.' || ch == '-' {
            out.push(ch);
        } else if !out.is_empty() {
            break;
        }
    }
    if out.is_empty() { None } else { Some(out) }
}

fn read_abundance(path: PathBuf) -> HashMap<String, (u32, f64, f64, f64)> {
    let mut map = HashMap::new();
    let data = fs::read_to_string(path).expect("read abundance");
    for (line_no, line) in data.lines().enumerate() {
        if line_no == 0 {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }
        let target = parts[0].to_string();
        let length: u32 = parts[1].parse().unwrap_or(0);
        let eff_len: f64 = parts[2].parse().unwrap_or(0.0);
        let est: f64 = parts[3].parse().unwrap_or(0.0);
        let tpm: f64 = parts[4].parse().unwrap_or(0.0);
        map.insert(target, (length, eff_len, est, tpm));
    }
    map
}

fn rel_error(a: f64, b: f64) -> f64 {
    let denom = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() / denom
}

fn percentile(values: &[f64], pct: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let clamped = pct.clamp(0.0, 100.0);
    let idx = ((clamped / 100.0) * (values.len() - 1) as f64).round() as usize;
    values[idx]
}
