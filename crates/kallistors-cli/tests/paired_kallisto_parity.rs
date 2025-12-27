use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

#[test]
fn paired_outputs_match_kallisto() {
    if !kallisto_available() {
        return;
    }
    let tempdir = TempDir::new().expect("tempdir");
    let fasta = tempdir.path().join("transcripts.fa");
    let index = tempdir.path().join("index.idx");
    let reads1 = tempdir.path().join("reads_1.fq");
    let reads2 = tempdir.path().join("reads_2.fq");
    let kallisto_out = tempdir.path().join("kallisto_out");
    let ours_out = tempdir.path().join("kallistors_out");

    let transcripts = synthetic_transcripts();
    write_fasta(&fasta, &transcripts).expect("write fasta");
    write_fastq_pairs(&reads1, &reads2, &transcripts, 50, 90, 16).expect("write fastq");

    let status = Command::new("kallisto")
        .arg("index")
        .arg("-i")
        .arg(&index)
        .arg(&fasta)
        .status()
        .expect("kallisto index");
    assert!(status.success(), "kallisto index failed");

    let status = Command::new("kallisto")
        .arg("quant")
        .arg("-i")
        .arg(&index)
        .arg("-o")
        .arg(&kallisto_out)
        .arg(&reads1)
        .arg(&reads2)
        .status()
        .expect("kallisto quant");
    assert!(status.success(), "kallisto quant failed");

    let status =
        run_kallistors_quant(&index, &ours_out, &reads1, &reads2).expect("kallistors quant");
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
    assert_eq!(
        json_u64(&kallisto_info, "index_version"),
        json_u64(&ours_info, "index_version")
    );
    assert_eq!(
        json_u64(&kallisto_info, "k-mer length"),
        json_u64(&ours_info, "k-mer length")
    );

    let p_kallisto = json_f64(&kallisto_info, "p_pseudoaligned").unwrap_or(0.0);
    let p_ours = json_f64(&ours_info, "p_pseudoaligned").unwrap_or(0.0);
    assert!((p_kallisto - p_ours).abs() < 0.1);

    let pu_kallisto = json_f64(&kallisto_info, "p_unique").unwrap_or(0.0);
    let pu_ours = json_f64(&ours_info, "p_unique").unwrap_or(0.0);
    assert!((pu_kallisto - pu_ours).abs() < 0.1);

    let kallisto_abund = read_abundance(kallisto_out.join("abundance.tsv"));
    let ours_abund = read_abundance(ours_out.join("abundance.tsv"));
    assert_eq!(kallisto_abund.len(), ours_abund.len());
    for (target, (len_k, eff_k, est_k, tpm_k)) in kallisto_abund {
        let (len_o, eff_o, est_o, tpm_o) = ours_abund
            .get(&target)
            .copied()
            .unwrap_or((0, 0.0, 0.0, 0.0));
        assert_eq!(len_k, len_o, "length mismatch for {target}");
        assert!(
            approx_eq(eff_k, eff_o, 1e-2),
            "eff_length mismatch for {target}"
        );
        assert!(
            approx_eq(est_k, est_o, 1e-2),
            "est_counts mismatch for {target}"
        );
        assert!(approx_eq(tpm_k, tpm_o, 1e-2), "tpm mismatch for {target}");
    }

    let kallisto_flens = fs::read_to_string(kallisto_out.join("flens.txt")).ok();
    let ours_flens = fs::read_to_string(ours_out.join("flens.txt")).ok();
    if let Some(kallisto_flens) = kallisto_flens {
        let ours_flens = ours_flens.expect("expected flens.txt from kallistors");
        assert_eq!(kallisto_flens.trim(), ours_flens.trim());
    } else {
        assert!(ours_flens.is_none(), "unexpected flens.txt from kallistors");
    }
}

fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn synthetic_transcripts() -> Vec<(String, Vec<u8>)> {
    let mut transcripts = Vec::new();
    for i in 0..4 {
        let name = format!("tx{i}");
        let seq = make_seq(i as usize, 240);
        transcripts.push((name, seq));
    }
    transcripts
}

fn write_fasta(path: &Path, transcripts: &[(String, Vec<u8>)]) -> std::io::Result<()> {
    let mut out = String::new();
    for (name, seq) in transcripts {
        out.push('>');
        out.push_str(name);
        out.push('\n');
        out.push_str(std::str::from_utf8(seq).unwrap_or(""));
        out.push('\n');
    }
    std::fs::write(path, out)
}

fn write_fastq_pairs(
    reads1: &Path,
    reads2: &Path,
    transcripts: &[(String, Vec<u8>)],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> std::io::Result<()> {
    let mut fq1 = Vec::new();
    let mut fq2 = Vec::new();
    let pairs = synthesize_paired_reads(transcripts, read_len, frag_len, count);
    for (i, (r1, r2)) in pairs.iter().enumerate() {
        fq1.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        fq1.extend_from_slice(r1);
        fq1.extend_from_slice(b"\n+\n");
        fq1.extend_from_slice(vec![b'I'; r1.len()].as_slice());
        fq1.extend_from_slice(b"\n");

        fq2.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        fq2.extend_from_slice(r2);
        fq2.extend_from_slice(b"\n+\n");
        fq2.extend_from_slice(vec![b'I'; r2.len()].as_slice());
        fq2.extend_from_slice(b"\n");
    }
    fs::write(reads1, fq1)?;
    fs::write(reads2, fq2)?;
    Ok(())
}

fn synthesize_paired_reads(
    transcripts: &[(String, Vec<u8>)],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut out = Vec::with_capacity(count);
    let mut idx = 0usize;
    let mut attempts = 0usize;
    while out.len() < count && attempts < count * 10 {
        let seq = &transcripts[idx % transcripts.len()].1;
        idx += 1;
        attempts += 1;
        if seq.len() < frag_len {
            continue;
        }
        let start = (idx * 7) % (seq.len() - frag_len + 1);
        let frag = &seq[start..start + frag_len];
        let read1 = frag[..read_len].to_vec();
        let read2 = revcomp(&frag[frag_len - read_len..]);
        out.push((read1, read2));
    }
    out
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        });
    }
    out
}

fn make_seq(seed: usize, len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(len);
    let mut state = 0x9e37_79b9_7f4a_7c15u64 ^ (seed as u64);
    for _ in 0..len {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let idx = (state >> 32) as usize & 3;
        seq.push(bases[idx]);
    }
    seq
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

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

fn run_kallistors_quant(
    index: &Path,
    out_dir: &Path,
    reads1: &Path,
    reads2: &Path,
) -> std::io::Result<std::process::ExitStatus> {
    if let Ok(bin) = std::env::var("CARGO_BIN_EXE_kallistors-cli") {
        Command::new(bin)
            .arg("quant")
            .arg("-i")
            .arg(index)
            .arg("-o")
            .arg(out_dir)
            .arg(reads1)
            .arg(reads2)
            .status()
    } else {
        Command::new("cargo")
            .arg("run")
            .arg("-p")
            .arg("kallistors-cli")
            .arg("--quiet")
            .arg("--")
            .arg("quant")
            .arg("-i")
            .arg(index)
            .arg("-o")
            .arg(out_dir)
            .arg(reads1)
            .arg(reads2)
            .status()
    }
}
