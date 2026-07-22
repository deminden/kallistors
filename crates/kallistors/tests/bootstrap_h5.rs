use std::ffi::OsStr;
use std::fs;
use std::path::Path;
use std::process::{Command, ExitStatus};

use tempfile::TempDir;

#[test]
fn quant_writes_bootstrap_h5() {
    let tempdir = TempDir::new().expect("tempdir");
    let fasta = tempdir.path().join("transcripts.fa");
    let index = tempdir.path().join("index.idx");
    let reads1 = tempdir.path().join("reads_1.fq");
    let reads2 = tempdir.path().join("reads_2.fq");
    let out = tempdir.path().join("quant_out");

    let transcripts = synthetic_transcripts();
    write_fasta(&fasta, &transcripts).expect("write fasta");
    write_fastq_pairs(&reads1, &reads2, &transcripts, 50, 90, 24).expect("write reads");

    let status = run_kallistors(&[
        OsStr::new("index"),
        OsStr::new("-i"),
        index.as_os_str(),
        fasta.as_os_str(),
    ])
    .expect("kallistors index");
    assert!(status.success(), "kallistors index failed");

    let status = run_kallistors(&[
        OsStr::new("quant"),
        OsStr::new("-i"),
        index.as_os_str(),
        OsStr::new("-o"),
        out.as_os_str(),
        OsStr::new("-b"),
        OsStr::new("3"),
        OsStr::new("--seed"),
        OsStr::new("11"),
        reads1.as_os_str(),
        reads2.as_os_str(),
    ])
    .expect("kallistors quant");
    assert!(status.success(), "kallistors quant failed");

    let run_info = fs::read_to_string(out.join("run_info.json")).expect("run_info");
    assert_eq!(json_u64(&run_info, "n_bootstraps"), Some(3));

    let h5 = hdf5_pure::File::open(out.join("abundance.h5")).expect("open h5");
    let ids = h5
        .dataset("aux/ids")
        .expect("ids dataset")
        .read_string()
        .expect("read ids");
    let expected_ids = transcripts
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    assert_eq!(ids, expected_ids);

    let num_bootstrap = h5
        .dataset("aux/num_bootstrap")
        .expect("num_bootstrap dataset")
        .read_i32()
        .expect("read num_bootstrap");
    assert_eq!(num_bootstrap, vec![3]);

    let est_counts = h5
        .dataset("est_counts")
        .expect("est_counts dataset")
        .read_f64()
        .expect("read est_counts");
    assert_eq!(est_counts.len(), transcripts.len());
    assert!(est_counts.iter().any(|value| *value > 0.0));

    let bias_observed = h5.dataset("aux/bias_observed").expect("bias_observed");
    assert_eq!(
        bias_observed.shape().expect("bias_observed shape"),
        vec![4096]
    );
    assert_eq!(
        bias_observed.read_i32().expect("read bias_observed"),
        vec![1; 4096]
    );

    let bias_normalized = h5.dataset("aux/bias_normalized").expect("bias_normalized");
    assert_eq!(
        bias_normalized.shape().expect("bias_normalized shape"),
        vec![4096]
    );
    assert!(
        bias_normalized
            .read_f64()
            .expect("read bias_normalized")
            .iter()
            .all(|value| *value == 1.0)
    );

    for idx in 0..3 {
        let counts = h5
            .dataset(&format!("bootstrap/bs{idx}"))
            .expect("bootstrap dataset")
            .read_f64()
            .expect("read bootstrap counts");
        assert_eq!(counts.len(), transcripts.len());
    }
}

#[test]
fn quant_plaintext_bootstrap_skips_h5() {
    let tempdir = TempDir::new().expect("tempdir");
    let fasta = tempdir.path().join("transcripts.fa");
    let index = tempdir.path().join("index.idx");
    let reads1 = tempdir.path().join("reads_1.fq");
    let reads2 = tempdir.path().join("reads_2.fq");
    let out = tempdir.path().join("quant_out");

    let transcripts = synthetic_transcripts();
    write_fasta(&fasta, &transcripts).expect("write fasta");
    write_fastq_pairs(&reads1, &reads2, &transcripts, 50, 90, 12).expect("write reads");

    let status = run_kallistors(&[
        OsStr::new("index"),
        OsStr::new("-i"),
        index.as_os_str(),
        fasta.as_os_str(),
    ])
    .expect("kallistors index");
    assert!(status.success(), "kallistors index failed");

    let status = run_kallistors(&[
        OsStr::new("quant"),
        OsStr::new("-i"),
        index.as_os_str(),
        OsStr::new("-o"),
        out.as_os_str(),
        OsStr::new("-b"),
        OsStr::new("2"),
        OsStr::new("--plaintext"),
        reads1.as_os_str(),
        reads2.as_os_str(),
    ])
    .expect("kallistors quant");
    assert!(status.success(), "kallistors quant failed");

    assert!(out.join("bs_abundance_0.tsv").exists());
    assert!(out.join("bs_abundance_1.tsv").exists());
    assert!(!out.join("abundance.h5").exists());
}

fn run_kallistors(args: &[&OsStr]) -> std::io::Result<ExitStatus> {
    if let Ok(bin) = std::env::var("CARGO_BIN_EXE_kallistors") {
        Command::new(bin).args(args).status()
    } else {
        Command::new("cargo")
            .arg("run")
            .arg("-p")
            .arg("kallistors")
            .arg("--quiet")
            .arg("--")
            .args(args)
            .status()
    }
}

fn synthetic_transcripts() -> Vec<(String, Vec<u8>)> {
    let mut transcripts = Vec::new();
    for i in 0..4 {
        transcripts.push((format!("tx{i}"), make_seq(i, 240)));
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
    fs::write(path, out)
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
    while out.len() < count {
        let seq = &transcripts[idx % transcripts.len()].1;
        idx += 1;
        if seq.len() < frag_len {
            continue;
        }
        let start = (idx * 7) % (seq.len() - frag_len + 1);
        let frag = &seq[start..start + frag_len];
        out.push((
            frag[..read_len].to_vec(),
            revcomp(&frag[frag_len - read_len..]),
        ));
    }
    out
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

fn make_seq(seed: usize, len: usize) -> Vec<u8> {
    let bases = *b"ACGT";
    let mut seq = Vec::with_capacity(len);
    let mut state = 0x9e37_79b9_7f4a_7c15u64 ^ (seed as u64);
    for _ in 0..len {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seq.push(bases[(state >> 32) as usize & 3]);
    }
    seq
}

fn json_u64(input: &str, key: &str) -> Option<u64> {
    let needle = format!("\"{}\"", key);
    let idx = input.find(&needle)?;
    let after = &input[idx + needle.len()..];
    let colon = after.find(':')?;
    let rest = after[colon + 1..].trim_start();
    let mut digits = String::new();
    for ch in rest.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
        } else if !digits.is_empty() {
            break;
        }
    }
    digits.parse().ok()
}
