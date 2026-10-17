use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus};
use std::thread;
use std::time::{Duration, Instant};

use tempfile::TempDir;

pub struct SyntheticIndex {
    _tempdir: TempDir,
    pub index_path: PathBuf,
    pub transcripts: Vec<Vec<u8>>,
}

pub fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .output()
        .is_ok_and(|output| output.status.success())
}

pub fn synthetic_transcripts() -> Vec<(String, Vec<u8>)> {
    let mut transcripts = Vec::new();
    for i in 0..4 {
        let name = format!("tx{i}");
        let seq = make_seq(i as usize, 240);
        transcripts.push((name, seq));
    }
    // Shared sequence creates ambiguous ECs; the flanks remain transcript-specific.
    let shared = transcripts[0].1[70..170].to_vec();
    transcripts[1].1[70..170].copy_from_slice(&shared);
    transcripts
}

pub fn build_synthetic_index() -> SyntheticIndex {
    build_fixture(false)
}

pub fn build_upstream_synthetic_index() -> SyntheticIndex {
    build_fixture(true)
}

fn build_fixture(upstream: bool) -> SyntheticIndex {
    let tempdir = TempDir::new().expect("synthetic index directory");
    let fasta_path = tempdir.path().join("synthetic.fa");
    let index_path = tempdir.path().join("synthetic.idx");
    let transcripts = synthetic_transcripts();
    write_fasta(&fasta_path, &transcripts).expect("write synthetic FASTA");

    if upstream {
        let output = Command::new("kallisto")
            .args(["index", "-i"])
            .arg(&index_path)
            .arg(&fasta_path)
            .output()
            .expect("run kallisto index");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    } else {
        kallistors::index::build_index(&index_path, &[fasta_path], Default::default())
            .expect("build synthetic index");
    }
    SyntheticIndex {
        _tempdir: tempdir,
        index_path,
        transcripts: transcripts.into_iter().map(|(_, seq)| seq).collect(),
    }
}

pub fn write_fasta(path: &Path, transcripts: &[(String, Vec<u8>)]) -> std::io::Result<()> {
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

fn make_seq(seed: usize, len: usize) -> Vec<u8> {
    let bases = *b"ACGT";
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

pub fn synthesize_paired_reads(
    transcripts: &[Vec<u8>],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut out = Vec::with_capacity(count);
    let mut idx = 0usize;
    let mut attempts = 0usize;
    while out.len() < count && attempts < count * 10 {
        let seq = &transcripts[idx % transcripts.len()];
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

pub fn to_fastq_pairs(pairs: &[(Vec<u8>, Vec<u8>)]) -> (Vec<u8>, Vec<u8>) {
    let mut fq1 = Vec::new();
    let mut fq2 = Vec::new();
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
    (fq1, fq2)
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
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

pub fn json_u64(input: &str, key: &str) -> Option<u64> {
    serde_json::from_str::<serde_json::Value>(input)
        .expect("valid run_info JSON")
        .get(key)?
        .as_u64()
}

pub fn json_f64(input: &str, key: &str) -> Option<f64> {
    serde_json::from_str::<serde_json::Value>(input)
        .expect("valid run_info JSON")
        .get(key)?
        .as_f64()
}

pub fn read_abundance(path: PathBuf) -> HashMap<String, (u32, f64, f64, f64)> {
    let mut map = HashMap::new();
    let data = fs::read_to_string(path).expect("read abundance");
    let mut lines = data.lines();
    assert_eq!(
        lines.next(),
        Some("target_id\tlength\teff_length\test_counts\ttpm")
    );
    for (line_no, line) in lines.enumerate() {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 5, "malformed abundance row {line_no}");
        let target = parts[0].to_string();
        let length: u32 = parts[1].parse().expect("target length");
        let eff_len: f64 = parts[2].parse().expect("effective length");
        let est: f64 = parts[3].parse().expect("estimated count");
        let tpm: f64 = parts[4].parse().expect("TPM");
        assert!(
            [eff_len, est, tpm]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0)
        );
        assert!(
            map.insert(target, (length, eff_len, est, tpm)).is_none(),
            "duplicate transcript ID"
        );
    }
    map
}

pub fn write_fastq_pairs(
    reads1: &Path,
    reads2: &Path,
    transcripts: &[(String, Vec<u8>)],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> std::io::Result<()> {
    let sequences: Vec<_> = transcripts.iter().map(|(_, seq)| seq.clone()).collect();
    let pairs = synthesize_paired_reads(&sequences, read_len, frag_len, count);
    let (left, right) = to_fastq_pairs(&pairs);
    fs::write(reads1, left)?;
    fs::write(reads2, right)
}

pub fn run_kallistors(args: &[&OsStr]) -> std::io::Result<std::process::ExitStatus> {
    Command::new(env!("CARGO_BIN_EXE_kallistors"))
        .args(args)
        .status()
}

pub fn command_status_with_timeout(
    command: &mut Command,
    timeout: Duration,
) -> io::Result<ExitStatus> {
    let description = format!("{command:?}");
    let mut child = command.spawn()?;
    let start = Instant::now();

    loop {
        if let Some(status) = child.try_wait()? {
            return Ok(status);
        }
        if start.elapsed() >= timeout {
            if let Some(status) = child.try_wait()? {
                return Ok(status);
            }
            let kill_result = child.kill();
            let wait_result = child.wait();
            if let Err(error) = kill_result
                && error.kind() != io::ErrorKind::InvalidInput
            {
                return Err(error);
            }
            wait_result?;
            return Err(io::Error::new(
                io::ErrorKind::TimedOut,
                format!("command exceeded {timeout:?} and was terminated: {description}"),
            ));
        }
        thread::sleep(Duration::from_millis(10));
    }
}

pub fn synthesize_reads(transcripts: &[Vec<u8>], read_len: usize, count: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(count);
    let mut idx = 0usize;
    let mut attempts = 0usize;
    while out.len() < count && attempts < count * 10 {
        let seq = &transcripts[idx % transcripts.len()];
        idx += 1;
        attempts += 1;
        if seq.len() < read_len {
            continue;
        }
        let start = (idx * 7) % (seq.len() - read_len + 1);
        out.push(seq[start..start + read_len].to_vec());
    }
    out
}

pub fn to_fastq(reads: &[Vec<u8>]) -> Vec<u8> {
    let mut fq = Vec::new();
    for (i, read) in reads.iter().enumerate() {
        fq.extend_from_slice(format!("@r{}\n", i + 1).as_bytes());
        fq.extend_from_slice(read);
        fq.extend_from_slice(b"\n+\n");
        fq.extend_from_slice(vec![b'I'; read.len()].as_slice());
        fq.extend_from_slice(b"\n");
    }
    fq
}
