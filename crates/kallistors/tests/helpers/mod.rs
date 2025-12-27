use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

pub struct SyntheticIndex {
    _tempdir: TempDir,
    pub index_path: PathBuf,
    pub transcripts: Vec<Vec<u8>>,
}

pub fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

pub fn synthetic_transcripts() -> Vec<(String, Vec<u8>)> {
    let mut transcripts = Vec::new();
    for i in 0..4 {
        let name = format!("tx{i}");
        let seq = make_seq(i as usize, 240);
        transcripts.push((name, seq));
    }
    transcripts
}

pub fn build_synthetic_index() -> Option<SyntheticIndex> {
    if !kallisto_available() {
        return None;
    }
    let tempdir = TempDir::new().ok()?;
    let fasta_path = tempdir.path().join("synthetic.fa");
    let index_path = tempdir.path().join("synthetic.idx");
    let transcripts = synthetic_transcripts();
    write_fasta(&fasta_path, &transcripts).ok()?;

    let status = Command::new("kallisto")
        .arg("index")
        .arg("-i")
        .arg(&index_path)
        .arg(&fasta_path)
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }

    Some(SyntheticIndex {
        _tempdir: tempdir,
        index_path,
        transcripts: transcripts.into_iter().map(|(_, seq)| seq).collect(),
    })
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
