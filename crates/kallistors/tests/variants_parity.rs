use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    PseudoalignOptions, Strand, build_bifrost_index, pseudoalign_paired_bifrost_with_options,
    pseudoalign_single_end_bifrost,
};

mod helpers;

#[test]
fn parity_on_synthetic_variants() {
    let tmp = temp_dir("kallistors-variants");
    let _ = fs::create_dir_all(&tmp);

    let read1 = tmp.join("reads_1.fq");
    let read2 = tmp.join("reads_2.fq");
    let gc_reads = tmp.join("reads_gc.fq");
    let dataset = match helpers::build_synthetic_index() {
        Some(dataset) => dataset,
        None => return,
    };
    let transcripts = dataset.transcripts;

    let mut rng = SimpleRng::new(0xA5A5_1234_1111_2222);
    let pairs = synthesize_paired_reads(&mut rng, &transcripts, 40, 80, 12);
    write_fastq(&read1, &pairs.iter().map(|p| &p.0).collect::<Vec<_>>()).expect("write read1");
    write_fastq(&read2, &pairs.iter().map(|p| &p.1).collect::<Vec<_>>()).expect("write read2");

    let gc = synthesize_reads_gc_biased(&mut rng, &transcripts, 40, 20, 0.6);
    write_fastq(&gc_reads, &gc.iter().collect::<Vec<_>>()).expect("write gc reads");

    assert_kallisto_paired_parity(&dataset.index_path, &read1, &read2);
    assert_kallisto_single_overhang_parity(&dataset.index_path, &gc_reads);

    let _ = fs::remove_dir_all(&tmp);
}

fn assert_kallisto_paired_parity(index: &Path, reads1: &Path, reads2: &Path) {
    let bifrost_index = build_bifrost_index(index).expect("bifrost index");

    let bifrost_counts = {
        let file1 = File::open(reads1).expect("reads1 file");
        let file2 = File::open(reads2).expect("reads2 file");
        let mut reader1 = FastqReader::new(BufReader::new(file1));
        let mut reader2 = FastqReader::new(BufReader::new(file2));
        pseudoalign_paired_bifrost_with_options(
            &bifrost_index,
            &mut reader1,
            &mut reader2,
            Strand::Unstranded,
            PseudoalignOptions::default(),
        )
        .expect("paired bifrost counts")
    };

    let kallisto_aligned =
        run_kallisto_quant_paired(index, reads1, reads2).expect("kallisto synthetic paired count");
    assert_eq!(
        bifrost_counts.reads_aligned,
        kallisto_aligned,
        "bifrost/kallisto synthetic paired parity mismatch for {} and {}",
        reads1.display(),
        reads2.display()
    );
}

fn assert_kallisto_single_overhang_parity(index: &Path, reads: &Path) {
    let bifrost_index = build_bifrost_index(index).expect("bifrost index");

    let bifrost_counts = {
        let file = File::open(reads).expect("reads file");
        let mut reader = FastqReader::new(BufReader::new(file));
        pseudoalign_single_end_bifrost(&bifrost_index, &mut reader).expect("bifrost counts")
    };

    let kallisto_aligned = run_kallisto_quant_single_overhang(index, reads)
        .expect("kallisto synthetic single-overhang count");
    assert_eq!(
        bifrost_counts.reads_aligned,
        kallisto_aligned,
        "bifrost/kallisto synthetic single-overhang parity mismatch for {}",
        reads.display()
    );
}

fn run_kallisto_quant_single_overhang(index: &Path, reads: &Path) -> Option<u64> {
    let out_dir = temp_dir("kallistors-variants-kallisto");
    let _ = fs::remove_dir_all(&out_dir);
    fs::create_dir_all(&out_dir).ok()?;
    let status = Command::new("kallisto")
        .arg("quant")
        .arg("--single")
        .arg("-l")
        .arg("100")
        .arg("-s")
        .arg("20")
        .arg("--single-overhang")
        .arg("-i")
        .arg(index)
        .arg("-o")
        .arg(&out_dir)
        .arg(reads)
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).ok()?;
    let key = "\"n_pseudoaligned\"";
    let idx = run_info.find(key)?;
    let after = &run_info[idx + key.len()..];
    let colon = after.find(':')?;
    let rest = after[colon + 1..].trim_start();
    let mut digits = String::new();
    for c in rest.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else if !digits.is_empty() {
            break;
        }
    }
    digits.parse().ok()
}

fn run_kallisto_quant_paired(index: &Path, reads1: &Path, reads2: &Path) -> Option<u64> {
    let out_dir = temp_dir("kallistors-variants-kallisto-paired");
    let _ = fs::remove_dir_all(&out_dir);
    fs::create_dir_all(&out_dir).ok()?;
    let status = Command::new("kallisto")
        .arg("quant")
        .arg("-i")
        .arg(index)
        .arg("-o")
        .arg(&out_dir)
        .arg(reads1)
        .arg(reads2)
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).ok()?;
    let key = "\"n_pseudoaligned\"";
    let idx = run_info.find(key)?;
    let after = &run_info[idx + key.len()..];
    let colon = after.find(':')?;
    let rest = after[colon + 1..].trim_start();
    let mut digits = String::new();
    for c in rest.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else if !digits.is_empty() {
            break;
        }
    }
    digits.parse().ok()
}

fn temp_dir(prefix: &str) -> PathBuf {
    let pid = std::process::id();
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!("{}-{}-{}", prefix, pid, now))
}

fn write_fastq(path: &Path, reads: &[&Vec<u8>]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (i, seq) in reads.iter().enumerate() {
        writeln!(file, "@r{}", i + 1)?;
        writeln!(file, "{}", String::from_utf8_lossy(seq))?;
        writeln!(file, "+")?;
        writeln!(file, "{}", "I".repeat(seq.len()))?;
    }
    Ok(())
}

fn synthesize_paired_reads(
    rng: &mut SimpleRng,
    transcripts: &[Vec<u8>],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let seq = &transcripts[rng.next_usize(transcripts.len())];
        if seq.len() < frag_len {
            continue;
        }
        let start = rng.next_usize(seq.len() - frag_len + 1);
        let frag = &seq[start..start + frag_len];
        let read1 = frag[..read_len].to_vec();
        let read2 = revcomp(&frag[frag_len - read_len..]);
        out.push((read1, read2));
    }
    out
}

fn synthesize_reads_gc_biased(
    rng: &mut SimpleRng,
    transcripts: &[Vec<u8>],
    read_len: usize,
    count: usize,
    gc_min: f64,
) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(count);
    let candidates: Vec<&Vec<u8>> = transcripts
        .iter()
        .filter(|seq| gc_fraction(seq) >= gc_min)
        .collect();
    let pool: Vec<&Vec<u8>> = if candidates.is_empty() {
        transcripts.iter().collect()
    } else {
        candidates
    };
    for _ in 0..count {
        let seq = pool[rng.next_usize(pool.len())];
        if seq.len() < read_len {
            continue;
        }
        let start = rng.next_usize(seq.len() - read_len + 1);
        out.push(seq[start..start + read_len].to_vec());
    }
    out
}

fn gc_fraction(seq: &[u8]) -> f64 {
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
        .count();
    if seq.is_empty() {
        0.0
    } else {
        gc as f64 / seq.len() as f64
    }
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

struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn next_usize(&mut self, bound: usize) -> usize {
        if bound == 0 {
            return 0;
        }
        (self.next_u64() as usize) % bound
    }
}
