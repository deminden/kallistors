use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    build_bifrost_index_with_positions, pseudoalign_single_end_bifrost_with_options,
    FragmentFilter, PseudoalignOptions, Strand, StrandSpecific,
};

fn kallisto_available() -> bool {
    Command::new("kallisto")
        .arg("version")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn run_kallisto_quant(index: &Path, reads: &Path, out_dir: &Path, strand: Option<&str>) -> u64 {
    let mut cmd = Command::new("kallisto");
    cmd.arg("quant")
        .arg("--single")
        .arg("-l")
        .arg("100")
        .arg("-s")
        .arg("20")
        .arg("-i")
        .arg(index)
        .arg("-o")
        .arg(out_dir)
        .arg(reads);
    if let Some(mode) = strand {
        cmd.arg(mode);
    }
    let status = cmd.status().expect("run kallisto quant");
    assert!(status.success(), "kallisto quant failed");
    let run_info = fs::read_to_string(out_dir.join("run_info.json")).expect("read run_info");
    parse_n_pseudoaligned(&run_info).expect("parse n_pseudoaligned")
}

fn parse_n_pseudoaligned(run_info: &str) -> Option<u64> {
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

#[test]
fn shade_strand_parity_with_kallisto() {
    if !kallisto_available() {
        return;
    }

    let tmp = temp_dir("kallistors-shade-parity");
    let _ = fs::remove_dir_all(&tmp);
    fs::create_dir_all(&tmp).expect("create temp dir");

    let fasta = tmp.join("shade_transcripts.fa");
    let reads = tmp.join("shade_reads.fq");
    let index = tmp.join("shade.idx");

    let transcripts = synthetic_transcripts();
    write_shade_fasta(&fasta, &transcripts).expect("write shade fasta");
    write_reads(&reads, &transcripts).expect("write reads");

    let status = Command::new("kallisto")
        .arg("index")
        .arg("-i")
        .arg(&index)
        .arg(&fasta)
        .status()
        .expect("run kallisto index");
    assert!(status.success(), "kallisto index failed");

    let k_unstranded = run_kallisto_quant(&index, &reads, &tmp.join("k_unstranded"), None);
    let k_fr = run_kallisto_quant(&index, &reads, &tmp.join("k_fr"), Some("--fr-stranded"));
    let k_rf = run_kallisto_quant(&index, &reads, &tmp.join("k_rf"), Some("--rf-stranded"));

    let index = build_bifrost_index_with_positions(&index, true).expect("load index");
    let filter = Some(FragmentFilter {
        fragment_length: 100,
        single_overhang: false,
    });

    let mut reader = FastqReader::new(BufReader::new(File::open(&reads).unwrap()));
    let ours = pseudoalign_single_end_bifrost_with_options(
        &index,
        &mut reader,
        Strand::Unstranded,
        filter,
        PseudoalignOptions::default(),
    )
    .expect("pseudoalign unstranded");
    assert_eq!(ours.reads_aligned, k_unstranded);

    let mut reader = FastqReader::new(BufReader::new(File::open(&reads).unwrap()));
    let ours = pseudoalign_single_end_bifrost_with_options(
        &index,
        &mut reader,
        Strand::Unstranded,
        filter,
        PseudoalignOptions {
            strand_specific: Some(StrandSpecific::FR),
            ..PseudoalignOptions::default()
        },
    )
    .expect("pseudoalign fr");
    assert_eq!(ours.reads_aligned, k_fr);

    let mut reader = FastqReader::new(BufReader::new(File::open(&reads).unwrap()));
    let ours = pseudoalign_single_end_bifrost_with_options(
        &index,
        &mut reader,
        Strand::Unstranded,
        filter,
        PseudoalignOptions {
            strand_specific: Some(StrandSpecific::RF),
            ..PseudoalignOptions::default()
        },
    )
    .expect("pseudoalign rf");
    assert_eq!(ours.reads_aligned, k_rf);

    let _ = fs::remove_dir_all(&tmp);
}

fn write_shade_fasta(path: &Path, transcripts: &[(String, Vec<u8>)]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (name, seq) in transcripts {
        writeln!(file, ">{}", name)?;
        writeln!(file, "{}", String::from_utf8_lossy(seq))?;
        writeln!(file, ">{}_shade_alt", name)?;
        writeln!(file, "{}", String::from_utf8_lossy(seq))?;
    }
    Ok(())
}

fn write_reads(path: &Path, transcripts: &[(String, Vec<u8>)]) -> std::io::Result<()> {
    let mut rng = SimpleRng::new(0x1234_5678_9abc_def0);
    let mut reads = Vec::new();

    for i in 0..120 {
        let seq = &transcripts[rng.next_usize(transcripts.len())].1;
        if seq.len() < 90 {
            continue;
        }
        let start = rng.next_usize(seq.len() - 90 + 1);
        reads.push((
            format!("fwd_{i}"),
            mutate(&mut rng, &seq[start..start + 90], 0.02),
        ));
    }

    for i in 0..60 {
        let seq = &transcripts[rng.next_usize(transcripts.len())].1;
        if seq.len() < 80 {
            continue;
        }
        let start = rng.next_usize(seq.len() - 80 + 1);
        reads.push((format!("rev_{i}"), revcomp(&seq[start..start + 80])));
    }

    for i in 0..30 {
        reads.push((format!("rand_{i}"), random_seq(&mut rng, 70)));
    }

    let mut file = File::create(path)?;
    for (name, seq) in reads {
        writeln!(file, "@{name}")?;
        writeln!(file, "{}", String::from_utf8_lossy(&seq))?;
        writeln!(file, "+")?;
        writeln!(file, "{}", "I".repeat(seq.len()))?;
    }
    Ok(())
}

fn synthetic_transcripts() -> Vec<(String, Vec<u8>)> {
    let mut rng = SimpleRng::new(0x9e37_79b9_7f4a_7c15);
    let mut out = Vec::new();
    for i in 0..6 {
        let seq = random_seq(&mut rng, 300);
        out.push((format!("tx{i}"), seq));
    }
    out
}

fn mutate(rng: &mut SimpleRng, seq: &[u8], rate: f64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        if rng.next_f64() < rate {
            let mut next = bases[rng.next_usize(bases.len())];
            while next == b {
                next = bases[rng.next_usize(bases.len())];
            }
            out.push(next);
        } else {
            out.push(b);
        }
    }
    out
}

fn random_seq(rng: &mut SimpleRng, len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        out.push(bases[rng.next_usize(bases.len())]);
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
            0
        } else {
            ((self.next_u64() >> 32) as usize) % bound
        }
    }

    fn next_f64(&mut self) -> f64 {
        let val = self.next_u64() >> 11;
        (val as f64) / ((1u64 << 53) as f64)
    }
}

fn temp_dir(prefix: &str) -> PathBuf {
    let pid = std::process::id();
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!("{}-{}-{}", prefix, pid, now))
}
