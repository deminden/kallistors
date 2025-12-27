use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};

use kallistors::io::FastqReader;
use kallistors::pseudoalign::{
    build_bifrost_index, build_kmer_ec_index, pseudoalign_single_end,
    pseudoalign_single_end_bifrost, EcCounts,
};

#[test]
fn parity_on_synthetic_variants() {
    let tmp = temp_dir("kallistors-variants");
    let _ = fs::create_dir_all(&tmp);

    let read1 = tmp.join("reads_1.fq");
    let read2 = tmp.join("reads_2.fq");
    let gc_reads = tmp.join("reads_gc.fq");
    let index = fixture_path("synthetic.idx");
    let fasta = fixture_path("synthetic_transcripts.fa");

    let transcripts = read_transcripts(&fasta);

    let mut rng = SimpleRng::new(0xA5A5_1234_1111_2222);
    let pairs = synthesize_paired_reads(&mut rng, &transcripts, 40, 80, 12);
    write_fastq(&read1, &pairs.iter().map(|p| &p.0).collect::<Vec<_>>()).expect("write read1");
    write_fastq(&read2, &pairs.iter().map(|p| &p.1).collect::<Vec<_>>()).expect("write read2");

    let gc = synthesize_reads_gc_biased(&mut rng, &transcripts, 40, 20, 0.6);
    write_fastq(&gc_reads, &gc.iter().collect::<Vec<_>>()).expect("write gc reads");

    assert_parity(&index, &read1);
    assert_parity(&index, &read2);
    assert_parity(&index, &gc_reads);

    let _ = fs::remove_dir_all(&tmp);
}

fn assert_parity(index: &Path, reads: &Path) {
    let naive_index = build_kmer_ec_index(index).expect("naive index");
    let bifrost_index = build_bifrost_index(index).expect("bifrost index");

    let naive_counts = {
        let file = File::open(reads).expect("reads file");
        let mut reader = FastqReader::new(BufReader::new(file));
        pseudoalign_single_end(&naive_index, &mut reader).expect("naive counts")
    };

    let bifrost_counts = {
        let file = File::open(reads).expect("reads file");
        let mut reader = FastqReader::new(BufReader::new(file));
        pseudoalign_single_end_bifrost(&bifrost_index, &mut reader).expect("bifrost counts")
    };

    assert_eq!(naive_counts.reads_processed, bifrost_counts.reads_processed);
    assert_eq!(naive_counts.reads_aligned, bifrost_counts.reads_aligned);
    assert_eq!(normalize_ec(&naive_counts), normalize_ec(&bifrost_counts));
}

fn normalize_ec(counts: &EcCounts) -> HashMap<Vec<u32>, u32> {
    let mut map = HashMap::new();
    for (ec, &count) in counts.ec_list.iter().zip(counts.counts.iter()) {
        let mut key = ec.clone();
        key.sort_unstable();
        map.insert(key, count);
    }
    map
}

fn temp_dir(prefix: &str) -> PathBuf {
    let pid = std::process::id();
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!("{}-{}-{}", prefix, pid, now))
}

fn fixture_path(name: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests");
    path.push("fixtures");
    path.push(name);
    path
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

fn read_transcripts(path: &Path) -> Vec<(String, Vec<u8>)> {
    let mut fasta = String::new();
    File::open(path)
        .expect("fasta file")
        .read_to_string(&mut fasta)
        .expect("read fasta");
    parse_fasta(&fasta)
}

fn synthesize_paired_reads(
    rng: &mut SimpleRng,
    transcripts: &[(String, Vec<u8>)],
    read_len: usize,
    frag_len: usize,
    count: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let seq = &transcripts[rng.next_usize(transcripts.len())].1;
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
    transcripts: &[(String, Vec<u8>)],
    read_len: usize,
    count: usize,
    gc_min: f64,
) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(count);
    let candidates: Vec<&Vec<u8>> = transcripts
        .iter()
        .map(|t| &t.1)
        .filter(|seq| gc_fraction(seq) >= gc_min)
        .collect();
    let pool: Vec<&Vec<u8>> = if candidates.is_empty() {
        transcripts.iter().map(|t| &t.1).collect()
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

fn parse_fasta(input: &str) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    let mut name = None;
    let mut seq = Vec::new();
    for line in input.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.push((n, seq));
                seq = Vec::new();
            }
            name = Some(rest.trim().to_string());
            continue;
        }
        seq.extend_from_slice(line.trim().as_bytes());
    }
    if let Some(n) = name {
        out.push((n, seq));
    }
    out
}
