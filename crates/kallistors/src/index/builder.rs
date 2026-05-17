use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use flate2::read::GzDecoder;

use super::graph_build::build_kmer_unitig_graph_with_report;
use super::writer::write_index_with_report;
use crate::{Error, Result};

const MAX_K: usize = 32;

#[derive(Debug, Clone)]
pub struct IndexBuildOptions {
    pub k: usize,
    pub g: Option<usize>,
    pub threads: usize,
    pub make_unique: bool,
    pub ec_max_size: i32,
}

impl Default for IndexBuildOptions {
    fn default() -> Self {
        Self {
            k: 31,
            g: None,
            threads: 1,
            make_unique: false,
            ec_max_size: -1,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct IndexBuildReport {
    pub total: Duration,
    pub fasta_parse: Duration,
    pub graph_build: Duration,
    pub minimizer_mphf: Duration,
    pub ec_build: Duration,
    pub write: Duration,
    pub transcripts: usize,
    pub unitigs: usize,
    pub kmers: usize,
    pub minimizers: usize,
}

#[derive(Clone)]
pub(super) struct Transcript {
    pub(super) name: String,
    pub(super) original_len: u32,
    pub(super) seq: Vec<u8>,
}

pub fn build_index(
    output: &Path,
    fasta_paths: &[PathBuf],
    options: IndexBuildOptions,
) -> Result<()> {
    build_index_with_report(output, fasta_paths, options).map(|_| ())
}

pub fn build_index_with_report(
    output: &Path,
    fasta_paths: &[PathBuf],
    options: IndexBuildOptions,
) -> Result<IndexBuildReport> {
    let total_start = Instant::now();
    validate_options(&options, fasta_paths)?;
    let k = options.k;
    let g = options.g.unwrap_or_else(|| default_minimizer_len(k));
    validate_k_g(k, g)?;

    let fasta_start = Instant::now();
    let transcripts = load_transcripts(fasta_paths, options.make_unique)?;
    let fasta_parse = fasta_start.elapsed();

    let (graph, graph_report) = build_kmer_unitig_graph_with_report(
        &transcripts,
        k,
        g,
        options.ec_max_size,
        options.threads,
    )?;
    let write_report = write_index_with_report(output, &graph)?;
    let unitigs = graph.unitigs.len() + graph.km_unitigs.len();
    let kmers = graph
        .unitigs
        .iter()
        .map(|seq| seq.len().saturating_sub(graph.k) + 1)
        .sum::<usize>()
        + graph.km_unitigs.len();
    let minimizers = graph.minimizer_keys.len();
    Ok(IndexBuildReport {
        total: total_start.elapsed(),
        fasta_parse,
        graph_build: graph_report.graph_build,
        minimizer_mphf: graph_report.minimizer_index + write_report.mphf,
        ec_build: graph_report.ec_build,
        write: write_report.total,
        transcripts: graph.transcript_names.len(),
        unitigs,
        kmers,
        minimizers,
    })
}

pub fn default_minimizer_len(k: usize) -> usize {
    if k <= 13 {
        k - 2
    } else if k <= 17 {
        k - 4
    } else if k <= 19 {
        k - 6
    } else {
        k - 8
    }
}

fn validate_options(options: &IndexBuildOptions, fasta_paths: &[PathBuf]) -> Result<()> {
    if options.threads == 0 {
        return Err(Error::InvalidFormat(
            "threads must be greater than zero".into(),
        ));
    }
    if fasta_paths.is_empty() {
        return Err(Error::InvalidFormat("no FASTA files specified".into()));
    }
    if options.ec_max_size < -1 {
        return Err(Error::InvalidFormat("invalid max EC size".into()));
    }
    for path in fasta_paths {
        if !path.exists() {
            return Err(Error::MissingFile(path.clone()));
        }
    }
    Ok(())
}

fn validate_k_g(k: usize, g: usize) -> Result<()> {
    if !(3..MAX_K).contains(&k) {
        return Err(Error::InvalidFormat(format!(
            "invalid k-mer length {k}, minimum is 3 and maximum is {}",
            MAX_K - 1
        )));
    }
    if k.is_multiple_of(2) {
        return Err(Error::InvalidFormat("k needs to be an odd number".into()));
    }
    if g <= 2 || g > k - 2 {
        return Err(Error::InvalidFormat(format!(
            "invalid minimizer size {g}, minimum is 3 and maximum is k - 2"
        )));
    }
    Ok(())
}

fn load_transcripts(fasta_paths: &[PathBuf], make_unique: bool) -> Result<Vec<Transcript>> {
    let mut out = Vec::new();
    let mut names = HashSet::new();
    let mut rng = Mt19937::new(42);

    for path in fasta_paths {
        let mut reader = open_fasta_reader(path)?;
        let mut line = String::new();
        let mut current_name: Option<String> = None;
        let mut current_seq = Vec::new();

        loop {
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if let Some(rest) = trimmed.strip_prefix('>') {
                flush_transcript(
                    &mut out,
                    &mut names,
                    &mut rng,
                    current_name.take(),
                    &mut current_seq,
                    make_unique,
                )?;
                let name = rest.split_whitespace().next().unwrap_or("").to_string();
                current_name = Some(name);
            } else {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
        }
        flush_transcript(
            &mut out,
            &mut names,
            &mut rng,
            current_name.take(),
            &mut current_seq,
            make_unique,
        )?;
    }

    Ok(out)
}

fn open_fasta_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let is_gz = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));
    if is_gz {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn flush_transcript(
    out: &mut Vec<Transcript>,
    names: &mut HashSet<String>,
    rng: &mut Mt19937,
    name: Option<String>,
    seq: &mut Vec<u8>,
    make_unique: bool,
) -> Result<()> {
    let Some(mut name) = name else {
        return Ok(());
    };
    if name.is_empty() {
        name = out.len().to_string();
    }
    if !names.insert(name.clone()) {
        if !make_unique {
            return Err(Error::InvalidFormat(format!(
                "repeated name in FASTA file: {name}; rerun with --make-unique"
            )));
        }
        let base = name;
        let mut suffix = 1usize;
        loop {
            let candidate = format!("{base}_{suffix}");
            if names.insert(candidate.clone()) {
                name = candidate;
                break;
            }
            suffix += 1;
        }
    }

    let original_len = seq.len() as u32;
    let mut normalized = Vec::with_capacity(seq.len());
    for &base in seq.iter() {
        let b = base.to_ascii_uppercase();
        match b {
            b'A' | b'C' | b'G' | b'T' => normalized.push(b),
            b'U' => normalized.push(b'T'),
            _ => normalized.push(dna_from_mt(rng.next_u32())),
        }
    }
    if normalized.len() >= 10
        && normalized[normalized.len() - 10..]
            .iter()
            .all(|&b| b == b'A')
    {
        while normalized.last() == Some(&b'A') {
            normalized.pop();
        }
    }
    out.push(Transcript {
        name,
        original_len,
        seq: normalized,
    });
    seq.clear();
    Ok(())
}

fn dna_from_mt(v: u32) -> u8 {
    b"ACGT"[(v & 0x03) as usize]
}

struct Mt19937 {
    mt: [u32; 624],
    index: usize,
}

impl Mt19937 {
    fn new(seed: u32) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = seed;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { mt, index: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            self.twist();
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    fn twist(&mut self) {
        for i in 0..624 {
            let y = (self.mt[i] & 0x8000_0000) | (self.mt[(i + 1) % 624] & 0x7fff_ffff);
            let mut next = self.mt[(i + 397) % 624] ^ (y >> 1);
            if (y & 1) != 0 {
                next ^= 0x9908_b0df;
            }
            self.mt[i] = next;
        }
        self.index = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizer_defaults_match_kallisto() {
        assert_eq!(default_minimizer_len(31), 23);
        assert_eq!(default_minimizer_len(19), 13);
        assert_eq!(default_minimizer_len(17), 13);
        assert_eq!(default_minimizer_len(13), 11);
    }

    #[test]
    fn mt19937_matches_kallisto_dna_masking_pattern() {
        let mut rng = Mt19937::new(42);
        let vals = (0..4)
            .map(|_| dna_from_mt(rng.next_u32()) as char)
            .collect::<String>();
        assert_eq!(vals, "GTAG");
    }
}
