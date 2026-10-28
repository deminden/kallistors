use super::workers::IndexWorkers;
use std::collections::HashMap;
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::hash::Hasher;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use flate2::read::GzDecoder;
use rayon::prelude::*;

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
    pub aa: bool,
}

impl Default for IndexBuildOptions {
    fn default() -> Self {
        Self {
            k: 31,
            g: None,
            threads: 1,
            make_unique: false,
            ec_max_size: -1,
            aa: false,
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
    /// Serialization and file I/O, excluding MPHF construction.
    pub write: Duration,
    pub transcripts: usize,
    pub unitigs: usize,
    pub kmers: usize,
    pub minimizers: usize,
    /// Logical edge-table scratch size, excluding allocator overhead and RSS.
    pub edge_scratch_bytes: usize,
}

#[derive(Clone)]
pub(super) struct Transcript {
    pub(super) name: String,
    pub(super) original_len: u32,
    pub(super) seq: Vec<u8>,
}

#[derive(Default)]
struct NameFingerprintHasher(u64);
impl Hasher for NameFingerprintHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write_u64(&mut self, value: u64) {
        self.0 = value;
    }
    fn write(&mut self, bytes: &[u8]) {
        self.0 = super::bifrost::wyhash(bytes, 0);
    }
}

#[derive(Default)]
struct TranscriptNames {
    first: HashMap<u64, usize, BuildHasherDefault<NameFingerprintHasher>>,
    collisions: HashMap<u64, Vec<usize>>,
}

impl TranscriptNames {
    fn insert(&mut self, name: &str, transcripts: &[Transcript]) -> bool {
        self.insert_hashed(
            super::bifrost::wyhash(name.as_bytes(), 0),
            name,
            transcripts,
        )
    }
    fn insert_hashed(&mut self, hash: u64, name: &str, transcripts: &[Transcript]) -> bool {
        match self.first.entry(hash) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(transcripts.len());
                true
            }
            std::collections::hash_map::Entry::Occupied(entry) => {
                if transcripts[*entry.get()].name == name
                    || self
                        .collisions
                        .get(&hash)
                        .is_some_and(|ids| ids.iter().any(|&id| transcripts[id].name == name))
                {
                    false
                } else {
                    self.collisions
                        .entry(hash)
                        .or_default()
                        .push(transcripts.len());
                    true
                }
            }
        }
    }
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

    let workers = IndexWorkers::new(options.threads)?;
    workers.install(|| {
        // Keep FASTA allocations on the same requested worker pool that owns
        // normalization, mapping and disposal, allowing allocator reuse there.
        let fasta_start = Instant::now();
        let transcripts = load_transcripts(
            fasta_paths,
            options.make_unique,
            options.aa,
            options.threads,
        )?;
        let fasta_parse = fasta_start.elapsed();

        let (mut graph, graph_report) = build_kmer_unitig_graph_with_report(
            transcripts,
            k,
            g,
            options.ec_max_size,
            options.threads,
        )?;
        let write_report = write_index_with_report(output, &mut graph, options.threads)?;
        let unitigs = graph.unitigs.len() + graph.km_unitigs.len();
        let kmers = graph
            .unitigs
            .iter()
            .map(|seq| seq.len().saturating_sub(graph.k) + 1)
            .sum::<usize>()
            + graph.km_unitigs.len();
        let minimizers = graph.minimizer_keys.len();
        let mut report = IndexBuildReport {
            total: total_start.elapsed(),
            fasta_parse,
            graph_build: graph_report.graph_build,
            minimizer_mphf: graph_report.minimizer_index + write_report.mphf,
            ec_build: graph_report.ec_build,
            write: write_report.total - write_report.mphf,
            transcripts: graph.transcript_names.len(),
            unitigs,
            kmers,
            minimizers,
            edge_scratch_bytes: graph_report.edge_scratch_bytes,
        };
        if options.threads > 1 {
            // Materialization and payload workers own many of these allocations.
            // Keep disposal bounded by the same requested worker pool.
            rayon::join(
                || graph.unitigs.into_par_iter().for_each(drop),
                || {
                    rayon::join(
                        || graph.nodes.into_par_iter().for_each(drop),
                        || graph.transcript_names.into_par_iter().for_each(drop),
                    )
                },
            );
            drop((
                graph.km_unitigs,
                graph.minimizer_keys,
                graph.unitig_bitmap_blocks,
                graph.km_bitmap_blocks,
                graph.transcript_lengths,
            ));
        } else {
            drop(graph);
        }
        report.total = total_start.elapsed();
        Ok(report)
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

fn load_transcripts(
    fasta_paths: &[PathBuf],
    make_unique: bool,
    aa: bool,
    threads: usize,
) -> Result<Vec<Transcript>> {
    let mut out = Vec::new();
    let mut names = TranscriptNames::default();
    let mut rng = Mt19937::new(42);

    for path in fasta_paths {
        let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
        let is_gz = path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));
        if is_gz {
            read_transcript_file(
                BufReader::with_capacity(256 * 1024, GzDecoder::new(file)),
                &mut out,
                &mut names,
                &mut rng,
                make_unique,
                aa,
                threads > 1,
            )?;
        } else {
            read_transcript_file(
                BufReader::with_capacity(256 * 1024, file),
                &mut out,
                &mut names,
                &mut rng,
                make_unique,
                aa,
                threads > 1,
            )?;
        }
    }
    if threads > 1 {
        let workers = IndexWorkers::new(threads)?;
        let needs_masking = workers.install(|| {
            out.par_iter_mut()
                .map(|transcript| {
                    if aa {
                        transcript.seq = crate::util::aa_to_comma_free(&transcript.seq);
                        transcript.original_len = transcript.seq.len() as u32;
                        trim_poly_a(&mut transcript.seq);
                        false
                    } else {
                        transcript.seq.make_ascii_uppercase();
                        let ambiguous = transcript
                            .seq
                            .chunks(32)
                            .any(|chunk| crate::dna::encode_kmer_pair(chunk).is_none());
                        if !ambiguous {
                            trim_poly_a(&mut transcript.seq);
                        }
                        ambiguous
                    }
                })
                .collect::<Vec<_>>()
        });
        // Random substitutions consume the same RNG in transcript/base order.
        for (transcript, mask) in out.iter_mut().zip(needs_masking) {
            if mask {
                for base in &mut transcript.seq {
                    match *base {
                        b'A' | b'C' | b'G' | b'T' => {}
                        b'U' => *base = b'T',
                        _ => *base = dna_from_mt(rng.next_u32()),
                    }
                }
                trim_poly_a(&mut transcript.seq);
            }
        }
    }
    Ok(out)
}

fn read_transcript_file<R: BufRead>(
    mut reader: R,
    out: &mut Vec<Transcript>,
    names: &mut TranscriptNames,
    rng: &mut Mt19937,
    make_unique: bool,
    aa: bool,
    defer_normalization: bool,
) -> Result<()> {
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();
    let mut process_line = |bytes: &[u8], known_ascii: bool| -> Result<()> {
        if known_ascii && !bytes.starts_with(b">") {
            let mut end = bytes.len();
            while end != 0 && matches!(bytes[end - 1], b'\n' | b'\r') {
                end -= 1;
            }
            current_seq.extend_from_slice(&bytes[..end]);
            return Ok(());
        }
        let line = if known_ascii {
            // SAFETY: every contributing input buffer was checked as ASCII.
            // Pending lines retain that proof only if all fragments did.
            unsafe { std::str::from_utf8_unchecked(bytes) }
        } else {
            std::str::from_utf8(bytes).map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "stream did not contain valid UTF-8",
                )
            })?
        };
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if let Some(rest) = trimmed.strip_prefix('>') {
            flush_transcript(
                out,
                names,
                rng,
                current_name.take(),
                &mut current_seq,
                make_unique,
                aa,
                defer_normalization,
            )?;
            let name = rest.split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
        } else {
            current_seq.extend_from_slice(trimmed.as_bytes());
        }
        Ok(())
    };
    let mut pending = Vec::new();
    let mut pending_ascii = true;
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            if !pending.is_empty() {
                process_line(&pending, pending_ascii)?;
            }
            break;
        }
        let ascii = buffer.is_ascii();
        let mut start = 0;
        for end in memchr::memchr_iter(b'\n', buffer) {
            let part = &buffer[start..=end];
            if pending.is_empty() {
                process_line(part, ascii)?;
            } else {
                pending.extend_from_slice(part);
                process_line(&pending, pending_ascii && ascii)?;
                pending.clear();
                pending_ascii = true;
            }
            start = end + 1;
        }
        pending.extend_from_slice(&buffer[start..]);
        if start < buffer.len() {
            pending_ascii &= ascii;
        }
        let consumed = buffer.len();
        reader.consume(consumed);
    }
    flush_transcript(
        out,
        names,
        rng,
        current_name.take(),
        &mut current_seq,
        make_unique,
        aa,
        defer_normalization,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn flush_transcript(
    out: &mut Vec<Transcript>,
    names: &mut TranscriptNames,
    rng: &mut Mt19937,
    name: Option<String>,
    seq: &mut Vec<u8>,
    make_unique: bool,
    aa: bool,
    defer_normalization: bool,
) -> Result<()> {
    let Some(mut name) = name else {
        return Ok(());
    };
    if name.is_empty() {
        name = out.len().to_string();
    }
    if !names.insert(&name, out) {
        if !make_unique {
            return Err(Error::InvalidFormat(format!(
                "repeated name in FASTA file: {name}; rerun with --make-unique"
            )));
        }
        let base = name;
        let mut suffix = 1usize;
        loop {
            let candidate = format!("{base}_{suffix}");
            if names.insert(&candidate, out) {
                name = candidate;
                break;
            }
            suffix += 1;
        }
    }

    let mut normalized = if defer_normalization {
        seq.clone()
    } else if aa {
        crate::util::aa_to_comma_free(seq)
    } else {
        let mut normalized = seq.clone();
        normalized.make_ascii_uppercase();
        for chunk in normalized.chunks_mut(32) {
            if crate::dna::encode_kmer_pair(chunk).is_some() {
                continue;
            }
            for base in chunk {
                match *base {
                    b'A' | b'C' | b'G' | b'T' => {}
                    b'U' => *base = b'T',
                    _ => *base = dna_from_mt(rng.next_u32()),
                }
            }
        }
        normalized
    };
    // Quantification retains the pre-trimming length even though graph keys
    // are collected from the normalized, poly-A-trimmed sequence.
    let original_len = normalized.len() as u32;
    if !defer_normalization {
        trim_poly_a(&mut normalized);
    }
    out.push(Transcript {
        name,
        original_len,
        seq: normalized,
    });
    seq.clear();
    Ok(())
}

fn trim_poly_a(seq: &mut Vec<u8>) {
    // Ten trailing As trigger removal of the entire tail, not just ten bases.
    // Shorter tails remain part of the reference, matching the index convention.
    if seq.len() >= 10 && seq[seq.len() - 10..].iter().all(|&b| b == b'A') {
        while seq.last() == Some(&b'A') {
            seq.pop();
        }
    }
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
    fn transcript_name_hash_collisions_retain_exact_membership() {
        let mut names = TranscriptNames::default();
        let mut transcripts = Vec::new();
        for name in ["alpha", "beta", "alpha_1"] {
            assert!(names.insert_hashed(7, name, &transcripts));
            transcripts.push(Transcript {
                name: name.into(),
                original_len: 0,
                seq: Vec::new(),
            });
        }
        for name in ["alpha", "beta", "alpha_1"] {
            assert!(!names.insert_hashed(7, name, &transcripts));
        }
        assert!(names.insert_hashed(7, "gamma", &transcripts));
    }

    #[test]
    fn minimizer_defaults_match_kallisto() {
        assert_eq!(default_minimizer_len(31), 23);
        assert_eq!(default_minimizer_len(19), 13);
        assert_eq!(default_minimizer_len(17), 13);
        assert_eq!(default_minimizer_len(13), 11);
    }

    #[test]
    fn buffered_fasta_lines_preserve_boundaries_and_utf8_errors() {
        let input =
            ">alpha\u{2003}description\r\nacguNN\r\n\r\n>beta\nCCGTA\r\r\n>\nTGCA".as_bytes();
        let mut expected = None;
        for capacity in [65536, 1, 2, 3, 7, 64] {
            let mut out = Vec::new();
            let mut names = TranscriptNames::default();
            let mut rng = Mt19937::new(42);
            read_transcript_file(
                BufReader::with_capacity(capacity, input),
                &mut out,
                &mut names,
                &mut rng,
                false,
                false,
                false,
            )
            .unwrap();
            let actual = out
                .into_iter()
                .map(|tr| (tr.name, tr.original_len, tr.seq))
                .collect::<Vec<_>>();
            assert_eq!(actual.len(), 3);
            assert_eq!(actual[0].0, "alpha");
            assert_eq!(actual[2].0, "2");
            if let Some(expected) = &expected {
                assert_eq!(&actual, expected, "capacity={capacity}");
            } else {
                expected = Some(actual);
            }
            for invalid in [b">bad\xff\nACGT".as_slice(), b">good\nAC\xffGT".as_slice()] {
                let mut out = Vec::new();
                let mut names = TranscriptNames::default();
                let mut rng = Mt19937::new(42);
                assert!(
                    read_transcript_file(
                        BufReader::with_capacity(capacity, invalid),
                        &mut out,
                        &mut names,
                        &mut rng,
                        false,
                        false,
                        false
                    )
                    .is_err()
                );
            }
        }
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
