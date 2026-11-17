# kallistors

A high-performance Rust implementation of kallisto-style pseudoalignment and RNA-seq transcript quantification.

Build kallisto v13-compatible indexes, quantify paired-end or single-end reads, and generate
single-cell BUS output. Use kallistors as a command-line tool or a Rust library.

In the published eight-thread benchmarks, full quantification was **6.1–15.4× faster on
GENCODE datasets** and **2.8–3.5× faster on mature-plus-nascent RNA** than original kallisto
0.52.0. Index construction was **9.2–10.9× faster**.

## Performance

The existing 0.6.0 benchmark reports remain relevant to 0.6.1. The results below use native
release builds on the same Apple Silicon machine, eight requested threads, one warmup and
the median of three measured runs. **Speedup = comparator wall time / kallistors wall time**.

### Full transcript quantification

Timings include index loading, FASTQ decoding, mapping, abundance inference and normal output.
Ranges cover the measured paired-end and single-end cases; single-end uses read 1.

| Reference and workload | Comparator | Measured speedup |
| --- | --- | ---: |
| GENCODE v50 ALL, real bulk RNA-seq | kallisto 0.52.0 | **6.07–15.40×** |
| Mature + nascent, TT-seq | kallisto 0.52.0 | **2.79–3.54×** |
| GENCODE v50 ALL, real bulk RNA-seq | Salmon 2.8.0 sketch | **2.07–2.42×** |
| GENCODE v50 ALL, real bulk RNA-seq | Salmon 2.8.0 selective alignment | **5.16–5.90×** |
| GENCODE v50 ALL, real bulk RNA-seq | piscem 0.23.0 + piscem-infer 0.7.0 | **2.95–4.03×** |
| Mature + nascent, TT-seq | Salmon 2.8.0 sketch | **1.00–1.37×** |
| Mature + nascent, TT-seq | Salmon 2.8.0 selective alignment | **3.77–4.32×** |
| Mature + nascent, TT-seq | piscem 0.23.0 + piscem-infer 0.7.0 | **1.23–1.31×** |

The GENCODE comparison with original kallisto covers muscle, airway and fibroblast libraries;
the Salmon/piscem ranges cover muscle and fibroblast. TT-seq single-end is effectively tied
with Salmon sketch. Salmon selective alignment uses transcript-only references without
genomic decoys. The piscem times include piscem-infer.

### Index construction

Both tools build from the same transcript target sets, with eight requested threads.
Times exclude downloads, reference preparation and compilation.

| Reference | kallistors, seconds | kallisto 0.52.0, seconds | Speedup |
| --- | ---: | ---: | ---: |
| GENCODE v50 ALL | 11.12 | 121.49 | **10.92×** |
| Mature + nascent | 98.69 | 907.84 | **9.20×** |

### Single-cell processing

The complete 10x PBMC FASTQ-to-gene/UMI-matrix pipeline was **2.26× faster** than original
kallisto with the same downstream bustools processing, and **1.33× faster** than the
piscem + alevin-fry pipeline at eight threads. BUS counts and corrected matrices match
the original-kallisto control exactly.

Performance depends on workload, hardware and build settings. Mapping policies and inference
models differ across tools; full reports include peak memory, numerical agreement, individual
runs and frozen input/tool provenance:

- [Original kallisto: quantification, single-cell processing and index construction](https://github.com/deminden/kallistors/blob/main/docs/kallisto-comparison-results.md)
- [Salmon and piscem: full workflows, native mapping and accuracy](https://github.com/deminden/kallistors/blob/main/docs/salmon-piscem-results.md)
- [Fulgor and Themisto: mapping-only comparisons](https://github.com/deminden/kallistors/blob/main/docs/mapping-comparison-results.md)
- [Benchmark preparation and methodology](https://github.com/deminden/kallistors/blob/main/docs/benchmarks.md)

## Features

- Build and load kallisto v13-compatible transcriptome indexes from transcript FASTA.
- Quantify paired-end and single-end reads, including stranded libraries, optional
  sequence-bias correction and bootstrap samples.
- Write kallisto-style `abundance.tsv`, `run_info.json` and `abundance.h5` output.
- Generate single-cell BUS output for supported technologies, including 10x Genomics,
  with supported FASTQ/BAM input and BAM output workflows.
- Run multi-threaded mapping and inference with automatic bulk match-cache sizing;
  override its memory budget with `--match-cache-size`.

## Install

Requires **Rust 1.97 or newer** to build from source.

```bash
cargo install kallistors --locked
kallistors --help
```

[Prebuilt binaries](https://github.com/deminden/kallistors/releases) are available for Linux
x86_64, macOS x86_64/arm64 and Windows x86_64, with SHA-256 checksums.

## Quick start

```bash
# Build an index from transcript sequences
kallistors index -i transcripts.idx -t 8 transcripts.fa.gz

# Quantify paired-end reads
kallistors quant -i transcripts.idx -o paired_out -t 8 reads_1.fq.gz reads_2.fq.gz

# Quantify single-end reads with a fragment mean and standard deviation
kallistors quant -i transcripts.idx -o single_out -t 8 \
  --single -l 200 -s 20 reads.fq.gz

# Generate BUS output for a 10x Genomics v3 library
kallistors bus -i transcripts.idx -o bus_out -t 8 -x 10XV3 reads_1.fq.gz reads_2.fq.gz
```

BUS output is processed downstream into a gene-count matrix, for example with bustools.
Add `-b 100 --seed 42` for bootstrap samples, or `--bias --transcripts transcripts.fa.gz`
for sequence-specific bias correction during quantification.

## Rust library

```toml
[dependencies]
kallistors = "0.6.1"
```

The library exposes index construction/loading, pseudoalignment and transcript quantification.
See the [Rust API examples](https://github.com/deminden/kallistors/blob/main/docs/development.md#using-the-rust-library)
for paired-end quantification, threading and fragment sampling.

## Compatibility and documentation

Existing kallisto v13 indexes can be reused. Abundance estimates and bootstrap values can
differ from original kallisto; general long-read quantification and fusion detection are
outside the supported surface. Long-read filtering is available for BUS workflows.

- [Quick start and performance comparisons](https://github.com/deminden/kallistors/blob/main/README.md)
- [Command-line guide](https://github.com/deminden/kallistors/blob/main/docs/usage.md)
- [Rust API examples](https://github.com/deminden/kallistors/blob/main/docs/development.md#using-the-rust-library)
- [Supported behavior and compatibility limits](https://github.com/deminden/kallistors/blob/main/docs/compatibility.md)
- [Match-cache memory settings](https://github.com/deminden/kallistors/blob/main/docs/match-cache.md)
- [0.6.1 release notes](https://github.com/deminden/kallistors/blob/main/docs/releases/0.6.1.md)

Licensed under the MIT License; see the included `LICENSE` file.
