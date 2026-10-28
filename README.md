# kallistors

A high-performance Rust implementation of kallisto-style pseudoalignment and RNA-seq transcript quantification.

## At a glance

- **Fastest measured full quantification**, compared with Salmon,
  piscem + piscem-infer and original kallisto. Nascent single-end at eight threads is
  effectively tied with Salmon sketch.
- Builds kallisto v13-compatible transcriptome indexes **9.2–10.9× faster** than original
  kallisto at eight threads in the benchmarks below. Also loads existing compatible indexes.
- Quantifies paired-end and single-end reads **6.1–15.4× faster** on the real GENCODE datasets
  and **2.8–3.5× faster** on mature-plus-nascent than original kallisto at eight threads.
  Writes kallisto-style `abundance.tsv`, `run_info.json`, and default `abundance.h5` output.
- Supports single-cell BUS workflows, including supported BAM input and output modes.
- Uses memory-efficient storage, exact-match caching and efficient parallel processing.

## Full quantification compared with other tools

Latest available measurements on the same Apple Silicon machine, using matched prepared
reads and reference target sets. Read results use one warmup and
the median of three runs. Ratios are **comparator time / kallistors time**; above 1 means
kallistors finished sooner. Peak RSS is shown as **kallistors / comparator** in GiB.

GENCODE quantification ranges cover the muscle and fibroblast libraries; mature-plus-nascent
quantification uses TT-seq. PE means paired-end; SE uses R1 only. Each row includes the complete
quantification or FASTQ-to-matrix workflow.

| Workload | Reference | Comparator | Time ratio · 1 thread | Time ratio · 8 threads | 8-thread peak RSS, GiB · current/comparator |
| --- | --- | --- | --- | --- | --- |
| PE quantification | GENCODE v50 | Salmon 2.8.0 sketch | 2.66–2.72× | 2.42× | 2.43–2.60 / 1.03–1.39 |
| PE quantification | GENCODE v50 | Salmon 2.8.0 selective alignment | 6.03–7.49× | 5.16–5.90× | 2.43–2.60 / 2.49–2.59 |
| PE quantification | GENCODE v50 | piscem 0.23.0 + piscem-infer 0.7.0 | 1.83–1.90× | 2.95–3.15× | 2.43–2.60 / 0.76–0.79 |
| SE quantification | GENCODE v50 | Salmon 2.8.0 sketch | 2.31–2.58× | 2.07–2.09× | 2.33–2.50 / 1.00–1.09 |
| SE quantification | GENCODE v50 | Salmon 2.8.0 selective alignment | 6.55–7.61× | 5.40–5.42× | 2.33–2.50 / 2.43–2.51 |
| SE quantification | GENCODE v50 | piscem 0.23.0 + piscem-infer 0.7.0 | 1.55–1.79× | 3.06–4.03× | 2.33–2.50 / 0.75–0.77 |
| PE quantification | Mature + nascent | Salmon 2.8.0 sketch | 1.86× | 1.37× | 9.21 / 3.30 |
| PE quantification | Mature + nascent | Salmon 2.8.0 selective alignment | 5.40× | 3.77× | 9.21 / 5.61 |
| PE quantification | Mature + nascent | piscem 0.23.0 + piscem-infer 0.7.0 | 1.37× | 1.31× | 9.21 / 3.22 |
| SE quantification | Mature + nascent | Salmon 2.8.0 sketch | 1.49× | 1.00× | 10.26 / 3.26 |
| SE quantification | Mature + nascent | Salmon 2.8.0 selective alignment | 6.41× | 4.32× | 10.26 / 5.57 |
| SE quantification | Mature + nascent | piscem 0.23.0 + piscem-infer 0.7.0 | 1.11× | 1.23× | 10.26 / 3.19 |
| PBMC FASTQ to matrix | GENCODE v50 | piscem 0.23.0 + alevin-fry 0.18.3 | — | 1.33× | 2.03 / 0.76 |

Full quantification includes mapping and native abundance inference/output. Salmon selective
alignment uses the transcript-only reference without genomic decoys; piscem bulk quantification
includes piscem-infer. PBMC matrix processing includes bustools for kallistors and alevin-fry
for piscem; the complete alevin-fry pipeline was measured at eight threads only.
Mapping policies, inference models and native outputs differ across tools.

See the [full Salmon/piscem comparison](docs/salmon-piscem-results.md#full-quantification) for
individual datasets, CPU/memory measurements, settings and reproducibility records.
Mapping-only results have separate tables in the
[Fulgor/Themisto report](docs/mapping-comparison-results.md#mapping-command-medians) and the
[Salmon/piscem report](docs/salmon-piscem-results.md#native-mapping).

## Comparison with original kallisto

Compared with original **kallisto 0.52.0**, using the same prepared reads and original-built
indexes. Times below are **kallistors / original** in seconds; RSS is peak memory usage.
Bulk timings cover the full quantification command, including loading and output.
PBMC includes the same downstream bustools matrix pipeline.

| Workload | 1 thread, s · current / original | 8 threads, s · current / original | 8-thread time ratio · original/current | 8-thread peak RSS, GiB · current / original |
| --- | --- | --- | --- | --- |
| Muscle · PE | 28.45 / 101.37 | 6.33 / 46.11 | 7.29× | 2.43 / 2.61 |
| Muscle · SE | 19.29 / 164.56 | 4.86 / 57.61 | 11.84× | 2.33 / 2.62 |
| Airway · PE | 102.20 / 315.29 | 22.15 / 180.99 | 8.17× | 2.42 / 2.89 |
| Airway · SE | 73.79 / 1011.02 | 17.11 / 263.42 | 15.40× | 2.50 / 2.87 |
| Fibroblast · PE | 166.32 / 536.59 | 34.20 / 207.47 | 6.07× | 2.60 / 2.84 |
| Fibroblast · SE | 99.79 / 1267.38 | 21.74 / 301.25 | 13.86× | 2.50 / 2.79 |
| TT-seq · mature + nascent · PE | 147.23 / 343.07 | 34.88 / 97.17 | 2.79× | 9.21 / 16.08 |
| TT-seq · mature + nascent · SE | 100.71 / 344.16 | 27.13 / 96.04 | 3.54× | 10.26 / 16.00 |
| PBMC · FASTQ to matrix | 407.15 / 940.30 | 80.81 / 182.37 | 2.26× | 2.03 / 2.59 |

Counts for each transcript group (EC) match in every ordinary bulk condition. Full PBMC BUS
counts and corrected gene/UMI matrices also match exactly. Abundance estimates can still
differ because of fragment sampling and statistical models. Matching output does not establish
biological accuracy. Timings depend on the workload, hardware and build settings.

The [full original-kallisto comparison](docs/kallisto-comparison-results.md) includes all
datasets, simulation accuracy, bootstrap modes, individual timings, CPU/RSS and numerical
agreement. Completed [sequence-bias results](docs/kallisto-comparison-results.md#completed-sequence-bias-observations)
have their own section because the abundance estimates differ materially; their timing ratios
are excluded from the table above. See the [benchmark guide](docs/benchmarks.md) for preparation
and methodology, and the [0.6.0 notes](docs/releases/0.6.0.md) for release changes.

## Index construction

Direct comparison of completed index construction on the same Apple Silicon machine,
using the same transcript target sets and eight requested threads. **Build-time ratio = tool
time / kallistors time**; for example, 4× means kallistors built the index four times faster.
Every successful entry is the median of three builds after one warmup.

| Reference | Tool | Build time, s | Build-time ratio | Peak RSS, GiB | Index size, GiB |
| --- | --- | ---: | ---: | ---: | ---: |
| GENCODE v50 | kallistors | 11.12 | 1.00× | 4.87 | 1.23 |
| GENCODE v50 | kallisto 0.52.0 | 121.49 | 10.92× | 7.69 | 1.19 |
| GENCODE v50 | Salmon 2.8.0 | 52.18 | 4.69× | 10.11 | 1.96 |
| GENCODE v50 | piscem 0.23.0 | 47.12 | 4.24× | 5.14 | 0.57 |
| GENCODE v50 | Fulgor 4.2.0 | 46.28 | 4.16× | 14.15 | 0.35 |
| GENCODE v50 | Themisto 3.2.2 | 201.08 | 18.08× | 8.97 | 0.48 |
| Mature + nascent | kallistors | 98.69 | 1.00× | 18.80 | 5.94 |
| Mature + nascent | kallisto 0.52.0 | 907.84 | 9.20× | 20.15 | 6.13 |
| Mature + nascent | Salmon 2.8.0 | 255.78 | 2.59× | 15.42 | 5.31 |
| Mature + nascent | piscem 0.23.0 | 240.47 | 2.44× | 11.94 | 3.08 |
| Mature + nascent | Fulgor 4.2.0 | Unavailable | — | — | — |
| Mature + nascent | Themisto 3.2.2 | 1504.83 | 15.25× | 12.37 | 3.23 |

Times exclude downloads, reference preparation and tool compilation; index sizes cover the
saved indexes. All completed runs are retained, including background-activity flags.
Fulgor's nascent build failed.
See the [original build details](docs/kallisto-comparison-results.md#index-construction),
[Salmon/piscem build details](docs/salmon-piscem-results.md#index-construction), and
[Fulgor/Themisto setup details](docs/mapping-comparison-results.md#setup-at-eight-requested-threads).

## Compatibility notes

- **Indexes:** build and load kallisto v13-compatible transcriptome indexes; amino-acid
  BUS indexes are also supported.
- **Quantification:** paired-end estimates fragment lengths from reads; single-end requires
  a mean and standard deviation. Optional sequence-bias correction requires the transcript FASTA.
- **Outputs:** kallisto-style abundance files, bootstrap samples and supported BUS/BAM workflows.
  Abundance estimates and bootstrap values can differ from original kallisto.
- **Limits:** general long-read quantification and fusion detection are not implemented.
  Long-read filtering is available for BUS workflows.

See [compatibility details](docs/compatibility.md) for supported options and output differences,
and the [numerical comparisons](docs/kallisto-comparison-results.md#native-abundance-agreement)
for measured agreement with original kallisto.

## Usage

### As a Binary

Download [prebuilt binaries](https://github.com/deminden/kallistors/releases) for Linux x86_64,
macOS x86_64/arm64 or Windows x86_64. Releases include SHA-256 checksums. To install from source,
use Rust 1.97 or newer:

```bash
cargo install --locked --git https://github.com/deminden/kallistors kallistors
```

Common commands:

```bash
# Build a transcriptome index
kallistors index -i transcripts.idx -t 8 transcripts.fa.gz

# Quantify paired-end reads
kallistors quant -i transcripts.idx -o paired_out -t 8 reads_1.fq.gz reads_2.fq.gz

# Quantify single-end reads: supply fragment mean and standard deviation in bases
kallistors quant -i transcripts.idx -o single_out -t 8 \
    --single -l 200 -s 20 reads.fq.gz

# Generate BUS output for 10x Genomics v3
kallistors bus -i transcripts.idx -o bus_out -t 8 -x 10XV3 reads_1.fq.gz reads_2.fq.gz
```

- Quantification writes `abundance.tsv`, `run_info.json` and `abundance.h5` by default.
  Use `--plaintext` to omit HDF5 output.
- Add `-b 100 --seed 42` to generate 100 bootstrap samples.
- Enable sequence-bias correction with `--bias --transcripts transcripts.fa.gz`.
- Bulk mapping chooses its cache size automatically. Override it with, for example,
  `--match-cache-size 4GiB`; see [memory settings](docs/match-cache.md).
- BUS output needs downstream processing to produce a gene-count matrix; see the
  [complete single-cell workflow](docs/benchmarks.md#single-cell-complete-10x-pbmc-1k-v3).

The [command-line guide](docs/usage.md) covers additional BUS technologies, BAM input/output,
batch input, amino-acid indexes, long-read BUS filtering and output options.

### As a Crate

Use kallistors from Rust to build indexes, map reads and estimate transcript abundance.
Add to `Cargo.toml`:

```toml
[dependencies]
kallistors = { git = "https://github.com/deminden/kallistors" }
```

See [Rust examples and configuration](docs/development.md#using-the-rust-library) for index
construction, paired-end quantification, threading and fragment sampling.

## Documentation

- [Command-line guide](docs/usage.md) and [compatibility details](docs/compatibility.md).
- [Benchmark methodology](docs/benchmarks.md) and full results for
  [Salmon/piscem](docs/salmon-piscem-results.md), [original kallisto](docs/kallisto-comparison-results.md)
  and [mapping-only tools](docs/mapping-comparison-results.md).
- [Development guide](docs/development.md): Rust API examples, testing, debugging and source layout.
- [Release notes](docs/releases/0.6.0.md).

## Contributing

Open an issue to discuss ideas, report bugs, or request features.

Pull requests are encouraged, especially for:

- performance improvements
- correctness / numerical stability fixes
- additional tests (including cross-validation vs original)
- documentation, examples, and benchmarking

Development workflow notes are in [docs/development.md](docs/development.md).

## License

This project is licensed under the MIT License. See `LICENSE`.
