# Changelog

## Unreleased

## 0.5.0 — 2026-10-14

### Performance

- Eight-thread full paired quantification measured **2.72x faster than original kallisto v0.52.0**
  (47.673 s to 17.532 s, 63.22% less elapsed time) on Apple Silicon with fastp-trimmed
  SRR13638690 and GENCODE v50 ALL transcripts, including index load and default HDF5 output.
  Peak RSS is **2.85 GiB upstream versus 3.01 GiB for kallistors**, or **5.7% higher**.
  This is one measured pair with warmed inputs; background CPU activity was recorded.

### Changes

- Remove the per-read scan of the entire transcript on-list in both scanners. Dummy transcript
  cardinality is calculated only when a synthetic D-list hit needs it.
- Add a D-list/on-list regression for both scanners and `scripts/compare_ec_counts.py` for exact
  full EC count comparison independent of EC numbering. All six documented real-read maps
  remain exact against upstream.
- Add paired-only benchmark selection with `--modes` and configurable measured repetitions with
  `--runs`; the defaults remain both read modes and three measured runs after one warmup.
- Document profiling, the main bottleneck, and validation in the
  [paired performance report](docs/paired_performance_fix_2026-10-14.md).
- Pack unitig bases and reuse exact read/minimizer encoding, preserving candidate and fallback order.
- Preserve buffered index reads and fill unitig minimizer buckets with ordered parallel lookup.
- Derive loader minimizers from packed bases instead of allocating full per-position tables.
  Use flat EC bounds for jumps and avoid duplicate EC graphs in ordinary quant/pseudoalign;
  retain positional, strand and diagnostic metadata when needed. Reuse count buffers for prefix
  offsets and release loading scratch after its last use.
- Restrict repeated EM work to active transcripts; use the requested quant thread count for
  parallel denominators and per-transcript updates with original contribution order.
- Reuse positive primary matches and worker MPHF scratch with explicit index/k/g cache ownership;
  retain full minimizer tie lists when they exceed inline cache capacity.
- Reduce packed FASTQ discard copies and use equality-checked numeric EC hashing.
- Add alternating whole-process variant benchmarks and numerical, encoding, cache, loader-order
  and FASTQ regressions; add portable peak RSS measurement with `--memory`. See the [performance report](docs/performance_0.5.0.md).

See [0.5.0 release notes](docs/releases/0.5.0.md) for benchmark provenance, compatibility
evidence, the threaded library API and migration guidance.

## 0.4.2 — 2026-10-14

### Fixed

- Paired-end quantification uses the empirical fragment histogram for transcript-specific
  effective lengths, including the means used by bias correction and the effective lengths used
  by bootstraps. This corrects substantial short-transcript TPM differences.
- Both read scanners complete the next valid k-mer's fallback probe before resuming jumps.
  Two previously ambiguous TTN reads now receive the upstream singleton EC.
- Genome BAM indexing handles both BAI and CSI results from the updated noodles API.
- HDF5 fixed ASCII types use the updated hdf5-pure public API.

### Added

- Fragment-model regression tests and small TTN/AGAP12P read fixtures; the full GENCODE v50
  scanner regression is opt-in and covers both scanner implementations.
- Reproducible ENA downloads with MD5/gzip checks, fastp trimming, GENCODE v50 ALL-transcript
  comparisons, and an optional R 4.6.1 / fastp 1.4.0 environment.
- Repeated index and paired/single-end benchmarks, retained final quantification outputs,
  and abundance audit summaries that distinguish output rounding from quantification changes.
- A tracked diagnostic patch for upstream kallisto v0.52.0 and reports covering the mismatches,
  their fixes, and the measured performance.

### Changed

- Refreshed Rust dependencies and `Cargo.lock`; modernized fixed-size slice iteration for current
  Clippy. Rust 1.97 remains the declared minimum.
- Updated the README and library example to describe empirical paired quantification, exact EC
  evidence, residual fragment-sampling differences, and workload-dependent speed.
- Release automation uses version-specific notes, locked CI dependencies, and includes the changelog and documentation
  in binary archives. The crate package includes the license text.

See [0.4.2 release notes](docs/releases/0.4.2.md) for dependency versions, migration guidance,
validation evidence, and publication commands.

## 0.4.1 — 2026-07-22

- Required Rust 1.97 / edition 2024 and refreshed dependencies.
- Adapted HDF5 handling, hardened BUS parity coverage, and refreshed release documentation.

## 0.4.0

- Introduced native binary release packaging and refreshed dependencies.
- Added and hardened BUS processing, including the modular input and technology implementation.

Older implementation history is available in Git commits and the compatibility roadmap.
