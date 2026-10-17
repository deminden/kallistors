# Changelog

## 0.5.1

### Quantification and compatibility

- Associate abundant k-mers with explicit stored EC-node keys, preserving matches when
  Bifrost hash-table order changes during serialization.
- Apply paired fragment quotas to upstream-compatible 8 MiB sequence buffers, counting both
  mates and carrying quotas across parser chunks. The default goal is 10,000 observations;
  concurrent workers can overshoot. All reads contribute to EC counts.
- Expose `--fragment-length-samples` on `quant` and `pseudoalign`; zero collects all qualifying
  observations. Single-end quantification uses its supplied fragment model.

### Performance and memory

- Compact singleton minimizer buckets, sequence/EC metadata, loader scratch and EM working arrays.
- Reuse loader bucket lookups and decode packed sequence/EC data without unnecessary allocations.
- Use exact grouped hashing, precomputed MPHF divisors, compact caches and packed reference
  windows. AArch64 NEON encoding has a portable fallback; wide values and minimizer ties retain
  complete fallback handling.
- Parallelize independent transcript updates and convergence checks while retaining each
  transcript's arithmetic order.
- Bound paired input queues and release blocked producers after input errors.
- Overlap TSV/HDF5 output and bound HDF5 chunks to 16,384 elements with DEFLATE level 6.

### Structure and validation

- Separate pseudoalignment diagnostics, quantification output and BUS genome/BAM output into
  focused modules, preserving public entry points.
- Cover output formatting, encoding, caches, loader associations, fragment sampling, input
  errors, EM arithmetic and HDF5 chunk tails with regression tests.
- Provide manifest-driven ENA/fasterp preparation, GENCODE cDNA and published Pachter Lab
  mature-plus-nascent reference preparation, and a shared quantification runner with complete EC preflight, phase timing, memory measurements and resumable
  result records. Native HDF5 validation includes fragment statistics; reporting exports JSON/TSV.
  See the [benchmark guide](docs/benchmarks.md).
- Use `target-cpu=native` for source builds and `target-cpu=generic` in CI and release binaries.

### API and requirements

- Rust callers using exhaustive `PseudoalignOptions` literals must include
  `fragment_length_samples`; `..PseudoalignOptions::default()` selects the 10,000 goal.
- Existing v13 indexes and output dataset schemas remain supported.
- Add `strength_reduce 0.2.4` for exact precomputed divisors. The minimum Rust version is 1.97.

See [release notes](docs/releases/0.5.1.md), [sampling behavior](docs/fragment_length_sampling.md)
and [compatibility limits](docs/compatibility.md) for usage details.
