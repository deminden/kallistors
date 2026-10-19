# Changelog

## 0.5.1

### Correctness and compatibility

- Advance ambiguous jump/midpoint probes to the next valid k-mer as original kallisto does,
  preserving mature/nascent EC membership and paired alignment.
- Preserve paired strand filtering, reverse-unitig backoff matches and abundant-node EC
  associations when serialized Bifrost hash-table order changes.
- Reject empty ECs in BUS output when incompatible transcript hits collapse the read
  intersection. Preserve the rejection trace without emitting a mapping or inflating counts.
- Match original kallisto's paired fragment-sampling goal of 10,000 observations and 8 MiB
  combined-sequence buffers. Workers may overshoot; every read contributes to EC counts.
  `--fragment-length-samples 0` collects all qualifying observations.
- Preserve v13 index compatibility and kallisto-style output schemas. Bound HDF5 chunks to
  16,384 elements and validate complete output with native HDF5.

### Speed and memory

- Decode gzip with directly buffered Rust `zlib-rs`, retaining CRC, multi-member and truncation
  checks. Use synchronous sequence-only parsing for single-thread paired counting.
- Load packed graph/EC data with compact minimizer metadata and deferred transcript names.
  Reuse exact query hashes, conservative all-anchor fingerprints and collision-checked caches.
  Borrow unchanged ECs through paired counting after fragment sampling.
- Bound queued reads and share complete worker EC keys. Release graph storage before EM.
  Adapt EM scheduling to graph work: keep small components local and share dominant components
  across the same worker pool. Preserve scalar arithmetic order, stopping rules and the final
  update; skip only bitwise fixed work.
- Build indexes with byte-proven repeat-window exclusion, bounded packed collections, parallel
  radix sorting and batched exact edge lookups. Compact graph/EC scratch and retire completed
  allocations before later stages. Preserve ordered graph claims, all minimizer ties and
  deterministic serialization.
- Adapt large-reference construction with exact four-byte relative edge ranks, five-byte keys
  packed before topology, and a prefix-owned conservative rejection filter. Verify every admitted
  key exactly; preserve existing small-reference paths and full-domain fallbacks.
- Replace dense large-reference walk locations with exact full-word landmarks, including internal
  transcript starts. Batch graph probes, compact join records, and retire completed bitmap/MPHF
  scratch before later allocations. Preserve ordered claims and byte-identical indexes.
- Use one requested worker pool throughout construction. Report exclusive build stages,
  separate MPHF construction from writing and include final cleanup in the internal total.
- Format TSV in bounded ordered chunks and overlap TSV/HDF5 output when threaded.
- Process ordinary single-cDNA BUS layouts with bounded packed input, concurrent decoding,
  ordered parallel counting and the shared compact mapping index. Borrow unchanged ECs and
  omit unused placement traces. Overlap mapping with ordered BUS output and reuse full-key
  cDNA answers across workers and lanes using bounded storage. Stream transcript names and
  allocate length statistics only for long reads. Preserve extended-mode processing, exact
  record counts and read-limit boundaries even when unread trailing FASTQ data is malformed.

### Maintenance and requirements

- Extract diagnostics, quantification writers and BUS genome/BAM output into focused modules.
  Consolidate overlapping tests and benchmark scripts; retain distinct regressions and exact
  EC checks. Documentation describes current behavior and current release measurements.
- Provide checksum-verified ENA retrieval, fastp trimming, GENCODE v50 ALL/protein-coding
  preparation and official Pachter Lab mature-plus-nascent references. Shared benchmark tools
  record complete runs, phases, peak memory, EC parity, native HDF5 and abundance checks.
- Add complete 10x PBMC v3 retrieval and barcode-preserving transcript trimming. Compare
  complete BUS counts and shared bustools-corrected gene/UMI matrices independently of timing.
- Require Rust 1.97 or newer. Source builds use `target-cpu=native`; CI and release binaries
  use `target-cpu=generic`. Add `strength_reduce 0.2.4` for exact precomputed divisors.
- Rust callers with exhaustive `PseudoalignOptions` literals must include
  `fragment_length_samples`; `..PseudoalignOptions::default()` selects the sampling goal.

See the [release notes](docs/releases/0.5.1.md), [current benchmarks](docs/benchmarks.md),
[sampling behavior](docs/fragment_length_sampling.md) and [compatibility limits](docs/compatibility.md).
