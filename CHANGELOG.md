# Changelog

## 0.6.0

### Correctness and compatibility

- Prevent Windows MSVC debug CLI startup from overflowing the default main-thread stack.
  Reserve 8 MiB for the executable without changing worker stacks or library callers.

- Load transcript orientation annotations for explicit BUS strand filtering, including
  single-cDNA FASTQ, paired technologies and BAM input. Preserve the compact unstranded path.
  Apply the same annotation requirement in read diagnostics, and correct reverse-complement
  lookup in the optional direct-kmer diagnostic scanner.
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

- Share requested bootstrap workers between independent samples and large within-sample
  components, including a single sample. Retain the original seed sequence, result order,
  full-target prior and exact ordered equations. Reuse compact IDs in large scalar samples,
  omit unobservable convergence counts and sequence TSV/bootstrap work within the worker budget.
- Sample exact SplitMix chunks in the existing pool with a serial replay on any rejected
  candidate. Serial callers can share a bounded exact inverse-CDF table and batched draws;
  larger libraries retain the compact search directory. Preserve every draw and RNG state.
- Accelerate large-reference loading with batched exact lookups, compact fingerprint and
  annotation storage, and bounded ordered parallel positional decoding. Reuse complete
  immutable EC-pair intersections within mapping invocations, retaining exact multiset
  semantics, serial fallbacks, index compatibility and worker limits.
- Accelerate single-cDNA BUS counting with incremental exact minimizer hashes,
  bounded conservative candidate directories, compact rank words and inline read-cache keys.
  Preserve ordered verification and rescue, complete cache-key checks and index serialization.
- Extend ordered parallel BUS counting to compound/synthetic barcodes, tagged
  layouts, paired FASTQ and ordinary BAM input. Retain complete paired BUS
  fragment histograms and omit unused diagnostic traces and BAM quality arrays.
- Keep ordinary stranded ECs borrowed through orientation filtering, use packed
  positional graphs for stranded counting, and switch to lighter counting when
  bias sampling is complete. Preserve comprehensive strand filters and exact
  fragment/bias sampling semantics.
- Decode gzip with directly buffered Rust `zlib-rs`, retaining CRC, multi-member and truncation
  checks. Use synchronous sequence-only parsing for single-thread paired counting.
- Load packed graph/EC data with compact minimizer metadata and deferred transcript names.
  Reuse exact query hashes, conservative all-anchor fingerprints and collision-checked caches.
  Borrow unchanged ECs through paired counting after fragment sampling.
- Reduce single-end memory with packed positional loading, lossless compact position groups,
  shared positional/mapping EC lists and block bounds, tightly owned block arrays, bounded
  prefetch and shared worker EC keys. Derive strand summaries from raw positions and select
  narrower descriptors only when exact bounds prove they fit. Preserve the supplied fragment
  model, complete coordinate range, exact EC membership, strand-aware loading and bias correction.
- Dispatch compact position packing once per block. Borrow annotation-block ranges and stop
  fragment-position searches at the first valid coordinate. Filter owned ECs in place and store
  shared immutable EC keys without spare capacity.
- Fold ordinary single-end ECs while scanning, retaining first/best placement without a full
  hit list. Keep unchanged ECs borrowed through fragment filtering and copy new counted keys.
  Use smaller single-end batches to bound live input storage and improve worker balance.
- Bound queued reads and share complete worker EC keys. Release graph storage before EM.
  Adapt EM scheduling to graph work: keep small components local and share dominant components
  across the same worker pool. Preserve scalar arithmetic order, stopping rules and the final
  update. Share compact components and exact witness scheduling with eligible bias-corrected
  runs, resetting certificates at model refreshes. Skip repeated bitwise fixed work and incoming
  sums whose positive-zero result is proved from finite operands; retain every denominator.
- Compute independent bias effective-length rescales through the existing requested pool,
  sharing one equation implementation across serial and parallel scheduling. Keep histograms,
  normalization and each transcript's contribution order unchanged.
- Release completed mapping storage before loading bias sequences. Borrow requested header
  keys and assign sequence buffers directly to output slots, preserving duplicate/missing-name
  behavior and complete UTF-8 validation without retaining a second name table.
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
- Format TSV in bounded ordered chunks and overlap TSV/HDF5 output when threaded. Use
  lower DEFLATE effort for HDF5 output, trading larger compressed files for less writing work
  while preserving every native dataset, type, shape and numeric value.
- Process ordinary single-cDNA BUS layouts with bounded packed input, concurrent decoding,
  ordered parallel counting and the shared compact mapping index. Borrow unchanged ECs and
  omit unused placement traces. Overlap mapping with ordered BUS output and reuse full-key
  cDNA answers across workers and lanes using bounded storage. Stream transcript names and
  allocate length statistics only for long reads. Preserve extended-mode processing, exact
  record counts and read-limit boundaries even when unread trailing FASTQ data is malformed.

### Maintenance and requirements

- Separate graph construction, packed index loading, fast matching, EM work, fragment/bias
  models, bootstraps and BUS configuration/input into focused private modules. Keep public APIs,
  arithmetic order, scanner candidate order and storage layouts unchanged. Share single-end
  fragment/strand filtering across serial and batched counting.
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
- Remove inactive workspace lint tables and a stale dead-code annotation. Retain the legacy
  `--kallisto-sparse-hits` option as an explicit compatibility no-op. Strict workspace Clippy
  checks remain part of the development workflow.
- Rust callers with exhaustive `PseudoalignOptions` literals must include
  `fragment_length_samples` and `match_cache_size`; `..PseudoalignOptions::default()` selects
  the sampling goal and automatic cache budget.
- Keep ordinary CLI help focused while retaining developer diagnostics by name. Add a compact
  package README, inherited author metadata, portable mapping commands and automated
  documentation/script checks. Expand CI to supported operating systems and the minimum Rust
  version, and exercise release binaries before upload.
- Retire stale benchmark publishers and reject duplicate, missing or undeclared trial records
  before a dependent comparison starts. Preserve completed benchmark evidence separately.

See the [release notes](docs/releases/0.6.0.md), [current benchmarks](docs/benchmarks.md),
[sampling behavior](docs/fragment_length_sampling.md) and [compatibility limits](docs/compatibility.md).
