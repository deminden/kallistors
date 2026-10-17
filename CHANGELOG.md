# Changelog

## 0.5.1

### Quantification and compatibility

- Retain strand annotations for paired `quant`/`pseudoalign` strand flags and reject empty
  strand-filter results instead of counting them as unstranded matches.
- Let validated backoff probes confirm jumps on reverse-oriented unitigs, preserving SVIL
  read ECs with Rust-built indexes in both scanners.
- Associate abundant k-mers with explicit stored EC-node keys, preserving matches when
  Bifrost hash-table order changes during serialization.
- Apply paired fragment quotas to upstream-compatible 8 MiB sequence buffers, counting both
  mates and carrying quotas across parser chunks. The default goal is 10,000 observations;
  concurrent workers can overshoot. All reads contribute to EC counts.
- Expose `--fragment-length-samples` on `quant` and `pseudoalign`; zero collects all qualifying
  observations. Single-end quantification uses its supplied fragment model.

### Performance and memory

- Exclude repeated transcript windows only after full byte comparison, using shared ends and
  sparse internal anchors to find candidates. Preserve every seam window and all transcript
  metadata. Deduplicate bounded groups before parallel integer radix sorting; append completed
  shards immediately and avoid copying the entire distinct-key table merely to shrink capacity.
- Parse buffered FASTA lines directly and normalize independent sequences in parallel. Ambiguity
  substitutions retain the same RNG and transcript/base order; original lengths and polyA trimming
  are unchanged. Share the exact vector DNA encoder with mapping and minimizer construction.
- Use exact singleton descriptors and small-bucket rejection guards, verifying full keys on
  guard matches. Batch canonical-edge queries, update both endpoints with idempotent CAS and
  publish the common first edge directly. AArch64 cache hints overlap exact range-key reads and positive edge updates;
  larger values retain full-width fallbacks.
- Losslessly pack each completed parallel collection group into two-bit scratch before sorting its words;
  release the original ASCII buffers and restore packed records to transcript order. Reuse worker
  collection buffers and preserve the same group boundaries and complete k-mer set.
- Keep identifier bytes in one arena during topology construction and reconstruct owned names
  exactly after mapping. Name-hash collisions still require complete string equality.
- Compress exact prefix ranks with dense-block/full-width fallbacks and use narrow location,
  state and join records only when their complete domains fit.
- Reuse the requested worker pool across construction stages. Derive continuation eligibility
  from stable rank fields through atomic loads, eliminating a separate degree snapshot. Finalize
  edge records in parallel and reuse their allocation after workers join, with an alignment fallback.
- Keep parallel path preparation separate from serial ordered claims. Use narrow prepared steps
  when ranks fit, materialize accepted branch paths in parallel, and retain bounded residual
  buffers and serial long-path fallback. Radix-sort compact oriented starts; preserve ordered
  joins without a redundant candidate sort and reuse sequence ownership while joining. Fill one
  exact oriented-start table through counted disjoint slices; bounded metadata chunks avoid large
  geometric reallocations.
- Encode whole-unitig transcript runs without repeating their known bounds, retaining complete
  partial-run and large-ID fallbacks. Flatten runs into stable ranges and store each EC block's
  grouped positions contiguously, eliminating per-transcript position allocations. Release graph
  lookup buffers before building payloads. Extend transcript runs through proven nonbranching unitig
  interiors. EC payloads, minimizer work and node serialization retain ordered parallel output.
- Build MPHF presence/collision bitsets per worker, merge collisions exactly and preserve survivor
  order. Keep minimizer keys numeric through hashing and serialize their original byte order;
  use 32-bit survivor/bucket scratch only when all values fit, with full-width fallbacks.
  Stream minimizer hashes through a bounded queue retaining all ties. Serialize tiny
  bitmaps with a native byte oracle, pack graph sequences directly into output, and use the same
  bytes for checksums. Group node buffers and retain portable packing fallbacks.
- Release wholly unused scratch pages on macOS/Linux before later allocations overlap them:
  deduplication tails, completed graph/edge/join buffers and serialization scratch. Retain virtual
  capacity where needed and preserve live values and normal allocator ownership.
- Report non-overlapping index stages, separate MPHF construction from writing and include final
  allocation cleanup in the internal total. Dispose independent allocations in the requested pool.
- Compact singleton minimizer buckets, sequence/EC metadata, loader scratch and EM working arrays.
- Decode graph sequences directly into packed words, avoiding duplicate ordinary/fixed/abundant
  ASCII storage. Share adjacent identical EC payloads, keep narrow block descriptors with wide
  fallbacks and retain abundant-node key associations and candidate order.
- Build eligible minimizer tables in one lookup pass, retaining ordered repeated/special events
  and aligned full-width runtime descriptors. Release redundant rank storage and use bounded,
  exact gap-encoded bitmap scratch with a raw fallback.
- Share complete paired-worker EC keys through collision-checked sharded tables and keep worker
  counters independent. Rendezvous work dispatch bounds queued read storage while preserving
  upstream sampling buffers and all reads.
- Release graph storage before EM, consume unneeded EC input after preparation and narrow
  transcript IDs in place. Bootstrap input remains available when needed.
- Defer ordinary transcript-name allocation until output after scanning all names for shades.
  Retain the same open index file, complete target IDs and eager compatibility loaders.
- Use exact grouped hashing, precomputed MPHF divisors, compact caches and packed reference
  windows. AArch64 NEON encoding has a portable fallback; wide values and minimizer ties retain
  complete fallback handling.
- Group independent EM components in stable order, fuse their denominator/update work and
  skip only bitwise fixed components. Preserve mandatory iterations, stopping rules, final
  threshold updates and every transcript's arithmetic order; cache entirely fixed branches without
  changing component partitions. Bias correction retains its general path.
- Reuse scan scratch and exact mode-aware matching caches. Filter ordinary minimizer buckets
  against every possible read anchor, bypassing conflicting/special buckets. Use compact MPHF
  rank records with full-width fallbacks and exact reciprocal reduction. Bounded AArch64 cache
  hints accelerate ordered loader writes without changing candidate order.
- Fold paired hits directly into ECs after fragment sampling finishes when placement metadata
  is unused. Reject an already hard-rejected pair before scanning its second mate.
- Bound paired input queues and release blocked producers after input errors.
- Format TSV rows in bounded ordered chunks, overlap TSV/HDF5 output and bound HDF5 chunks
  to 16,384 elements with DEFLATE level 6.

### Structure and validation

- Support checksum-verified GENCODE v50 ALL and protein-coding reference preparation.
- Separate pseudoalignment diagnostics, quantification output and BUS genome/BAM output into
  focused modules, preserving public entry points.
- Cover output formatting, encoding, caches, loader associations, fragment sampling, input
  errors, EM arithmetic and HDF5 chunk tails with regression tests.
- Provide manifest-driven ENA/fasterp preparation, GENCODE cDNA and published Pachter Lab
  mature-plus-nascent reference preparation, and a shared quantification runner with complete EC preflight, phase timing, memory measurements and resumable
  result records. Native HDF5 validation includes fragment statistics; reporting exports JSON/TSV.
  Benchmark cancellation also terminates its timed process group. See the [benchmark guide](docs/benchmarks.md).
- Use `target-cpu=native` for source builds and `target-cpu=generic` in CI and release binaries.

### API and requirements

- Rust callers using exhaustive `PseudoalignOptions` literals must include
  `fragment_length_samples`; `..PseudoalignOptions::default()` selects the 10,000 goal.
- Existing v13 indexes and output dataset schemas remain supported.
- Add `strength_reduce 0.2.4` for exact precomputed divisors. The minimum Rust version is 1.97.

See [release notes](docs/releases/0.5.1.md), [sampling behavior](docs/fragment_length_sampling.md)
and [compatibility limits](docs/compatibility.md) for usage details.
