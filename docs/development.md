# Development notes

This page collects developer-facing workflows for checks, parity debugging, and
abundance investigation. For the patched upstream `kallisto` build and its extra
tracing flags, see [kallisto_debug_build.md](kallisto_debug_build.md).

## Toolchain

Development requires Rust 1.97 or newer. The workspace uses Rust edition 2024 and declares the
minimum supported Rust version in the root `Cargo.toml`. Source builds use
`target-cpu=native` from `.cargo/config.toml`; CI overrides `RUSTFLAGS` with
`-C target-cpu=generic` for portable release artifacts and validation.

## Required checks

Before committing changes, run:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

## Test strategy and prerequisites

The default suite is deterministic and self-contained. Small mapping fixtures always use the
Rust index builder, regardless of installed tools. They include shared transcript sequence,
unique sequence and unmapped reads; threaded mapping cases span multiple work batches.

Integration tests are grouped into six executable targets to limit repeated compilation and
linking. Run just the affected target during development:

| Target | Contract checked |
| --- | --- |
| `index_builder` | Index structure, normalization, empty graphs and repeated names |
| `mapping` | Complete EC maps against the naive scanner, strand behavior, fragment sampling and multi-batch threading |
| `quantification` | EM against an independent dense reference, serial/threaded bit parity, count conservation, active bias correction and TSV/HDF5/bootstrap values |
| `bus` | Technology layouts, input validation, BUS records and decoded BAM/index outputs |
| `upstream` | Opt-in comparisons with original kallisto, including both builders crossed with both quantifiers |
| `full_index` | Opt-in GENCODE v50 scanner/abundant-node/fallback regressions and official NAC ambiguous-jump parity |

For example, `cargo test --locked -p kallistors --test mapping` runs only mapping integration
checks. Individual modules and cases can be selected with the usual test-name filter.

Keep tests tied to observable behavior or a documented algorithm/format boundary. JSON and TSV
checks parse schemas and values without fixing whitespace, JSON field order or version strings.
Parallel tests must prove their fixture reaches multiple batches; bias tests must show a real
correction. Comparisons against a shared implementation also need independent expected values
or invariants. A test that compares a function to itself adds no protection. Deadlock regressions
run in subprocesses with a watchdog that terminates and reaps a stalled process.

Small focused unit tests remain beside the implementation for packing/width boundaries,
minimizer selection, encoding and sampling quotas. Integration helpers live in
`crates/kallistors/tests/helpers/mod.rs`. Avoid duplicate fixtures and success-only checks when
an existing case can verify the same option together with its output.

BUS validation, input inference, barcode geometry and corrected-tag variants use labelled case
tables with shared tiny fixtures. Keep separate cases for distinct parser branches and error
conditions; every matrix assertion should identify its case. Reuse output checks for verbose
reporting and basic record decoding. The paired batch regression checks complete fragment
histograms with facing mates, including a different expected fragment length for each batch.

Tests that directly compare against original kallisto are ignored by default. With kallisto
v0.52.0 on `PATH`, enable them with:

```bash
cargo test --locked -p kallistors --test upstream --test bus -- --ignored
```

Full-index regressions require prepared GENCODE v50 indexes and are also ignored by default.
The reverse-unitig SVIL case additionally uses `data/srr-small/gencode.v50.kallistors.idx`,
built with `kallistors index` from the same ALL transcript FASTA; `KALLISTORS_TEST_RUST_INDEX`
can select its path. The other full-index cases use the original builder's index.
The paired NAC regression also requires the original-built Pachter Lab mature-plus-nascent index
at `data/official-human-nac/human_nac.mature-plus-nascent.idx`; `KALLISTORS_TEST_NAC_INDEX`
can select its path. It checks two real TT-seq pairs in both scanners and at one/eight threads.
Enabling an external test without its prerequisite fails clearly. Ordinary tests do not download
data or measure speed: use the [benchmark workflow](benchmarks.md) to evaluate performance.
HDF5 readback through the Rust library checks values; native-library validation below supplies
an independent format check. Passing synthetic tests does not establish general real-data or
biological accuracy.

## Native HDF5 output validation

For changes to HDF5 output, also read a generated quantification directory with the native
HDF5 library. This catches format errors that a writer and reader from the same implementation
can both accept. Install `h5py` and `numpy` in a Python environment, then run:

```bash
python scripts/check_hdf5_output.py path/to/quant-output --json hdf5-readback.json
```

The check reads every dataset, compares all transcript IDs, estimated counts, effective lengths
and target lengths with `abundance.tsv`, and validates processed-read and bootstrap metadata.
The Rust `quantification` integration tests include a 70,003-target case with a partial final chunk.

## Release workflow

The crate version lives in `crates/kallistors/Cargo.toml`. Version-specific release notes
live in `docs/releases/`; the binary publication workflow reads the file matching its tag.

For a release:

1. Bump the crate version and refresh its package entry in `Cargo.lock`. Update `CHANGELOG.md`,
   the README, and `docs/releases/<version>.md`. Document supported behavior and compatibility
   requirements. Keep the maintained benchmark evidence compact and portable; retain raw execution outputs separately.
2. Run the required fmt/clippy/test checks above with the locked dependencies. Run any explicitly
   enabled external-data regressions needed for the changed behavior.
3. Build and check the publication package:

   ```bash
   cargo build --locked --release --package kallistors --bin kallistors
   target/release/kallistors version
   cargo publish --package kallistors --locked --dry-run
   ```

   Use `--allow-dirty` when reviewing uncommitted package changes. The dry run compiles
   the crate package without publishing. Check that the package includes the license
   text and read fixtures, and excludes downloaded data and generated benchmark outputs.
   Keep `crates/kallistors/LICENSE` identical to the root `LICENSE`.
4. Commit the reviewed release contents and tag that commit `v<version>`. Pushing the tag starts
   publication; the release notes are consumed directly by the binary workflow.

Two GitHub workflows handle release publication:

- `Publish to crates.io` runs on version tags (`v*` and `[0-9]*`). It verifies
  that the tag matches the crate version, runs fmt/clippy/tests, performs a
  locked dry run, then publishes when `CARGO_REGISTRY_TOKEN` is configured.
- `Release binaries` runs on version-tag pushes, when a GitHub Release is published, and can also
  be dispatched manually for an existing tag. It checks out the tag, verifies the
  tag/version match, builds `kallistors` with `cargo build --locked --release`
  for Linux x86_64, macOS x86_64, macOS arm64, and Windows x86_64, then uploads
  archives plus `.sha256` checksum files to the release. On a tag push or manual dispatch, it
  creates a missing GitHub Release using `docs/releases/<version>.md` when available, otherwise
  generated GitHub notes. Archives include README, license, changelog, documentation, and benchmark scripts/manifests.

## Kallistors CLI debugging tools

- `trace-reads` supports EC traces, per-hit dumps, intersection dumps, and
  positions visited.
- Minimizer debug includes MPH lookup info, unitig match decisions, and
  special/D-list flags.
- Debug helpers include `minimizer-lookup`, `minimizer-bitmap-scan`,
  `minimizer-unitig-mphf-check`, and `minimizer-mphf-keys` for inspecting
  MPHF/index consistency.
- Use `target/release/kallistors <command> --help` for the current CLI flag
  list; debug flags change more often than the user-facing quant interface.

## Abundance diffs

Compare `abundance.tsv` outputs and rank transcripts by relative error:

```bash
scripts/compare_abundance.py \
  --kallisto /path/to/kallisto/abundance.tsv \
  --kallistors /path/to/kallistors/abundance.tsv \
  --top 50
```

For an audit that includes effective lengths, common six-digit precision, missing targets,
expressed-transcript errors, and total TPM mass difference:

```bash
python3 scripts/compare_abundance.py \
  --kallisto /path/to/kallisto/abundance.tsv \
  --kallistors /path/to/kallistors/abundance.tsv \
  --summary-json target/validation/abundance-audit.json
```

The official kallisto TSV uses six significant digits. Compare with high-precision diagnostic
output when distinguishing floating-point differences from substantive abundance mismatches.
Exact `run_info.json` counters alone do not establish full EC or abundance parity.

Compare complete EC count tables from `kallistors pseudoalign` and the diagnostic upstream build:

```bash
python3 scripts/compare_ec_counts.py \
  --reference /path/to/upstream/diagnostic.ec_counts.tsv \
  --candidate /path/to/kallistors.ec.tsv \
  --report target/validation/ec-comparison.json
```

Both runs must use the same index: transcript IDs are index-specific. The script ignores EC
numbering, compares counted transcript sets, and exits nonzero on a mismatch. The
[benchmark guide](benchmarks.md) gives the complete reproduction commands and phase boundaries.

## Full-reference and performance checks

Prepare GENCODE v50 ALL and the public libraries using the [benchmark guide](benchmarks.md).
The full-index regressions require both original- and Rust-built indexes of that same reference:

```bash
cargo build --locked --release
target/release/kallistors index -t 8 \
  -i data/srr-small/gencode.v50.kallistors.idx \
  data/srr-small/gencode.v50.transcripts.fa.gz
KALLISTORS_TEST_BIN="$PWD/target/release/kallistors" \
  cargo test --locked -p kallistors --test full_index -- --ignored
```

`KALLISTORS_TEST_RUST_INDEX` selects another Rust-built reference. Use complete EC comparisons
across both builders and quantifiers; summary counters alone cannot establish compatibility.
For a mismatch, create a seeded paired subset with `subsample_benchmark_pairs.py`, then use
`trace-reads` and the [upstream diagnostic build](kallisto_debug_build.md) to inspect affected reads.

Use `benchmark_reads.py --stages` for complete paired/single runs, phase timing and peak RSS.
It checks every counted EC before timing. `benchmark_index.py --threads 8 --keep-last-index`
measures construction and retains the final indexes for cross-loading checks. Both tools use
one warmup and three measured runs by default. Keep input, executable and output hashes;
compare identical reads and references, retain CPU-activity flags and use a new result directory
after changing an executable. Native HDF5 checks and abundance comparisons follow timing.
See the benchmark guide for runnable commands, manifests, trimming and reporting.

## Source layout

| Responsibility | Modules |
| --- | --- |
| CLI orchestration | `commands/quant.rs`, `commands/pseudoalign.rs`, `commands/bus.rs` |
| Index loading and building | `index/`, `pseudoalign/bifrost.rs` |
| Read matching | `pseudoalign/core/`, `pseudoalign/single.rs`, `pseudoalign/paired.rs` |
| Threaded input and EC aggregation | `pseudoalign/threaded.rs` |
| EC dumps, minimizer lookup and read traces | `pseudoalign/diagnostics/{ec_dump,minimizer,tracing}.rs` |
| Fragment models, EM and bootstraps | `quant/mod.rs`, `quant/components.rs` |
| Abundance TSV/HDF5 and run metadata | `quant/output/{tsv,hdf5,metadata}.rs` |
| BUS input, records and technology definitions | `commands/bus/{input,output,technology}.rs` |
| BUS genome projection and BAM output | `commands/bus/{genome,pseudobam}.rs` |
| BUS unit/integration tests | `commands/bus/tests.rs`, `tests/bus.rs` |

Public diagnostic and quantification output APIs are re-exported from their parent modules.
Keep module boundaries separate from algorithm changes. Changes to mapping must preserve
candidate/fallback order, EC contribution order, cache ownership and wide-value handling.
Changes to EM must preserve per-transcript addition order and convergence rules.

Validate structural changes with complete EC maps, matched-model abundances, native HDF5
readback and diagnostic/BAM output checks. Compare one- and eight-thread whole-process timing
and memory with alternating controls. Use the same executable for all validation of a candidate.

## Quantification

- **Input:** direct `flate2::bufread::MultiGzDecoder` with the Rust `zlib-rs` backend and
  1 MiB compressed/decoded buffers. CRC, concatenated-member and truncation checks remain active.
  Single-thread paired counting parses sequence-only batches synchronously. Diagnostic paths
  retain full records, and sampling follows the same record order and sequence-buffer quotas.
- **Loading:** packed graph words, shared EC lists, compact MPHF ranks and bounded bitmap
  scratch reduce allocations. Deferred transcript names use the same open index file and preserve
  every output row. Shade/bias/positional loaders keep eager names. Abundant-node keys remain
  explicitly associated with their EC payloads; source candidate order is preserved.
- **Matching:** jump and midpoint probes advance across ambiguous bases to the next valid k-mer,
  preserving upstream iterator positions and end-of-read behavior. Caches verify complete keys and mode/index ownership. Conservative 12-bit
  fingerprints check all nine possible g-mer anchors for k=31/g=23; conflicting or special
  buckets bypass rejection. Full keys still prove a hit. Encoded fallback queries, exact overlap
  reuse and bounded AArch64 vector paths preserve ties, strand priority and scalar fallbacks.
- **Paired counting:** once sampling finishes, unused placements are omitted. Unchanged ECs
  remain borrowed through mate counting; filtering or divergent intersections materialize owned
  vectors. Per-worker counts share collision-checked EC keys. Bounded queues and rendezvous
  dispatch retain all reads and release blocked producers on input failure.
- **EM:** graph storage is released first; bootstrap input is retained when needed. Independent
  components run in parallel without changing each transcript's contribution order. Only bitwise
  fixed components skip repeated work. Mandatory iterations, convergence and the final threshold
  update remain intact; bias correction uses the general path. Work is partitioned by connected
  component, so a large component can limit EM parallelism.
- **Output:** bounded ordered TSV chunks preserve serial bytes. Threaded output overlaps TSV
  and HDF5. DEFLATE level 6 and chunks of at most 16,384 elements preserve dataset names,
  types and shapes; native HDF5 readback independently verifies interoperability.

## Index construction

The builder produces deterministic v13-compatible nucleotide indexes with exact key checks and
full-width/scalar fallbacks. Important boundaries:

1. Buffered FASTA parsing preserves Unicode/error handling, input order, original lengths,
   ambiguity RNG order and polyA normalization. Independent normalization runs in parallel.
2. Repeat-window exclusions require full byte equality. Shared-end/internal anchors only find
   candidates; seam windows and all transcript metadata remain. Bounded parallel collections
   losslessly pack sequences and deduplicate words before radix sorting.
3. Canonical-edge symmetry, exact bucket descriptors and batched probes avoid duplicate lookup
   work. Idempotent endpoint updates preserve both ends. Parallel path preparation feeds ordered
   serial claims; accepted branch sequences are materialized in parallel.
4. Compact ranks, oriented starts, joins and locations retain complete-domain fallbacks. Counted
   disjoint fills and contiguous transcript/EC position storage bound temporary allocations.
   Identifier bytes share an arena until owned names can be restored exactly.
5. Minimizer hashing shares exact REP arithmetic with mapping and retains every tie. Parallel
   MPHF collision merging preserves survivor order. Small bitmap serialization is byte-checked
   against CRoaring, with the general serializer as fallback; node output remains ordered.
6. Completed scratch is released before later allocations overlap. Page retirement on macOS/Linux
   touches only wholly unused page-aligned ranges and preserves allocator ownership/live bytes.
   The requested worker pool handles construction and independent cleanup.

Stage clocks are exclusive: FASTA, graph, EC, minimizer/MPHF and write. The internal total includes
final allocation cleanup; whole-process timing also includes startup and pool teardown. Complete
index bytes and cross-loaded EC maps protect compatibility; speed or summary counts alone do not.
