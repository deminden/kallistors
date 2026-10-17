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
| `full_index` | Opt-in GENCODE v50 scanner, abundant-node, paired fallback and reverse-unitig SVIL regressions |

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
   requirements. Keep generated benchmark outputs outside the maintained documentation.
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

## GENCODE v50 comparison and regression workflow

The [benchmark guide](benchmarks.md) documents public ENA retrieval, checksums, fasterp
trimming, shared references and timing commands. [Fragment sampling](fragment_length_sampling.md)
describes the paired default. Regression tests exercise the current library and CLI, including
regressions for fragment models, abundant-node associations and fallback probes.

```bash
# Prepare reads and the reference independently.
cargo install --locked fasterp --version 0.2.1
python3 scripts/prepare_benchmark.py --samples SRR13638690
python3 scripts/prepare_reference.py --reference gencode50

# Prepare the shared reference with an upstream v0.52.0 executable.
target/tools/kallisto-0.52.0 index -t 8 \
  -i data/srr-small/gencode.v50.kallisto.idx \
  data/srr-small/gencode.v50.transcripts.fa.gz

# Enable full-reference scanner and abundant-node regressions explicitly.
cargo build --locked --release
target/release/kallistors index -t 8 \
  -i data/srr-small/gencode.v50.kallistors.idx \
  data/srr-small/gencode.v50.transcripts.fa.gz
KALLISTORS_TEST_BIN="$PWD/target/release/kallistors" \
  cargo test --locked -p kallistors --test full_index -- --ignored

# Optional fresh eight-thread comparison; final outputs are retained for validation.
python3 scripts/benchmark_reads.py \
  --samples SRR13638690 --thread-counts 8 \
  --diagnostic-kallisto target/tools/kallisto-src/build/src/kallisto \
  --hdf-python target/hdf5-native-check/bin/python \
  --out target/benchmarks/gencode50-release
```

The upstream executable is downloaded separately, as described in the guide. The normal suite
includes synthetic empirical-fragment tests; downloaded-index tests are ignored by default.
The benchmark runner saves completed run records as it progresses. Abundance outputs from the final
measured repetition are retained for native HDF5 and abundance validation.

## Per-read EC traces

Trace per-read EC decisions for a small set of single-end reads:

```bash
printf "READ_ID_1\nREAD_ID_2\n" > read_list.txt
target/release/kallistors trace-reads \
  --index data/srr-small/gencode.v50.kallisto.idx \
  --reads data/srr13638690/trimmed-fasterp/SRR13638690_1.trimmed.fastq.gz \
  --read-list read_list.txt \
  --out read_traces.tsv \
  --fragment-length 200
```

## Real-data EC parity

Use the complete EC comparison in [the benchmark guide](benchmarks.md#correctness-and-interpretation).
The quantification runner exports original ECs with the tracked diagnostic
patch, runs current `pseudoalign`, and fails on any counted transcript-set difference.
For faster investigation, construct paired subsets with `scripts/subsample_benchmark_pairs.py`
and pass the same subset/index to both tools. Use `trace-reads` on selected reads after a
mismatch is found.

## Index builder validation

Build and inspect an index from a supplied small transcript FASTA:

```bash
cargo build --release
mkdir -p target/validation
target/release/kallistors index \
  -i target/validation/reduced.idx \
  -t 4 \
  --timings \
  /path/to/small-transcripts.fa

target/release/kallistors index-info --index target/validation/reduced.idx
target/tools/kallisto-0.52.0 inspect target/validation/reduced.idx
```

Compare index builds with one warmup and three alternating measured runs:

```bash
python3 scripts/benchmark_index.py \
  --fasta data/srr-small/gencode.v50.transcripts.fa.gz \
  --threads 8 --out target/validation/index-build
```

For compatibility validation, run paired quant on the same reads with:

- an upstream `kallisto index` output
- a `kallistors index` output loaded by `kallistors quant`
- a `kallistors index` output loaded by upstream `kallisto quant`

Compare complete EC maps as well as `run_info.json` counters. Use the abundance-diff workflow
above with matched fragment models to check effective lengths and estimated counts.

## Performance comparison

The [benchmark guide](benchmarks.md) describes current-versus-original controls, repeated
stage timings, CPU/power context, RSS/footprint and complete EC/HDF5 checks. Use
`scripts/benchmark_reads.py --stages` for optional phase timing alongside paired/single
total time and peak RSS on macOS/Linux. It verifies complete EC parity before any timing.
Use a new output directory after changing either executable.

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

The threaded loader retains packed sequence metadata and compact MPHF rank records, with
full-width fallbacks. On AArch64, bounded cache hints fetch upcoming inline descriptors and
overflow slots during ordered loading; all writes retain their original order. Exact reciprocal division and peeled initial hash levels change lookup
work without changing index bytes or candidate order. Per-worker caches verify complete keys,
matching flags and index ownership; scratch allocations are reused between reads. Intersections
retain sorted-set semantics, including singleton and unbalanced-size paths.

The ordinary quantification loader reads two-bit graph sequences directly, without duplicate
ASCII unitigs or k-mers. Cold sequence helpers decode packed bases on demand; full and positional
loaders retain their existing fields. Flat EC blocks share adjacent identical transcript lists and
use narrow metadata only when every bound fits. Canonical fallback entries retain complete keys,
original candidate order and explicit abundant-node associations.

Eligible k=31, g=23 graphs build the aligned runtime minimizer table in one lookup pass. Temporary
first-position/ordinal records retain singletons; repeated and special buckets retain ordered
events. Finalization restores source candidate order and the same full-width runtime descriptors.
Larger layouts use the general loader. Gap-encoded bitmap scratch has exact-sized allocation and
raw fallback, and redundant MPHF rank arrays are released only after their packed copy is complete.

Paired workers keep independent counters while sharing identical complete EC keys through sharded
intern tables. Compact 32-bit fingerprints and head IDs reduce per-worker lookup storage, with
full-width fallback IDs. Every hash collision is resolved by comparing transcript vectors. Pool roots are
released before merging, and public EC lists retain their owned-vector interface. Paired work uses
rendezvous dispatch and bounded mate queues; parser chunks still assemble into upstream's unchanged
8 MiB sampling buffers. Queue sizes do not alter the sampling quota or omit reads.

Ordinary quantification scans every transcript name for shade annotations, using a reusable
buffer, then defers name allocation until output. The deferred reader retains the same open
index file; full IDs and UTF-8 replacement behavior are preserved. Shade references, bias,
positional and diagnostic loaders keep eager names. Every transcript still appears in output.

Mapping releases graph storage before EM. Owned EC input is released once the EM representation
has been prepared; bootstrap runs retain the input they need. EM transcript IDs narrow in place
when their full domain fits. Wholly unused scratch pages can be retired on macOS/Linux while
normal allocator ownership and live boundary bytes are preserved.

For ordinary k=31, g=23 indexes, an optional byte fingerprint rejects a bucket only when none
of the read k-mer's nine canonical g-mer anchors can match it. Checking only the selected
minimizer would lose relaxed matches. Conflicting fingerprints and special/KM buckets always
bypass rejection. Fingerprints travel through the loader's sorted count batches and do not
change the serialized v13 index.

After the paired fragment quota is exhausted, unstranded counting without bias correction folds
hits directly into ECs, retaining online jump/rejection behavior and the separate final all-hit
intersection. Placement records are still collected whenever sampling or another consumer needs
them. A hard rejection of mate 1 can reject the pair before scanning mate 2; both records have
already been parsed and the processed-pair count is retained.

Without bias correction, EM groups independent transcript/EC components while retaining each
transcript's original contribution order. Each component computes its denominators before
updating transcripts. Mandatory initial iterations omit convergence bookkeeping; subsequent
iterations use the original stopping predicate. Only exact equality of every alpha bit permits
skipping a component. Entirely fixed branches cache that fact and skip scheduling with their
original component boundaries intact. Thresholding resets component and branch state before
the final update. Bias correction
retains the general EM path. Regression tests compare serial and threaded floating-point bits
against the independent reference implementation.

TSV formatting uses bounded, ordered chunks and retains the serial writer's exact bytes.
Default output still includes HDF5 with DEFLATE level 6 and bounded 16,384-element chunks.

## Index construction

Buffered FASTA parsing checks ASCII once per input buffer and directly appends proven-ASCII
sequence lines. Pending fragments retain the conjunction of their buffer proofs; Unicode headers
and invalid UTF-8 retain full validation, whitespace and CR/LF behavior. Parallel
normalization uppercases and validates independent sequences; ambiguity substitutions still
consume MT19937 in transcript/base order. Original lengths precede identical polyA trimming.
Single-threaded construction retains the serial path. The [shared DNA encoder](../crates/kallistors/src/dna.rs)
serves mapping, graph initialization and minimizer keys, with byte/position/orientation oracles.

Collection excludes windows only after complete byte equality proves they already occur in an
earlier transcript in the same bounded group. Shared prefix/suffix anchors and sparse internal
word anchors select candidates. Seam windows and all transcript metadata remain. Collection uses
bounded 2,048-transcript groups with sparse internal anchors selected at 1/16 density. Portable word
comparisons locate the exact differing byte after equal 32-byte blocks. Group deduplication precedes
parallel radix sorting. Completed shards append immediately; final ordering is independent of worker
completion order. Large parallel
collections losslessly pack each group before sorting its words, after borrowed anchors have
finished. Disjoint mutable groups are dispatched to the requested workers; packed records are restored to original transcript
order. The distinct table keeps its virtual allocation until lookup mapping finishes, avoiding a
second full-table copy. Fully unused initialized tail pages can be discarded without moving live keys.

[Packed transcript scratch](../crates/kallistors/src/index/packed_transcripts.rs) uses two bits per
normalized base and extracts exact k-mer words directly, including unaligned byte boundaries and
tails. Shared group buffers retain every sequence until its mapping completes. After packing consumes
the ASCII sequences, only original transcript lengths remain; the emptied transcript-record
allocation is retired before topology construction. Identifier bytes
share one UTF-8 arena during topology construction; owned names are reconstructed unchanged after
mapping. The name table hashes each identifier once, uses one entry probe and verifies full
equality for every hash match, including complete collisions. Bounded DNA anchor maps use fast word hashing
with complete-key comparisons.

Lookup descriptors store complete singleton suffixes/ranks where they fit. Two-, three- and
four-key guards reject definite misses and verify full keys on positives. Larger graphs use exact
rank ranges: shared block bases and packed 12-bit deltas, with exact u32 dense-block and
large-rank fallbacks. Active query prefixes avoid dummy probes; AArch64 nonblocking cache hints overlap
independent key reads and positive atomic updates without changing equality checks.
Canonical batches cover up to 32 sources with 16-bit source/orientation metadata. Active query and
metadata buffers are reused within bounded worker chunks. Other architectures retain ordinary reads.
Canonical extension selection queries one reverse-complement representative and updates both
endpoints through idempotent CAS, including duplicate self-edge cases. Compact records retain
rank/orientation/degree information and exact continuation eligibility, with full-width fallback.
Finalization first gives nonunique ranks stable degree sentinels. Atomic neighbor reads then
update only continuation flags; rank fields remain unchanged, so no separate degree snapshot is
needed. Bounded prefetch batches overlap those reads. All workers join before immutable ownership
transfer; matching size/alignment and native byte order are required, with ordinary collection as
fallback. Possible join successors are collected in source/orientation order while the edge
storage is live. The edge allocation is then retired before oriented-start sorting, link
allocation and joined-sequence assembly. The owned query words preserve the same candidate order
and complete-word duplicate checks; an empty query set leaves the paths unchanged.

Parallel branch/residual preparation never claims vertices. Serial ordered claims truncate at the
first earlier canonical claim, preserving orientation, unitig order and reverse collisions. Prepared
groups retain their exact claimed prefixes without a separate count array. The continuation layout
is selected once per group. Path prefixes use bounded SIMD decoding on little-endian AArch64 with
a scalar fallback. Accepted branch sequences materialize in parallel; residual work retains bounded batches and a long-path
serial fallback. Prepared steps use u32 when ranks fit, otherwise u64. Compact oriented starts sort
by exact numeric keys. A counting pass divides exact spare capacity into disjoint worker slices;
length is published only after every worker succeeds. Records implement Copy; on an error the
vector retains length zero.
All-active, all-inactive and sparse masks retain the same IDs and keys. Exact prefix ranges and
equality scans of small ranges retain complete word and duplicate checks. The paths retain exact first/last canonical ranks, avoiding repeated endpoint membership
searches during query collection. Byte-only internal constructors retain the checked lookup
fallback. Ordered query resolution preserves join priority without another sort. Freed path
preparation scratch is collected on its existing worker owners before joining. Accepted path metadata uses bounded ordered chunks and exact
final assembly. Joined sequences reuse owned buffers. Missing join links use an out-of-domain sentinel, preserving
state zero and cycles. Transcript runs advance through proven nonbranching unitig interiors; the
independent byte oracle checks every skipped continuation in test builds. Parallel location lookup
batches span short unitigs and preserve complete canonical membership, offsets and orientations;
the serial path retains bounded per-unitig batches.

[Transcript-run scratch](../crates/kallistors/src/index/transcript_runs.rs) stores exact whole-unitig
runs in three words and partial runs in five. Large unitig IDs retain the full representation.
Known bounds are reconstructed from unitig lengths before stable per-unitig assembly.
[Block positions](../crates/kallistors/src/index/block_positions.rs) retain sorted unique positions
and group boundaries in contiguous arrays, avoiding an allocation for every transcript in an EC.
Lookup buffers are released before EC/minimizer payloads. ECs, minimizer work, MPHF marking and grouped
node serialization run in parallel. Worker-private seen/collision sets merge both within-worker
and cross-worker collisions exactly. Seen/collision pairs share one local record, and construction
reuses the mapper's exact remainder from a floor reciprocal with boundary/modulo oracles. Survivor and
serialized bitvector order are unchanged.
Survivor IDs and temporary buckets use u32 only when their complete domains fit; usize/u64
fallbacks serialize identically. Numeric minimizer keys are converted to their original big-endian
bytes at hashing/serialization boundaries, avoiding a full key-layout copy.
Index construction and mapping share the exact grouped REP-state initializer. Rolling minimizer
hashes feed a fixed monotonic queue retaining every tied minimum. Numeric
radix sorting preserves key byte order. Tiny bitmap and packed-sequence oracles protect serialization;
Native small-array and single run bitmap paths use the locked CRoaring representation only when
its choice is proven; general cases use the library serializer. Tests compare every emitted byte
through cardinality, container, adjacency and run-size boundaries. AArch64 sequence packing has a
scalar fallback. Independent allocations are released in the same
requested worker pool. [Scratch-page retirement](../crates/kallistors/src/index/scratch_pages.rs)
releases only page-aligned, wholly unused initialized ranges on macOS/Linux. Live bytes, allocator
headers and ownership are preserved; spare capacity remains reusable. Dead starts, candidates,
join links, consumed path metadata and the completed graph-section buffer are retired at their last
use. Owner-pool allocator collection skips allocations supplied by system-allocator library hosts.
Stage clocks do not overlap; MPHF is separate from writing and the internal
total includes final cleanup. Whole-process time additionally includes startup and pool teardown.

The shared-window union oracle protects prefix/suffix and internal seams, changed flanks and
mutations. Serial/parallel graph checks cover canonical claims, joins, narrow/wide steps, reverse
collisions and long cycles; lookup and MPHF oracles protect exact keys and serialized bytes.
Real-reference validation checks complete index bytes, then complete paired/single EC maps in both
quantifiers across both builders. Summary counters alone are insufficient. See the
[benchmark guide](benchmarks.md) for reference preparation, trimming and phase/memory measurement.
