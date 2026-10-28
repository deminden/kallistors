# Development notes

This page collects developer-facing workflows for checks, parity debugging, and
abundance investigation. For the patched upstream `kallisto` build and its extra
tracing flags, see [kallisto_debug_build.md](kallisto_debug_build.md).

## Toolchain

Development requires Rust 1.97 or newer. The workspace uses Rust edition 2024 and declares the
minimum supported Rust version in the root `Cargo.toml`. Source builds use
`target-cpu=native` from `.cargo/config.toml`; CI overrides `RUSTFLAGS` with
`-C target-cpu=generic` for portable release artifacts and validation.

## Using the Rust library

Add the dependency shown in the [README](../README.md#as-a-crate), then use the examples
below to build an index or quantify paired-end reads.

Serial quantification APIs use one thread. For parallel EM, call
`em_quantify_with_fragment_histogram_threaded` with an explicit thread count. The CLI
`quant -t N` uses N threads for mapping and EM in successive phases. HDF5 output uses
bounded chunks.

Exhaustive `PseudoalignOptions` literals must include `fragment_length_samples` and
`match_cache_size`. Use `..PseudoalignOptions::default()` for the 10,000-observation goal
and automatic [match-cache budget](match-cache.md). Set `fragment_length_samples` to zero
to collect all qualifying fragment observations. This option controls fragment estimation;
every read pair still contributes to EC counting. See the
[fragment-sampling policy](fragment_length_sampling.md) for quotas and threaded behavior.

### Build an index

```rust
use std::error::Error;
use std::path::{Path, PathBuf};

use kallistors::index::{build_index, IndexBuildOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let fasta = vec![PathBuf::from("transcripts.fa.gz")];
    build_index(
        Path::new("transcripts.idx"),
        &fasta,
        IndexBuildOptions::default(),
    )?;
    Ok(())
}
```

### Quantify paired-end reads

```rust
use std::error::Error;
use std::path::Path;

use kallistors::index::Index;
use kallistors::io::open_fastq_reader;
use kallistors::pseudoalign::{
    build_bifrost_index, pseudoalign_paired_bifrost_with_options, PseudoalignOptions, Strand,
};
use kallistors::quant::{em_quantify_with_fragment_histogram, EcCountsInput, QuantOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let index_path = Path::new("path/to/index.idx");
    let mut reads_1 = open_fastq_reader(Path::new("reads_1.fq.gz"))?;
    let mut reads_2 = open_fastq_reader(Path::new("reads_2.fq.gz"))?;
    let index = build_bifrost_index(index_path)?;
    let ec = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut reads_1,
        &mut reads_2,
        Strand::Unstranded,
        PseudoalignOptions::default(),
    )?;

    let meta = Index::load(index_path)?;
    let input = EcCountsInput {
        ec_list: kallistors::ec::EcList {
            classes: ec.ec_list,
        },
        counts: ec.counts,
    };
    let lengths: Vec<u32> = meta.transcripts.iter().map(|t| t.length).collect();
    let quant = em_quantify_with_fragment_histogram(
        &input,
        &lengths,
        None,
        None,
        ec.fragment_length_hist.as_deref(),
        QuantOptions::default(),
    )?;

    println!("reads processed: {}", ec.reads_processed);
    println!("reads aligned: {}", ec.reads_aligned);
    println!("targets: {}", quant.est_counts.len());
    Ok(())
}
```

## Required checks

Before committing changes, run:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

Repository checks use Python and Bash without downloading data or running benchmarks:

```bash
python3 -B scripts/check_repository.py
python3 -B -m unittest discover -s scripts/tests -v
```

These check Markdown links/anchors, shell examples, Python syntax and corrupt benchmark-record
fixtures. The benchmark handoff gate requires each declared condition's warmup and all three
distinct measurements, unchanged executables and uninterrupted native/probe timings.

CI runs formatting and Clippy on Linux, tests stable Rust on Linux/macOS/Windows, and tests the
declared minimum Rust 1.97 on Linux. The lightweight repository job is separate from native
benchmarks. Release builds also run `python3 -B scripts/smoke_cli.py PATH_TO_BINARY`: a tiny
index, mapping and quantification workflow with expected counts and default HDF5 output.

The Windows MSVC CLI reserves an 8 MiB main-thread stack because the unoptimized generated
argument parser exceeds the linker's 1 MiB default. The package build script applies this
to the executable, including installed builds; worker stacks and library callers are unchanged.
The initial committed stack size remains the linker default.

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

Arithmetic oracles compare non-NaN results bit for bit, including subnormals, infinities and
signed zero. Check NaN propagation by classification: Rust does not guarantee arithmetic NaN
payloads or signs, even for identical expressions. Use runtime operands when testing optimized
kernels. See [Rust's NaN contract](https://doc.rust-lang.org/std/primitive.f32.html#nan-bit-patterns).

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

## Mapping route selection

Select routes by required metadata. Ordinary stranded mapping needs an orientation
hit; fragment filtering and active bias sampling also need a best placement.
Union, no-jump and Shade strand filters retain their complete hit collection.
Diagnostic traces remain separate from counting. Bias sampling can finish while
the accumulator remains available to quantification; its presence alone must
not force full mapping metadata for subsequent reads.

BUS uses packed batches for ordinary single-cDNA layouts, including compound
barcode/UMI slices. Synthetic barcodes, header/tag extraction and paired layouts
use the existing extraction rules to prepare bounded batches for parallel mapping.
Ordinary BAM follows the same mapping/output layer after tag decoding. Results
are emitted in input order, and each pair contributes its eligible fragment length
to the complete BUS histogram. Do not substitute bulk's sampled histogram here.
Answer caches are scoped to the strand policy, including tag-disabled reads.

Packed positional loading retains transcript coordinates and orientation bits;
it changes in-memory graph storage without changing v13 serialization. Alternate
scanners and requested diagnostics retain their established representations.

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

The crate version lives in `crates/kallistors/Cargo.toml`. Its compact package README lives in
`crates/kallistors/README.md` and uses repository URLs so documentation links also work outside
a checkout. Version-specific release notes
live in `docs/releases/`; the binary publication workflow reads the file matching its tag.

For a release:

1. Bump the crate version and refresh its package entry in `Cargo.lock`. Update `CHANGELOG.md`,
   both READMEs, and `docs/releases/<version>.md`. Document supported behavior and compatibility
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
  for Linux x86_64, macOS x86_64, macOS arm64, and Windows x86_64, exercises each binary on a
  tiny index-to-quantification workflow, then uploads
  archives plus `.sha256` checksum files to the release. On a tag push or manual dispatch, it
  creates a missing GitHub Release using `docs/releases/<version>.md` when available, otherwise
  generated GitHub notes. Archives include README, license, changelog, documentation, and benchmark scripts/manifests.

## Kallistors CLI debugging tools

Developer commands are hidden from the ordinary command list and remain callable by name.
Run `kallistors COMMAND --help` for their options. Existing diagnostic names and accepted
compatibility switches remain available; the legacy no-op and unsupported quant BAM flags
are hidden from normal help as well.

- `index-ec-dump` and `ec-from-index` export index EC data; `ec-stats`, `ec-export` and
  `ec-compare` inspect and canonicalize EC text. `transcript-lookup` resolves target names.
- `unitig-window` inspects graph sequence windows and minimizers.
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
one warmup and three measured runs by default; use that protocol and report medians for the
bulk and original-kallisto construction panels. Keep input, executable and output hashes;
compare identical reads and references, retain CPU-activity flags and use a new result directory
after changing an executable. Native HDF5 checks and abundance comparisons follow timing.
See the benchmark guide for runnable commands, manifests, trimming and reporting.

## Source layout

| Responsibility | Modules |
| --- | --- |
| CLI orchestration | `commands/{quant,pseudoalign,bus}.rs` |
| Graph construction | `index/graph_build.rs` coordinates `graph_build/{collection,compaction,joins,locations,payloads,minimizers}.rs` |
| Exact graph lookup and edge storage | `index/graph_build/{lookup,edges,keys,relative_edges,sparse_locations}.rs` |
| Index loading | `pseudoalign/bifrost.rs` exposes constructors; `bifrost/{loader,encoded,ec_storage,buckets,decoding}.rs` own loading, storage and decoding |
| Native index codecs | `index/{binary,bitmap,hashing}.rs` share integer encoding, bitmap decoding and exact minimizer hashing; native and portable-compatible bitmap policies remain explicit |
| MPHF storage and construction | `index/bifrost.rs` owns lookup storage; `index/writer/mphf.rs` owns native MPHF construction and serialization |
| Read matching | `pseudoalign/core/fast/{matching,rescue,scanner}.rs`; `core/read.rs` dispatches to the fast path or `read/{probe,general,finish}.rs` for ordered general probes, scanning and final EC folding |
| Worker caches and ownership | `pseudoalign/core/runtime.rs` owns thread-local context; `runtime/{match_cache,minimizer_cache,block_cache,budget}.rs` own cache records, lookup policies and scoped allocation |
| Single/paired counting and sampling | `pseudoalign/{single,paired}.rs`; single-end fragment/strand filtering is shared across serial and batched counting |
| Threaded input and EC aggregation | `pseudoalign/threaded.rs` |
| EC dumps, minimizer lookup and read traces | `pseudoalign/diagnostics/{ec_dump,minimizer,tracing,paired}.rs` |
| EM scheduling and arithmetic | `quant/mod.rs` coordinates `quant/work.rs`; `quant/components.rs` groups independent components and `components/{scalar,parallel}.rs` schedule their exact updates |
| Fragment/bias models, input and bootstraps | `quant/{models,input,bootstrap}.rs`; `quant/bootstrap/sampling.rs` owns the shared exact sampling plan |
| Abundance TSV/HDF5 and run metadata | `quant/output/{tsv,hdf5,metadata}.rs` |
| BUS configuration, tags and general alignment | `commands/bus/{configuration,tags,alignment}.rs` |
| BUS FASTQ/BAM input, records and technology definitions | `commands/bus/{input,bam_input,output,technology}.rs` |
| BUS packed input, ordered parallel counting and bounded reader shutdown | `commands/bus/fast.rs`; `bus/{normalized,bam_counting,counting}.rs` share mapping for compound FASTQ and BAM input |
| Exact shared cDNA reuse and collision/replacement regression | `commands/bus/fast/query_cache.rs` |
| BUS genome projection and BAM output | `commands/bus/{genome,pseudobam}.rs` |
| BUS unit/integration tests | `commands/bus/tests.rs`; `tests/bus.rs` contains shared fixtures and `tests/bus/` groups distinct regressions |
| Graph regressions | `index/graph_build/tests/` groups collection, compaction, lookup/edges, minimizers and payload checks |
| Codec, cache and EM regressions | Existing tests live beside their implementations in `index/{bifrost,writer}/`, `index/skipped_position_tests.rs`, `core/runtime/` and `quant/components/tests.rs` |

Public index-loading, diagnostic and quantification APIs are re-exported from their parent modules.
Implementation modules stay private to their subsystem. Keep independent naive/fallback mapping
paths as correctness oracles; their lack of production fast-path traffic does not make them unused.
Keep module boundaries separate from algorithm changes. Changes to mapping must preserve
candidate/fallback order, EC contribution order, cache ownership and wide-value handling.
Changes to EM must preserve per-transcript addition order and convergence rules.

### Mapping routes

| Workload | Route |
| --- | --- |
| Ordinary single-end and paired-end bulk, unstranded or FR/RF | Encoded fast scanner at one or multiple threads; counting batches reuse worker scratch |
| Fragment and bias sampling | Retain placement/hit observations while sampling; use count-only shortcuts when the required quotas are complete and filters permit |
| Standard single-cDNA BUS FASTQ layouts | Packed sequence input and exact shared cDNA-answer cache |
| Compound, tagged and paired BUS FASTQ layouts; ordinary BAM counting | Normalize required input fields, then use shared counting; single-cDNA reads reuse the exact answer cache and paired reads use the paired counting API |
| Shade indexes, union and alternate compatibility policies | General scanner preserves the requested filtering and candidate order |
| Traces, BAM projection, long-read, amino-acid and unmapped diagnostics | Retain the metadata and processing required by the requested output; scanner eligibility is checked independently |

Bulk quantification, non-debug pseudoalignment and ordinary BUS counting use
`prepare_read_mapping`, including compound FASTQ and BAM counting. It checks the same
eligibility policy as the fast scanner before preparing bounded selectors. Selectors currently
support k=31/g=23, retain their existing 48 MiB payload ceiling and only narrow the candidate
search: every selected candidate still requires exact verification. Shade indexes and explicit
general-scanner policies skip preparation. Legacy `prepare_single_cell_candidate_cache*`
APIs retain their unconditional preparation behaviour and share the same selector builder.

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
  every output row. Ordinary unstranded single-end loading keeps the graph packed while
  retaining all positions needed by the fragment filter. Grouped positions use inline storage
  or one lossless coordinate arena per block. Exact value and offset bounds select 16-, 24- or 32-bit
  descriptors; repeated coordinates remain full-width, with wide offset/length fallback.
  Packing dispatches on descriptor width once per block.
  Native singleton position bitmaps decode directly to their complete coordinate and sign,
  avoiding a temporary one-element allocation. Consecutive native singleton groups pack
  directly into the same lossless descriptors without an intermediate position-set vector;
  On 64-bit targets, one complete 32-bit coordinate, two 24-bit descriptors or three 16-bit
  descriptors fit in the existing group descriptor without allocation. Each narrower descriptor
  retains its own strand bit and is selected only after all coordinate bounds are proved.
  Larger coordinates, repeated positions and 32-bit targets retain allocated storage.
  Compact arenas proved to contain only singleton positions omit indirect-entry flags.
  Native singleton-group validation also computes the maximum coordinate in the same pass.
  Other encodings retain the general decoder.
  Positional bitmap payloads are decoded directly from the buffered node bytes, without
  intermediate serialized buffers.
  Large packed references batch graph lookups and use bounded, ordered parallel annotation
  decoding when workers are available. Small references and single-worker loading retain
  the serial decoder. Paged annotation ownership and narrow node directories avoid millions
  of separate allocations; wide spans, exceptional encodings and ordered errors retain
  their complete fallback behavior. One- and two-transcript ECs fit inline where lossless.
  Parallel unitig counting advances eight independent MPHF keys together and gathers compact
  rank words before testing their hit bits. Lane order, scalar hash/rank results and full-width
  or final-map fallbacks remain exact; scratch storage is bounded to one batch.
  Positional and mapping views share the same owned transcript IDs and block bounds.
  Strand summaries are derived from raw positions instead of cached per-block arrays. Bias,
  strand-aware and diagnostic loaders retain their required metadata and sequence access.
  Shade/bias loaders keep eager names. Abundant-node keys remain
  explicitly associated with their EC payloads; source candidate order is preserved.
- **Matching:** jump and midpoint probes advance across ambiguous bases to the next valid k-mer,
  preserving upstream iterator positions and end-of-read behavior. Caches verify complete keys and mode/index ownership. Conservative 12-bit
  fingerprints check all nine possible g-mer anchors for k=31/g=23; conflicting or special
  buckets bypass rejection. Full keys still prove a hit. Encoded fallback queries, exact overlap
  reuse and bounded AArch64 vector paths preserve ties, strand priority and scalar fallbacks.
  Read minimizer initialization shares exact forward and reverse-complement group lookups
  for 23-base windows; other lengths retain the general initializer.
  Ordinary primary/tail probes and Bifrost rescue retain their distinct candidate families.
  Encoded rescue preserves their ordering and delegates special/overcrowded markers to the
  full matcher when required.
  Invocation-scoped worker caches reuse complete intersections of immutable annotation pairs.
  Full identities and values prove reuse; bounded answer storage and index lifetimes prevent
  stale references. Sorted multiset fallbacks preserve duplicate multiplicity and ordering.
- **Paired counting:** once sampling finishes, unused placements are omitted. Unchanged ECs
  remain borrowed through mate counting; filtering or divergent intersections materialize owned
  vectors. Per-worker counts share collision-checked EC keys. Bounded queues and rendezvous
  dispatch retain all reads and release blocked producers on input failure.
- **Single-end counting:** threaded workers share complete, collision-checked EC keys through
  the same pool used for paired reads; immutable keys store exactly their transcript payload.
  Ordinary reads fold ECs during scanning and retain first/best placement without
  collecting a full hit list. Both synchronous and threaded counting keep unchanged ECs borrowed
  through exact fragment filtering; filtering materializes a borrowed EC only after its first
  rejected transcript, while owned intersections filter in place. Only new counted keys are
  copied. Borrowed annotation ranges and early coordinate acceptance avoid temporary position
  vectors while preserving boundary padding, transcript order and the supplied fragment model.
  When fragment and strand filters address the same annotation block, they share its decoded
  positions while retaining independent acceptance rules, including mixed-strand groups.
  Threaded dispatch uses 4,096-read batches and one prefetched batch with rendezvous delivery;
  receivers are released before joining blocked producers on failure. Bias, comprehensive
  strand filters and diagnostics retain their required hit metadata.
- **BUS strand loading:** explicit FR/RF filtering loads transcript orientation annotations
  for single-cDNA FASTQ, paired layouts and BAM input. The unstranded compact loader stays
  separate. Diagnostic commands use the shared strand-aware mapping loader even when no
  fragment filter is requested. Reverse-complement direct-kmer entries are keyed by the
  query sequence, with orientation recorded in the entry flag.
- **BUS counting:** simple single-cDNA layouts use the compact count-only scanner. A bounded
  mapping/output queue preserves input and first-seen EC order. The shared sequence cache
  uses two-way replacement and verifies both a hash and full bytes under a short shard lock;
  mapping runs outside locks. Its answers belong
  to one immutable index/options invocation, borrowed ECs retain index lifetimes, and owned
  intersection storage is capped. Equal-length unequal ECs skip the impossible subset proof;
  count-only scanning stops once the final EC fold has irreversibly rejected a read. Full
  diagnostics and placement paths retain their existing behavior. Positive/negative collisions
  and concurrent replacement are
  covered by a focused native-index regression. Decoders honor the remaining input limit before
  forming each batch; interleaved limits count complete record groups.
- **EM:** graph storage is released first; bootstrap input is retained when needed. Independent
  components run in parallel without changing each transcript's contribution order. Scheduling
  adapts to the counted EC graph: a component larger than one worker's fair share and 262,144
  edges shares its independent EC denominators and transcript updates across the same pool.
  All component denominators finish before alpha updates; scalar additions stay ordered, and only
  integer change counts and fixed-point flags are reduced. Smaller components keep local
  iterations for cache reuse and lower synchronization cost. Eligible ordinary and bias-corrected
  runs use the same compact component engine and ordered update kernels. Bias refreshes project
  active alpha values to the original target order, refresh weights and invalidate fixed-point
  certificates at their original iteration boundaries. Iteration-zero, absent-model and other
  exceptional cases retain the general path.
  A changing component can witness that convergence is still impossible; other components then
  catch up in bounded tiles using the same updates. This preserves the first global stop,
  mandatory iterations and final threshold update. Only bitwise fixed components skip repeated
  updates. Incoming sums can also omit a proven positive-zero result when the previous alpha
  and singleton are both positive zero and every actual weight and current scale is finite.
  All denominators still execute before this decision; signed-zero and nonfinite cases retain
  literal arithmetic. Work is partitioned by connected component and within dominant components,
  so a large connected graph can use the requested threads without approximating its ECs.
- **Bias model:** observed/expected histograms and their normalization remain ordered. Independent
  per-transcript effective-length rescales can use the existing requested EM pool. One shared
  equation implementation has serial and parallel scheduling specializations; small adapters
  preserve scalar execution without maintaining a second model implementation.
- **Bias FASTA:** release completed mapping storage before allocating transcript sequences.
  Header lookups borrow existing eager names and assign sequence buffers directly to requested
  slots. Preserve requested order, independent copies for duplicate requests, last-record wins
  for duplicate FASTA headers and UTF-8 validation of the complete source. Ordinary quantification
  keeps deferred-name loading; the eager-name view does not allocate a second name table.
- **Bootstraps:** independent samples use bounded requested workers and a conservative
  2 GiB budget for parallel worker scratch. Large single samples reuse the ordered incoming
  component engine; when several samples leave spare workers, their components share the same
  pool. Tiny or oversized models retain scalar scheduling. Precomputed root seeds and indexed
  collection preserve sample order. Compact active coordinates retain the full-target uniform
  prior and final scatter; eligible scalar samples narrow IDs in their existing allocation.
  Shared kernels preserve contribution order, first convergence stop and final pruning.
  Mandatory rounds and intermediate component updates omit only unobservable change counts.
  Exact parallel sampling addresses independent SplitMix chunks and commits them only when
  every original rejection check succeeds; otherwise it replays serially from the untouched
  RNG. Serial callers can share an exact inverse-CDF table of at most 64 MiB; larger libraries
  and multiworker sampling use the bounded search directory. Main TSV writing finishes before
  bootstrap work starts, keeping worker pools within the request.
- **Output:** bounded ordered TSV chunks preserve serial bytes. Threaded output overlaps TSV
  and HDF5. DEFLATE level 1 and chunks of at most 16,384 elements preserve dataset names,
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
   Location scratch selects an exact representation from graph bounds: the existing direct
   32-bit layout covers at most 2^22 unitigs and 2^30−1 positions. When that layout no longer fits,
   graphs with at most 2^31−1 positions use a global position plus strand bit in 32 bits. Ordered
   starts and one unitig bucket per 64 positions recover the full coordinate. Remaining graphs
   use the full-width layout. These are transient build layouts; graph membership, ordered claims,
   decoded coordinates and serialized v13 data are preserved.
   Above 2^28 keys, bounded collection avoids long-transcript worker buffers growing without
   limit. When prefix and rank bounds permit, four-byte edge pairs encode exact 14-bit relative
   ranks. Ten-base roots split only overfull ranges through a sparse prefix trie; the oriented
   word and outgoing base recover strand and complete rank. Wider layouts remain available.
   Directions, outgoing codes and degree/non-self sentinels remain exact. Smaller graphs keep
   their existing formats. `--timings` reports logical edge-scratch bytes separately from RSS.

   Large graphs with suffixes fitting 40 bits pack keys in place from eight to five bytes
   **before topology**. Disjoint parallel waves preserve unread words, original allocation
   alignment and all live bytes. Complete prefix, suffix and padding checks establish membership.
   A temporary prefix-owned filter is budgeted at 12 bits/key, rounded to power-of-two blocks;
   three hash-derived bit positions share one word, reducing false admissions with one random
   read. It rejects definite misses, and every positive still undergoes exact key verification.
   Compact boxed blocks retain one owner per prefix and use a zero sentinel for empty ranges.
   Completed worker scratch is collected before filter allocation; the filter is retired before
   path construction.

   Reference walks use sparse full-word landmarks when their conservative temporary memory
   bound beats the dense layout. Unitig endpoints and every transcript-start k-mer cover forward,
   reverse and interior-start walks. Prefix-owned hash tables retain full key equality and full
   coordinates; dense direct/global layouts remain for other graph shapes. This allows full
   topology keys to be retired before transcript-run mapping.

   Batched exact probes reuse two complete prefix bounds and use compact rank results.
   Residual path preparation uses larger batches for large graphs while serial claims retain original seed order. Join queries retain complete
   words and 40-bit source identifiers in 13 bytes. Word-reconstruction hints are retired before
   join probes, whose membership checks no longer need them.

5. Minimizer hashing shares exact REP arithmetic with mapping and retains every tie. Parallel
   MPHF collision merging preserves survivor order, generates initial IDs implicitly on large
   builds, compacts later survivors in place, reserves exact survivor capacity and retires
   consumed shards. Serialized graph bitmaps are released before MPHF construction. Small bitmap
   serialization is byte-checked against CRoaring, with the general serializer as fallback; node output remains ordered.
6. Completed scratch is released before later allocations overlap. Page retirement on macOS/Linux
   touches only wholly unused page-aligned ranges and preserves allocator ownership/live bytes.
   The requested worker pool handles construction and independent cleanup.

Stage clocks are exclusive: FASTA, graph, EC, minimizer/MPHF and write. The internal total includes
final allocation cleanup; whole-process timing also includes startup and pool teardown. Complete
index bytes and cross-loaded EC maps protect compatibility; speed or summary counts alone do not.
