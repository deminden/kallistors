# Development notes

This page collects developer-facing workflows for checks, parity debugging, and
abundance investigation. For the patched upstream `kallisto` build and its extra
tracing flags, see [kallisto_debug_build.md](kallisto_debug_build.md).

## Toolchain

Development requires Rust 1.97 or newer. The workspace uses Rust edition 2024 and declares the
minimum supported Rust version in the root `Cargo.toml`.

## Required checks

Before committing changes, run:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

## Release workflow

The crate version lives in `crates/kallistors/Cargo.toml`. Version 0.5.0 is prepared in this
checkout. Its [release notes](releases/0.5.0.md) describe the performance changes,
validation evidence, and remaining compatibility differences.

For a release:

1. Bump the crate version and refresh its package entry in `Cargo.lock`. Update `CHANGELOG.md`,
   the README, and `docs/releases/<version>.md`; retain dates and versions on historical reports.
2. Run the required fmt/clippy/test checks above with the locked dependencies. Run any explicitly
   enabled external-data regressions needed for the changed behavior.
3. Build and check the publication package:

   ```bash
   cargo build --locked --release --package kallistors --bin kallistors
   target/release/kallistors version
   cargo publish --package kallistors --locked --dry-run
   ```

   During preparation before committing, add `--allow-dirty` to the dry run. It creates and
   compiles the crate package without publishing. Check that the package includes the license
   text and read fixtures, and excludes downloaded data and generated benchmark outputs.
   Keep `crates/kallistors/LICENSE` identical to the root `LICENSE`.
4. Commit the reviewed release contents and tag that commit `v0.5.0`. Pushing the tag starts
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
  generated GitHub notes. Archives include README, license, changelog, and documentation.

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
[paired performance report](paired_performance_fix_2026-10-14.md) includes the complete reproduction
commands and the profiling procedure's interpretation.

## GENCODE v50 comparison and regression workflow

The [read comparison guide](read_comparison_2026-10-14.md) documents the three public ENA
runs, download checksums, fastp trimming, index building, and multicore timing commands. The
[fix report](quantification_fixes_2026-10-14.md) records the latest complete EC validation and
the remaining fragment-sampling differences.

```bash
# Optional R/fastp environment; R is not needed to run kallistors.
mamba env create -f environment.yml
mamba run -n kallistors-r bash scripts/prepare_benchmark.sh

# Prepare the shared reference with an upstream v0.52.0 executable.
target/tools/kallisto-0.52.0 index -t 8 \
  -i data/srr-small/gencode.v50.kallisto.idx \
  data/srr-small/gencode.v50.transcripts.fa.gz

# Enable the full-reference TTN/AGAP12P scanner regression explicitly.
cargo build --locked --release
KALLISTORS_TEST_BIN="$PWD/target/release/kallistors" \
  cargo test --locked -p kallistors --test single_ttn_regression -- --ignored

# Optional fresh eight-thread comparison; final outputs are copied outside timing.
python3 scripts/benchmark_reads.py \
  --samples SRR13638690_larger --thread-counts 8 --keep-last-output \
  --out-prefix target/benchmarks/gencode50-release
```

The upstream executable is downloaded separately, as described in the guide. The normal suite
includes synthetic empirical-fragment tests; the downloaded-index test remains ignored by default.
The benchmark runner saves completed run records as it progresses. Abundance outputs are retained
only for the final measured repetition when `--keep-last-output` is enabled.

## Per-read EC traces

Trace per-read EC decisions for a small set of single-end reads:

```bash
printf "READ_ID_1\nREAD_ID_2\n" > read_list.txt
target/release/kallistors trace-reads \
  --index data/gencode.v49_kallisto.idx \
  --reads data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  --read-list read_list.txt \
  --out read_traces.tsv \
  --fragment-length 200
```

## Real-data subset parity

Run deterministic paired-prefix parity against the real-read dataset:

```bash
python3 scripts/real_subset_parity.py \
  --mode paired \
  --threads 32 \
  --sizes 64,256,1024,4096,16384,65536,131072,262144 \
  --kallisto-bin ./kallisto_src/build/src/kallisto \
  --kallistors-bin ./target/release/kallistors \
  --report /tmp/kallistors-parity.tsv
```

When investigating an experimental fast-path branch, pass `--fast-env` to run
the candidate under the relevant environment flag and emit `trace-compare`
output for the first divergent read:

```bash
python3 scripts/real_subset_parity.py \
  --mode paired \
  --threads 32 \
  --sizes 64,256,1024,4096 \
  --fast-env KALLISTORS_FAST_DELTA_BOUNDED_INCREMENTAL_SCAN=1 \
  --trace-out /tmp/kallistors-trace-compare.tsv
```

## Index builder validation

Build and inspect a reduced real-transcript index:

```bash
cargo build --release
target/release/kallistors index \
  -i /tmp/reduced.kallistors.idx \
  -t 4 \
  --timings \
  data/subsets/reduced_64_transcripts.fa.gz

target/release/kallistors index-info --index /tmp/reduced.kallistors.idx
./kallisto_src/build/src/kallisto inspect /tmp/reduced.kallistors.idx
```

Run a large GENCODE build benchmark:

```bash
/usr/bin/time -o /tmp/kallistors-index.time -f 'elapsed=%E maxrss=%MKB' \
  target/release/kallistors index \
  -i /tmp/gencode.kallistors.idx \
  -t 8 \
  --timings \
  data/gencode.v49.transcripts.fa.gz
```

Compare against upstream and write a Markdown report:

```bash
python3 scripts/bench_index_build.py \
  --fasta data/gencode.v49.transcripts.fa.gz \
  --threads 8 \
  --kallisto-bin ./kallisto_src/build/src/kallisto \
  --kallistors-bin ./target/release/kallistors \
  --out target/validation/index_build_benchmark.md
```

For compatibility validation, run paired quant on the same reads with:
- an upstream `kallisto index` output
- a `kallistors index` output loaded by `kallistors quant`
- a `kallistors index` output loaded by upstream `kallisto quant`

The expected first-order check is exact `run_info.json` count parity. `abundance.tsv` comparisons
use the same tolerance-based diff workflow described above.

## Performance comparison

The [0.5.0 performance report](performance_0.5.0.md) records eight-thread comparisons
against original kallisto, CPU contention, complete EC validation and the
`scripts/benchmark_quant_variants.py` workflow.
