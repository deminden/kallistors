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

The crate version lives in `crates/kallistors/Cargo.toml`. For a release, bump
that version, refresh `Cargo.lock`, update docs that mention the release state,
run the required checks above, then tag the checked commit with the matching
version (`v0.4.1` for crate version `0.4.1`).

Two GitHub workflows handle release publication:
- `Publish to crates.io` runs on version tags (`v*` and `[0-9]*`). It verifies
  that the tag matches the crate version, runs fmt/clippy/tests, performs a
  locked dry run, then publishes when `CARGO_REGISTRY_TOKEN` is configured.
- `Release binaries` runs when a GitHub Release is published, and can also be
  dispatched manually for an existing tag. It checks out the tag, verifies the
  tag/version match, builds `kallistors` with `cargo build --locked --release`
  for Linux x86_64, macOS x86_64, macOS arm64, and Windows x86_64, then uploads
  archives plus `.sha256` checksum files to the release.

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

Run deterministic paired-prefix parity against the local real dataset:

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
