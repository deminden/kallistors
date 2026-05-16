# Development notes

This page collects developer-facing workflows for checks, parity debugging, and
abundance investigation. For the patched upstream `kallisto` build and its extra
tracing flags, see [kallisto_debug_build.md](kallisto_debug_build.md).

## Required checks

Before committing changes, run:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

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

Run deterministic paired-prefix parity against the checked-in real dataset:

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
