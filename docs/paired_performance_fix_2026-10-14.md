# Eight-thread paired quantification bottleneck — 2026-10-14

This report covers an unreleased optimization after v0.4.2, based on commit
`a8795c2165048578c6a72479f83f692cb0dedc79`. The published v0.4.2 binaries do not contain it.
The input is the same fastp-trimmed SRR13638690 sample (4,411,003 pairs), shared upstream-built
GENCODE v50 ALL-transcript index (670,670 targets), and Apple Silicon host used in the
[earlier comparison](local_read_comparison_2026-10-14.md). Every full-read run here uses eight threads.

## Diagnosis

The v0.4.2 scanners begin each read by computing
`index.onlist.as_deref().map(onlist_cardinality)`. The cardinality helper scans every flag in the
on-list to obtain the dummy transcript ID used for synthetic D-list hits. Ordinary reads never
need that dummy ID, but still pay for the complete reference scan. Across both mates, this input
causes approximately 5.9 trillion flag visits.

Two macOS `sample` captures found the scanner as the dominant active leaf symbol. Disassembly of
the hottest scanner offsets shows vector byte loads and widened addition reductions, matching
the on-list counting loop. The FASTQ producers were mostly waiting on full worker queues;
decompression was not the main bottleneck. These stack samples were taken during mapping, not EM.

The instrumented upstream v0.52.0 build was also run with all optional per-read diagnostics
disabled. Timestamped stderr locates the boundaries between loading, mapping, and EM. Upstream
stores its D-list dummy unitig in the loaded index and uses it when a D-list hit requires it;
ordinary matching does not perform this Rust reference-wide flag count.

| Stage | kallistors v0.4.2 profile | Fixed kallistors first run | Instrumented upstream |
| --- | ---: | ---: | ---: |
| Index graph decode + minimizer count/fill | 7.416 s | 7.609 s | About 5.585 s through startup |
| Mapping / EC merge | 62.286 s | 16.758 s | About 11.269 s |
| EM | 9.017 s | 9.015 s | About 31.211 s |
| Full process wall time | 79.687 s | 34.275 s | 48.733 s |

The Rust stage timers for threaded decompression and pseudoalignment overlap: producer time
includes queue waits, and the pseudoalignment timer includes worker joining and EC merge. Do not
sum those two timers. Upstream's boundaries are inferred from stderr timestamps and are not
identical instrumentation. The baseline wall time includes two three-second sampling captures;
it is diagnostic evidence, not an unprofiled timing median.

## Fix and preservation of behavior

Both [fast](../crates/kallistors/src/pseudoalign/core/fast.rs) and
[general](../crates/kallistors/src/pseudoalign/core/read.rs) scanners retain an optional reference
to the on-list, use its presence for the same branch decisions, and calculate cardinality only
when actually appending a synthetic D-list hit. The on-list remains immutable during a read, so
this delays the same pure calculation without changing its value.

There is no index-format or library-structure change, no stored count that could become stale if
library callers modify the public on-list, and no change to search candidates, jumps, EC
intersections, fragment sampling, or EM. The new synthetic
[regression](../crates/kallistors/tests/onlist_dummy.rs) checks normal and D-list reads with both
scanners, plus a matching D-list with no on-list. Existing TTN/AGAP12P coverage checks the prior
fallback-probe correction.

## Complete EC validation

The complete EC count maps are compared independently of per-run EC numbering. The reference
maps come directly from the diagnostic upstream outputs preserved during the
[v0.4.2 audit](quantification_fixes_2026-10-14.md). The new
[comparison script](../scripts/compare_ec_counts.py) reports different classes and total absolute
count error and returns a failure exit status on any mismatch.

All six complete maps remain exact against upstream:

| Dataset | Mode | Nonzero ECs | Different classes | Absolute count difference |
| --- | --- | ---: | ---: | ---: |
| ERR2260014 | paired | 3,604 | 0 | 0 |
| ERR2260014 | single | 3,534 | 0 | 0 |
| DRR023395 | paired | 181,044 | 0 | 0 |
| DRR023395 | single | 146,676 | 0 | 0 |
| SRR13638690 | paired | 190,711 | 0 | 0 |
| SRR13638690 | single | 192,645 | 0 | 0 |

The full large paired run has 4,411,003 processed fragments, 4,246,701 aligned, and 209,744 unique,
matching upstream. The comparison JSON is `ec/upstream-summary.json` under the artifact directory.

Paired abundance also remains numerically unchanged from the preserved v0.4.2 binary in the
first and final measured fixed outputs: all 670,670 effective lengths are identical, maximum
estimated-count difference is 1.63e-9, and maximum TPM difference is 1.09e-11. Tiny differences reflect floating-point summation order
with threaded EC collection; the fragment model and observation histogram are unchanged.
The residual abundance differences from upstream documented in the earlier audit therefore
remain. This optimization fixes speed, not the separate fragment-sampling policy difference.

Validation passed:

```bash
cargo fmt --all
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
cargo test --locked --workspace --all-features
KALLISTORS_TEST_BIN="$PWD/target/release/kallistors" \
  cargo test --locked -p kallistors --test single_ttn_regression -- --ignored
```

## Reproduce EC comparison

With the downloaded inputs and diagnostic upstream binary prepared:

```bash
mkdir -p target/validation
KALLISTO_DIAGNOSTIC_ECS=1 target/tools/kallisto-src/build/src/kallisto quant \
  -i data/srr-small/gencode.v50.kallisto.idx -t 8 \
  -o target/validation/paired-upstream \
  data/srr13638690/trimmed/SRR13638690_1.trimmed.fastq.gz \
  data/srr13638690/trimmed/SRR13638690_2.trimmed.fastq.gz

target/release/kallistors pseudoalign --threads 8 \
  --index data/srr-small/gencode.v50.kallisto.idx \
  --reads data/srr13638690/trimmed/SRR13638690_1.trimmed.fastq.gz \
  --reads2 data/srr13638690/trimmed/SRR13638690_2.trimmed.fastq.gz \
  --out target/validation/paired-fixed.ec.tsv

python3 scripts/compare_ec_counts.py \
  --reference target/validation/paired-upstream/diagnostic.ec_counts.tsv \
  --candidate target/validation/paired-fixed.ec.tsv \
  --report target/validation/paired-ec-comparison.json
```

## Focused timing comparison

The official upstream release binary is used for end-to-end timing. This quick check uses one
warmup and one measured run per tool, at eight threads and paired-end only, with outputs retained
outside timing. It is a single-run check, not a repeated median or a scaling benchmark:

```bash
python3 scripts/benchmark_local_reads.py \
  --samples SRR13638690_larger --modes paired --thread-counts 8 --runs 1 \
  --keep-last-output --out-prefix target/paired-performance-2026-10-14/quant-comparison
```

The runner retains its default of one warmup and three measured runs, while `--modes` and `--runs`
allow a focused check. No profiler, compiler, tests, or EC analysis ran concurrently with these
end-to-end timings:

| Tool | Warmup | Measured paired wall time |
| --- | ---: | ---: |
| Official kallisto v0.52.0 | 47.350541 s | 48.939757 s |
| Fixed kallistors development build | 33.761472 s | 33.955053 s |

The fixed run takes 30.62% less elapsed time than upstream, or a 1.4413x speedup in this check.
Both tools report exactly 4,411,003 processed, 4,246,701 aligned, and 209,744 unique fragments.
Raw records and retained quantification outputs share the `quant-comparison` prefix under the
artifact directory. The first stage-timed fixed run was also close at 34.275 seconds, but this
focused comparison is not a repeated-median estimate and does not establish performance for
other datasets or thread counts.

Local profiling artifacts, stage logs, EC maps, comparison JSON, and the preserved v0.4.2 binary
are under ignored `target/paired-performance-2026-10-14/`. The upstream diagnostic build and its
tracked patch are documented in [the build guide](kallisto_debug_build.md).
