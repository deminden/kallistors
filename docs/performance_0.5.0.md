# Quantification performance — kallistors 0.5.0

Version 0.5.0 accelerates index loading, matching and EM while preserving complete EC count parity
with original kallisto. This report describes the measured upstream comparison and correctness
checks for the 0.5.0 implementation.

## Comparison with original kallisto

Apple Silicon macOS, **eight threads**, **4,411,003 fastp-trimmed SRR13638690 read pairs**, and the
same upstream-built **GENCODE v50 ALL-transcript index with 670,670 targets**. Full process timing
includes index loading, every read, EM, default HDF5/TSV output and shutdown.

| Build | Elapsed time | Peak RSS |
| --- | ---: | ---: |
| Original kallisto v0.52.0 | 47.673446 s | 2.85 GiB |
| kallistors 0.5.0 | 17.532188 s | 3.01 GiB |

**2.72x faster total quantification; 63.22% less elapsed time; 5.7% more peak resident memory.**

This is one measured pair with already warmed inputs and no additional warmup, not a three-run
median. Both processes were measured using macOS `/usr/bin/time -l`; optional Rust stage timings
were omitted. No compiler or correctness job ran concurrently. The before-run CPU snapshot
recorded WindowServer at 45.0% before original kallisto, and no other process above the harness's
25% threshold before kallistors. Background activity and the single-pair method limit the precision
of the speed estimate. Results are specific to this workload and host; no new single-end or
index-build speed claim follows from this measurement.

Both tools report 4,411,003 processed fragments, 4,246,701 pseudoaligned and 209,744 unique.
Complete EC parity is checked independently below.

## Memory allocation diagnosis and reduction

An allocation-capacity audit of the GENCODE v50 ALL index found **1,831,997,688 bytes** in
precomputed minimizer words and their sequence offsets. The loader materialized an eight-byte
canonical word at every g-window, including positions absent from the stored bitmaps, and retained
these tables during mapping and EM. The existing packed-base stores occupied **116,427,136 bytes**
and already contained the information needed to reconstruct each word. These are allocation sizes,
not RSS measurements.

The loader now extracts only bitmap-selected words from packed bases using constant-time shifts,
bit reversal and two-bit complement. It retains the same canonical representation, MPHF lookup
and ordered bucket append. No full per-position minimizer table is allocated.

Ordinary quantification also retained nested EC blocks alongside the flat EC table. Jump distances
read bounds from the nested copy, preventing its removal. Jumps now use identical flat bounds,
and ordinary `quant` and `pseudoalign` use the compact counting loader. Positional fragment filters,
strand filtering, detailed diagnostics and alternate scanners retain the nested metadata they need.

The minimizer count buffer becomes prefix offsets in place, avoiding another **400,817,368-byte**
allocation on this reference. Serialized MPHF bytes, node scratch, bitmaps and fill cursors are
released after their last use. Packed matching, parallel MPHF resolution and parallel EM remain
in place.

### Peak memory measurement

Peak RSS is the maximum resident set size over the complete process lifetime, including index
loading, mapping, EM and output writing. The macOS report expresses it in bytes: **3,057,057,792**
for original kallisto and **3,230,564,352** for kallistors. Kallistors uses **5.7% more peak memory**
(**0.16 GiB more**, **1.057x**) in the measurement above. GiB means 2^30 bytes. This is process
peak resident memory, not average usage, allocated bytes or total machine memory. It covers all
threads sharing the process address space. Positional and diagnostic modes retain additional
metadata; their memory usage is not represented by this paired counting measurement.

## Performance changes

- Calculate on-list cardinality only when a synthetic D-list hit needs it. Ordinary reads avoid
  scanning all 670,670 transcript flags.
- Pack unitig bases into u64 words and extract candidate k-mers with shifts and one/two word loads.
  Reconstruct loader minimizers from these same words and avoid duplicate EC graphs during counting.
- Reuse validated read encodings for minimizer representation and adjacent read words. Roll the
  exact REP hash states between adjacent g-windows, retaining bounds, ties and candidate order.
- Preserve buffered index bytes when skipping node k-mers. Resolve MPHF lookups in parallel, then
  append events in original order so bucket ordering is preserved.
- Restrict repeated EM work to active transcripts after preserving the uniform initial prior.
  Compute EC denominators independently and update transcripts in parallel, preserving each
  transcript's original EC/edge contribution order and floating-point expression.
- Refresh transposed EM weights when bias changes effective lengths. Preserve zero-iteration
  behavior, convergence thresholds, iteration limits and bootstrap sampling.
- Reuse positive primary matches and retain MPHF scratch per worker. Cached primary misses still
  populate fallback candidates; unique index/k/g ownership prevents stale cache reuse. Oversized
  minimizer tie lists are recomputed in full.
- Skip discarded FASTQ lines without copying their bytes. Hash numeric EC keys efficiently while
  comparing complete vectors for equality on collisions.

`quant -t N` uses N threads for mapping and EM in successive phases. Existing library quantification
entry points remain serial; `em_quantify_with_fragment_histogram_threaded` accepts an explicit
thread count. Index format, CLI syntax and output layout are unchanged.

## Exact EC validation against upstream

All checks use the same GENCODE v50 ALL index and fastp-trimmed inputs. Complete counted transcript
sets are compared independently of per-run EC numbering, using diagnostic original kallisto output.

| Dataset | Paired nonzero ECs | Single-end nonzero ECs | Different classes / count L1 |
| --- | ---: | ---: | --- |
| ERR2260014 | 3,604 | 3,534 | 0 / 0 |
| DRR023395 | 181,044 | 146,676 | 0 / 0 |
| SRR13638690 | 190,711 | 192,645 | 0 / 0 |

Exact EC parity does not imply byte-identical upstream abundance. Kallistors retains all qualifying
paired fragment observations; upstream samples near a bounded goal. The remaining paired abundance
differences, bootstrap RNG differences and bias-enabled real-data validation limits are described
in the [quantification audit](quantification_fixes_2026-10-14.md).

Numerical regressions compare serial and threaded EM exactly with sparse targets, zero iterations,
Gaussian/empirical fragment models and bias-weight refresh. Encoding tests cover word boundaries,
every supported k/g size, tied minimizers, mixed case, invalid bases and scanner jumps. Compact/full loader regressions compare bucket order, EC values, bounds and complete read hits
in both scanners; positional loading retains strand annotations. Prefix-offset tests verify buffer
reuse, including empty buckets. Further regressions cover cache ownership, complete tie lists, loader bucket order, small FASTQ buffers,
CRLF, truncated records and concatenated gzip members. The full-index TTN/AGAP12P regression checks
both scanners. Formatting, strict Clippy and workspace tests passed.

## Reproduce the upstream comparison

Prepare the inputs using the [read comparison guide](read_comparison_2026-10-14.md), then
build kallistors and compare the two tools from the repository root:

```bash
cargo build --locked --release
python3 scripts/benchmark_quant_variants.py \
  --control target/tools/kallisto-0.52.0 \
  --candidate target/release/kallistors \
  --threads 8 --warmups 0 --runs 1 --no-timings --memory \
  --out target/benchmarks/kallisto-comparison
```

`--memory` records peak RSS using `/usr/bin/time` on macOS or Linux. The harness saves the raw
memory report per run and normalizes macOS bytes / Linux KiB to bytes in its JSON output. It reports
the median of measured per-process peaks when more than one run is requested. Omit `--memory` for
a timing-only comparison. Use `--warmups 1 --runs 3` for a repeated comparison. `--index`, `--reads` and `--reads2` select
another paired workload. `--candidate-first` reverses starting order. The output directory must be
new or empty. The harness saves commands, before-run CPU snapshots, wall times and final outputs;
cleanup happens outside timing. Equal run-info counters alone do not establish exact EC parity:
use `scripts/compare_ec_counts.py` with diagnostic upstream EC output as the reference.
