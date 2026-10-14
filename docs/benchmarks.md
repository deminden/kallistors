# Benchmarks

This page indexes benchmark results and correctness evidence. The inputs and raw artifacts are not
part of the repository. Results are machine- and cache-sensitive; compare only runs measured with
the same inputs, binaries, and timing method.

The paired/single-end comparison using GENCODE v50 and three ENA human RNA-seq runs is recorded
in [read_comparison_2026-10-14.md](read_comparison_2026-10-14.md). Its download,
trimming, repeated index-build, and multicore quantification workflows are implemented by
`scripts/prepare_benchmark.sh`, `scripts/benchmark_index.py`, and
`scripts/benchmark_reads.py`.

The [quantification mismatch audit](quantification_mismatch_audit_2026-10-14.md) identifies
paired-end effective-length differences and a minimal single-end TTN scanner reproducer.
The [fix report](quantification_fixes_2026-10-14.md) records both corrections, full EC validation,
and the repeated eight-thread comparison after the fixes. The original October 14 timing tables
are historical measurements from before these corrections.

## 0.5.0 performance evidence

The [performance report](performance_0.5.0.md) records the implementation included in
0.5.0: full fastp-trimmed SRR13638690 paired reads, GENCODE v50 ALL transcripts, Apple Silicon,
eight threads, full process wall time including load, default HDF5 output and shutdown.

| Tool | Full elapsed time | Peak RSS |
| --- | ---: | ---: |
| Original kallisto v0.52.0 | 47.673 s | 2.85 GiB |
| kallistors 0.5.0 | 17.532 s | 3.01 GiB |

**2.72x faster; 63.22% less elapsed time; 5.7% more peak memory (0.16 GiB extra).**
One measured pair with already warmed inputs and no additional warmup, using macOS
`/usr/bin/time -l`. RSS covers the complete process lifetime; GiB = 2^30 bytes. CPU contention
is recorded in the [performance report](performance_0.5.0.md), alongside exact byte counts,
allocation diagnosis, limitations and reproduction. Complete EC maps match original kallisto in
all six cases. Earlier diagnosis is in the
[on-list bottleneck report](paired_performance_fix_2026-10-14.md).

## Historical release evidence — 0.4.2

GENCODE v50 ALL transcripts, shared upstream-built index, fastp-trimmed SRR13638690, Apple Silicon
macOS, eight threads, one warmup and the median of three measured runs:

| Mode | kallisto v0.52.0 | kallistors with 0.4.2 fixes | EC count map |
| --- | ---: | ---: | --- |
| paired | 46.694 s | 78.123 s | Exact |

In this historical 0.4.2 measurement, the large paired workload was slower in kallistors.
Complete EC count maps match for all three October datasets in both read modes. Remaining paired abundance differences arise from fragment
observation sampling; the largest paired TPM mass difference fell from 3.94% to 0.080%.
The large single-end post-fix timing was stopped after warmups, so there is no new single-end
median. See the [fix report](quantification_fixes_2026-10-14.md) for provenance and limitations.

## Historical May 19, 2026 benchmark

The sections below describe a different host, reference, release snapshot, and timing protocol.

### Environment

- Date: 2026-05-19
- OS: Linux x86_64
- CPU: AMD Ryzen 9 7950X3D
- Thread counts: `1, 2, 4, 8, 16, 32`
- Timing method: `/usr/bin/time -f '%e %M'`
- Builds: release binaries
- Upstream kallisto: `kallisto_src/build/src/kallisto`, version `0.52.0`
- kallistors: `target/release/kallistors`, benchmark snapshot from the pre-`0.4.0`
  development line

## Paired-End Quant

### Full Dataset, 32 Threads

Input:
- Reads:
  - `data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz`
  - `data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz`
- Read pairs: `4,408,640`
- Index: `data/gencode.v49_kallisto.idx`
- Flags: `-t 32`

Method:
- one warmup run for each binary
- five measured runs for each binary
- alternating kallisto and kallistors runs
- table reports median elapsed time and median peak RSS
- count columns report the latest post-fix full-file validation on the same input

| tool | median elapsed | mean elapsed | min | max | median RSS | n_processed | n_pseudoaligned | n_unique |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| kallisto | `55.81s` | `56.49s` | `55.44s` | `57.76s` | `2784 MiB` | `4,408,640` | `4,244,771` | `276,251` |
| kallistors | `32.79s` | `32.38s` | `31.53s` | `33.09s` | `5504 MiB` | `4,408,640` | `4,244,771` | `276,251` |

Result:
- kallistors median wall time is `1.70x` faster than upstream kallisto at `-t 32`.
- kallistors uses about `2.7 GiB` more peak RSS on this run.
- Post-fix full-file paired validation matches upstream kallisto counts exactly:
  `n_processed=4,408,640`, `n_pseudoaligned=4,244,771`, `n_unique=276,251`. The deterministic
  `1,048,576`-pair prefix in the multicore table below is exact at all measured thread counts.

Raw artifacts:
- `target/benchmarks/paired_after_fix_2026-05-19/full_paired_after_fix_runs.tsv`
- `target/benchmarks/paired_after_fix_2026-05-19/full_paired_after_fix_summary.json`
- `KALLISTORS_RUN_REAL_PAIRED_REGRESSION=1 cargo test -p kallistors --test real_paired_regression`
  runs the isolated full-index pair regression; normal CI skips the expensive index load.

### Multicore Scaling, 1,048,576 Pairs

Input:
- Reads:
  - `data/subsets/real_mate1_n1048576.fastq.gz`
  - `data/subsets/real_mate2_n1048576.fastq.gz`
- Read pairs: `1,048,576`
- Index: `data/gencode.v49_kallisto.idx`

Method:
- one warmup run per binary per thread count
- five measured runs per binary per thread count
- alternating kallisto and kallistors runs
- table reports median elapsed time and median peak RSS

| threads | kallisto median | kallistors median | speedup | kallisto RSS | kallistors RSS | run_info parity |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | `56.87s` | `76.81s` | `0.74x` | `2228 MiB` | `4180 MiB` | exact |
| 2 | `41.63s` | `45.65s` | `0.91x` | `2360 MiB` | `4580 MiB` | exact |
| 4 | `31.41s` | `30.21s` | `1.04x` | `2380 MiB` | `4658 MiB` | exact |
| 8 | `27.08s` | `22.42s` | `1.21x` | `2435 MiB` | `4749 MiB` | exact |
| 16 | `25.76s` | `19.19s` | `1.34x` | `2542 MiB` | `4947 MiB` | exact |
| 32 | `25.47s` | `17.56s` | `1.45x` | `2763 MiB` | `5005 MiB` | exact |

Observed scaling:
- kallistors is slower at one and two threads on this subset.
- kallistors becomes faster from four threads onward.
- The largest measured benefit on this subset is at 32 threads: `25.47s -> 17.56s`.

Raw artifacts:
- `target/benchmarks/paired_2026-05-19/paired_1m_multicore_median5.tsv`
- `target/benchmarks/paired_2026-05-19/paired_1m_multicore_summary.tsv`

## Single-End Quant

Measured after the positional fragment-filter parity fixes.

Input:
- Index: `data/gencode.v49_kallisto.idx`
- Reads:
  - `data/subsets/real_mate1_n1048576.fastq.gz`
  - `data/subsets/real_mate1_n2097152.fastq.gz`
- Flags: `--single -l 200 -s 20 -t 32`
- Method: one warmup per binary, alternating measured runs

### Patched Rust vs Upstream Kallisto

| subset | runs | upstream kallisto mean | patched Rust mean | elapsed delta | kallisto RSS | patched RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,048,576 reads | 3 | `27.17s` | `18.69s` | `31%` faster | `2723 MiB` | `5687 MiB` |
| 2,097,152 reads | 2 | `36.35s` | `21.70s` | `40%` faster | `2918 MiB` | `5721 MiB` |

## Index Build

Input:
- Reference transcripts: `data/gencode.v49.transcripts.fa.gz`
- `kallistors index` threads: `-t 8`

| tool | elapsed | peak RSS | index size |
| --- | ---: | ---: | ---: |
| kallisto | `6:24.81` | `12175612KB` | `878M` |
| kallistors | `2:02.23` | `9975360KB` | `894M` |

Result:
- kallistors built the GENCODE index about `3.15x` faster.
- kallistors peak RSS was about `18%` lower for this index-build run.

Use `scripts/bench_index_build.py` to regenerate an index-build report for another machine, FASTA,
or thread count.

## Performance methodology

The [0.5.0 performance report](performance_0.5.0.md) records eight-thread comparisons
against original kallisto, CPU contention, exact EC evidence and the reproducible
`scripts/benchmark_quant_variants.py` workflow.
