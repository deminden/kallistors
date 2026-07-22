# Benchmarks

This page records benchmark results for a local real dataset. The inputs and raw artifacts are not
part of the repository. Results are machine- and cache-sensitive; compare only runs measured with
the same inputs, binaries, and timing method.

## Environment

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

### Patched Rust vs Original Rust

The original Rust binary was a clean release build from commit
`82008709d72d5975d997bf05935767d2d63d0308`.

| subset | runs | original Rust mean | patched Rust mean | elapsed delta | original RSS | patched RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,048,576 reads | 5 | `20.10s` | `20.61s` | `+2.5%` slower | `5004 MiB` | `5693 MiB` |
| 2,097,152 reads | 3 | `23.40s` | `23.37s` | `0.1%` faster | `5020 MiB` | `5702 MiB` |

The speed impact is effectively neutral at larger prefixes, with about `680-690 MiB` more peak RSS
from loading full positional sets.

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
