# Latest Benchmark Summary

The detailed benchmark record now lives in [docs/benchmarks.md](docs/benchmarks.md).

Latest paired-end headline, measured on 2026-05-19 with one warmup and median of five measured
runs on the checked-in full paired dataset:

| tool | threads | median elapsed | median RSS |
| --- | ---: | ---: | ---: |
| kallisto | 32 | `55.81s` | `2784 MiB` |
| kallistors | 32 | `32.79s` | `5504 MiB` |

On the deterministic `1,048,576`-pair subset, kallistors reaches `1.45x` faster than kallisto at
32 threads with exact `run_info.json` parity. The latest full-file paired count check is also exact:
`4,408,640` processed, `4,244,771` pseudoaligned, `276,251` unique. Full multicore tables,
single-end results, index-build results, methodology, and raw artifact paths are in
[docs/benchmarks.md](docs/benchmarks.md).
