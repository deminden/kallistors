# Fresh comparison with original kallisto

Apple Silicon; original kallisto 0.52.0 and the same frozen kallistors executable used in the modern-tool panels. All primary values are medians of three runs after one warmup. Builds use eight workers only; paired-end, R1-only and single-cell processing use one and eight workers.

| Reference | Dataset | Mode | Workers | Original, s | kallistors, s | Original/current time | Original RSS, GiB | kallistors RSS, GiB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR1039508 | paired | 1 | 316.463 | 144.017 | 2.197× | 2.449 | 1.277 |
| gencode | SRR1039508 | paired | 8 | 178.462 | 26.207 | 6.810× | 2.890 | 1.341 |
| gencode | SRR1039508 | single | 1 | 1015.560 | 95.291 | 10.657× | 2.401 | 1.468 |
| gencode | SRR1039508 | single | 8 | 269.022 | 19.706 | 13.652× | 2.870 | 1.501 |
| gencode | SRR13638690 | paired | 1 | 102.153 | 40.028 | 2.552× | 2.256 | 1.468 |
| gencode | SRR13638690 | paired | 8 | 46.077 | 8.084 | 5.700× | 2.610 | 1.433 |
| gencode | SRR13638690 | single | 1 | 168.754 | 30.958 | 5.451× | 2.208 | 1.468 |
| gencode | SRR13638690 | single | 8 | 60.969 | 6.492 | 9.392× | 2.622 | 1.368 |
| gencode | SRR307903 | paired | 1 | 532.582 | 233.848 | 2.277× | 2.375 | 1.468 |
| gencode | SRR307903 | paired | 8 | 209.235 | 44.036 | 4.751× | 2.847 | 1.563 |
| gencode | SRR307903 | single | 1 | 1268.719 | 139.010 | 9.127× | 2.324 | 1.468 |
| gencode | SRR307903 | single | 8 | 301.096 | 26.689 | 11.282× | 2.789 | 1.475 |
| gencode | fragment-truth | paired | 1 | 94.581 | 41.113 | 2.301× | 2.313 | 1.277 |
| gencode | fragment-truth | paired | 8 | 58.540 | 8.160 | 7.174× | 2.678 | 1.246 |
| gencode | fragment-truth | single | 1 | 426.244 | 33.701 | 12.648× | 2.327 | 1.468 |
| gencode | fragment-truth | single | 8 | 113.636 | 7.039 | 16.143× | 2.733 | 1.391 |
| gencode | index | index | 8 | 121.797 | 10.986 | 11.087× | 7.608 | 5.253 |
| gencode | pbmc_1k_v3 | single_cell | 1 | 940.636 | 507.248 | 1.854× | 2.165 | 2.033 |
| gencode | pbmc_1k_v3 | single_cell | 8 | 181.735 | 96.163 | 1.890× | 2.617 | 2.033 |
| nascent | SRR10919051 | paired | 1 | 339.392 | 188.524 | 1.800× | 15.225 | 9.102 |
| nascent | SRR10919051 | paired | 8 | 96.170 | 41.547 | 2.315× | 16.215 | 9.154 |
| nascent | SRR10919051 | single | 1 | 340.765 | 137.456 | 2.479× | 15.269 | 10.630 |
| nascent | SRR10919051 | single | 8 | 96.083 | 33.479 | 2.870× | 16.317 | 10.651 |
| nascent | index | index | 8 | 900.543 | 94.799 | 9.499× | 20.214 | 18.818 |

GENCODE v50 ALL is used for cDNA and PBMC; TT-seq uses only the official Pachter lab mature-plus-nascent FASTA assets. Strand settings follow metadata, and both tools receive identical empirical single-end fragment models and fastp-prepared inputs. The independent simulation is adapter-free and untrimmed.

Full EC maps are checked before bulk timing. The primary read comparison uses identical original-built index bytes; both programs additionally quantify matched prefixes using the freshly kallistors-built indexes to check index compatibility. Prefix compatibility timings are not primary performance observations.

The official release executable supplies original timing controls. The separate instrumented original is used only for full EC correctness exports. Every full single-cell run has exact weighted BUS and canonical gene/UMI matrix parity after the same bustools pipeline.

Native HDF5 abundance comparisons include maximum absolute differences, p99, p99.9 and p99.99 at transcript and gene level for every measured replicate. Results include the complete feature universe and the nonzero union, so zeros cannot conceal differences among expressed features. The largest transcript differences are labelled in a separate table. Native fragment sampling and EM can produce abundance differences even when EC counts are exactly equal; realized fragment histograms, within-tool replicate variation and one-versus-eight-worker comparisons are retained separately. Real-data agreement does not establish which estimator is more accurate; the independent simulation supplies stated-model truth.

Single-end inputs use the same supplied Gaussian mean and SD. Kallistors does not serialize that model into `aux/fld`, while original kallisto stores an integerized Gaussian histogram. Its absence is reported as unavailable; no zero histogram is substituted and no histogram equality or L1 difference is asserted. Native count, effective-length and TPM comparisons remain complete. Empirical paired-end histograms are compared on their common one-base bin coordinates.

See the [complete numeric tables](kallisto-comparison-tables/README.md), [machine-readable results](kallisto-comparison-results.json), and [methods](modern-benchmarks.md).
