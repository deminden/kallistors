# Current kallistors versus original kallisto

Apple Silicon, native kallistors release build and original kallisto 0.52.0.

Fresh alternating controls, one warmup and three measured runs. Quantification and single-cell pipelines: 1 and 8 requested workers. Index construction: 8 only. All completed timing observations retained, with contention in raw records.

## Published evidence

This page contains only the latest completed comparison. The frozen measurement commit is
`ff74fddf928322d40568ae858ac245859eb4f499`; after the unpublished history was reorganized,
the equivalent code checkpoint is `2054ff7`. The executable hashes remain authoritative.

- [Complete numeric results](kallisto-comparison-results.json): all 76 performance records,
  individual wall times, CPU/RSS, contention, native abundance comparisons, fragment models,
  simulation accuracy, bootstrap summaries, within-tool variation, thread agreement and largest differences.
- [Protocols and provenance](benchmark-provenance.json): original command configuration,
  prepared input/tool hashes, compiler settings, host information, EC checks, full PBMC parity,
  index compatibility, completion records and recovery amendment.
- [Preparation and runner guide](benchmarks.md).

Unlinked evidence filenames below refer to the frozen local workspace
`target/benchmarks/original-ff74fdd.noindex/`. Large native outputs and logs remain there;
the published JSON contains the complete numerical summaries. Absolute paths in metadata
identify those records and are not download links.

## Inputs and modes

GENCODE v50 ALL supplies the bulk cDNA and PBMC reference; TT-seq uses the official
Pachter Lab human mature-plus-nascent reference. Real bulk libraries were prepared with
fastp 1.4.0. PE processes both mates; SE processes R1 with the identical supplied Gaussian
model in both tools. The table rounds calibration parameters for readability; the protocol
retains their full precision. Paired fragment distributions are estimated during each run.

| Reference | Dataset | Prepared pairs | Strand | Supplied SE mean / SD, bases |
| --- | --- | ---: | --- | --- |
| gencode | SRR13638690 | 4,411,003 | reverse | 144.2279 / 50.1407 |
| gencode | SRR1039508 | 22,176,801 | unstranded | 144.5484 / 55.5962 |
| gencode | SRR307903 | 32,460,465 | reverse | 190.1364 / 22.0902 |
| nascent | SRR10919051 | 10,172,656 | unstranded | 135.6360 / 39.3825 |

The independent `fragment-truth` simulation contains 4,000,000 adapter-free pairs and is
untrimmed. Its generation parameters and truth files are pinned in the protocol. The full
10x PBMC 1k v3 input has 65,893,590 pairs across two lanes; retained R1 barcode/UMI sequences
are unchanged after cDNA filtering. BUS mapping and the complete FASTQ-to-matrix pipeline are
reported separately, with identical bustools processing for both tools.

Bias and bootstrap panels use the full muscle library in paired and single-end modes.
Quantification includes index loading, decompression, mapping, fragment/model work, EM and
normal TSV/HDF5 output. Reference preparation and read trimming are outside these timers.
The bulk exact-match cache uses the automatic budget; its memory is included in peak RSS.

## Quantification, single-cell and index performance

Paired values are **kallistors / original**. CPU is user plus system time. All rows use three measured runs after one warmup.

| Reference | Dataset | Mode | Threads | Scope | Wall s · current/original | CPU s · current/original | Peak RSS GiB · current/original | Original/current wall |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR13638690 | paired | 1 | default_quantification | 28.450 / 101.367 | 28.420 / 101.270 | 1.724 / 2.219 | 3.56× |
| gencode | SRR13638690 | paired | 8 | default_quantification | 6.327 / 46.112 | 48.180 / 141.970 | 2.429 / 2.615 | 7.29× |
| gencode | SRR13638690 | single | 1 | default_quantification | 19.288 / 164.560 | 19.260 / 164.460 | 1.719 / 2.294 | 8.53× |
| gencode | SRR13638690 | single | 8 | default_quantification | 4.864 / 57.608 | 33.500 / 231.450 | 2.333 / 2.622 | 11.84× |
| gencode | SRR13638690 | paired | 1 | bootstrap-1 | 32.301 / 124.643 | 32.270 / 124.560 | 1.724 / 2.280 | 3.86× |
| gencode | SRR13638690 | paired | 8 | bootstrap-1 | 7.398 / 70.087 | 55.940 / 166.810 | 2.436 / 2.663 | 9.47× |
| gencode | SRR13638690 | single | 1 | bootstrap-1 | 23.885 / 189.723 | 23.850 / 189.600 | 1.719 / 2.273 | 7.94× |
| gencode | SRR13638690 | single | 8 | bootstrap-1 | 6.091 / 83.080 | 41.890 / 256.320 | 2.327 / 2.680 | 13.64× |
| gencode | SRR13638690 | paired | 1 | bootstrap-30 | 154.086 / 879.608 | 153.980 / 878.990 | 1.724 / 2.281 | 5.71× |
| gencode | SRR13638690 | paired | 8 | bootstrap-30 | 35.079 / 303.662 | 228.050 / 2069.200 | 2.423 / 3.130 | 8.66× |
| gencode | SRR13638690 | single | 1 | bootstrap-30 | 147.977 / 934.416 | 147.870 / 933.870 | 1.719 / 2.294 | 6.31× |
| gencode | SRR13638690 | single | 8 | bootstrap-30 | 34.355 / 313.827 | 218.040 / 2159.990 | 2.336 / 3.086 | 9.13× |
| gencode | SRR1039508 | paired | 1 | default_quantification | 102.197 / 315.289 | 102.130 / 315.100 | 1.647 / 2.445 | 3.09× |
| gencode | SRR1039508 | paired | 8 | default_quantification | 22.149 / 180.992 | 179.440 / 408.000 | 2.418 / 2.890 | 8.17× |
| gencode | SRR1039508 | single | 1 | default_quantification | 73.795 / 1011.023 | 73.740 / 1010.470 | 1.800 / 2.399 | 13.70× |
| gencode | SRR1039508 | single | 8 | default_quantification | 17.105 / 263.418 | 133.230 / 1361.930 | 2.503 / 2.870 | 15.40× |
| gencode | SRR307903 | paired | 1 | default_quantification | 166.315 / 536.588 | 166.200 / 536.110 | 1.815 / 2.417 | 3.23× |
| gencode | SRR307903 | paired | 8 | default_quantification | 34.198 / 207.465 | 285.970 / 751.780 | 2.596 / 2.843 | 6.07× |
| gencode | SRR307903 | single | 1 | default_quantification | 99.791 / 1267.382 | 99.720 / 1266.500 | 1.789 / 2.322 | 12.70× |
| gencode | SRR307903 | single | 8 | default_quantification | 21.738 / 301.252 | 176.130 / 1777.380 | 2.497 / 2.791 | 13.86× |
| gencode | fragment-truth | paired | 1 | default_quantification | 29.207 / 94.228 | 29.160 / 94.140 | 1.577 / 2.346 | 3.23× |
| gencode | fragment-truth | paired | 8 | default_quantification | 6.947 / 57.670 | 48.800 / 125.490 | 2.313 / 2.675 | 8.30× |
| gencode | fragment-truth | single | 1 | default_quantification | 24.095 / 427.155 | 24.050 / 426.920 | 1.761 / 2.286 | 17.73× |
| gencode | fragment-truth | single | 8 | default_quantification | 5.852 / 114.000 | 41.440 / 568.960 | 2.397 / 2.672 | 19.48× |
| nascent | SRR10919051 | paired | 1 | default_quantification | 147.233 / 343.073 | 147.070 / 342.520 | 9.157 / 14.849 | 2.33× |
| nascent | SRR10919051 | paired | 8 | default_quantification | 34.878 / 97.173 | 227.170 / 521.120 | 9.214 / 16.085 | 2.79× |
| nascent | SRR10919051 | single | 1 | default_quantification | 100.709 / 344.160 | 100.590 / 343.840 | 10.205 / 15.147 | 3.42× |
| nascent | SRR10919051 | single | 8 | default_quantification | 27.130 / 96.043 | 160.990 / 519.020 | 10.263 / 15.999 | 3.54× |
| gencode | pbmc_1k_v3 | single_cell | 1 | FASTQ_to_matrix | 407.147 / 940.299 | 405.730 / 938.980 | 2.033 / 2.158 | 2.31× |
| gencode | pbmc_1k_v3 | single_cell | 1 | BUS_mapping | 396.088 / 929.294 | 394.870 / 928.180 | 1.515 / 2.158 | 2.35× |
| gencode | pbmc_1k_v3 | single_cell | 8 | FASTQ_to_matrix | 80.807 / 182.367 | 605.160 / 1368.490 | 2.033 / 2.585 | 2.26× |
| gencode | pbmc_1k_v3 | single_cell | 8 | BUS_mapping | 72.644 / 174.286 | 593.420 / 1356.800 | 1.683 / 2.585 | 2.40× |
| gencode | index | index | 8 | index_build | 11.124 / 121.488 | 68.700 / 609.130 | 4.869 / 7.689 | 10.92× |
| nascent | index | index | 8 | index_build | 98.687 / 907.842 | 586.860 / 4303.310 | 18.795 / 20.152 | 9.20× |

Native HDF5 comparisons report full-target and nonzero-union max, p99, p99.9, p99.99, transcript and gene counts/TPM. Real-data agreement does not establish biological accuracy. Fragment models are recorded separately; paired empirical samples may vary with concurrency. Bootstrap draws use different RNGs across programs and are compared as distributions.

Timing uses the uninstrumented original v0.52.0 release. The separate instrumented original exports exact EC evidence only. Both timing tools consume identical original-built v13 index bytes and identical prepared FASTQs. Newly built candidate indexes receive a separate cross-tool compatibility check. Bias correction uses the same FASTA; general bias-estimate equivalence is not assumed.

Inputs, tool hashes, compiler settings, source, workflow and analysis dependency versions are pinned in `protocol.json`. Lossless compressed TSV outputs retain full SHA-verified content; HDF5 remains directly readable.

## Index construction

Eight requested threads only, one warmup and three complete measured builds per tool and
reference. CPU is user plus system time; disk size counts the serialized index, excluding
inputs, build scratch and tool binaries. The main performance table reports median wall time.

| Reference | Tool | Measured build times, s | Median CPU, s | Median peak RSS, GiB | Index size, GiB |
| --- | --- | --- | ---: | ---: | ---: |
| gencode | kallisto | 120.924, 121.734, 121.488 | 609.130 | 7.689 | 1.190 |
| gencode | kallistors | 11.124, 10.862, 11.136 | 68.700 | 4.869 | 1.228 |
| nascent | kallisto | 905.942, 909.478, 907.842 | 4303.310 | 20.152 | 6.132 |
| nascent | kallistors | 97.333, 98.687, 99.509 | 586.860 | 18.795 | 5.944 |

Original/kallistors median build-time ratios are **10.92× on GENCODE** and **9.20× on
mature + nascent**. The [combined construction table](../README.md#index-construction)
also compares Salmon, piscem, Fulgor and Themisto on the same host, reference target sets and
eight-thread budget. Salmon and piscem use one warmup and three builds; Fulgor and Themisto
have one successful setup observation per available reference. See the
[Salmon/piscem build details](salmon-piscem-results.md#index-construction) for their measurements.

Fresh kallistors-built indexes passed separate cross-tool EC compatibility checks. Exact file
bytes need not match between builders. Primary read timings use identical original-built
indexes in both quantifiers. Fulgor and Themisto setup observations and their distinct index
formats are documented in the [mapping report](mapping-comparison-results.md#setup-at-eight-requested-threads).

## Native abundance agreement

Every ordinary bulk EC preflight is exact, and every full PBMC trial has exact raw weighted
BUS and corrected gene/UMI matrix parity. That does not imply equal native bulk abundance
estimates. Below are the largest absolute differences and worst p99 across the three measured
comparisons in each default condition, using the union of expressed targets. Values are
absolute counts or TPM, not percentages. Simulation agreement here is between programs;
accuracy against known truth is recorded separately in `simulation_accuracy` in the numeric results.

| Reference | Dataset | Mode | Threads | Transcript count max | Transcript count p99 | Transcript TPM max | Gene count max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR1039508 | paired | 1 | 1.2551e-10 | 5.12372e-13 | 5.57066e-12 | 5.82077e-11 |
| gencode | SRR1039508 | paired | 8 | 5.27716 | 0.0384023 | 108.807 | 1.40002 |
| gencode | SRR1039508 | single | 1 | 1.09139e-10 | 7.10543e-13 | 5.22959e-12 | 2.91038e-11 |
| gencode | SRR1039508 | single | 8 | 1.45519e-10 | 1.02318e-12 | 1.09139e-11 | 8.73115e-11 |
| gencode | SRR13638690 | paired | 1 | 1.62981e-09 | 1.13687e-13 | 2.18279e-11 | 2.91038e-11 |
| gencode | SRR13638690 | paired | 8 | 1.42852 | 0.0309496 | 65.3055 | 0.455951 |
| gencode | SRR13638690 | single | 1 | 1.28057e-09 | 1.13687e-13 | 4.72937e-11 | 2.32831e-10 |
| gencode | SRR13638690 | single | 8 | 4.36557e-10 | 1.98952e-13 | 3.27418e-11 | 2.32831e-10 |
| gencode | SRR307903 | paired | 1 | 2.83762e-10 | 6.82121e-13 | 8.6402e-12 | 6.91216e-11 |
| gencode | SRR307903 | paired | 8 | 5.87953 | 0.0704155 | 7.88487 | 1.55943 |
| gencode | SRR307903 | single | 1 | 1.14596e-10 | 9.09495e-13 | 1.31877e-11 | 1.16415e-10 |
| gencode | SRR307903 | single | 8 | 2.76486e-10 | 1.19371e-12 | 1.04592e-11 | 1.74623e-10 |
| gencode | fragment-truth | paired | 1 | 4.77485e-11 | 2.27374e-13 | 6.0254e-12 | 1.45519e-11 |
| gencode | fragment-truth | paired | 8 | 0.0145634 | 0.000926203 | 0.0976727 | 0.0071329 |
| gencode | fragment-truth | single | 1 | 5.36602e-11 | 2.27374e-13 | 6.59384e-12 | 7.27596e-12 |
| gencode | fragment-truth | single | 8 | 4.54747e-11 | 3.69482e-13 | 7.27596e-12 | 1.45519e-11 |
| nascent | SRR10919051 | paired | 1 | 4.65661e-10 | 5.68434e-14 | 1.45519e-11 | 4.65661e-10 |
| nascent | SRR10919051 | paired | 8 | 0.775836 | 0.00134368 | 22.5027 | 0.775836 |
| nascent | SRR10919051 | single | 1 | 4.65661e-10 | 5.68434e-14 | 1.45519e-11 | 4.65661e-10 |
| nascent | SRR10919051 | single | 8 | 4.65661e-10 | 1.13687e-13 | 1.45519e-11 | 4.65661e-10 |

The [complete results](kallisto-comparison-results.json) retain all-feature and nonzero-union
tails through p99.99, effective lengths, gene TPM, labelled largest differences, bootstrap
comparisons, and separate within-program and 1-versus-8-thread variation. No difference is
removed by averaging replicates or omitting zero-count targets. Real-data agreement alone
does not identify the more accurate estimator.

Single-end runs use the same supplied Gaussian mean and SD. kallistors does not serialize
this model into `aux/fld`; original stores an integerized histogram. The missing serialized
histogram is recorded as unavailable, not replaced with zeros or treated as a mismatch.
Paired empirical histograms are compared on shared one-base bin coordinates.

## Completed sequence-bias observations

These completed runs use materially different native abundance estimates. Timing ratios below describe the measured commands, do not establish equivalent-output acceleration, and are excluded from ordinary-mode headline comparisons. No additional original-kallisto bias run was performed during recovery. Elapsed time alone does not establish the cause of the cost.

| Mode | Requested threads | Tool | Median wall s | Wall range s | Median CPU s | Median peak RSS GiB | Background activity |
| --- | ---: | --- | ---: | --- | ---: | ---: | --- |
| paired | 1 | kallisto | 1843.882 | 1829.282–1852.891 | 1842.880 | 3.843 | flagged |
| paired | 1 | kallistors | 57.730 | 57.704–57.852 | 57.670 | 2.508 | none recorded |
| paired | 8 | kallisto | 1785.191 | 1769.350–1794.826 | 1881.540 | 3.976 | flagged |
| paired | 8 | kallistors | 13.673 | 13.640–13.778 | 88.640 | 2.835 | flagged |
| single | 1 | kallisto | 1893.624 | 1888.041–1895.051 | 1892.600 | 3.852 | flagged |
| single | 1 | kallistors | 44.976 | 44.870–44.989 | 44.930 | 2.490 | none recorded |
| single | 8 | kallisto | 1787.291 | 1771.874–1793.900 | 1958.150 | 4.246 | flagged |
| single | 8 | kallistors | 11.952 | 11.888–12.006 | 70.430 | 2.744 | none recorded |

One warmup and three measured runs per condition. Raw warmup and measurement observations remain in the original extra-mode `runs.json` files. All completed pairs are retained.

- paired, 1 requested threads: measured wall-time ratio 31.94×.
- paired, 8 requested threads: measured wall-time ratio 130.56×.
- single, 1 requested threads: measured wall-time ratio 42.10×.
- single, 8 requested threads: measured wall-time ratio 149.53×.

### Native abundance differences

The table summarizes the largest observed absolute difference and the worst replicate percentile across all three paired comparisons, using the nonzero union. Complete transcript and gene tails for both all-feature and nonzero-union universes are retained in `sequence-bias-native.tsv`; no native value or completed observation has been changed.

| Mode | Threads | Metric | Max absolute | p99 | p99.9 | p99.99 |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| paired | 1 | effective_length | 3811.86513 | 1.80257479 | 10.3087114 | 17.9773842 |
| paired | 1 | gene_counts | 3056.51507 | 0.0861496489 | 2.17150471 | 993.099576 |
| paired | 1 | gene_tpm | 2593.93237 | 0.321126294 | 2.76719433 | 850.513475 |
| paired | 1 | transcript_counts | 3057.37159 | 0.92442571 | 10.6698745 | 348.163351 |
| paired | 1 | transcript_tpm | 2594.86522 | 0.306359997 | 2.83868668 | 31.3005895 |
| paired | 8 | effective_length | 3811.1475 | 15.6850867 | 43.7271861 | 155.451483 |
| paired | 8 | gene_counts | 2988.90712 | 0.57422241 | 8.34983725 | 983.271722 |
| paired | 8 | gene_tpm | 2544.23353 | 2.40948655 | 17.3353236 | 897.81329 |
| paired | 8 | transcript_counts | 2989.01189 | 4.5031259 | 37.5352337 | 597.479889 |
| paired | 8 | transcript_tpm | 2544.34772 | 1.73703941 | 14.0927382 | 78.0646706 |
| single | 1 | effective_length | 1422.75727 | 2.05692218 | 10.1999003 | 16.2991224 |
| single | 1 | gene_counts | 3.39241672 | 0.083678568 | 0.926942169 | 2.93694739 |
| single | 1 | gene_tpm | 6.05991386 | 0.289462726 | 2.06897896 | 5.46758142 |
| single | 1 | transcript_counts | 164.974724 | 0.586400403 | 6.89992041 | 66.7487182 |
| single | 1 | transcript_tpm | 15.9413266 | 0.26360678 | 2.22833714 | 6.48333964 |
| single | 8 | effective_length | 4256.43007 | 13.4055656 | 36.1298888 | 125.76639 |
| single | 8 | gene_counts | 88.367629 | 0.44352946 | 6.52174314 | 46.2096154 |
| single | 8 | gene_tpm | 449.458072 | 1.75438966 | 12.7868688 | 139.433296 |
| single | 8 | transcript_counts | 1513.07836 | 2.94766186 | 26.7063601 | 329.016468 |
| single | 8 | transcript_tpm | 449.458072 | 1.3470597 | 8.95427132 | 61.67355 |

## Storage interruption and recovery

The original coordinator stopped at its 8 GiB free-space guard before the first index build. Bulk and full PBMC observations were complete. Duplicate SHA-256-identical historical canonical BUS files were consolidated using hardlinks; every evidence path and byte was retained. The frozen coordinator, protocol, executables, source and completed native observations were not edited. The recovery runner restores the saved analysis arrays and invokes only the original index panels, followed by the original final provenance checks. See `recovery/amendment.json` and `recovery/evidence-storage.json`.
