# Salmon and piscem comparison

Latest available measurements on the same Apple Silicon machine. Bulk cases use identical
prepared FASTQs, reference FASTA hashes, strand settings and supplied single-end fragment
models. PBMC uses identical prepared lanes, barcode/UMI reads and the common on-list.
The saved metadata comparison is in the [complete numeric evidence](salmon-piscem-results.json).

Tools: Salmon **2.8.0**, piscem **0.23.0**, piscem-infer **0.7.0**, alevin-fry **0.18.3**.
kallistors uses the same frozen executable as the [original-kallisto report](kallisto-comparison-results.md).
Read conditions use one warmup and three measured runs at one and eight requested threads.
Full alevin-fry matrix processing uses eight threads only. Index construction uses eight
threads, one warmup and three builds. Each table reports medians; the numeric evidence
includes individual observations, ranges, CPU/RSS, exact commands and executable hashes.

Paired values below are **kallistors / comparator**. Ratios are comparator wall time divided
by kallistors wall time. Complete pipelines sum sequential stage wall/CPU times within each
repetition before taking the median; pipeline peak RSS is the maximum stage peak.

## Full quantification

GENCODE v50 ALL supplies the muscle, fibroblast and independent simulation reference.
TT-seq uses the official Pachter Lab mature-plus-nascent reference. PE processes both mates;
SE processes R1 with the same supplied Gaussian parameters in each tool. Bias correction
and bootstraps are off. Times include loading, decoding, mapping, native abundance inference
and normal output, including intermediate RAD output for piscem.

Salmon sketch and transcriptome-only selective alignment are separate methods; selective
alignment has no genomic decoys in this comparison. Salmon uses native VB inference;
kallistors and piscem-infer use their native EM/model procedures. These complete workflow
timings do not imply identical mapping or abundance estimates.

### Paired-end

| Reference | Dataset | Comparator | Threads | Wall s | CPU s | Peak RSS GiB | Time ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR13638690 | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 28.450 / 54.141 | 28.420 / 54.070 | 1.724 / 0.659 | 1.90× |
| gencode | SRR13638690 | Salmon 2.8.0 selective alignment | 1 | 28.450 / 212.987 | 28.420 / 214.730 | 1.724 / 2.211 | 7.49× |
| gencode | SRR13638690 | Salmon 2.8.0 sketch | 1 | 28.450 / 77.414 | 28.420 / 79.210 | 1.724 / 0.833 | 2.72× |
| gencode | SRR13638690 | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 6.327 / 19.900 | 48.180 / 126.680 | 2.429 / 0.791 | 3.15× |
| gencode | SRR13638690 | Salmon 2.8.0 selective alignment | 8 | 6.327 / 37.350 | 48.180 / 288.730 | 2.429 / 2.487 | 5.90× |
| gencode | SRR13638690 | Salmon 2.8.0 sketch | 8 | 6.327 / 15.295 | 48.180 / 113.050 | 2.429 / 1.031 | 2.42× |
| gencode | SRR307903 | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 166.315 / 305.066 | 166.200 / 304.780 | 1.815 / 0.708 | 1.83× |
| gencode | SRR307903 | Salmon 2.8.0 selective alignment | 1 | 166.315 / 1003.281 | 166.200 / 1021.910 | 1.815 / 2.403 | 6.03× |
| gencode | SRR307903 | Salmon 2.8.0 sketch | 1 | 166.315 / 441.606 | 166.200 / 459.760 | 1.815 / 1.283 | 2.66× |
| gencode | SRR307903 | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 34.198 / 100.788 | 285.970 / 677.350 | 2.596 / 0.762 | 2.95× |
| gencode | SRR307903 | Salmon 2.8.0 selective alignment | 8 | 34.198 / 176.307 | 285.970 / 1388.250 | 2.596 / 2.594 | 5.16× |
| gencode | SRR307903 | Salmon 2.8.0 sketch | 8 | 34.198 / 82.920 | 285.970 / 644.360 | 2.596 / 1.389 | 2.42× |
| gencode | fragment-truth | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 29.207 / 64.615 | 29.160 / 64.530 | 1.577 / 0.658 | 2.21× |
| gencode | fragment-truth | Salmon 2.8.0 selective alignment | 1 | 29.207 / 128.005 | 29.160 / 130.310 | 1.577 / 2.223 | 4.38× |
| gencode | fragment-truth | Salmon 2.8.0 sketch | 1 | 29.207 / 73.673 | 29.160 / 75.970 | 1.577 / 0.862 | 2.52× |
| gencode | fragment-truth | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 6.947 / 32.648 | 48.800 / 224.900 | 2.313 / 0.772 | 4.70× |
| gencode | fragment-truth | Salmon 2.8.0 selective alignment | 8 | 6.947 / 23.527 | 48.800 / 178.570 | 2.313 / 2.524 | 3.39× |
| gencode | fragment-truth | Salmon 2.8.0 sketch | 8 | 6.947 / 14.595 | 48.800 / 108.000 | 2.313 / 1.081 | 2.10× |
| nascent | SRR10919051 | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 147.233 / 201.880 | 147.070 / 201.750 | 9.157 / 3.132 | 1.37× |
| nascent | SRR10919051 | Salmon 2.8.0 selective alignment | 1 | 147.233 / 794.562 | 147.070 / 795.170 | 9.157 / 5.430 | 5.40× |
| nascent | SRR10919051 | Salmon 2.8.0 sketch | 1 | 147.233 / 274.321 | 147.070 / 275.110 | 9.157 / 3.182 | 1.86× |
| nascent | SRR10919051 | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 34.878 / 45.693 | 227.170 / 342.180 | 9.214 / 3.217 | 1.31× |
| nascent | SRR10919051 | Salmon 2.8.0 selective alignment | 8 | 34.878 / 131.501 | 227.170 / 1039.930 | 9.214 / 5.606 | 3.77× |
| nascent | SRR10919051 | Salmon 2.8.0 sketch | 8 | 34.878 / 47.781 | 227.170 / 372.730 | 9.214 / 3.298 | 1.37× |

### Single-end

| Reference | Dataset | Comparator | Threads | Wall s | CPU s | Peak RSS GiB | Time ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR13638690 | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 19.288 / 34.530 | 19.260 / 34.450 | 1.719 / 0.653 | 1.79× |
| gencode | SRR13638690 | Salmon 2.8.0 selective alignment | 1 | 19.288 / 146.846 | 19.260 / 148.550 | 1.719 / 2.218 | 7.61× |
| gencode | SRR13638690 | Salmon 2.8.0 sketch | 1 | 19.288 / 49.726 | 19.260 / 51.800 | 1.719 / 0.832 | 2.58× |
| gencode | SRR13638690 | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 4.864 / 19.593 | 33.500 / 125.200 | 2.333 / 0.768 | 4.03× |
| gencode | SRR13638690 | Salmon 2.8.0 selective alignment | 8 | 4.864 / 26.348 | 33.500 / 201.100 | 2.333 / 2.425 | 5.42× |
| gencode | SRR13638690 | Salmon 2.8.0 sketch | 8 | 4.864 / 10.170 | 33.500 / 73.440 | 2.333 / 1.001 | 2.09× |
| gencode | SRR307903 | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 99.791 / 154.527 | 99.720 / 154.220 | 1.789 / 0.653 | 1.55× |
| gencode | SRR307903 | Salmon 2.8.0 selective alignment | 1 | 99.791 / 654.130 | 99.720 / 673.500 | 1.789 / 2.352 | 6.55× |
| gencode | SRR307903 | Salmon 2.8.0 sketch | 1 | 99.791 / 230.212 | 99.720 / 253.470 | 1.789 / 0.969 | 2.31× |
| gencode | SRR307903 | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 21.738 / 66.422 | 176.130 / 395.480 | 2.497 / 0.746 | 3.06× |
| gencode | SRR307903 | Salmon 2.8.0 selective alignment | 8 | 21.738 / 117.467 | 176.130 / 922.650 | 2.497 / 2.514 | 5.40× |
| gencode | SRR307903 | Salmon 2.8.0 sketch | 8 | 21.738 / 45.099 | 176.130 / 349.720 | 2.497 / 1.092 | 2.07× |
| gencode | fragment-truth | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 24.095 / 46.091 | 24.050 / 46.010 | 1.761 / 0.659 | 1.91× |
| gencode | fragment-truth | Salmon 2.8.0 selective alignment | 1 | 24.095 / 121.126 | 24.050 / 123.580 | 1.761 / 2.240 | 5.03× |
| gencode | fragment-truth | Salmon 2.8.0 sketch | 1 | 24.095 / 59.339 | 24.050 / 62.390 | 1.761 / 0.841 | 2.46× |
| gencode | fragment-truth | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 5.852 / 31.139 | 41.440 / 214.410 | 2.397 / 0.995 | 5.32× |
| gencode | fragment-truth | Salmon 2.8.0 selective alignment | 8 | 5.852 / 23.784 | 41.440 / 178.610 | 2.397 / 2.481 | 4.06× |
| gencode | fragment-truth | Salmon 2.8.0 sketch | 8 | 5.852 / 12.248 | 41.440 / 88.700 | 2.397 / 1.053 | 2.09× |
| nascent | SRR10919051 | piscem 0.23.0 + piscem-infer 0.7.0 | 1 | 100.709 / 111.429 | 100.590 / 111.310 | 10.205 / 3.125 | 1.11× |
| nascent | SRR10919051 | Salmon 2.8.0 selective alignment | 1 | 100.709 / 645.750 | 100.590 / 646.250 | 10.205 / 5.425 | 6.41× |
| nascent | SRR10919051 | Salmon 2.8.0 sketch | 1 | 100.709 / 150.426 | 100.590 / 151.380 | 10.205 / 3.178 | 1.49× |
| nascent | SRR10919051 | piscem 0.23.0 + piscem-infer 0.7.0 | 8 | 27.130 / 33.322 | 160.990 / 228.290 | 10.263 / 3.189 | 1.23× |
| nascent | SRR10919051 | Salmon 2.8.0 selective alignment | 8 | 27.130 / 117.293 | 160.990 / 927.840 | 10.263 / 5.571 | 4.32× |
| nascent | SRR10919051 | Salmon 2.8.0 sketch | 8 | 27.130 / 27.164 | 160.990 / 208.330 | 10.263 / 3.256 | 1.00× |

The TT-seq single-end eight-thread medians are effectively tied for kallistors and Salmon
sketch (27.130 s and 27.164 s). The README rounds that ratio to 1.00×.

## Native mapping

These commands include index loading, decoding and RAD intermediate output. Bulk Salmon
mapping was measured separately from its quantification command; piscem mapping is the first
stage of its quantification pipeline. PE and SE retain their native fragment/strand policies.
These are not the unstranded independent-read Fulgor/Themisto workloads. The table reports
the mapping commands themselves; full quantification ratios above include inference too.

| Reference | Dataset | Mapper | Threads | Wall s | CPU s | Peak RSS GiB |
| --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR13638690 | piscem 0.23.0 | 1 | 39.924 | 39.880 | 0.659 |
| gencode | SRR13638690 | piscem 0.23.0 | 8 | 6.988 | 54.170 | 0.791 |
| gencode | SRR13638690 | Salmon 2.8.0 selective alignment | 1 | 192.897 | 192.780 | 2.089 |
| gencode | SRR13638690 | Salmon 2.8.0 selective alignment | 8 | 32.587 | 255.860 | 2.398 |
| gencode | SRR13638690 | Salmon 2.8.0 sketch | 1 | 46.849 | 46.800 | 0.704 |
| gencode | SRR13638690 | Salmon 2.8.0 sketch | 8 | 8.477 | 65.510 | 0.884 |
| gencode | SRR13638690-se | piscem 0.23.0 | 1 | 19.584 | 19.540 | 0.653 |
| gencode | SRR13638690-se | piscem 0.23.0 | 8 | 3.999 | 25.210 | 0.768 |
| gencode | SRR13638690-se | Salmon 2.8.0 selective alignment | 1 | 122.340 | 122.250 | 2.082 |
| gencode | SRR13638690-se | Salmon 2.8.0 selective alignment | 8 | 20.784 | 162.710 | 2.344 |
| gencode | SRR13638690-se | Salmon 2.8.0 sketch | 1 | 24.127 | 24.080 | 0.698 |
| gencode | SRR13638690-se | Salmon 2.8.0 sketch | 8 | 4.466 | 33.830 | 0.852 |
| gencode | SRR307903 | piscem 0.23.0 | 1 | 238.658 | 238.460 | 0.657 |
| gencode | SRR307903 | piscem 0.23.0 | 8 | 41.503 | 329.270 | 0.762 |
| gencode | SRR307903 | Salmon 2.8.0 selective alignment | 1 | 825.579 | 825.100 | 2.098 |
| gencode | SRR307903 | Salmon 2.8.0 selective alignment | 8 | 139.795 | 1112.700 | 2.405 |
| gencode | SRR307903 | Salmon 2.8.0 sketch | 1 | 275.002 | 274.800 | 0.703 |
| gencode | SRR307903 | Salmon 2.8.0 sketch | 8 | 48.881 | 387.590 | 0.900 |
| gencode | SRR307903-se | piscem 0.23.0 | 1 | 108.290 | 108.050 | 0.653 |
| gencode | SRR307903-se | piscem 0.23.0 | 8 | 20.021 | 146.490 | 0.746 |
| gencode | SRR307903-se | Salmon 2.8.0 selective alignment | 1 | 516.577 | 516.240 | 2.083 |
| gencode | SRR307903-se | Salmon 2.8.0 selective alignment | 8 | 89.166 | 708.950 | 2.376 |
| gencode | SRR307903-se | Salmon 2.8.0 sketch | 1 | 135.583 | 135.430 | 0.697 |
| gencode | SRR307903-se | Salmon 2.8.0 sketch | 8 | 24.782 | 194.720 | 0.885 |
| gencode | fragment-truth | piscem 0.23.0 | 1 | 34.043 | 34.000 | 0.658 |
| gencode | fragment-truth | piscem 0.23.0 | 8 | 6.152 | 47.490 | 0.772 |
| gencode | fragment-truth | Salmon 2.8.0 selective alignment | 1 | 102.011 | 101.930 | 2.089 |
| gencode | fragment-truth | Salmon 2.8.0 selective alignment | 8 | 17.835 | 138.660 | 2.413 |
| gencode | fragment-truth | Salmon 2.8.0 sketch | 1 | 41.892 | 41.840 | 0.702 |
| gencode | fragment-truth | Salmon 2.8.0 sketch | 8 | 7.726 | 59.560 | 0.870 |
| gencode | fragment-truth-se | piscem 0.23.0 | 1 | 15.580 | 15.540 | 0.659 |
| gencode | fragment-truth-se | piscem 0.23.0 | 8 | 3.215 | 23.340 | 0.995 |
| gencode | fragment-truth-se | Salmon 2.8.0 selective alignment | 1 | 78.827 | 78.750 | 2.083 |
| gencode | fragment-truth-se | Salmon 2.8.0 selective alignment | 8 | 14.620 | 113.450 | 2.358 |
| gencode | fragment-truth-se | Salmon 2.8.0 sketch | 1 | 20.850 | 20.800 | 0.702 |
| gencode | fragment-truth-se | Salmon 2.8.0 sketch | 8 | 3.946 | 29.910 | 0.909 |
| gencode | pbmc_1k_v3 | piscem 0.23.0 | 1 | 549.038 | 548.710 | 0.671 |
| gencode | pbmc_1k_v3 | piscem 0.23.0 | 8 | 101.264 | 798.110 | 0.761 |
| nascent | SRR10919051 | piscem 0.23.0 | 1 | 194.102 | 193.990 | 3.132 |
| nascent | SRR10919051 | piscem 0.23.0 | 8 | 33.942 | 266.680 | 3.217 |
| nascent | SRR10919051 | Salmon 2.8.0 selective alignment | 1 | 786.441 | 786.000 | 5.430 |
| nascent | SRR10919051 | Salmon 2.8.0 selective alignment | 8 | 128.986 | 1023.820 | 5.607 |
| nascent | SRR10919051 | Salmon 2.8.0 sketch | 1 | 267.943 | 267.780 | 3.159 |
| nascent | SRR10919051 | Salmon 2.8.0 sketch | 8 | 45.377 | 357.840 | 3.259 |
| nascent | SRR10919051-se | piscem 0.23.0 | 1 | 103.758 | 103.660 | 3.125 |
| nascent | SRR10919051-se | piscem 0.23.0 | 8 | 19.594 | 135.780 | 3.189 |
| nascent | SRR10919051-se | Salmon 2.8.0 selective alignment | 1 | 638.218 | 637.840 | 5.425 |
| nascent | SRR10919051-se | Salmon 2.8.0 selective alignment | 8 | 114.732 | 910.720 | 5.565 |
| nascent | SRR10919051-se | Salmon 2.8.0 sketch | 1 | 141.855 | 141.760 | 3.154 |
| nascent | SRR10919051-se | Salmon 2.8.0 sketch | 8 | 24.188 | 189.090 | 3.218 |

## Single-cell mapping and matrix processing

All 65,893,590 prepared PBMC pairs are processed. The native mapping comparison uses BUS
output for kallistors and RAD output for piscem. Complete processing adds bustools or
alevin-fry. Barcode correction, multimapper handling and UMI resolution differ; runtime
comparison alone does not establish matrix equality. alevin-fry's downstream stages require
at least two threads, so there is no complete one-thread alevin-fry result.

### FASTQ to matrix

| Comparator | Threads | Wall s | CPU s | Peak RSS GiB | Time ratio |
| --- | --- | --- | --- | --- | --- |
| piscem 0.23.0 + alevin-fry 0.18.3 | 8 | 80.807 / 107.660 | 605.160 / 829.710 | 2.033 / 0.761 | 1.33× |

### Mapping only

| Comparator | Threads | Wall s | CPU s | Peak RSS GiB | Time ratio |
| --- | --- | --- | --- | --- | --- |
| piscem 0.23.0 | 1 | 396.088 / 549.038 | 394.870 / 548.710 | 1.515 / 0.671 | 1.39× |
| piscem 0.23.0 | 8 | 72.644 / 101.264 | 593.420 / 798.110 | 1.683 / 0.761 | 1.39× |

## Index construction

Same normalized transcript reference and eight requested threads. Index formats and native
indexing algorithms differ. Times exclude reference preparation, downloads and compilation;
disk size counts the serialized index rather than scratch or peak temporary disk use.
Each row uses one warmup and three complete builds.

| Reference | Comparator | Wall s | CPU s | Peak RSS GiB | Index GiB | Build-time ratio |
| --- | --- | --- | --- | --- | --- | --- |
| gencode | piscem 0.23.0 | 11.124 / 47.115 | 68.700 / 271.480 | 4.869 / 5.142 | 1.228 / 0.570 | 4.24× |
| gencode | Salmon 2.8.0 | 11.124 / 52.183 | 68.700 / 281.990 | 4.869 / 10.110 | 1.228 / 1.955 | 4.69× |
| nascent | piscem 0.23.0 | 98.687 / 240.466 | 586.860 / 1122.430 | 18.795 / 11.939 | 5.944 / 3.079 | 2.44× |
| nascent | Salmon 2.8.0 | 98.687 / 255.781 | 586.860 / 1129.540 | 18.795 / 15.424 | 5.944 / 5.311 | 2.59× |

## Evidence and reproduction

- [Full numeric evidence and commands](salmon-piscem-results.json): complete quantification,
  native mapping, index and single-cell observations; input matching audit; simulation
  accuracy and thread agreement for Salmon/piscem; source result paths and hashes.
- [kallistors numerical evidence](kallisto-comparison-results.json): current native abundance,
  simulation accuracy, fragment models, bootstrap and thread/replicate comparisons.
- [Tool pins and protocols](benchmark-provenance.json) and
  [tool version lock](../scripts/benchmarks/modern/tools.lock.json).
- [Data preparation and benchmark guide](benchmarks.md).

The `source_records` and `observations` fields identify each measurement's frozen inputs,
executable, method and command. Source paths refer to local evidence workspaces; they are
not download links. Native outputs, raw logs and frozen executables remain in those workspaces.
The independent simulation has known fragment origins under its stated model; real-library
agreement does not establish biological accuracy. No exact cross-tool EC or matrix equivalence
is assumed for Salmon/piscem.
