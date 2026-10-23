# Modern quantifier benchmark results

Apple Silicon; pinned native binaries; one warmup and three measured runs per condition. All times and peak RSS values below are medians. Index construction uses eight workers only. Bulk and single-cell mapping use one and eight workers; alevin-fry complete processing uses eight only. Cell Ranger execution is deferred.

Read times cover prepared FASTQ to native abundance or count matrix. Downloading, trimming and reference preparation are outside the timers. See the [methods](modern-benchmarks.md) for exact settings, algorithm differences and limitations.

## Index construction

| Reference | Builder | Wall, s | Range, s | CPU, s | Peak RSS, GiB | Index, GiB |
| --- | --- | --- | --- | --- | --- | --- |
| GENCODE v50 ALL | kallistors | 11.845 | 11.540–13.583 | 71.76 | 4.862 | 1.228 |
| GENCODE v50 ALL | piscem | 47.115 | 46.437–47.202 | 271.48 | 5.142 | 0.570 |
| GENCODE v50 ALL | Salmon | 52.183 | 50.580–52.999 | 281.99 | 10.110 | 1.955 |
| Official mature + nascent | kallistors | 98.732 | 97.845–99.528 | 590.42 | 18.013 | 5.944 |
| Official mature + nascent | piscem | 240.466 | 240.300–241.950 | 1122.43 | 11.939 | 3.079 |
| Official mature + nascent | Salmon | 255.781 | 255.683–256.310 | 1129.54 | 15.424 | 5.311 |

## Paired-end quantification

| Dataset | Method | Wall at 1, s | Wall at 8, s | RSS at 1, GiB | RSS at 8, GiB |
| --- | --- | --- | --- | --- | --- |
| SRR13638690 | kallistors | 44.814 | 9.215 | 1.423 | 1.446 |
| SRR13638690 | Salmon sketch | 77.414 | 15.295 | 0.833 | 1.031 |
| SRR13638690 | Salmon selective alignment | 212.987 | 37.350 | 2.211 | 2.487 |
| SRR13638690 | piscem + piscem-infer | 54.141 | 19.900 | 0.659 | 0.791 |
| SRR307903 | kallistors | 275.460 | 50.990 | 1.423 | 1.501 |
| SRR307903 | Salmon sketch | 441.606 | 82.920 | 1.283 | 1.389 |
| SRR307903 | Salmon selective alignment | 1003.281 | 176.307 | 2.403 | 2.594 |
| SRR307903 | piscem + piscem-infer | 305.066 | 100.788 | 0.708 | 0.762 |
| fragment-truth | kallistors | 52.009 | 10.515 | 1.248 | 1.244 |
| fragment-truth | Salmon sketch | 73.673 | 14.595 | 0.862 | 1.081 |
| fragment-truth | Salmon selective alignment | 128.005 | 23.527 | 2.223 | 2.524 |
| fragment-truth | piscem + piscem-infer | 64.615 | 32.648 | 0.658 | 0.772 |
| SRR10919051 | kallistors | 445.947 | 88.798 | 8.735 | 8.742 |
| SRR10919051 | Salmon sketch | 274.321 | 47.781 | 3.182 | 3.298 |
| SRR10919051 | Salmon selective alignment | 794.562 | 131.501 | 5.430 | 5.606 |
| SRR10919051 | piscem + piscem-infer | 201.880 | 45.693 | 3.132 | 3.217 |

## Single-end quantification

| Dataset | Method | Wall at 1, s | Wall at 8, s | RSS at 1, GiB | RSS at 8, GiB |
| --- | --- | --- | --- | --- | --- |
| SRR13638690-se | kallistors | 33.899 | 7.271 | 1.423 | 1.445 |
| SRR13638690-se | Salmon sketch | 49.726 | 10.170 | 0.832 | 1.001 |
| SRR13638690-se | Salmon selective alignment | 146.846 | 26.348 | 2.218 | 2.425 |
| SRR13638690-se | piscem + piscem-infer | 34.530 | 19.593 | 0.653 | 0.768 |
| SRR307903-se | kallistors | 159.997 | 30.474 | 1.423 | 1.446 |
| SRR307903-se | Salmon sketch | 230.212 | 45.099 | 0.969 | 1.092 |
| SRR307903-se | Salmon selective alignment | 654.130 | 117.467 | 2.352 | 2.514 |
| SRR307903-se | piscem + piscem-infer | 154.527 | 66.422 | 0.653 | 0.746 |
| fragment-truth-se | kallistors | 39.081 | 8.403 | 1.423 | 1.446 |
| fragment-truth-se | Salmon sketch | 59.339 | 12.248 | 0.841 | 1.053 |
| fragment-truth-se | Salmon selective alignment | 121.126 | 23.784 | 2.240 | 2.481 |
| fragment-truth-se | piscem + piscem-infer | 46.091 | 31.139 | 0.659 | 0.995 |
| SRR10919051-se | kallistors | 293.369 | 60.960 | 9.888 | 9.910 |
| SRR10919051-se | Salmon sketch | 150.426 | 27.164 | 3.178 | 3.256 |
| SRR10919051-se | Salmon selective alignment | 645.750 | 117.293 | 5.425 | 5.571 |
| SRR10919051-se | piscem + piscem-infer | 111.429 | 33.322 | 3.125 | 3.189 |

Single-end observations use R1 from the paired library. All methods receive the same frozen paired-end fragment mean and SD for each real library. Simulation fragment lengths are known by construction. These layouts are related inputs, not biological replicates.

## Single-cell processing

| Pipeline | Workers | Wall, s | CPU, s | Peak RSS, GiB | Cells passing gate | Median UMIs | Median detected genes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kallistors + bustools | 1 | 811.331 | 810.74 | 2.033 | 1202 | 6957.5 | 2376.5 |
| kallistors + bustools | 8 | 151.380 | 1166.02 | 2.033 | 1202 | 6957.5 | 2376.5 |
| piscem + alevin-fry | 8 | 107.660 | 829.71 | 0.761 | 1202 | 7056.0 | 2390.5 |

The common descriptive gate is at least 500 UMIs and 200 detected genes. Native retained barcodes are not a biological cell count. Barcode correction, multimapper handling and UMI resolution differ between bustools and alevin-fry.

At eight workers, **1,202 cells** pass the gate in both workflows; called-cell Jaccard is **1.00000**. Matrix differences are reported over the full barcode/gene union and the common-called subset.

| Matrix entries | p99 | p99.9 | p99.99 | Maximum absolute difference |
| --- | --- | --- | --- | --- |
| Full nonzero union | 1 | 4 | 11 | 176 |
| Common-called nonzero union | 1 | 4 | 12 | 176 |

## Independent simulation accuracy

The following summaries use measured replicate 1 at eight workers, selected before observing performance. Known fragment origins and nominal molecule TPM come from the independent simulation. Full percentile, correlation and one-worker results are in the supplementary tables.

| Input | Method | Level | Count MAE | Count max abs. | TPM MAE | TPM log1p Pearson |
| --- | --- | --- | --- | --- | --- | --- |
| fragment-truth | kallistors | transcript | 0.938146 | 2188.55 | 0.303404 | 0.890632 |
| fragment-truth | kallistors | gene | 1.21728 | 1728.52 | 1.18659 | 0.97233 |
| fragment-truth | Salmon sketch | transcript | 0.873082 | 2779.21 | 0.292254 | 0.899683 |
| fragment-truth | Salmon sketch | gene | 1.242 | 3249 | 1.18983 | 0.973828 |
| fragment-truth | Salmon selective alignment | transcript | 0.852321 | 4278.19 | 0.282035 | 0.903633 |
| fragment-truth | Salmon selective alignment | gene | 1.226 | 4464.95 | 1.18464 | 0.9741 |
| fragment-truth | piscem + piscem-infer | transcript | 0.904422 | 1987.53 | 0.300359 | 0.892544 |
| fragment-truth | piscem + piscem-infer | gene | 1.24242 | 1725.69 | 1.20099 | 0.971922 |
| fragment-truth-se | kallistors | transcript | 1.44222 | 5621.5 | 0.439983 | 0.814379 |
| fragment-truth-se | kallistors | gene | 1.37705 | 1738.56 | 1.33327 | 0.967726 |
| fragment-truth-se | Salmon sketch | transcript | 1.74601 | 6859.58 | 1.09037 | 0.729722 |
| fragment-truth-se | Salmon sketch | gene | 1.84966 | 5439.87 | 4.93367 | 0.939328 |
| fragment-truth-se | Salmon selective alignment | transcript | 1.71575 | 7751.47 | 1.07541 | 0.734051 |
| fragment-truth-se | Salmon selective alignment | gene | 1.81935 | 7976.65 | 4.90725 | 0.94053 |
| fragment-truth-se | piscem + piscem-infer | transcript | 1.369 | 6871.12 | 0.658199 | 0.79511 |
| fragment-truth-se | piscem + piscem-infer | gene | 1.35632 | 1791.79 | 2.91013 | 0.956326 |

This simulation omits indels, polymorphism, contamination, PCR and sequence biases. It cannot establish a universal biological accuracy ranking. Real-library and single-cell comparisons measure agreement only.

## Mapping, agreement and reproducibility

* Separate native mapping tables include index loading, decompression and intermediate output. Aggregate EC, RAD and BUS outputs have different workloads; these timings are not pure mapping-kernel timings.
* Full pairwise agreement, replicate variability and within-method one-versus-eight-worker agreement include maximum absolute differences, p99, p99.9 and p99.99. Gene-by-state and mature/nascent subsets are included for TT-seq.
* Wall-time and RSS ranges, CPU time, pipeline stages, retained-barcode and UMI summaries, and processed-input audits are in the [publication tables](modern-benchmark-tables/README.md).
* Native processed-record counters match each frozen input count. These are technical triplicates on one host and two distinct references; no cross-host or biological replication is implied.
* Salmon selective alignment uses a transcript-only reference without genomic decoys. The TT-seq reference is the official Ensembl 108 mature-plus-nascent set, not an augmented GENCODE v50 set.
