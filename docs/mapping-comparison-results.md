# Mapping comparison: kallistors, Fulgor and Themisto

Apple Silicon using native release builds. Frozen kallistors commit `ff74fddf928322d40568ae858ac245859eb4f499`.

**Coverage:** all five GENCODE read files completed with all three tools. Both official mature-plus-nascent read files completed with kallistors and Themisto. Fulgor nascent is unavailable: its pinned GGCAT builder aborted before producing an index. This limitation is retained as a setup failure, not a fabricated mapping time.

## Published evidence

This page contains only the latest completed comparison. The frozen measurement commit is
`ff74fddf928322d40568ae858ac245859eb4f499`; after the unpublished history was reorganized,
the equivalent code checkpoint is `2054ff7`. The executable hashes remain authoritative.

- [Complete numeric results](mapping-comparison-results.json): all 48 tool conditions,
  individual wall times, CPU/RSS, observed threads, contention, throughput, all agreement
  comparisons, setup time/disk and completion checks.
- [Exact native commands](mapping-comparison-commands.json).
- [Portable command templates](mapping-reproduction.md): configure tool, index, read and
  scratch paths for another checkout while retaining the native command settings.
- [Protocols and provenance](benchmark-provenance.json): input/reference hashes, recursive
  source pins, dependency/compiler settings, executable hashes and recorded setup failures.

Unlinked evidence filenames below refer to the frozen local workspace
`target/benchmarks/mapping-fulgor-themisto-ff74fdd.noindex/`. Large native outputs and logs
remain there; the published JSON contains the complete numerical summaries. Absolute paths
in metadata identify those records and are not download links.

## What was measured

One warmup plus median of three measurements at 1 and 8 requested threads, serial execution under the shared native-measurement lock. Each competitor has fresh alternating kallistors controls. All completed trials are retained, including contention flags; there is no rejection or best-run selection. There are 192 timing commands (144 measured, 48 warmups), 38 full-output correctness commands, 48 tool-condition summaries and seven complete read-file panels.

Wall time includes startup, index loading, decompression and unavoidable native output work, using `/dev/null` in every timing command. kallistors aggregates weighted equivalence classes; Fulgor and Themisto enumerate and format per-read transcript sets even to the null sink. These are descriptive native CLI timings, **not pure mapping-kernel or equivalent-output speedups**. No EM, bias model, native paired-fragment concordance, BUS, barcode/UMI processing or count matrix is timed in this extension.

R1 and R2 are independent single-read passes. Summing mates is not a native paired-end result. PBMC is the full 65,893,590 cDNA R2 reads from both prepared lanes, combined as unchanged concatenated gzip members. The fixture verified concatenated-member consumption; all retained full Fulgor/Themisto outputs were checked for a complete unique zero-based read-ID permutation. Bulk muscle has 4,411,003 reads per mate; fibroblast 32,460,465 per mate; nascent TT-seq 10,172,656 per mate.

## Mapping command medians

Each ratio is competitor wall time / its fresh kallistors control. The algorithms and native outputs differ. CPU/RSS, exact samples, ranges, reads/s and contention counts for every tool condition are in `mapping-performance.tsv`; all raw records and exact commands are in `native-observations.json` and `commands.json`.

Paired values are **kallistors / comparator**. CPU is user plus system time.

| Reference | Read file | Comparator | Requested threads | Wall s · current/comparator | CPU s · current/comparator | Peak RSS GiB · current/comparator | Time ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR13638690-R1 | fulgor | 1 | 12.903 / 28.709 | 12.870 / 43.460 | 1.572 / 0.366 | 2.22× |
| gencode | SRR13638690-R1 | fulgor | 8 | 3.387 / 7.207 | 22.090 / 54.980 | 2.152 / 0.402 | 2.13× |
| gencode | SRR13638690-R1 | themisto | 1 | 12.919 / 85.853 | 12.900 / 89.740 | 1.572 / 1.290 | 6.65× |
| gencode | SRR13638690-R1 | themisto | 8 | 3.389 / 17.217 | 22.100 / 128.330 | 2.164 / 1.228 | 5.08× |
| gencode | SRR13638690-R2 | fulgor | 1 | 13.325 / 29.435 | 13.300 / 44.620 | 1.573 / 0.367 | 2.21× |
| gencode | SRR13638690-R2 | fulgor | 8 | 3.495 / 7.275 | 22.540 / 55.570 | 2.154 / 0.403 | 2.08× |
| gencode | SRR13638690-R2 | themisto | 1 | 13.365 / 90.937 | 13.340 / 94.520 | 1.573 / 1.273 | 6.80× |
| gencode | SRR13638690-R2 | themisto | 8 | 3.499 / 17.234 | 22.450 / 134.240 | 2.162 / 1.211 | 4.93× |
| gencode | SRR307903-R1 | fulgor | 1 | 68.551 / 148.034 | 68.500 / 225.380 | 1.637 / 0.358 | 2.16× |
| gencode | SRR307903-R1 | fulgor | 8 | 14.228 / 41.052 | 116.760 / 315.530 | 2.326 / 0.367 | 2.89× |
| gencode | SRR307903-R1 | themisto | 1 | 68.486 / 366.329 | 68.430 / 383.500 | 1.637 / 3.258 | 5.35× |
| gencode | SRR307903-R1 | themisto | 8 | 14.244 / 66.528 | 116.760 / 543.340 | 2.322 / 2.520 | 4.67× |
| gencode | SRR307903-R2 | fulgor | 1 | 67.224 / 148.509 | 67.170 / 226.100 | 1.652 / 0.358 | 2.21× |
| gencode | SRR307903-R2 | fulgor | 8 | 14.024 / 41.146 | 114.720 / 315.540 | 2.338 / 0.367 | 2.93× |
| gencode | SRR307903-R2 | themisto | 1 | 67.103 / 363.457 | 67.040 / 380.890 | 1.652 / 3.248 | 5.42× |
| gencode | SRR307903-R2 | themisto | 8 | 14.031 / 66.507 | 114.760 / 538.910 | 2.326 / 2.533 | 4.74× |
| gencode | pbmc_1k_v3-R2 | fulgor | 1 | 394.847 / 697.889 | 394.580 / 1040.350 | 1.658 / 0.358 | 1.77× |
| gencode | pbmc_1k_v3-R2 | fulgor | 8 | 76.674 / 147.158 | 624.740 / 1116.260 | 2.358 / 0.369 | 1.92× |
| gencode | pbmc_1k_v3-R2 | themisto | 1 | 394.476 / 1314.084 | 394.240 / 1351.420 | 1.658 / 7.337 | 3.33× |
| gencode | pbmc_1k_v3-R2 | themisto | 8 | 76.748 / 237.640 | 625.140 / 1933.500 | 2.359 / 6.153 | 3.10× |
| nascent | SRR10919051-R1 | themisto | 1 | 94.998 / 391.924 | 94.920 / 398.530 | 9.160 / 4.654 | 4.13× |
| nascent | SRR10919051-R1 | themisto | 8 | 25.252 / 69.860 | 147.860 / 549.230 | 9.217 / 4.565 | 2.77× |
| nascent | SRR10919051-R2 | themisto | 1 | 98.119 / 399.985 | 98.040 / 406.530 | 9.160 / 4.654 | 4.08× |
| nascent | SRR10919051-R2 | themisto | 8 | 25.821 / 70.753 | 152.340 / 557.520 | 9.217 / 4.582 | 2.74× |

kallistors completed every available mapping condition sooner in this experiment. Fulgor had the smallest GENCODE mapping RSS. On nascent reads kallistors used about 9.2 GiB versus Themisto about 4.6 GiB; on full PBMC it used about 2.36 GiB versus Themisto about 6.15 GiB and Fulgor about 0.37 GiB at eight requested threads. Memory ranking therefore depends on the reference and query. These observations do not establish universal superiority over published algorithms.

## Thread accounting

Requested threads are CLI budgets, not a strict physical-core cap. Fulgor `-t 1` runs a mapping worker plus a parser; `-t 8` uses seven mapping workers plus a parser. Themisto runs the requested number of mapping callbacks plus a producer and main thread. Observed OS maxima were Fulgor 3/9 and Themisto 3/10 at requested 1/8. kallistors observed 1/10 on GENCODE and 2/10 on nascent (auxiliary startup/loading work can add an OS thread). The Fulgor one-thread CPU totals exceed wall time, confirming parser overlap. Full observations are retained separately from worker counts.

## Reference and algorithm conventions

GENCODE v50 ALL retains 670,670 transcript targets; the official Pachter mature-plus-nascent reference retains 267,211 targets (227,665 mature and 39,546 nascent). IDs, input order, sequence duplicates and short targets are preserved. Every Fulgor file color is remapped through its actual `print-filenames` order to one transcript; Themisto uses `--sequence-colors` in input sequence order. Both exact target universes are pinned in the reference manifests.

The existing prepared reference already uses uppercase ACGT, including the original deterministic replacement of ambiguity bases with A/C/G/T by position modulo four. That upstream preparation is shared. The derived competitor reference additionally matches the native index builder: if the last ten bases are A, remove the entire trailing A run. GENCODE: 11,659 affected targets, 180,582 bases removed, 65 targets shorter than k=31 retained. Nascent: 2,319 affected targets, 36,386 bases removed, 62 targets shorter than k retained. The derived FASTA bytes differ from the source; source and derived hashes are reported separately in `references/*/reference.json`. This is not a claim to benchmark unnormalized genomic sequence.

All commands are unstranded and use k=31. Fulgor stable default hybrid `.fur` uses full intersection by omitting `-r`; `-r 1` would select a different threshold-union path in this release. Themisto uses `--threshold 1`, selecting its specialized intersection path and excluding unknown k-mers. Known positive k-mer colors are combined across both orientations. The independent fixture covers duplicates, short targets/queries, reverse complements, unknown bases, unmapped reads, polyA handling and IDs. Neither tool is run with per-query deduplication or a cherry-picked index variant.

Pinned primary implementations are `source/fulgor/src/ps_full_intersection.cpp`, its query parser in `external/FQFeeder`, and `source/themisto/include/pseudoalign.hh` (intersection worker and dispatcher). They intersect positive k-mer colors; streaming lookup/color reuse can accelerate this without changing the declared constraint set. kallistors follows its frozen native graph-skipping/verification/rescue path. Measured differences below are real; this experiment has not localized the mechanism of each differing read and does not ascribe every difference to graph skipping.

## Complete mapping agreement

Weighted EC means the complete transcript target tuple and its integer read count. The L1 columns sum absolute count differences over the complete EC universe, optionally including the empty class; they are not counts of differing individual reads. kallistors CLI exports aggregate ECs only, so no kallistors per-read precision/recall/Jaccard or read disagreement count is claimed.

| Reference | Read file | Comparison | Mapped kallistors | Mapped comparator | Differing weighted ECs | Mapped EC count L1 | L1 incl. unmapped |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | SRR13638690-R1 | kallistors_vs_fulgor | 4276153 | 4256145 | 3230 | 21414 | 41422 |
| gencode | SRR13638690-R1 | kallistors_vs_themisto | 4276153 | 4256145 | 3230 | 21414 | 41422 |
| gencode | SRR13638690-R2 | kallistors_vs_fulgor | 4289613 | 4265151 | 3450 | 25070 | 49532 |
| gencode | SRR13638690-R2 | kallistors_vs_themisto | 4289613 | 4265151 | 3450 | 25070 | 49532 |
| gencode | SRR307903-R1 | kallistors_vs_fulgor | 29331097 | 29314980 | 6396 | 16675 | 32792 |
| gencode | SRR307903-R1 | kallistors_vs_themisto | 29331097 | 29314980 | 6396 | 16675 | 32792 |
| gencode | SRR307903-R2 | kallistors_vs_fulgor | 29356733 | 29341756 | 6647 | 15699 | 30676 |
| gencode | SRR307903-R2 | kallistors_vs_themisto | 29356733 | 29341756 | 6647 | 15699 | 30676 |
| gencode | pbmc_1k_v3-R2 | kallistors_vs_fulgor | 43165816 | 42936597 | 34172 | 232609 | 461828 |
| gencode | pbmc_1k_v3-R2 | kallistors_vs_themisto | 43165816 | 42936597 | 34172 | 232609 | 461828 |
| nascent | SRR10919051-R1 | kallistors_vs_themisto | 9202826 | 9122291 | 17300 | 80921 | 161456 |
| nascent | SRR10919051-R2 | kallistors_vs_themisto | 9181806 | 9091452 | 18335 | 90618 | 180972 |

The eight-thread aggregate comparisons are identical to the one-thread comparisons above. Every tool passed its complete weighted EC 1-versus-8 check. Fulgor and Themisto matched every individual read target set on all five GENCODE files at both thread settings, including unmapped reads: 139,636,526 input reads across those five files. Themisto also matched every read target set between 1 and 8 threads on both nascent mates. Full per-read Jaccard/precision/recall is 1 for the available Fulgor/Themisto comparisons, excluding both-unmapped pairs from those means. No cross-tool per-read claim is available for nascent or for kallistors.

For PBMC, kallistors mapped 43,165,816 versus 42,936,597 for Fulgor/Themisto. For nascent R1 it mapped 9,202,826 versus 9,122,291 for Themisto. More mapped reads is not evidence of biological accuracy. Complete aggregate and available per-read comparison records are in `mapping-agreement.json` and `mapping-agreement.tsv`; every full native correctness output is preserved as a lossless, SHA-verified gzip archive with canonical counters and read vectors.

## Setup at eight requested threads

These are single successful setup observations, not median-of-three index comparisons. Times exclude tool compilation and reference-adapter preparation. Disk columns count serialized indexes only, not sources, references, per-transcript shards, manifests, correctness archives or temporary peak disk use. Builder memory arguments of 8 GiB are resource hints, not RSS limits. The frozen original-built kallistors indexes were reused; no fresh kallistors setup time is mixed into this table.

| Reference | Tool | Wall s | User CPU s | System CPU s | Peak RSS GiB | Logical index GiB | Allocated index GiB | Contention |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gencode | fulgor | 50.926 | 157.09 | 40.39 | 14.170 | 0.350 | 0.350 | False |
| gencode | themisto | 202.120 | 506.08 | 5.83 | 8.945 | 0.480 | 0.495 | False |
| nascent | themisto | 1492.039 | 2438.74 | 33.29 | 12.477 | 3.233 | 3.233 | True |

### Build-time comparison with kallistors

The same campaign measured kallistors index construction on this host and the same reference
target sets at eight requested threads: median **11.123517 s** for GENCODE and **98.687386 s**
for mature + nascent, each after one warmup and three builds. Dividing the competitor setup
wall times above by those kallistors medians gives **4.58× for Fulgor** and **18.17× for
Themisto** on GENCODE, and **15.12× for Themisto** on mature + nascent. These directly compare
observed construction cost; Fulgor/Themisto contribute one successful setup observation each.
The [README index table](../README.md#index-construction) includes original kallisto too;
[all kallisto/kallistors build observations](kallisto-comparison-results.md#index-construction)
retain the three-run measurements and resource use. The mapping runs themselves reused
original-built kallistors indexes, as described above.

### Recorded setup failures and resource recovery

Fulgor nascent aborted after 73.58 OS wall seconds, 555.88 user + 17.12 system CPU seconds, peak RSS about 9.57 GiB. Its bundled GGCAT `build_unitigs.rs:176` unwrapped a missing unitig lookup. No index or nascent mapping time was produced. Root cause remains unconfirmed; there was no source patch or parameter sweep. The failure is in `nascent-setup-failure.json` and its raw build log.

The first Themisto nascent attempt exited after 32.02 OS wall seconds, 149.06 user + 6.17 system CPU seconds, peak RSS about 9.27 GiB, with a KMC LUT-prefix read error. The preserved 217,164-byte prefix declared a signature map requiring over 1 MiB, proving that file was structurally invalid. About 22 GiB of failed suffix scratch had left only 2.4 GiB free. Logs do not prove the original cause. After preserving the prefix and suffix boundary samples, failed scratch was retired and SHA-identical immutable historical evidence files were consolidated with hardlinks: every original path and byte was retained. One retry with the same executable, reference, parameters and threads succeeded after roughly 80 GiB became free. Recovery supports resource pressure as an explanation but does not establish causality.

The resource retry and Themisto-only nascent workflow were saved as new scripts and provenance amendments. The frozen GENCODE workflow and completed observations were not edited. See `nascent-themisto-resource-amendment.json`, `nascent-themisto-only-amendment.json`, and `storage-recovery/remaining-historical-evidence-sharing.json`. Shared historical aliases must be copied to a fresh inode before future mutation. Earlier storage guards and regenerable shard/build-cache retirement are documented in `storage-recovery/amendment.json` and companion records. Canonical mapping evidence, frozen executables and referenced inputs are retained.

## Reproducibility and environment

- fulgor: v4.2.0, commit `ef6541952c1153c1b7b4bfa00aadbf741a340eff`, executable SHA-256 `a058593dab54d5d5191593e54100302d63219e310f4ccb2ffbb85020b32ce352`.
- themisto: v3.2.2, commit `08685f75674634205c120333e991a1ef889b826e`, executable SHA-256 `01bb654abd0a7995aa5325f3a2386e77aa515e5e720dba8e8a548e9f578ae627`.
- kallistors executable SHA-256 `6cd5e9f4e9a9ae2e2242947d73fed623822fbc88a1e9a0afafc365f5b7f6668d`.
- Exact native commands: `commands.json`, `native-observations.json`, and per-run logs/OS resource files. Frozen protocol files: `mapping-protocol-gencode.json` and `mapping-protocol-nascent-themisto-only.json`.
- Recursive source/submodule pins: `fulgor-source.json`, `themisto-source.json`, `dependency-provenance-supplement.json`; compiler flags, runtime dependencies, configuration/build commands and executable hashes: `*-build.json` and `tool-compatibility-amendment.json`.
- Fulgor: Apple Clang 21.0.0 with native CPU flags and native Rust. Themisto: local GCC 16.2.0/OpenMP with native CPU flags, scoped Rust 1.77.2 for its pinned older Rust dependencies, and `-Wno-error=template-body` for a GCC diagnostic severity change. Local compiler/runtime bottles were relocated and re-signed as packaging; no competitor algorithm source was modified.
- Themisto display banner captured the parent repository git description because of upstream CMake directory handling. It is not a trustworthy version identifier; the actual source commit and executable hashes above are authoritative.
- Derived Fulgor transcript files are regenerable from pinned normalized FASTA and `colors.tsv`. GENCODE shards were retired during storage recovery and have an exact restoration helper. Build caches may need re-creation for rebuilding, but runtime binaries and required local GCC libraries are preserved.

## Original-kallisto comparison

See the [separate original-kallisto report](kallisto-comparison-results.md) for full quantification,
bootstrap, single-cell and index measurements. Its dedicated
[sequence-bias section](kallisto-comparison-results.md#completed-sequence-bias-observations)
retains completed timings and material native abundance differences. Bias timings are outside
this mapping-only comparison.

## Evidence status

All available panels and complete exact within-tool checks passed. Fulgor nascent remains an explicit coverage limitation. `final-verification.json` verifies the pinned evidence closure, including archived mapping output hashes, without rerunning native work. `COMPLETE.json` records terminal completion of this comparison with the documented unavailable condition; it does not claim that missing measurement was performed.
