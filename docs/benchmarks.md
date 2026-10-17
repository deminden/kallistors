# Benchmarks

Compare current kallistors with official kallisto v0.52.0 using the same index and trimmed
FASTQs. The quantification runner requires complete EC parity before timing, saves executable
and input hashes, and validates retained abundance/HDF5 outputs after timing. Results depend on
the input, reference, hardware, power mode and compiler settings.

## Tools

Python 3.9 or newer runs the benchmark tools. Native HDF5 validation also needs h5py and NumPy:

```bash
python3 -m venv target/hdf5-native-check
target/hdf5-native-check/bin/python -m pip install h5py numpy
```

Put the official timing executable at `target/tools/kallisto-0.52.0`, or pass `--kallisto`.
Build a separate instrumented original as described in the [debugging guide](kallisto_debug_build.md).
Pass it with `--diagnostic-kallisto`; it supplies EC evidence and is never the timing control.

| Command | Purpose |
| --- | --- |
| `prepare_benchmark.py` | ENA download verification and fasterp trimming |
| `prepare_reference.py` | GENCODE cDNA or Pachter Lab NAC FASTA preparation |
| `subsample_benchmark_pairs.py` | Exact seeded paired subsets |
| `benchmark_reads.py` | Quantification, thread scaling, optional stage timing and EC preflight |
| `benchmark_index.py` | Index build time, size and peak RSS |
| `report_benchmarks.py` | Portable JSON/TSV export from completed runs |
| `compare_ec_counts.py` | Independent comparison of complete counted EC maps |
| `compare_abundance.py` | Transcript abundance and effective-length comparisons |
| `check_hdf5_output.py` | Independent native HDF5 validation and fragment-model statistics |
| `diagnostics/prepare_quant_replay.py` | Fixed-EC, fixed-fragment-model inputs for EM investigation |

`benchmark_common.py` supplies shared execution, memory measurement, hashes, downloads and
result writing. The [ENA manifest](../scripts/benchmarks/human_datasets.json) defines four
public datasets and their directories. The [suite manifest](../scripts/benchmarks/quantification.json)
defines larger-library/reference combinations; the runner contains no accession-specific plans.

## Retrieve and trim reads

Use fasterp 0.2.1 with the included [bounds fix](../scripts/benchmarks/fasterp-0.2.1-trim-bounds.patch).
The released version can panic when front-quality trimming and overlap adapter clipping leave
disjoint intervals (observed on SRR1039508). The patch keeps that read empty so the existing
minimum-length filter rejects its pair. It changes no kallistors code or quantification settings.

Build the pinned source outside this Rust workspace (requires Rust, a C/C++ toolchain and CMake):

```bash
FASTERP_SOURCE=$(mktemp -d)
curl -fL https://static.crates.io/crates/fasterp/fasterp-0.2.1.crate \
  -o "$FASTERP_SOURCE/fasterp.crate"
python3 -c 'import hashlib,pathlib,sys; assert hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest() == "71ded328b10ccbff9cf55e57ec0a359e28b6fc4d71977a46dc48a73d869077e4"' \
  "$FASTERP_SOURCE/fasterp.crate"
tar -xzf "$FASTERP_SOURCE/fasterp.crate" -C "$FASTERP_SOURCE"
patch -d "$FASTERP_SOURCE/fasterp-0.2.1" -p1 \
  < scripts/benchmarks/fasterp-0.2.1-trim-bounds.patch
cargo install --locked --path "$FASTERP_SOURCE/fasterp-0.2.1" \
  --root "$PWD/target/tools/fasterp-fixed"
python3 scripts/prepare_benchmark.py --samples all \
  --fasterp target/tools/fasterp-fixed/bin/fasterp
```

`--samples all` selects the four current manifest entries: the small sample, two larger libraries and a TT-seq sample.
`--samples SRR1039508,SRR307903` selects those two. The manifest pins mate URLs and MD5s;
compressed sizes are checked when supplied. Downloads support resumption and gzip integrity is
checked. `--data DIRECTORY` changes the data root; `--fasterp PATH` selects the trimmer.

| Run | Source | Deposited pairs |
| --- | --- | ---: |
| SRR13638690 | Temporal muscle, PRJNA699558 | 4,451,399 |
| SRR1039508 | Untreated airway smooth muscle N61311, PRJNA229998 | 22,935,521 |
| SRR307903 | ENCODE BJ fibroblast poly(A)+ RNA, PRJNA30709 | 36,010,036 |
| SRR10919051 | HEK293 wild-type control 2, 1-hour 4sU TT-seq, PRJNA602480 | 10,566,118 |

Preparation enables fasterp adapter detection, front/tail trimming with a four-base Q20 window,
minimum length 50, qualified-base Q15, at most 40% unqualified bases and five Ns. Defaults are
four threads and gzip level 4. `trimmed-fasterp/preparation.json` records settings, trimmer
version/hash, input checksums, output hashes and actual before/after pair counts. Cached outputs
must match those records. The benchmark verifies prepared-read hashes before measuring.

A pair is two FASTQ records. Retrieval and trimming are outside quantification timing.

The nascent sample is [GSM4277828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4277828),
a human HEK293 control labelled with 4-thiouridine for one hour. TT-seq enriches newly synthesized
RNA and includes unspliced transcripts; it is not a pure-unspliced preparation. Its reference is
only the official Pachter Lab mature-plus-nascent FASTA below. Mature-only/nascent-only/mixed
EC partitions describe compatible target sets, not a validated biological splicing classification.

## Prepare references

Reference preparation writes FASTA provenance and target counts in `reference.json`.
Index construction is a separate step using the same official builder for both quantifiers.

### GENCODE v50 ALL cDNAs

```bash
python3 scripts/prepare_reference.py --reference gencode50
target/tools/kallisto-0.52.0 index -t 8 \
  -i data/srr-small/gencode.v50.kallisto.idx \
  data/srr-small/gencode.v50.transcripts.fa.gz
```

The ALL transcript FASTA has 670,670 targets and MD5 `4d5bd41b247480e12c0b370490f52ddd`.

For construction measurements that exclude gzip decoding, decompress the same verified file:

```bash
mkdir -p data/gencode50-all
gzip -dc data/srr-small/gencode.v50.transcripts.fa.gz \
  > data/gencode50-all/gencode.v50.transcripts.fa
```

The uncompressed ALL FASTA is 1,575,713,375 bytes, with SHA-256
`179bbe798a0e455d05b9456ab5f7ef51bffdcab6118d4697ea3219181c5e7c93`.
Pass this file to `benchmark_index.py --fasta` when reproducing the uncompressed measurements.

### Smaller index-building target: protein-coding transcripts

```bash
python3 scripts/prepare_reference.py --reference gencode50-protein-coding \
  --transcripts data/srr-small/gencode.v50.transcripts.fa.gz
python3 scripts/benchmark_index.py \
  --fasta data/gencode50-protein-coding/gencode.v50.protein_coding.fa \
  --threads 8 --warmups 1 --runs 3 --stages --keep-last-index \
  --out target/benchmarks/protein-coding-index
```

This selects 288,650 transcripts with biotype `protein_coding`, the eighth pipe-delimited
GENCODE header field, from the checksum-verified ALL FASTA. IDs, sequences and input order are
preserved. Transcripts with other biotypes are excluded even when their gene is protein
coding. Both builders receive the same uncompressed FASTA. The preparation record includes
source and output hashes and target/base counts. This optional smaller target does not replace
full-reference compatibility checks.

### Intronic reference: Pachter Lab human NAC

```bash
python3 scripts/prepare_reference.py --reference official-nac
target/tools/kallisto-0.52.0 index -t 8 -k 31 -m 23 \
  -i data/official-human-nac/human_nac.mature-plus-nascent.idx \
  data/official-human-nac/human_nac.mature-plus-nascent.fa
```

Intronic benchmarks use only Pachter Lab's published
[mature](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_nac.mature.fa.gz)
and [nascent](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_nac.nascent.fa.gz)
FASTAs, concatenated in that order with IDs and sequences preserved. They use Ensembl 108
with selected gene biotypes:
227,665 mature and 39,546 nascent targets. Asset sizes and gzip integrity are checked;
retrieved SHA-256 hashes are recorded, without claiming publisher-provided checksums.

`--out DIRECTORY` changes the reference directory. `--transcripts FASTA.gz` can reuse the
verified GENCODE cDNA file for `--reference gencode50`. Reference preparation does not invoke
platform-specific timing tools.

## Exact paired subsets

To select four million airway pairs after trimming:

```bash
python3 scripts/subsample_benchmark_pairs.py \
  --reads data/srr1039508/trimmed-fasterp/SRR1039508_1.trimmed.fastq.gz \
  --reads2 data/srr1039508/trimmed-fasterp/SRR1039508_2.trimmed.fastq.gz \
  --total-pairs "$(python3 -c 'import json; print(json.load(open("data/srr1039508/trimmed-fasterp/preparation.json"))["pairs_after"])')" \
  --pairs 4000000 --seed 42 --out data/srr1039508/trimmed-fasterp/subset
```

Use SRR307903's paths and recorded pair count for its subset. Selection is uniform without
replacement, checks mate names and FASTQ structure, preserves order, and writes deterministic
compressed mates plus `subset.json`. Subset compression is separate from fasterp compression.

## Quantification timing

Build the candidate with locked dependencies. Source builds use `target-cpu=native`; CI uses
`target-cpu=generic` for portable artifacts. Record compiler/CPU-target settings with results.

```bash
cargo build --locked --release --package kallistors --bin kallistors
python3 scripts/benchmark_reads.py \
  --samples SRR13638690 --modes paired,single --thread-counts 1,8 \
  --diagnostic-kallisto target/tools/kallisto-src/build/src/kallisto \
  --hdf-python target/hdf5-native-check/bin/python \
  --warmups 1 --runs 3 --stages --out target/benchmarks/read-scaling
```

All requested configurations must pass complete EC preflight before any timed quantification.
Each uses the same index bytes, mode, fragment parameters and thread count in both tools.
`--reads FASTQ --reads2 FASTQ --index INDEX` supplies custom paired inputs. For single-end-only
input, use `--reads FASTQ --modes single`. Single-end defaults are mean 200 and SD 20, configurable
with `--fragment-length` and `--fragment-sd`. Use `--reference-metadata reference.json` to attach
reference provenance to a direct run.

Defaults are eight threads, one full warmup and three measured runs per tool/configuration.
Measured order alternates. Every execution includes index loading, all reads, model fitting,
EM, normal TSV/HDF5 output and shutdown. Last measured outputs are retained for validation;
other outputs are removed after their records are saved. Counter differences and disagreement
with preflight EC totals fail the run. Equal counters alone never establish EC parity.

Peak RSS is enabled by default on macOS/Linux; `--no-memory` disables it. macOS uses
`/usr/bin/time -l` and also records peak footprint. Linux uses GNU time and converts KiB to
bytes. `--stages` additionally records live progress-marker arrival times on either platform.
CPU contention is sampled outside the benchmark process tree; available power/swap context
and OS elapsed/CPU times are retained. Apple Silicon hardware is identified by family only. An OS-real versus stage-clock gap over two seconds is
flagged as interrupted timing. Interrupted records remain saved and cannot be exported as
an uninterrupted result; use a fresh series for replacement measurements.

| Comparable phase | Boundary |
| --- | --- |
| Loading/setup | Process start to read-processing marker, including index/input setup |
| Read processing | Read-processing start to end; decompression, matching and EC counting overlap |
| Model + EM | Read-processing end to EM completion |
| Output/shutdown | EM completion to exit, including writing and cleanup |

Individual stage medians need not sum to the median total. TSV phase columns use the phases
of the median-total execution; for even repetition counts, they average the two central runs.
Progress-marker timing is an external approximation, not an isolated kernel measurement.

`--pause-file PATH` exits with status 75 after a complete execution if that file exists.
Remove it and pass `--resume` to continue the same output directory. Resume checks frozen
input/executable/script hashes and protocol, preserves completed records, and adds one warmup
per tool/configuration (`--resume-warmups` controls that count). Incomplete attempts are archived.
Changing inputs or tools requires a new directory. `protocol.json`, `runs.json` and `summary.json`
use the same schema across direct and suite quantification runs.

## Configured suites and reports

Prepare the full libraries and references first, then select a manifest group:

```bash
python3 scripts/benchmark_reads.py \
  --suite scripts/benchmarks/quantification.json --group gencode50 \
  --diagnostic-kallisto target/tools/kallisto-src/build/src/kallisto \
  --hdf-python target/hdf5-native-check/bin/python \
  --thread-counts 8 --runs 3 --warmups 1 --stages \
  --out target/benchmarks/gencode50-suite
python3 scripts/report_benchmarks.py \
  --benchmark target/benchmarks/gencode50-suite \
  --json target/benchmarks/gencode50-results.json \
  --tsv target/benchmarks/gencode50-results.tsv
```

`--group gencode50` selects cDNA-only cases; `--group official-nac` selects the nascent-RNA
case using the Pachter Lab mature-plus-nascent reference. A suite case names its accession,
index, reference metadata, modes and optional subset; read paths come from the ENA manifest.
Custom suites use the same schema. Reporting rechecks EC evidence, validated output hashes, run completeness, warmups,
counters and interruption status, then exports portable JSON and TSV without rerunning quant.
The reports include HDF5 validation, abundance differences, fragment statistics and reference
provenance. NAC ECs are partitioned as `mature_only`, `nascent_only` or `mixed` using the
recorded mature-target boundary.

For index build comparisons:

```bash
python3 scripts/benchmark_index.py \
  --fasta data/srr-small/gencode.v50.transcripts.fa.gz --threads 8 \
  --runs 3 --warmups 1 --stages --keep-last-index \
  --out target/benchmarks/index-build
```

Index benchmarking shares execution, provenance, memory and summary code with quantification.
Its outputs can be exported with `report_benchmarks.py` using the same command structure.
`--keep-last-index` retains the last measured index from each builder for compatibility checks;
other generated indexes are removed after measurement. Without that option all generated indexes
are removed. `--stages` records kallistors' internal build stages in JSON and TSV:

| Build stage | Work included |
| --- | --- |
| `fasta_parse` | FASTA reading and transcript normalization |
| `graph_build` | K-mer collection/sorting, topology, compaction and k-mer locations |
| `ec_build` | Transcript mapping and EC/position blocks |
| `minimizer_mphf` | Selected minimizer indexing and MPHF construction |
| `write` | Serialization and file writing, excluding MPHF construction |

These stages do not overlap. Unattributed setup and cleanup can leave their sum below the
internal total; independently calculated stage medians need not sum to the median total.
Use external whole-process elapsed time for speed comparisons. Original kallisto's normal
progress log is retained, without claiming equivalent internal stage measurements.

Index time/size comparisons do not by themselves establish index compatibility. Load both
builders' indexes in both quantifiers and compare complete counted EC maps on identical reads,
including paired/single and relevant strand options. Index bytes may differ between builders.

### Eight-thread quantification measurement

Current release measurements cover four complete human libraries, including two larger bulk
libraries and a genuine nascent-enriched TT-seq library. All use fasterp 0.2.1, identical trimmed
mates/index bytes per comparison, eight requested threads, one warmup and three measured
complete runs. All current measurements use the fixed binary. The larger cDNA comparisons retain their
unchanged control observations, with fresh complete original checks interleaved with current
runs. TT-seq rotates all three programs. The small-library controls were refreshed after a
single-thread check showed timing drift; the single-thread comparison uses current/original
only. CPU activity, power, swap and interruption observations are retained. Small original
checks also serve as the first observations in the refreshed original medians.

The first three libraries use GENCODE v50 ALL (670,670 targets). TT-seq uses
only the official Pachter Lab mature-plus-nascent reference (267,211 targets; 6.58 GB index).
All indexes use k=31/g=23. Retrieval/trimming are outside timing.

| Library | Trimmed pairs | Current s | v0.5.0 s | Original s | Current vs original |
| --- | ---: | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle | 4,410,990 | 8.065 | 17.341 | 45.158 | 5.60× |
| SRR1039508 · airway | 22,297,617 | 53.288 | 44.770 | 178.942 | 3.36× |
| SRR307903 · fibroblast | 32,878,740 | 71.653 | 85.372 | 212.451 | 2.97× |
| SRR10919051 · TT-seq | 10,261,737 | 47.921 | 91.817 | 98.344 | 2.05× |

Peak RSS is the whole-process high-water mark; GiB means 2^30 bytes. macOS peak footprint is
recorded separately in the [current measurement evidence](benchmark-results.json).

| Library | Current GiB | v0.5.0 GiB | Original GiB | RSS reduction vs original |
| --- | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle | 1.218 | 3.064 | 2.847 | 57.2% |
| SRR1039508 · airway | 1.338 | 3.801 | 3.129 | 57.2% |
| SRR307903 · fibroblast | 1.326 | 3.702 | 3.091 | 57.1% |
| SRR10919051 · TT-seq | 9.174 | 14.960 | 14.973 | 38.7% |

#### Where time is spent

| Library / tool | Startup/setup s | Read processing s | Model + EM s | Output/shutdown s |
| --- | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle / original | 5.074 | 9.561 | 29.436 | 1.041 |
| SRR13638690 · muscle / current | 1.255 | 4.487 | 1.769 | 0.529 |
| SRR1039508 · airway / original | 5.172 | 31.025 | 141.934 | 1.107 |
| SRR1039508 · airway / current | 1.256 | 16.456 | 35.029 | 0.544 |
| SRR307903 · fibroblast / original | 5.201 | 60.765 | 145.147 | 1.106 |
| SRR307903 · fibroblast / current | 1.251 | 33.578 | 36.240 | 0.542 |
| SRR10919051 · TT-seq / original | 45.676 | 36.970 | 13.561 | 2.403 |
| SRR10919051 · TT-seq / current | 15.671 | 27.202 | 4.981 | 0.103 |

On the larger cDNA libraries, original model/EM dominates total time; index setup remains
about five seconds. Read processing also improves, so the observed gains extend beyond loading.
For TT-seq against the larger NAC reference, loading/setup is
45.68 → 15.67 s,
read processing 36.97 → 27.20 s, and model/EM
13.56 → 4.98 s (original → current). Loading contributes
substantially to this case's gain; the read and model phases also improve. These are phase
observations, not isolated measurements of individual optimizations.

The read phase includes gzip input, parsing, matching and EC counting. Startup includes index
loading and input setup; output includes deferred names, complete default TSV/HDF5 and shutdown.
Model/EM timings use each program's default fragment sample. Progress-marker arrival is an
external approximation; independent phase medians need not sum to the median total. The
released v0.5.0 binary lacks phase boundaries, so its stages are not inferred. The recorded
matched-native v0.5.0 control on the small library also supports a source-code comparison.
These phase tables distinguish read work from index-loading and EM overhead; no universal
speed factor or precise attribution to an individual optimization is claimed.

#### Complete EC and output validation

All current/original full-library EC maps are exact, independent of EC numbering:

| Library | Nonzero ECs | Aligned pairs | Differing ECs | Count L1 |
| --- | ---: | ---: | ---: | ---: |
| SRR13638690 | 190,693 | 4,246,399 | 0 | 0 |
| SRR1039508 | 742,488 | 21,137,312 | 0 | 0 |
| SRR307903 | 744,744 | 29,637,796 | 0 | 0 |
| SRR10919051 | 212,719 | 8,890,219 | 0 | 0 |

The published v0.5.0 binary does **not** preserve the same complete EC map on every larger case
(SRR1039508: 10 differing ECs, count L1 21; SRR307903: 18 differing ECs, count L1 31; SRR10919051: 5
differing ECs, count L1 5). Its rows are delivered-binary time/RSS observations, not equivalent-work
speedup claims. Current/original parity is mandatory and passes. Fixes in 0.5.1 include explicit
abundant-node associations and paired fallback/strand behavior.

Observed total-time regression: SRR1039508 · airway: current takes 19.0% longer than the published
v0.5.0 download, with 64.8% lower peak RSS. Current model/EM takes 35.03 s. EC/model differences prevent assigning the total regression
to one code change from these timings alone.

For TT-seq, the shared counted map partitions by Pachter's mature-target boundary as follows:
mature only: **251,972 pairs (2.83%)**; nascent only: **6,676,534 pairs (75.10%)**; mixed:
**1,961,713 pairs (22.07%)**. These are target-compatibility partitions, not a measurement of
biological splicing accuracy. A mixed EC contains both target types. The assay enriches newly synthesized
RNA and is not exclusively unspliced.

Exact EC parity does not imply equal abundances. Concurrent scheduling can change the
near-10,000-observation fragment sample, effective lengths and EM stopping. v0.5.0 instead
collects all qualifying observations. The retained eight-thread current/original outputs show:

| Library | Fragment observations, current / original | Estimated-count L1 / original estimate mass | Maximum count difference |
| --- | ---: | ---: | ---: |
| SRR13638690 | 11,828 / 11,828 | 0.000000% | 8.73115e-10 |
| SRR1039508 | 17,383 / 17,383 | 0.000000% | 1.74623e-10 |
| SRR307903 | 16,111 / 16,111 | 0.000000% | 1.07684e-09 |
| SRR10919051 | 59,241 / 59,241 | 0.000000% | 1.16415e-10 |

Every current/original retained output passes native HDF5 readback across all transcript rows. The
v0.5.0 download fails native large-vector readback on SRR1039508, SRR307903, SRR10919051
(chunk-dimension encoding); its valid 17-digit TSV supplies abundance comparisons. The small
full-reference comparison independently found the same v0.5.0 HDF5 encoding limitation; current
bounded chunks pass.

#### Single-thread control and build scope

On the 4,410,990-pair library, the current median is **38.798 s**, versus
**93.967 s** for freshly measured original controls (**2.42× faster**).
Read processing is 25.313 versus 49.658 s;
model/EM is 8.999 versus 28.598 s.
Fragment histograms are identical and effective lengths are
identical; maximum estimated-count difference
is 2.62e-09. One-thread runs cover only this library. Index construction is measured
at eight threads only.

Current source and the small matched-native v0.5.0 control use Rust 1.99.0 with
`target-cpu=native`. The control differs from its release tag only by phase messages. The
checksum-verified [published v0.5.0 download](https://github.com/deminden/kallistors/releases/tag/v0.5.0)
was built with native settings on its CI runner; its effective CPU/compiler settings are not
assumed identical to this host. Original timing uses the unmodified official kallisto v0.52.0
binary; the separately instrumented original supplies EC evidence only. Current release CI
uses `target-cpu=generic`, so source-build factors are not promises for prebuilt binaries.

All hardware descriptions are Apple Silicon. Power, interruption and sampled background CPU
flags are recorded per execution; all observations remain in the evidence. Three-run medians
across these libraries characterize these inputs/builds, not every tissue, assay or machine.
Original control checks (small-library quantification checks also form the first refreshed original observations):

| Fresh original control | Seconds | Earlier median s | Difference |
| --- | ---: | ---: | ---: |
| SRR1039508-cdna-full--paired--t8 / 8 threads | 177.604 | 178.942 | -0.75% |
| SRR307903-cdna-full--paired--t8 / 8 threads | 210.692 | 212.451 | -0.83% |
| quantification / 8 threads | 45.158 | 45.632 | -1.04% |
| quantification / 1 thread | 93.967 | 99.090 | -5.17% |
| index / 8 threads | 121.392 | 123.813 | -1.96% |

See [fragment sampling](fragment_length_sampling.md) and [compatibility](compatibility.md).

### Eight-thread index build measurement

Full uncompressed GENCODE v50 ALL: **670,670 transcripts and 207,228,052 distinct k-mers**,
k=31/g=23, eight requested threads, one warmup and three complete measured builds per binary.
The shared verified FASTA is 1,575,713,375 bytes. The same build/source controls described above
apply; all current measurements use the corrected release source.

| Builder | Median seconds | Peak RSS GiB | Peak footprint GiB |
| --- | ---: | ---: | ---: |
| Current kallistors | 11.115 | 5.765 | 4.576 |
| Published v0.5.0 download | 111.160 | 12.921 | 14.207 |
| v0.5.0, matched native build | 110.400 | 13.156 | 14.207 |
| Original kallisto v0.52.0 | 123.813 | 7.845 | 7.152 |

Current is **10.00× faster with 55.4% lower RSS** than the download,
**9.93× with 56.2% lower RSS** than matched-native v0.5.0, and
**11.14× with 26.5% lower RSS** than original. Decompression is excluded;
writing, allocation cleanup and shutdown are included.

| Current internal stage | Median seconds |
| --- | ---: |
| FASTA parsing | 0.458 |
| Graph construction | 8.372 |
| EC construction | 0.868 |
| Minimizers and MPHF | 0.869 |
| Serialization and writing | 0.465 |

Current stages are exclusive; independent stage medians need not sum to the median total.
v0.5.0's reported write timer overlaps its MPHF timer; those raw phases are not added.
Original's normal progress log has no equivalent complete construction-stage breakdown.
RSS and macOS footprint are separate whole-process high-water measurements.

Every final index reproduces the complete upstream EC map on all 4,410,990 pairs:
190,693 nonzero classes and 4,246,399 aligned pairs, with zero differences. Original kallisto
also loads the current-built index with exact parity. Current, published v0.5.0 and
matched-native v0.5.0 produce byte-identical Rust indexes. Hashes and sizes are recorded.

All observations and control rechecks are retained with CPU/power context in the
[measurement evidence](benchmark-results.json). These factors describe this reference and build.
See [construction details](development.md#index-construction) for the current implementation.

## Correctness and interpretation

Compare complete EC count maps with `compare_ec_counts.py`, independent of EC numbering.
Match index bytes and transcript ordering. Missing, malformed or empty comparison inputs fail
explicitly. The native HDF5 checker verifies IDs, lengths, counts, bootstrap metadata/values and
fragment histograms. Candidate TSV values must match HDF5 exactly; official TSV values are
checked at its six-significant-digit precision. `compare_abundance.py` reports numeric
differences and rejects missing/empty inputs, duplicate IDs and mismatched transcript sets.

[Fragment sampling](fragment_length_sampling.md) can change effective lengths, EM convergence
and ambiguous isoform estimates even with identical ECs. Compare saved histograms when
interpreting model/EM speed. Fixed EC/model replays and matched iteration settings isolate EM
behavior. Bulk poly(A) libraries can assess reference cost and compatibility; biological
nascent-RNA accuracy requires suitable assays and independent evidence.
