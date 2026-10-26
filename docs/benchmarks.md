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
| `prepare_benchmark.py` | Verified ENA/10x retrieval and fastp trimming |
| `prepare_reference.py` | GENCODE cDNA or Pachter Lab NAC FASTA preparation |
| `subsample_benchmark_pairs.py` | Exact seeded paired subsets |
| `benchmark_reads.py` | Quantification, thread scaling, optional stage timing and EC preflight |
| `benchmark_bus.py` | Complete single-cell BUS timing, counted-record parity and corrected gene/UMI matrices |
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

Use [fastp 1.4.0](https://github.com/OpenGene/fastp/releases/tag/v1.4.0), installed using
its official build instructions or a package manager. All new preparation uses fastp; pass its
executable explicitly when it is not on the search path:

```bash
python3 scripts/prepare_benchmark.py --samples all --fastp /path/to/fastp
```

`--samples all` selects the four current manifest entries: the small sample, two larger libraries and a TT-seq sample.
`--samples SRR1039508,SRR307903` selects those two. The manifest pins mate URLs and MD5s;
compressed sizes are checked when supplied. Downloads support resumption and gzip integrity is
checked. `--data DIRECTORY` changes the data root; `--fastp PATH` selects the trimmer.

| Run | Source | Deposited pairs |
| --- | --- | ---: |
| SRR13638690 | Temporal muscle, PRJNA699558 | 4,451,399 |
| SRR1039508 | Untreated airway smooth muscle N61311, PRJNA229998 | 22,935,521 |
| SRR307903 | ENCODE BJ fibroblast poly(A)+ RNA, PRJNA30709 | 36,010,036 |
| SRR10919051 | HEK293 wild-type control 2, 1-hour 4sU TT-seq, PRJNA602480 | 10,566,118 |

Preparation enables fastp paired-end adapter detection, front/tail trimming with a four-base Q20 window,
minimum length 50, qualified-base Q15, at most 40% unqualified bases and five Ns. Defaults are
four threads and gzip level 4. `trimmed-fastp/preparation.json` records settings, trimmer
version/hash, input checksums, output hashes and actual before/after pair counts. Cached outputs
must match those records. Earlier measured bulk-RNA results retain their recorded fasterp
0.2.1 provenance; they are not relabelled as fastp measurements. The benchmark verifies prepared-read hashes before measuring.

A pair is two FASTQ records. Retrieval and trimming are outside quantification timing.

### Single-cell: complete 10x PBMC 1k v3

The [10x PBMC 1k v3 dataset](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)
contains two sequencing lanes from a healthy donor, with 1,222 vendor-reported cells.
Read 1 contains a 16-base barcode and 12-base UMI; read 2 contains the transcript sequence.
Preparation downloads the complete 5,549,312,000-byte FASTQ archive, verifies SHA-256
`873f6511a7cb732c0a709d2fb9201ca4f1230a70f3a04e46c0198a9e4be896fd`, and safely extracts
the four R1/R2 files. Sample-index I1 reads are not BUS inputs.

```bash
python3 scripts/prepare_benchmark.py --tenx-pbmc-v3 --threads 8 --fastp /path/to/fastp
```

Fastp 1.4.0 trims and filters only the transcript read, using the same Q20 four-base
front/tail windows and Q15/40%/five-N/50-base filters described above, with poly-G trimming
enabled and duplication evaluation disabled. Adapter detection is single-end because the
barcode read is not an opposing transcript mate. Every passing transcript ID is matched to
its original R1 record; barcode and UMI sequences and qualities remain untouched. Full
record-order/count checks and gzip CRC validation cover both lanes, including all concatenated
gzip members. `preparation.json` records both raw and retained pair counts, commands and hashes.

The preparation also retrieves the standard v3 barcode on-list from the
[official BUStools getting-started release](https://github.com/BUStools/getting_started/releases/tag/species_mixing),
with SHA-256 `843a6f7038db8cb3c06f3dc21cc69d04139ffa10689780518d7a9e42dc2e819b`.
Install [bustools 0.45.1](https://github.com/BUStools/bustools/releases/tag/v0.45.1) and NumPy
for the shared independent validation:

```bash
python3 -m pip install numpy
python3 scripts/benchmark_bus.py \
  --bustools /path/to/bustools \
  --index data/srr-small/gencode.v50.kallisto.idx \
  --reads data/tenx-pbmc-1k-v3/trimmed-fastp/pbmc_1k_v3_S1_L00{1,2}_R{1,2}_001.fastq.gz \
  --preparation data/tenx-pbmc-1k-v3/preparation.json \
  --onlist data/tenx-pbmc-1k-v3/10xv3_whitelist.txt \
  --reference-fasta data/gencode50-all/gencode.v50.transcripts.fa \
  --threads 8 --warmups 1 --runs 3 --out target/benchmarks/pbmc-v3
```

The runner first checks 100,000 retained pairs, then measures the complete library with
one warmup and three timed runs per program, alternating program order between repetitions.
Report medians of whole-process time and peak RSS. Both programs use `bus -x 10xv3` and the
same full GENCODE v50 index. Retrieval, trimming and bustools are outside these timings.
Ordinary single-cDNA layouts with one barcode, UMI and transcript slice use compact counting
indexes, packed sequence batches and ordered parallel mapping. At multiple requested threads,
each input has a bounded prefetch queue and its own decoder; one-thread runs decode synchronously.
Mapping overlaps ordered BUS writing through a bounded queue. A bounded shared cDNA cache
uses two-way replacement, compares complete query bytes and reuses exact mapped or unmapped
answers across workers and lanes.
Counting borrows unchanged graph ECs, materializes only changed intersections and omits unused
placement metadata. Count-only scanning stops after irreversible EC rejection. Transcript names stream from the same index file; long-read statistics are
allocated only for long-read runs. Read limits also bound decoding, including interleaved records.
Tagged, paired,
BAM, explicit-strand, long-read and diagnostic modes retain their full processing paths.
Measure this workflow separately from `quant`.

Every warmup and measured BUS output is compared as a complete weighted multiset of
`(barcode, UMI, exact transcript set, flags)` records. Canonical EC IDs use full transcript
tuples; bustools sorts and aggregates records before complete byte comparison. Record order,
EC numbering, header text and unused padding are not biological differences. Processed counts
must also match the complete retained library. The first measured output from each program
then undergoes identical bustools barcode correction, sorting, UMI collapsing and gene counting;
the sparse matrix and its gene/barcode labels must match exactly. This checks compatibility
with the same downstream pipeline; it does not establish agreement with Cell Ranger or perform
cell calling. Matrix barcode counts include background barcodes and are not detected-cell counts.

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
  --threads 8 --warmups 0 --runs 2 --stages --keep-last-index \
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
  --reads data/srr1039508/trimmed-fastp/SRR1039508_1.trimmed.fastq.gz \
  --reads2 data/srr1039508/trimmed-fastp/SRR1039508_2.trimmed.fastq.gz \
  --total-pairs "$(python3 -c 'import json; print(json.load(open("data/srr1039508/trimmed-fastp/preparation.json"))["pairs_after"])')" \
  --pairs 4000000 --seed 42 --out data/srr1039508/trimmed-fastp/subset
```

Use SRR307903's paths and recorded pair count for its subset. Selection is uniform without
replacement, checks mate names and FASTQ structure, preserves order, and writes deterministic
compressed mates plus `subset.json`. Subset compression is separate from fastp compression.

## Quantification timing

Build the candidate with locked dependencies. Source builds use `target-cpu=native`; CI uses
`target-cpu=generic` for portable artifacts. Record compiler/CPU-target settings with results.
Build the timing CLI as its own target, then freeze its SHA-256 before preflight and timing.
Build replay examples in a separate invocation; rebuild the CLI afterward because development
dependencies can change transitive feature resolution and the executable hash.

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
  --runs 2 --warmups 0 --stages --keep-last-index \
  --out target/benchmarks/index-build
```

Use two measured builds per tool and no warmup for the current index protocol; report their
arithmetic mean. With two observations, the shared summary statistic equals this mean.
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
complete runs.

All current quantification runs use the normal CLI build: one warmup and three complete
measurements per configuration. The unchanged pre-adaptive 0.5.1 baseline observations are
retained, with fresh complete baseline checks interleaved after the first current measurement
in each configuration. Original kallisto and published v0.5.0 controls retain earlier verified
observations and original control rechecks. Their ratios are contextual comparisons across
sessions. All reported measurements and control checks are retained.

The first three libraries use GENCODE v50 ALL (670,670 targets). TT-seq uses
only the official Pachter Lab mature-plus-nascent reference (267,211 targets; 6.58 GB index).
All indexes use k=31/g=23. Retrieval/trimming are outside timing.

| Library | Trimmed pairs | Current s | v0.5.0 s | Original s | Current vs original |
| --- | ---: | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle | 4,410,990 | 8.119 | 17.341 | 45.158 | 5.56× |
| SRR1039508 · airway | 22,297,617 | 27.183 | 44.770 | 178.942 | 6.58× |
| SRR307903 · fibroblast | 32,878,740 | 44.249 | 85.372 | 212.451 | 4.80× |
| SRR10919051 · TT-seq | 10,261,737 | 44.307 | 91.817 | 98.344 | 2.22× |

Peak RSS is the whole-process high-water mark; GiB means 2^30 bytes. macOS peak footprint is
recorded separately in the [current measurement evidence](benchmark-results.json).

| Library | Current GiB | v0.5.0 GiB | Original GiB | RSS reduction vs original |
| --- | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle | 1.229 | 3.064 | 2.847 | 56.8% |
| SRR1039508 · airway | 1.336 | 3.801 | 3.129 | 57.3% |
| SRR307903 · fibroblast | 1.321 | 3.702 | 3.091 | 57.3% |
| SRR10919051 · TT-seq | 9.174 | 14.960 | 14.973 | 38.7% |

#### Where time is spent

| Library / tool | Startup/setup s | Read processing s | Model + EM s | Output/shutdown s |
| --- | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle / original | 5.074 | 9.561 | 29.436 | 1.041 |
| SRR13638690 · muscle / current | 1.265 | 4.531 | 1.793 | 0.534 |
| SRR1039508 · airway / original | 5.172 | 31.025 | 141.934 | 1.107 |
| SRR1039508 · airway / current | 1.271 | 16.523 | 8.828 | 0.548 |
| SRR307903 · fibroblast / original | 5.201 | 60.765 | 145.147 | 1.106 |
| SRR307903 · fibroblast / current | 1.262 | 33.640 | 8.741 | 0.554 |
| SRR10919051 · TT-seq / original | 45.676 | 36.970 | 13.561 | 2.403 |
| SRR10919051 · TT-seq / current | 15.713 | 27.260 | 1.181 | 0.112 |

On the larger cDNA libraries, original model/EM dominates total time; index setup remains
about five seconds. Read processing also improves, so the observed gains extend beyond loading.
For TT-seq against the larger NAC reference, loading/setup is
45.68 → 15.71 s,
read processing 36.97 → 27.26 s, and model/EM
13.56 → 1.18 s (original → current). Loading contributes
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

#### Adaptive EM scheduling

Scheduling uses the counted EC graph and requested thread count. Components exceeding both
262,144 edges and one worker's fair share of all multi-EC edges share their ordered denominator
and transcript work across the same pool. Smaller components retain local rounds and cache
reuse. Each scalar sum stays in its original order; only integer changes and fixed-point flags
are reduced. Eligible single-thread and bias-corrected runs share the compact component engine;
bias weight refreshes retain their original boundaries and reset fixed-work certificates.
Zero-iteration and exceptional model cases keep their general fallbacks. Graph topology and
requested workers select component scheduling, independently of tissue or accession.

The largest component accounts for 6.0% of multi-EC edges in the small library, 84.8% in airway,
85.4% in fibroblast and 91.6% in TT-seq. This topology explains why scheduling across components
alone cannot use eight workers effectively on all these inputs.

| Library | Component-local s | Adaptive s | Total speedup | Component-local / adaptive RSS GiB |
| --- | ---: | ---: | ---: | ---: |
| SRR13638690 · muscle | 8.079 | 8.119 | 1.00× | 1.223 / 1.229 |
| SRR1039508 · airway | 53.160 | 27.183 | 1.96× | 1.337 / 1.336 |
| SRR307903 · fibroblast | 71.464 | 44.249 | 1.62× | 1.320 / 1.321 |
| SRR10919051 · TT-seq | 48.079 | 44.307 | 1.09× | 9.175 / 9.174 |

These complete-process medians use the verified standard CLI builds and control protocol above. The small
case is interpreted as unchanged within timing noise. Fixed-EC/fixed-histogram full-convergence
replays additionally produce byte-identical output for all four libraries, independent of fragment
sampling variation. Default-run fragment histograms differ on SRR13638690 · muscle, SRR307903 ·
fibroblast, SRR10919051 · TT-seq, as concurrent sampling permits; complete abundance audits are
retained. Replay times are diagnostic single executions, not repeated speed estimates.
Executable/input hashes, all observations and output audits are in the measurement evidence.
Fresh baseline checks against the retained three-run medians:

| Configuration | Fresh baseline s | Retained median s | Difference |
| --- | ---: | ---: | ---: |
| SRR13638690 | 8.118 | 8.079 | +0.48% |
| SRR10919051-official-nac-full--paired--t8 | 48.007 | 48.079 | -0.15% |
| SRR1039508-cdna-full--paired--t8 | 53.095 | 53.160 | -0.12% |
| SRR307903-cdna-full--paired--t8 | 71.954 | 71.464 | +0.69% |
| SRR13638690--t1 | 39.042 | 38.872 | +0.44% |


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

Adaptive EM shares the dominant airway component across threads, removing the observed
regression against published v0.5.0. Complete fixed-EC/fixed-histogram replays preserve every
output byte, including full convergence and the final threshold update.

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
| SRR307903 | 16,205 / 16,111 | 0.000731% | 4.27427 |
| SRR10919051 | 59,355 / 59,241 | 0.000160% | 2.03495 |

Every current/original retained output passes native HDF5 readback across all transcript rows. The
v0.5.0 download fails native large-vector readback on SRR1039508, SRR307903, SRR10919051
(chunk-dimension encoding); its valid 17-digit TSV supplies abundance comparisons. The small
full-reference comparison independently found the same v0.5.0 HDF5 encoding limitation; current
bounded chunks pass.

#### Single-thread control and build scope

On the 4,410,990-pair library, the current median is **38.855 s**, versus
**93.967 s** for original controls (**2.42× faster**).
The retained pre-adaptive baseline is 38.872 s.
Observed ranges are 38.736–39.393 s baseline and
38.652–38.964 s current; interpret differences within these ranges as timing variation.
Read processing is 25.292 versus 49.658 s;
model/EM is 8.963 versus 28.598 s.
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
k=31/g=23, eight requested threads. The verified FASTA is 1,575,713,375 bytes.
Current construction uses **two complete measured builds, arithmetic mean, no warmup**.
Fresh preceding-builder controls bracket the candidates in control/current/current/control order.
Original and published v0.5.0 retain their verified three-run medians; ratios to them compare
sessions and are contextual. Retrieval and decompression are excluded; writing, cleanup and
shutdown are included.

| Builder | Seconds | Peak RSS GiB | Peak footprint GiB |
| --- | ---: | ---: | ---: |
| Current kallistors, two-run mean | 11.449 | 5.323 | 4.558 |
| Published v0.5.0 download | 111.160 | 12.921 | 14.207 |
| v0.5.0, matched native build | 110.400 | 13.156 | 14.207 |
| Original kallisto v0.52.0 | 123.813 | 7.845 | 7.152 |

Current observations are **11.563 / 11.336 s**.
The fresh preceding-builder mean is **11.280 s**;
the current mean is 1.5% higher, within the 2% guard, with all four indexes byte-identical. Large-reference scratch
layouts are not selected for this smaller graph.

| Current internal stage | Mean seconds |
| --- | ---: |
| FASTA parsing | 0.428 |
| Graph construction | 8.621 |
| EC construction | 0.894 |
| Minimizers and MPHF | 0.919 |
| Serialization and writing | 0.507 |

Current stages are exclusive. Whole-process timing also includes startup and shutdown;
v0.5.0's write timer overlaps MPHF and must not be added to its MPHF timer.
Peak RSS and macOS footprint are separate whole-process high-water measurements, averaged
across the same two runs.

All current GENCODE outputs have the previously cross-loaded index SHA-256
`ca6427c35fe6ea971122caf69cd7e364b6a9fc90fb7a81af6429d47e4c5e5738`.
The complete EC proof therefore applies to these same bytes: **4,410,990 pairs,
190,693 nonzero ECs and 4,246,399 aligned pairs**, with zero differences from original.

### Mature-plus-nascent index construction

Apple Silicon with 32 GiB RAM, eight requested threads, k=31/g=23. The uncompressed official
Pachter Lab mature-plus-nascent FASTA (Ensembl 108) contains **267,211 targets and
2,394,099,698 bases**. Its graph has **1,563,276,809 k-mers and 21,707,203 unitigs**.
Use the two official FASTA files unchanged; nascent sequences are not reconstructed.

Current uses **two complete measured builds, arithmetic mean, no warmup**:
**96.955 / 95.850 s**.
Fresh preceding-builder controls bracket both candidates. A qualifying change requires at least
2% faster mean total construction, no mean peak RSS increase, exact complete index hashes,
and the full GENCODE eight-thread guard. Individual observations and contention are retained.

| Builder | Total seconds | Peak RSS GiB | Peak footprint GiB | Index GB |
| --- | ---: | ---: | ---: | ---: |
| Current kallistors, two-run mean | 96.403 | 18.397 | 18.087 | 6.382 |
| Original kallisto v0.52.0, retained single control | 903.511 | 19.729 | 33.390 | 6.583 |

The contextual ratio to the retained original control is **9.37× faster**, with
**6.8% lower peak RSS**. Original was measured in an earlier session; this is a
preliminary comparison, not a fresh matched original benchmark. Memory pressure, filesystem
cache and machine capacity can affect these measurements. No warmup or cache flushing was used
for the new builds. System-wide swap counters do not measure process-specific swapping.

| Current internal stage | Mean seconds |
| --- | ---: |
| FASTA parsing | 0.370 |
| Graph construction | 75.945 |
| EC construction | 7.012 |
| Minimizers and MPHF | 7.675 |
| Serialization and writing | 3.920 |

Construction combines bounded collections, exact four-byte adaptive edges, early five-byte
key packing, a temporary prefix-owned rejection filter, and exact sparse full-word walk
landmarks. Full-key verification, ordered claims and complete fallback layouts preserve graph
membership. Compact queries and earlier bitmap/MPHF scratch retirement reduce overlapping
allocations. See [construction details](development.md#index-construction).

Every complete current index has SHA-256
`2a5157b5ea55cf8fbb9c5ac4e7dc9d1e7217026fcfa0b9895d691ba3ac3f40a9`.
The existing four-way cross-loading proof applies unchanged: both programs load both builders'
indexes and reproduce **212,719 nonzero ECs and 8,890,219 aligned pairs** on all
**10,261,737 TT-seq pairs**, with zero differences. Quantification was not rerun for these
scratch-only changes; the new proof is complete index-byte identity.

The [measurement evidence](benchmark-results.json) identifies the current construction binary
and source separately from the frozen quantification build, and records inputs, individual
construction times, exclusive stages, peak memory and existing cross-loading checks.

## Correctness and interpretation

Compare complete EC count maps with `compare_ec_counts.py`, independent of EC numbering.
Match index bytes and transcript ordering. Missing, malformed or empty comparison inputs fail
explicitly. The native HDF5 checker verifies IDs, lengths, counts, bootstrap metadata/values and
fragment histograms. Candidate TSV values must match HDF5 exactly; official TSV values are
checked at its six-significant-digit precision. `compare_abundance.py` reports numeric
differences and rejects missing/empty inputs, duplicate IDs and mismatched transcript sets.

Every measured bulk replicate is retained and audited separately. Absolute errors for estimated
counts, TPM and effective length include the maximum, p99, p99.9 and p99.99, using linear
interpolation at `(N-1)*p`. Reports include all transcripts (with zeros) and transcripts expressed
in either output. `compare_abundance.py --hdf5` compares native f64 counts/effective lengths and
derives TPM before text rounding; TSV comparisons remain separate. Replicate 1 is also compared
with replicates 2 and 3 within each program to expose fragment-sampling and numerical variation.
Exact corrected single-cell matrix/label equality establishes zero absolute error at every
percentile, both over the full matrix and its union of nonzero entries.

[Fragment sampling](fragment_length_sampling.md) can change effective lengths, EM convergence
and ambiguous isoform estimates even with identical ECs. Compare saved histograms when
interpreting model/EM speed. Fixed EC/model replays and matched iteration settings isolate EM
behavior. Bulk poly(A) libraries can assess reference cost and compatibility; biological
nascent-RNA accuracy requires suitable assays and independent evidence.

### Quantification model audit

The audit uses three retained eight-thread paired runs per tool on each complete fastp-trimmed
library: SRR13638690 (4,411,003 pairs), SRR1039508 (22,176,801), SRR307903 (32,460,465) and
SRR10919051 (10,172,656). The first three use GENCODE v50 ALL; the last uses the official
Pachter Lab mature-plus-nascent reference. It compares all nine cross-tool run pairings and all
three within-tool pairings per program, using native HDF5 values. These are correctness checks;
the fresh full-library runs do not supply repeated performance estimates.

| Library | Largest cross-tool count difference | Largest original-within-original count difference | Largest cross-tool / original-within-original TPM difference |
| --- | ---: | ---: | ---: |
| SRR13638690 | 1.080 | 1.080 | 68.983 / 68.983 |
| SRR1039508 | 5.277 | 5.277 | 108.807 / 108.807 |
| SRR307903 | 22.754 | 22.754 | 28.259 / 28.259 |
| SRR10919051 | 2.996 | 2.996 | 65.558 / 65.558 |

Within each cDNA library, the three retained kallistors runs agree to numerical precision;
the nascent runs have modest sample-dependent variation. Replaying all ten distinct fragment
histograms with the current quantifier reproduces every saved original/current output. All
effective lengths and iteration counts match. Maximum residuals across all transcript rows are
`1.921e-9` estimated counts and `4.221e-10` TPM. Fresh full fibroblast and nascent mappings also
preserve complete counted EC maps, and their quantifications independently match a retained
original run when their sampled histograms match.

The largest relative count difference among transcripts reaching 10 estimated counts in either
output is 1.64%, using absolute difference divided by the larger estimate. The 15.46% relative
TPM outlier has only about 0.05 estimated counts. Large absolute TPM outliers are mainly explained
by short-transcript effective lengths and normalization. These maxima describe these observations;
three repeats do not establish distributional equivalence or bound other inputs.

CALR-201 and CALR-233 have identical observed EC incidence in airway and fibroblast. Controlled
fixed-update replays separate the small model change from six additional EM updates in the
fibroblast case. The [sampling guide](fragment_length_sampling.md#interpreting-transcript-level-differences)
explains this ambiguity and the existing all-observation option. Portable model fingerprints,
error percentiles and the fresh checks are recorded under `quantification_variation_audit` in
[the measurement evidence](benchmark-results.json). The audited mode is default unbiased
short-read paired quantification without bootstrap sampling.
