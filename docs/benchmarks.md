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
result writing. The [ENA manifest](../scripts/benchmarks/human_datasets.json) defines all five
public datasets and their directories. The [suite manifest](../scripts/benchmarks/quantification.json)
defines larger-library/reference combinations; the runner contains no accession-specific plans.

## Retrieve and trim reads

```bash
cargo install --locked fasterp --version 0.2.1
python3 scripts/prepare_benchmark.py --samples SRR13638690
```

`--samples all` selects all five manifest entries, including the two larger libraries.
`--samples SRR1039508,SRR307903` selects those two. The manifest pins mate URLs and MD5s;
compressed sizes are checked when supplied. Downloads support resumption and gzip integrity is
checked. `--data DIRECTORY` changes the data root; `--fasterp PATH` selects the trimmer.

| Run | Source | Deposited pairs |
| --- | --- | ---: |
| SRR13638690 | Temporal muscle, PRJNA699558 | 4,451,399 |
| SRR1039508 | Untreated airway smooth muscle N61311, PRJNA229998 | 22,935,521 |
| SRR307903 | ENCODE BJ fibroblast poly(A)+ RNA, PRJNA30709 | 36,010,036 |

Preparation enables fasterp adapter detection, front/tail trimming with a four-base Q20 window,
minimum length 50, qualified-base Q15, at most 40% unqualified bases and five Ns. Defaults are
four threads and gzip level 4. `trimmed-fasterp/preparation.json` records settings, trimmer
version/hash, input checksums, output hashes and actual before/after pair counts. Cached outputs
must match those records. The benchmark verifies prepared-read hashes before measuring.

A pair is two FASTQ records. Retrieval and trimming are outside quantification timing.

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
source and output hashes and target/base counts. This smaller reference supports iterative
index-build measurements: k=31 and g=23 produce 109,981,201 distinct k-mers and
1,017,525 unitigs. It does not replace full-reference compatibility checks.

### Intronic reference: Pachter Lab human NAC

```bash
python3 scripts/prepare_reference.py --reference official-nac
target/tools/kallisto-0.52.0 index -t 8 -k 31 -m 23 \
  -i data/official-human-nac/human_nac.mature-plus-nascent.idx \
  data/official-human-nac/human_nac.mature-plus-nascent.fa
```

Intronic benchmarks use only Pachter Lab's published [mature](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_nac.mature.fa.gz)
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

Prepare the full libraries, subsets and references first, then select a manifest group:

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

`--group gencode50` selects cDNA-only cases; `--group official-nac` selects all intronic
cases using the Pachter Lab mature-plus-nascent reference. A suite case names its accession,
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

On Apple Silicon, quantifying **4,410,990 paired fragments** from SRR13638690 against the complete
GENCODE v50 ALL index (670,670 transcripts, k=31, g=23) took a median **8.80 s**, with **1.15 GiB peak RSS**, and
**1.10 GiB peak footprint**. Reads use the documented ENA/fasterp preparation. One full
warmup per binary preceded three measured runs at eight threads, alternating with the preceding
committed Rust baseline on identical inputs and index. Source builds used `target-cpu=native`.
All measurements were retained; the CPU monitor flagged activity in 0 candidate and
0 control runs. Such sampling cannot exclude all background contention.

| Measurement | Preceding baseline | Current implementation |
| --- | ---: | ---: |
| Whole-process median | 9.272 s | 8.799 s |
| Median peak RSS | 2.069 GiB | 1.146 GiB |

The fresh comparison reduced median peak RSS by **44.6%** with a
**1.054×** total speedup. Ten accepted memory optimization rounds each required at least
2% lower median peak RSS and no higher median total time against the immediately preceding
accepted binary. Partial improvements were combined; failed candidates did not count.

| Phase | Current median seconds |
| --- | ---: |
| Index loading/setup | 1.561 |
| Read processing | 4.900 |
| Model + EM | 1.803 |
| Output/shutdown | 0.533 |
| Complete process | 8.799 |

Totals include index loading, every read pair, fragment estimation, EM, full default TSV/HDF5
output and cleanup. Deferred transcript-name loading is included in output/shutdown. Stage
medians are independent and need not sum to the total median. Peak RSS is the operating system's
whole-process high-water mark, not the live size of the final graph. macOS footprint is a separate
measurement.

Complete paired EC comparison preserved **190,693 nonzero classes and 4,246,399 counts** against
both the preceding Rust baseline and original kallisto's diagnostic export. A full-mate single-end
comparison preserved **191,262 classes and 4,276,198 counts** against the Rust baseline. Native HDF5
readback validated every one of the 670,670 full transcript IDs, lengths and output rows.

Default fragment sampling retains upstream buffer/quota behavior. Histograms can differ with
scheduling, so exact EC parity does not imply exact default abundances. With
`--fragment-length-samples 0`, the matched all-observation model produced identical histograms,
target IDs, target lengths and effective lengths. The maximum estimated-count difference was
1.1e-09. Fixed-input serial/threaded regressions also check floating-point bits. These are
measurements of this workload and build, not a universal memory or speed estimate.

### Eight-thread index build measurement

On Apple Silicon, building the complete GENCODE v50 ALL transcript FASTA with eight threads:
670,670 transcripts and 207,228,052 distinct k-mers; k=31 and g=23. The input was exactly
decompressed before timing, with no transcript filtering.

| Measurement | Current kallistors |
| --- | ---: |
| Median whole-process elapsed time | 12.39 s |
| Median peak RSS | 5.08 GiB |
| Median peak footprint | 4.59 GiB |

One warmup preceded three measured builds. Source builds used `target-cpu=native`.
All measured runs were uninterrupted. Five-second sampling recorded a background process at
26% CPU in one measured build; alternating controls and all three observations were retained.
RSS and macOS
peak footprint are separate OS measurements. These values describe this reference and
machine family; they do not establish a universal timing estimate.

Median internal stages were 0.65 s for FASTA parsing, 9.19 s for graph construction,
0.89 s for EC construction, 0.97 s for minimizer/MPHF work and 0.61 s for writing.
Independent stage medians need not sum to the median whole-process time.

Construction excludes only byte-proven repeat windows, preserving seams and complete metadata.
Shared ends and sparse internal anchors select candidate repeats; full equality proves each skip.
Bounded deduplication and radix sorting limit collection work; completed parallel groups are
losslessly packed before their word sort and restored to transcript order. Canonical-edge symmetry, exact
small-bucket guards, active query batches and AArch64 cache hints reduce graph lookup work.
The requested worker pool is reused across stages; edge finalization, accepted
branch sequences, ECs, minimizers, collision marking and grouped serialization run in parallel.
Serial ordered claims/joins retain deterministic graph order; long residual paths keep bounded
fallback handling. Compact path/link/location scratch, sequence ownership reuse and early lookup
release limit temporary allocations. Exact packed rank bounds, counted start-table filling,
contiguous grouped EC positions and narrow MPHF scratch avoid overlapping full-size buffers.
Wholly unused scratch pages are reclaimed on macOS/Linux. All tied minimizers and the v13 representation are preserved.
The internal total includes final cleanup; process startup and pool teardown remain in whole-process
time. See [construction details](development.md#index-construction) for exactness guarantees.

Compatibility validation used 99,029 fasterp-trimmed pairs retained from the first 100,000
SRR13638690 pairs. All 24 complete EC maps matched within their mode/strand case:
paired/single × unstranded/FR/RF × both index builders × both quantifiers. Normal quantification
counters and native HDF5 readback passed for both indexes in paired and single-end modes.
Complete index-byte regression checks passed for both the protein-coding subset and full
GENCODE v50 ALL reference (670,670 transcripts and 207,228,052 distinct k-mers). Serial/threaded
construction also matched across the nucleotide k/g and amino-acid validation matrix.


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
