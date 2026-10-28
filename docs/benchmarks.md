# Benchmarks

Compare kallistors with original kallisto, Salmon, piscem, Fulgor and Themisto on prepared
reads and matched reference target sets. The original-kallisto quantification runner uses
identical index bytes and requires complete EC parity before timing. Comparisons with other
tools report their native mapping and inference differences. Runners save executable/input
hashes and validate retained outputs; results depend on workload, hardware and build settings.

## Latest completed results

Current benchmark evidence is organized by workload, with frozen tool versions and input
provenance identified in each report:

1. [Salmon and piscem comparison](salmon-piscem-results.md): full PE/SE quantification,
   native mapping, single-cell processing and index construction.
2. [Original kallisto comparison](kallisto-comparison-results.md): default bulk quantification,
   bootstrap modes, full PBMC BUS/matrices, index construction and separate sequence-bias results.
3. [Fulgor and Themisto mapping comparison](mapping-comparison-results.md): independent
   single-read mapping, native output differences, actual worker/parser counts and setup costs.

The reports link complete numeric exports and [frozen protocols](benchmark-provenance.json).
The recipes below explain preparation and individual runner panels; the frozen protocols and
[mapping command record](mapping-comparison-commands.json) specify the measured configurations.

## Tools

The individual preparation and timing scripts use Python 3.9 or newer. Native HDF5 validation
also needs h5py and NumPy; choose versions compatible with the interpreter:

```bash
python3 -m venv target/hdf5-native-check
target/hdf5-native-check/bin/python -m pip install h5py numpy
```

Put the official timing executable at `target/tools/kallisto-0.52.0`, or pass `--kallisto`.
Build a separate instrumented original as described in the [debugging guide](kallisto_debug_build.md).
Pass it with `--diagnostic-kallisto`; it supplies EC evidence and is never the timing control.

The full comparison suite has stricter requirements. `modern_benchmarks.py`,
`run_kallisto_comparison.py` and `run_full_comparison.py` enforce the
[analysis environment lock](../scripts/benchmarks/modern/analysis.lock.json): Python **3.14.7**,
NumPy 2.5.4, SciPy 1.18.1, h5py 3.16.0 and psutil 7.2.2. Use that Python patch version
to create a separate environment, then install the checksum-pinned dependencies:

```bash
python3.14 -m venv target/benchmark-analysis
target/benchmark-analysis/bin/python -m pip install --require-hashes \
  -r scripts/benchmarks/modern/analysis-requirements.lock
```

The suite verifies the interpreter and package versions before executing a panel. Tool versions
and supported platform downloads are pinned in the
[tool lock](../scripts/benchmarks/modern/tools.lock.json). Bootstrap reuses explicitly supplied
executables or retrieves pinned distributions; building fastp or bustools from source needs
their build dependencies, including a C/C++ toolchain and make/CMake.

| Command | Purpose |
| --- | --- |
| `prepare_benchmark.py` | Verified ENA/10x retrieval and fastp trimming |
| `prepare_reference.py` | GENCODE cDNA or Pachter Lab NAC FASTA preparation |
| `subsample_benchmark_pairs.py` | Exact seeded paired subsets |
| `benchmark_reads.py` | Quantification, thread scaling, optional stage timing and EC preflight |
| `benchmark_bus.py` | Complete single-cell BUS timing, counted-record parity and corrected gene/UMI matrices |
| `benchmark_index.py` | Index build time, size and peak RSS |
| `report_benchmarks.py` | Portable JSON/TSV export from completed runs |
| `modern_benchmarks.py` | Prepare, run and report the pinned Salmon/piscem comparison |
| `run_kallisto_comparison.py` | Compare original kallisto using the prepared reference panels |
| `run_full_comparison.py` | Run the alternatives and original-kallisto panels sequentially |
| `compare_ec_counts.py` | Independent comparison of complete counted EC maps |
| `compare_abundance.py` | Transcript abundance and effective-length comparisons |
| `check_hdf5_output.py` | Independent native HDF5 validation and fragment-model statistics |
| `diagnostics/prepare_quant_replay.py` | Fixed-EC, fixed-fragment-model inputs for EM investigation |

`benchmark_common.py` supplies shared execution, memory measurement, hashes, downloads and
result writing. The [ENA manifest](../scripts/benchmarks/human_datasets.json) defines four
public datasets and their directories. The [suite manifest](../scripts/benchmarks/quantification.json)
defines larger-library/reference combinations; the runner contains no accession-specific plans.

`modern_benchmark/` contains shared tool adapters, measurement, validation and analysis code.
`original_comparison/` reuses those helpers for original-kallisto controls. The diagnostic
patches and replay helper support the independent correctness checks described in the
[debugging guide](kallisto_debug_build.md).

### Full comparison suite

From the repository root, build the native candidate and prepare the additional airway library
used by the original-kallisto panel. Supply fastp 1.4.0 and the two original-kallisto executables
described above:

```bash
cargo build --locked --release --package kallistors --bin kallistors
target/benchmark-analysis/bin/python scripts/prepare_benchmark.py \
  --samples SRR1039508 --fastp /path/to/fastp
target/benchmark-analysis/bin/python scripts/run_full_comparison.py \
  --workspace target/benchmarks/comparison \
  --data data \
  --kallistors target/release/kallistors \
  --kallisto target/tools/kallisto-0.52.0 \
  --diagnostic-kallisto target/tools/kallisto-src/build/src/kallisto
```

The coordinator retrieves and prepares the other libraries and both references, runs the
Salmon/piscem panels sequentially, then runs the original-kallisto comparison. It includes
index construction, bulk quantification, single-cell processing and numerical checks;
the original panel also includes bootstraps and sequence-bias correction. Use a new workspace
for a new comparison and leave frozen inputs, code and executables unchanged during execution.
`status.json` records the current phase or a failure; final summaries record completeness.
Run one native measurement campaign at a time, without competing builds or profiling.

For a single reference, use `modern_benchmarks.py --workspace PATH --reference gencode50 all`
with the same analysis interpreter; `--reference official-nac` selects mature-plus-nascent.
Its individual actions are `bootstrap`, `retrieve`, `prepare`, `plan`, `pilot`, `run` and `report`;
`--tool NAME=PATH` supplies an executable to bootstrap. The Fulgor/Themisto extension is
documented separately in its [mapping report](mapping-comparison-results.md) and
[command record](mapping-comparison-commands.json). Use the
[portable templates](mapping-reproduction.md) to configure local tools, indexes and reads.

### Generated reports

Results stay in the selected benchmark workspace. `modern_benchmarks.py report` writes its
numeric summary and detailed tables under `tables/`; the original-kallisto runner writes
`summary.json`, `README.md` and `tables/` in its workspace. The full coordinator keeps these
under its `gencode/`, `nascent/` and `original/` subdirectories.

Repository documentation is updated separately from reviewed evidence. Running a benchmark
or generating its report does not copy results into `docs/`. Completed runs retain their
frozen workflow snapshots; use those snapshots to reproduce an existing campaign, and a
new workspace for changed scripts or executables.

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
must match those records. The benchmark verifies prepared-read hashes before measuring.

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
Ordinary BUS counting uses compact indexes and ordered parallel mapping for single-cDNA,
compound, tagged and paired FASTQ layouts, including explicit strand selection. Packed
single-cDNA input overlaps decoding, mapping and ordered output; other layouts and ordinary
BAM input normalize their required fields into bounded batches for the shared counting workers.
Single-cDNA reads can reuse exact full-query answers across workers and lanes. Long-read,
amino-acid, unmapped-read and BAM-output modes retain the metadata their features require.
See [mapping routes](development.md#mapping-routes) for implementation details.
Measure BUS separately from `quant`.

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
EM, normal TSV/HDF5 output and shutdown. Every measured output is retained and validated;
warmup outputs are removed after their records are saved. Counter differences and disagreement
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
  --thread-counts 1,8 --runs 3 --warmups 1 --stages \
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

Use one warmup and three measured builds per tool at eight threads; report the median.
Keep each observation and record background activity without selecting only the fastest run.
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
