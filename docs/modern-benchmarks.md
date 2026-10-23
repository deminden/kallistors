# Modern RNA-seq quantifier comparison

This protocol compares kallistors with current Salmon and piscem workflows.
It measures performance and quantification separately. Existing original-kallisto
controls are documented in [the benchmark report](benchmarks.md).
No kallistors implementation changes are part of this experiment.
Completed measurements are in the [current results](modern-benchmark-results.md).
A fresh original-kallisto control cohort follows both modern reference panels;
its design is described under [Original-kallisto refresh](#original-kallisto-refresh).

## Complete serial comparison

Run the entire panel with the locked analysis interpreter after building the
release executable:

```sh
python scripts/run_full_comparison.py --workspace target/benchmarks/full-comparison.noindex
```

The coordinator freezes the implementation and executable, collects the GENCODE
and official mature-plus-nascent panels, exports their tables, then runs the full
original-kallisto comparison. Index construction uses eight workers only;
paired-end, single-end and single-cell mapping use one and eight. Every condition
has one warmup and three measured observations. The complete alevin-fry pipeline
uses eight workers because its downstream stages require at least two.
Cell Ranger remains deferred.

`--tools-from WORKSPACE` reuses verified pinned installations.
`--reuse-gencode-inputs WORKSPACE` reuses checksum-verified reference sequences,
prepared reads and independent simulation truth. These options reuse inputs and
software; every performance observation is newly collected. Completed trials
retain their original immutable identities when the same command is resumed.
The coordinator records its commands, protocol and progress in the workspace.

The original-kallisto analysis compares native HDF5 counts and derived TPM for
every measured repetition at transcript and gene level. It reports maximum
absolute differences and p99, p99.9 and p99.99 over both all features and the
nonzero union. Separate tables label the largest transcript differences, record
realized fragment histograms, compare repetitions within each tool, and compare
one versus eight workers. Nascent analyses additionally separate mature,
nascent and gene-by-state estimates. Exact EC and weighted BUS comparisons remain
independent checks; real-data abundance agreement is not an accuracy ranking.

## Pinned software

| Component | Version | Role |
| --- | --- | --- |
| kallistors | 0.5.1 development snapshot | Exact executable SHA-256 and complete source snapshot recorded |
| [Salmon](https://github.com/COMBINE-lab/salmon/releases/tag/v2.8.0) | 2.8.0 | Sketch and selective alignment, followed by native inference |
| [piscem](https://github.com/COMBINE-lab/piscem/releases/tag/v0.23.0) | 0.23.0 | Bulk and single-cell mapping |
| [piscem-infer](https://github.com/COMBINE-lab/piscem-infer/releases/tag/v0.7.0) | 0.7.0 | Bulk abundance inference |
| [alevin-fry](https://github.com/COMBINE-lab/alevin-fry/releases/tag/v0.18.3) | 0.18.3 | Barcode correction, collation and UMI resolution |
| [bustools](https://github.com/BUStools/bustools/releases/tag/v0.45.1) | 0.45.1 | Complete kallistors BUS-to-matrix workflow |
| fastp | 1.4.0 | Shared read preparation |
| [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-release-notes) | 10.1.0 | Execution deferred at the user's request |

[`tools.lock.json`](../scripts/benchmarks/modern/tools.lock.json) pins release
assets by SHA-256, architecture and source revision. Native Apple Silicon and
Linux distributions are used, with no emulation. The analysis environment pins
Python 3.14.7, NumPy, SciPy, h5py and psutil, including wheel hashes. Results must
identify the development snapshot by its source and executable hashes; the
version string alone is insufficient. Kallistors uses the repository release
build settings, including `target-cpu=native`; comparator executables are the
pinned official native distributions. Their microarchitecture compiler settings
are not asserted to match. This is a comparison of these recorded executable
builds, not an experiment isolating compiler effects.

## Reference and input panels

The primary reference is **GENCODE v50 ALL**, with **670,670 transcripts** and
**86,411 gene identifiers**. All tools receive the same FASTA, target order and
transcript-to-gene table. Transcript versions are retained; gene versions are
removed after checking that this creates no collision. No biotype selection,
sequence deduplication or poly-A clipping is permitted. Salmon uses
`--keepDuplicates --no-clip`.

Thirteen sequences contain ambiguous bases. They are normalized once, before
any index construction, by replacing each non-ACGT character with
`ACGT[position mod 4]`. This avoids different tools making different random
replacements. The source and normalized FASTA hashes are recorded. All builders
use k=31; Salmon and piscem use m=19, and kallistors uses its native minimizer
setting. Different index structures are intrinsic to the methods.

| Panel | Input | Retained read pairs | Purpose |
| --- | --- | ---: | --- |
| Bulk, smaller library | SRR13638690 | 4,411,003 | End-to-end performance and agreement |
| Bulk, larger library | SRR307903 | 32,460,465 | Scaling with library size |
| Independent fragment-origin control | Seed 48321 | 4,000,000 | Accuracy under a stated generative model |
| Standard 10x Chromium 3′ v3 | PBMC 1k, both sequencing lanes | 65,893,590 | Mapping and complete gene-by-cell matrices |
| Nascent-enriched TT-seq, separate reference panel | SRR10919051 | 10,172,656 | Mature-plus-nascent reference scaling and agreement |

Each bulk library is tested as both paired-end and single-end at 1 and 8 threads.
Single-end inputs use only the exact prepared R1 files from the corresponding
paired-end library; they are not independent biological replicates. The source
records do not report both fragment mean and SD. Use one frozen empirical estimate
per real library, shared by every method:

| Library | Mean, bases | Population SD, bases | Observations |
| --- | ---: | ---: | ---: |
| SRR13638690 | 144.2279 | 50.140709624 | 10,000 |
| SRR307903 | 190.1364 | 22.090246604 | 10,000 |
| SRR10919051, nascent-enriched TT-seq | 135.636 | 39.382474579 | 10,000 |
| Fragment-origin simulation | 250 | 30 | Known generative parameters |

Real-library moments come from the stored paired-fragment histogram of the first
single-worker original-kallisto control using the identical fastp-prepared inputs.
These existing observations are reused for calibration, separately from the
subsequent original-kallisto performance refresh. The histogram, estimator version, source HDF5 checksum,
input checksums and reproduction command are pinned in
[`library-models.json`](../scripts/benchmarks/modern/library-models.json). The workflow
recomputes the moments from those counts and rejects mismatched prepared inputs.
They are empirical estimates, not measurements reported by the sample submitter
or independent biological truth. Original sampling/1000-bin limits apply, and
reducing the histogram to a Gaussian mean/SD is an approximation. Real-data
abundance differences remain agreement results rather than accuracy claims.
Native Gaussian truncation and binning remain tool-specific; supplied parameters
are common, while reported realized moments can differ. Native model metadata
are retained for audit. Simulation R1 retains the paired control's known
fragment-origin counts and nominal molecule TPM.

### Nascent-enriched reference panel

The TT-seq extension uses only the official Pachter lab
[mature](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_nac.mature.fa.gz)
and [nascent](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_nac.nascent.fa.gz)
FASTA assets from release v1: **227,665 mature and 39,546 nascent targets**
from Ensembl 108. These are a separate reference, not introns added to GENCODE v50.
No intronic sequences are reconstructed. Asset sizes and SHA-256 checksums are
pinned; every tool receives the same canonicalized target sequences and order.
The shared ambiguity policy above applies, with any replacements counted in the
reference manifest. Gene identifiers and mature/nascent labels come from the
official FASTA headers.

[GSM4277828 / SRR10919051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4277828)
is a HEK293 control labelled with 4-thiouridine for one hour. It enriches newly
synthesized RNA and contains both mature and unspliced transcripts. The
[TT-seq methods, section 2.13](https://pmc.ncbi.nlm.nih.gov/articles/PMC8780909/)
name the NEBNext Ultra RNA Library Prep Kit, which the
[manufacturer describes as non-directional](https://www.neb.com/en-gb/products/e7530-nebnext-ultra-rna-library-prep-kit-for-illumina).
The separate stranded mRNA assay in section 2.12 does not define this TT-seq library.
The common primary settings are therefore unstranded: `IU`/`U` and kallistors'
unstranded default. The study used an S2 RNA spike-in; all primary tools receive
the same human-only reference, so this panel does not quantify the spike-in.

Both paired-end and matching R1-only inputs are included at 1 and 8 workers,
with one warmup and three measured runs. Index construction uses 8 workers only.
The empirical single-end calibration is reused from the first stored one-worker
paired-end original-kallisto observation on the identical fastp-prepared files,
as for the cDNA libraries. It is pinned in
[`nascent-library-models.json`](../scripts/benchmarks/modern/nascent-library-models.json).
The nascent panel runs serially in its own frozen workspace, after the cDNA and
single-cell cohort; it cannot change their measurement protocol.

Agreement includes transcript estimates, gene sums, gene-by-state sums, and
mature and nascent target subsets. Model-based mature/nascent abundance totals
are reported separately. These are not known splicing truth, RNA velocity or
transcription-rate measurements; no biological accuracy ranking follows from
their agreement.

### Shared input preparation

ENA URLs, checksums and retrieval metadata are reused from
[`human_datasets.json`](../scripts/benchmarks/human_datasets.json). The 10x library
comes from the [manufacturer's PBMC 1k v3 dataset](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0).
The manufacturer reports 1,222 cells; this is contextual information, **not known
cell-count truth**. The 100,000-pair prefixes are format and command pilots only;
they must not enter the main performance tables.

All real libraries use the existing fastp preparation. For 10x data, R2 is
trimmed and filtered while the matched R1 barcode and UMI sequence is retained
unchanged. Both pipelines receive the same surviving pairs and lane order.
Shared downloading, reference normalization and trimming are preparation steps
outside the quantifier timing. Thus the complete-pipeline scope is **prepared
FASTQ to count matrix**, not an unqualified raw-data-to-analysis runtime.

The independent control samples 50,000 transcripts of at least 400 bases from
the complete reference without selecting for unique mappability. Molecule
abundances are lognormal; fragment lengths follow a clipped normal distribution
centered at 250 bases. PE100 reads have random fragment orientation and 0.1%
substitution errors. PCG64 seed 48321, exact fragment-origin counts, nominal
molecule TPM and the full model are saved. These adapter-free Q30 reads do not
need trimming. This is a deliberately limited accuracy control: it omits indels,
polymorphism, contamination, PCR and sequence biases. It cannot establish a
universal accuracy ranking. Real-data agreement is never labelled accuracy.

## Execution design

* Run bulk quantification methods **serially on the same host**, at **1 and 8
  requested worker threads**, with **one warmup followed by three measured runs**.
  **Index construction uses 8 threads only**, with the same warmup and repetitions.
* Rotate method order deterministically between repetitions. Retain all three
  measurements; do not discard slow runs or stop early based on results.
* Record OS, architecture, power/swap context, competing processes, executable
  hashes, source snapshot, complete commands, logs and workflow hashes.
* Use the final successful 8-thread index for quantification at both thread
  counts. Salmon shares one index between its two mapping modes; piscem shares
  one index between bulk and single-cell mapping. Do not charge a shared index
  twice.
* Report medians, observed minima and maxima, and n=3. Three repetitions quantify
  technical variability on this host, not variability across biological
  samples, machines or transcriptomes. Avoid universal speedup claims.

**Single-cell mapping uses 1 and 8 worker threads**, as do the complete
kallistors/bustools and original kallisto/bustools pipelines. Each condition
receives one warmup and three measured runs. Piscem's one-worker result covers
mapping only; its complete alevin-fry pipeline uses eight workers.
Ordinary stranded BUS input uses parallel counting with the strand policy applied
before caching transcript sets. Diagnostic traces and placement output are not
constructed for this counting-only workload.
Alevin-fry 0.18.3 enforces a minimum of two threads in permit-list generation,
collation and quantification. Bustools correction and counting are single-threaded stages
even in a pipeline with an 8-thread budget.


On Apple Silicon these settings specify tool worker counts, not CPU affinity.
Auxiliary I/O and runtime threads may exist. CPU seconds and CPU/wall ratios
should accompany any discussion of parallel scaling.

### Performance accounting

Native `/usr/bin/time` records user CPU, system CPU and peak RSS for each command.
Complete sequential pipeline wall time is the sum of the individual command wall
times; CPU time is summed and peak RSS is the maximum stage peak. The wrapper's
elapsed time, including monitoring and metadata collection between commands, is
also retained separately. Do not sum stage medians to estimate pipeline time:
the primary median is calculated from each complete repetition's total.

A separate 50 ms sampler records simultaneous process-tree RSS. It can miss brief
peaks and double count shared pages, so it is labelled separately from native
peak RSS. These measures must not be confused with proportional set size.
After validated read outputs are available, temporary BUS and RAD products are
hashed and removed outside the timer. Count matrices, abundance tables, labels,
native metadata and logs are retained. This bounds disk use without changing
what is timed.

Final index logical and allocated disk sizes are recorded after construction,
outside the timer. They are not a measurement of peak temporary disk demand.

The mapping table measures native mapping **including index loading, input
reading/decompression and intermediate output**. Kallistors writes aggregate EC
counts; Salmon and piscem write per-fragment RAD representations. These are
useful native workflow measurements, not identical output workloads or a pure
mapping-kernel microbenchmark. For bulk kallistors and Salmon, mapping-only
commands are measured separately from full quantification; piscem's isolated
mapping command is already the first stage of its full pipeline. For single
cell, native BUS and RAD mapping are the first pipeline stages.

### Algorithms and correction policies

Bulk strand settings follow the library metadata and are identical across tools.
The [ENCODE BJ library](https://www.encodeproject.org/experiments/ENCSR000COP/)
explicitly reports reverse strand specificity. SRR13638690 uses TruSeq RNA Exome,
according to its [experiment record](https://www.ebi.ac.uk/ena/browser/view/SRX10035917)
and [publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC8067552/); reverse
orientation is inferred from that first-strand/second-strand-marking protocol and
checked against the read prefix. Both cDNA libraries therefore use
`--rf-stranded` in kallistors, `ISR` for paired-end and `SR` for single-end in
Salmon and piscem-infer. The random-orientation simulation uses `IU`/`U` and
kallistors' unstranded default. Protocol sources and the distinction between
reported and protocol-inferred orientation are retained. Library-type warnings
are kept; no tool may auto-detect its strand setting while another is forced to
a different setting.

Salmon's primary selective-alignment panel uses the same transcript-only
reference as the other methods, **without genomic decoys**. It must be described
as transcriptome-only selective alignment. It does not evaluate Salmon's
recommended decoy-aware reference against contamination or unannotated genomic
transcription. Bias correction and bootstraps are off in all primary bulk runs.
Salmon's native deterministic VB inference, kallistors' EM and fragment sampling,
and piscem-infer's EM, range factorization and fragment estimation retain their
native defaults. Their statistical objectives and stopping criteria differ;
output equality is not expected. Salmon writes 15 significant digits for abundance
comparison, avoiding its default three-digit count rounding as a source of apparent
algorithmic disagreement. Full command and native metadata retention
makes these differences auditable.

For single cell, both workflows use the same 10x v3 barcode on-list and unchanged
16-base barcode/12-base UMI fields. Alevin-fry uses unique barcode correction in
a Hamming-1 neighborhood, unfiltered-list minimum 10 reads, and supported
`cr-like` resolution. Its pinned release rejects `--umi-edit-dist 1` with this
resolver, so that unsupported combination is not used. Tiny-cell resolution
behavior is retained and recorded in native `quant.json`. Bustools performs
sort, on-list correction, a second sort and gene counting. Its UMI and multimapper
rules differ from alevin-fry's; do not describe these as identical correction
algorithms.

The primary single-cell analysis uses explicit forward-transcript strand filtering:
`--fr-stranded` for kallistors and `-d fw` for alevin-fry. Kallistors loads the
required orientation annotations from the existing v13 index; unstranded BUS
keeps the compact loader. Piscem stores orientation-labelled mappings before
alevin-fry applies the forward-orientation policy in its permit-list/collation
workflow. Kallistors applies that policy during BUS mapping. Mapping-stage output
counts have different filtering scopes; the complete matrix panel uses the same
forward policy. The CLI exposes `--single-cell-strand fw|both`; a
separate unstranded panel must use `both` consistently in both workflows and be
labelled explicitly. Validate commands on the prefix pilot before full runs.

Strand-filtering verification on 100,000 prepared PBMC pairs gives exact weighted
barcode/UMI/transcript-set/flag/count parity with original kallisto: 62,033 forward,
10,294 reverse and 64,899 unstranded mappings. This verifies software compatibility,
not biological accuracy or full-library performance. Toy controls additionally
cover paired BUS and BAM, both index builders, and 1/8-thread processing. No index
rebuild is required for the strand fix. A changed executable requires a new frozen
benchmark workspace; do not combine read measurements from different snapshots.

### Agreement and accuracy

Bulk abundance tables must contain exactly the shared transcript universe.
Compare transcript counts and TPM, then counts and summed transcript TPM at gene
level. Report mean absolute difference, maximum absolute difference, p99, p99.9,
p99.99, log1p Pearson correlation and Spearman correlation. Agreement uses full
feature sets, including zero estimates. Simulation tables compare against known
fragment origins and nominal expression, separately from pairwise agreement.

Canonicalize both matrices to genes × barcodes, checking label uniqueness and
matrix dimensions. Report total UMIs, detected genes and raw barcode counts. A
common descriptive cell gate requires **at least 500 UMIs and 200 genes** in both
workflows; report called-cell overlap/Jaccard and per-cell UMI/gene summaries.
This gate is not an EmptyDrops significance test or a substitute for biological
cell calling. Report differences over the full union of genes and barcodes;
include nonzero-union percentiles because all-entry percentiles are dominated by
implicit zeros. Also compare the matrix within cells passing the common gate in
both pipelines, and report per-barcode UMI and detected-gene differences there. Native and common-gate cell counts must remain distinct.

Measured replicate 1 is selected in advance for cross-tool agreement; it is not
chosen for favorable speed or agreement. Replicate-to-replicate variation is
reported separately. Biological truth is unavailable for the PBMC library, so
UMI and cell differences cannot establish which method is more accurate.

## Reproduce

Build the kallistors snapshot intended for comparison before freezing it:

```sh
cargo build --release --locked
python3.14 -m venv target/modern-analysis
target/modern-analysis/bin/python -m pip install --only-binary=:all: --require-hashes \
  -r scripts/benchmarks/modern/analysis-requirements.lock
target/modern-analysis/bin/python scripts/modern_benchmarks.py all
```

Use the exact Python patch version in `analysis.lock.json`. Native distributions
are downloaded and checked by `bootstrap`; existing pinned fastp and bustools
builds are reused when present. Source builds require a C/C++ compiler, CMake,
zlib and the fastp build dependencies. An existing installation can be supplied
with `--tool NAME=/path/to/executable`; the version is checked and an immutable
copy is frozen. All commands accept `--workspace` and `--data` for portable data
and result locations.

The `all` action freezes software and inputs, collects the 8-thread index panel,
validates prefix pilots, collects the full read panel, and exports tables. Pilots
reuse checksum-verified primary indexes when their tools and reference match,
so index construction is not repeated merely for command validation.

Run the separate mature-plus-nascent panel with the same software pins:

```sh
target/modern-analysis/bin/python scripts/modern_benchmarks.py \
  --workspace target/modern-nascent --reference official-nac all
```

Run the two reference panels serially. The nascent command selects only the TT-seq
library and its R1 counterpart, rather than timing unrelated libraries again.

After both panels finish, generate the repository report and publication tables:

```sh
target/modern-analysis/bin/python scripts/publish_modern_benchmarks.py \
  --gencode target/benchmarks/modern-comparison.noindex \
  --nascent target/modern-nascent
```

Publication rejects incomplete warmup/triplicate groups and interrupted timers.
It exports CPU, wall-time and RSS tables, index sizes, agreement and accuracy,
and processed-record audits. Analysis sources are retained separately from the
frozen measurement workflow.

Individual actions are `bootstrap`, `retrieve`, `prepare`, `plan`, `pilot`, `run`
and `report`. `retrieve` reuses verified existing preparation records. `pilot`
runs separate prefix fixtures with no publication claims. `run --only index`
executes the index panel; `run --only SRR13638690,SRR307903,fragment-truth` executes
the paired-end bulk panel against the completed indexes. Append `-se` to each
bulk case ID to select its single-end counterpart; omitting `--only` includes all
bulk modes and the full single-cell library. Workflow snapshots are retained by
hash. An explicitly audited orchestration revision can be applied after index collection and before any read measurements
with `--amend-unstarted-workflow`; the timer, tools and reference must remain
unchanged. If read adapters are extended before read
collection, every completed index command and worker count must still match
exactly; this proof and the earlier specification are retained.

A successful run resumes from saved records;
changed tools, workflow or completed-case inputs are rejected. Interrupted output
is archived separately. An unstarted case can be amended after a pilot, with the
previous specification retained for audit.

For an executable correction, recollect every measurement of the changed tool.
Unchanged competitor controls may be retained only through an explicit audit of
their commands, executable hashes, inputs and measurement code. Retained records
keep their original identities and workflow snapshots; corrected measurements
use a distinct frozen executable. Publication checks reject commands that use
an executable outside the final experiment manifest.

The result directory contains frozen binaries and sources, reference and input
manifests, per-command logs/accounting, individual run records and tables.
`report` exports TSV, Markdown and booktabs LaTeX tables. Incomplete triplicates
are excluded from primary tables and listed explicitly. A missing/deferred tool
is never assigned zero runtime or zero memory.

## Deferred Cell Ranger panel

Cell Ranger 10.1 requires x86-64 Linux. Its execution is deferred, with no
emulation or cross-host speed ratio. The extension must build a matched GENCODE
v50 genome/GTF reference, use 8 cores, and compare exon-only counting with
`--include-introns=false`, without BAM or secondary-analysis work in the measured
counting scope. Genome alignment remains an algorithmic difference from
transcriptome mapping. Installation/archive checksums, exact reference artifacts,
index construction and both native and common-gate matrices must be captured on
that host. This extension has **not** been executed or validated here.

## Original-kallisto refresh

After both modern reference panels and their analysis finish, run fresh controls
against the pinned official [kallisto 0.52.0 release](https://github.com/pachterlab/kallisto/releases/tag/v0.52.0).
Use the identical frozen kallistors executable, references and prepared reads.
The refresh also includes the existing airway library, SRR1039508: **22,176,801
fastp-retained pairs**. Its ENA record names TruSeq RNA Sample Prep Kit v2, so the
common setting is unstranded. Its matching paired-end calibration is pinned in
[`kallisto-extra-library.json`](../scripts/benchmarks/kallisto-extra-library.json).
All four real bulk libraries and the independent simulation are tested as
paired-end and R1-only input at 1 and 8 workers.

Fresh index builds use 8 workers only, with one warmup and three measurements
for each builder on each reference. Primary read controls share the final
original-built index bytes. Separate matched-prefix checks use each freshly
kallistors-built index, checking that both programs can read it with identical
EC counts. These prefixes are compatibility fixtures, not primary speed results.
Full-library EC preflight is mandatory for every primary bulk configuration.

The complete PBMC library is tested at 1 and 8 workers for both mapping and the
same prepared-FASTQ-to-matrix bustools pipeline. Every warmup and measured output
must have exact weighted barcode/UMI/transcript-set/flag/count parity and exact
canonical gene/UMI matrices. Raw intermediates are retired only after full checks.

Every primary condition uses one warmup and three measured runs, alternating
tool order. Record wall and CPU time, peak RSS and index disk size. Abundance
agreement uses native HDF5 f64 values, rather than original kallisto's rounded
TSV counts, and includes transcript, gene and mature/nascent state comparisons,
maximum absolute differences, p99, p99.9 and p99.99. Replicate variation,
one-versus-eight-worker agreement and independent simulation accuracy remain
distinct. Exact EC parity does not imply identical sampled fragment models.
Bulk runs additionally record CLI-message phase boundaries, enabling estimates
of startup/index loading, input setup, read processing, model setup, EM and output.
These are wall-time boundary estimates, not isolated mapping-kernel measurements
or per-stage CPU/RSS. The model-and-EM total is a derived aggregate and must not
be added to its components.

```sh
target/modern-analysis/bin/python scripts/run_kallisto_comparison.py \
  --workspace target/original-comparison \
  --gencode target/benchmarks/modern-comparison.noindex \
  --nascent target/modern-nascent \
  --kallisto target/tools/kallisto-0.52.0 \
  --diagnostic-kallisto target/tools/kallisto-src/build/src/kallisto
```

The official release supplies timing controls; the separate instrumented build
supplies EC correctness exports only. See [diagnostic build instructions](kallisto_debug_build.md).
The optional `--wait-for` accepts a preceding pipeline status JSON and waits for
`complete`; a preceding failure holds the new cohort for inspection. The runner
freezes its own tools and workflow and publishes `kallisto-comparison-results.md`,
machine-readable results, and TSV/LaTeX tables only after the full panel passes.
