# Local real-read comparison

This report records reproducible local comparisons of kallistors 0.4.1 and upstream kallisto
0.52.0 on three paired-end human RNA-seq runs. The compact run is useful for a quick parity check;
the larger runs provide progressively more work for multicore timing. All quantification comparisons
use the same upstream-built index, so the speed measurements focus on quantification.

The subsequent [quantification audit](quantification_mismatch_audit_2026-10-14.md) found material
paired-end effective-length differences despite matching EC count maps, plus the reproducible
two-read TTN discrepancy in the larger single-end sample. The count-parity columns below refer to
run counters, and do not establish transcript abundance parity.
Both defects are now corrected; the [fix report](quantification_fixes_2026-10-14.md) contains
full EC validation, abundance checks, and the new eight-thread timing results. Timing and abundance
values in this original report describe the earlier build.

## Software and reference

- R 4.6.1 in the `kallistors-r` Conda environment.
- fastp 1.4.0 for adapter and quality trimming.
- Rust dependencies were updated to their latest stable releases before the comparison; the
  resolved versions are pinned in `Cargo.lock`.
- Apple Silicon macOS host with 10 logical CPUs and 32 GiB RAM.
- Upstream kallisto v0.52.0 from the [official release](https://github.com/pachterlab/kallisto/releases/tag/v0.52.0).
- kallistors built in release mode from this checkout using the versions in `Cargo.lock`.
- GENCODE human release 50 (GRCh38.p14), the full `gencode.v50.transcripts.fa.gz` FASTA (ALL
  regions; not the protein-coding-only file). The official [release page](https://www.gencodegenes.org/human/)
  describes the transcript FASTA. The built index contains 670,670 targets.

## Datasets

| ENA run | Description | Raw pairs | Read length | Compressed mates |
| --- | --- | ---: | ---: | ---: |
| [ERR2260014](https://www.ebi.ac.uk/ena/browser/view/ERR2260014) | Small human paired-end RNA-seq run | 15,738 | 99 bp | 1.8 MB |
| [DRR023395](https://www.ebi.ac.uk/ena/browser/view/DRR023395) | Human pulmonary artery endothelial cells (HPAEC), study PRJDB1532 | 913,953 | 100 bp | 145 MB |
| [SRR13638690](https://www.ebi.ac.uk/ena/browser/view/SRR13638690) | Human temporal muscle RNA-seq, study PRJNA699558 | 4,451,399 | about 171 bp per mate | 871 MB |

ENA supplies per-file MD5s; `scripts/prepare_local_benchmark.sh` checks them after download and
checks each gzip stream. SRR13638690 is about five times the trimmed pair count of DRR023395 and
adds a materially larger workload while keeping its paired raw FASTQs below 1 GiB combined.

fastp was run on all three samples with paired adapter detection, Q20 sliding-window trimming at both
ends, a 50 bp minimum length, and four threads. It retained 9,045 pairs from ERR2260014,
867,962 pairs from DRR023395, and 4,411,003 pairs from SRR13638690. Reports are written beside
each trimmed dataset under `data/`.

## Prerequisites

Run these procedures from a source checkout. Build kallistors with
`cargo build --locked --release`. Download the executable for your platform from the
[official kallisto v0.52.0 release](https://github.com/pachterlab/kallisto/releases/tag/v0.52.0),
extract it, and place it at `target/tools/kallisto-0.52.0` with executable permission.
Both benchmark scripts also accept `--kallisto /path/to/kallisto` and
`--kallistors /path/to/kallistors` to use different binary locations.

The optional environment is created with `mamba env create -f environment.yml`. Data preparation
also works directly with fastp 1.4.0 on `PATH` and a working Python 3 installation; R is not used
by the comparison scripts.

## Reproduce data preparation

From the repository root, after activating or using the `kallistors-r` environment:

```bash
mamba run -n kallistors-r bash scripts/prepare_local_benchmark.sh
```

The script downloads all three ENA runs and the GENCODE v50 full transcript FASTA into ignored `data/`
paths, verifies the ENA MD5s and gzip streams, then writes fastp-trimmed paired FASTQs and reports.
The metadata and FASTQ locations can also be retrieved from the [ENA run file report API](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=DRR023395&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes,read_count,library_layout).

Build the shared index used for read comparisons separately; the index benchmark deletes its
temporary indexes after timing:

```bash
target/tools/kallisto-0.52.0 index -t 8 \
  -i data/srr-small/gencode.v50.kallisto.idx \
  data/srr-small/gencode.v50.transcripts.fa.gz
```

## Reproduce index-build comparison

The script uses the same transcript FASTA, k=31, and four threads for both tools. It performs one
warmup per tool and three measured builds per tool, alternating the tool order between measured
runs. Each build goes into a temporary directory; the raw timings and median summary are saved in
`target/benchmarks/local-real-2026-10-14/`.

```bash
python3 scripts/benchmark_local_index.py --threads 4
```

| Tool | Median build time | Median index size |
| --- | ---: | ---: |
| kallisto 0.52.0 | 139.03 s | 1,336,070,019 bytes |
| kallistors 0.4.1 | 123.16 s | 1,376,743,621 bytes |

Kallistors was 1.13x faster on the median of three builds, with an index about 3.0% larger.

## Reproduce paired-end and single-end quantification

Quantification uses one shared index built by upstream kallisto and four thread counts. For each
dataset, mode, and thread count, the script runs one warmup per tool followed by three measured
runs; tool order alternates across replicates. It reports median wall time, speedup, and parity for
`n_processed`, `n_pseudoaligned`, and `n_unique`. Single-end quantification uses the trimmed mate 1
file with `--single -l 200 -s 20` for both tools. The compact sample uses four threads; the larger
sample uses one, two, four, and eight threads. The speedup column is kallisto median time divided by
kallistors median time, so values above 1 indicate kallistors was faster.

```bash
python3 scripts/benchmark_local_reads.py \
  --samples ERR2260014_compact,DRR023395_larger \
  --thread-counts 1,2,4,8 \
  --compact-thread-counts 4 --larger-thread-counts 1,2,4,8
```

Raw run records and median summaries are written to `target/benchmarks/local-real-2026-10-14/`.
Those build and quant result files are ignored by Git so the repository does not acquire large
data or machine-specific timing artifacts. The saved `quant.runs.tsv` contains each warmup and
measured repetition; `quant.summary.tsv` contains the medians below.

| Dataset | Mode | Threads | kallisto median | kallistors median | kallisto / kallistors | Processed inputs* | Pseudoaligned | Unique | Counts match |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| ERR2260014 | Paired | 4 | 8.235 s | 9.263 s | 0.89x | 9,045 | 6,211 | 1,138 | Yes |
| ERR2260014 | Single | 4 | 15.987 s | 12.492 s | 1.28x | 9,045 | 6,587 | 993 | Yes |
| DRR023395 | Paired | 1 | 46.725 s | 84.056 s | 0.56x | 867,962 | 825,231 | 54,366 | Yes |
| DRR023395 | Paired | 2 | 38.248 s | 50.165 s | 0.76x | 867,962 | 825,231 | 54,366 | Yes |
| DRR023395 | Paired | 4 | 31.325 s | 32.049 s | 0.98x | 867,962 | 825,231 | 54,366 | Yes |
| DRR023395 | Paired | 8 | 31.747 s | 26.920 s | 1.18x | 867,962 | 825,231 | 54,366 | Yes |
| DRR023395 | Single | 1 | 58.196 s | 55.804 s | 1.04x | 867,962 | 818,562 | 46,441 | Yes |
| DRR023395 | Single | 2 | 40.811 s | 41.758 s | 0.98x | 867,962 | 818,562 | 46,441 | Yes |
| DRR023395 | Single | 4 | 33.579 s | 28.239 s | 1.19x | 867,962 | 818,562 | 46,441 | Yes |
| DRR023395 | Single | 8 | 29.407 s | 23.756 s | 1.24x | 867,962 | 818,562 | 46,441 | Yes |

*Paired mode reports processed read pairs; single-end mode reports processed mate-1 reads.

For DRR023395, both tools improved from one to eight threads in paired mode (46.73 to 31.75 s
for kallisto; 84.06 to 26.92 s for kallistors). In single-end mode, the medians fell from 58.20 to
29.41 s and from 55.80 to 23.76 s, respectively. Kallistors was slower at low thread counts in
paired mode, then faster at eight threads; the single-end results also varied by thread count. These
are local wall-time measurements on one host, with a shared upstream index, and should not be read
as a general performance claim.

## Larger SRR13638690 comparison at eight threads

Only eight-thread runs were requested for this additional dataset. For paired and single-end modes,
the command below performs one warmup and three measured runs per tool, then writes raw records and
medians to an ignored `target/` path. Single-end mode uses trimmed mate 1 with the same fragment
length settings documented above.

```bash
python3 scripts/benchmark_local_reads.py \
  --samples SRR13638690_larger --thread-counts 8 \
  --out-prefix target/benchmarks/local-real-2026-10-14/srr13638690-8thread
```

fastp retained 4,411,003 pairs from the ENA-reported 4,451,399. At eight threads:

| Mode | kallisto median | kallistors median | kallisto / kallistors | Processed inputs* | Pseudoaligned | Unique kallisto / kallistors |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Paired | 48.433 s | 80.051 s | 0.61x | 4,411,003 | 4,246,701 | 209,744 / 209,744 |
| Single | 59.122 s | 51.537 s | 1.15x | 4,411,003 | 4,246,922 | 205,821 / 205,819 |

*Paired mode reports read pairs; single-end mode reports mate-1 reads. Processed and pseudoaligned
counts matched in both modes. Unique counts matched in paired mode; in single-end mode kallistors
reported two fewer uniquely pseudoaligned reads (205,819 vs 205,821).

An additional untimed eight-thread single-end output comparison covered all 670,670 transcripts.
The TPM values were exactly equal for 604,161 transcripts. Among transcripts with positive kallisto
TPM, relative-error median, 90th percentile, and 99th percentile were 6.74e-7, 2.34e-6, and
3.97e-6; maximum relative error was 0.00145 (0.145%). Thus the outputs are very close by this
summary, but not bit-for-bit identical. The one-off abundance files are kept under the ignored
`target/benchmarks/local-real-2026-10-14/srr13638690-abundance-check/` directory.

The original ERR2260014 abundance summary used percentiles over the full transcript table, most
of which has zero estimated abundance. The subsequent audit found substantial differences for
expressed short transcripts, including 64,342.7 vs 49,999.4 TPM for ENST00000427161.1. Percentiles
over all targets and matching run counters concealed this quantification discrepancy.
