# Command-line guide

Start with the [README quick start](../README.md#as-a-binary) for installation, index
construction, paired/single-end quantification and standard 10x BUS output. This guide
covers additional options and protocols. Use `kallistors <command> --help` for the full
option list, and check the [compatibility guide](compatibility.md) for supported combinations.

## Quantification options

`quant -t N` sets the thread count for mapping and abundance inference. Paired-end
quantification estimates fragment lengths from the reads. Single-end quantification
requires the fragment mean (`-l`) and standard deviation (`-s`) in bases. The paired
histogram determines transcript-specific conditional means and effective lengths;
single-end uses a Gaussian model. See [fragment sampling](fragment_length_sampling.md).

Quantification writes `abundance.tsv`, `run_info.json` and, by default, `abundance.h5`
in the output directory. Field names and HDF5 structure follow kallisto conventions.
Use `--plaintext` to omit HDF5 output and write any bootstrap results as TSV files.

### Strand selection

Bulk quantification is unstranded by default. Choose the library's orientation with
`--fr-stranded` (first mate in transcript orientation) or `--rf-stranded` (first mate in
reverse-complement orientation). Paired libraries require the opposite orientation for
the second mate. For example, quantify a reverse-stranded paired library with:

```bash
kallistors quant -i transcripts.idx -o counts -t 8 --rf-stranded \
    reads_1.fq.gz reads_2.fq.gz
```

The same strand flags apply to single-end reads together with `--single -l MEAN -s SD`.
Choose only one strand flag; omit both for an unstranded library.

### Bootstrap samples

Use `-b/--bootstrap-samples` to request bootstrap samples and `--seed` to select the seed:

```bash
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
    -t 8 \
    -b 100 \
    --seed 42 \
    reads_1.fq reads_2.fq
```

HDF5 stores bootstrap count vectors under `/bootstrap/bs*`; `--plaintext` writes
`bs_abundance_*.tsv`. The HDF5 layout is compatible with original kallisto, but exact
bootstrap values can differ. See [output compatibility](compatibility.md#output).

### Sequence-bias correction

Enable sequence-specific bias correction with `--bias --transcripts transcripts.fa.gz`.
The transcript FASTA is required. The [completed bias comparison](kallisto-comparison-results.md#completed-sequence-bias-observations)
reports the observed differences from original kallisto.

## Mapping without abundance inference

`pseudoalign` writes counted equivalence classes (ECs): groups of transcripts compatible
with a read or paired fragment. It skips abundance inference and uses a TSV file as its
output, rather than a quantification directory:

```bash
kallistors pseudoalign --index transcripts.idx \
    --reads reads_1.fq.gz --reads2 reads_2.fq.gz \
    --threads 8 --out ec_counts.tsv
```

Omit `--reads2` for single-end mapping. The TSV columns are `ec_id`, `count` and
`transcripts`; the last column contains comma-separated, zero-based transcript IDs from
the index. EC numbers can change with processing order, so compare complete transcript
sets and their counts rather than row numbers. This is aggregated mapping output, not a
per-read alignment file.

For single-end EC comparisons with `quant`, also supply the same fragment mean using
`--fragment-length MEAN`, and match any strand or `--single-overhang` options. Without
that mean, `pseudoalign` does not apply the single-end fragment-length filter.

## Index options

`kallistors index` defaults to `k = 31` and one thread. When `-m/--min-size` is omitted,
it selects a kallisto-compatible minimizer length. Pass `-t N` to set the build thread
count and `--timings` to print timing information. Multiple FASTA files can be supplied
in input order. Repeated transcript names are rejected by default; `--make-unique`
appends numeric suffixes to make them distinct. Keep those renamed identifiers consistent
with downstream annotations.

Inspect an existing index's version, k-mer/minimizer lengths and transcript/graph counts with:

```bash
kallistors index-info --index transcripts.idx
```

This reads metadata; mapping compatibility requires a separate mapping check. The default
builder input is nucleotide transcript FASTA; `--aa` builds an amino-acid BUS index from
protein FASTA:

```bash
kallistors index \
    -i proteins.idx \
    --aa \
    proteins.fa.gz
```

## Single-cell protocols

The `-x` technology setting describes barcode, UMI and cDNA locations. List the supported
presets without supplying an index or reads:

```bash
kallistors bus --list
```

The README shows standard `10XV3` input. For a complete FASTQ-to-matrix workflow, see the
[PBMC benchmark recipe](benchmarks.md#single-cell-complete-10x-pbmc-1k-v3).

### Strand selection

Nucleotide BUS defaults to unstranded processing. Choose a strand flag from the library
protocol; the technology preset alone does not infer strandedness. For example, retain
10x v3 cDNA reads in transcript orientation with:

```bash
kallistors bus -i path/to/index.idx -o bus_fw -x 10XV3 --fr-stranded reads_1.fq reads_2.fq
```

### Amino-acid BUS mapping

Translate cDNA reads against the protein index built with `index --aa`:

```bash
kallistors bus \
    -i proteins.idx \
    -o bus_out \
    -x 10XV3 \
    --aa \
    reads_1.fq reads_2.fq
```

Amino-acid mode supports single-cDNA technologies and reports `n_frame_clashes` in
`run_info.json`. Paired BUS technologies are rejected in this mode.

### Transcriptome and genome BAM output

These output modes are supported by `bus`; `quant --pseudobam` and `quant --genomebam`
currently return an unsupported-feature error.

Write transcriptome pseudobam output for supported single-cDNA or paired technologies:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --pseudobam \
    reads_1.fq reads_2.fq
```

Project placements to genome coordinates and write a BAI or CSI index. Supply matching
GTF annotations and chromosome lengths:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o genome_bus_out \
    -x 10XV3 \
    --genomebam \
    --gtf annotations.gtf \
    --chromosomes chromosomes.txt \
    reads_1.fq reads_2.fq
```

Unpaired multi-sequence custom technologies are rejected for BAM output and strand-specific
processing. See [BUS compatibility](compatibility.md#bus).

### Custom technology and UMI tags

Custom technology triples can extract an RX UMI from FASTQ comments:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 0,0,16:RX:1,0,0 \
    reads_1.fq reads_2.fq
```

### Batch and interleaved input

A batch file lists a sample ID followed by the FASTQ files for that sample:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    -B samples.tsv \
    --batch-barcodes
```

Interleaved FASTQ input uses kallisto's `--inleaved` spelling:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --inleaved \
    interleaved.fq
```

### SmartSeq3

SmartSeq3 uses kallisto's default UMI tag sequence; override it with `-T` when needed:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x SMARTSEQ3 \
    -T ATTGCGCAATG \
    sample_bc1.fq sample_bc2.fq reads_1.fq reads_2.fq
```

### BAM input

Read barcode tags CR/CB and UMI tags UR/RX from BAM input:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV2 \
    --bam \
    alignments.bam
```

### Unmapped-read diagnostics and long reads

Write per-read unmapped k-mer ratios alongside BUS output:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --unmapped \
    reads_1.fq reads_2.fq
```

Long-read BUS mode filters reads above the unmapped k-mer threshold:

```bash
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --long \
    --platform ONT \
    --error-rate 0.01 \
    --threshold 0.8 \
    --unmapped \
    reads_1.fq reads_2.fq
```

General long-read quantification is not implemented; these options apply to BUS processing.

## Match-cache memory

Bulk `quant` and `pseudoalign` accept `--match-cache-size auto` (the default) or an
aggregate budget such as `256MiB` or `4GiB`. Automatic budgets scale with the number of
mapping workers. See [match-cache sizing](match-cache.md) for defaults, rounding and
memory accounting.
