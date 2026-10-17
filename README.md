# kallistors

A Rust implementation of kallisto-style pseudoalignment and RNA-seq quantification.

## At a glance

- Builds and loads kallisto v13-compatible transcriptome indexes.
- Quantifies paired-end and single-end reads, with kallisto-style `abundance.tsv`, `run_info.json`,
  and default `abundance.h5` output.
- Supports single-cell BUS workflows, including supported BAM input and output modes.
- Uses compact index/EM storage, cached exact matching and parallel quantification.
- Paired fragment estimation follows kallisto's sequence-buffer sampling rules; all reads
  contribute to EC counts. See [sampling behavior](docs/fragment_length_sampling.md).

## Performance and benchmarking

Quantification uses the requested thread count for mapping and EM. Compact index storage,
packed matching, bounded input queues and buffered output reduce work and memory use.
Source builds use `target-cpu=native`; release binaries use a portable CPU target.

The [benchmark guide](docs/benchmarks.md) documents public datasets, read trimming,
GENCODE cDNA and Pachter Lab mature-plus-nascent references, repeated comparisons with original
kallisto, phase timing, memory measurement and complete EC validation. Performance depends on the
input, reference, thread count, build and hardware; compare complete quantification runs.

Current release details are in the [0.5.1 notes](docs/releases/0.5.1.md) and
[changelog](CHANGELOG.md).

## Compatibility notes

The notes below describe the supported interface. See [compatibility limits](docs/compatibility.md)
for differences from original kallisto.

| Area | Status |
| --- | --- |
| Existing `kallisto` index loading | Supported |
| Pure Rust index building | Supported for normal nucleotide transcript FASTA, v13-compatible output; amino-acid BUS indexes are supported with `index --aa` |
| Paired-end quant | Uses an empirical fragment histogram for effective lengths |
| Single-end quant | Uses the supplied Gaussian fragment model |
| Single-cell BUS output | Fixed presets and custom `-x` triples; batch/interleaved FASTQ, BAM input tags, SmartSeq3 tags, supported amino-acid and long-read BUS modes, pseudobam, and projected/sorted genome BAM with BAI or CSI index |
| Sequence-specific bias correction | Optional with `--bias` and `--transcripts` |
| Bootstrap / H5 output | Supported for quant; `abundance.h5` is written by default, plaintext bootstrap TSVs with `--plaintext` |
| Long-read and fusion detection | Long-read BUS filtering/output is implemented; general long-read quant and fusion detection are not implemented |
| CLI option coverage | Partial |

Current validation and output compatibility:

- Complete EC comparisons use counted transcript sets, independent of EC numbering.
- Regression coverage includes paired/single-end matching, TTN/AGAP12P fallback probes,
  abundant-node associations, empirical fragment models, and serial/threaded EM.
  Full-index regression tests are explicitly enabled when the downloaded reference is present.
- Paired fragment-length estimation follows original's 10,000-observation goal and 8 MiB
  sequence-buffer quotas, including possible concurrent overshoot. Every read contributes
  to EC counting. See [fragment sampling](docs/fragment_length_sampling.md).
- Pure Rust indexes are checked against standard `MAX_KMER_SIZE=32` upstream builds.
  Format compatibility does not require byte-identical index files.
- Native HDF5 readback checks target IDs, estimated counts, effective lengths, target lengths
  and run metadata against TSV. Bootstrap output is supported, but its RNG/sample values
  are not expected to match upstream exactly.
- Exact abundance equality and real-data sequence-bias parity are not general claims.

Sequence-specific bias correction is optional and enabled only with `--bias`.
The index builder targets nucleotide transcript FASTA by default. Use `index --aa` for BUS indexes
built from protein FASTA; unsupported builder features still include distinguish mode and
d-list-specific behavior.

## Usage

### As a Binary

Building or installing from source requires Rust 1.97 or newer. The project uses Rust edition 2024.

```bash
# Download prebuilt binaries from GitHub Releases, or install the CLI from source.
#
# Install the CLI from git
cargo install --locked --git https://github.com/deminden/kallistors kallistors

# Build a kallisto-compatible transcriptome index
kallistors index \
    -i transcripts.idx \
    -t 8 \
    --timings \
    transcripts.fa.gz

# Build an amino-acid BUS index from protein FASTA
kallistors index \
    -i proteins.idx \
    --aa \
    proteins.fa.gz

# Quantify (paired-end)
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
    reads_1.fq reads_2.fq

# Quantify with bootstrap samples written to abundance.h5, the default binary output
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
    -b 100 \
    --seed 42 \
    reads_1.fq reads_2.fq

# Generate BUS output for a 10x Genomics v3-style single-cell run
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    reads_1.fq reads_2.fq

# Generate BUS output by translating cDNA reads against an amino-acid index
kallistors bus \
    -i proteins.idx \
    -o bus_out \
    -x 10XV3 \
    --aa \
    reads_1.fq reads_2.fq

# Write transcriptome pseudobam output for supported single-cDNA or paired BUS technologies
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --pseudobam \
    reads_1.fq reads_2.fq

# Project pseudobam placements to genome coordinates and write a BAI or CSI index
kallistors bus \
    -i path/to/index.idx \
    -o genome_bus_out \
    -x 10XV3 \
    --genomebam \
    --gtf annotations.gtf \
    --chromosomes chromosomes.txt \
    reads_1.fq reads_2.fq

# Custom BUS technology triples support RX UMI extraction from FASTQ comments
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 0,0,16:RX:1,0,0 \
    reads_1.fq reads_2.fq

# Generate BUS output from a batch file: sample_id followed by the FASTQ files for that sample
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    -B samples.tsv \
    --batch-barcodes

# Interleaved FASTQ input is accepted with kallisto's --inleaved spelling
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --inleaved \
    interleaved.fq

# SmartSeq3 uses kallisto's default UMI tag sequence; override it with -T if needed
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x SMARTSEQ3 \
    -T ATTGCGCAATG \
    sample_bc1.fq sample_bc2.fq reads_1.fq reads_2.fq

# Read a BAM file with CR/CB barcode tags and UR/RX UMI tags
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV2 \
    --bam \
    alignments.bam

# Write per-read unmapped k-mer ratios alongside BUS output
kallistors bus \
    -i path/to/index.idx \
    -o bus_out \
    -x 10XV3 \
    --unmapped \
    reads_1.fq reads_2.fq

# Long-read BUS mode filters reads above the unmapped k-mer threshold
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

# Quantify (single-end)
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
    --single \
    -l 200 \
    -s 20 \
    reads.fq.gz
```

Notes:

- Versioned [GitHub Releases](https://github.com/deminden/kallistors/releases) attach native release archives for Linux x86_64, macOS x86_64,
  macOS arm64, and Windows x86_64, plus SHA-256 checksum files.
- Paired-end quant uses the observed fragment histogram for transcript-specific conditional means
  and effective lengths. Single-end uses `-l` and `-s` for its Gaussian fragment model.
- Quant writes `abundance.tsv`, `run_info.json`, and, unless `--plaintext` is set, `abundance.h5`
  in `out_dir` (matching kallisto field names and upstream-readable HDF5 structure).
- Bootstrap samples use `-b/--bootstrap-samples`; H5 stores bootstrap count vectors under
  `/bootstrap/bs*`, while `--plaintext` writes `bs_abundance_*.tsv`. Bootstrap H5 layout is
  compatible, but exact upstream bootstrap sample values are not currently claimed.
- `--bias` requires `--transcripts` to provide the transcript FASTA.
- `kallistors index` defaults to `k = 31`, `threads = 1`, and kallisto-compatible minimizer length
  selection when `-m/--min-size` is omitted.
- `bus --aa` currently supports single-cDNA BUS technologies and reports `n_frame_clashes` in
  `run_info.json`; paired BUS technologies are rejected in amino-acid mode.
- BUS pseudobam/genomebam output is supported for single-cDNA and paired BUS technologies. Unpaired
  multi-sequence custom technologies are rejected for BAM output and strand-specific processing
  until their placement semantics are made explicit.

### Performance and validation

- Data preparation, benchmark commands and measurement methodology are in
  [docs/benchmarks.md](docs/benchmarks.md).
- Developer-facing trace workflows, abundance diffs, and parity-debug commands are in
  [docs/development.md](docs/development.md).
- Instrumented upstream kallisto build notes are in
  [docs/kallisto_debug_build.md](docs/kallisto_debug_build.md).
- Reproducible downloads, MD5/gzip checks, fasterp trimming, and repeated comparisons are documented
  in [the benchmark guide](docs/benchmarks.md). Install fasterp separately with Cargo.
  Quantification benchmarks require complete EC parity before timing and export JSON/TSV reports.

### As a Crate

The serial quantification APIs use one thread. Use `em_quantify_with_fragment_histogram_threaded`
with an explicit thread count for parallel EM. The CLI `quant -t N` uses N threads for mapping
and EM in successive phases. HDF5 output uses bounded chunks. Set the
`PseudoalignOptions::fragment_length_samples` field when using exhaustive struct literals, or use
`..PseudoalignOptions::default()` for the 10,000-observation goal. Use zero to collect all
qualifying fragment observations.

Add to `Cargo.toml`:
```toml
[dependencies]
kallistors = { git = "https://github.com/deminden/kallistors" }
```

Build an index from Rust:
```rust
use std::error::Error;
use std::path::{Path, PathBuf};

use kallistors::index::{build_index, IndexBuildOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let fasta = vec![PathBuf::from("transcripts.fa.gz")];
    build_index(
        Path::new("transcripts.idx"),
        &fasta,
        IndexBuildOptions::default(),
    )?;
    Ok(())
}
```

Minimal library usage sketch:
```rust
use std::error::Error;
use std::path::Path;

use kallistors::index::Index;
use kallistors::io::open_fastq_reader;
use kallistors::pseudoalign::{
    build_bifrost_index, pseudoalign_paired_bifrost_with_options, PseudoalignOptions, Strand,
};
use kallistors::quant::{em_quantify_with_fragment_histogram, EcCountsInput, QuantOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let index_path = Path::new("path/to/index.idx");
    let mut reads_1 = open_fastq_reader(Path::new("reads_1.fq.gz"))?;
    let mut reads_2 = open_fastq_reader(Path::new("reads_2.fq.gz"))?;
    let index = build_bifrost_index(index_path)?;
    let ec = pseudoalign_paired_bifrost_with_options(
        &index,
        &mut reads_1,
        &mut reads_2,
        Strand::Unstranded,
        PseudoalignOptions::default(),
    )?;

    let meta = Index::load(index_path)?;
    let input = EcCountsInput {
        ec_list: kallistors::ec::EcList {
            classes: ec.ec_list,
        },
        counts: ec.counts,
    };
    let lengths: Vec<u32> = meta.transcripts.iter().map(|t| t.length).collect();
    let quant = em_quantify_with_fragment_histogram(
        &input,
        &lengths,
        None,
        None,
        ec.fragment_length_hist.as_deref(),
        QuantOptions::default(),
    )?;

    println!("reads processed: {}", ec.reads_processed);
    println!("reads aligned: {}", ec.reads_aligned);
    println!("targets: {}", quant.est_counts.len());
    Ok(())
}
```

### Development docs

Developer-facing checks, trace workflows, abundance-diff commands, and real-data
parity workflows are in [docs/development.md](docs/development.md). The
[source layout](docs/development.md#source-layout) describes mapping, quantification,
diagnostics and BUS output modules.

## Contributing

Open an issue to discuss ideas, report bugs, or request features.

Pull requests are encouraged, especially for:

- performance improvements
- correctness / numerical stability fixes
- additional tests (including cross-validation vs original)
- documentation, examples, and benchmarking

Development workflow notes are in [docs/development.md](docs/development.md).

## License

This project is licensed under the MIT License. See `LICENSE`.
