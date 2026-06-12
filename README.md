# kallistors
kallistors: a Rust implementation of kallisto-style pseudoalignment and quantification.

## At A Glance
- Paired-end quant is faster than local upstream `kallisto` on the checked-in real dataset. Latest
  median-of-5 full-file benchmark on `7950X3D`, `-t 32`: `kallisto 55.81s`, `kallistors 32.79s`
  (`1.70x` faster). See [docs/benchmarks.md](docs/benchmarks.md).
- Builds pure Rust kallisto v13-compatible transcriptome indexes for normal nucleotide FASTA.
  Last measured GENCODE build: `kallisto 6:24.81`, `kallistors 2:02.23` (`3.15x` faster,
  about `18%` lower peak RSS).
- Quant writes kallisto-style `abundance.tsv`, `run_info.json`, and default `abundance.h5`.

## Compatibility notes (v0.4.0)

This is a focused reimplementation at the current stage, not a drop-in replacement for `kallisto`.

| Area | Status |
| --- | --- |
| Existing `kallisto` index loading | Supported |
| Pure Rust index building | Supported for normal nucleotide transcript FASTA, v13-compatible output; amino-acid BUS indexes are supported with `index --aa` |
| Paired-end quant | Implemented; deterministic prefix and latest full-file checks have exact `run_info.json` count parity |
| Single-end quant | Implemented, with synthetic and selected parity coverage |
| Single-cell BUS output | Implemented for fixed technology presets and kallisto-style custom `-x` triples, including batch/interleaved input, batch barcode sidecars, BAM input tags and missing-tag skips, SmartSeq3 tags, `--num`, `--union`, `--no-jump`, `--unmapped` ratio output, amino-acid `--aa` mode for single-cDNA technologies, first-pass `--long` threshold filtering with `--platform`/`--error-rate` and valid `novel.fastq`, long/paired BUS sidecars, pseudobam output, and projected/sorted genome BAM output with BAI index |
| Sequence-specific bias correction | Optional with `--bias` and `--transcripts` |
| Bootstrap / H5 output | Supported for quant; `abundance.h5` is written by default, plaintext bootstrap TSVs with `--plaintext` |
| Long-read and fusion detection | Long-read BUS filtering/output is implemented; general long-read quant and fusion detection are not implemented |
| CLI option coverage | Partial |

Current real-data status:
- Paired-end parity is exact on the deterministic `1,048,576`-pair subset across `1, 2, 4, 8, 16,
  32` threads.
- The active working tree's latest full-file paired check matches kallisto counts exactly:
  `n_processed=4,408,640`, `n_pseudoaligned=4,244,771`, `n_unique=276,251`.
- Full-file runtime is currently faster than local upstream `kallisto`.
- Pure Rust index output is accepted by both `kallistors` and upstream `kallisto inspect/quant` on
  the reduced real subset and the GENCODE validation path. Generated indexes are format-compatible,
  not byte-identical to upstream indexes.
- `abundance.tsv` is validated with floating-point tolerances where tests compare estimates. The
  README does not claim bit-for-bit abundance parity.
- H5 output is readable by upstream `kallisto h5dump` and uses kallisto-compatible default bias
  datasets when sequence bias is disabled. H5 writing uses a pure Rust writer, so building
  `kallistors` does not require system HDF5 development headers/libraries. Bootstrap dataset
  layout is compatible, but bootstrap sample values are not yet expected to be byte-identical to
  upstream kallisto because the RNG/sampling path still differs.

Sequence-specific bias correction is optional and enabled only with `--bias`.
The index builder targets nucleotide transcript FASTA by default. Use `index --aa` for BUS indexes
built from protein FASTA; unsupported builder features still include distinguish mode and
d-list-specific behavior.

## Usage

### As a Binary

```bash
# Download prebuilt binaries from GitHub Releases, or install the CLI from source.
#
# Install the CLI from git
cargo install --git https://github.com/deminden/kallistors kallistors

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

# Project pseudobam placements to genome coordinates and write pseudoalignments.bam.bai
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
- Versioned GitHub Releases attach native release archives for Linux x86_64, macOS x86_64,
  macOS arm64, and Windows x86_64, plus SHA-256 checksum files.
- Paired-end quant estimates fragment length mean/sd from pseudoaligned pairs.
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
- Current benchmark tables, raw artifact paths, and methodology are in
  [docs/benchmarks.md](docs/benchmarks.md).
- Developer-facing trace workflows, abundance diffs, and parity-debug commands are in
  [docs/development.md](docs/development.md).
- Instrumented upstream kallisto build notes are in
  [docs/kallisto_debug_build.md](docs/kallisto_debug_build.md).

### Parity expectations
- Synthetic parity: `variants_parity` requires exact aligned-count parity with `kallisto`.
- Real paired-prefix parity requires exact `n_processed`, `n_pseudoaligned`, and `n_unique` counts at
  every checked prefix.
- The isolated full-index paired regression test is opt-in so CI is not forced to load the large
  local GENCODE index; run it with `KALLISTORS_RUN_REAL_PAIRED_REGRESSION=1 cargo test -p kallistors --test real_paired_regression`.
- Tests that compare read-level ECs require exact EC equality for aligned reads.
- Abundance and percentage comparisons use floating-point tolerances where appropriate; exact
  bit-for-bit `abundance.tsv` parity is not currently claimed.

### As a Crate

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
use kallistors::quant::{em_quantify, EcCountsInput, QuantOptions};

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
    let quant = em_quantify(&input, &lengths, None, None, QuantOptions::default())?;

    println!("reads processed: {}", ec.reads_processed);
    println!("reads aligned: {}", ec.reads_aligned);
    println!("targets: {}", quant.est_counts.len());
    Ok(())
}
```

### Development docs

Developer-facing checks, trace workflows, abundance-diff commands, and real-data
parity workflows are in [docs/development.md](docs/development.md).

## Contributing

Contributions are very welcome! 
If you’d like to help improve `kallistors`, feel free to open an issue to discuss ideas, report bugs, or request features.

Pull requests are encouraged, especially for:
- performance improvements
- correctness / numerical stability fixes
- additional tests (including cross-validation vs original)
- documentation, examples, and benchmarking

Development workflow notes are in [docs/development.md](docs/development.md).

## License

This project is licensed under the MIT License. See `LICENSE`.
