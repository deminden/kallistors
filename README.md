# kallistors
kallistors: a Rust implementation of kallisto-style pseudoalignment and quantification.

## At A Glance
- Exact full-file paired-end `run_info.json` count parity with `kallisto` on the checked-in real
  dataset: `n_processed = 4408640`, `n_pseudoaligned = 4244771`, `n_unique = 276251`.
- Last measured full-file benchmark on `7950X3D`, `-t 32`:
  `kallisto 59.12s`, `kallistors 30.76s` (`1.92x` faster, about `48%` less wall time). See
  [bench_latest.md](bench_latest.md) for the benchmark record.
- Builds pure Rust kallisto v13-compatible transcriptome indexes for normal nucleotide FASTA.
  Last measured GENCODE build: `kallisto 6:24.81`, `kallistors 2:02.23` (`3.15x` faster,
  about `18%` lower peak RSS).
- Uses the `zlib-rs` gzip backend for faster compressed FASTQ input.
- Uses packed/reusable FASTQ batches plus long-lived worker `EcCounts` state to reduce allocation
  and threaded handoff overhead.

## Compatibility notes (v0.3.1)

This is a focused reimplementation at the current stage, not a drop-in replacement for `kallisto`.

| Area | Status |
| --- | --- |
| Existing `kallisto` index loading | Supported |
| Pure Rust index building | Supported for normal nucleotide transcript FASTA, v13-compatible output |
| Paired-end quant | Exact `run_info.json` count parity on the checked-in paired dataset |
| Single-end quant | Implemented, with synthetic and selected parity coverage |
| Sequence-specific bias correction | Optional with `--bias` and `--transcripts` |
| Bootstrap / H5 output | Supported for quant; `abundance.h5` is written by default, plaintext bootstrap TSVs with `--plaintext` |
| Long-read, UMI/BUS/technology modes, fusion detection | Not implemented |
| CLI option coverage | Partial |

Current real-data status:
- Paired-end parity is exact on deterministic prefixes through `262144` pairs from the checked-in
  real dataset under `data/` with `--threads 32`.
- Full-file paired-end quant has exact `run_info.json` count parity on the same checked-in index and
  reads.
- Full-file runtime is currently faster than local upstream `kallisto`; current timings and stage
  measurements are in [bench_latest.md](bench_latest.md).
- Pure Rust index output is accepted by both `kallistors` and upstream `kallisto inspect/quant` on
  the reduced real subset and the GENCODE validation path. Generated indexes are format-compatible,
  not byte-identical to upstream indexes.
- `abundance.tsv` is validated with floating-point tolerances where tests compare estimates. The
  README does not claim bit-for-bit abundance parity.
- H5 output is readable by upstream `kallisto h5dump` and uses kallisto-compatible default bias
  datasets when sequence bias is disabled. Bootstrap dataset layout is compatible, but bootstrap
  sample values are not yet expected to be byte-identical to upstream kallisto because the
  RNG/sampling path still differs.

Sequence-specific bias correction is optional and enabled only with `--bias`.
The index builder currently targets nucleotide transcript FASTA. It intentionally rejects kallisto
features outside that supported builder surface, including amino-acid mode, distinguish mode,
d-list-specific behavior, and technology-specific modes.

## Usage

### As a Binary

```bash
# Install the CLI from git
cargo install --git https://github.com/deminden/kallistors kallistors

# Build a kallisto-compatible transcriptome index
kallistors index \
    -i transcripts.idx \
    -t 8 \
    --timings \
    transcripts.fa.gz

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

# Quantify with plaintext bootstrap outputs instead of H5
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
    -b 100 \
    --plaintext \
    reads_1.fq reads_2.fq

# Quantify (single-end)
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
    --single \
    -l 200 \
    -s 20 \
    reads.fq.gz

# Pseudoalign (single-end; emits EC counts)
kallistors pseudoalign \
    --index path/to/index.idx \
    --reads reads.fq.gz \
    --fragment-length 200 \
    --out ec_counts.tsv
```

Notes:
- Paired-end quant estimates fragment length mean/sd from pseudoaligned pairs.
- Quant writes `abundance.tsv`, `run_info.json`, and, unless `--plaintext` is set, `abundance.h5`
  in `out_dir` (matching kallisto field names and upstream-readable HDF5 structure).
- Bootstrap samples use `-b/--bootstrap-samples`; H5 stores bootstrap count vectors under
  `/bootstrap/bs*`, while `--plaintext` writes `bs_abundance_*.tsv`. Bootstrap H5 layout is
  compatible, but exact upstream bootstrap sample values are not currently claimed.
- `--bias` requires `--transcripts` to provide the transcript FASTA.
- `kallistors index` defaults to `k = 31`, `threads = 1`, and kallisto-compatible minimizer length
  selection when `-m/--min-size` is omitted.

### Recent Compatibility And Performance Work
- Added a pure Rust `kallistors index` builder and `kallistors::index::build_index` library entry
  point for kallisto v13-compatible nucleotide transcript indexes.
- Added quant bootstrap sampling with kallisto-style H5 output, upstream-readable HDF5 metadata,
  kallisto-compatible default bias datasets, and optional plaintext `bs_abundance_*.tsv` files.
- The builder reads plain or gzipped FASTA, normalizes sequences with kallisto-compatible ambiguous
  base replacement and poly-A clipping, rejects duplicate transcript names unless `--make-unique`
  is used, and preserves original transcript lengths in index metadata.
- Index building now keeps normalized transcripts in memory, builds compacted dBG/EC payloads
  directly, writes BooPHF-compatible minimizer indexes, and reports per-stage timings with
  `--timings`.
- The latest GENCODE index benchmark is faster than local upstream `kallisto index` while preserving
  upstream load and quant parity. Details are in [bench_latest.md](bench_latest.md).
- Full-file paired-end `run_info.json` count parity now matches `kallisto` exactly on the checked-in
  real dataset.
- The last full-file mismatch was fixed with a narrow Bifrost-style retry on probe/backoff misses
  after prior evidence exists.
- The loader now uses a much cheaper minimizer count/fill path, which removed most of the old
  startup penalty.
- `flate2` now uses the `zlib-rs` backend.
- The threaded path now transports reads as packed/reusable FASTQ batches with one contiguous
  backing buffer plus per-record offsets instead of allocating owned FASTQ payloads per read.
- Threaded workers accumulate directly into long-lived `EcCounts` state instead of rebuilding and
  merging fresh per-batch count maps.
- The common fast pseudoalignment path reuses encoded k-mer codes through minimizer candidate
  lookup, match-cache keying, and jump/middle/scan probes.
- Hot pseudoalignment env flags are cached once per process instead of being read from the
  environment in per-read/per-k-mer paths.
- EM now pre-splits singleton/nonzero multi-transcript ECs and stores multi-EC transcript and weight
  metadata in contiguous arrays.
- Small direct-mapped hot lookup caches for MPHF minimizer lookup and EC block lookup were enlarged
  to reduce collision misses in the common path.
- Quant now reuses transcript metadata already loaded by the Bifrost index path instead of reloading
  index metadata and cloning EC tables before EM.
- The Bifrost loader skips redundant graph pre-scanning on the supported `k <= 32` path, lazily
  allocates shade metadata, reuses graph-node buffers, and avoids paired positional payload loading
  when paired fragment estimation only needs flat EC block bounds.

### Real-data benchmark (current checked-in dataset)
The detailed benchmark record lives in [bench_latest.md](bench_latest.md). The commands below are
the benchmark inputs used for the last measured result.

Dataset + index (in `data/` in this repo):
- Paired FASTQs:
  `data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz`
  `data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz`
- Reference transcripts: `data/gencode.v49.transcripts.fa.gz`
- Index: `data/gencode.v49_kallisto.idx`

```bash
# Run kallisto
kallisto_src/build/src/kallisto quant \
  -i data/gencode.v49_kallisto.idx -o /tmp/kallisto_full -t 32 \
  data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  "data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz"

# Run kallistors
./target/release/kallistors quant \
  -i data/gencode.v49_kallisto.idx -o /tmp/kallistors_full -t 32 \
  data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  "data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz"
```

Last measured result: `kallisto 59.12s`, `kallistors 30.76s` on 2026-05-17, Linux x86_64,
`7950X3D`, `-t 32`, with exact `run_info.json` counts. Stage timings in
[bench_latest.md](bench_latest.md) are instrumentation timings for individual phases and should
not be treated as additive wall-clock components.

### Index-build benchmark
The latest large-transcriptome index benchmark used `data/gencode.v49.transcripts.fa.gz`:

```bash
/usr/bin/time -o /tmp/kallisto-index.time -f 'elapsed=%E maxrss=%MKB' \
  kallisto_src/build/src/kallisto index \
  -i /tmp/gencode.kallisto.idx \
  data/gencode.v49.transcripts.fa.gz

/usr/bin/time -o /tmp/kallistors-index.time -f 'elapsed=%E maxrss=%MKB' \
  ./target/release/kallistors index \
  -i /tmp/gencode.kallistors.idx \
  -t 8 \
  --timings \
  data/gencode.v49.transcripts.fa.gz
```

Last measured result: `kallisto 6:24.81`, peak RSS `12175612KB`, index size `878M`;
`kallistors 2:02.23`, peak RSS `9975360KB`, index size `894M`. The kallistors-built GENCODE index
loads in upstream `kallisto inspect`, loads in `kallistors index-info`, and gives exact synthetic
paired-read count parity through both `kallistors quant` and upstream `kallisto quant`.
Use `scripts/bench_index_build.py` to regenerate a Markdown report for another machine or FASTA.

### Parity tests
- Synthetic parity: `variants_parity` requires exact aligned-count parity with `kallisto`.
- Real paired-prefix parity requires exact `n_processed`, `n_pseudoaligned`, and `n_unique` counts at
  every checked prefix.
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
