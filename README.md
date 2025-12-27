# kallistors

kallistors: a Rust implementation of kallisto-style pseudoalignment and quantification.

## Compatibility notes (current gaps)

This is a minimal reimplementation (at current stage). Known differences vs kallisto today:
- No index builder yet (load-only).
- No bootstrap or H5 output.
- No long-read, UMI/BUS/technology modes, or fusion detection.
- Minimal CLI surface; some kallisto options are missing.
- Not a drop-in replacement yet; use for experimentation and validation, not as a production substitute.

Sequence-specific bias correction is optional and enabled only with `--bias`.

## Usage

### As a Binary

```bash
# Build
git clone https://github.com/deminden/kallistors
cd kallistors
cargo build --workspace --release

# Convenience wrapper (similar to kallisto-style usage)
kallistors() { ./target/release/kallistors-cli "$@"; }

# Quantify (paired-end)
kallistors quant \
    -i path/to/index.idx \
    -o out_dir \
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


Notes:
- Paired-end quant estimates fragment length mean/sd from pseudoaligned pairs.
- Quant writes `abundance.tsv` and `run_info.json` in `out_dir` (matching kallisto field names).
- `--bias` requires `--transcripts` to provide the transcript FASTA.
```

### As a Crate

Add to `Cargo.toml`:
```toml
[dependencies]
kallistors = { git = "https://github.com/deminden/kallistors" }
```

Use in code (simpler, CLI-like flow):
```rust
use kallistors::index::Index;
use kallistors::io::open_fastq_reader;
use kallistors::pseudoalign::{
    build_bifrost_index, PseudoalignOptions, Strand, pseudoalign_paired_bifrost_with_options,
};
use kallistors::quant::{em_quantify, QuantOptions};

let index_path = "path/to/index.idx";
let index = build_bifrost_index(index_path)?;
let mut r1 = open_fastq_reader("reads_1.fq.gz".as_ref())?;
let mut r2 = open_fastq_reader("reads_2.fq.gz".as_ref())?;

let ec_counts = pseudoalign_paired_bifrost_with_options(
    &index,
    &mut r1,
    &mut r2,
    Strand::Unstranded,
    PseudoalignOptions::default(),
)?;

let meta = Index::load(index_path)?;
let lengths = meta.transcripts.iter().map(|t| t.length).collect::<Vec<_>>();

let result = em_quantify(
    &kallistors::quant::EcCountsInput {
        ec_list: kallistors::ec::EcList {
            classes: ec_counts.ec_list,
        },
        counts: ec_counts.counts,
    },
    &lengths,
    None,
    None,
    QuantOptions::default(),
)?;
```

Use in code (explicit EC counts, lower-level flow):
```rust
use kallistors::index::Index;
use kallistors::quant::{em_quantify, load_ec_counts, QuantOptions};

let index = Index::load("path/to/index.idx")?;
let ec_input = load_ec_counts("path/to/ec_counts.tsv")?;
let lengths = index.transcripts.iter().map(|t| t.length).collect::<Vec<_>>();

let result = em_quantify(
    &ec_input,
    &lengths,
    None,
    None,
    QuantOptions::default(),
)?;
```



### Real-data parity and benchmarks

Dataset + index (already in `data/` in this repo):
- Reads: `data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed.fastq.gz`
- Reference transcripts: `data/Homo_sapiens.GRCh38.cdna.all.fa.gz`
- Index: `data/Human_kallisto_index` (built from the reference transcript FASTA)

Generate a ~1/100 subset (single-end; deterministic sampling):
```bash
scripts/sample_fastq.py \
  --input1 data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed.fastq.gz \
  --output1 data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz \
  --fraction 0.01 \
  --seed 42
```

Run the real-data parity test (skips if files or `kallisto` are missing):
```bash
cargo test -p kallistors-cli --test real_kallisto_parity -- --nocapture
```

Run benchmarks vs kallisto and write a report:
```bash
scripts/bench_real_kallisto.py \
  --index data/Human_kallisto_index \
  --reads data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz \
  --threads 8
```

Latest benchmark output is stored in `bench_latest.md`.
Latest benchmark summary (2025-12-27, macOS arm64, 8 threads, ~43k reads):
- kallisto wall 11.106s / cpu 30.802s; kallistors wall 36.013s / cpu 36.408s
- run_info: n_pseudoaligned 41444 vs 41113; n_unique 2420 vs 2111
- abundance parity: p99 rel error tpm 0.00882; est_counts 1.43e-07

## Contributing

Contributions are very welcome! 
If you’d like to help improve `kallistors`, feel free to open an issue to discuss ideas, report bugs, or request features.

Pull requests are encouraged — especially for:
- performance improvements
- correctness / numerical stability fixes
- additional tests (including cross-validation vs original)
- documentation, examples, and benchmarking

### Development notes

- Please run formatting and linting before submitting:
  ```bash
  cargo fmt --all
  cargo clippy --workspace --all-targets --all-features
  cargo test --workspace --all-features
  ```

## License

This project is licensed under the BSD-2-Clause (Simplified BSD) License. See `LICENSE`.
