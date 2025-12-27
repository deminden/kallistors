# kallistors

kallistors: a Rust implementation of kallisto-style pseudoalignment and quantification.

## Compatibility notes (current gaps)

This is a minimal reimplementation (at current stage). Known differences vs kallisto today:
- No index builder yet (load-only).
- No bootstrap or H5 output.
- No long-read, UMI/BUS/technology modes, or fusion detection.
- Minimal CLI surface; some kallisto options are missing.

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

## Development

Recommended checks:
- cargo fmt --check
- cargo clippy --all-targets -- -D warnings
- cargo test --workspace

Regen note: synthetic datasets and local outputs can be regenerated via `scripts/generate_simple_dataset.py` and `scripts/compare_kallisto_simple.py`.

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
