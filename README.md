# kallistors

kallistors: a Rust implementation of kallisto-style pseudoalignment and quantification.

## Compatibility notes (current gaps)

This is a minimal reimplementation (at current stage). Known differences vs kallisto today:
- No index builder yet (load-only).
- No bootstrap, H5 output, or full `run_info.json` parity.
- No long-read, UMI/BUS/technology modes, or fusion detection.
- No gzip FASTQ support yet.
- Minimal CLI surface; some kallisto options are missing.
- No multi-threaded pseudoalign/quant paths yet.

Sequence-specific bias correction is optional and enabled only with `--bias`.

## Usage

### As a Binary

```bash
# Build
git clone https://github.com/deminden/kallistors
cd kallistors
cargo build --workspace --release


# Pseudoalign (single-end)
./target/release/kallistors-cli pseudoalign \
    --index path/to/index.idx \
    --reads reads.fq \
    --fragment-length 200 \
    --out ec_counts.tsv

# Inspect index
./target/release/kallistors-cli index-info --index path/to/index.idx

# Quantify
./target/release/kallistors-cli quant \
    --index path/to/index.idx \
    --ec ec_counts.tsv \
    --out abund.tsv \
    --fragment-length 200 \
    --fragment-length-sd 20
```

### As a Crate

Add to `Cargo.toml`:
```toml
[dependencies]
kallistors = { git = "https://github.com/deminden/kallistors" }
```

Use in code:
```rust
use kallistors::index::Index;
use kallistors::quant::{em_quantify, load_ec_counts, QuantOptions};

let index = Index::load("path/to/index.idx")?;
let ec_input = load_ec_counts("ec_counts.tsv")?;
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
