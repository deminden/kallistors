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

### Kallisto alignment parity (recent changes)
- Bifrost path mirrors kallisto minimizer traversal with minhash candidates and canonicalized minimizers.
- Orientation-aware minimizer offsets (forward/rev) and relaxed window checks around minimizer positions.
- Online intersection tracking with `intersection_empty` skip parity.
- Jump logic records a hit before skipping to match kallisto’s jump behavior.
- Special/overcrowded minimizers handled via D-list / special unitig checks.

### Debugging tools
- `trace-reads` supports EC traces, per-hit dumps, intersection dumps, and positions visited.
- Minimizer debug includes MPH lookup info, unitig match decisions, and special/D-list flags.
- Debug helpers: `minimizer-lookup`, `minimizer-bitmap-scan`, `minimizer-unitig-mphf-check`,
  and `minimizer-mphf-keys` for inspecting MPHF/index consistency.

Debug flags (CLI)
- Pseudoalign/quant parity switches: `--kallisto-enum`, `--kallisto-strict`, `--kallisto-local-fallback`,
  `--kallisto-fallback`, `--kallisto-direct-kmer`, `--kallisto-bifrost-find`,
  `--discard-special-only`, `--skip-overcrowded-minimizer`, `--no-jump`, `--union`,
  `--min-range`, `--dfk-onlist`.
- Trace reads: `--hits-out`, `--hits-intersection-out`, `--positions-visited-out`,
  `--minimizer-positions`, `--minimizer-out`, `--local-kmer-out`, `--gene-map`,
  plus the parity switches above and strand/fragment options (`--strand`, `--fragment-length`,
  `--single-overhang`, `--fr-stranded`, `--rf-stranded`).

### Real-data benchmark (latest, no debug)
Dataset + index (in `data/` in this repo):
- Reads trimmed with fastp and 1/100 selected: `data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz`
- Reference transcripts: `data/Homo_sapiens.GRCh38.cdna.all.fa.gz`
- Index (rebuilt): `data/kallisto_index`
- Mean read length: 168.90 (used `-l 169 -s 20`)

```bash

# Run kallisto (no debug)
kallisto quant \
  -i data/kallisto_index -o data/kallisto_bench_run \
  --single -l 169 -s 20 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz

# Run kallistors (no debug)
./target/release/kallistors-cli quant \
  -i data/kallisto_index -o data/kallistors_bench_run \
  --single -l 169 -s 20 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
```

Latest results (2025-12-29, macOS arm64, no debug):
- Speed: kallisto real 19.26s; kallistors real 36.38s
- n_pseudoaligned: kallisto 41626 / 43817; kallistors 41666 / 43817
- TPM Pearson 0.974175, TPM MAE 11.206 (TPM filter max(k_tpm, o_tpm) > 1; n = 7824)
- est_counts Pearson 0.998049, est_counts MAE 0.303

### Parity tests
- Synthetic parity: `variants_parity` allows a 1-read drift in aligned count; EC sets must match when aligned.
- Other parity tests keep strict equality on aligned counts and EC sets.

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



### Debugging abundance diffs

Compare `abundance.tsv` outputs and rank transcripts by relative error:
```bash
scripts/compare_abundance.py \
  --kallisto /path/to/kallisto/abundance.tsv \
  --kallistors /path/to/kallistors/abundance.tsv \
  --top 50
```

Trace per-read EC decisions for a small set of reads (single-end):
```bash
printf "READ_ID_1\nREAD_ID_2\n" > read_list.txt
target/release/kallistors-cli trace-reads \
  --index data/Human_kallisto_index \
  --reads data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz \
  --read-list read_list.txt \
  --out read_traces.tsv \
  --fragment-length 200
```

## Contributing

Contributions are very welcome! 
If you’d like to help improve `kallistors`, feel free to open an issue to discuss ideas, report bugs, or request features.

Pull requests are encouraged, especially for:
- performance improvements
- correctness / numerical stability fixes
- additional tests (including cross-validation vs original)
- documentation, examples, and benchmarking

### Development notes

- Please run formatting and linting before submitting:
  ```bash
  cargo fmt --all
  cargo clippy --workspace --all-targets --all-features -- -D warnings
  cargo test --workspace --all-features
  ```

## License

This project is licensed under the BSD-2-Clause (Simplified BSD) License. See `LICENSE`.
