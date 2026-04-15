# kallistors

kallistors: a Rust implementation of kallisto-style pseudoalignment and quantification.

## Compatibility notes (v0.2.1)

This is a minimal reimplementation (at current stage). Known differences vs kallisto today:
- No index builder yet (load-only).
- No bootstrap or H5 output.
- No long-read, UMI/BUS/technology modes, or fusion detection.
- Minimal CLI surface; some kallisto options are missing.
- Not a drop-in replacement yet; use for experimentation and validation, not as a production substitute.

Current real-data status:
- Paired-end parity is exact on deterministic prefixes through `262144` pairs from the checked-in real dataset under `data/`.
- Full-file paired-end quant is still slightly off from `kallisto` on the same index and reads.
- Full-file startup/runtime is still materially slower than `kallisto`; the remaining gap is mostly in index loading/materialization.

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
- If exact minimizer anchoring fails, match candidates are scanned across the valid minimizer window
  (`rel_pos - (k-g) .. rel_pos`) before rejection.
- Online intersection tracking with `intersection_empty` skip parity.
- Jump logic records a hit before skipping to match kallisto’s jump behavior.
- Special/overcrowded minimizers handled via D-list / special unitig checks, with upstream-style
  Bifrost fallback for overcrowded minimizers instead of a global exact-kmer fallback.


### Real-data benchmark (current checked-in dataset)
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
./target/release/kallistors-cli quant \
  -i data/gencode.v49_kallisto.idx -o /tmp/kallistors_full -t 32 \
  data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  "data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz"
```

Latest results (2026-04-15, Linux x86_64, `7950X3D`, `-t 32`):
- Full-file paired quant:
  `kallisto` wall `59.21s`
  `kallistors` wall `190.76s`
- Full-file `run_info.json`:
  `kallisto` `n_pseudoaligned = 4244771`, `n_unique = 276251`
  `kallistors` `n_pseudoaligned = 4244773`, `n_unique = 276251`
- Deterministic paired-prefix parity:
  exact through `262144` pairs on `data/subsets/`

### Parity tests
- Synthetic parity: `variants_parity` allows a 1-read drift in aligned count; EC sets must match when aligned.
- Other parity tests keep strict equality on aligned counts and EC sets.

### As a Crate

Add to `Cargo.toml`:
```toml
[dependencies]
kallistors = { git = "https://github.com/deminden/kallistors" }
```

Use in code:
```rust
use kallistors::index::Index;
use kallistors::io::open_fastq_reader;
use kallistors::pseudoalign::{
    build_bifrost_index, pseudoalign_paired_bifrost_with_options, PseudoalignOptions, Strand,
};
use kallistors::quant::{em_quantify, EcCountsInput, QuantOptions};

let index_path = "path/to/index.idx";
let mut reads_1 = open_fastq_reader("reads_1.fq.gz".as_ref())?;
let mut reads_2 = open_fastq_reader("reads_2.fq.gz".as_ref())?;
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
let quant = em_quantify(
    &input,
    &lengths,
    None,
    None,
    QuantOptions::default(),
)?;

println!("reads processed: {}", ec.reads_processed);
println!("reads aligned: {}", ec.reads_aligned);
println!("targets: {}", quant.est_counts.len());
```

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
  --index data/gencode.v49_kallisto.idx \
  --reads data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
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
