# Real-data benchmark vs kallisto
Generated: 2026-04-18

## Environment
- OS: Linux x86_64
- CPU: AMD Ryzen 9 7950X3D
- Threads used: 32

## Inputs
- Reads:
  - `data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz`
  - `data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz`
- Read pairs: 4,408,640
- Reference transcripts: `data/gencode.v49.transcripts.fa.gz`
- Index: `data/gencode.v49_kallisto.idx`
- Threads: `-t 32`

## Benchmark (no debug flags)
Commands:
```bash
kallisto_src/build/src/kallisto quant \
  -i data/gencode.v49_kallisto.idx -o /tmp/kallisto_full -t 32 \
  data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  "data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz"

./target/release/kallistors-cli quant \
  -i data/gencode.v49_kallisto.idx -o /tmp/kallistors_full -t 32 \
  data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  "data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz"
```

Timings:
- kallisto: real 56.24s
- kallistors: real 68.63s, user 1006.42s, sys 17.45s, peak RSS 5,653,304 KB

run_info:
- kallisto `n_pseudoaligned`: 4244771 / 4408640
- kallistors `n_pseudoaligned`: 4244771 / 4408640
- kallisto `n_unique`: 276251
- kallistors `n_unique`: 276251

Current `kallistors` stage timings:
- `index_header_parse 2.113s`
- `graph_decode 6.859s`
- `minimizer_count_pass 1.390s`
- `minimizer_fill_pass 9.002s`
- `fastq_read_decompress 34.613s`
- `pseudoalign 47.39s`

## Read-level parity
- Full-file paired parity is exact against `kallisto` on the checked-in real dataset.
- Deterministic paired prefixes are exact through `262144` pairs in `data/subsets/`.

## Notable improvements behind this result
- Bifrost-style retry on probe/backoff misses after prior evidence exists, fixing the last full-file paired mismatch.
- `flate2` switched to the `zlib-rs` backend.
- Threaded workers now accumulate directly into long-lived `EcCounts`.
- Threaded FASTQ transport now uses packed/reusable batches with one contiguous byte buffer plus per-record offsets instead of owned `FastqRecord` payloads.
