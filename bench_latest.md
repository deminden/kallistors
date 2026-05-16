# Real-data benchmark vs kallisto
Generated: 2026-05-16

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

./target/release/kallistors quant \
  -i data/gencode.v49_kallisto.idx -o /tmp/kallistors_full -t 32 \
  data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade.gz \
  "data/SRR13638690_RNA_seq_of_homo_sapiens_temporal_muscle_of_low_grade (2).gz"
```

Timings:
- kallisto: real 57.14s
- kallistors: real 70.73s, user 1163.61s, sys 16.25s

run_info:
- kallisto `n_pseudoaligned`: 4244771 / 4408640
- kallistors `n_pseudoaligned`: 4244771 / 4408640
- kallisto `n_unique`: 276251
- kallistors `n_unique`: 276251

Current `kallistors` stage timings:
- `index_header_parse 2.183s`
- `graph_decode 6.855s`
- `minimizer_count_pass 1.358s`
- `minimizer_fill_pass 9.062s`
- `fastq_read_decompress 29.023s`
- `pseudoalign 36.610s`
- `ec_merge 0.192s`
- `em 8.075s`

Delta vs the pre-optimization `kallistors` control from this round of work:
- Wall time: `77.60s -> 70.73s`, `6.87s` faster, about `8.9%`.
- Pseudoalign stage: `41.488s -> 36.610s`, about `11.8%`.
- EM stage: `9.673s -> 8.075s`, about `16.5%`.

## Read-level parity
- Full-file paired parity is exact against `kallisto` on the checked-in real dataset.
- Deterministic paired prefixes are exact through `262144` pairs in `data/subsets/` with
  `--threads 32`.

## Notable improvements behind this result
- Bifrost-style retry on probe/backoff misses after prior evidence exists, fixing the last full-file paired mismatch.
- `flate2` switched to the `zlib-rs` backend.
- Threaded workers now accumulate directly into long-lived `EcCounts`.
- Threaded FASTQ transport now uses packed/reusable batches with one contiguous byte buffer plus per-record offsets instead of owned `FastqRecord` payloads.
- Fast pseudoalignment reuses encoded k-mer codes through minimizer candidate lookup, match-cache
  keying, and jump/middle/scan probes.
- Hot pseudoalignment environment flags are cached once per process instead of being looked up in
  per-read/per-k-mer paths.
- The EM loop pre-splits singleton/nonzero multi-transcript ECs and stores multi-EC transcript and
  weight metadata in contiguous arrays.
- Tiny hot direct-mapped lookup caches for MPHF minimizer lookup and EC block lookup were enlarged
  to reduce collisions in the common path.
