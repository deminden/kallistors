# Real-data benchmark vs kallisto
Generated: 2026-05-17

## Environment
- OS: Linux x86_64
- CPU: AMD Ryzen 9 7950X3D
- Threads used: 32
- Build: release profile with fat LTO and `.cargo/config.toml` `target-cpu=native`

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
- kallisto: real 59.12s, user 161.25s, sys 7.12s
- kallistors: real 30.76s, user 404.40s, sys 9.31s

run_info:
- kallisto `n_pseudoaligned`: 4244771 / 4408640
- kallistors `n_pseudoaligned`: 4244771 / 4408640
- kallisto `n_unique`: 276251
- kallistors `n_unique`: 276251

Current `kallistors` stage timings:
- `index_header_parse 0.000s`
- `graph_decode 3.872s`
- `minimizer_count_pass 0.842s`
- `minimizer_fill_pass 5.454s`
- `fastq_read_decompress 9.568s`
- `pseudoalign 12.678s`
- `ec_merge 0.125s`
- `em 7.731s`

Delta vs local upstream `kallisto` on this run:
- Wall time: `59.12s -> 30.76s`, `28.36s` faster for `kallistors`.
- Speedup: about `1.92x` faster, or about `48%` less wall time.

Delta vs the previous public `kallistors` benchmark:
- Wall time: `70.73s -> 30.76s`, `39.97s` faster, about `56.5%`.

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
- Quant reuses transcript metadata already loaded by the Bifrost index path and moves EC
  classes/counts into EM instead of cloning them.
- The Bifrost loader skips the redundant graph pre-scan on the supported `k <= 32` path, lazily
  allocates shade metadata, reuses graph-node payload buffers, and avoids paired positional payload
  loading when paired fragment estimation only needs flat EC block bounds.
