# Real-data benchmark vs kallisto
Generated: 2026-02-14

## Environment
- OS: macOS arm64
- CPU cores: 8

## Inputs
- Reads: `data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz`
- Read count: 43,817
- Mean read length: 168.904 (used `-l 169 -s 20`)
- Index: `data/Human_kallisto_index`
- Threads: `-t 8`

## Benchmark (no debug flags)
Commands:
```bash
kallisto_src/build/src/kallisto quant \
  -i data/Human_kallisto_index -o data/kallisto_bench_run_20260214_rerun1 \
  --single -l 169 -s 20 -t 8 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz

./target/release/kallistors-cli quant \
  -i data/Human_kallisto_index -o data/kallistors_bench_run_20260214_rerun1 \
  --single -l 169 -s 20 -t 8 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
```

Timings:
- kallisto: real 10.32s, user 28.17s, sys 2.00s
- kallistors: real 41.01s, user 32.77s, sys 11.04s

run_info:
- kallisto `n_pseudoaligned`: 41626 / 43817
- kallistors `n_pseudoaligned`: 41620 / 43817

Accuracy (TPM filter `max(k_tpm, o_tpm) > 1`; n = 7699):
- TPM Pearson: 0.999999
- TPM MAE: 0.071
- est_counts Pearson: 1.000000
- est_counts MAE: 0.0019

## Read-level parity (25k sample, same index)
- Before relaxed minimizer-window matching: 79 mismatches (40 `no_hits_ok`, 39 `ok_no_hits`)
- After jump/backoff and special-overcrowded alignment work: 4 mismatches (4 `no_hits_ok`, 0 `ok_no_hits`)
- Artifacts:
  - `data/kallisto_vs_kallistors_read_diff_25000_step14.tsv`
  - `data/kallisto_vs_kallistors_read_diff_no_hits_ok_25000_step14.txt`
  - `data/kallisto_vs_kallistors_read_diff_ok_no_hits_25000_step14.txt`
