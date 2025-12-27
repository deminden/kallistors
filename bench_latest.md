# Real-data benchmark vs kallisto
Generated: 2025-12-27T20:01:07.370567

## Environment
- OS: macOS-15.7.2-arm64-arm-64bit
- Machine: arm64
- Processor: arm
- CPU cores: 8

## Inputs
- Index: data/Human_kallisto_index
- Reads: data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
- Fragment length: 200 +- 20
- Threads: 8

## Commands
```bash
kallisto quant --single -l 200 -s 20 -i data/Human_kallisto_index -o /var/folders/vh/5r90mw8d7y99ysqy4r6hd0gr0000gn/T/kallistors-bench-htnonw32/kallisto -t 8 data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
target/release/kallistors-cli quant --single -l 200 -s 20 -i data/Human_kallisto_index -o /var/folders/vh/5r90mw8d7y99ysqy4r6hd0gr0000gn/T/kallistors-bench-htnonw32/kallistors -t 8 data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
```

## Timings
- kallisto wall: 11.106s, cpu: 30.802s
- kallistors wall: 36.013s, cpu: 36.408s

## run_info parity
| field | kallisto | kallistors |
| --- | --- | --- |
| n_targets | 328868 | 328868 |
| n_processed | 43817 | 43817 |
| n_pseudoaligned | 41444 | 41113 |
| n_unique | 2420 | 2111 |
| p_pseudoaligned | 94.6 | 93.8 |
| p_unique | 5.5 | 4.8 |

## abundance.tsv parity
- Missing transcripts: 0
- Max eff_length abs diff: 4.94999e-05

Percentiles for relative error (abs(a-b)/max(1,a,b)):

| metric | p50 | p90 | p99 | max |
| --- | --- | --- | --- | --- |
| est_counts | 0 | 0 | 1.43e-07 | 1 |
| tpm | 0 | 0 | 0.00882 | 1 |
