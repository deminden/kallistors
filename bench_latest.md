# Real-data benchmark vs kallisto
Generated: 2025-12-27T20:01:07.370567

## Environment
- OS: macOS-15.7.2-arm64-arm-64bit
- Machine: arm64
- Processor: arm
- CPU cores: 8

## 2025-12-29 benchmark
Inputs:
- Reads: data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
- Mean read length: 168.90 (used -l 169 -s 20)

Commands:
```bash
kallisto_src/build/src/kallisto quant \
  -i data/kallisto_index -o data/kallisto_bench_run \
  --single -l 169 -s 20 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz

./target/release/kallistors-cli quant \
  -i data/kallisto_index -o data/kallistors_bench_run \
  --single -l 169 -s 20 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
```
Timings:
- kallisto real 19.26s, user 17.41s, sys 1.47s
- kallistors real 36.38s, user 28.25s, sys 7.13s
run_info:
- kallisto n_pseudoaligned 41626 / 43817
- kallistors n_pseudoaligned 41666 / 43817
Accuracy (TPM filter max(k_tpm, o_tpm) > 1; n = 7824):
- TPM Pearson 0.974175, TPM MAE 11.206
- est_counts Pearson 0.998049, est_counts MAE 0.303

### Original index (data/Human_kallisto_index)
Commands:
```bash
kallisto_src/build/src/kallisto quant \
  -i data/Human_kallisto_index -o data/kallisto_bench_run_old \
  --single -l 169 -s 20 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz

./target/release/kallistors-cli quant \
  -i data/Human_kallisto_index -o data/kallistors_bench_run_old \
  --single -l 169 -s 20 \
  data/SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz
```
Timings:
- kallisto real 19.50s, user 17.50s, sys 1.50s
- kallistors real 39.48s, user 29.00s, sys 8.80s
run_info:
- kallisto n_pseudoaligned 41626 / 43817
- kallistors n_pseudoaligned 41665 / 43817
Accuracy (TPM filter max(k_tpm, o_tpm) > 1; n = 7840):
- TPM Pearson 0.974031, TPM MAE 11.614
- est_counts Pearson 0.997615, est_counts MAE 0.325
