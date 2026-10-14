# Quantification fixes — 2026-10-14

This follow-up records the fixes identified by the
[quantification mismatch audit](quantification_mismatch_audit_2026-10-14.md).

Paired-end EM now uses the observed fragment histogram to calculate transcript-specific conditional
fragment means and effective lengths. The histogram is also used consistently by bias correction;
bootstrap quantification continues to use the resulting effective lengths. Single-end quantification
retains its supplied Gaussian fragment model.

The fast and general read scanners now complete the next valid k-mer's fallback probe without
starting another jump. They then resume ordinary jumping, following the actual upstream v0.52.0
loop. A scan that suppressed jumping over the entire interval was rejected during full-read
validation, because it lost medium-sample alignments.

Regression coverage includes empirical conditional means at histogram boundaries, command-level
paired quantification on two observed fragment lengths, the two TTN splice reads, and an AGAP12P
read that catches the over-broad fallback scan. The full-reference scanner regression runs both
scanner implementations against the downloaded GENCODE v50 index:

```bash
cargo test --locked --workspace --all-features
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
cargo build --release --locked
KALLISTORS_TEST_BIN="$PWD/target/release/kallistors" \
  cargo test --locked --test single_ttn_regression -- --ignored
```

The workspace test suite, strict Clippy command, and the explicit full-reference scanner regression
all passed. The latter uses both fast and general scanners and checks three processed, aligned,
and unique reads: two TTN-206 reads and one AGAP12P-205 read. Normal CI leaves this downloaded-index
regression ignored; the command above enables it explicitly.

## Complete EC count validation

All nonzero transcript-set/count entries were compared with the diagnostic upstream v0.52.0
build. This checks the entire mapping evidence passed into EM, beyond the headline counters.
All six cases match exactly:

| Dataset | Mode | Nonzero ECs | Differing EC counts | Sum of absolute count differences |
| --- | --- | ---: | ---: | ---: |
| ERR2260014 | paired | 3,604 | 0 | 0 |
| ERR2260014 | single | 3,534 | 0 | 0 |
| DRR023395 | paired | 181,044 | 0 | 0 |
| DRR023395 | single | 146,676 | 0 | 0 |
| SRR13638690 | paired | 190,711 | 0 | 0 |
| SRR13638690 | single | 192,645 | 0 | 0 |

The SRR13638690 single-end singleton TTN count is restored from 1,426 to upstream's 1,428;
the corresponding ambiguous class drops from 197,661 to 197,659. No other EC counts change.
Full EC maps and the zero-difference summary are preserved under ignored
`target/quant-fixes-2026-10-14/ec/`.

## Abundance validation

Checks use the same trimmed reads and shared GENCODE v50 ALL-transcript index (670,670 targets),
eight threads, no bias correction, and no bootstraps. Single-end uses mate 1 and `-l 200 -s 20`.
The diagnostic upstream writer uses 17 significant digits to separate ordinary TSV rounding from
algorithmic differences. Paired fragment histograms can vary with upstream worker scheduling.

TPM mass shift is half the sum of absolute transcript TPM differences divided by upstream's
total TPM, expressed as a percentage. Results are:

| Dataset / mode | TPM mass shift after fix | Maximum absolute TPM difference |
| --- | ---: | ---: |
| ERR2260014 paired | 2.37e-15% | 1.36e-12 |
| ERR2260014 single | 2.56e-15% | 7.28e-12 |
| DRR023395 paired | 0.0100% | 4.6633 |
| DRR023395 single | 3.12e-14% | 1.50e-11 |
| SRR13638690 paired | 0.08024% | 69.3149 |

The paired discrepancy was about 4.46% for ERR2260014, 1.05% for DRR023395, and 3.94% for
SRR13638690 before the fix. The large-sample result compares retained official-release outputs from
the final measured paired repetition; its six-digit text rounding also contributes a small error.
For example, SNORD89 now has 7,794.93 TPM upstream versus 7,862.31 in kallistors, compared with
7,828.39 versus 4,358.17 in the earlier check. Upstream's sampled histogram changes between runs.
The remaining paired-sample discrepancy comes from the fragment histogram: upstream stops near
its observation goal and can overshoot with multiple workers; kallistors collects all qualifying
observations. This fix applies the empirical calculation to the existing histogram. It does not
change observation sampling. Exact EC parity therefore does not imply exact abundance parity for
every paired-end dataset.

Sampling is primarily a compatibility decision here. The large sample collects 46,496 fragment
observations from 4,411,003 pairs. Collection uses the mates' existing best matches, EC-block
lookups, arithmetic, and a histogram increment; it does not remap the reads. Capping observations
would leave FASTQ reading and pseudoalignment unchanged. A speed benefit has not been measured,
but inspection suggests it would be small. The current implementation retains all observations.

Audit JSON files are preserved under `target/quant-fixes-2026-10-14/`. The earlier short-RNA TPM
differences and the diagnostic source patch are documented in the linked mismatch audit.

## Eight-thread repeated benchmark

The larger SRR13638690 benchmark uses one warmup and three measured runs per binary and mode,
with alternating binary order and the median reported. The final measured abundance outputs are
retained outside the timed region:

```bash
python3 scripts/benchmark_local_reads.py \
  --samples SRR13638690_larger --thread-counts 8 --keep-last-output \
  --out-prefix target/benchmarks/local-real-2026-10-14/srr13638690-8thread-fixed
```

At the user's request, the remaining single-end repetitions were stopped. The paired configuration
had already completed its warmups and all three measured repetitions:

| Mode | Threads | kallisto median | kallistors median | kallisto / kallistors |
| --- | ---: | ---: | ---: | ---: |
| paired | 8 | 46.694 s | 78.123 s | 0.598x |

Measured paired elapsed times were 46.694, 47.252, and 46.436 seconds for kallisto, and 78.123,
77.373, and 82.259 seconds for kallistors. Kallistors remains slower on this paired workload.
These fixes address correctness; the change from the earlier 80.05-second kallistors median is
not sufficient evidence of a speed improvement.

Both single-end warmups completed with identical counters: 4,411,003 processed, 4,246,922
pseudoaligned, and 205,821 unique. No new single-end timing median or full large-sample single-end
abundance comparison is reported. The complete single-end EC comparison and isolated TTN
quantification regression described above passed.

The final paired quantification outputs are under
`target/benchmarks/local-real-2026-10-14/srr13638690-8thread-fixed.outputs/`.
Paired run and summary TSVs share the `srr13638690-8thread-fixed` prefix. Because the benchmark was
terminated before its end-of-run write, these TSV timings were recovered from console output,
rounded to 0.001 seconds; `.status.json` records that limitation and the single-end stop. The runner
now saves each completed run and configuration outside the timed region so future interruptions
preserve completed results automatically.
