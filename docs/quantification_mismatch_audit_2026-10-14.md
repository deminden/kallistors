# Quantification mismatch audit — 2026-10-14

This report records two compatibility defects before their fixes. The paired-end effective-length
model has the larger impact. A premature jump during a single-end fallback probe explains the
two-read TTN discrepancy. Matching `run_info.json` counters does not establish
abundance parity; the earlier comparison report only checked those counters for paired-end runs.
The subsequent implementation and validation are recorded in
[quantification_fixes_2026-10-14.md](quantification_fixes_2026-10-14.md).

This audit uses the previously downloaded, fastp-trimmed human datasets and the same upstream-built
GENCODE v50 ALL-transcript index, containing 670,670 targets. Quantification checks use eight
threads; single-end checks use trimmed mate 1 with `-l 200 -s 20`. No bias correction or bootstrap
sampling is enabled. These additional checks diagnose correctness, and are not new speed medians.

## Findings

| Priority | Finding | Evidence | Recommended change |
| --- | --- | --- | --- |
| High | Paired-end EM uses a fitted Gaussian instead of the observed fragment histogram | All three paired-end EC count maps agree exactly, but short-transcript TPMs differ by tens of percent | Pass empirical, transcript-length-truncated fragment means into EM and effective-length calculation |
| Medium | A single-end fallback probe initiates another jump | Two identical reads miss a unique TTN k-mer at read position 63; disabling jumps restores the unique EC | Complete the next valid k-mer's fallback probe without starting a jump, in both scanner implementations |
| No corrective action | Tiny differences in the two smaller single-end samples | High-precision output differs by at most about 1e-11 TPM; rounded estimates agree | Use a numerical tolerance and distinguish text precision from algorithmic changes |

`n_unique` counts reads assigned to an EC containing one transcript (paired fragments in paired
mode). It does not count distinct transcripts.

## Paired-end effective lengths: material differences

Representative differences from the official kallisto release and current kallistors:

| Dataset | Transcript | kallisto TPM | kallistors TPM | Relative change |
| --- | --- | ---: | ---: | ---: |
| ERR2260014 | ENST00000427161.1, short lncRNA | 64,342.7 | 49,999.4 | -22.3% |
| DRR023395 | MT-TP / ENST00000387461.2 | 3,675.96 | 1,742.44 | -52.6% |
| SRR13638690 | SNORD89 / ENST00000390981.1 | 7,828.39 | 4,358.17 | -44.3% |

These are not only near-zero transcript differences: the examples include approximately 22, 24,
and 292 estimated fragments, respectively. The short-transcript differences also change the global
TPM normalization and therefore affect other transcripts.

The total TPM mass shift, defined as half the sum of absolute TPM differences divided by the
reference TPM sum, is 4.46%, 1.05%, and 3.94% for the compact, medium, and larger paired datasets.
This is a compatibility discrepancy, not a measurement of which tool is biologically more accurate.

The underlying nonzero EC count maps match exactly for all paired datasets:

| Dataset | Nonzero ECs | Differing EC counts |
| --- | ---: | ---: |
| ERR2260014 | 3,604 | 0 |
| DRR023395 | 181,044 | 0 |
| SRR13638690 | 190,711 | 0 |

The cause is visible in [the quant command](../crates/kallistors/src/commands/quant.rs): paired
observations are reduced to mean and SD before calling EM. [EM](../crates/kallistors/src/quant/mod.rs)
then calls `mean_frag_lens_by_transcript`, fitting and truncating a Gaussian. The existing empirical
`fragment_length_hist` is written to HDF5 metadata but is not used by EM.

Upstream [MinCollector::compute_mean_frag_lens_trunc](https://github.com/pachterlab/kallisto/blob/v0.52.0/src/MinCollector.cpp)
uses cumulative empirical mass and counts. For a transcript of length L, its conditional fragment
mean is `sum(f * histogram[f], f <= L) / sum(histogram[f], f <= L)`; effective length is then
`L - conditional_mean + 1`, with upstream's short-length fallback.

The compact dataset isolates the model difference particularly well. Both stored HDF5 histograms
are exactly equal: 286 fragment observations, mean 170.7937062937063 bp. For the 157 bp transcript
ENST00000427161.1, the empirical conditional mean is 121.59848484848484 bp, giving effective length
36.401515151515156. Kallistors instead produces effective length 48.16748837645359. Reconstructing
its Gaussian with the same observed mean and sample SD reproduces that effective length.

There is also a sampling-policy difference to consider after fixing the model. Kallistors collects
all qualifying fragment observations; upstream stops near its 10,000-observation goal and can
overshoot with multiple workers. In these eight-thread checks, stored histogram observation counts
were 15,361 vs 15,653 for DRR023395 and 11,644 vs 46,496 for SRR13638690 (kallisto vs kallistors).
More observations are not inherently incorrect, but exact upstream compatibility requires a
deliberate sampling policy. Upstream's worker scheduling can also change the sampled histogram.

## Single-end TTN discrepancy: a small but real scanner defect

For SRR13638690 single-end input, the entire EC count comparison differs in just two classes:

| Transcript IDs in shared index | kallisto count | kallistors count |
| --- | ---: | ---: |
| `{96635}` | 1,428 | 1,426 |
| `{96633,96635,96636,96637,96638,96641,96642,96643}` | 197,659 | 197,661 |

The reads are `SRR13638690.2809610` and `SRR13638690.3196574`. Their trimmed sequences are identical,
126 bases long. They are preserved in the small [reproduction fixture](../crates/kallistors/tests/fixtures/single_ttn_gencode50.fixture).
Quality characters are replaced with `I`; pseudoalignment does not use them. The fixture requires
the shared GENCODE v50 index because the numeric transcript and unitig IDs are reference-specific.

Both the official kallisto binary and a diagnostic build assign the two reads uniquely to transcript
96635, ENST00000425332.4 / TTN-206. Current kallistors assigns them to eight TTN transcripts. Running
its diagnostic pseudoaligner with `--no-jump` restores both unique assignments.

The Rust trace visits positions 0 through 33 and then accepts a jump to the end. Upstream's accepted
hit stream includes position 63, where unitig 377463 has the distinguishing singleton EC `{96635}`.
The Rust [general scanner](../crates/kallistors/src/pseudoalign/core/read.rs) and
[fast scanner](../crates/kallistors/src/pseudoalign/core/fast.rs) retain `backoff_until`, but still
evaluate a new jump on the immediate fallback probe. The actual v0.52.0
[upstream loop](https://github.com/pachterlab/kallisto/blob/v0.52.0/src/KmerIndex.cpp) performs one
fallback probe, then returns to its outer loop: an inner `nextPos` shadows the outer stop variable.
Its comment suggests a longer bounded scan, but that is not its executed behavior. Suppressing
every jump until the endpoint was tested and rejected because it lost medium-sample alignments.
The fix completes the immediate fallback probe before allowing the next outer-loop jump.

The two-read difference changes the most affected TTN isoform by about 269 estimated counts and
1.11 TPM after EM redistribution; this is not 269 differently mapped reads. TTN gene-level TPM
changes by only about 0.016. The impact is much smaller than the paired-end effective-length issue,
but the deterministic reproducer makes this a useful compatibility regression to fix.

Reproduce the count discrepancy with the existing index:

```bash
target/tools/kallisto-0.52.0 quant --single -l 200 -s 20 -t 8 \
  -i data/srr-small/gencode.v50.kallisto.idx \
  -o target/ttn-reproducer-kallisto \
  crates/kallistors/tests/fixtures/single_ttn_gencode50.fixture

target/release/kallistors quant --single -l 200 -s 20 -t 8 \
  -i data/srr-small/gencode.v50.kallisto.idx \
  -o target/ttn-reproducer-kallistors \
  crates/kallistors/tests/fixtures/single_ttn_gencode50.fixture
```

Before the fix, `n_processed / n_pseudoaligned / n_unique` values were `2 / 2 / 2` for kallisto and
`2 / 2 / 0` for kallistors. Both tools now produce `2 / 2 / 2`. The diagnostic
`pseudoalign --no-jump --fragment-length 200 --threads 8` check also produces the singleton class
with count 2.

## Harmless output differences and audit artifacts

The released kallisto TSV writer uses six significant digits, while kallistors writes full float
precision. Counting unequal TSV values alone therefore overstates discrepancies. A diagnostic
upstream build using 17 significant digits gives maximum absolute TPM differences of 7.28e-12 and
9.49e-12 for ERR2260014 and DRR023395 single-end outputs, respectively. These need no quantification
fix; small changes in floating-point accumulation order are expected.

[compare_abundance.py](../scripts/compare_abundance.py) now supports `--summary-json`. It reports
effective-length differences, discrepancies after rounding to a common precision, absolute-error
quantiles, relative TPM errors among expressed transcripts, missing targets, and TPM mass shift.
Using only medians over all transcripts hides discrepancies when most reference targets have zero
abundance.

```bash
python3 scripts/compare_abundance.py \
  --kallisto target/quant-mismatch-investigation/ERR2260014-paired-kallisto/abundance.tsv \
  --kallistors target/quant-mismatch-investigation/ERR2260014-paired-kallistors/abundance.tsv \
  --summary-json target/quant-mismatch-investigation/compact-paired-summary.json
```

Audit summaries, full EC maps, observed fragment histograms, and per-read traces are kept under
ignored `target/quant-mismatch-investigation/`. The source diagnostic patch is preserved in
[kallisto_v0.52.0_quant_audit.patch](../scripts/diagnostics/kallisto_v0.52.0_quant_audit.patch).
It applies to upstream tag v0.52.0, adds optional diagnostic output and high-precision TSV output,
and includes build fixes for the current Apple Silicon toolchain. EC and histogram output is
enabled with `KALLISTO_DIAGNOSTIC_ECS=1`. Per-read hit tracing uses
`KALLISTO_DIAGNOSTIC_HITS=<path>`; singleton read capture also needs
`KALLISTO_DIAGNOSTIC_READS=<path>` and `KALLISTO_DIAGNOSTIC_TID=96635` for this index.

The pre-fix artifacts remain available for comparison with the fixed results.
