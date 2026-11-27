# Paired fragment-length sampling

## Behavior

Paired quantification estimates its fragment-length distribution from uniquely assigned pairs
whose existing graph matches support a valid fragment length. The histogram determines
transcript-specific conditional means and effective lengths.

The default follows ordinary short-read quantification in original kallisto v0.52.0:

1. Aim for 10,000 qualifying fragment observations.
2. Form logical sequence buffers of 8 MiB, counting both mate sequences and one
   terminating byte per mate, with upstream's strict `< buffer_size` fit rule.
3. Number the logical buffers in input order.
4. Give buffer `i` a quota of the goal minus the observations kept by buffers
   `1..i`. The first buffer receives the full goal.
5. Keep each buffer's first qualifying observations, in read order, up to its quota.

Original kallisto snapshots a shared remaining quota at the start of each buffer and
publishes the buffer's count when it finishes. Concurrent buffers can therefore start
with the same quota and overshoot the goal, and their scheduling can change the final
histogram and, consequently, abundance estimates. The quota above is the one that shared
counter gives when one worker runs, and kallistors applies it for any number of workers,
so the sample is the first qualifying observations in input order: the same in every run
for a fixed input and goal, at every thread count, and never above the goal. This matches
original kallisto on one thread. With several threads original kallisto's sample depends
on timing, so identical EC counts alone do not guarantee identical abundance.

Before 0.7.0, buffer `i`'s quota followed buffers `1..=i - W`, where `W` was the worker
count: the quota a shared counter gives when `W` workers run in step. That sample was fixed
for a given thread count but differed between thread counts.
Parser chunks retain their smaller size for input throughput. A logical sampling buffer
can span several parser chunks or split one; its quota is carried across all those chunks.
Once sampling finishes, mapping uses the usual small work batches. No sequence copies
are needed for split chunks. `--fragment-length-samples 0` retains ordinary work batches.

Every read pair is processed for EC counts after estimation reaches its goal. Single-end
quantification continues to use its supplied `-l` / `-s` fragment model.

The upstream rules are in
[ReadProcessor::processBuffer](https://github.com/pachterlab/kallisto/blob/4e9f29cf3b021260415430c057a22469ca081391/src/ProcessReads.cpp),
including `findFragmentLength`, `flengoal`, and the fragment histogram update.

## Configuration

Both `quant` and `pseudoalign` accept `--fragment-length-samples N`. The default is `10000`;
`0` collects all qualifying observations, which is useful for comparing the
counting and EM implementations with an identical observation model. For the Rust API, set
`PseudoalignOptions::fragment_length_samples`.

Paired `quant` reports the actual observation count, mean, and standard deviation.

## Validation

Regression tests cover byte-boundary splits, cross-chunk quotas, variable mate lengths,
oversized pairs, remaining batch quotas, ordered quotas under shuffled completion order,
short inputs, exact EC
invariance between sampled and all-observation runs, and serial/multithreaded execution using
self-contained Rust-built synthetic indexes. Opt-in upstream comparisons use a separately
built original-kallisto index. The complete EC map remains the correctness gate for speed
experiments. To isolate EM arithmetic, replay the same ECs and histogram with both implementations.
Within kallistors, the default sampled model is fixed for a given input at every thread
count, and `--fragment-length-samples 0` provides an all-observation model; original
kallisto requires diagnostic controls for a matched model.

Sampling is retained for upstream compatibility. Its speed effect is assessed with complete
quantification runs.

## Interpreting transcript-level differences

Identical counted ECs can yield different abundance estimates when fragment histograms differ.
Sampled effective lengths affect transcript counts and TPM normalization, especially for short
transcripts. Ambiguous isoforms can also respond to model changes and EM convergence; inspect
gene totals alongside transcript-level differences.

The [latest numerical comparison](kallisto-comparison-results.md#native-abundance-agreement)
reports complete transcript/gene distributions and links every measured replicate, realized
fragment model, within-program variation and one-versus-eight-thread comparison. Inspect counts,
effective lengths, expressed-target percentiles and within-program variation together.
Agreement with one upstream estimate does not establish biological accuracy.

`--fragment-length-samples 0` uses all qualifying observations. It can add processing work
and produces a different model from the default sampled run. Replaying an identical saved histogram and EC map with matched EM
settings isolates implementation differences from sampling differences.
