# Paired fragment-length sampling

## Behavior

Paired quantification estimates its fragment-length distribution from uniquely assigned pairs
whose existing graph matches support a valid fragment length. The histogram determines
transcript-specific conditional means and effective lengths.

The default follows ordinary short-read quantification in original kallisto v0.52.0:

1. Aim for 10,000 qualifying fragment observations.
2. Form logical sequence buffers of 8 MiB, counting both mate sequences and one
   terminating byte per mate, with upstream's strict `< buffer_size` fit rule.
3. At the start of each logical buffer, snapshot the remaining quota.
4. Collect qualifying observations up to that buffer's quota.
5. Publish the buffer's observation count when it finishes.

Concurrent batches can start with the same remaining quota and overshoot the goal. Their
scheduling can change the final histogram and, consequently, abundance estimates. Original
kallisto uses this behavior too; identical EC counts alone do not guarantee identical abundance.
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

Paired `quant` reports the actual observation count, mean, and standard deviation. The actual
count may exceed the requested goal with multiple threads.

## Validation

Regression tests cover byte-boundary splits, cross-chunk quotas, variable mate lengths,
oversized pairs, remaining batch quotas, concurrent overshoot, short inputs, exact EC
invariance between sampled and all-observation runs, and serial/multithreaded execution using a
real kallisto-built synthetic index. The complete EC map remains the correctness gate for speed
experiments. To isolate EM arithmetic, replay the same ECs and histogram with both implementations.
Within kallistors, `--fragment-length-samples 0` provides a deterministic all-observation
fragment model for a fixed input; original kallisto requires diagnostic controls for a matched
model.

Sampling is retained for upstream compatibility. Its speed effect is assessed with complete
quantification runs.
