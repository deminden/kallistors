# Kallisto compatibility TODOs

This checklist is ordered. Do not proceed to the next item until the current item is validated against kallisto output.

## Phase A: Index parsing parity (metadata)
- [x] Parse and report index version (should be 13)
- [x] Parse and report k-mer length (match `kallisto inspect`)
- [x] Parse and report minimizer length (match `kallisto inspect`)
- [x] Parse and report number of unitigs (match `kallisto inspect`)
- [x] Parse and report number of k-mers (match `kallisto inspect`)
- [x] Parse and report number of targets and total length (match `kallisto inspect` + transcript FASTA)

Acceptance test:
- `kallisto inspect data/synthetic.idx` matches `kallistors-cli index-info` for all reported fields.

## Phase B: EC data parity
- [x] Parse EC mapping and emit list in kallisto-compatible order (debug-only)
- [x] Support loading kallisto EC text files (`output.ec.txt`) and validate ordering
- [x] Export EC list to text and compare with kallisto EC list (if available)
- [x] Add debug-only EC compare utility (remove later)

Notes:
- kallisto does not expose an EC list export for a plain index; comparison requires an external tool or custom build.

Acceptance test:
- EC list checksum matches a known reference for a synthetic index.

## Phase C: Pseudoalignment parity (single-end)
- [x] Implement naive k-mer lookup against kallisto index structures (debug path)
- [x] Compute EC per read and count (debug path)
- [ ] Validate EC counts for toy reads against kallisto pseudoalignment output

Acceptance test:
- EC counts identical for a fixed synthetic reads set.

Notes:
- kallisto does not expose EC counts directly; validation will require a custom tool or patched build.

## Phase D: Quant parity (basic EM)
- [x] Implement EM on EC counts
- [x] Validate against kallisto quant on toy data (within tolerance)

Acceptance test:
- TPMs/est_counts within tolerance for synthetic data.

Notes:
- `kallisto quant --bias` segfaults on the current simple synthetic dataset; bias parity remains blocked upstream.

## Phase E: Real-data subset parity
- [x] Normalize the current real-data inputs into deterministic small subsets under `data/subsets/`
- [x] Add a scripted parity runner for `kallisto` vs `kallistors` on those subsets
- [x] Validate exact `run_info.json` parity on the smallest single-end subset
- [x] Promote to larger single-end subsets and capture the first size that diverges
- [x] Promote to paired-end subsets after single-end parity is stable
- [x] If real-data parity diverges, save the minimal failing subset and collect per-read traces

Acceptance test:
- `scripts/make_real_subsets.py` produces the tracked subset layout and
  `scripts/real_subset_parity.py --mode single` reports exact `n_processed`,
  `n_pseudoaligned`, and `n_unique` parity for the smallest subset before scaling up.

Notes:
- Minimal failing paired real-data subset was `n=64` on a reduced transcriptome index derived
  from the sample itself.
- Root cause: paired pseudoalignment dropped the whole pair when one mate had no EC, while
  instrumented `kallisto` preserves the informative mate EC and still pseudoaligns the pair.
- Fixed in `crates/kallistors/src/pseudoalign/paired.rs`; regression test added in
  `crates/kallistors/tests/paired_parity.rs`.
- Later real-data paired mismatches came from overcrowded minimizer handling.
- Current status:
  exact on paired deterministic prefixes through `262144` pairs on `data/subsets/`
  full-file paired quant still differs by `2` pseudoaligned pairs

## Phase F: Loader performance parity with original `kallisto`
- [ ] Stop rebuilding the minimizer map as `Vec<Vec<u64>>`; switch to a search-oriented static structure closer to Bifrost `MinimizerIndex`
- [ ] Eliminate the intermediate `(idx, pos_id)` collection step in the loader; insert directly into the final minimizer structure during decode
- [ ] Match upstream unitig loading more closely: keep the linear `unitig_id / total_unitig_len / curr_unitig_len` walk and avoid any extra remapping work per minimizer
- [ ] Match upstream threaded unitig insertion more closely: parallelize writes into the final minimizer structure instead of parallel compute plus serial merge
- [ ] Keep the `km_unitigs` and special minimizer paths as direct inserts into the final minimizer structure, not into temporary vectors
- [ ] Reduce allocator churn in the loader by avoiding millions of tiny `Vec` growth operations and by using contiguous storage where possible
- [ ] Remove the large shutdown/drop cost caused by eagerly materialized per-minimizer vectors
- [ ] Add explicit stage timing for index load, minimizer reconstruction, pseudoalignment, and EM so loader regressions are visible immediately
- [ ] Benchmark full GENCODE index startup against `kallisto` after each structural change; do not accept changes that only reshuffle overhead

Acceptance test:
- On `data/gencode.v49_kallisto.idx`, `kallistors-cli quant` should spend most of its wall time in
  pseudoalignment/EM rather than loader startup, and the gap to `kallisto quant` on the same
  1024-read subset should shrink materially without using an on-disk cache.

Notes:
- Upstream `kallisto`/Bifrost loads minimizer occurrences directly into `MinimizerIndex` during
  `readBinaryIndex(...)` in `kallisto_src/ext/bifrost/src/IO.tcc` instead of first constructing a
  Rust-side `Vec<Vec<u64>>`.
- The upstream threaded path uses `hmap_min_unitigs.init_threads()`,
  `add_unitig_p(minz_rep, pos_id_unitig)`, and `release_threads()` to write unitig minimizer
  positions directly into the final structure.
- Our current `rayon` change improves scheduling, but it does not address the bigger structural gap:
  we still build, merge, and later drop a large eager `minz_positions` representation.
- Current full-file benchmark on the checked-in real dataset (`-t 32`):
  `kallisto` wall `59.21s`
  `kallistors` wall `190.76s`
