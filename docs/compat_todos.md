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
- [x] Stop rebuilding the minimizer map as `Vec<Vec<u64>>`; switch to a flat contiguous minimizer store
- [x] Eliminate the intermediate `(idx, pos_id)` collection step in the loader; write into the final minimizer storage during the fill pass
- [ ] Match upstream unitig loading more closely: keep the linear `unitig_id / total_unitig_len / curr_unitig_len` walk and avoid any extra remapping work per minimizer
- [x] Match upstream threaded unitig insertion more closely: parallelize writes into the final minimizer structure instead of parallel compute plus serial merge
- [ ] Keep the `km_unitigs` and special minimizer paths as direct inserts into the final minimizer structure, not into temporary vectors
- [x] Reduce allocator churn in the loader by avoiding millions of tiny `Vec` growth operations and by using contiguous storage where possible
- [x] Remove the large shutdown/drop cost caused by eagerly materialized per-minimizer vectors
- [x] Add explicit stage timing for index load, minimizer reconstruction, pseudoalignment, and EM so loader regressions are visible immediately
- [ ] Benchmark full GENCODE index startup against `kallisto` after each structural change; do not accept changes that only reshuffle overhead
- [x] Replace per-occurrence atomic minimizer counting with buffered local batch counting that mirrors upstream `readBinaryIndex(...)` batching while preserving the existing exact storage write path

Acceptance test:
- On `data/gencode.v49_kallisto.idx`, `kallistors-cli quant` should spend most of its wall time in
  pseudoalignment/EM rather than loader startup, and the gap to `kallisto quant` on the same
  1024-read subset should shrink materially without using an on-disk cache.

Notes:
- Upstream `kallisto`/Bifrost loads minimizer occurrences directly into `MinimizerIndex` during
  `readBinaryIndex(...)` in `kallisto_src/ext/bifrost/src/IO.tcc` instead of first constructing a
  Rust-side vector-of-vectors representation.
- The upstream threaded path uses `hmap_min_unitigs.init_threads()`,
  `add_unitig_p(minz_rep, pos_id_unitig)`, and `release_threads()` to write unitig minimizer
  positions directly into the final structure.
- Current common-case search optimizations mirror upstream more closely:
  direct-mapped MPHF minimizer caching, reusable minimizer candidate buffers, and a fast path that
  stays out of shade/debug/fallback logic.
- Common-case candidate validation now uses precomputed encoded unitig / km-unitig base stores
  rather than repeated byte-slice equality in the hot path.
- Pathological and repeated k-mer reuse now extend beyond a single read via a thread-local fast
  match cache in the worker hot path. That cache memoizes exact common-path and overcrowded/special
  fallback outcomes for recurring k-mers without changing default behavior.
- EC block lookup now uses a small thread-local block-range cache, and fast-path minimizer
  enumeration uses a thread-local minimizer-candidate cache to avoid redoing the same
  `minimizers_for_kmer` work for recurring windows.
- Pathological selective fallback now stores pre-encoded forward k-mer codes, so fallback matching
  no longer re-encodes the candidate unitig sequence for each overcrowded/special hit.
- Tried a more original-shaped fast matcher over a flatter typed minimizer bucket index
  (regular / short / special flags split out of raw `pos_id` words). It preserved parity on the
  real paired ladder through `4096`, but it regressed the full-file benchmark (`131.31s` wall), so
  it was reverted in favor of the simpler cached matcher.
- The current loader win comes from mirroring the upstream static-MPHF read path more closely on
  the count side: unitig/km minimizer positions are buffered into large ordered batches, each batch
  is processed in parallel into thread-local sorted minimizer-id runs, and those local runs are
  merged into the global counts without a relaxed atomic increment on every occurrence.
- Tried extending that same batched/local approach to the fill/write pass. It was much faster, but
  it changed large-subset/full-file pseudoalignment counts, so it was reverted. The exact write
  path still uses the previous storage fill logic, while only the count pass takes the new batched
  route.
- Tried a quant-only loader shortcut that omitted nested `ec_blocks` and kept only the flattened EC
  representation. It looked fine on tiny subsets but regressed larger paired parity, so normal
  quant/pseudoalign commands were switched back to the standard builder.
- Remaining search-side gap is still in candidate validation and pathological fallback handling,
  not in loader timing or EC merge.
- Current full-file benchmark on the checked-in real dataset (`-t 32`):
  `kallisto` wall `56.61s`
  `kallistors` wall `81.85s`
  `kallistors` stage timings:
  `index_header_parse 2.113s`
  `graph_decode 6.859s`
  `minimizer_count_pass 1.390s`
  `minimizer_fill_pass 9.002s`
  `fastq_read_decompress 34.613s`
  `pseudoalign 45.433s`
  `em 10.100s`

## Phase G: Mirror remaining upstream search-loop optimizations
- [ ] Replace per-k-mer minimizer recomputation with a rolling `minHashIterator`-style state for the common quant path
- [ ] Preserve iterator-local minimizer state across adjacent read windows so the next k-mer search can reuse the current minimizer set instead of rebuilding it from slices
- [ ] Mirror upstream overcrowded-minimizer handling: when a bucket is overcrowded, advance to the next minimizer candidate first and only fall back to direct search if no alternate minimizer proves the hit
- [ ] Add a fast “minimizer presence changed” gate like upstream `Search.tcc` so repeated adjacent windows do not redo MPHF lookups when the primary minimizer position/hash is unchanged
- [ ] Mirror the original quant jump loop more literally: maintain one rolling read iterator and one rolling search state, and use the next/middle-hit checks to skip windows before falling back to incremental scanning
- [ ] Mirror upstream incremental fallback more literally after a failed jump: scan only the intended checkpoint windows up to the next stop, rather than dropping immediately into a broader compatibility-heavy path
- [ ] Rework common-case bucket probing to stay as a single compact late-branch scan over raw entries, matching `packed_tiny_vector` semantics without introducing a slower typed side index
- [ ] Delay `match_kmer_direct(...)` in the common path so it behaves like a true rare fallback instead of an eager branch on every overcrowded marker
- [ ] After each search-loop change, rerun the deterministic real-data ladder through `262144` plus the full file; reject changes that improve wall time but move full-file parity away from exactness

Acceptance test:
- On the checked-in paired real dataset, the common search path should stop rebuilding minimizer
  candidates for every adjacent k-mer window, full-file `pseudoalign` time should drop materially
  from the current `83.624s`, and the deterministic paired ladder should remain exact through
  `262144` while the full-file drift does not worsen.

Strategy:
- Start with iterator-state reuse, because upstream Bifrost search is built around `minHashIterator`
  and `getNewMin(...)`, while `kallistors` still rebuilds minimizer candidates from slices in the
  hot loop. The first implementation should introduce a rolling minimizer state object for one read
  window that can advance by one base and expose:
  `primary hash`, `primary position`, `alternate minimizers`, and `next minimizer after hash`.
- Then change the fast matcher to consume that rolling state directly instead of calling
  `minimizers_for_kmer_into`, `minimizers_ranked_for_kmer`, and `minhash_next_after_hash` on raw
  slices for each checked k-mer.
- After that, mirror the original overcrowded branch from `CompactedDBG.tcc`: keep scanning the
  current compact bucket, and when an overcrowded marker is hit, ask the rolling iterator for the
  next minimizer before escalating to direct fallback. This should reduce expensive fallback calls
  on repetitive reads.
- Next, rewrite the jump/backoff loop to match the `KmerIndex.cpp` control flow more closely:
  one rolling read iterator, one current match state, explicit `nextPos` / `middlePos` probes, and
  a narrow incremental scan window when the jump cannot be proven. The intent is to reduce both
  branchiness and repeated k-mer reconstruction in the common paired quant path.
- Keep the bucket representation compact. The previous typed candidate experiment was slower, so any
  further “original-like” work should preserve the current raw compact storage and only change the
  probing/control flow around it.

Notes:
- Upstream Bifrost search in `kallisto_src/ext/bifrost/src/CompactedDBG.tcc` does not immediately
  take a heavyweight fallback on overcrowded minimizers. It calls `it_min.getNewMin(mhr)` and
  retries the next minimizer first.
- Upstream search in `kallisto_src/ext/bifrost/src/Search.tcc` keeps `minz_pres` and reuses
  minimizer-presence state across adjacent windows so repeated windows do not always rescan the same
  minimizers.
- The kallisto quant loop in `kallisto_src/src/KmerIndex.cpp` keeps one rolling `KmerIterator`,
  performs explicit next/middle-hit probes, and only falls back to incremental scanning for a
  bounded region. `kallistors` has similar logic now, but it still pays more to rebuild minimizer
  state and fallback candidates inside that loop.
- The previous typed minimizer bucket experiment showed that “copy the representation” is not
  enough; the missing win is in the iterator/probing semantics of the original loop, not in adding
  another layer of candidate structs.
- Tried a naive rolling minimizer window in the fast path and an eager upstream-inspired
  overcrowded retry. It reduced full-file wall time sharply, but it also regressed full-file parity
  (`aligned 4,245,700` / `4,245,946` vs the current correct `4,244,773`), so it was reverted.
  The next attempt must mirror upstream iterator semantics more exactly instead of approximating
  them with a simple shifted minimizer cache.
- Tried swapping the fast path over to exact `minHashKmer`-style primary/next minimizer
  selection without changing the rest of the search-loop control flow. That also regressed
  full-file parity (`aligned 4,245,294`), so the remaining mismatch is not just candidate
  selection; the surrounding iterator/probing semantics must change together with it.
- Tried a wider integrated rewrite of the fast quant loop to be more `KmerIndex.cpp`-like:
  immediate bounded backoff scans after failed jump probes, plus tighter shared match-state
  updates. That improved end-to-end wall time sharply (`85.05s` total, `pseudoalign 40.489s`),
  but it regressed full-file parity even harder (`aligned 4,245,888`), so it was reverted.
  The salvage path is to reintroduce that speedup as isolated deltas only after each delta is
  proven on the deterministic paired ladder and the full file.
- Next salvage sequence for the faster state:
  1. Keep the current correct baseline loop intact as the control branch.
  2. Reapply only the bounded incremental-scan change while preserving the old hit/intersection
     bookkeeping and compare the first divergent pair against the baseline.
  3. Reapply only the faster probe-state update rules and compare again.
  4. Keep only the deltas that preserve exact full-file parity; reject any speedup that changes
     `n_pseudoaligned`, even if the wall time improvement is large.
- Salvage tooling now exists for this workflow:
  - internal investigation flags in `PseudoalignOptions` driven by environment variables
    (`KALLISTORS_TRACE_FAST_PATH`, `KALLISTORS_FAST_DELTA_*`) so experimental fast-path deltas
    can be enabled without changing default CLI behavior
  - `kallistors-cli trace-compare` to run baseline vs fast-path traces side-by-side on selected
    reads or read pairs and emit structured diffs for probe positions, minimizer candidates, jump
    decisions, and final ECs
  - `scripts/real_subset_parity.py --fast-env ...` to reproduce a mismatched subset under an
    experimental fast branch, isolate the first divergent paired read, emit a `trace-compare`
    report, and dump a `kallisto --ec-trace` for that first divergent pair
- Resolution of the remaining full-file `+1`:
  - the offending pair was `SRR13638690.2703620`
  - the root cause was not threading; it was a default matcher divergence from Bifrost minimizer
    choice on the left mate at `read_pos=98`
  - the default ranked-minimizer path chose `0000d4257d17944b` (`1110991:1459:u`) and missed
    completely, then advanced to `read_pos=99`
  - `--kallisto-bifrost-find` instead chose `0000545397f45d50`
    (`748115:102:u,1395444:95:u`), accepted unitig `748115`, and collapsed the running EC
    intersection to empty
  - the fix is a narrow Bifrost-style retry on probe/backoff misses after prior evidence exists:
    when the default matcher misses despite prior accepted hits, retry the same k-mer through
    `match_kmer_at_pos(..., kallisto_bifrost_find = true)` before advancing
  - the retry is now wired into both:
    - `ec_for_read_bifrost_fast(...)`
    - `ec_for_read_bifrost(...)`
- validation after the fix:
  - isolated pair `SRR13638690.2703620`: `1 -> 0` pseudoaligned
  - full file: exact parity restored at `4,244,771` pseudoaligned and `276,251` unique

Input-path work:
- switched `flate2` to the `zlib-rs` backend in `crates/kallistors/Cargo.toml`
  after refreshing the registry state so Cargo could resolve `zlib-rs 0.6.3`
- measured full-file paired quant with `zlib-rs` at `83.11s` wall, down from
  `84.02s` on the previous exact build
- removed one layer of threaded worker churn by accumulating directly into
  long-lived worker `EcCounts` state instead of rebuilding and merging a fresh
  `EcCounts` per batch:
  - `crates/kallistors/src/pseudoalign/threaded.rs`
  - `crates/kallistors/src/pseudoalign/single.rs`
  - `crates/kallistors/src/pseudoalign/paired.rs`
  - `crates/kallistors/src/pseudoalign/utils.rs`
- measured full-file paired quant after that change at `82.85s` wall; this is
  a small but real improvement (`-0.26s`, about `-0.3%`)
- the remaining input-path gap is now mostly the per-record FASTQ object model
  in `crates/kallistors/src/io/mod.rs`: `FastqReader::next_record()` still
  allocates four `Vec<u8>` objects per record
- replaced the threaded owned-record transport with packed/reusable FASTQ
  batches built around:
  - one contiguous backing buffer for all record bytes in a batch
  - per-record offsets and lengths instead of four owned `Vec<u8>` fields
  - direct worker consumption from lightweight record views
- this touched:
  - `crates/kallistors/src/io/mod.rs`
  - `crates/kallistors/src/pseudoalign/threaded.rs`
  - `crates/kallistors/src/pseudoalign/single.rs`
  - `crates/kallistors/src/pseudoalign/paired.rs`
- measured full-file paired quant after the packed-batch change at `68.63s`
  wall with exact real-data parity (`4,244,771` pseudoaligned, `276,251`
  unique); this is a major improvement over `82.85s` (`-14.22s`, about
  `-17.2%`)
- current remaining gap to upstream `kallisto -t 32` (`56.24s`) is now mostly
  in pseudoalignment proper rather than loader startup or per-record FASTQ
  allocation/copying
