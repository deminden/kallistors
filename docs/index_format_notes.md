# kallisto index format notes (from kallisto_src)

Sources consulted:
- `kallisto_src/src/KmerIndex.cpp`
- `kallisto_src/src/KmerIndex.h`
- `kallisto_src/ext/bifrost/src/IO.tcc`
- `kallisto_src/ext/bifrost/src/Kmer.hpp`

## High-level layout (index built via `kallisto index`)

The index is written in two stages:
1) `KmerIndex::BuildDeBruijnGraph(...)` writes the header + dBG + MPHFs.
2) `KmerIndex::write(std::ofstream&)` appends node data and transcript metadata.

Observed binary layout (little-endian on x86_64):

1. `size_t INDEX_VERSION`
2. `size_t dbg_size` (bytes for `dbg.writeBinary` output)
3. `dbg_size` bytes: `CompactedDBG::writeBinary` (graph + index)
4. `size_t mphf_size`
5. `mphf_size` bytes: boophf MPHFs (`mphf->save`)
6. `uint64_t dlist_size`
7. `uint64_t dlist_overhang`
8. `dlist_size` entries of `Kmer` (binary)
9. `size_t node_count`
10. For each node:
    - head k-mer string, length `k` bytes (ASCII)
    - `uint32_t node_size`
    - `node_size` bytes: serialized `Node`
11. `int num_trans`
12. `num_trans` x `int` transcript lengths
13. `num_trans` x (`size_t name_len` + `name_len` bytes) transcript names
14. `size_t onlist_size`
15. `onlist_size` bytes: Roaring bitmap (`onlist_sequences.write`)

## Extracting k-mer length

`CompactedDBG::readBinaryGraph` reads a small graph header:
- `size_t file_format_version`
- `int rk` (k-mer length)
- `int rg` (minimizer length)

This header is at the beginning of the `dbg_size` bytes.

## Kmer binary size

`Kmer` uses `MAX_KMER_SIZE` (default 32). Storage is `MAX_K/4` bytes, so:
- default `Kmer` binary size = 32/4 = 8 bytes

If kallisto is compiled with a different `MAX_KMER_SIZE`, the d-list entry size will differ.

## Notes / assumptions for kallistors

- We currently assume:
  - `size_t` = 8 bytes
  - little-endian host
  - `MAX_KMER_SIZE` = 32 (8-byte `Kmer`)
- The MVP parser only needs `k`, transcript names, and lengths.
- Full node/EC table parsing is deferred.
