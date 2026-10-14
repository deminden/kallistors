# Modified kallisto build and debug usage

Two diagnostic workflows are documented here. The October 14 v0.52.0 quantification audit uses
the tracked patch below. The older `kallisto_src` workflow later on this page uses a separate local
patch whose CLI tracing flags are not in the tracked audit patch or stock upstream kallisto.

## Reproducible v0.52.0 quantification audit build

The tracked [quant audit patch](../scripts/diagnostics/kallisto_v0.52.0_quant_audit.patch) applies
to upstream tag v0.52.0 (`4e9f29cf3b021260415430c057a22469ca081391`). It adds EC counts,
fragment histograms, per-read diagnostics, and 17-digit abundance output, plus CMake/Apple Silicon
build compatibility changes. It was used for the [mismatch audit](quantification_mismatch_audit_2026-10-14.md).

From the repository root, with Git and CMake available:

```bash
git clone --branch v0.52.0 --depth 1 https://github.com/pachterlab/kallisto.git \
  target/tools/kallisto-src
git -C target/tools/kallisto-src apply \
  "$PWD/scripts/diagnostics/kallisto_v0.52.0_quant_audit.patch"
cmake -S target/tools/kallisto-src -B target/tools/kallisto-src/build \
  -DUSE_BAM=OFF -DUSE_HDF5=OFF -DENABLE_AVX2=OFF -DCOMPILATION_ARCH=OFF
cmake --build target/tools/kallisto-src/build --parallel 8

KALLISTO_DIAGNOSTIC_ECS=1 target/tools/kallisto-src/build/src/kallisto quant \
  -i data/srr-small/gencode.v50.kallisto.idx -t 8 \
  -o target/validation/upstream-ec \
  data/srr13638690/trimmed/SRR13638690_1.trimmed.fastq.gz \
  data/srr13638690/trimmed/SRR13638690_2.trimmed.fastq.gz
```

This emits `diagnostic.ec_counts.tsv` and `diagnostic.fragment_lengths.tsv` in the quant output
directory. Optional environment controls:

| Variable | Output |
| --- | --- |
| `KALLISTO_DIAGNOSTIC_HITS=<path>` | Per-read EC/hit-unitig and position trace |
| `KALLISTO_DIAGNOSTIC_READS=<path>` | Singleton read IDs and sequences; also set the TID below |
| `KALLISTO_DIAGNOSTIC_TID=96635` | Select TTN-206 in the documented GENCODE v50 index |

Numeric transcript and unitig IDs depend on the index. The no-HDF5/no-BAM diagnostic binary is
for correctness investigation; use the official release executable for benchmark timing.

## Historical local tracing build

## Build (macOS arm64 + BAM support)

Requirements:
- CMake 3.x (use `/Applications/CMake.app/Contents/bin/cmake`)
- autotools (autoconf, automake, libtool)
- htslib (built from `kallisto_src/ext/htslib` when `USE_BAM=ON`)

Commands:
```sh
cd kallisto_src
rm -rf build
mkdir build
cd build
/Applications/CMake.app/Contents/bin/cmake .. -DUSE_BAM=ON -DENABLE_AVX2=OFF -DCOMPILATION_ARCH=OFF
make
```

Binary output:
`kallisto_src/build/src/kallisto`

## Debug outputs added to kallisto

These flags are added to the `kallisto quant` and `kallisto inspect` CLIs.

### Per-read EC trace (quant)
Flags:
- `--ec-trace=FILE` write per-read trace rows
- `--ec-trace-max-reads=INT` limit to first N reads (0 = all)
- `--hit-dump=FILE` unitig/block hit dump per read

Example:
```sh
./kallisto quant \
  -i INDEX \
  -o OUTDIR \
  --ec-trace OUTDIR/ec_trace.tsv \
  --ec-trace-max-reads 1000 \
  --hit-dump OUTDIR/hits.tsv \
  READS_1.fq READS_2.fq
```

`ec_trace.tsv` columns:
```
type	read_id	read_name1	read_name2	part	hit_idx	unitig_id	unitig_pos	read_pos	block_id	block_lb	block_ub	ec	running_ec	action
```
Row types:
- `READ` metadata about the read and hit counts
- `HIT` per-hit candidate EC + running intersection state
- `FINAL` final EC or empty reason

`hits.tsv` columns:
```
read_id	read_name	part	hit_idx	unitig_id	unitig_pos	read_pos	block_id	block_lb	block_ub
```

### Index-level k-mer -> EC dump (inspect)
Flags:
- `--ec-dump=FILE` write k-mer/EC mapping dump
- `--ec-dump-limit=INT` maximum k-mers to dump (default 1000)

Example:
```sh
./kallisto inspect INDEX --ec-dump index_ec_dump.tsv --ec-dump-limit 500
```

`index_ec_dump.tsv` columns:
```
unitig_id	unitig_pos	kmer	block_id	block_lb	block_ub	ec
```
