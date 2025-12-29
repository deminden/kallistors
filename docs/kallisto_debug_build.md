# Modified kallisto build and debug usage

This project vendors upstream kallisto under `kallisto_src` and includes
local debug instrumentation for EC tracing and hit dumps.

## Build (macOS arm64 + BAM support)

Requirements:
- CMake 3.x (use `/Applications/CMake.app/Contents/bin/cmake`)
- autotools (autoconf, automake, libtool)
- htslib (built from `kallisto_src/ext/htslib` when `USE_BAM=ON`)

Commands:
```sh
cd /Users/denisdemin/Code/kallistors/kallisto_src
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
