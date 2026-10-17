# Instrumented kallisto for correctness checks

Use the official kallisto v0.52.0 release executable as the performance control. A separate
instrumented build supplies complete EC counts and targeted traces for correctness work.
The patches target upstream commit `4e9f29cf3b021260415430c057a22469ca081391` (tag v0.52.0).

## Build

The [build compatibility patch](../scripts/diagnostics/kallisto_v0.52.0_build_compat.patch)
contains CMake/Apple Silicon adjustments and the Bifrost copy-constructor build fix.
The [diagnostics patch](../scripts/diagnostics/kallisto_v0.52.0_diagnostics.patch) provides EC
exports, mapping-only execution, selected reads, fragment-sampling traces and 17-digit TSVs.

```bash
git clone --branch v0.52.0 --depth 1 https://github.com/pachterlab/kallisto.git \
  target/tools/kallisto-src
git -C target/tools/kallisto-src apply \
  "$PWD/scripts/diagnostics/kallisto_v0.52.0_build_compat.patch"
git -C target/tools/kallisto-src apply \
  "$PWD/scripts/diagnostics/kallisto_v0.52.0_diagnostics.patch"
cmake -S target/tools/kallisto-src -B target/tools/kallisto-src/build \
  -DUSE_BAM=OFF -DUSE_HDF5=OFF -DENABLE_AVX2=OFF -DCOMPILATION_ARCH=OFF
cmake --build target/tools/kallisto-src/build --parallel 8
```

Use a fresh checkout for these commands. The diagnostic executable is passed separately as
`--diagnostic-kallisto` to the [benchmark runner](benchmarks.md); the runner rejects using the
same executable bytes for diagnostic work and the timing control.

## Complete EC exports

```bash
KALLISTO_DIAGNOSTIC_ECS=1 KALLISTO_DIAGNOSTIC_MAPPING_ONLY=1 \
  target/tools/kallisto-src/build/src/kallisto quant \
  -i data/srr-small/gencode.v50.kallisto.idx -t 8 \
  -o target/validation/upstream-ec \
  data/srr13638690/trimmed-fasterp/SRR13638690_1.trimmed.fastq.gz \
  data/srr13638690/trimmed-fasterp/SRR13638690_2.trimmed.fastq.gz
```

This writes `diagnostic.ec_counts.tsv` and exits before fragment-model fitting and EM.
Omit mapping-only mode to run full quantification and also export
`diagnostic.fragment_lengths.tsv`. Numeric transcript/unitig IDs depend on the index;
compare counted transcript sets using identical index bytes.

## Targeted diagnostics

| Variable | Behavior |
| --- | --- |
| `KALLISTO_DIAGNOSTIC_HITS=PATH` | Read ECs and hit-unitig/position traces |
| `KALLISTO_DIAGNOSTIC_READS=PATH` | Singleton read IDs/sequences; also set `KALLISTO_DIAGNOSTIC_TID=ID` |
| `KALLISTO_DIAGNOSTIC_SELECTED_PAIRS=PATH` | Read name, both sequences and final EC |
| `KALLISTO_DIAGNOSTIC_EC_FILTER='0,1;2,3'` | Select counted transcript sets for paired-read export |
| `KALLISTO_DIAGNOSTIC_SAMPLING=PATH` | Native buffer quotas and accepted pair ordinals/fragment lengths |
| `KALLISTO_DIAGNOSTIC_STOP_AFTER_SAMPLING=1` | Stop after the sampling goal, for prefix investigation |

Selected-pair export requests names and synchronizes writes across workers. Sampling traces
retain ordinal information without requesting names, preserving sequence-buffer boundaries.
The sampling-prefix shortcut can be combined with mapping-only mode. These diagnostic
shortcuts are excluded from timing; the runner clears diagnostic environment controls in
measured processes.

## Fixed-EC/model replay

`scripts/diagnostics/prepare_quant_replay.py` reads complete EC exports and native HDF5
histograms to prepare all transcript lengths, ECs and fragment bins. The `replay_ec_quant`
Rust example runs the existing quantifier without rereading FASTQs or an index:

```bash
python scripts/diagnostics/prepare_quant_replay.py \
  --ecs path/to/diagnostic.ec_counts.tsv \
  --original path/to/original/abundance.h5 --candidate path/to/candidate/abundance.h5 \
  --out target/validation/replay
cargo run --locked --release --example replay_ec_quant -- \
  target/validation/replay/ecs.tsv target/validation/replay/lengths \
  target/validation/replay/original.histogram target/validation/replay/quant.tsv 8 default
```

Use a Python interpreter with h5py and NumPy for replay preparation. Target-name filtering
selects metadata for inspection; it does not remove transcripts or ECs from the replay.
Matched ECs, models and iteration settings separate EM behavior from mapping and sampling.
