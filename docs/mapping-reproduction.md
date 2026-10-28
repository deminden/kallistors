# Portable mapping commands

These templates use configurable paths for the native commands in the
[Fulgor/Themisto comparison](mapping-comparison-results.md). The
[original argument arrays](mapping-comparison-commands.json) remain the record of the
published runs. New runs produce new observations; keep their inputs, tool versions and
build settings with their results.

Use a POSIX shell on Linux or macOS. Build the pinned Fulgor and Themisto versions listed in
the report and prepare the reference as described under
[reference semantics](mapping-comparison-results.md#reference-and-algorithm-conventions).
The tools require their own index formats. For transcript-level colors, Fulgor needs one
transcript per FASTA file, while Themisto uses one color per sequence in a combined FASTA.
Preserve duplicate sequences, short transcripts, target IDs and input order. Apply the
documented ambiguity normalization and trailing-polyA clipping before building competitor
indexes; an arbitrary transcript FASTA does not reproduce the published target semantics.

## Local paths

Set these paths to existing tools, prepared inputs and indexes. Choose a fresh working
directory for new outputs. The read file is one bulk mate or the prepared PBMC cDNA reads;
mates are queried independently in this mapping-only comparison.

```bash
KALLISTORS_BIN=/path/to/kallistors
FULGOR_BIN=/path/to/fulgor
THEMISTO_BIN=/path/to/themisto
KALLISTORS_INDEX=/path/to/transcripts.kallisto.idx
COMPETITOR_FASTA=/path/to/transcripts-clipped.fa
FULGOR_FILE_LIST=/path/to/filenames.txt
READS_FASTQ=/path/to/reads.fastq.gz
MAPPING_WORK=/path/to/new-mapping-workspace
mkdir -p "$MAPPING_WORK/fulgor/temporary" "$MAPPING_WORK/themisto/temporary"
```

`FULGOR_FILE_LIST` contains the local path of each per-transcript FASTA on its own line.
Rewrite those entries when relocating prepared inputs. After building, recover the actual
Fulgor file/color order using its `print-filenames` command and map colors back to transcript
IDs; do not assume the builder preserved list order. Themisto colors follow FASTA sequence
order. The original-built kallisto v13 index is reused by kallistors.

## Index construction

Build each competitor index once before read measurements. These commands use eight requested
threads and the measured index settings; record setup resources separately from query resources.

```bash
"$FULGOR_BIN" build \
  -l "$FULGOR_FILE_LIST" -o "$MAPPING_WORK/fulgor/transcripts" \
  -k 31 -m 19 -d "$MAPPING_WORK/fulgor/temporary" -g 8 -t 8

"$THEMISTO_BIN" build \
  -k 31 -i "$COMPETITOR_FASTA" \
  --index-prefix "$MAPPING_WORK/themisto/transcripts" \
  --temp-dir "$MAPPING_WORK/themisto/temporary" \
  --mem-gigas 8 --n-threads 8 --sequence-colors
```

The published Fulgor mature-plus-nascent build failed; there is no completed Fulgor nascent
index or mapping result to reproduce. The successful nascent comparison uses Themisto and
kallistors. The eight-GiB builder memory argument is a resource hint, not a peak-RSS limit.

## Query commands

Choose `MAPPING_THREADS=1` or `8`. Keep output at `/dev/null` for the native timing commands:

```bash
MAPPING_THREADS=8
MAPPING_OUTPUT=/dev/null

"$KALLISTORS_BIN" pseudoalign \
  --index "$KALLISTORS_INDEX" --reads "$READS_FASTQ" \
  --out "$MAPPING_OUTPUT" --threads "$MAPPING_THREADS"

"$FULGOR_BIN" pseudoalign \
  -i "$MAPPING_WORK/fulgor/transcripts.fur" -q "$READS_FASTQ" \
  -o "$MAPPING_OUTPUT" -t "$MAPPING_THREADS" --format ascii --verbose

"$THEMISTO_BIN" pseudoalign \
  --index-prefix "$MAPPING_WORK/themisto/transcripts" --query-file "$READS_FASTQ" \
  --out-file "$MAPPING_OUTPUT" --temp-dir "$MAPPING_WORK/themisto/temporary" \
  --n-threads "$MAPPING_THREADS" --threshold 1
```

All three commands are unstranded. Fulgor's omitted `-r` selects full intersection;
`-r 1` would choose a different path. Themisto's threshold of one uses its intersection
path and excludes unknown k-mers. kallistors uses its native pseudoalignment semantics.
Fulgor/Themisto still enumerate and format per-read target sets to the null sink; kallistors
aggregates EC counts. Keep those output costs explicit when comparing elapsed time.

## Measurement and agreement

For each reference, read file, thread count and competitor, run one warmup per tool and
three measured pairs with a fresh kallistors control. Alternate control/competitor order
between repetitions and measure one command at a time. Retain every completed observation,
including background-contention flags. Record wall time, user/system CPU, peak RSS, actual
worker/parser counts, command arguments and executable/input hashes. Use whole-command
medians; the templates above do not collect resource measurements automatically.

Run separate full-output correctness commands using a different output filename for each
tool, replacing `MAPPING_OUTPUT` for that command. Preserve the numeric read-ID convention
and check complete read coverage. Canonicalize each tool's transcript IDs before comparing
weighted EC counts; Fulgor/Themisto per-read sets and kallistors aggregate ECs are different
output formats. A matched processed-read total alone is insufficient. See the report's
[agreement definitions](mapping-comparison-results.md#complete-mapping-agreement) and the
[benchmark guide](benchmarks.md) for interpretation and preparation.
