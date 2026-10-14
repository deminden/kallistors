#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

FASTP_BIN="${FASTP_BIN:-fastp}"
THREADS="${THREADS:-4}"
SMALL_DIR="data/srr-small"
LARGE_DIR="data/drr023395"
MUSCLE_DIR="data/srr13638690"
REF="$SMALL_DIR/gencode.v50.transcripts.fa.gz"

fetch_gzip() {
  local url="$1"
  local path="$2"
  local expected_md5="${3:-}"

  if [[ ! -s "$path" ]] || ! gzip -t "$path" >/dev/null 2>&1; then
    mkdir -p "$(dirname "$path")"
    curl -fL --retry 3 "$url" -o "$path"
  fi

  gzip -t "$path"
  if [[ -n "$expected_md5" ]]; then
    python3 - "$path" "$expected_md5" <<'PY'
import hashlib
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
expected = sys.argv[2].lower()
hasher = hashlib.md5()
with path.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        hasher.update(block)
digest = hasher.hexdigest()
if digest != expected:
    raise SystemExit(f"MD5 mismatch for {path}: expected {expected}, got {digest}")
print(f"MD5 OK: {path}")
PY
  fi
}

fetch_gzip \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR226/004/ERR2260014/ERR2260014_1.fastq.gz" \
  "$SMALL_DIR/ERR2260014_1.fastq.gz" \
  "695aee1b6358b4514b488d1008f9c290"
fetch_gzip \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR226/004/ERR2260014/ERR2260014_2.fastq.gz" \
  "$SMALL_DIR/ERR2260014_2.fastq.gz" \
  "b92a4181c93bbde2228891735e0e26b0"
fetch_gzip \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR023/DRR023395/DRR023395_1.fastq.gz" \
  "$LARGE_DIR/DRR023395_1.fastq.gz" \
  "56fab394ca1e9ddcb077c7064c430d1f"
fetch_gzip \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR023/DRR023395/DRR023395_2.fastq.gz" \
  "$LARGE_DIR/DRR023395_2.fastq.gz" \
  "4ee867f3bcd04c2754a8ca65a7b83d91"
fetch_gzip \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/090/SRR13638690/SRR13638690_1.fastq.gz" \
  "$MUSCLE_DIR/SRR13638690_1.fastq.gz" \
  "b6282c8f9f618dea68cc92c887e74616"
fetch_gzip \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/090/SRR13638690/SRR13638690_2.fastq.gz" \
  "$MUSCLE_DIR/SRR13638690_2.fastq.gz" \
  "13e9ea4645d46a6c0b3a72b40afedbe5"
fetch_gzip \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz" \
  "$REF"

mkdir -p "$SMALL_DIR/trimmed" "$LARGE_DIR/trimmed" "$MUSCLE_DIR/trimmed"
"$FASTP_BIN" \
  --in1 "$SMALL_DIR/ERR2260014_1.fastq.gz" \
  --in2 "$SMALL_DIR/ERR2260014_2.fastq.gz" \
  --out1 "$SMALL_DIR/trimmed/ERR2260014_1.trimmed.fastq.gz" \
  --out2 "$SMALL_DIR/trimmed/ERR2260014_2.trimmed.fastq.gz" \
  --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
  --length_required 50 --thread "$THREADS" \
  --json "$SMALL_DIR/trimmed/fastp.json" --html "$SMALL_DIR/trimmed/fastp.html"
"$FASTP_BIN" \
  --in1 "$LARGE_DIR/DRR023395_1.fastq.gz" \
  --in2 "$LARGE_DIR/DRR023395_2.fastq.gz" \
  --out1 "$LARGE_DIR/trimmed/DRR023395_1.trimmed.fastq.gz" \
  --out2 "$LARGE_DIR/trimmed/DRR023395_2.trimmed.fastq.gz" \
  --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
  --length_required 50 --thread "$THREADS" \
  --json "$LARGE_DIR/trimmed/fastp.json" --html "$LARGE_DIR/trimmed/fastp.html"
"$FASTP_BIN" \
  --in1 "$MUSCLE_DIR/SRR13638690_1.fastq.gz" \
  --in2 "$MUSCLE_DIR/SRR13638690_2.fastq.gz" \
  --out1 "$MUSCLE_DIR/trimmed/SRR13638690_1.trimmed.fastq.gz" \
  --out2 "$MUSCLE_DIR/trimmed/SRR13638690_2.trimmed.fastq.gz" \
  --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
  --length_required 50 --thread "$THREADS" \
  --json "$MUSCLE_DIR/trimmed/fastp.json" --html "$MUSCLE_DIR/trimmed/fastp.html"

printf 'Prepared reads and reference under data/srr-small, data/drr023395, and data/srr13638690.\n'
