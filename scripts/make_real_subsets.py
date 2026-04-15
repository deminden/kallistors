#!/usr/bin/env python3
"""Create deterministic paired FASTQ subsets from the real-data inputs in data/."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterable, List, Tuple


def parse_sizes(value: str) -> List[int]:
    sizes = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        size = int(item)
        if size <= 0:
            raise argparse.ArgumentTypeError("subset sizes must be positive")
        sizes.append(size)
    if not sizes:
        raise argparse.ArgumentTypeError("at least one subset size is required")
    return sorted(set(sizes))


def infer_inputs(data_dir: Path) -> Tuple[Path, Path]:
    candidates = sorted(
        path for path in data_dir.glob("*.gz") if "gencode" not in path.name.lower()
    )
    mate1 = []
    mate2 = []
    for path in candidates:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            header = handle.readline().strip()
        if header.endswith("/1"):
            mate1.append(path)
        elif header.endswith("/2"):
            mate2.append(path)
    if not mate1 or not mate2:
        raise SystemExit("could not infer mate 1 / mate 2 FASTQ files in data/")
    return max(mate1, key=lambda path: path.stat().st_size), max(
        mate2, key=lambda path: path.stat().st_size
    )


def read_record(handle: gzip.GzipFile) -> List[str]:
    lines = []
    for _ in range(4):
        line = handle.readline()
        if not line:
            return []
        lines.append(line)
    return lines


def copy_prefix(src1: Path, src2: Path, dst1: Path, dst2: Path, n_reads: int) -> None:
    with gzip.open(src1, "rt", encoding="utf-8", errors="ignore") as in1, gzip.open(
        src2, "rt", encoding="utf-8", errors="ignore"
    ) as in2, gzip.open(dst1, "wt", encoding="utf-8") as out1, gzip.open(
        dst2, "wt", encoding="utf-8"
    ) as out2:
        for _ in range(n_reads):
            rec1 = read_record(in1)
            rec2 = read_record(in2)
            if not rec1 or not rec2:
                raise SystemExit(f"input ended before {n_reads} reads were copied")
            out1.writelines(rec1)
            out2.writelines(rec2)


def write_manifest(
    out_dir: Path,
    src1: Path,
    src2: Path,
    outputs: Iterable[Tuple[int, Path, Path]],
) -> None:
    manifest = out_dir / "manifest.tsv"
    with manifest.open("w", encoding="utf-8") as handle:
        handle.write("reads\tmate1\tmate2\tsource_mate1\tsource_mate2\n")
        for reads, mate1, mate2 in outputs:
            handle.write(
                f"{reads}\t{mate1}\t{mate2}\t{src1.name}\t{src2.name}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", help="Directory containing the real FASTQs")
    parser.add_argument(
        "--sizes",
        default="64,256,1024",
        type=parse_sizes,
        help="Comma-separated read counts for each subset",
    )
    parser.add_argument(
        "--out-dir",
        default="data/subsets",
        help="Output directory for subset FASTQs",
    )
    parser.add_argument("--mate1", help="Optional explicit mate 1 FASTQ(.gz)")
    parser.add_argument("--mate2", help="Optional explicit mate 2 FASTQ(.gz)")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.mate1 and args.mate2:
        src1 = Path(args.mate1)
        src2 = Path(args.mate2)
    else:
        src1, src2 = infer_inputs(data_dir)

    outputs = []
    for reads in args.sizes:
        dst1 = out_dir / f"real_mate1_n{reads}.fastq.gz"
        dst2 = out_dir / f"real_mate2_n{reads}.fastq.gz"
        copy_prefix(src1, src2, dst1, dst2, reads)
        outputs.append((reads, dst1, dst2))
        print(f"wrote {dst1} and {dst2}")

    write_manifest(out_dir, src1, src2, outputs)
    print(f"wrote {out_dir / 'manifest.tsv'}")


if __name__ == "__main__":
    main()
