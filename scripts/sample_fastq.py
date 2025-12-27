#!/usr/bin/env python3
"""Deterministically sample reads from FASTQ (optionally paired)."""

import argparse
import gzip
import random
import sys
from typing import IO, Optional, Tuple


def open_fastq(path: str, mode: str) -> IO[str]:
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def read_record(handle: IO[str]) -> Optional[Tuple[str, str, str, str]]:
    lines = [handle.readline() for _ in range(4)]
    if lines[0] == "":
        return None
    if any(line == "" for line in lines):
        raise ValueError("truncated FASTQ record")
    return (lines[0], lines[1], lines[2], lines[3])


def write_record(handle: IO[str], record: Tuple[str, str, str, str]) -> None:
    handle.write(record[0])
    handle.write(record[1])
    handle.write(record[2])
    handle.write(record[3])


def sample_single(input_path: str, output_path: str, fraction: float, rng: random.Random) -> None:
    total = 0
    kept = 0
    with open_fastq(input_path, "rt") as reader, open_fastq(output_path, "wt") as writer:
        while True:
            record = read_record(reader)
            if record is None:
                break
            total += 1
            if rng.random() < fraction:
                write_record(writer, record)
                kept += 1
    print(f"kept {kept} of {total} reads ({kept / max(total, 1):.4f})")


def sample_paired(
    input1: str,
    input2: str,
    output1: str,
    output2: str,
    fraction: float,
    rng: random.Random,
) -> None:
    total = 0
    kept = 0
    with (
        open_fastq(input1, "rt") as reader1,
        open_fastq(input2, "rt") as reader2,
        open_fastq(output1, "wt") as writer1,
        open_fastq(output2, "wt") as writer2,
    ):
        while True:
            rec1 = read_record(reader1)
            rec2 = read_record(reader2)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError("paired FASTQ files have different lengths")
            total += 1
            if rng.random() < fraction:
                write_record(writer1, rec1)
                write_record(writer2, rec2)
                kept += 1
    print(f"kept {kept} of {total} read pairs ({kept / max(total, 1):.4f})")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input1", required=True, help="FASTQ(.gz) for read 1")
    parser.add_argument("--output1", required=True, help="Output FASTQ(.gz) for read 1")
    parser.add_argument("--input2", help="FASTQ(.gz) for read 2")
    parser.add_argument("--output2", help="Output FASTQ(.gz) for read 2")
    parser.add_argument("--fraction", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.fraction <= 0.0 or args.fraction > 1.0:
        print("--fraction must be in (0, 1]", file=sys.stderr)
        sys.exit(2)

    rng = random.Random(args.seed)
    if args.input2 or args.output2:
        if not args.input2 or not args.output2:
            print("--input2 and --output2 must be provided together", file=sys.stderr)
            sys.exit(2)
        sample_paired(args.input1, args.input2, args.output1, args.output2, args.fraction, rng)
    else:
        sample_single(args.input1, args.output1, args.fraction, rng)


if __name__ == "__main__":
    main()
