#!/usr/bin/env python3
"""Compare kallisto and kallistors index build time on the same transcript FASTA."""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


def build_once(
    tool: str,
    binary: Path,
    fasta: Path,
    threads: int,
    temp_root: Path,
) -> tuple[float, int]:
    with tempfile.TemporaryDirectory(prefix="local-index-", dir=temp_root) as tmp:
        output = Path(tmp) / "transcripts.idx"
        command = [str(binary), "index", "-i", str(output), "-t", str(threads)]
        if tool == "kallisto":
            command.append(str(fasta))
        else:
            command.extend(["--timings", str(fasta)])

        start = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        elapsed = time.perf_counter() - start
        if result.returncode != 0:
            raise RuntimeError(
                f"{tool} failed ({result.returncode})\n"
                f"command: {' '.join(command)}\n{result.stdout}\n{result.stderr}"
            )
        if not output.is_file():
            raise RuntimeError(f"{tool} did not write {output}")
        return elapsed, output.stat().st_size


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kallisto", type=Path, default=Path("target/tools/kallisto-0.52.0"))
    parser.add_argument("--kallistors", type=Path, default=Path("target/release/kallistors"))
    parser.add_argument(
        "--fasta", type=Path, default=Path("data/srr-small/gencode.v50.transcripts.fa.gz")
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument(
        "--out-prefix", type=Path, default=Path("target/benchmarks/local-real-2026-10-14/index")
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    binaries = {
        "kallisto": (repo / args.kallisto).resolve(),
        "kallistors": (repo / args.kallistors).resolve(),
    }
    fasta = (repo / args.fasta).resolve()
    if args.warmups != 1 or args.runs != 3:
        parser.error("this benchmark uses one warmup and three measured runs")
    if not fasta.is_file():
        parser.error(f"missing transcript FASTA: {fasta}")
    for tool, binary in binaries.items():
        if not binary.is_file():
            parser.error(f"missing {tool} executable: {binary}")

    temp_root = (repo / "target/benchmarks/local-real-2026-10-14/tmp").resolve()
    temp_root.mkdir(parents=True, exist_ok=True)
    raw_rows: list[dict[str, Any]] = []
    measured: dict[str, list[tuple[float, int]]] = {tool: [] for tool in binaries}

    for tool, binary in binaries.items():
        print(f"{tool} threads={args.threads}: warmup", flush=True)
        elapsed, size = build_once(tool, binary, fasta, args.threads, temp_root)
        raw_rows.append(
            {
                "tool": tool,
                "threads": args.threads,
                "phase": "warmup",
                "replicate": 0,
                "wall_seconds": f"{elapsed:.6f}",
                "index_bytes": size,
            }
        )

    for replicate in range(1, args.runs + 1):
        order = list(binaries)
        if replicate % 2 == 0:
            order.reverse()
        for tool in order:
            elapsed, size = build_once(tool, binaries[tool], fasta, args.threads, temp_root)
            measured[tool].append((elapsed, size))
            raw_rows.append(
                {
                    "tool": tool,
                    "threads": args.threads,
                    "phase": "measured",
                    "replicate": replicate,
                    "wall_seconds": f"{elapsed:.6f}",
                    "index_bytes": size,
                }
            )
        print(
            f"replicate {replicate}/{args.runs}: "
            f"kallisto={measured['kallisto'][-1][0]:.2f}s, "
            f"kallistors={measured['kallistors'][-1][0]:.2f}s",
            flush=True,
        )

    medians = {
        tool: statistics.median(run[0] for run in tool_runs)
        for tool, tool_runs in measured.items()
    }
    median_sizes = {
        tool: int(statistics.median(run[1] for run in tool_runs))
        for tool, tool_runs in measured.items()
    }
    summary_rows = [
        {
            "threads": args.threads,
            "runs": args.runs,
            "warmups_per_tool": args.warmups,
            "kallisto_median_seconds": f"{medians['kallisto']:.6f}",
            "kallistors_median_seconds": f"{medians['kallistors']:.6f}",
            "kallisto_over_kallistors_speedup": f"{medians['kallisto'] / medians['kallistors']:.4f}",
            "kallisto_index_bytes_median": median_sizes["kallisto"],
            "kallistors_index_bytes_median": median_sizes["kallistors"],
        }
    ]
    write_tsv(
        args.out_prefix.with_suffix(".runs.tsv"),
        raw_rows,
        ["tool", "threads", "phase", "replicate", "wall_seconds", "index_bytes"],
    )
    write_tsv(
        args.out_prefix.with_suffix(".summary.tsv"),
        summary_rows,
        [
            "threads",
            "runs",
            "warmups_per_tool",
            "kallisto_median_seconds",
            "kallistors_median_seconds",
            "kallisto_over_kallistors_speedup",
            "kallisto_index_bytes_median",
            "kallistors_index_bytes_median",
        ],
    )
    print(f"wrote {args.out_prefix}.runs.tsv and {args.out_prefix}.summary.tsv")


if __name__ == "__main__":
    main()
