#!/usr/bin/env python3
"""Compare paired and single-end quantification across thread counts."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


def run_quant(
    tool: str,
    binary: Path,
    index: Path,
    reads: tuple[Path, Path],
    mode: str,
    threads: int,
    temp_root: Path,
    save_to: Path | None = None,
) -> tuple[float, dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix="local-quant-", dir=temp_root) as tmp:
        out_dir = Path(tmp) / "quant"
        command = [str(binary), "quant", "-i", str(index), "-o", str(out_dir), "-t", str(threads)]
        if mode == "single":
            command.extend(["--single", "-l", "200", "-s", "20", str(reads[0])])
        else:
            command.extend([str(reads[0]), str(reads[1])])

        start = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        elapsed = time.perf_counter() - start
        if result.returncode != 0:
            raise RuntimeError(
                f"{tool} failed ({result.returncode})\n"
                f"command: {' '.join(command)}\n{result.stdout}\n{result.stderr}"
            )
        info_path = out_dir / "run_info.json"
        if not info_path.is_file():
            raise RuntimeError(f"{tool} did not write {info_path}")
        info = json.loads(info_path.read_text(encoding="utf-8"))
        if save_to is not None:
            shutil.copytree(out_dir, save_to, dirs_exist_ok=True)
            (save_to / "command.json").write_text(json.dumps(command, indent=2) + "\n")
            (save_to / "stderr.log").write_text(result.stderr)
        return elapsed, info


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
    parser.add_argument("--index", type=Path, default=Path("data/srr-small/gencode.v50.kallisto.idx"))
    parser.add_argument("--thread-counts", default="1,2,4,8")
    parser.add_argument(
        "--modes", default="paired,single", help="comma-separated modes: paired,single (default: both)"
    )
    parser.add_argument(
        "--samples",
        default="all",
        help="comma-separated sample keys, or 'all' (default): "
        "ERR2260014_compact,DRR023395_larger,SRR13638690_larger",
    )
    parser.add_argument(
        "--compact-thread-counts",
        help="optional thread-count override for ERR2260014_compact",
    )
    parser.add_argument(
        "--larger-thread-counts",
        help="optional thread-count override for DRR023395_larger",
    )
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument(
        "--keep-last-output",
        action="store_true",
        help="retain each tool's final measured quant output for abundance comparison (copying is untimed)",
    )
    parser.add_argument(
        "--out-prefix", type=Path, default=Path("target/benchmarks/real-reads-2026-10-14/quant")
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    binaries = {
        "kallisto": (repo / args.kallisto).resolve(),
        "kallistors": (repo / args.kallistors).resolve(),
    }
    index = (repo / args.index).resolve()
    samples = {
        "ERR2260014_compact": (
            repo / "data/srr-small/trimmed/ERR2260014_1.trimmed.fastq.gz",
            repo / "data/srr-small/trimmed/ERR2260014_2.trimmed.fastq.gz",
        ),
        "DRR023395_larger": (
            repo / "data/drr023395/trimmed/DRR023395_1.trimmed.fastq.gz",
            repo / "data/drr023395/trimmed/DRR023395_2.trimmed.fastq.gz",
        ),
        "SRR13638690_larger": (
            repo / "data/srr13638690/trimmed/SRR13638690_1.trimmed.fastq.gz",
            repo / "data/srr13638690/trimmed/SRR13638690_2.trimmed.fastq.gz",
        ),
    }
    if args.samples == "all":
        selected_samples = list(samples)
    else:
        selected_samples = [sample.strip() for sample in args.samples.split(",") if sample.strip()]
        unknown = sorted(set(selected_samples) - set(samples))
        if unknown:
            parser.error(f"unknown sample key(s): {', '.join(unknown)}")
        if not selected_samples:
            parser.error("--samples must name at least one sample or be 'all'")
    selected_modes = [mode.strip() for mode in args.modes.split(",") if mode.strip()]
    if not selected_modes or set(selected_modes) - {"paired", "single"}:
        parser.error("--modes must name paired, single, or paired,single")
    thread_counts_by_sample = {
        "ERR2260014_compact": [
            int(value)
            for value in (args.compact_thread_counts or args.thread_counts).split(",")
        ],
        "DRR023395_larger": [
            int(value)
            for value in (args.larger_thread_counts or args.thread_counts).split(",")
        ],
        "SRR13638690_larger": [int(value) for value in args.thread_counts.split(",")],
    }
    if args.warmups != 1:
        parser.error("this benchmark records exactly one warmup per tool and configuration")
    if args.runs < 1:
        parser.error("--runs must be positive (default: three measured runs)")
    if not index.is_file():
        parser.error(f"missing shared kallisto index: {index}")
    for tool, binary in binaries.items():
        if not binary.is_file():
            parser.error(f"missing {tool} executable: {binary}")
    for sample in selected_samples:
        read_pair = samples[sample]
        for read in read_pair:
            if not read.is_file():
                parser.error(f"missing {sample} trimmed FASTQ: {read}; run scripts/prepare_benchmark.sh")

    temp_root = (repo / "target/benchmarks/real-reads-2026-10-14/tmp").resolve()
    temp_root.mkdir(parents=True, exist_ok=True)
    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    def record_run(row: dict[str, Any]) -> None:
        # Persist outside the timed region so completed work survives interruption.
        raw_rows.append(row)
        write_tsv(args.out_prefix.with_suffix(".runs.tsv"), raw_rows, list(row))

    for sample in selected_samples:
        read_pair = samples[sample]
        for mode in selected_modes:
            for threads in thread_counts_by_sample[sample]:
                print(f"{sample} {mode} threads={threads}: warmup", flush=True)
                for tool, binary in binaries.items():
                    elapsed, info = run_quant(
                        tool, binary, index, read_pair, mode, threads, temp_root
                    )
                    print(
                        f"  {tool} warmup counts: "
                        f"processed={info.get('n_processed')}, "
                        f"pseudoaligned={info.get('n_pseudoaligned')}, "
                        f"unique={info.get('n_unique')}",
                        flush=True,
                    )
                    record_run(
                        {
                            "sample": sample,
                            "mode": mode,
                            "threads": threads,
                            "tool": tool,
                            "phase": "warmup",
                            "replicate": 0,
                            "wall_seconds": f"{elapsed:.6f}",
                            "n_processed": info.get("n_processed"),
                            "n_pseudoaligned": info.get("n_pseudoaligned"),
                            "n_unique": info.get("n_unique"),
                        }
                    )

                measured: dict[str, list[tuple[float, dict[str, Any]]]] = {
                    tool: [] for tool in binaries
                }
                for replicate in range(1, args.runs + 1):
                    order = list(binaries)
                    if replicate % 2 == 0:
                        order.reverse()
                    for tool in order:
                        save_to = None
                        if args.keep_last_output and replicate == args.runs:
                            save_to = (
                                args.out_prefix.parent
                                / (args.out_prefix.name + ".outputs")
                                / sample
                                / mode
                                / f"threads{threads}"
                                / tool
                            )
                        elapsed, info = run_quant(
                            tool, binaries[tool], index, read_pair, mode, threads, temp_root, save_to
                        )
                        measured[tool].append((elapsed, info))
                        record_run(
                            {
                                "sample": sample,
                                "mode": mode,
                                "threads": threads,
                                "tool": tool,
                                "phase": "measured",
                                "replicate": replicate,
                                "wall_seconds": f"{elapsed:.6f}",
                                "n_processed": info.get("n_processed"),
                                "n_pseudoaligned": info.get("n_pseudoaligned"),
                                "n_unique": info.get("n_unique"),
                            }
                        )
                    print(
                        f"  replicate {replicate}/{args.runs}: "
                        f"kallisto={measured['kallisto'][-1][0]:.3f}s, "
                        f"kallistors={measured['kallistors'][-1][0]:.3f}s",
                        flush=True,
                    )

                medians = {
                    tool: statistics.median(run[0] for run in tool_runs)
                    for tool, tool_runs in measured.items()
                }
                infos = {
                    tool: measured[tool][-1][1]
                    for tool in binaries
                }
                fields = ("n_processed", "n_pseudoaligned", "n_unique")
                parity = all(
                    infos["kallisto"].get(field) == infos["kallistors"].get(field)
                    for field in fields
                )
                summary_rows.append(
                    {
                        "sample": sample,
                        "mode": mode,
                        "threads": threads,
                        "runs": args.runs,
                        "warmups_per_tool": args.warmups,
                        "kallisto_median_seconds": f"{medians['kallisto']:.6f}",
                        "kallistors_median_seconds": f"{medians['kallistors']:.6f}",
                        "kallisto_over_kallistors_speedup": f"{medians['kallisto'] / medians['kallistors']:.4f}",
                        "run_info_counts_match": parity,
                        "n_processed": infos["kallisto"].get("n_processed"),
                        "n_pseudoaligned": infos["kallisto"].get("n_pseudoaligned"),
                        "n_unique": infos["kallisto"].get("n_unique"),
                    }
                )
                write_tsv(
                    args.out_prefix.with_suffix(".summary.tsv"),
                    summary_rows,
                    list(summary_rows[-1]),
                )

    raw_fields = [
        "sample",
        "mode",
        "threads",
        "tool",
        "phase",
        "replicate",
        "wall_seconds",
        "n_processed",
        "n_pseudoaligned",
        "n_unique",
    ]
    summary_fields = [
        "sample",
        "mode",
        "threads",
        "runs",
        "warmups_per_tool",
        "kallisto_median_seconds",
        "kallistors_median_seconds",
        "kallisto_over_kallistors_speedup",
        "run_info_counts_match",
        "n_processed",
        "n_pseudoaligned",
        "n_unique",
    ]
    write_tsv(args.out_prefix.with_suffix(".runs.tsv"), raw_rows, raw_fields)
    write_tsv(args.out_prefix.with_suffix(".summary.tsv"), summary_rows, summary_fields)
    print(f"wrote {args.out_prefix}.runs.tsv and {args.out_prefix}.summary.tsv")


if __name__ == "__main__":
    main()
