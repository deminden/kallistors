#!/usr/bin/env python3
"""Benchmark kallisto index vs kallistors index and write a Markdown report."""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import List, Optional


def run_command(cmd: List[str]) -> tuple[subprocess.CompletedProcess[str], float]:
    start = time.perf_counter()
    result = subprocess.run(cmd, text=True, capture_output=True, check=False)
    return result, time.perf_counter() - start


def tail(text: str, lines: int = 40) -> str:
    parts = text.strip().splitlines()
    return "\n".join(parts[-lines:])


def read_index_info(kallistors_bin: str, index: Path) -> dict[str, str]:
    result, _ = run_command([kallistors_bin, "index-info", "--index", str(index)])
    if result.returncode != 0:
        return {"error": tail(result.stderr or result.stdout)}
    out: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        out[key.strip()] = value.strip()
    return out


def inspect_with_kallisto(kallisto_bin: str, index: Path) -> tuple[int, str]:
    result, _ = run_command([kallisto_bin, "inspect", str(index)])
    return result.returncode, tail((result.stdout or "") + (result.stderr or ""))


def command_available(path: str) -> bool:
    if os.path.sep in path:
        return Path(path).exists()
    return shutil.which(path) is not None


def format_size(path: Optional[Path]) -> str:
    if path is None or not path.exists():
        return "missing"
    size = path.stat().st_size
    units = ["B", "KiB", "MiB", "GiB"]
    value = float(size)
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return str(size)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", default="data/gencode.v49.transcripts.fa.gz")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--kallisto-bin", default="./kallisto_src/build/src/kallisto")
    parser.add_argument("--kallistors-bin", default="./target/release/kallistors")
    parser.add_argument("--out", default="target/validation/index_build_benchmark.md")
    parser.add_argument(
        "--work-dir",
        default="",
        help="Directory for generated indexes. Defaults to a temporary directory.",
    )
    parser.add_argument("--keep", action="store_true", help="Keep temporary indexes.")
    parser.add_argument("--skip-kallisto", action="store_true")
    parser.add_argument("--skip-kallistors", action="store_true")
    args = parser.parse_args()

    fasta = Path(args.fasta)
    if not fasta.exists():
        raise SystemExit(f"missing FASTA: {fasta}")
    if not args.skip_kallisto and not command_available(args.kallisto_bin):
        raise SystemExit(f"kallisto binary not found: {args.kallisto_bin}")
    if not args.skip_kallistors and not command_available(args.kallistors_bin):
        raise SystemExit(f"kallistors binary not found: {args.kallistors_bin}")

    temp_context: tempfile.TemporaryDirectory[str] | None = None
    if args.work_dir:
        work_dir = Path(args.work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
    elif args.keep:
        work_dir = Path(tempfile.mkdtemp(prefix="kallistors-index-bench-"))
    else:
        temp_context = tempfile.TemporaryDirectory(prefix="kallistors-index-bench-")
        work_dir = Path(temp_context.name)

    kallisto_idx = work_dir / "kallisto.idx"
    kallistors_idx = work_dir / "kallistors.idx"
    rows = []

    if not args.skip_kallisto:
        cmd = [args.kallisto_bin, "index", "-i", str(kallisto_idx), "-k", str(args.k), str(fasta)]
        result, elapsed = run_command(cmd)
        rows.append(("kallisto", cmd, result, elapsed, kallisto_idx))

    if not args.skip_kallistors:
        cmd = [
            args.kallistors_bin,
            "index",
            "-i",
            str(kallistors_idx),
            "-k",
            str(args.k),
            "-t",
            str(args.threads),
            "--timings",
            str(fasta),
        ]
        result, elapsed = run_command(cmd)
        rows.append(("kallistors", cmd, result, elapsed, kallistors_idx))

    report = []
    report.append("# Index build benchmark\n\n")
    report.append(f"Generated: {datetime.now().isoformat()}\n\n")
    report.append("## Environment\n")
    report.append(f"- OS: {platform.platform()}\n")
    report.append(f"- Machine: {platform.machine()}\n")
    report.append(f"- CPU cores: {os.cpu_count()}\n")
    report.append(f"- FASTA: {fasta}\n")
    report.append(f"- k: {args.k}\n")
    report.append(f"- kallistors threads: {args.threads}\n")
    report.append(f"- work dir: {work_dir}\n\n")
    report.append("## Results\n")
    report.append("| tool | exit | wall_s | index_size |\n")
    report.append("| --- | ---: | ---: | --- |\n")
    for tool, _cmd, result, elapsed, index_path in rows:
        report.append(
            f"| {tool} | {result.returncode} | {elapsed:.3f} | {format_size(index_path)} |\n"
        )
    report.append("\n## Commands\n")
    report.append("```bash\n")
    for _tool, cmd, _result, _elapsed, _index_path in rows:
        report.append(" ".join(cmd) + "\n")
    report.append("```\n\n")

    report.append("## Metadata\n")
    for tool, _cmd, result, _elapsed, index_path in rows:
        report.append(f"### {tool}\n")
        if result.returncode != 0 or not index_path.exists():
            report.append("Index was not produced.\n\n")
            continue
        if command_available(args.kallistors_bin):
            info = read_index_info(args.kallistors_bin, index_path)
            report.append("```json\n")
            report.append(json.dumps(info, indent=2, sort_keys=True) + "\n")
            report.append("```\n")
        else:
            report.append("kallistors index-info skipped: binary unavailable.\n\n")
        if command_available(args.kallisto_bin):
            inspect_status, inspect_output = inspect_with_kallisto(args.kallisto_bin, index_path)
            report.append(f"kallisto inspect exit: {inspect_status}\n\n")
            report.append("```text\n")
            report.append(inspect_output + "\n")
            report.append("```\n")
        report.append("\n")

    report.append("## Logs\n")
    for tool, _cmd, result, _elapsed, _index_path in rows:
        report.append(f"### {tool}\n")
        report.append("```text\n")
        combined = (result.stdout or "") + (result.stderr or "")
        report.append(tail(combined) + "\n")
        report.append("```\n\n")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(report), encoding="utf-8")
    print(f"wrote {out}")
    if args.keep and not args.work_dir:
        print(f"temporary indexes kept in {work_dir}")


if __name__ == "__main__":
    main()
