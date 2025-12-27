#!/usr/bin/env python3
"""Benchmark kallisto vs kallistors on a real FASTQ subset."""

import argparse
import json
import os
import platform
import resource
import shutil
import subprocess
import tempfile
import time
from datetime import datetime
from typing import Dict, Iterable, List, Tuple


def read_abundance(path: str) -> Dict[str, Tuple[int, float, float, float]]:
    out: Dict[str, Tuple[int, float, float, float]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            if i == 0:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            target = parts[0]
            length = int(parts[1])
            eff_len = float(parts[2])
            est = float(parts[3])
            tpm = float(parts[4])
            out[target] = (length, eff_len, est, tpm)
    return out


def percentile(sorted_vals: List[float], pct: float) -> float:
    if not sorted_vals:
        return 0.0
    if pct <= 0:
        return sorted_vals[0]
    if pct >= 100:
        return sorted_vals[-1]
    idx = int(round((pct / 100.0) * (len(sorted_vals) - 1)))
    return sorted_vals[idx]


def rel_error(a: float, b: float) -> float:
    denom = max(1.0, abs(a), abs(b))
    return abs(a - b) / denom


def child_cpu_time() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def run_command(cmd: List[str]) -> Tuple[subprocess.CompletedProcess, float, float]:
    cpu_before = child_cpu_time()
    start = time.perf_counter()
    result = subprocess.run(cmd, check=False, text=True, capture_output=True)
    elapsed = time.perf_counter() - start
    cpu_after = child_cpu_time()
    return result, elapsed, max(0.0, cpu_after - cpu_before)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--index",
        default="data/Human_kallisto_index",
        help="Path to kallisto index",
    )
    parser.add_argument(
        "--reads",
        default=(
            "data/"
            "SRR13638690_RNA-seq_of_homo_sapiens_temporal_muscle_of_low_grade_migraine_1_trimmed_subset.fastq.gz"
        ),
        help="Subset FASTQ(.gz)",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--fragment-length", type=int, default=200)
    parser.add_argument("--fragment-length-sd", type=int, default=20)
    parser.add_argument("--kallisto-bin", default="kallisto")
    parser.add_argument("--kallistors-bin", default="kallistors-cli")
    parser.add_argument("--out", default="bench_latest.md")
    args = parser.parse_args()

    if not os.path.exists(args.index):
        raise SystemExit(f"missing index: {args.index}")
    if not os.path.exists(args.reads):
        raise SystemExit(f"missing reads: {args.reads}")

    if shutil.which(args.kallisto_bin) is None:
        raise SystemExit(f"kallisto not found: {args.kallisto_bin}")
    if shutil.which(args.kallistors_bin) is None:
        raise SystemExit(f"kallistors-cli not found: {args.kallistors_bin}")

    with tempfile.TemporaryDirectory(prefix="kallistors-bench-") as tmp:
        k_out = os.path.join(tmp, "kallisto")
        ours_out = os.path.join(tmp, "kallistors")

        kallisto_cmd = [
            args.kallisto_bin,
            "quant",
            "--single",
            "-l",
            str(args.fragment_length),
            "-s",
            str(args.fragment_length_sd),
            "-i",
            args.index,
            "-o",
            k_out,
            "-t",
            str(args.threads),
            args.reads,
        ]

        kallistors_cmd = [
            args.kallistors_bin,
            "quant",
            "--single",
            "-l",
            str(args.fragment_length),
            "-s",
            str(args.fragment_length_sd),
            "-i",
            args.index,
            "-o",
            ours_out,
            "-t",
            str(args.threads),
            args.reads,
        ]

        k_result, k_wall, k_cpu = run_command(kallisto_cmd)
        if k_result.returncode != 0:
            raise SystemExit(f"kallisto failed: {k_result.stderr.strip()}")

        o_result, o_wall, o_cpu = run_command(kallistors_cmd)
        if o_result.returncode != 0:
            raise SystemExit(f"kallistors-cli failed: {o_result.stderr.strip()}")

        with open(os.path.join(k_out, "run_info.json"), "r", encoding="utf-8") as handle:
            k_info = json.load(handle)
        with open(os.path.join(ours_out, "run_info.json"), "r", encoding="utf-8") as handle:
            o_info = json.load(handle)

        k_abund = read_abundance(os.path.join(k_out, "abundance.tsv"))
        o_abund = read_abundance(os.path.join(ours_out, "abundance.tsv"))

        errors_est: List[float] = []
        errors_tpm: List[float] = []
        max_eff_len_diff = 0.0
        missing = 0
        for target, (length_k, eff_k, est_k, tpm_k) in k_abund.items():
            ours = o_abund.get(target)
            if ours is None:
                missing += 1
                continue
            length_o, eff_o, est_o, tpm_o = ours
            if length_k != length_o:
                raise SystemExit(f"length mismatch for {target}: {length_k} vs {length_o}")
            max_eff_len_diff = max(max_eff_len_diff, abs(eff_k - eff_o))
            errors_est.append(rel_error(est_k, est_o))
            errors_tpm.append(rel_error(tpm_k, tpm_o))

        errors_est.sort()
        errors_tpm.sort()

        report = []
        report.append("# Real-data benchmark vs kallisto\n")
        report.append(f"Generated: {datetime.now().isoformat()}\n")
        report.append("\n## Environment\n")
        report.append(f"- OS: {platform.platform()}\n")
        report.append(f"- Machine: {platform.machine()}\n")
        report.append(f"- Processor: {platform.processor()}\n")
        report.append(f"- CPU cores: {os.cpu_count()}\n")
        report.append("\n## Inputs\n")
        report.append(f"- Index: {args.index}\n")
        report.append(f"- Reads: {args.reads}\n")
        report.append(
            f"- Fragment length: {args.fragment_length} +- {args.fragment_length_sd}\n"
        )
        report.append(f"- Threads: {args.threads}\n")
        report.append("\n## Commands\n")
        report.append("```bash\n")
        report.append(" ".join(kallisto_cmd) + "\n")
        report.append(" ".join(kallistors_cmd) + "\n")
        report.append("```\n")
        report.append("\n## Timings\n")
        report.append(f"- kallisto wall: {k_wall:.3f}s, cpu: {k_cpu:.3f}s\n")
        report.append(f"- kallistors wall: {o_wall:.3f}s, cpu: {o_cpu:.3f}s\n")
        report.append("\n## run_info parity\n")
        report.append(
            "| field | kallisto | kallistors |\n| --- | --- | --- |\n"
        )
        fields = [
            "n_targets",
            "n_processed",
            "n_pseudoaligned",
            "n_unique",
            "p_pseudoaligned",
            "p_unique",
        ]
        for field in fields:
            report.append(f"| {field} | {k_info.get(field)} | {o_info.get(field)} |\n")
        report.append("\n## abundance.tsv parity\n")
        report.append(f"- Missing transcripts: {missing}\n")
        report.append(f"- Max eff_length abs diff: {max_eff_len_diff:.6g}\n")
        report.append("\nPercentiles for relative error (abs(a-b)/max(1,a,b)):\n")
        report.append("\n| metric | p50 | p90 | p99 | max |\n| --- | --- | --- | --- | --- |\n")
        report.append(
            "| est_counts |"
            f" {percentile(errors_est, 50):.3g} |"
            f" {percentile(errors_est, 90):.3g} |"
            f" {percentile(errors_est, 99):.3g} |"
            f" {errors_est[-1] if errors_est else 0.0:.3g} |\n"
        )
        report.append(
            "| tpm |"
            f" {percentile(errors_tpm, 50):.3g} |"
            f" {percentile(errors_tpm, 90):.3g} |"
            f" {percentile(errors_tpm, 99):.3g} |"
            f" {errors_tpm[-1] if errors_tpm else 0.0:.3g} |\n"
        )

        with open(args.out, "w", encoding="utf-8") as handle:
            handle.writelines(report)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
