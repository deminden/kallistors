#!/usr/bin/env python3
"""Run kallisto vs kallistors on a small synthetic dataset and report diffs."""

import argparse
import json
import os
import re
import subprocess
from pathlib import Path


def run(cmd, cwd=None, capture=False):
    if capture:
        return subprocess.run(cmd, cwd=cwd, check=True, text=True, capture_output=True)
    subprocess.run(cmd, cwd=cwd, check=True)
    return None


def parse_kallisto_pseudoaligned(run_info_path: Path) -> int:
    with run_info_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return int(data.get("n_pseudoaligned", 0))


def parse_kallistors_aligned(stdout: str) -> int:
    match = re.search(r"aligned:\s*(\d+)", stdout)
    if not match:
        return 0
    return int(match.group(1))


def load_abundance(path: Path):
    out = {}
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            parts = line.strip().split("\t")
            if "name" in idx:
                name = parts[idx["name"]]
            elif "target_id" in idx:
                name = parts[idx["target_id"]]
            else:
                name = parts[idx.get("transcript_id", 0)]
            est_counts = float(parts[idx.get("est_counts", 3)])
            tpm = float(parts[idx.get("tpm", 4)])
            out[name] = (est_counts, tpm)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", help="Dataset directory (default: data)")
    parser.add_argument("--single-overhang", action="store_true")
    parser.add_argument("--bias", action="store_true")
    parser.add_argument("--mean", type=float, default=100.0)
    parser.add_argument("--sd", type=float, default=20.0)
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    fasta = data_dir / "simple_transcripts.fa"
    reads = data_dir / "simple_reads.fq"
    index = data_dir / "simple.idx"
    k_out = data_dir / "kallisto_simple"
    ours_ec = data_dir / "simple.ec.txt"
    ours_bias = data_dir / "simple.bias.txt"
    ours_out = data_dir / "simple.abund.tsv"

    if not fasta.exists() or not reads.exists():
        raise SystemExit("Missing dataset; run scripts/generate_simple_dataset.py first.")

    if k_out.exists():
        for p in k_out.glob("*"):
            p.unlink()
    k_out.mkdir(parents=True, exist_ok=True)

    run(["kallisto", "index", "-i", str(index), str(fasta)])

    k_cmd = [
        "kallisto",
        "quant",
        "--single",
        "-l",
        str(args.mean),
        "-s",
        str(args.sd),
        "-i",
        str(index),
        "-o",
        str(k_out),
        str(reads),
    ]
    if args.single_overhang:
        k_cmd.append("--single-overhang")
    if args.bias:
        k_cmd.append("--bias")
    try:
        run(k_cmd)
    except subprocess.CalledProcessError as exc:
        print(f"kallisto failed with return code {exc.returncode}")
        return

    kallisto_aligned = parse_kallisto_pseudoaligned(k_out / "run_info.json")

    p_cmd = [
        "cargo",
        "run",
        "-p",
        "kallistors-cli",
        "--",
        "pseudoalign",
        "--index",
        str(index),
        "--reads",
        str(reads),
        "--fragment-length",
        str(args.mean),
        "--out",
        str(ours_ec),
    ]
    if args.single_overhang:
        p_cmd.append("--single-overhang")
    if args.bias:
        p_cmd.append("--bias")
        p_cmd.extend(["--bias-out", str(ours_bias)])

    result = run(p_cmd, capture=True)
    ours_aligned = parse_kallistors_aligned(result.stdout)

    q_cmd = [
        "cargo",
        "run",
        "-p",
        "kallistors-cli",
        "--",
        "quant",
        "--index",
        str(index),
        "--ec",
        str(ours_ec),
        "--out",
        str(ours_out),
        "--fragment-length",
        str(args.mean),
        "--fragment-length-sd",
        str(args.sd),
    ]
    if args.bias:
        q_cmd.append("--bias")
        q_cmd.extend(["--bias-counts", str(ours_bias)])
        q_cmd.extend(["--transcripts", str(fasta)])

    run(q_cmd)

    ours_abund = load_abundance(ours_out)
    kallisto_abund = load_abundance(k_out / "abundance.tsv")

    max_est = 0.0
    max_tpm = 0.0
    worst = None
    for name, (k_est, k_tpm) in kallisto_abund.items():
        o_est, o_tpm = ours_abund.get(name, (0.0, 0.0))
        diff_est = abs(k_est - o_est)
        diff_tpm = abs(k_tpm - o_tpm)
        if diff_est > max_est or diff_tpm > max_tpm:
            max_est = max(max_est, diff_est)
            max_tpm = max(max_tpm, diff_tpm)
            worst = (name, k_est, o_est, k_tpm, o_tpm)

    print(f"kallisto aligned: {kallisto_aligned}")
    print(f"kallistors aligned: {ours_aligned}")
    if worst:
        name, k_est, o_est, k_tpm, o_tpm = worst
        print(f"max est_counts diff: {max_est:.6f}")
        print(f"max tpm diff: {max_tpm:.6f}")
        print(
            f"worst transcript {name}: est_counts {k_est:.6f} vs {o_est:.6f}, "
            f"tpm {k_tpm:.6f} vs {o_tpm:.6f}"
        )
    else:
        print("no transcripts found in abundance files")


if __name__ == "__main__":
    main()
