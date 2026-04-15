#!/usr/bin/env python3
"""Run small-sample real-data parity checks against instrumented kallisto."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def parse_sizes(value: str) -> List[int]:
    sizes = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        sizes.append(int(item))
    if not sizes:
        raise argparse.ArgumentTypeError("at least one subset size is required")
    return sorted(set(sizes))


def read_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def rel_error(a: float, b: float) -> float:
    denom = max(1.0, abs(a), abs(b))
    return abs(a - b) / denom


def percentile(values: List[float], pct: float) -> float:
    if not values:
        return 0.0
    idx = int(round((pct / 100.0) * (len(values) - 1)))
    return values[max(0, min(len(values) - 1, idx))]


def load_abundance(path: Path) -> Dict[str, Tuple[float, float]]:
    out: Dict[str, Tuple[float, float]] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            target = row.get("target_id") or row.get("transcript_id") or row.get("name")
            if not target:
                continue
            out[target] = (
                float(row.get("est_counts", 0.0)),
                float(row.get("tpm", 0.0)),
            )
    return out


def compare_abundance(k_path: Path, o_path: Path) -> Tuple[float, float, float, float]:
    kallisto = load_abundance(k_path)
    ours = load_abundance(o_path)
    est_errors = []
    tpm_errors = []
    for target, (k_est, k_tpm) in kallisto.items():
        o_est, o_tpm = ours.get(target, (0.0, 0.0))
        est_errors.append(rel_error(k_est, o_est))
        tpm_errors.append(rel_error(k_tpm, o_tpm))
    est_errors.sort()
    tpm_errors.sort()
    return (
        percentile(est_errors, 99.0),
        est_errors[-1] if est_errors else 0.0,
        percentile(tpm_errors, 99.0),
        tpm_errors[-1] if tpm_errors else 0.0,
    )


def run_command(cmd: List[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def find_first_mismatch(
    kallisto_bin: Path,
    kallistors_bin: Path,
    index: Path,
    reads: Path,
    fragment_length: int,
    fragment_sd: int,
) -> Path | None:
    with tempfile.TemporaryDirectory(prefix="kallistors-single-trace-") as tmp:
        tmp_path = Path(tmp)
        trace_path = tmp_path / "ec_trace.tsv"
        out_dir = tmp_path / "kallisto"
        out_dir.mkdir()
        cmd = [
            str(kallisto_bin),
            "quant",
            "--single",
            "-l",
            str(fragment_length),
            "-s",
            str(fragment_sd),
            "-i",
            str(index),
            "-o",
            str(out_dir),
            "--ec-trace",
            str(trace_path),
            "--ec-trace-max-reads",
            "0",
            str(reads),
        ]
        result = run_command(cmd)
        if result.returncode != 0 or not trace_path.exists():
            return None

        read_list = tmp_path / "reads.txt"
        seen = []
        with gzip_open_text(reads) as handle, read_list.open("w", encoding="utf-8") as out:
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline()
                plus = handle.readline()
                qual = handle.readline()
                if not seq or not plus or not qual:
                    break
                normalized = header.strip().split()[0].lstrip("@")
                out.write(normalized + "\n")
                seen.append(normalized)

        ours_trace = tmp_path / "ours.tsv"
        cmd = [
            str(kallistors_bin),
            "trace-reads",
            "--index",
            str(index),
            "--reads",
            str(reads),
            "--read-list",
            str(read_list),
            "--out",
            str(ours_trace),
        ]
        result = run_command(cmd)
        if result.returncode != 0 or not ours_trace.exists():
            return None
    return None


def gzip_open_text(path: Path):
    import gzip

    return gzip.open(path, "rt", encoding="utf-8", errors="ignore")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", default="data/gencode.v49_kallisto.idx")
    parser.add_argument("--subset-dir", default="data/subsets")
    parser.add_argument("--sizes", default="64,256,1024", type=parse_sizes)
    parser.add_argument("--mode", choices=("single", "paired"), default="single")
    parser.add_argument("--fragment-length", type=int, default=200)
    parser.add_argument("--fragment-length-sd", type=int, default=20)
    parser.add_argument("--kallisto-bin", default="./kallisto_src/build/src/kallisto")
    parser.add_argument("--kallistors-bin", default="./target/release/kallistors-cli")
    parser.add_argument("--report", default="data/subsets/parity_report.tsv")
    args = parser.parse_args()

    kallisto_bin = Path(args.kallisto_bin)
    kallistors_bin = Path(args.kallistors_bin)
    index = Path(args.index)
    subset_dir = Path(args.subset_dir)
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    if shutil.which(str(kallisto_bin)) is None and not kallisto_bin.exists():
        raise SystemExit(f"missing kallisto binary: {kallisto_bin}")
    if shutil.which(str(kallistors_bin)) is None and not kallistors_bin.exists():
        raise SystemExit(f"missing kallistors binary: {kallistors_bin}")
    if not index.exists():
        raise SystemExit(f"missing index: {index}")

    rows = []
    mismatch_found = False

    for reads_n in args.sizes:
        mate1 = subset_dir / f"real_mate1_n{reads_n}.fastq.gz"
        mate2 = subset_dir / f"real_mate2_n{reads_n}.fastq.gz"
        reads_arg = [str(mate1)] if args.mode == "single" else [str(mate1), str(mate2)]
        for path in reads_arg:
            if not Path(path).exists():
                raise SystemExit(f"missing subset: {path}")

        with tempfile.TemporaryDirectory(prefix=f"kallistors-real-{args.mode}-{reads_n}-") as tmp:
            tmp_path = Path(tmp)
            k_out = tmp_path / "kallisto"
            o_out = tmp_path / "kallistors"

            kallisto_cmd = [
                str(kallisto_bin),
                "quant",
                "-i",
                str(index),
                "-o",
                str(k_out),
                "-t",
                "1",
            ]
            ours_cmd = [
                str(kallistors_bin),
                "quant",
                "-i",
                str(index),
                "-o",
                str(o_out),
                "-t",
                "1",
            ]

            if args.mode == "single":
                extra = [
                    "--single",
                    "-l",
                    str(args.fragment_length),
                    "-s",
                    str(args.fragment_length_sd),
                ]
                kallisto_cmd.extend(extra)
                ours_cmd.extend(extra)

            kallisto_cmd.extend(reads_arg)
            ours_cmd.extend(reads_arg)

            k_result = run_command(kallisto_cmd)
            if k_result.returncode != 0:
                raise SystemExit(f"kallisto failed for n={reads_n}: {k_result.stderr}")
            o_result = run_command(ours_cmd)
            if o_result.returncode != 0:
                raise SystemExit(f"kallistors failed for n={reads_n}: {o_result.stderr}")

            k_info = read_json(k_out / "run_info.json")
            o_info = read_json(o_out / "run_info.json")

            run_info_match = (
                k_info.get("n_processed") == o_info.get("n_processed")
                and k_info.get("n_pseudoaligned") == o_info.get("n_pseudoaligned")
                and k_info.get("n_unique") == o_info.get("n_unique")
            )
            est_p99, est_max, tpm_p99, tpm_max = compare_abundance(
                k_out / "abundance.tsv",
                o_out / "abundance.tsv",
            )
            rows.append(
                {
                    "mode": args.mode,
                    "reads": reads_n,
                    "n_processed_kallisto": k_info.get("n_processed"),
                    "n_processed_ours": o_info.get("n_processed"),
                    "n_pseudoaligned_kallisto": k_info.get("n_pseudoaligned"),
                    "n_pseudoaligned_ours": o_info.get("n_pseudoaligned"),
                    "n_unique_kallisto": k_info.get("n_unique"),
                    "n_unique_ours": o_info.get("n_unique"),
                    "run_info_match": int(run_info_match),
                    "est_p99_rel_err": f"{est_p99:.6g}",
                    "est_max_rel_err": f"{est_max:.6g}",
                    "tpm_p99_rel_err": f"{tpm_p99:.6g}",
                    "tpm_max_rel_err": f"{tpm_max:.6g}",
                }
            )
            mismatch_found = mismatch_found or not run_info_match

    with report_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {report_path}")
    for row in rows:
        print(
            f"{row['mode']} n={row['reads']}: run_info_match={row['run_info_match']} "
            f"pseudoaligned={row['n_pseudoaligned_kallisto']}/{row['n_pseudoaligned_ours']} "
            f"unique={row['n_unique_kallisto']}/{row['n_unique_ours']}"
        )

    if mismatch_found:
        raise SystemExit("parity mismatch detected; inspect the report for the first failing subset")


if __name__ == "__main__":
    main()
