#!/usr/bin/env python3
"""Run small-sample real-data parity checks against instrumented kallisto."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
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


def run_command_env(cmd: List[str], env: Dict[str, str]) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    merged.update(env)
    return subprocess.run(cmd, text=True, capture_output=True, check=False, env=merged)


def parse_env_assignments(values: List[str]) -> Dict[str, str]:
    env: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise argparse.ArgumentTypeError(f"invalid --fast-env assignment: {value}")
        key, raw = value.split("=", 1)
        env[key.strip()] = raw.strip()
    return env


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
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore")


def collect_headers(path: Path) -> List[str]:
    headers: List[str] = []
    with gzip_open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not seq or not plus or not qual:
                break
            headers.append(header.strip().split()[0].lstrip("@"))
    return headers


def extract_pair_subset(src1: Path, src2: Path, header: str, out1: Path, out2: Path) -> None:
    with gzip_open_text(src1) as handle1, gzip_open_text(src2) as handle2:
        with gzip.open(out1, "wt", encoding="utf-8") as dst1, gzip.open(
            out2, "wt", encoding="utf-8"
        ) as dst2:
            while True:
                h1 = handle1.readline()
                h2 = handle2.readline()
                if not h1 or not h2:
                    break
                s1, p1, q1 = handle1.readline(), handle1.readline(), handle1.readline()
                s2, p2, q2 = handle2.readline(), handle2.readline(), handle2.readline()
                if not all([s1, p1, q1, s2, p2, q2]):
                    break
                key1 = h1.strip().split()[0].lstrip("@")
                key2 = h2.strip().split()[0].lstrip("@")
                if key1 == header and key2 == header:
                    dst1.write(h1)
                    dst1.write(s1)
                    dst1.write(p1)
                    dst1.write(q1)
                    dst2.write(h2)
                    dst2.write(s2)
                    dst2.write(p2)
                    dst2.write(q2)
                    return
    raise SystemExit(f"failed to extract paired read {header}")


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
    parser.add_argument(
        "--fast-env",
        action="append",
        default=[],
        help="environment assignment for fast branch, e.g. KALLISTORS_FAST_DELTA_BOUNDED_INCREMENTAL_SCAN=1",
    )
    parser.add_argument(
        "--compare-out",
        default="data/subsets/fast_trace_compare.tsv",
        help="trace-compare report for the first mismatched subset/read",
    )
    parser.add_argument(
        "--kallisto-trace-out",
        default="data/subsets/fast_kallisto_ec_trace.tsv",
        help="instrumented kallisto ec-trace for the first mismatched pair",
    )
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
    first_mismatch_inputs: Tuple[int, Path, Path] | None = None
    fast_env = parse_env_assignments(args.fast_env)

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
            o_result = run_command_env(ours_cmd, fast_env) if fast_env else run_command(ours_cmd)
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
            if not run_info_match and first_mismatch_inputs is None:
                first_mismatch_inputs = (reads_n, mate1, mate2)

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
        if first_mismatch_inputs and args.mode == "paired":
            reads_n, mate1, mate2 = first_mismatch_inputs
            compare_out = Path(args.compare_out)
            compare_out.parent.mkdir(parents=True, exist_ok=True)
            headers = collect_headers(mate1)
            with tempfile.TemporaryDirectory(prefix="kallistors-fast-compare-") as tmp:
                tmp_path = Path(tmp)
                read_list = tmp_path / "reads.txt"
                read_list.write_text("\n".join(headers) + "\n", encoding="utf-8")
                compare_cmd = [
                    str(kallistors_bin),
                    "trace-compare",
                    "--index",
                    str(index),
                    "--reads",
                    str(mate1),
                    "--reads2",
                    str(mate2),
                    "--read-list",
                    str(read_list),
                    "--out",
                    str(compare_out),
                ]
                compare_result = (
                    run_command_env(compare_cmd, fast_env) if fast_env else run_command(compare_cmd)
                )
                if compare_result.returncode != 0:
                    raise SystemExit(
                        f"trace-compare failed for n={reads_n}: {compare_result.stderr}"
                    )
                first_diff = None
                with compare_out.open("r", encoding="utf-8") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        if (
                            row["aligned_baseline"] != row["aligned_fast"]
                            or row["reason_baseline"] != row["reason_fast"]
                            or row["ec_baseline"] != row["ec_fast"]
                        ):
                            first_diff = row
                            break
                if first_diff:
                    print(
                        "first divergent paired read:",
                        first_diff["read"],
                        first_diff["reason_baseline"],
                        first_diff["reason_fast"],
                    )
                    with tempfile.TemporaryDirectory(prefix="kallisto-fast-first-pair-") as tmp2:
                        tmp_path = Path(tmp2)
                        one1 = tmp_path / "pair_1.fastq.gz"
                        one2 = tmp_path / "pair_2.fastq.gz"
                        extract_pair_subset(mate1, mate2, first_diff["read"], one1, one2)
                        k_out = tmp_path / "kallisto"
                        k_out.mkdir()
                        kallisto_trace = Path(args.kallisto_trace_out)
                        kallisto_trace.parent.mkdir(parents=True, exist_ok=True)
                        kallisto_cmd = [
                            str(kallisto_bin),
                            "quant",
                            "-i",
                            str(index),
                            "-o",
                            str(k_out),
                            "-t",
                            "1",
                            "--ec-trace",
                            str(kallisto_trace),
                            "--ec-trace-max-reads",
                            "0",
                            str(one1),
                            str(one2),
                        ]
                        k_result = run_command(kallisto_cmd)
                        if k_result.returncode != 0:
                            raise SystemExit(
                                f"kallisto ec-trace failed for {first_diff['read']}: {k_result.stderr}"
                            )
        raise SystemExit("parity mismatch detected; inspect the report for the first failing subset")


if __name__ == "__main__":
    main()
