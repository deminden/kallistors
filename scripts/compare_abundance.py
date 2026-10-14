#!/usr/bin/env python3
"""Compare kallisto vs kallistors abundance.tsv outputs."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Dict, Iterable, Tuple


def load_abundance(path: Path) -> Dict[str, Dict[str, float]]:
    data: Dict[str, Dict[str, float]] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return data
        idx = {name: i for i, name in enumerate(header)}
        name_key = "name" if "name" in idx else "target_id" if "target_id" in idx else "transcript_id"
        for row in reader:
            if not row:
                continue
            name = row[idx.get(name_key, 0)]
            est_counts = float(row[idx.get("est_counts", 3)])
            tpm = float(row[idx.get("tpm", 4)])
            data[name] = {
                "est_counts": est_counts,
                "tpm": tpm,
                "eff_length": float(row[idx["eff_length"]]) if "eff_length" in idx else 0.0,
            }
    return data


def rel_error(a: float, b: float) -> float:
    denom = max(1.0, a, b)
    return abs(a - b) / denom


def audit_summary(
    kallisto: Dict[str, Dict[str, float]],
    kallistors: Dict[str, Dict[str, float]],
    significant_digits: int,
) -> dict:
    """Separate output rounding from differences in expressed transcripts and effective lengths."""
    shared = sorted(kallisto.keys() & kallistors.keys())
    if not shared:
        raise ValueError("abundance files have no shared transcript IDs")

    def rounded(value):
        return float(format(value, f".{significant_digits}g"))

    def quantiles(values):
        ordered = sorted(values)
        if not ordered:
            return {"median": None, "p90": None, "p99": None, "max": None}
        return {
            "median": statistics.median(ordered),
            "p90": ordered[int(0.9 * (len(ordered) - 1))],
            "p99": ordered[int(0.99 * (len(ordered) - 1))],
            "max": ordered[-1],
        }

    result = {
        "shared_transcripts": len(shared),
        "missing_from_kallistors": sorted(kallisto.keys() - kallistors.keys()),
        "missing_from_kallisto": sorted(kallistors.keys() - kallisto.keys()),
        "comparison_significant_digits": significant_digits,
        "metrics": {},
    }
    for metric in ("est_counts", "tpm", "eff_length"):
        errors = [abs(kallisto[t][metric] - kallistors[t][metric]) for t in shared]
        worst = sorted(
            shared,
            key=lambda t: abs(kallisto[t][metric] - kallistors[t][metric]),
            reverse=True,
        )[:10]
        result["metrics"][metric] = {
            "absolute_errors": quantiles(errors),
            "sum_absolute_error": sum(errors),
            "unequal_after_rounding": sum(
                rounded(kallisto[t][metric]) != rounded(kallistors[t][metric])
                for t in shared
            ),
            "largest_differences": [
                {
                    "transcript": t,
                    "kallisto": kallisto[t][metric],
                    "kallistors": kallistors[t][metric],
                }
                for t in worst
            ],
        }
    positive = [t for t in shared if kallisto[t]["tpm"] > 0.0]
    result["positive_kallisto_tpm_transcripts"] = len(positive)
    result["positive_tpm_relative_errors"] = quantiles(
        [
            abs(kallisto[t]["tpm"] - kallistors[t]["tpm"]) / kallisto[t]["tpm"]
            for t in positive
        ]
    )
    reference_tpm_sum = sum(kallisto[t]["tpm"] for t in shared)
    result["tpm_mass_shift_percent"] = (
        50.0 * result["metrics"]["tpm"]["sum_absolute_error"] / reference_tpm_sum
        if reference_tpm_sum else None
    )
    return result


def build_rows(
    kallisto: Dict[str, Dict[str, float]],
    kallistors: Dict[str, Dict[str, float]],
    metric: str,
) -> Iterable[Tuple[float, str, float, float, float]]:
    for name, k_vals in kallisto.items():
        k_val = k_vals.get(metric, 0.0)
        o_val = kallistors.get(name, {}).get(metric, 0.0)
        abs_diff = abs(k_val - o_val)
        rel_diff = rel_error(k_val, o_val)
        yield (rel_diff, name, abs_diff, k_val, o_val)


def report_top(rows: Iterable[Tuple[float, str, float, float, float]], top: int) -> None:
    print("rank\ttranscript\trel_diff\tabs_diff\tkallisto\tkallistors")
    for idx, (rel_diff, name, abs_diff, k_val, o_val) in enumerate(rows, start=1):
        print(
            f"{idx}\t{name}\t{rel_diff:.6g}\t{abs_diff:.6g}\t{k_val:.6g}\t{o_val:.6g}"
        )
        if idx >= top:
            break


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kallisto", required=True, help="Path to kallisto abundance.tsv")
    parser.add_argument("--kallistors", required=True, help="Path to kallistors abundance.tsv")
    parser.add_argument("--top", type=int, default=50, help="Rows per ranking (default: 50)")
    parser.add_argument(
        "--summary-json",
        type=Path,
        help="Write a quantification audit summary instead of rankings",
    )
    parser.add_argument(
        "--significant-digits",
        type=int,
        default=6,
        help="Precision used to separate text rounding from other differences (default: kallisto's 6)",
    )
    parser.add_argument(
        "--sort-by",
        choices=("tpm", "est_counts", "max"),
        default="max",
        help="Ranking metric (default: max)",
    )
    args = parser.parse_args()
    if args.significant_digits < 1:
        parser.error("--significant-digits must be positive")

    k_path = Path(args.kallisto)
    o_path = Path(args.kallistors)
    if not k_path.exists() or not o_path.exists():
        missing = [str(p) for p in (k_path, o_path) if not p.exists()]
        print(f"Skipping: missing file(s): {', '.join(missing)}")
        return

    kallisto = load_abundance(k_path)
    kallistors = load_abundance(o_path)
    if not kallisto or not kallistors:
        print("Skipping: empty abundance file(s).")
        return

    if args.summary_json:
        summary = audit_summary(kallisto, kallistors, args.significant_digits)
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote audit summary: {args.summary_json}")
        return

    metrics = ["est_counts", "tpm"] if args.sort_by == "max" else [args.sort_by]
    for metric in metrics:
        print(f"# top {args.top} by {metric}")
        rows = list(build_rows(kallisto, kallistors, metric))
        rows.sort(key=lambda item: (-item[0], item[1]))
        report_top(rows, args.top)
        print()

    if args.sort_by == "max":
        combined = []
        for name, k_vals in kallisto.items():
            k_est = k_vals.get("est_counts", 0.0)
            k_tpm = k_vals.get("tpm", 0.0)
            o_vals = kallistors.get(name, {})
            o_est = o_vals.get("est_counts", 0.0)
            o_tpm = o_vals.get("tpm", 0.0)
            rel = max(rel_error(k_est, o_est), rel_error(k_tpm, o_tpm))
            abs_diff = max(abs(k_est - o_est), abs(k_tpm - o_tpm))
            combined.append((rel, name, abs_diff, k_est, o_est, k_tpm, o_tpm))
        combined.sort(key=lambda item: (-item[0], item[1]))
        print(f"# top {args.top} by max(est_counts, tpm)")
        print(
            "rank\ttranscript\trel_diff\tabs_diff\tkallisto_est\tkallistors_est\t"
            "kallisto_tpm\tkallistors_tpm"
        )
        for idx, (rel, name, abs_diff, k_est, o_est, k_tpm, o_tpm) in enumerate(
            combined, start=1
        ):
            print(
                f"{idx}\t{name}\t{rel:.6g}\t{abs_diff:.6g}\t"
                f"{k_est:.6g}\t{o_est:.6g}\t{k_tpm:.6g}\t{o_tpm:.6g}"
            )
            if idx >= args.top:
                break


if __name__ == "__main__":
    main()
