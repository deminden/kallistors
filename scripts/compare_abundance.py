#!/usr/bin/env python3
"""Compare kallisto vs kallistors abundance.tsv outputs."""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple


def load_abundance(path: Path) -> Dict[str, Dict[str, float]]:
    data = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_id", "length", "eff_length", "est_counts", "tpm"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(f"{path}: missing abundance columns")
        for line_number, row in enumerate(reader, 2):
            name = row["target_id"]
            if not name or name in data or None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{line_number}: malformed row or duplicate target")
            length = int(row["length"])
            values = {key: float(row[key]) for key in ("eff_length", "est_counts", "tpm")}
            if length < 0 or any(not math.isfinite(v) or v < 0 for v in values.values()):
                raise ValueError(f"{path}:{line_number}: invalid abundance value")
            data[name] = values
    if not data:
        raise ValueError(f"{path}: empty abundance table")
    return data


def rel_error(a: float, b: float) -> float:
    denom = max(1.0, a, b)
    return abs(a - b) / denom


def quantiles(values):
    """Linear interpolation at (N-1)*p, matching NumPy's default percentile rule."""
    ordered = sorted(values)
    probabilities = {"median": .5, "p90": .9, "p99": .99,
                     "p99_9": .999, "p99_99": .9999, "max": 1.0}
    if not ordered:
        return dict.fromkeys(probabilities)
    result = {}
    for name, probability in probabilities.items():
        position = probability * (len(ordered) - 1)
        lower, upper = math.floor(position), math.ceil(position)
        weight = position - lower
        result[name] = ordered[lower] * (1 - weight) + ordered[upper] * weight
    return result


def load_hdf5_abundance(path):
    """Use saved f64 counts/effective lengths; derive TPM before TSV rounding."""
    import h5py
    import numpy as np
    with h5py.File(path, "r") as handle:
        names = handle["aux/ids"].asstr()[()]
        counts = handle["est_counts"][()]
        lengths = handle["aux/eff_lengths"][()]
    if (counts.shape != names.shape or lengths.shape != names.shape or not len(names)
            or len(set(names)) != len(names) or not all(names)
            or not np.isfinite(counts).all() or not np.isfinite(lengths).all()
            or (counts < 0).any() or (lengths < 0).any()
            or ((lengths == 0) & (counts > 0)).any()):
        raise ValueError(f"{path}: invalid abundance arrays or transcript IDs")
    rates = np.divide(counts, lengths, out=np.zeros_like(counts), where=lengths > 0)
    total = rates.sum()
    tpm = rates * (1_000_000 / total) if total else rates
    return {name: {"est_counts": float(count), "eff_length": float(length), "tpm": float(value)}
            for name, count, length, value in zip(names, counts, lengths, tpm)}


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

    result = {
        "shared_transcripts": len(shared),
        "missing_from_kallistors": sorted(kallisto.keys() - kallistors.keys()),
        "missing_from_kallisto": sorted(kallistors.keys() - kallisto.keys()),
        "comparison_significant_digits": significant_digits,
        "percentile_method": "linear interpolation at (N-1)*p; zeros included",
        "metrics": {},
    }
    expressed = [t for t in shared if kallisto[t]["est_counts"] > 0 or kallistors[t]["est_counts"] > 0]
    for metric in ("est_counts", "tpm", "eff_length"):
        errors = [abs(kallisto[t][metric] - kallistors[t][metric]) for t in shared]
        worst = heapq.nlargest(
            10,
            shared,
            key=lambda t: abs(kallisto[t][metric] - kallistors[t][metric]),
        )
        result["metrics"][metric] = {
            "absolute_errors": quantiles(errors),
            "expressed_in_either_transcripts": len(expressed),
            "expressed_in_either_absolute_errors": quantiles(
                abs(kallisto[t][metric] - kallistors[t][metric]) for t in expressed),
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
    parser.add_argument("--hdf5", action="store_true", help="Compare native HDF5 f64 values instead of TSV")
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
    try:
        loader = load_hdf5_abundance if args.hdf5 else load_abundance
        kallisto = loader(k_path)
        kallistors = loader(o_path)
        if kallisto.keys() != kallistors.keys():
            raise ValueError("Transcript sets differ; the outputs do not share the same reference")
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Abundance comparison failed: {error}\n")

    if args.summary_json:
        summary = audit_summary(kallisto, kallistors, args.significant_digits)
        summary["source_precision"] = "native HDF5 f64; TPM derived from counts/effective lengths" if args.hdf5 else "TSV text"
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
