#!/usr/bin/env python3
"""Compare complete EC count maps independently of per-run EC numbering."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path


def load_counts(path: Path) -> Counter[tuple[int, ...]]:
    counts: Counter[tuple[int, ...]] = Counter()
    with path.open() as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected EC ID, count, transcripts")
            count = int(fields[1])
            transcripts = tuple(sorted(int(t) for t in fields[2].split(",") if t))
            if count < 0 or (count and not transcripts):
                raise ValueError(f"{path}:{line_number}: invalid EC count or empty counted class")
            if count:
                counts[transcripts] += count
    return counts


def compare(reference: Counter, candidate: Counter) -> dict:
    different = sorted(ec for ec in reference.keys() | candidate.keys() if reference[ec] != candidate[ec])
    return {
        "exact": not different,
        "reference_nonzero_ecs": len(reference),
        "candidate_nonzero_ecs": len(candidate),
        "reference_total_counts": sum(reference.values()),
        "candidate_total_counts": sum(candidate.values()),
        "different_classes": len(different),
        "count_l1_difference": sum(abs(reference[ec] - candidate[ec]) for ec in different),
        "first_differences": [
            {"transcripts": ec, "reference": reference[ec], "candidate": candidate[ec]}
            for ec in different[:20]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    summary = compare(load_counts(args.reference), load_counts(args.candidate))
    result = json.dumps(summary, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(result)
    print(result, end="")
    raise SystemExit(0 if summary["exact"] else 1)


if __name__ == "__main__":
    main()
