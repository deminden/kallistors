#!/usr/bin/env python3
"""Compare kallisto accepted k-mer positions with kallistors accepted hits."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict


def norm_read_id(text: str) -> str:
    text = text.strip()
    if text.startswith("@"):
        text = text[1:]
    return text.split()[0]


def load_read_list(path: str) -> list[str]:
    out: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            rid = norm_read_id(line)
            if rid:
                out.append(rid)
    return sorted(set(out))


def parse_int_set(text: str) -> set[int]:
    text = text.strip()
    if not text:
        return set()
    return {int(v) for v in text.split(",") if v}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kallisto-kmer-dump", required=True)
    parser.add_argument("--kallistors-hits", required=True)
    parser.add_argument("--kallistors-positions-visited", required=True)
    parser.add_argument("--read-list", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--detail-out")
    args = parser.parse_args()

    reads = load_read_list(args.read_list)
    read_set = set(reads)

    kallisto_rows: dict[str, dict[int, list[dict[str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with open(args.kallisto_kmer_dump, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read_name"])
            if rid not in read_set:
                continue
            pos = int(row["read_pos"])
            kallisto_rows[rid][pos].append(row)

    ours_hits: dict[str, dict[int, set[str]]] = defaultdict(lambda: defaultdict(set))
    with open(args.kallistors_hits, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read"])
            if rid not in read_set:
                continue
            pos = int(row["read_pos"])
            ours_hits[rid][pos].add(row["ec"])

    visited_by_read: dict[str, set[int]] = {}
    with open(args.kallistors_positions_visited, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read"])
            if rid not in read_set:
                continue
            visited_by_read[rid] = parse_int_set(row["positions_visited"])

    with open(args.summary_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "read",
                "k_accept_count",
                "ours_hit_positions_all",
                "ours_hit_positions_visited",
                "ours_synthetic_hit_positions",
                "extra_positions_visited",
                "missing_positions_visited",
                "first_extra_pos",
                "first_extra_ec_ids",
                "first_missing_pos",
                "first_missing_k_is_hit",
                "first_missing_k_special",
                "first_missing_k_overcrowded",
                "first_missing_k_jump_used",
                "k_accept_special_count",
                "k_accept_overcrowded_count",
                "k_accept_jump_count",
            ]
        )

        detail_writer = None
        detail_handle = None
        if args.detail_out:
            detail_handle = open(args.detail_out, "w", encoding="utf-8", newline="")
            detail_writer = csv.writer(detail_handle, delimiter="\t")
            detail_writer.writerow(
                [
                    "read",
                    "pos",
                    "k_accept",
                    "ours_hit",
                    "ours_hit_visited",
                    "ours_ec_ids",
                    "k_is_hit",
                    "k_hit_action",
                    "k_skip_reason",
                    "k_special",
                    "k_overcrowded",
                    "k_jump_used",
                ]
            )

        for rid in reads:
            by_pos = kallisto_rows.get(rid, {})
            visited = visited_by_read.get(rid, set())
            ours_by_pos = ours_hits.get(rid, {})

            k_accept_positions = sorted(
                pos
                for pos, rows in by_pos.items()
                if any(r.get("hit_action") == "accept" for r in rows)
            )
            k_accept_set = set(k_accept_positions)

            ours_all_positions = sorted(ours_by_pos)
            ours_visited_positions = sorted(pos for pos in ours_all_positions if pos in visited)
            ours_visited_set = set(ours_visited_positions)

            synthetic_positions = sorted(pos for pos in ours_all_positions if pos not in visited)
            extra_positions = sorted(pos for pos in ours_visited_positions if pos not in k_accept_set)
            missing_positions = sorted(pos for pos in k_accept_positions if pos not in ours_visited_set)

            first_extra = extra_positions[0] if extra_positions else None
            first_extra_ecs = "|".join(sorted(ours_by_pos.get(first_extra, set()))) if first_extra is not None else "-"

            first_missing = missing_positions[0] if missing_positions else None
            miss_row = None
            if first_missing is not None:
                rows = by_pos.get(first_missing, [])
                accept_rows = [r for r in rows if r.get("hit_action") == "accept"]
                miss_row = accept_rows[0] if accept_rows else (rows[0] if rows else None)

            k_accept_special = 0
            k_accept_overcrowded = 0
            k_accept_jump = 0
            for pos in k_accept_positions:
                rows = by_pos.get(pos, [])
                if any(r.get("minimizer_is_special") == "1" for r in rows):
                    k_accept_special += 1
                if any(r.get("minimizer_is_overcrowded") == "1" for r in rows):
                    k_accept_overcrowded += 1
                if any(r.get("jump_used") == "1" for r in rows):
                    k_accept_jump += 1

            writer.writerow(
                [
                    rid,
                    len(k_accept_positions),
                    len(ours_all_positions),
                    len(ours_visited_positions),
                    len(synthetic_positions),
                    ",".join(str(v) for v in extra_positions) if extra_positions else "-",
                    ",".join(str(v) for v in missing_positions) if missing_positions else "-",
                    first_extra if first_extra is not None else "-",
                    first_extra_ecs,
                    first_missing if first_missing is not None else "-",
                    miss_row.get("is_hit", "-") if miss_row else "-",
                    miss_row.get("minimizer_is_special", "-") if miss_row else "-",
                    miss_row.get("minimizer_is_overcrowded", "-") if miss_row else "-",
                    miss_row.get("jump_used", "-") if miss_row else "-",
                    k_accept_special,
                    k_accept_overcrowded,
                    k_accept_jump,
                ]
            )

            if detail_writer is not None:
                all_positions = sorted(set(by_pos) | set(ours_by_pos))
                for pos in all_positions:
                    rows = by_pos.get(pos, [])
                    accept = any(r.get("hit_action") == "accept" for r in rows)
                    row = None
                    if rows:
                        accept_rows = [r for r in rows if r.get("hit_action") == "accept"]
                        row = accept_rows[0] if accept_rows else rows[0]
                    ours_hit = pos in ours_by_pos
                    ours_hit_visited = ours_hit and pos in visited
                    detail_writer.writerow(
                        [
                            rid,
                            pos,
                            int(accept),
                            int(ours_hit),
                            int(ours_hit_visited),
                            "|".join(sorted(ours_by_pos.get(pos, set()))) if ours_hit else "-",
                            row.get("is_hit", "-") if row else "-",
                            row.get("hit_action", "-") if row else "-",
                            row.get("skip_reason", "-") if row else "-",
                            row.get("minimizer_is_special", "-") if row else "-",
                            row.get("minimizer_is_overcrowded", "-") if row else "-",
                            row.get("jump_used", "-") if row else "-",
                        ]
                    )

        if detail_handle is not None:
            detail_handle.close()


if __name__ == "__main__":
    main()
