#!/usr/bin/env python3
"""Summarize first divergence positions for kallisto vs kallistors traces."""

import argparse
import csv
from collections import defaultdict


def norm_read_id(text: str) -> str:
    text = text.strip()
    if text.startswith("@"):
        text = text[1:]
    return text.split()[0]


def read_id_set(path: str) -> set[str]:
    out: set[str] = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value:
                out.add(value)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kmer-dump", required=True, help="kallisto kmer_dump.tsv")
    parser.add_argument("--trace", required=True, help="kallistors trace TSV (read,aligned,reason,...)")
    parser.add_argument("--candidates", required=True, help="kallistors minimizer-candidates TSV")
    parser.add_argument("--visited", required=True, help="kallistors positions-visited TSV")
    parser.add_argument("--no-hits-ok", required=True, help="read IDs where kallisto=no-hit and ours=ok")
    parser.add_argument("--ok-no-hits", required=True, help="read IDs where kallisto=ok and ours=no-hit")
    parser.add_argument("--out", required=True, help="output TSV path")
    args = parser.parse_args()

    no_hits_ok = read_id_set(args.no_hits_ok)
    ok_no_hits = read_id_set(args.ok_no_hits)
    all_ids = sorted(no_hits_ok | ok_no_hits)

    kallisto_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open(args.kmer_dump, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read_name"])
            if rid not in no_hits_ok and rid not in ok_no_hits:
                continue
            row["read_pos"] = str(int(row["read_pos"]))
            kallisto_rows[rid].append(row)

    ours_by_pos: dict[str, dict[int, dict[str, int]]] = defaultdict(
        lambda: defaultdict(
            lambda: {"matched": 0, "mphf_hit": 0, "special": 0, "overcrowded": 0}
        )
    )
    with open(args.candidates, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read"])
            pos = int(row["read_pos"])
            slot = ours_by_pos[rid][pos]
            if row["matched"] == "1":
                slot["matched"] = 1
            if row["mphf_hit"] == "1":
                slot["mphf_hit"] = 1
            if row["has_special"] == "1":
                slot["special"] = 1
            if row["overcrowded"] == "1":
                slot["overcrowded"] = 1

    visited: dict[str, set[int]] = {}
    with open(args.visited, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read"])
            raw = row["positions_visited"].strip()
            visited[rid] = set(int(v) for v in raw.split(",") if v)

    reason_by_read: dict[str, str] = {}
    with open(args.trace, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reason_by_read[norm_read_id(row["read"])] = row["reason"]

    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "read",
                "class",
                "trace_reason",
                "divergence_pos",
                "k_hit_action",
                "k_skip_reason",
                "k_is_hit",
                "k_special",
                "k_overcrowded",
                "k_jump_used",
                "k_jump_from",
                "k_jump_to",
                "ours_visited",
                "ours_matched",
                "ours_mphf_hit",
                "ours_special",
                "ours_overcrowded",
            ]
        )

        for rid in all_ids:
            cls = "no_hits_ok" if rid in no_hits_ok else "ok_no_hits"
            krows = sorted(kallisto_rows[rid], key=lambda r: int(r["read_pos"]))
            ours = ours_by_pos[rid]
            vis = visited.get(rid, set())
            divergence_pos = None
            krow: dict[str, str] | None = None

            if cls == "ok_no_hits":
                for row in krows:
                    if row.get("hit_action") != "accept":
                        continue
                    pos = int(row["read_pos"])
                    if ours[pos]["matched"] == 0:
                        divergence_pos = pos
                        krow = row
                        break
            else:
                for pos in sorted(p for p, state in ours.items() if state["matched"] == 1):
                    same_pos = [row for row in krows if int(row["read_pos"]) == pos]
                    if any(row.get("hit_action") == "accept" for row in same_pos):
                        continue
                    divergence_pos = pos
                    krow = same_pos[0] if same_pos else None
                    break

            if divergence_pos is None:
                if cls == "ok_no_hits":
                    for row in krows:
                        if row.get("hit_action") == "accept":
                            divergence_pos = int(row["read_pos"])
                            krow = row
                            break
                else:
                    for pos in sorted(ours):
                        if ours[pos]["matched"] == 1:
                            divergence_pos = pos
                            break

            if krow is None:
                krow = {
                    "hit_action": "-",
                    "skip_reason": "-",
                    "is_hit": "-",
                    "minimizer_is_special": "-",
                    "minimizer_is_overcrowded": "-",
                    "jump_used": "-",
                    "jump_from_pos": "-",
                    "jump_to_pos": "-",
                }

            if divergence_pos is None:
                ours_state = {"matched": 0, "mphf_hit": 0, "special": 0, "overcrowded": 0}
            else:
                ours_state = ours[divergence_pos]

            writer.writerow(
                [
                    rid,
                    cls,
                    reason_by_read.get(rid, "-"),
                    divergence_pos if divergence_pos is not None else "-",
                    krow.get("hit_action", "-"),
                    krow.get("skip_reason", "-"),
                    krow.get("is_hit", "-"),
                    krow.get("minimizer_is_special", "-"),
                    krow.get("minimizer_is_overcrowded", "-"),
                    krow.get("jump_used", "-"),
                    krow.get("jump_from_pos", "-"),
                    krow.get("jump_to_pos", "-"),
                    int(divergence_pos in vis) if divergence_pos is not None else 0,
                    ours_state["matched"],
                    ours_state["mphf_hit"],
                    ours_state["special"],
                    ours_state["overcrowded"],
                ]
            )


if __name__ == "__main__":
    main()
