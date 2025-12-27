#!/usr/bin/env python3
"""Generate a small synthetic transcriptome and reads for quick kallisto parity checks."""

import argparse
import os
import random
from typing import List

BASES = "ACGT"


def rand_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def mutate_seq(rng: random.Random, seq: str, error_rate: float) -> str:
    if error_rate <= 0:
        return seq
    out = list(seq)
    for i, c in enumerate(out):
        if rng.random() < error_rate:
            choices = [b for b in BASES if b != c]
            out[i] = rng.choice(choices)
    return "".join(out)


def write_fasta(path: str, names: List[str], seqs: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for name, seq in zip(names, seqs):
            handle.write(f">{name}\n")
            handle.write(seq)
            handle.write("\n")


def write_fastq(path: str, reads: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for i, seq in enumerate(reads):
            handle.write(f"@read_{i}\n")
            handle.write(seq)
            handle.write("\n+\n")
            handle.write("I" * len(seq))
            handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="data", help="Output directory (default: data)")
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--num-transcripts", type=int, default=4)
    parser.add_argument("--transcript-length", type=int, default=320)
    parser.add_argument("--shared-length", type=int, default=60)
    parser.add_argument("--read-length", type=int, default=80)
    parser.add_argument("--num-reads", type=int, default=200)
    parser.add_argument("--num-random", type=int, default=20)
    parser.add_argument("--error-rate", type=float, default=0.01)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)

    shared = rand_seq(rng, args.shared_length)
    names = [f"tx{i}" for i in range(args.num_transcripts)]
    seqs = []
    for _ in names:
        seq = list(rand_seq(rng, args.transcript_length))
        if args.shared_length > 0 and args.shared_length < args.transcript_length:
            insert_at = rng.randint(10, args.transcript_length - args.shared_length - 10)
            seq[insert_at : insert_at + args.shared_length] = list(shared)
        seqs.append("".join(seq))

    fasta_path = os.path.join(args.out_dir, "simple_transcripts.fa")
    write_fasta(fasta_path, names, seqs)

    reads = []
    for _ in range(args.num_reads):
        t_idx = rng.randrange(len(seqs))
        seq = seqs[t_idx]
        if args.read_length >= len(seq):
            start = 0
        else:
            start = rng.randint(0, len(seq) - args.read_length)
        read = seq[start : start + args.read_length]
        reads.append(mutate_seq(rng, read, args.error_rate))

    for _ in range(args.num_random):
        reads.append(rand_seq(rng, args.read_length))

    fastq_path = os.path.join(args.out_dir, "simple_reads.fq")
    write_fastq(fastq_path, reads)

    print(f"wrote {fasta_path}")
    print(f"wrote {fastq_path}")


if __name__ == "__main__":
    main()
