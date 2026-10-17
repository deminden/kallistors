#!/usr/bin/env python3
"""Select an exact, seeded uniform subset of paired FASTQ records after trimming."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import random
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', type=Path, required=True)
    parser.add_argument('--reads2', type=Path, required=True)
    parser.add_argument('--total-pairs', type=int, required=True)
    parser.add_argument('--pairs', type=int, default=4_000_000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.pairs <= args.total_pairs:
        parser.error('require 0 < pairs <= total-pairs')
    if args.out.exists() and any(args.out.iterdir()):
        parser.error('--out must be empty/new')
    args.out.mkdir(parents=True, exist_ok=True)
    destinations = [args.out / f'mate{mate}.fastq.gz' for mate in (1, 2)]
    rng = random.Random(args.seed)
    needed = args.pairs
    start = time.perf_counter()
    with gzip.open(args.reads, 'rb') as left, gzip.open(args.reads2, 'rb') as right, \
            destinations[0].open('wb') as raw1, destinations[1].open('wb') as raw2, \
            gzip.GzipFile(filename='', mode='wb', fileobj=raw1, compresslevel=4, mtime=0) as out1, \
            gzip.GzipFile(filename='', mode='wb', fileobj=raw2, compresslevel=4, mtime=0) as out2:
        for ordinal in range(args.total_pairs):
            records = [[source.readline() for _ in range(4)] for source in (left, right)]
            if any(not record[-1] or not record[0].startswith(b'@') or not record[0][1:].strip()
                   or not record[2].startswith(b'+')
                   or len(record[1].rstrip(b'\r\n')) != len(record[3].rstrip(b'\r\n'))
                   for record in records):
                raise RuntimeError(f'Truncated/invalid FASTQ pair {ordinal}')
            names = [record[0][1:].split(None, 1)[0] for record in records]
            names = [name[:-2] if name.endswith((b'/1', b'/2')) else name for name in names]
            if not names[0] or names[0] != names[1]:
                raise RuntimeError(f'FASTQ mate names differ at pair {ordinal}')
            remaining = args.total_pairs - ordinal
            if rng.randrange(remaining) < needed:
                out1.writelines(records[0])
                out2.writelines(records[1])
                needed -= 1
        if left.read(1) or right.read(1) or needed:
            raise RuntimeError('Input pair count does not match --total-pairs')
    files = []
    for path in destinations:
        h = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
                h.update(block)
        files.append({'path': str(path.resolve()), 'sha256': h.hexdigest(), 'bytes': path.stat().st_size})
    result = {'inputs': [str(args.reads.resolve()), str(args.reads2.resolve())],
              'total_pairs': args.total_pairs, 'selected_pairs': args.pairs, 'seed': args.seed,
              'method': 'sequential selection without replacement: randrange(remaining) < still_needed',
              'files': files, 'elapsed_seconds': time.perf_counter() - start}
    (args.out / 'subset.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
