#!/usr/bin/env python3
"""Compare index build time, file size and optional peak RSS with alternating controls."""

import argparse
from pathlib import Path
import sys

from benchmark_common import digest, executable, fingerprint, run_process, summarize, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kallisto', default='target/tools/kallisto-0.52.0')
    parser.add_argument('--kallistors', default='target/release/kallistors')
    parser.add_argument('--fasta', type=Path, default=Path('data/srr-small/gencode.v50.transcripts.fa.gz'))
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--memory', action=argparse.BooleanOptionalAction, default=sys.platform in ('darwin', 'linux'))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.threads < 1 or args.runs < 1 or args.warmups < 0:
            raise ValueError('Require positive threads/runs and nonnegative warmups')
        if args.out.exists() and any(args.out.iterdir()):
            raise ValueError('--out must be new or empty')
        args.out = args.out.resolve()
        tools = {key: executable(getattr(args, key)) for key in ('kallisto', 'kallistors')}
        config = {'fasta': fingerprint(args.fasta), 'binaries': {key: fingerprint(p) for key, p in tools.items()},
                  'threads': args.threads, 'runs': args.runs, 'warmups': args.warmups, 'memory': args.memory,
                  'scripts': {p.name: digest(p) for p in [Path(__file__), Path(__file__).with_name('benchmark_common.py')]}}
        write_json(args.out / 'protocol.json', {'schema_version': 1, 'kind': 'index', 'configuration': config})
        rows = []
        for phase, count in [('warmup', args.warmups), ('measured', args.runs)]:
            for rep in range(1, count + 1):
                for tool in list(tools) if rep % 2 else list(reversed(tools)):
                    if fingerprint(tools[tool]) != config['binaries'][tool]:
                        raise ValueError('An executable changed during measurement')
                    prefix = args.out / f'{phase}-{rep}-{tool}'
                    index = prefix.with_suffix('.idx')
                    row = run_process([str(tools[tool]), 'index', '-t', str(args.threads),
                                       '-i', str(index), config['fasta']['path']], prefix, args.memory)
                    if not index.is_file() or not index.stat().st_size:
                        raise ValueError(f'{tool} did not write an index')
                    row.update(tool=tool, phase=phase, replicate=rep, index_bytes=index.stat().st_size)
                    rows.append(row)
                    write_json(args.out / 'runs.json', rows)
                    index.unlink()
                    print(f'{phase} {rep} {tool}: {row["wall_seconds"]:.3f}s', flush=True)
        if fingerprint(args.fasta) != config['fasta'] or any(fingerprint(p) != config['binaries'][key] for key, p in tools.items()):
            raise ValueError('An input or executable changed during measurement')
        write_json(args.out / 'summary.json', {'schema_version': 1, 'kind': 'index', 'complete': True, **summarize(rows)})
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f'Index benchmark failed: {error}\n')


if __name__ == '__main__':
    main()
