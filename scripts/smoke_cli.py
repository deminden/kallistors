#!/usr/bin/env python3
"""Exercise a built binary on a tiny, self-contained index and quantification workflow."""

import argparse
import csv
import json
from pathlib import Path
import random
import subprocess
import tempfile


def smoke(binary):
    binary = Path(binary).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix='kallistors-smoke-') as directory:
        work = Path(directory)

        def run(*arguments):
            command = [str(binary), *map(str, arguments)]
            result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=60)
            if result.returncode:
                raise RuntimeError(f'{command} failed:\n{result.stdout}\n{result.stderr}')
            return result.stdout

        run('version')
        run('--help')
        if '10xv3' not in run('bus', '--list').lower():
            raise RuntimeError('Technology list is missing 10xv3')

        rng = random.Random(7)
        transcripts = [''.join(rng.choice('ACGT') for _ in range(240)) for _ in range(2)]
        (work / 'transcripts.fa').write_text(''.join(
            f'>tx{i}\n{sequence}\n' for i, sequence in enumerate(transcripts)), encoding='ascii')
        reads = [transcripts[0][20:80], transcripts[0][20:80], transcripts[1][60:120]]
        (work / 'reads.fq').write_text(''.join(
            f'@read{i}\n{sequence}\n+\n{"I" * len(sequence)}\n'
            for i, sequence in enumerate(reads)), encoding='ascii')
        run('index', '-i', 'transcripts.idx', '-t', '2', 'transcripts.fa')
        if 'transcripts: 2' not in run('index-info', '--index', 'transcripts.idx'):
            raise RuntimeError('Index inspection did not retain both transcripts')
        run('pseudoalign', '--index', 'transcripts.idx', '--reads', 'reads.fq',
            '--out', 'ecs.tsv', '--threads', '2', '--match-cache-size', '1MiB')
        with (work / 'ecs.tsv').open() as handle:
            counts = {row['transcripts']: int(row['count']) for row in csv.DictReader(handle, delimiter='\t')}
        if counts != {'0': 2, '1': 1}:
            raise RuntimeError(f'Unexpected mapped EC counts: {counts}')
        run('quant', '-i', 'transcripts.idx', '-o', 'quant', '-t', '2',
            '--single', '-l', '100', '-s', '10', '--match-cache-size', '1MiB', 'reads.fq')
        info = json.loads((work / 'quant/run_info.json').read_text())
        if info['n_processed'] != 3 or info['n_pseudoaligned'] != 3:
            raise RuntimeError(f'Unexpected quantification counts: {info}')
        with (work / 'quant/abundance.tsv').open() as handle:
            abundance = {row['target_id']: float(row['est_counts'])
                         for row in csv.DictReader(handle, delimiter='\t')}
        if abundance != {'tx0': 2.0, 'tx1': 1.0}:
            raise RuntimeError(f'Unexpected transcript abundances: {abundance}')
        with (work / 'quant/abundance.h5').open('rb') as handle:
            if handle.read(8) != b'\x89HDF\r\n\x1a\n':
                raise RuntimeError('Default quantification did not produce HDF5 output')
    print('Release smoke check passed: index, metadata, exact EC counts, quantification and HDF5 output.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    smoke(parser.parse_args().binary)
