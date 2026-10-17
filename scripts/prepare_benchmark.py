#!/usr/bin/env python3
"""Retrieve checksum-pinned RNA-seq libraries and prepare paired reads with fasterp."""

import argparse
import concurrent.futures
import json
from pathlib import Path
import shutil
import subprocess
import time

from benchmark_common import DATASETS, datasets, digest, fetch

FASTERP_VERSION = '0.2.1'


def prepare(row, data, binary, trimmer, threads):
    run = row['run_accession']
    root = data / row.get('data_directory', run.lower())
    root.mkdir(parents=True, exist_ok=True)
    (root / 'ena-metadata.json').write_text(json.dumps(row, indent=2) + '\n')
    urls, hashes = row['fastq_ftp'].split(';'), row['fastq_md5'].split(';')
    sizes = row['fastq_bytes'].split(';') if row.get('fastq_bytes') else [None, None]
    if len(urls) != 2 or len(hashes) != 2 or len(sizes) != 2:
        raise RuntimeError(f'Manifest must supply two mates: {run}')
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        files = list(pool.map(lambda item: fetch(item[0], root / Path(item[0]).name, item[1], item[2]),
                              zip(urls, hashes, sizes)))
    trimmed = root / 'trimmed-fasterp'
    outputs = [trimmed / f'{run}_{mate}.trimmed.fastq.gz' for mate in (1, 2)]
    command = [binary, '-i', files[0]['path'], '-I', files[1]['path'],
               '-o', str(outputs[0]), '-O', str(outputs[1]),
               '--cut-front', '--cut-tail', '--cut-mean-quality', '20',
               '--cut-window-size', '4', '--length-required', '50',
               '--qualified-quality-phred', '15', '--unqualified-percent-limit', '40',
               '--n-base-limit', '5', '--compression-level', '4', '--threads', str(threads),
               '--json', str(trimmed / 'fasterp.json'), '--html', str(trimmed / 'fasterp.html')]
    record = trimmed / 'preparation.json'
    if record.exists():
        saved = json.loads(record.read_text())
        if saved['command'] != command or saved['trimmer'] != trimmer:
            raise RuntimeError(f'Trimming settings changed; use a fresh --data directory: {record}')
        if any(not Path(item['path']).is_file() or digest(Path(item['path'])) != item['sha256']
               for item in saved['outputs']):
            raise RuntimeError(f'Trimmed output is missing or changed: {record}')
        return saved
    if trimmed.exists() and any(trimmed.iterdir()):
        raise RuntimeError(f'Unrecorded trimming outputs exist: {trimmed}')
    trimmed.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with (trimmed / 'fasterp.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)
    report = json.loads((trimmed / 'fasterp.json').read_text())
    result = {'files': files, 'trimmer': trimmer, 'command': command,
              'trim_seconds': time.perf_counter() - start,
              'pairs_before': report['summary']['before_filtering']['total_reads'] // 2,
              'pairs_after': report['summary']['after_filtering']['total_reads'] // 2,
              'outputs': [{'path': str(path), 'bytes': path.stat().st_size,
                           'sha256': digest(path)} for path in outputs]}
    record.write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=DATASETS)
    parser.add_argument('--samples', default='all', help='comma-separated ENA accessions, or all')
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--fasterp', default='fasterp', help='fasterp executable or path')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    available = datasets(args.manifest)
    selected = list(available) if args.samples == 'all' else args.samples.split(',')
    if not selected or len(set(selected)) != len(selected) or set(selected) - available.keys():
        parser.error('select unique accessions from the dataset manifest, or all')
    if args.threads < 1:
        parser.error('--threads must be positive')
    binary = shutil.which(args.fasterp)
    if not binary and args.fasterp == 'fasterp':
        tools = Path(__file__).resolve().parents[1] / 'target/tools'
        for name in ('fasterp-fixed', 'fasterp'):
            binary = shutil.which('fasterp', path=str(tools / name / 'bin'))
            if binary:
                break
    if not binary:
        parser.error('fasterp is missing; see docs/benchmarks.md#retrieve-and-trim-reads for the pinned build')
    binary = str(Path(binary).resolve())
    version = subprocess.check_output([binary, '--version'], text=True).strip()
    if version != 'fasterp ' + FASTERP_VERSION:
        parser.error(f'expected fasterp {FASTERP_VERSION}, got {version!r}')
    trimmer = {'name': 'fasterp', 'version': FASTERP_VERSION, 'binary_sha256': digest(Path(binary))}
    data = args.data.resolve()
    for run in selected:
        result = prepare(available[run], data, binary, trimmer, args.threads)
        print(json.dumps({'run': run, **result}, indent=2), flush=True)


if __name__ == '__main__':
    main()
