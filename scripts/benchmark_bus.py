#!/usr/bin/env python3
"""Benchmark complete BUS generation and compare exact barcode/UMI/EC counts."""

import argparse
import gzip
import json
from pathlib import Path
import struct
import subprocess

import numpy as np

from benchmark_common import digest, executable, fingerprint, run_process, summarize, write_json

RECORD = np.dtype([('barcode', '<u8'), ('umi', '<u8'), ('ec', '<i4'),
                   ('count', '<u4'), ('flags', '<u4'), ('padding', '<u4')])


def canonical_bus(directory, out, classes, transcripts, bustools, threads):
    """Normalize complete transcript sets, never just EC numbers or marginal totals."""
    names = (directory / 'transcripts.txt').read_text().splitlines()
    if names != transcripts:
        raise ValueError('Transcript lists differ for the shared index')
    mapping = {}
    for line in (directory / 'matrix.ec').read_text().splitlines():
        number, targets = line.split('\t')
        number = int(number)
        targets = tuple(sorted(map(int, targets.split(','))))
        if number < 0 or number in mapping or not targets or len(set(targets)) != len(targets):
            raise ValueError('Invalid EC definition')
        if targets[0] < 0 or targets[-1] >= len(names):
            raise ValueError('EC target outside transcript list')
        mapping[number] = classes.setdefault(targets, len(classes))
    translation = np.full(max(mapping, default=-1) + 1, -1, dtype='<i4')
    for key, value in mapping.items():
        translation[key] = value
    temporary = out.with_suffix('.unsorted.bus')
    with (directory / 'output.bus').open('rb') as src, temporary.open('wb') as dst:
        header = src.read(20)
        magic, version, barcode_len, umi_len, text_len = struct.unpack('<4sIIII', header)
        if magic != b'BUS\0' or version != 1 or barcode_len != 16 or umi_len != 12:
            raise ValueError('Expected BUS v1 with 10xv3 barcode/UMI lengths')
        if len(src.read(text_len)) != text_len:
            raise ValueError('Truncated BUS header')
        dst.write(struct.pack('<4sIIII', magic, version, barcode_len, umi_len, 0))
        while block := src.read(RECORD.itemsize * 1_000_000):
            if len(block) % RECORD.itemsize:
                raise ValueError('Truncated BUS record')
            rows = np.frombuffer(block, dtype=RECORD).copy()
            if np.any(rows['ec'] < 0) or np.any(rows['ec'] >= len(translation)) or np.any(rows['count'] == 0):
                raise ValueError('Invalid BUS EC or count')
            rows['ec'] = translation[rows['ec']]
            if np.any(rows['ec'] < 0):
                raise ValueError('BUS references an undefined EC')
            rows['padding'] = 0
            dst.write(rows.tobytes())
    # Default bustools sorting retains flags and sums identical weighted records.
    subprocess.run([str(bustools), 'sort', '-t', str(threads), '-m', '2G',
                    '-o', str(out), str(temporary)], check=True)
    temporary.unlink()
    count, records = 0, 0
    with out.open('rb') as src:
        src.seek(20)
        while block := src.read(RECORD.itemsize * 1_000_000):
            rows = np.frombuffer(block, dtype=RECORD)
            count += int(rows['count'].sum(dtype=np.uint64))
            records += len(rows)
    info = json.loads((directory / 'run_info.json').read_text())
    if info['n_pseudoaligned'] != count or info['n_targets'] != len(transcripts):
        raise ValueError('BUS records disagree with run_info counts')
    return {'sha256': digest(out), 'mapped_reads': count, 'distinct_records': records,
            'run_info': info}


def compare_outputs(outputs, root, bustools, threads):
    root.mkdir(parents=True, exist_ok=True)
    classes = {}
    transcripts = (outputs[0][1] / 'transcripts.txt').read_text().splitlines()
    evidence = {}
    for label, directory in outputs:
        print('Validating complete BUS counts:', label, flush=True)
        canonical = root / (label + '.canonical.bus')
        evidence[label] = canonical_bus(directory, canonical, classes, transcripts, bustools, threads)
        if evidence[label]['sha256'] != next(iter(evidence.values()))['sha256']:
            write_json(root / 'parity-failure.json', evidence)
            raise ValueError(f'Exact barcode/UMI/transcript-set/flag/count mismatch: {label}')
    ec = root / 'canonical.matrix.ec'
    with ec.open('w') as handle:
        for targets, number in classes.items():
            handle.write(f'{number}\t' + ','.join(map(str, targets)) + '\n')
    write_json(root / 'parity.json', evidence)
    return evidence, ec


def prefix_reads(reads, out, pairs):
    result = []
    for mate, path in enumerate(reads[:2], 1):
        dest = out / f'preflight-{mate}.fastq.gz'
        with gzip.open(path, 'rb') as src, gzip.open(dest, 'wb', compresslevel=1) as dst:
            for _ in range(pairs * 4):
                line = src.readline()
                if not line:
                    raise ValueError('Not enough records for preflight')
                dst.write(line)
        result.append(dest)
    return result


def gene_map(fasta, transcripts, out):
    wanted = set(transcripts.read_text().splitlines())
    found = {}
    opener = gzip.open if fasta.suffix == '.gz' else open
    with opener(fasta, 'rt') as src:
        for line in src:
            if line.startswith('>'):
                name = line[1:].split()[0]
                fields = name.split('|')
                if len(fields) < 2:
                    raise ValueError('Expected GENCODE pipe-delimited transcript/gene headers')
                for alias in (name, fields[0]):
                    if alias in wanted:
                        found[alias] = fields[1]
    if set(found) != wanted:
        raise ValueError('Reference FASTA does not map every indexed transcript to a gene')
    with out.open('w') as dst:
        for name in transcripts.read_text().splitlines():
            dst.write(f'{name}\t{found[name]}\n')


def matrices(outputs, root, bustools, threads, onlist, fasta, ec):
    t2g = root / 'transcripts-to-genes.tsv'
    gene_map(fasta, outputs[0][1] / 'transcripts.txt', t2g)
    result = {}
    for label, directory in outputs:
        prefix = root / (label + '.gene')
        corrected = root / (label + '.corrected.unsorted.bus')
        sorted_bus = root / (label + '.corrected.bus')
        commands = [
            [bustools, 'correct', '-w', onlist, '-o', corrected, root / (label + '.canonical.bus')],
            [bustools, 'sort', '-t', str(threads), '-m', '2G', '-o', sorted_bus, corrected],
            [bustools, 'count', '--genecounts', '-g', t2g, '-e', ec,
             '-t', directory / 'transcripts.txt', '-o', prefix, sorted_bus],
        ]
        for command in commands:
            subprocess.run(list(map(str, command)), check=True)
        corrected.unlink()
        files = {suffix: fingerprint(Path(str(prefix) + suffix))
                 for suffix in ('.mtx', '.genes.txt', '.barcodes.txt')}
        result[label] = files
        if {key: value['sha256'] for key, value in files.items()} != {
                key: value['sha256'] for key, value in next(iter(result.values())).items()}:
            raise ValueError('Corrected gene/UMI matrices differ')
    write_json(root / 'matrix-parity.json', result)
    # Exact matrix and label equality proves zero error, including implicit zeros.
    matrix_path = Path(next(iter(result.values()))['.mtx']['path'])
    with matrix_path.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith('%'):
                rows, columns, nonzero = map(int, line.split())
                break
        else:
            raise ValueError('Missing sparse matrix dimensions')
    zeros = dict.fromkeys(('median', 'p90', 'p99', 'p99_9', 'p99_99', 'max'), 0)
    comparison = {'exact': True, 'shape': [rows, columns], 'nonzero_entries': nonzero,
                  'absolute_errors_all_entries': zeros,
                  'absolute_errors_union_nonzero_entries': zeros,
                  'percentile_method': 'all absolute differences are exactly zero by complete matrix/label byte equality'}
    write_json(root / 'matrix-differences.json', comparison)
    return {'files': result, 'differences': comparison}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kallistors', default='target/release/kallistors')
    parser.add_argument('--kallisto', default='target/tools/kallisto-0.52.0')
    parser.add_argument('--bustools', required=True)
    parser.add_argument('--index', type=Path, required=True)
    parser.add_argument('--reads', type=Path, nargs='+', required=True, help='R1 R2 for each lane')
    parser.add_argument('--preparation', type=Path, required=True)
    parser.add_argument('--onlist', type=Path)
    parser.add_argument('--reference-fasta', type=Path)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--preflight-pairs', type=int, default=100000)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if (len(args.reads) % 2 or args.threads < 1 or args.runs < 1 or args.warmups < 0
            or args.preflight_pairs < 1 or bool(args.onlist) != bool(args.reference_fasta)):
        parser.error('Supply paired lanes, positive sizes and both --onlist/--reference-fasta or neither')
    if args.out.exists() and any(args.out.iterdir()):
        parser.error('--out must be new or empty')
    args.out.mkdir(parents=True, exist_ok=True)
    binaries = {tool: executable(getattr(args, tool)) for tool in ('kallisto', 'kallistors')}
    bustools = executable(args.bustools)
    preparation = json.loads(args.preparation.read_text())
    reads = [fingerprint(path) for path in args.reads]
    expected = [lane[key] for lane in preparation['lanes'] for key in ('read1', 'read2')]
    if reads != expected:
        raise ValueError('Inputs do not match complete prepared paired-library provenance')
    provenance = {'executables': {key: fingerprint(path) for key, path in binaries.items()},
                  'bustools': fingerprint(bustools), 'index': fingerprint(args.index),
                  'reads': reads, 'preparation': preparation,
                  'onlist': fingerprint(args.onlist) if args.onlist else None,
                  'reference': fingerprint(args.reference_fasta) if args.reference_fasta else None,
                  'technology': '10xv3', 'threads': args.threads,
                  'warmups': args.warmups, 'runs': args.runs,
                  'scope': 'Complete BUS generation, including index loading and output; preprocessing and bustools excluded'}
    write_json(args.out / 'provenance.json', provenance)

    def command(tool, directory, inputs):
        return [binaries[tool], 'bus', '-x', '10xv3', '-t', str(args.threads),
                '-i', args.index, '-o', directory, *inputs]

    preflight = args.out / 'preflight'
    preflight.mkdir()
    inputs = prefix_reads(args.reads, preflight, args.preflight_pairs)
    preflight_outputs = []
    for tool in binaries:
        directory = preflight / tool
        run_process(command(tool, directory, inputs), preflight / (tool + '-run'))
        preflight_outputs.append((tool, directory))
    preflight_proof, _ = compare_outputs(preflight_outputs, preflight, bustools, args.threads)
    if any(item['run_info']['n_processed'] != args.preflight_pairs for item in preflight_proof.values()):
        raise ValueError('Incomplete preflight processing')
    print('Complete preflight parity passed; beginning full-library timing', flush=True)
    rows, outputs = [], []
    iteration = 0
    for phase, repetitions in [('warmup', args.warmups), ('measured', args.runs)]:
        for replicate in range(1, repetitions + 1):
            order = list(binaries) if iteration % 2 == 0 else list(reversed(binaries))
            iteration += 1
            for tool in order:
                label = f'{phase}-{replicate}-{tool}'
                directory = args.out / label
                print('Running', label, flush=True)
                row = run_process(command(tool, directory, args.reads), args.out / (label + '-run'))
                row.update(tool=tool, phase=phase, replicate=replicate)
                rows.append(row)
                outputs.append((label, directory))
                write_json(args.out / 'runs.json', rows)
                print(f'{label}: {row["wall_seconds"]:.3f} s; peak RSS {row["peak_rss_bytes"] / 2**30:.3f} GiB', flush=True)
    # Freeze checks outside timing, before validating every retained full-library output.
    for tool, binary in binaries.items():
        if fingerprint(binary) != provenance['executables'][tool]:
            raise ValueError(f'Executable changed during benchmark: {tool}')
    if fingerprint(args.index) != provenance['index'] or [fingerprint(p) for p in args.reads] != reads:
        raise ValueError('Benchmark inputs changed during execution')
    if args.onlist and (fingerprint(args.onlist) != provenance['onlist']
                        or fingerprint(args.reference_fasta) != provenance['reference']):
        raise ValueError('Barcode/reference inputs changed during execution')
    parity, ec = compare_outputs(outputs, args.out, bustools, args.threads)
    for label, item in parity.items():
        if item['run_info']['n_processed'] != preparation['total_pairs']:
            raise ValueError(f'Incomplete library processing: {label}')
    matrix = None
    if args.onlist:
        selected = [(label, directory) for label, directory in outputs if label.startswith('measured-1-')]
        matrix = matrices(selected, args.out, bustools, args.threads, args.onlist, args.reference_fasta, ec)
    result = {**summarize(rows), 'exact_bus_parity': True, 'parity': parity, 'matrix_parity': matrix}
    write_json(args.out / 'summary.json', result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
