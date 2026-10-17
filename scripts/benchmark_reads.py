#!/usr/bin/env python3
"""Compare quantification with mandatory complete EC preflight and optional phase timing."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from benchmark_common import (DATASETS, ORIGINAL_MARKERS, RUST_MARKERS, datasets, digest,
                              executable, fingerprint, prepared_reads, run_process,
                              summarize, timing_environment, write_json)
from compare_ec_counts import compare, load_counts


def configurations(args, rows):
    if args.suite:
        if args.reads or args.reads2 or args.samples != 'all' or args.index or args.reference_metadata:
            raise ValueError('--suite defines samples and references; do not combine with direct inputs')
        specs = json.loads(args.suite.read_text())['cases']
        if args.group:
            specs = [spec for spec in specs if spec.get('group') == args.group]
    elif args.reads:
        if args.samples != 'all':
            raise ValueError('Choose --samples or --reads')
        specs = [{'name': 'custom', 'reads': str(args.reads), 'reads2': str(args.reads2) if args.reads2 else None,
                  'index': str((args.index or args.data / 'srr-small/gencode.v50.kallisto.idx').resolve()),
                  'reference_metadata': str(args.reference_metadata.resolve()) if args.reference_metadata else None}]
    else:
        selected = list(rows) if args.samples == 'all' else args.samples.split(',')
        if len(set(selected)) != len(selected) or not selected or set(selected) - rows.keys():
            raise ValueError('Select unique accessions from the dataset manifest, or all')
        specs = [{'name': name, 'sample': name,
                  'index': str((args.index or args.data / 'srr-small/gencode.v50.kallisto.idx').resolve()),
                  'reference_metadata': str(args.reference_metadata.resolve()) if args.reference_metadata else None}
                 for name in selected]
    if args.group and not args.suite:
        raise ValueError('--group requires --suite')
    configs = []
    for spec in specs:
        unknown = set(spec) - {'name', 'group', 'sample', 'subset', 'reads', 'reads2', 'index', 'reference_metadata', 'modes'}
        if unknown:
            raise ValueError(f'Unknown suite fields: {sorted(unknown)}')
        if type(spec.get('subset', False)) is not bool:
            raise ValueError('Suite subset must be a boolean')
        name = spec['name']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
            raise ValueError(f'Invalid case name: {name}')
        provenance = None
        if 'sample' in spec:
            row = rows[spec['sample']]
            reads = prepared_reads(row, args.data, spec.get('subset', False))
            record = reads[0].parent / ('subset.json' if spec.get('subset') else 'preparation.json')
            provenance = {'dataset': row, 'preparation': json.loads(record.read_text())}
        else:
            reads = [Path(spec[key]) for key in ('reads', 'reads2') if spec.get(key)]
        index = Path(spec['index'])
        if not index.is_absolute():
            index = args.data / index
        reference = spec.get('reference_metadata')
        if reference:
            reference = Path(reference)
            if not reference.is_absolute():
                reference = args.data / reference
            reference = json.loads(reference.read_text())
            if 'total_targets' not in reference:
                reference['total_targets'] = reference.get('mature_targets', 0) + reference.get('nascent_targets', 0)
        modes = spec.get('modes', args.modes.split(','))
        if not modes or len(set(modes)) != len(modes) or set(modes) - {'paired', 'single'}:
            raise ValueError(f'Invalid modes for {name}')
        for mode in modes:
            if len(reads) != 2 and mode == 'paired' or not reads:
                raise ValueError(f'{name}: paired mode requires two FASTQs')
            for threads in args.thread_counts:
                configs.append({'id': f'{name}--{mode}--t{threads}', 'name': name,
                                'mode': mode, 'threads': threads, 'index': str(index.resolve()),
                                'reads': [str(p.resolve()) for p in reads[:2 if mode == 'paired' else 1]],
                                'fragment_length': args.fragment_length, 'fragment_sd': args.fragment_sd,
                                'reference': reference, 'provenance': provenance})
    if not configs or len({c['id'] for c in configs}) != len(configs):
        raise ValueError('Empty or duplicate benchmark configurations')
    return configs


def quant_command(binary, config, out):
    command = [str(binary), 'quant', '-t', str(config['threads']), '-i', config['index'], '-o', str(out)]
    if config['mode'] == 'single':
        command += ['--single', '-l', str(config['fragment_length']), '-s', str(config['fragment_sd'])]
    return command + config['reads']


def ec_preflight(config, args, tools):
    root = args.out / 'correctness' / config['id']
    report = root / 'ec.json'
    original, candidate = root / 'original/diagnostic.ec_counts.tsv', root / 'candidate.ec.tsv'
    if report.exists():
        saved = json.loads(report.read_text())
        if saved['files'] != {'original': fingerprint(original), 'candidate': fingerprint(candidate)}:
            raise ValueError(f'Changed EC evidence: {root}')
        result = compare(load_counts(original), load_counts(candidate))
        if not result['exact'] or result != saved['comparison']:
            raise ValueError(f'Invalid EC evidence: {root}')
        return saved
    if root.exists():
        archive = args.out / 'incomplete'
        archive.mkdir(exist_ok=True)
        shutil.move(str(root), str(archive / f'{config["id"]}-{len(list(archive.iterdir()))}'))
    root.mkdir(parents=True)
    env = timing_environment()
    env.update(KALLISTO_DIAGNOSTIC_ECS='1', KALLISTO_DIAGNOSTIC_MAPPING_ONLY='1')
    command = quant_command(args.diagnostic_kallisto, config, root / 'original')
    with (root / 'original.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=log, env=env, check=True)
    command = [str(tools['kallistors']), 'pseudoalign', '--threads', str(config['threads']),
               '--index', config['index'], '--reads', config['reads'][0], '--out', str(candidate)]
    if config['mode'] == 'paired':
        command += ['--reads2', config['reads'][1]]
    else:
        command += ['--fragment-length', str(config['fragment_length'])]
    with (root / 'candidate.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=log, env=timing_environment(), check=True)
    result = compare(load_counts(original), load_counts(candidate))
    saved = {'comparison': result, 'files': {'original': fingerprint(original), 'candidate': fingerprint(candidate)}}
    write_json(report, saved)
    if not result['exact']:
        raise RuntimeError(f'Complete EC parity failed: {config["id"]}; see {report}')
    print(f'{config["id"]}: complete EC parity passed', flush=True)
    return saved


def archive_incomplete(prefix, output):
    files = [output, prefix.with_suffix('.log'), prefix.with_suffix('.memory.log')]
    existing = [p for p in files if p.exists()]
    if existing:
        archive = prefix.parent / 'incomplete'
        archive.mkdir(exist_ok=True)
        archive = archive / f'{prefix.name}-{len(list(archive.iterdir()))}'
        archive.mkdir()
        for path in existing:
            shutil.move(str(path), str(archive / path.name))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kallisto', default='target/tools/kallisto-0.52.0')
    parser.add_argument('--kallistors', default='target/release/kallistors')
    parser.add_argument('--diagnostic-kallisto', required=True)
    parser.add_argument('--hdf-python', default=sys.executable)
    parser.add_argument('--manifest', type=Path, default=DATASETS)
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--samples', default='all')
    parser.add_argument('--suite', type=Path)
    parser.add_argument('--group')
    parser.add_argument('--index', type=Path)
    parser.add_argument('--reference-metadata', type=Path)
    parser.add_argument('--reads', type=Path)
    parser.add_argument('--reads2', type=Path)
    parser.add_argument('--modes', default='paired,single')
    parser.add_argument('--thread-counts', default='8')
    parser.add_argument('--fragment-length', type=float, default=200)
    parser.add_argument('--fragment-sd', type=float, default=20)
    parser.add_argument('--memory', action=argparse.BooleanOptionalAction, default=sys.platform in ('darwin', 'linux'))
    parser.add_argument('--stages', action='store_true')
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--resume-warmups', type=int, default=1)
    parser.add_argument('--pause-file', type=Path)
    args = parser.parse_args()
    try:
        args.thread_counts = [int(t) for t in args.thread_counts.split(',')]
        if (not args.thread_counts or len(set(args.thread_counts)) != len(args.thread_counts)
                or min(args.thread_counts) < 1 or args.runs < 1 or min(args.warmups, args.resume_warmups) < 0
                or not 0 < args.fragment_length < float('inf') or not 0 < args.fragment_sd < float('inf')):
            raise ValueError('Require positive unique thread counts, runs and fragment parameters; nonnegative warmups')
        if args.reads2 and not args.reads:
            raise ValueError('--reads2 requires --reads')
        tools = {key: executable(getattr(args, key)) for key in ('kallisto', 'kallistors')}
        args.diagnostic_kallisto = executable(args.diagnostic_kallisto)
        args.hdf_python = executable(args.hdf_python)
        if digest(args.diagnostic_kallisto) == digest(tools['kallisto']):
            raise ValueError('Use separate official timing and instrumented diagnostic kallisto executables')
        if not args.resume and args.out.exists() and any(args.out.iterdir()):
            raise ValueError('--out must be new or empty')
        args.out = args.out.resolve()
        configs = configurations(args, datasets(args.manifest))
        hdf_versions = json.loads(subprocess.check_output([str(args.hdf_python), '-c',
            'import h5py,numpy,json; print(json.dumps({"h5py":h5py.__version__,"numpy":numpy.__version__,"hdf5":h5py.version.hdf5_version}))'], text=True))
        files = sorted({p for c in configs for p in [c['index'], *c['reads']]})
        file_records = {p: fingerprint(p) for p in files}
        for case in configs:
            provenance = case.get('provenance')
            if not provenance:
                continue
            prep = provenance['preparation']
            outputs = {Path(m['path']).name: m for m in prep.get('outputs', prep.get('files', []))}
            for path in case['reads']:
                if Path(path).name not in outputs or outputs[Path(path).name]['sha256'] != file_records[path]['sha256']:
                    raise ValueError(f'Prepared FASTQ changed or missing provenance: {path}')
            case['expected_processed'] = prep.get('pairs_after', prep.get('selected_pairs'))
            if type(case['expected_processed']) is not int or case['expected_processed'] < 1:
                raise ValueError(f'Invalid prepared pair count: {case["id"]}')
        scripts = Path(__file__).parent
        configuration = {'cases': configs, 'files': list(file_records.values()),
                         'binaries': {key: fingerprint(p) for key, p in tools.items()},
                         'diagnostic': fingerprint(args.diagnostic_kallisto), 'hdf_python': {'invocation': str(args.hdf_python), 'binary': fingerprint(args.hdf_python), 'packages': hdf_versions},
                         'scripts': {p.name: digest(p) for p in scripts.glob('*.py')},
                         'warmups': args.warmups, 'runs': args.runs, 'memory': args.memory, 'stages': args.stages,
                         'environment': {k: v for k, v in timing_environment().items() if k.startswith('KALLISTORS_')},
                         'python': sys.version}
        args.out.mkdir(parents=True, exist_ok=True)
        protocol_path, runs_path = args.out / 'protocol.json', args.out / 'runs.json'
        rows, resume_events = [], []
        if args.resume:
            prior = json.loads(protocol_path.read_text())
            if prior['configuration'] != configuration:
                raise ValueError('Resume changes frozen inputs, binaries, scripts, environment or protocol; use a new --out')
            rows = json.loads(runs_path.read_text())
            resume_events = prior['resume_events'] + [{'at_utc': datetime.now(timezone.utc).isoformat(),
                                                      'prior_runs_sha256': digest(runs_path),
                                                      'warmups_per_tool_case': args.resume_warmups}]
        protocol = {'schema_version': 1, 'kind': 'quantification', 'configuration': configuration,
                    'resume_events': resume_events}
        write_json(protocol_path, protocol)
        write_json(runs_path, rows)
        completed = {(r['case'], r['tool'], r['phase'], r['replicate']) for r in rows}
        if len(completed) != len(rows) or any(not r['counters_match'] for r in rows):
            raise ValueError('Duplicate or failed saved runs')
        # Every requested case must pass before the first timing execution.
        preflight = {c['id']: ec_preflight(c, args, tools) for c in configs}
        def verify_tools():
            if any(fingerprint(binary) != configuration['binaries'][name] for name, binary in tools.items()):
                raise ValueError('A timing executable changed; use a new --out')
        verify_tools()
        if fingerprint(args.diagnostic_kallisto) != configuration['diagnostic']:
            raise ValueError('Diagnostic executable changed during EC preflight')
        expected = {}
        for row in rows:
            expected.setdefault(row['case'], row['counters'])

        def run(config, tool, phase, replicate):
            key = (config['id'], tool, phase, replicate)
            if key in completed:
                return
            verify_tools()
            folder = args.out / 'cases' / config['id']
            folder.mkdir(parents=True, exist_ok=True)
            prefix = folder / f'{phase}-{replicate}-{tool}'
            output = prefix.with_suffix('.output')
            archive_incomplete(prefix, output)
            command = quant_command(tools[tool], config, output)
            if args.stages and tool == 'kallistors':
                command.append('--timings')
            markers = (ORIGINAL_MARKERS if tool == 'kallisto' else RUST_MARKERS) if args.stages else None
            row = run_process(command, prefix, args.memory, markers)
            info = json.loads((output / 'run_info.json').read_text())
            counters = {k: info[k] for k in ('n_processed', 'n_pseudoaligned', 'n_unique')}
            if any(type(v) is not int or v < 0 for v in counters.values()):
                raise ValueError(f'Invalid quant counters: {config["id"]}')
            control = expected.setdefault(config['id'], counters)
            ec_total = preflight[config['id']]['comparison']['reference_total_counts']
            log_text = prefix.with_suffix('.log').read_text()
            rounds = re.search(r'algorithm ran for ([\d,]+) rounds', log_text)
            row['reported_em_rounds'] = int(rounds[1].replace(',', '')) if rounds else None
            row['versions'] = {key: info.get(key) for key in ('kallisto_version', 'index_version', 'k-mer length', 'n_targets')}
            row.update(case=config['id'], tool=tool, phase=phase, replicate=replicate, counters=counters,
                       counters_match=(counters == control and counters['n_pseudoaligned'] == ec_total
                                       and counters['n_processed'] == config.get('expected_processed', counters['n_processed'])
                                       and (not config.get('reference') or info['n_targets'] == config['reference']['total_targets'])))
            rows.append(row)
            write_json(runs_path, rows)
            if not row['counters_match']:
                raise RuntimeError(f'Quant counters disagree with the control or complete EC preflight: {key}')
            if phase != 'measured' or replicate != args.runs:
                shutil.rmtree(output)
            print(f'{config["id"]} {phase} {replicate} {tool}: {row["wall_seconds"]:.3f}s', flush=True)
            if args.pause_file and args.pause_file.exists():
                raise SystemExit(75)

        for config in configs:
            if args.resume:
                for rep in range(args.resume_warmups):
                    for tool in tools:
                        run(config, tool, f'resume-warmup-{len(resume_events)}', rep)
            for rep in range(args.warmups):
                for tool in tools:
                    run(config, tool, 'warmup', rep)
            for rep in range(1, args.runs + 1):
                for tool in list(tools) if rep % 2 else list(reversed(tools)):
                    run(config, tool, 'measured', rep)
        audits = {}
        for config in configs:
            folder = args.out / 'cases' / config['id']
            validation = {}
            for tool in tools:
                output = folder / f'measured-{args.runs}-{tool}.output'
                destination = folder / f'{tool}.hdf5.json'
                subprocess.run([str(args.hdf_python), str(scripts / 'check_hdf5_output.py'), str(output),
                                '--json', str(destination), '--significant-digits', '6' if tool == 'kallisto' else '17'], check=True)
                validation[tool] = json.loads(destination.read_text())
            subprocess.run([sys.executable, str(scripts / 'compare_abundance.py'), '--kallisto',
                            str(folder / f'measured-{args.runs}-kallisto.output/abundance.tsv'), '--kallistors',
                            str(folder / f'measured-{args.runs}-kallistors.output/abundance.tsv'),
                            '--summary-json', str(folder / 'abundance.json')], check=True)
            audits[config['id']] = {'hdf5': validation, 'abundance': json.loads((folder / 'abundance.json').read_text())}
        summary = {'schema_version': 1, 'kind': 'quantification', 'complete': True, 'cases': []}
        for config in configs:
            summary['cases'].append({'id': config['id'], **summarize([r for r in rows if r['case'] == config['id']]),
                                     'ec': preflight[config['id']], 'validation': audits[config['id']]})
        verify_tools()
        if {p.name: digest(p) for p in scripts.glob('*.py')} != configuration['scripts']:
            raise ValueError('Benchmark scripts changed during measurement')
        # Verify content again outside timing before marking the series complete.
        if any(fingerprint(p) != meta for p, meta in file_records.items()):
            raise ValueError('An input or index changed during measurement')
        write_json(args.out / 'summary.json', summary)
    except (ValueError, KeyError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Benchmark failed: {error}\n')


if __name__ == '__main__':
    main()
