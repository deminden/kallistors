#!/usr/bin/env python3
"""Collect fresh, serial modern-tool and original-kallisto publication panels."""
import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import traceback

from benchmark_common import fingerprint, write_json
from modern_benchmark.environment import check_analysis, load_environment, read_json, source_state


def reuse_gencode_inputs(source, destination, data):
    """Reuse prepared sequences and truth, never measurements or executables."""
    from modern_benchmark.inputs import bulk, single_cell, single_end
    from modern_benchmark.simulation import simulate

    previous = read_json(source / 'inputs.json')
    reference = previous['reference']
    if reference['kind'] != 'gencode50':
        raise ValueError('Expected the complete GENCODE v50 reference')
    for key in ('source', 'fasta', 'transcripts', 'gene_map'):
        if fingerprint(reference[key]['path']) != reference[key]:
            raise ValueError('Prepared reference changed: ' + key)
    cases = [bulk(data, accession) for accession in ('SRR13638690', 'SRR307903')]
    cases.append(simulate(reference, source / 'truth'))
    cases += [single_end(case) for case in cases.copy()]
    cell = single_cell(data)
    cell['strand'] = 'fw'
    cases.append(cell)
    inputs = {'reference': reference, 'cases': cases}
    if inputs != previous:
        raise ValueError('Prepared panel differs from its recorded inputs')
    write_json(destination / 'inputs.json', inputs)
    write_json(destination / 'input-reuse.json', {
        'source': fingerprint(source / 'inputs.json'),
        'scope': 'Verified reference, prepared FASTQs and independent simulation only; all timings are fresh.'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--kallistors', type=Path, default=Path('target/release/kallistors'))
    parser.add_argument('--kallisto', type=Path, default=Path('target/tools/kallisto-0.52.0'))
    parser.add_argument('--diagnostic-kallisto', type=Path,
                        default=Path('target/tools/kallisto-src/build/src/kallisto'))
    parser.add_argument('--tools-from', type=Path, help='Reuse verified pinned tool installations')
    parser.add_argument('--reuse-gencode-inputs', type=Path,
                        help='Reuse a verified prepared panel, including its independent simulation')
    args = parser.parse_args()
    check_analysis()
    root = args.workspace.resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / 'executor.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    scripts = Path(__file__).resolve().parent
    data = args.data.resolve()
    gencode, nascent, original = (root / name for name in ('gencode', 'nascent', 'original'))
    protocol = {
        'index_threads': [8], 'quantification_threads': [1, 8],
        'single_cell_threads': [1, 8], 'alevin_fry_complete_pipeline_threads': [8],
        'warmups': 1, 'measured_runs': 3, 'cellranger': 'deferred by user',
        'fresh_measurements': True, 'serial_execution': True,
        'kallistors': fingerprint(args.kallistors.resolve()),
        'kallisto': fingerprint(args.kallisto.resolve()),
        'diagnostic_kallisto': fingerprint(args.diagnostic_kallisto.resolve()),
        'source': source_state(), 'runner': fingerprint(Path(__file__)),
        'data': str(data),
        'tools_from': str(args.tools_from.resolve()) if args.tools_from else None,
        'reuse_gencode_inputs': str(args.reuse_gencode_inputs.resolve()) if args.reuse_gencode_inputs else None,
    }
    manifest = root / 'protocol.json'
    if manifest.exists() and read_json(manifest) != protocol:
        raise ValueError('Frozen full-comparison protocol changed; use a new workspace')
    write_json(manifest, protocol)

    def status(phase, **extra):
        write_json(root / 'status.json', {'phase': phase, **extra})
        print(phase, flush=True)

    def invoke(label, command):
        if source_state() != protocol['source']:
            raise ValueError('Implementation changed during the frozen comparison')
        status(label)
        with (root / 'commands.jsonl').open('a') as handle:
            handle.write(json.dumps(command) + '\n')
        subprocess.run(command, check=True)

    def modern(workspace, reference, action, *extra):
        invoke(workspace.name + ': ' + action, [sys.executable, str(scripts / 'modern_benchmarks.py'),
            '--workspace', str(workspace), '--data', str(data), '--reference', reference,
            action, *map(str, extra)])

    try:
        overrides = {'kallistors': str(args.kallistors.resolve())}
        if args.tools_from:
            environment = load_environment(args.tools_from.resolve())
            overrides.update({name: item['path'] for name, item in environment['tools'].items()
                              if name != 'kallistors' and 'path' in item})
        tool_args = [value for name, path in overrides.items() for value in ('--tool', name + '=' + path)]
        modern(gencode, 'gencode50', 'bootstrap', *tool_args)
        if args.reuse_gencode_inputs:
            status('verifying prepared GENCODE inputs and truth')
            reuse_gencode_inputs(args.reuse_gencode_inputs.resolve(), gencode, data)
        else:
            modern(gencode, 'gencode50', 'retrieve')
            modern(gencode, 'gencode50', 'prepare')
        for workspace, reference in ((gencode, 'gencode50'), (nascent, 'official-nac')):
            if workspace == nascent:
                environment = load_environment(gencode)
                tool_args = [value for name, item in environment['tools'].items() if 'path' in item
                             for value in ('--tool', name + '=' + item['path'])]
                modern(nascent, reference, 'bootstrap', *tool_args)
                modern(nascent, reference, 'retrieve')
                modern(nascent, reference, 'prepare')
            modern(workspace, reference, 'run', '--only', 'index')
            modern(workspace, reference, 'pilot')
            selected = ','.join(case['id'] for case in read_json(workspace / 'inputs.json')['cases'])
            modern(workspace, reference, 'run', '--only', selected)
            modern(workspace, reference, 'report')
        invoke('publishing modern comparisons', [sys.executable, str(scripts / 'publish_modern_benchmarks.py'),
            '--gencode', str(gencode), '--nascent', str(nascent)])
        invoke('original kallisto: full comparison', [sys.executable, str(scripts / 'run_kallisto_comparison.py'),
            '--workspace', str(original), '--gencode', str(gencode), '--nascent', str(nascent),
            '--data', str(data), '--kallisto', protocol['kallisto']['path'],
            '--diagnostic-kallisto', protocol['diagnostic_kallisto']['path']])
        if not read_json(original / 'summary.json')['complete']:
            raise ValueError('Original-kallisto comparison is incomplete')
        status('complete', modern_gencode=str(gencode), modern_nascent=str(nascent), original=str(original))
    except Exception:
        status('failed', error=traceback.format_exc())
        raise


if __name__ == '__main__':
    main()
