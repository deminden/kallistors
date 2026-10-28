"""Freeze provenance, reuse native runners, and execute the two-tool panel serially."""
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

from benchmark_common import fingerprint, write_json
from modern_benchmark.environment import check_analysis, content_id, load_environment, read_json
from modern_benchmark.inputs import prefix_case
from modern_benchmark.validation import validate_completed_panel

SCRIPTS = Path(__file__).resolve().parents[1]


def freeze_tools(args):
    current = load_environment(args.gencode)
    other = load_environment(args.nascent)
    if current['source_sha256'] != other['source_sha256']:
        raise ValueError('Both preceding panels must use the same kallistors source snapshot')
    for name in ('kallistors', 'bustools'):
        if current['tools'][name]['sha256'] != other['tools'][name]['sha256']:
            raise ValueError('Preceding panels used different executable bytes: ' + name)
    record = args.workspace / 'frozen-tools.json'
    if record.exists():
        tools = read_json(record)
        for item in tools.values():
            if fingerprint(item['path']) != item:
                raise ValueError('Frozen original-comparison executable changed')
    else:
        paths = {'kallisto': args.kallisto, 'diagnostic-kallisto': args.diagnostic_kallisto,
                 **{name: Path(current['tools'][name]['path']) for name in ('kallistors', 'bustools')}}
        tools = {}
        for name, path in paths.items():
            dest = args.workspace / 'bin' / name
            dest.parent.mkdir(exist_ok=True)
            shutil.copy2(path, dest)
            tools[name] = fingerprint(dest)
        write_json(record, tools)
    for name in ('kallistors', 'bustools'):
        if tools[name]['sha256'] != current['tools'][name]['sha256']:
            raise ValueError('Refresh must reuse the preceding panels\' executable: ' + name)
    if tools['kallisto']['sha256'] == tools['diagnostic-kallisto']['sha256']:
        raise ValueError('Official performance and instrumented correctness controls must be separate')
    for name in ('kallisto', 'diagnostic-kallisto'):
        version = subprocess.check_output([tools[name]['path'], 'version'], text=True).strip()
        if version != 'kallisto, version 0.52.0':
            raise ValueError('Expected the pinned kallisto 0.52.0 control: ' + version)
    return tools, current


def airway(data):
    model = read_json(SCRIPTS / 'benchmarks/kallisto-extra-library.json')
    directory = data / 'srr1039508/trimmed-fastp'
    prep = read_json(directory / 'preparation.json')
    if prep['trimmer']['name'] != 'fastp' or prep['trimmer']['version'] != '1.4.0':
        raise ValueError('Airway library requires the same pinned fastp preparation')
    if [(x['bytes'], x['sha256']) for x in prep['outputs']] != [(x['bytes'], x['sha256']) for x in model['prepared_inputs']]:
        raise ValueError('Airway calibration inputs differ')
    hist = model['calibration']['histogram_counts']
    n = sum(hist)
    mean = math.fsum(i * count for i, count in enumerate(hist)) / n
    sd = math.sqrt(math.fsum(count * (i - mean)**2 for i, count in enumerate(hist)) / n)
    return {'id': 'SRR1039508', 'kind': 'bulk', 'pairs': prep['pairs_after'], 'strand': model['strand'],
        'read1': [str((directory / Path(prep['outputs'][0]['path']).name).resolve())],
        'read2': [str((directory / Path(prep['outputs'][1]['path']).name).resolve())],
        'single_end_fragment_model': {'mean': mean, 'sd': sd, 'samples': n, 'source': 'frozen matching paired-end calibration'},
        'provenance': prep, 'library_metadata': {k: model[k] for k in ('strand', 'strand_basis', 'sources')}}


def suite_case(case, index, reference, name=None, prefix=False):
    model = {'mean': 250, 'sd': 30} if case.get('truth') else case['single_end_fragment_model']
    result = {'name': name or case['id'], 'index': str(index), 'reference_metadata': str(reference),
        'modes': ['paired', 'single'], 'strand': case.get('strand', 'unstranded'),
        'fragment_length': model['mean'], 'fragment_sd': model['sd'], 'expected_processed': case['pairs']}
    if not prefix and case['id'].startswith('SRR'):
        result['sample'] = case['id']
    else:
        result.update(reads=case['read1'][0], reads2=case['read2'][0])
    return result


def execute(args, status):
    check_analysis()
    args.gencode, args.nascent, args.data = (p.resolve() for p in (args.gencode, args.nascent, args.data))
    for root in (args.gencode, args.nascent):
        _, previous = validate_completed_panel(root)
        if read_json(root / 'inputs.json') != previous['inputs']:
            raise ValueError('Inputs changed after the preceding comparison')
    tools, environment = freeze_tools(args)
    paths = {name: record['path'] for name, record in tools.items()}
    panels = {}
    for name, root in [('gencode', args.gencode), ('nascent', args.nascent)]:
        inputs = read_json(root / 'inputs.json')
        cases = [c for c in inputs['cases'] if c['kind'] == 'bulk' and c.get('layout', 'paired') != 'single']
        if name == 'gencode':
            cases.insert(1, airway(args.data))
        panels[name] = {'reference': inputs['reference'], 'cases': cases}
    sc = next(c for c in read_json(args.gencode / 'inputs.json')['cases'] if c['kind'] == 'single_cell')
    sources = [*SCRIPTS.glob('*.py'), *SCRIPTS.glob('modern_benchmark/*.py'),
               *Path(__file__).parent.glob('*.py'), SCRIPTS / 'benchmarks/kallisto-extra-library.json',
               *[p for p in SCRIPTS.glob('benchmarks/modern/*') if p.is_file()],
               *SCRIPTS.glob('diagnostics/kallisto_v0.52.0*.patch')]
    workflow = {str(p.relative_to(SCRIPTS)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    spec = {'tools': tools, 'analysis': check_analysis(), 'kallistors_source_sha256': environment['source_sha256'],
        'panels': panels, 'single_cell': sc, 'workflow': workflow,
        'design': {'read_workers': [1, 8], 'single_cell_workers': [1, 8], 'index_workers': [8], 'warmups': 1, 'measured_runs': 3,
                   'bulk_layouts': ['paired', 'single'], 'scope': 'Fresh sequential two-tool controls, after alternatives',
                   'single_cell': 'Full mapping and prepared-FASTQ-to-matrix with identical bustools settings',
                   'bulk_stages': 'CLI-message wall-time boundaries; not isolated mapping kernels or per-stage CPU/RSS',
                   'correctness': 'Full EC counts; complete weighted barcode/UMI/transcript-set/flag/count and canonical matrices'}}
    manifest = args.workspace / 'experiment.json'
    if manifest.exists() and read_json(manifest) != spec:
        raise ValueError('Original-comparison protocol changed; use a new workspace')
    write_json(manifest, spec)
    revision = args.workspace / 'workflow' / content_id(workflow)
    for relative in workflow:
        dest = revision / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(SCRIPTS / relative, dest)
    commands_file = args.workspace / 'commands.json'
    commands = read_json(commands_file) if commands_file.exists() else []
    def invoke(label, command, summary):
        if summary.exists() and read_json(summary).get('complete'):
            return
        command = list(map(str, command))
        status(label)
        commands.append(command)
        write_json(commands_file, commands)
        subprocess.run(command, check=True)
    for name, panel in panels.items():
        ref = panel['reference']
        status('verifying ' + name + ' inputs')
        for item in [ref['fasta'], ref['transcripts'], ref['gene_map']]:
            if fingerprint(item['path']) != item:
                raise ValueError('Reference changed before original comparison')
        index_out = args.workspace / name / 'index'
        invoke('index construction: ' + name, [sys.executable, SCRIPTS / 'benchmark_index.py',
            '--kallisto', paths['kallisto'], '--kallistors', paths['kallistors'], '--fasta', ref['fasta']['path'],
            '--threads', '8', '--warmups', '1', '--runs', '3', '--memory', '--keep-last-index', '--out', index_out],
            index_out / 'summary.json')
        indexes = {tool: index_out / f'measured-3-{tool}.idx' for tool in ('kallisto', 'kallistors')}
        metadata = args.workspace / name / 'reference.json'
        write_json(metadata, {**ref, 'total_targets': ref['targets']})
        for compatibility in (True, False):
            selected = panel['cases']
            if compatibility:
                selected = [prefix_case(c, args.workspace / 'prefixes' / c['id']) for c in selected]
            filename = args.workspace / name / ('compatibility-suite.json' if compatibility else 'bulk-suite.json')
            write_json(filename, {'cases': [suite_case(c, indexes['kallistors' if compatibility else 'kallisto'], metadata,
                name=c['id'], prefix=compatibility) for c in selected]})
            out = args.workspace / name / ('compatibility' if compatibility else 'bulk')
            command = [sys.executable, SCRIPTS / 'benchmark_reads.py', '--kallisto', paths['kallisto'],
                '--kallistors', paths['kallistors'], '--diagnostic-kallisto', paths['diagnostic-kallisto'],
                '--hdf-python', sys.executable, '--data', args.data, '--suite', filename,
                '--thread-counts', '8' if compatibility else '1,8', '--warmups', '0' if compatibility else '1',
                '--runs', '1' if compatibility else '3', '--memory', '--out', out]
            if (out / 'protocol.json').exists() and not (out / 'summary.json').exists():
                command += ['--resume']
            if not compatibility:
                command += ['--stages']
            invoke(('index compatibility: ' if compatibility else 'bulk quantification: ') + name, command, out / 'summary.json')
    from .single_cell import run_single_cell
    run_single_cell(args.workspace, paths, sc, panels['gencode']['reference'],
                    args.workspace / 'gencode/index/measured-3-kallisto.idx', spec, status)
    from .report import report
    status('analysing original comparison')
    if {relative: hashlib.sha256((SCRIPTS / relative).read_bytes()).hexdigest() for relative in workflow} != workflow:
        raise ValueError('Comparison workflow changed during execution')
    for tool, item in tools.items():
        if fingerprint(item['path']) != item:
            raise ValueError('Comparison executable changed: ' + tool)
    report(args.workspace, spec)
