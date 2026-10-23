"""Serial, balanced warmup-plus-three experiments with immutable resume records."""
import shutil
from pathlib import Path
from benchmark_common import digest, fingerprint, write_json, system_context
from .environment import content_id, read_json
from .adapters import BUILDERS, BULK, SC, index_plan, bulk_plan, sc_plan, bulk_mapping_plan
from .measure import pipeline, disk_usage


def plan(inputs):
    return {'reference_targets': inputs['reference']['targets'], 'threads': [1, 8], 'index_threads': [8],
            'single_cell_threads': [1, 8],
            'warmups': 1, 'measured_runs': 3, 'builders': list(BUILDERS),
            'cases': [{'id': c['id'], 'kind': c['kind'], 'layout': c.get('layout', 'paired'),
                       'input_records': c.get('pairs', c.get('reads'))} for c in inputs['cases']],
            'bulk_methods': list(BULK), 'single_cell_methods': list(SC),
            'not_run': [{'method': 'piscem-alevin-fry', 'threads': 1,
                         'scope': 'complete_pipeline', 'reason': 'One-worker mapping only; alevin-fry enforces at least two workers.'},
                        {'method': 'cellranger', 'threads': 8, 'reason': 'Deferred by user; requires x86-64 Linux.'}]}


def workflow_state():
    root = Path(__file__).resolve().parents[1]
    paths = [*Path(__file__).parent.glob('*.py'), root / 'modern_benchmarks.py', root / 'benchmark_common.py', root / 'prepare_benchmark.py', root / 'prepare_reference.py', root / 'publish_modern_benchmarks.py']
    paths += list((root / 'benchmarks/modern').glob('*'))
    return {str(p.relative_to(root)): digest(p) for p in sorted(paths) if p.is_file()}


def run_one(root, key, case, method, threads, phase, repetition, commands, spec):
    destination = root / key
    record = destination / 'result.json'
    identity = content_id({'commands': commands, 'spec': spec})
    if record.exists():
        saved = read_json(record)
        if saved['identity'] != identity:
            # An explicit amendment may retain an unchanged competitor control
            # while replacing another executable. Never relabel its provenance.
            amendment = root / 'amendments/retained-trials.json'
            approved = read_json(amendment).get(key) if amendment.exists() else None
            actual = [(s['stage'], s['requested_workers'], s['command']) for s in saved['stages']]
            if not (approved and approved['previous_identity'] == saved['identity']
                    and approved['new_identity'] == identity
                    and approved['record_sha256'] == fingerprint(record)['sha256']
                    and actual == commands):
                raise ValueError(f'Changed specification for completed run: {key}')
        return saved
    if destination.exists():
        archive = root / 'incomplete' / key
        archive.parent.mkdir(parents=True, exist_ok=True)
        if archive.exists():
            raise ValueError(f'Inspect previous incomplete run: {archive}')
        shutil.move(destination, archive)
    destination.mkdir(parents=True)
    if case == 'index':
        (destination / 'artifacts').mkdir(parents=True, exist_ok=True)
        if method == 'piscem':
            (destination / 'artifacts/piscem').mkdir(parents=True, exist_ok=True)
    result = pipeline(commands, destination)
    result.update(identity=identity, case=case, method=method, threads=threads,
                  phase=phase, repetition=repetition, output=str(destination))
    write_json(record, result)
    print(f"{phase} {case} {method} t{threads} #{repetition}: "
          f"{result['wall_seconds']:.3f} s; {result['peak_rss_bytes'] / 2**30:.3f} GiB", flush=True)
    return result



def reuse_primary_indexes(workspace, root, spec):
    """Validate frozen primary indexes before reusing them for untimed pilots."""
    primary = workspace / 'runs'
    manifest = primary / 'experiment.json'
    if not manifest.exists():
        return False
    previous = read_json(manifest)
    if previous['tools'] != spec['tools'] or previous['inputs']['reference'] != spec['inputs']['reference']:
        raise ValueError('Pilot index reuse requires identical frozen tools and reference')
    pointers = [primary / f'index-current-{builder}.json' for builder in BUILDERS]
    if not all(p.exists() for p in pointers):
        return False
    evidence = {'experiment': fingerprint(manifest), 'indexes': {}}
    for builder, pointer in zip(BUILDERS, pointers, strict=True):
        saved = read_json(pointer)
        record = read_json(primary / saved['run'] / 'result.json')
        if not record.get('index_files'):
            raise ValueError(f'Missing index verification records for {builder}')
        for item in record['index_files']:
            if fingerprint(item['path']) != item:
                raise ValueError(f'Index changed before pilot: {item["path"]}')
        evidence['indexes'][builder] = {'pointer': saved, 'files': record['index_files']}
        write_json(root / pointer.name, saved)
    write_json(root / 'index-reuse.json', evidence)
    return True


def execute(workspace, tools, inputs, pilot=False, only=None, amend=False):
    root = workspace / ('pilot' if pilot else 'runs')
    root.mkdir(parents=True, exist_ok=True)
    selection = set(only.split(',')) if only else None
    if selection and selection - {'index', *[c['id'] for c in inputs['cases']]}:
        raise ValueError('Unknown --only case: ' + str(selection))
    spec = {'design': plan(inputs), 'inputs': inputs, 'tools': {k: fingerprint(v) for k, v in tools.items()},
            'workflow': workflow_state(), 'pilot': pilot}
    manifest = root / 'experiment.json'
    if manifest.exists() and read_json(manifest) != spec:
        old = read_json(manifest)
        # An unstarted biological case may be amended after its pilot. Completed
        # cases, measurement code, tools and reference are immutable. Keep the audit.
        completed = {read_json(p)['case'] for p in root.glob('**/result.json')}
        if any(old[k] != spec[k] for k in ('tools', 'pilot')) or old['inputs']['reference'] != inputs['reference']:
            raise ValueError(f'Frozen tools or reference changed: {manifest}')
        if old['workflow'] != spec['workflow']:
            timing_files = ('benchmark_common.py', 'modern_benchmark/measure.py')
            if not (amend and completed <= {'index'} and selection and 'index' not in selection
                    and all(old['workflow'][p] == spec['workflow'][p] for p in timing_files)):
                raise ValueError('Workflow changed; only an audited revision before read measurements is permitted')
            # Read-only extensions may add SE adapters after the index panel.
            # Prove every completed build still has identical commands/workers.
            evidence = []
            for record in root.glob('index/**/result.json'):
                row = read_json(record)
                expected = index_plan(tools, row['method'], inputs['reference'],
                                      Path(row['output']) / 'artifacts', row['threads'])
                actual = [(s['stage'], s['requested_workers'], s['command']) for s in row['stages']]
                if actual != expected:
                    raise ValueError('Completed index commands changed; use a new workspace')
                evidence.append({'run': str(record), 'commands_unchanged': True})
            write_json(root / 'amendments' / (content_id(old) + '-index-commands.json'), evidence)
        old_cases = {c['id']: c for c in old['inputs']['cases']}
        new_cases = {c['id']: c for c in inputs['cases']}
        if any(old_cases.get(c) != new_cases.get(c) for c in completed if c != 'index'):
            raise ValueError('A completed case changed; use a new workspace')
        write_json(root / 'amendments' / (content_id(old) + '.json'), old)
    write_json(manifest, spec)
    for relative in spec['workflow']:
        source = Path(__file__).resolve().parents[1] / relative
        dest = root / 'workflow-revisions' / content_id(spec['workflow']) / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(source, dest)
    if not (root / 'host.json').exists():
        write_json(root / 'host.json', system_context())
    repeats = [('pilot', 0)] if pilot else [('warmup', 0), *[('measured', i) for i in range(1, 4)]]
    # A separate directory per builder and trial prevents overwrites. Each subsequent
    # trial replaces the previous successful index after its size/hash evidence is saved.
    reused = pilot and reuse_primary_indexes(workspace, root, spec)
    if not reused and (selection is None or 'index' in selection):
        for threads in (8,):
            for phase, rep in repeats:
                order = BUILDERS[rep % len(BUILDERS):] + BUILDERS[:rep % len(BUILDERS)]
                for method in order:
                    key = f'index/t{threads}/{phase}-{rep}/{method}'
                    dest = root / key
                    commands = index_plan(tools, method, inputs['reference'], dest / 'artifacts', threads)
                    index_spec = {k: spec[k] for k in ('tools', 'workflow', 'pilot')}
                    index_spec['reference'] = inputs['reference']
                    row = run_one(root, key, 'index', method, threads, phase, rep, commands, index_spec)
                    artifact = dest / 'artifacts'
                    if not artifact.exists():
                        continue  # A later successful build superseded these recorded artifacts.
                    if not row.get('index_files'):
                        row['index_disk'] = disk_usage(artifact)
                        row['index_files'] = [fingerprint(p) for p in sorted(artifact.rglob('*')) if p.is_file()]
                        write_json(dest / 'result.json', row)
                    # Keep only the newest successful index of each builder.
                    pointer = root / f'index-current-{method}.json'
                    old = read_json(pointer) if pointer.exists() else None
                    if old and old['path'] != str(artifact) and Path(old['path']).exists():
                        shutil.rmtree(old['path'])
                    write_json(pointer, {'path': str(artifact), 'run': key})
    indexes = root / 'indexes'
    indexes.mkdir(exist_ok=True)
    for builder in BUILDERS:
        pointer = root / f'index-current-{builder}.json'
        if not pointer.exists():
            raise ValueError('Build indexes first, or include index in --only')
        src = Path(read_json(pointer)['path'])
        for child in src.iterdir():
            link = indexes / child.name
            if link.is_symlink():
                link.unlink()
            elif link.exists():
                raise ValueError(f'Unexpected non-symlink: {link}')
            link.symlink_to(child, target_is_directory=child.is_dir())
    for case in inputs['cases']:
        if selection is not None and case['id'] not in selection:
            continue
        methods = BULK if case['kind'] == 'bulk' else SC
        workers = spec['design']['single_cell_threads'] if case['kind'] == 'single_cell' else spec['design']['threads']
        for threads in ([8] if pilot else workers):
            for phase, rep in repeats:
                order = methods[rep % len(methods):] + methods[:rep % len(methods)]
                for method in order:
                    key = f"{case['id']}/t{threads}/{phase}-{rep}/{method}"
                    dest = root / key
                    adapter = bulk_plan if case['kind'] == 'bulk' else sc_plan
                    commands = adapter(tools, method, case, inputs['reference'], indexes, dest, threads)
                    case_spec = {k: spec[k] for k in ('tools', 'workflow', 'pilot')}
                    case_spec.update(reference=inputs['reference'], case=case)
                    row = run_one(root, key, case['id'], method, threads, phase, rep, commands, case_spec)
                    if case['kind'] == 'bulk' and 'mapping_probe' not in row:
                        mapping_out = dest / 'mapping-probe'
                        mapping_commands = bulk_mapping_plan(tools, method, case, inputs['reference'], indexes, mapping_out, threads)
                        if mapping_commands:
                            if mapping_out.exists():
                                archived = root / 'incomplete' / key / 'mapping-probe'
                                archived.parent.mkdir(parents=True, exist_ok=True)
                                shutil.move(mapping_out, archived)
                            row['mapping_probe'] = pipeline(mapping_commands, mapping_out)
                            write_json(dest / 'result.json', row)
                    # Matrix/abundance validation is deliberately outside the timer.
                    if case['kind'] == 'bulk':
                        from .metrics import load_bulk
                        values = load_bulk(dest, method, inputs['reference'])
                        row['mapped_fragments_estimated'] = float(values['transcript_counts'].sum())
                    elif not (method == 'piscem-alevin-fry' and threads == 1):
                        from .metrics import load_sc, sc_summary
                        row['matrix_summary'] = sc_summary(load_sc(dest, method))
                    else:
                        row['scope'] = 'mapping_only'
                    if 'retired_intermediates' not in row:
                        from .artifacts import retire_intermediates
                        row['retired_intermediates'] = retire_intermediates(dest)
                    write_json(dest / 'result.json', row)
