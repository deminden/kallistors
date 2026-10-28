"""Verify completed comparison panels before a dependent benchmark starts."""

from collections import Counter

from benchmark_common import fingerprint
from .environment import read_json


def validate_completed_panel(workspace):
    summary = read_json(workspace / 'tables/publication-summary.json')
    spec = read_json(workspace / 'runs/experiment.json')
    design = spec['design']
    if design.get('warmups', 1) != 1 or design.get('measured_runs', 3) != 3:
        raise ValueError('Publication requires one warmup and three measurements per condition')
    conditions = [('index', workers, method) for workers in design['index_threads']
                  for method in design['builders']]
    for case in spec['inputs']['cases']:
        single_cell = case['kind'] == 'single_cell'
        workers = design.get('single_cell_threads', [1, 8]) if single_cell else design['threads']
        methods = design['single_cell_methods' if single_cell else 'bulk_methods']
        conditions.extend((case['id'], count, method) for count in workers for method in methods)
    expected = len(conditions)
    rows = [read_json(p) for p in (workspace / 'runs').glob('**/result.json')]
    if summary['complete_groups'] != expected or len(rows) != expected * 4:
        raise ValueError(f'Publication requires a complete warmup-plus-triplicate panel: {workspace}')
    # Totals alone cannot distinguish a complete panel from duplicate or mislabeled records.
    expected_rows = Counter({(*condition, phase, repetition): 1 for condition in conditions
                             for phase, repetition in [('warmup', 0), ('measured', 1),
                                                       ('measured', 2), ('measured', 3)]})
    observed_rows = Counter((row['case'], row['threads'], row['method'], row['phase'], row['repetition'])
                            for row in rows)
    if observed_rows != expected_rows:
        raise ValueError('Missing, duplicate or undeclared observations in the experiment panel')
    executables = {item['path']: item for item in spec['tools'].values()}
    for item in executables.values():
        if fingerprint(item['path']) != item:
            raise ValueError('Publication executable differs from the frozen experiment')
    for row in rows:
        stages = row['stages'] + row.get('mapping_probe', {}).get('stages', [])
        if any(s['command'][0] not in executables for s in stages):
            raise ValueError('A result used an executable outside the frozen experiment')
        if any(s.get('timing_interrupted') for s in stages):
            raise ValueError('Interrupted timing is not a valid publication observation')
    return summary, spec
