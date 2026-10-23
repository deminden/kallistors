"""Export publication tables only from complete measured triplicates."""
import csv
from itertools import combinations
from pathlib import Path
import statistics
import numpy as np
from benchmark_common import write_json
from .environment import read_json
from .metrics import load_bulk, load_sc, differences, sc_agreement


def table(out, name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (out / (name + '.tsv')).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)
    def display(v):
        if v is None:
            return 'NA'
        return f'{v:.4g}' if isinstance(v, float) else str(v)
    lines = ['| ' + ' | '.join(keys) + ' |', '| ' + ' | '.join('---' for _ in keys) + ' |']
    lines += ['| ' + ' | '.join(display(r.get(k)) for k in keys) + ' |' for r in rows]
    (out / (name + '.md')).write_text('\n'.join(lines) + '\n')
    def tex(v):
        return display(v).replace('\\', r'\textbackslash{}').replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')
    lines = [r'\begin{tabular}{' + 'l' * len(keys) + '}', r'\toprule',
             ' & '.join(tex(k) for k in keys) + r' \\', r'\midrule']
    lines += [' & '.join(tex(r.get(k)) for k in keys) + r' \\' for r in rows]
    lines += [r'\bottomrule', r'\end{tabular}']
    (out / (name + '.tex')).write_text('\n'.join(lines) + '\n')


def report(workspace):
    root, out = workspace / 'runs', workspace / 'tables'
    out.mkdir(parents=True, exist_ok=True)
    experiment = read_json(root / 'experiment.json')
    reference = experiment['inputs']['reference']
    rows = [read_json(p) for p in sorted(root.glob('**/result.json'))]
    measured = [r for r in rows if r['phase'] == 'measured']
    groups = {}
    for r in measured:
        groups.setdefault((r['case'], r['threads'], r['method']), []).append(r)
    performance, stages, cells, agreement, accuracy, variability, mapping = [], [], [], [], [], [], []
    status = experiment['design']['not_run'].copy()
    cases = {c['id']: c for c in experiment['inputs']['cases']}
    for (case, threads, method), group in sorted(groups.items()):
        if sorted(r['repetition'] for r in group) != [1, 2, 3]:
            status.append({'case': case, 'method': method, 'threads': threads, 'reason': 'Incomplete triplicate; not published in primary table'})
            continue
        warmups = [r for r in rows if r['case'] == case and r['threads'] == threads and r['method'] == method and r['phase'] == 'warmup']
        if len(warmups) != 1:
            raise ValueError('Every published triplicate requires one warmup')
        group.sort(key=lambda r: r['repetition'])
        median = lambda key: statistics.median(r[key] for r in group)
        p = {'case': case, 'method': method, 'threads': threads,
             'layout': 'index' if case == 'index' else cases[case].get('layout', 'paired'),
             'scope': group[0].get('scope', 'index' if case == 'index' else 'complete_pipeline'),
             'n': 3, 'wall_s': median('wall_seconds'),
             'wall_min_s': min(r['wall_seconds'] for r in group), 'wall_max_s': max(r['wall_seconds'] for r in group),
             'cpu_s': median('cpu_seconds'), 'rss_gib': median('peak_rss_bytes') / 2**30,
             'sampled_tree_rss_gib': median('sampled_tree_peak_rss_bytes') / 2**30}
        if case == 'index':
            p['index_gib'] = statistics.median(r['index_disk']['logical_bytes'] for r in group) / 2**30
        performance.append(p)
        if case != 'index':
            probes = [r.get('mapping_probe') for r in group]
            if all(probes):
                mapping.append({'case': case, 'method': method, 'threads': threads,
                                'scope': 'native mapping and intermediate output; standalone process',
                                'wall_s': statistics.median(v['wall_seconds'] for v in probes),
                                'cpu_s': statistics.median(v['cpu_seconds'] for v in probes),
                                'rss_gib': statistics.median(v['peak_rss_bytes'] for v in probes) / 2**30})
            elif group[0].get('mapping_wall_seconds') is not None:
                seq = [next(s for s in r['stages'] if s['stage'] == 'mapping') for r in group]
                mapping.append({'case': case, 'method': method, 'threads': threads,
                                'scope': 'native mapping and intermediate output; first pipeline stage',
                                'wall_s': statistics.median(v['wall_seconds'] for v in seq),
                                'cpu_s': statistics.median(v['process_user_seconds'] + v['process_system_seconds'] for v in seq),
                                'rss_gib': statistics.median(v['peak_rss_bytes'] for v in seq) / 2**30})
        for i, stage in enumerate(group[0]['stages']):
            seq = [r['stages'][i] for r in group]
            stages.append({'case': case, 'method': method, 'threads': threads, 'stage': stage['stage'],
                           'stage_workers': stage['requested_workers'],
                           'wall_s': statistics.median(s['wall_seconds'] for s in seq),
                           'cpu_s': statistics.median(s['process_user_seconds'] + s['process_system_seconds'] for s in seq),
                           'rss_gib': statistics.median(s['peak_rss_bytes'] for s in seq) / 2**30})
        if case == 'index' or group[0].get('scope') == 'mapping_only':
            continue
        if cases[case]['kind'] == 'single_cell':
            cells.append({'case': case, 'method': method, 'threads': threads, **group[0]['matrix_summary']})
            # Expose run-to-run count variability separately from cross-tool disagreement.
            for r in group[1:]:
                variability.append({'case': case, 'method': method, 'threads': threads, 'rep_a': 1,
                    'rep_b': r['repetition'], 'metrics': sc_agreement(load_sc(group[0]['output'], method), load_sc(r['output'], method))})
        else:
            a = load_bulk(group[0]['output'], method, reference)
            for r in group[1:]:
                b = load_bulk(r['output'], method, reference)
                variability.append({'case': case, 'method': method, 'threads': threads, 'rep_a': 1,
                                    'rep_b': r['repetition'], 'metrics': {k: differences(a[k], b[k]) for k in a}})
            if cases[case].get('truth'):
                truth = np.load(cases[case]['truth'])
                for metric, truth_key in [('transcript_counts', 'fragment_counts'), ('transcript_tpm', 'molecule_tpm')]:
                    accuracy.append({'case': case, 'method': method, 'threads': threads, 'metric': metric,
                                     **differences(a[metric], truth[truth_key])})
    # Use measured replicate 1, selected before inspecting timing or agreement.
    for case in cases.values():
        for threads in (1, 8):
            keys = [k for k, g in groups.items() if k[0] == case['id'] and k[1] == threads and
                    len(g) == 3 and not g[0].get('scope') == 'mapping_only']
            for ka, kb in combinations(sorted(keys), 2):
                a = next(r for r in groups[ka] if r['repetition'] == 1)
                b = next(r for r in groups[kb] if r['repetition'] == 1)
                base = {'case': case['id'], 'threads': threads, 'method_a': ka[2], 'method_b': kb[2]}
                if case['kind'] == 'bulk':
                    x, y = load_bulk(a['output'], ka[2], reference), load_bulk(b['output'], kb[2], reference)
                    agreement.extend({**base, 'metric': k, **differences(x[k], y[k])} for k in x)
                else:
                    item = sc_agreement(load_sc(a['output'], ka[2]), load_sc(b['output'], kb[2]))
                    agreement.append({**base, 'metric': 'gene_umi', **{k: v for k, v in item.items() if not isinstance(v, dict)}})
                    write_json(out / f"{case['id']}-t{threads}-matrix-agreement.json", item)
    expected = [('index', t, m) for t in experiment['design']['index_threads'] for m in experiment['design']['builders']]
    for c in cases.values():
        methods = experiment['design']['bulk_methods' if c['kind'] == 'bulk' else 'single_cell_methods']
        workers = experiment['design'].get('single_cell_threads', [1, 8]) if c['kind'] == 'single_cell' else experiment['design']['threads']
        expected += [(c['id'], t, m) for t in workers for m in methods]
    excluded = {(r.get('case'), r.get('threads'), r.get('method')) for r in status}
    for key in expected:
        if key not in groups and key not in excluded:
            status.append({'case': key[0], 'threads': key[1], 'method': key[2],
                           'reason': 'Not measured; no performance or accuracy conclusion available.'})
    table(out, 'performance', performance)
    table(out, 'stages', stages)
    table(out, 'mapping', mapping)
    table(out, 'single-cell-counts', cells)
    table(out, 'agreement', agreement)
    table(out, 'accuracy', accuracy)
    table(out, 'not-run', status)
    write_json(out / 'repeat-variability.json', variability)
    from .publication import export
    export(workspace)
    from .experiment import workflow_state
    write_json(out / 'manifest.json', {'analysis_workflow': workflow_state(), 'experiment': experiment, 'rows': len(rows), 'measured': len(measured),
                                      'interpretation': 'Within-host medians of three after one warmup. Agreement is not accuracy. No Cell Ranger results; no one-thread alevin-fry pipeline.'})
    print(f'Tables written to {out}', flush=True)
