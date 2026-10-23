"""Publication views and prespecified thread/truth audits; no timed work."""
from functools import lru_cache
from pathlib import Path
import csv
import statistics
import shutil
import numpy as np
from benchmark_common import write_json
from .environment import read_json, content_id
from .metrics import differences, load_bulk, load_sc, sc_agreement
from .tables import table

LABELS = {
    'kallistors': 'kallistors', 'salmon-sketch': 'Salmon sketch',
    'salmon-sa': 'Salmon selective alignment', 'piscem-infer': 'piscem + piscem-infer',
    'kallistors-bustools': 'kallistors + bustools',
    'piscem-alevin-fry': 'piscem + alevin-fry',
    'salmon': 'Salmon', 'piscem': 'piscem',
}


def audit_inputs(workspace):
    """Validate native record counts, outside all performance measurements."""
    audited = []
    for scope in ('runs', 'pilot'):
        root = workspace / scope
        if not (root / 'experiment.json').exists():
            continue
        spec = read_json(root / 'experiment.json')
        cases = {c['id']: c for c in spec['inputs']['cases']}
        for path in sorted(root.glob('**/result.json')):
            row = read_json(path)
            if row['case'] == 'index':
                if row['threads'] != 8:
                    raise ValueError('Index trials must use 8 workers')
                continue
            case, dest, method = cases[row['case']], Path(row['output']), row['method']
            expected = case.get('pairs', case.get('reads'))
            if method in ('kallistors', 'kallistors-bustools'):
                source = dest / ('quant/run_info.json' if method == 'kallistors' else 'mapping/run_info.json')
                key = 'n_processed'
            elif method.startswith('salmon'):
                source, key = dest / 'quant/aux_info/meta_info.json', 'num_processed'
            else:
                source = dest / ('mapping.map_info.json' if method == 'piscem-infer' else 'mapping/map_info.json')
                key = 'num_reads'
            processed = read_json(source)[key]
            if processed != expected:
                raise ValueError(f'Processed-record mismatch: {path}: {processed} != {expected}')
            audited.append({'scope': scope, 'case': row['case'], 'method': method,
                'threads': row['threads'], 'phase': row['phase'], 'repetition': row['repetition'],
                'processed': processed, 'expected': expected, 'exact': True})
    write_json(workspace / 'tables/processed-input-audit.json', audited)
    return len(audited)


def _range(values, scale=1, digits=3):
    lo, hi = min(values) / scale, max(values) / scale
    return f'{lo:.{digits}f}–{hi:.{digits}f}'


def _paper_table(out, name, rows, caption, index=False):
    """Booktabs/longtable output with grouped datasets and readable numeric columns."""
    if not rows:
        return
    table(out, name, rows)
    def escape(value):
        return str(value).replace('\\', r'\textbackslash{}').replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('–', '--')
    headers = ['Method', 'Wall, s', 'CPU, s', 'Peak RSS, GiB', 'Index, GiB'] if index else [
        'Method', 'Workers', 'Wall, s', 'CPU, s', 'Peak RSS, GiB']
    lines = [r'% Requires \usepackage{booktabs,longtable}.', r'\begingroup', r'\small',
        r'\setlength{\tabcolsep}{4pt}', r'\begin{longtable}{@{}lrrrr@{}}',
        r'\caption{' + escape(caption) + r'}\label{tab:' + name + r'}\\', r'\toprule',
        ' & '.join(headers) + r' \\', r'\midrule', r'\endfirsthead', r'\toprule',
        ' & '.join(headers) + r' \\', r'\midrule', r'\endhead',
        r'\bottomrule', r'\endfoot']
    previous = None
    for row in rows:
        if not index and row['Dataset'] != previous:
            previous = row['Dataset']
            lines += [r'\addlinespace', r'\multicolumn{5}{@{}l}{\textit{' + escape(previous) + r'}} \\']
        wall = row['Wall, s'] + ' [' + row['Wall range, s'] + ']'
        rss = row['RSS, GiB'] + ' [' + row['RSS range, GiB'] + ']'
        values = [row['Method'], wall, row['CPU, s'], rss, row['Index, GiB']] if index else [
            row['Method'], row['Workers'], wall, row['CPU, s'], rss]
        lines.append(' & '.join(escape(v) for v in values) + r' \\')
    lines += [r'\end{longtable}', r'\noindent\footnotesize Medians of three measured runs after one warmup; brackets show observed minima and maxima. CPU is user plus system time. Pipeline peak RSS is the maximum stage peak.', r'\endgroup']
    (out / (name + '.tex')).write_text('\n'.join(lines) + '\n')


def export(workspace):
    """Derive tables from existing native results, retaining analysis provenance."""
    root, out = workspace / 'runs', workspace / 'tables'
    experiment = read_json(root / 'experiment.json')
    from .experiment import workflow_state
    analysis_workflow = workflow_state()
    source_root = Path(__file__).resolve().parents[1]
    revision = workspace / 'analysis-workflow-revisions' / content_id(analysis_workflow)
    for relative in analysis_workflow:
        dest = revision / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(source_root / relative, dest)
    audited_inputs = audit_inputs(workspace)
    reference = experiment['inputs']['reference']
    cases = {c['id']: c for c in experiment['inputs']['cases']}
    rows = [read_json(p) for p in sorted(root.glob('**/result.json'))]
    groups = {}
    for row in rows:
        if row['phase'] == 'measured':
            groups.setdefault((row['case'], row['threads'], row['method']), []).append(row)
    complete = {}
    for key, group in groups.items():
        if sorted(r['repetition'] for r in group) != [1, 2, 3]:
            continue
        warm = [r for r in rows if (r['case'], r['threads'], r['method']) == key and r['phase'] == 'warmup']
        if len(warm) != 1:
            raise ValueError(f'Publication triplicate without exactly one warmup: {key}')
        complete[key] = sorted(group, key=lambda r: r['repetition'])

    performance, bulk, single_cell, indexes, mapping = [], [], [], [], []
    for (case, threads, method), group in sorted(complete.items()):
        median = lambda field: statistics.median(r[field] for r in group)
        values = {
            'case': case, 'method': method, 'threads': threads,
            'scope': group[0].get('scope', 'index' if case == 'index' else 'complete_pipeline'),
            'layout': 'index' if case == 'index' else cases[case].get('layout', 'paired'),
            'n': 3, 'wall_s': median('wall_seconds'),
            'wall_min_s': min(r['wall_seconds'] for r in group),
            'wall_max_s': max(r['wall_seconds'] for r in group),
            'cpu_s': median('cpu_seconds'), 'rss_gib': median('peak_rss_bytes') / 2**30,
            'rss_min_gib': min(r['peak_rss_bytes'] for r in group) / 2**30,
            'rss_max_gib': max(r['peak_rss_bytes'] for r in group) / 2**30,
        }
        if case != 'index':
            records = cases[case].get('pairs', cases[case].get('reads'))
            values['input_records'] = records
            values['million_records_per_second'] = records / values['wall_s'] / 1e6
        else:
            values['index_gib'] = statistics.median(r['index_disk']['logical_bytes'] for r in group) / 2**30
        performance.append(values)
        paper = {
            'Dataset': case, 'Method': LABELS[method], 'Workers': threads,
            'Wall, s': f"{values['wall_s']:.3f}",
            'Wall range, s': _range([r['wall_seconds'] for r in group]),
            'CPU, s': f"{values['cpu_s']:.2f}",
            'RSS, GiB': f"{values['rss_gib']:.3f}",
            'RSS range, GiB': _range([r['peak_rss_bytes'] for r in group], 2**30),
        }
        if case == 'index':
            paper.pop('Dataset'); paper.pop('Workers')
            paper['Index, GiB'] = f"{values['index_gib']:.3f}"
            indexes.append(paper)
        elif cases[case]['kind'] == 'bulk':
            bulk.append(paper)
        elif values['scope'] != 'mapping_only':
            single_cell.append(paper)
        probes = [r.get('mapping_probe') for r in group]
        if all(probes):
            seq = probes
            cpu = lambda s: s['cpu_seconds']
        elif case != 'index':
            seq = [next(s for s in r['stages'] if s['stage'] == 'mapping') for r in group]
            cpu = lambda s: s['process_user_seconds'] + s['process_system_seconds']
        else:
            continue
        mapping.append({
            'Dataset': case, 'Method': LABELS[method], 'Workers': threads,
            'Wall, s': f"{statistics.median(s['wall_seconds'] for s in seq):.3f}",
            'Wall range, s': _range([s['wall_seconds'] for s in seq]),
            'CPU, s': f"{statistics.median(cpu(s) for s in seq):.2f}",
            'RSS, GiB': f"{statistics.median(s['peak_rss_bytes'] for s in seq) / 2**30:.3f}",
            'RSS range, GiB': _range([s['peak_rss_bytes'] for s in seq], 2**30),
        })

    @lru_cache(maxsize=8)
    def bulk_values(case, threads, method):
        return load_bulk(complete[(case, threads, method)][0]['output'], method, reference)

    thread_agreement = []
    for case in cases.values():
        methods = experiment['design']['bulk_methods' if case['kind'] == 'bulk' else 'single_cell_methods']
        for method in methods:
            keys = [(case['id'], t, method) for t in (1, 8)]
            if not all(k in complete and complete[k][0].get('scope') != 'mapping_only' for k in keys):
                continue
            base = {'case': case['id'], 'method': method, 'threads_a': 1, 'threads_b': 8, 'repetition': 1}
            if case['kind'] == 'bulk':
                a, b = (bulk_values(*k) for k in keys)
                thread_agreement.extend({**base, 'metric': k, **differences(a[k], b[k])} for k in a)
            else:
                a, b = (load_sc(complete[k][0]['output'], method) for k in keys)
                metrics = sc_agreement(a, b)
                thread_agreement.append({**base, 'metric': 'gene_umi',
                    **{k: v for k, v in metrics.items() if not isinstance(v, dict)}})
                write_json(out / f"{case['id']}-thread-agreement.json", metrics)

    with Path(reference['transcripts']['path']).open() as handle:
        universe = list(csv.DictReader(handle, delimiter='\t'))
    genes = sorted({r['gene_id'] for r in universe})
    gi = {name: i for i, name in enumerate(genes)}
    tx_to_gene = np.array([gi[r['gene_id']] for r in universe])
    accuracy = []
    for case in cases.values():
        if not case.get('truth'):
            continue
        with np.load(case['truth']) as truth:
            targets = {'transcript_counts': truth['fragment_counts'], 'transcript_tpm': truth['molecule_tpm']}
            targets.update({name.replace('transcript_', 'gene_'):
                np.bincount(tx_to_gene, weights=value, minlength=len(genes)) for name, value in list(targets.items())})
        for threads in (1, 8):
            for method in experiment['design']['bulk_methods']:
                if (case['id'], threads, method) not in complete:
                    continue
                estimates = bulk_values(case['id'], threads, method)
                accuracy.extend({'case': case['id'], 'method': method, 'threads': threads,
                    'metric': metric, **differences(estimates[metric], target)} for metric, target in targets.items())

    reference_name = 'GENCODE v50 ALL' if reference['kind'] == 'gencode50' else 'Official Pachter lab mature-plus-nascent (Ensembl 108)'
    allocation = []
    if reference['kind'] == 'official-nac':
        for (case, threads, method), group in sorted(complete.items()):
            if case == 'index':
                continue
            values = bulk_values(case, threads, method)
            mature, nascent = (float(values[state + '_counts'].sum()) for state in ('mature', 'nascent'))
            allocation.append({'case': case, 'method': method, 'threads': threads, 'repetition': 1,
                'mature_estimated_fragments': mature, 'nascent_estimated_fragments': nascent,
                'nascent_fraction_estimated_fragments': nascent / (mature + nascent) if mature + nascent else None})
    table(out, 'nascent-allocation', allocation)
    table(out, 'publication-performance', performance)
    _paper_table(out, 'paper-index', indexes, reference_name + ' index construction at eight workers on Apple Silicon.', index=True)
    _paper_table(out, 'paper-bulk-paired', [r for r in bulk if cases[r['Dataset']].get('layout', 'paired') != 'single'],
        'Prepared paired-end FASTQ to transcript abundance on Apple Silicon.')
    _paper_table(out, 'paper-bulk-single', [r for r in bulk if cases[r['Dataset']].get('layout') == 'single'],
        'Prepared single-end R1 FASTQ to transcript abundance on Apple Silicon; common empirical fragment models.')
    _paper_table(out, 'paper-single-cell', single_cell,
        'Prepared 10x PBMC FASTQ to gene-by-barcode count matrix on Apple Silicon.')
    _paper_table(out, 'paper-mapping', mapping,
        'Native mapping workflows including loading, decompression and intermediate output; output workloads differ.')
    table(out, 'thread-agreement', thread_agreement)
    table(out, 'accuracy', accuracy)
    write_json(out / 'publication-summary.json', {
        'performance': performance, 'thread_agreement': thread_agreement, 'accuracy': accuracy,
        'scope': 'Technical repetitions on one Apple Silicon host; agreement is not biological accuracy.',
        'nascent_allocation': allocation,
        'complete_groups': len(complete), 'measured_trials': sum(len(g) for g in complete.values()),
        'audited_read_trials': audited_inputs,
        'reference': {k: reference[k] for k in ('kind', 'targets', 'genes', 'bases')},
        'reference_sha256': reference['fasta']['sha256'],
        'tools': {k: {'sha256': v['sha256'], 'bytes': v['bytes']} for k, v in experiment['tools'].items()},
        'measurement_workflow': experiment['workflow'],
        'analysis_workflow': analysis_workflow,
    })
