"""Native-f64 agreement and publication tables for the fresh two-tool controls."""
import csv
from functools import lru_cache
from pathlib import Path
import shutil
import statistics

import h5py
import numpy as np

from benchmark_common import write_json
from benchmark_common import digest
from modern_benchmark.environment import read_json
from modern_benchmark.metrics import differences
from modern_benchmark.tables import table


@lru_cache(maxsize=2)
def reference_labels(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = np.array([r['transcript_id'] for r in rows])
    genes = sorted({r['gene_id'] for r in rows})
    index = {name: i for i, name in enumerate(genes)}
    return ids, np.array([index[r['gene_id']] for r in rows]), len(genes)


@lru_cache(maxsize=1)
def state_labels(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    states = sorted({(r['gene_id'], r['biotype']) for r in rows})
    index = {key: i for i, key in enumerate(states)}
    return (np.array([index[(r['gene_id'], r['biotype'])] for r in rows]), len(states),
            {state: np.array([r['biotype'] == state for r in rows]) for state in ('mature', 'nascent')})


def native_values(out, reference):
    ids, to_gene, genes = reference_labels(reference['transcripts']['path'])
    with h5py.File(out / 'abundance.h5') as handle:
        names = handle['aux/ids'].asstr()[:]
        counts = handle['est_counts'][:]
        lengths = handle['aux/eff_lengths'][:]
    if not np.array_equal(names, ids):
        lookup = {name: i for i, name in enumerate(names)}
        if len(lookup) != len(names) or set(lookup) != set(ids):
            raise ValueError('Native abundance target universe differs')
        order = [lookup[name] for name in ids]
        counts, lengths = counts[order], lengths[order]
    if (not np.isfinite(counts).all() or not np.isfinite(lengths).all()
            or (counts < 0).any() or (lengths < 0).any() or ((lengths == 0) & (counts > 0)).any()):
        raise ValueError('Invalid native abundance values')
    rates = np.divide(counts, lengths, out=np.zeros_like(counts), where=lengths > 0)
    tpm = rates * (1e6 / rates.sum()) if rates.sum() else rates
    result = {'transcript_counts': counts, 'transcript_tpm': tpm, 'effective_length': lengths,
            'gene_counts': np.bincount(to_gene, weights=counts, minlength=genes),
            'gene_tpm': np.bincount(to_gene, weights=tpm, minlength=genes)}
    if reference['kind'] == 'official-nac':
        to_state, states, masks = state_labels(reference['transcripts']['path'])
        result.update(gene_state_counts=np.bincount(to_state, weights=counts, minlength=states),
                      gene_state_tpm=np.bincount(to_state, weights=tpm, minlength=states))
        for state, mask in masks.items():
            result[state + '_counts'] = counts[mask]
            result[state + '_tpm'] = tpm[mask]
    return result


def fragment_model(out):
    """Read the realized model, separately from supplied single-end parameters."""
    with h5py.File(out / 'abundance.h5') as handle:
        histogram = np.asarray(handle['aux/fld'][:], dtype=float)
    if histogram.ndim != 1 or not np.isfinite(histogram).all() or (histogram < 0).any():
        raise ValueError('Invalid native fragment histogram')
    total = float(histogram.sum()) if histogram.size else None
    bins = np.arange(len(histogram), dtype=float)
    mean = float(np.dot(histogram, bins) / total) if total else None
    sd = float(np.sqrt(np.dot(histogram, (bins - mean) ** 2) / total)) if total else None
    return histogram, {'histogram_bins': len(histogram), 'histogram_serialized': bool(histogram.size),
                       'histogram_total': total, 'realized_mean': mean, 'realized_sd': sd}


def compare_histograms(a, b, layout):
    """An omitted single-end histogram is unavailable, not a zero distribution."""
    if not a.size or not b.size:
        if layout != 'single':
            raise ValueError('Paired-end empirical fragment histogram is missing')
        return {'histogram_comparison': 'not comparable: native histogram omitted',
                'histogram_exact': None, 'histogram_l1_difference': None}
    if a.shape != b.shape:
        raise ValueError('Unexpected nonempty fragment histogram binning difference')
    return {'histogram_comparison': 'same one-base bin coordinates',
            'histogram_exact': bool(np.array_equal(a, b)),
            'histogram_l1_difference': float(np.abs(a - b).sum())}


def report(workspace, spec):
    out = workspace / 'tables'
    out.mkdir(exist_ok=True)
    scripts = Path(__file__).resolve().parents[1]
    analysis_workflow = {relative: digest(scripts / relative) for relative in spec['workflow']}
    from modern_benchmark.environment import content_id
    revision_id = content_id(analysis_workflow)
    revision = workspace / 'analysis-workflow-revisions' / revision_id
    for relative in analysis_workflow:
        target = revision / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copy2(scripts / relative, target)
    write_json(out / 'analysis-manifest.json', {'measurement_workflow': spec['workflow'],
        'analysis_workflow': analysis_workflow, 'analysis_revision': revision_id,
        'native_measurements_recollected': False})
    performance, agreement, thread_agreement, accuracy, mapping, cell_counts, phases = [], [], [], [], [], [], []
    repeated_agreement, replicate_variation, models, largest_changes = [], [], [], []
    def summarize(group, reference, case, layout, threads, tool, scope):
        if len(group) != 3:
            raise ValueError('Publication requires exactly three measured runs')
        if any(r.get('timing_interrupted') for r in group):
            raise ValueError('Interrupted primary timing')
        cpu = lambda r: r.get('cpu_seconds', r.get('process_user_seconds', 0) + r.get('process_system_seconds', 0))
        return {'reference': reference, 'case': case, 'layout': layout, 'threads': threads, 'tool': tool,
            'scope': scope, 'n': 3, 'wall_s': statistics.median(r['wall_seconds'] for r in group),
            'wall_min_s': min(r['wall_seconds'] for r in group), 'wall_max_s': max(r['wall_seconds'] for r in group),
            'cpu_s': statistics.median(cpu(r) for r in group),
            'rss_gib': statistics.median(r['peak_rss_bytes'] for r in group) / 2**30,
            'rss_min_gib': min(r['peak_rss_bytes'] for r in group) / 2**30,
            'rss_max_gib': max(r['peak_rss_bytes'] for r in group) / 2**30}
    for name, panel in spec['panels'].items():
        index_rows = read_json(workspace / name / 'index/runs.json')
        for tool in ('kallisto', 'kallistors'):
            group = [r for r in index_rows if r['phase'] == 'measured' and r['tool'] == tool]
            row = summarize(group, name, 'index', 'index', 8, tool, 'index')
            row['index_gib'] = statistics.median(r['index_bytes'] for r in group) / 2**30
            performance.append(row)
        root = workspace / name / 'bulk'
        if not read_json(root / 'summary.json')['complete']:
            raise ValueError('Incomplete bulk panel')
        configs = read_json(root / 'protocol.json')['configuration']['cases']
        rows = read_json(root / 'runs.json')
        for c in configs:
            for tool in ('kallisto', 'kallistors'):
                group = [r for r in rows if r['case'] == c['id'] and r['tool'] == tool and r['phase'] == 'measured']
                performance.append(summarize(group, name, c['name'], c['mode'], c['threads'], tool, 'complete_quantification'))
                if all(r.get('stages_seconds') for r in group):
                    phases.extend({'reference': name, 'case': c['name'], 'layout': c['mode'], 'threads': c['threads'],
                        'tool': tool, 'phase': phase, 'wall_s': statistics.median(r['stages_seconds'][phase] for r in group),
                        'scope': 'CLI message boundary estimate; model_and_em is a derived aggregate'}
                        for phase in group[0]['stages_seconds'])
            original, candidate = (native_values(root / 'cases' / c['id'] / f'measured-1-{tool}.output', panel['reference'])
                                   for tool in ('kallisto', 'kallistors'))
            agreement.extend({'reference': name, 'case': c['name'], 'layout': c['mode'], 'threads': c['threads'],
                'metric': metric, 'precision': 'native HDF5 f64', **differences(original[metric], candidate[metric])}
                for metric in original)
            base = {'reference': name, 'case': c['name'], 'layout': c['mode'], 'threads': c['threads']}
            labels, _, _ = reference_labels(panel['reference']['transcripts']['path'])
            for replicate in (1, 2, 3):
                paths = {tool: root / 'cases' / c['id'] / f'measured-{replicate}-{tool}.output'
                         for tool in ('kallisto', 'kallistors')}
                estimates = {'kallisto': original, 'kallistors': candidate} if replicate == 1 else {
                    tool: native_values(path, panel['reference']) for tool, path in paths.items()}
                for metric, a in estimates['kallisto'].items():
                    b = estimates['kallistors'][metric]
                    for universe, mask in (('all_features', np.ones(len(a), dtype=bool)),
                                           ('nonzero_union', (a != 0) | (b != 0))):
                        if mask.any():
                            repeated_agreement.append({**base, 'replicate': replicate, 'metric': metric,
                                'universe': universe, 'precision': 'native HDF5 f64', **differences(a[mask], b[mask])})
                    if replicate == 1 and metric in ('transcript_counts', 'transcript_tpm'):
                        delta = np.abs(a - b)
                        order = np.argsort(-delta, kind='stable')[:25]
                        largest_changes.extend({**base, 'replicate': 1, 'metric': metric, 'rank': rank,
                            'transcript_id': str(labels[index]), 'original': float(a[index]),
                            'kallistors': float(b[index]), 'abs_difference': float(delta[index])}
                            for rank, index in enumerate(order, 1) if delta[index] > 0)
                if replicate != 1:
                    for tool, baseline in (('kallisto', original), ('kallistors', candidate)):
                        replicate_variation.extend({**base, 'tool': tool, 'replicate_a': 1,
                            'replicate_b': replicate, 'metric': metric,
                            **differences(baseline[metric], estimates[tool][metric])} for metric in baseline)
                histograms, descriptions = {}, {}
                for tool, path in paths.items():
                    histograms[tool], descriptions[tool] = fragment_model(path)
                a, b = (histograms[tool] for tool in ('kallisto', 'kallistors'))
                models.append({**base, 'replicate': replicate,
                    **compare_histograms(a, b, c['mode']),
                    'histogram_interpretation': 'empirical observations' if c['mode'] == 'paired' else
                        'configured Gaussian model; serialized bins are not empirical observations',
                    'supplied_gaussian_mean': c['fragment_length'] if c['mode'] == 'single' else None,
                    'supplied_gaussian_sd': c['fragment_sd'] if c['mode'] == 'single' else None,
                    **{tool + '_' + key: value for tool, description in descriptions.items()
                       for key, value in description.items()}})
            print('Analysed native abundance agreement:', name, c['id'], flush=True)
            source_case = next(case for case in panel['cases'] if case['id'] == c['name'])
            if source_case.get('truth'):
                _, to_gene, n_genes = reference_labels(panel['reference']['transcripts']['path'])
                with np.load(source_case['truth']) as truth:
                    targets = {'transcript_counts': truth['fragment_counts'], 'transcript_tpm': truth['molecule_tpm']}
                    targets.update({key.replace('transcript_', 'gene_'): np.bincount(to_gene, weights=value, minlength=n_genes)
                                    for key, value in list(targets.items())})
                for tool, estimates in [('kallisto', original), ('kallistors', candidate)]:
                    accuracy.extend({'case': c['name'], 'layout': c['mode'], 'threads': c['threads'], 'tool': tool,
                        'metric': metric, **differences(estimates[metric], target)} for metric, target in targets.items())
        by_key = {(c['name'], c['mode'], c['threads']): c for c in configs}
        for case in panel['cases']:
            for layout in ('paired', 'single'):
                for tool in ('kallisto', 'kallistors'):
                    values = [native_values(root / 'cases' / by_key[(case['id'], layout, t)]['id'] /
                              f'measured-1-{tool}.output', panel['reference']) for t in (1, 8)]
                    thread_agreement.extend({'reference': name, 'case': case['id'], 'layout': layout, 'tool': tool,
                        'threads_a': 1, 'threads_b': 8, 'metric': metric, **differences(values[0][metric], values[1][metric])}
                        for metric in values[0])
    parity = read_json(workspace / 'single-cell/parity.json')
    if not parity['complete'] or not parity['all_full_libraries_exact']:
        raise ValueError('Single-cell full parity not proven')
    sc_rows = [read_json(p) for p in (workspace / 'single-cell').glob('**/result.json')]
    for threads in spec['design']['single_cell_workers']:
        for tool in ('kallisto', 'kallistors'):
            group = sorted([r for r in sc_rows if r['phase'] == 'measured' and r['threads'] == threads and r['method'] == tool],
                           key=lambda r: r['repetition'])
            performance.append(summarize(group, 'gencode', 'pbmc_1k_v3', 'single_cell', threads, tool, 'prepared_FASTQ_to_matrix'))
            stages = [next(s for s in r['stages'] if s['stage'] == 'mapping') for r in group]
            mapping.append(summarize(stages, 'gencode', 'pbmc_1k_v3', 'single_cell', threads, tool, 'BUS_mapping'))
            cell_counts.append({'case': 'pbmc_1k_v3', 'threads': threads, 'tool': tool,
                **group[0]['matrix_summary'], 'max_abs': 0, 'p99': 0, 'p99_9': 0, 'p99_99': 0,
                'exact_canonical_matrix': True})
    for stem, values in [('performance', performance), ('agreement-native', agreement), ('thread-agreement-native', thread_agreement),
                         ('accuracy-native', accuracy), ('single-cell-mapping', mapping), ('single-cell-counts', cell_counts),
                         ('bulk-phases', phases), ('agreement-all-replicates-native', repeated_agreement),
                         ('replicate-variation-native', replicate_variation), ('fragment-models', models),
                         ('largest-transcript-differences', largest_changes)]:
        table(out, stem, values)
    from modern_benchmark.publication import _paper_table
    def paper_row(row):
        return {'Dataset': row['reference'] + ': ' + row['case'],
            'Method': 'original kallisto' if row['tool'] == 'kallisto' else 'kallistors',
            'Workers': row['threads'], 'Wall, s': f"{row['wall_s']:.3f}",
            'Wall range, s': f"{row['wall_min_s']:.3f}–{row['wall_max_s']:.3f}",
            'CPU, s': f"{row['cpu_s']:.2f}", 'RSS, GiB': f"{row['rss_gib']:.3f}",
            'RSS range, GiB': f"{row['rss_min_gib']:.3f}–{row['rss_max_gib']:.3f}"}
    for reference in spec['panels']:
        rows = []
        for row in performance:
            if row['layout'] == 'index' and row['reference'] == reference:
                item = paper_row(row)
                item['Index, GiB'] = f"{row['index_gib']:.3f}"
                rows.append(item)
        _paper_table(out, 'paper-index-' + reference, rows,
                     reference + ' index construction at eight workers on Apple Silicon.', index=True)
    for layout in ('paired', 'single', 'single_cell'):
        _paper_table(out, 'paper-' + layout.replace('_', '-'),
                     [paper_row(row) for row in performance if row['layout'] == layout],
                     layout.replace('_', ' ') + ' complete processing on Apple Silicon.')
    _paper_table(out, 'paper-single-cell-mapping', [paper_row(row) for row in mapping],
                 'Single-cell BUS mapping including index loading and intermediate output on Apple Silicon.')
    result = {'complete': True, 'source_precision': 'native HDF5 f64; TPM recomputed before text rounding',
        'measurement_workflow': spec['workflow'], 'analysis_workflow': analysis_workflow,
        'analysis_revision': revision_id,
        'design': spec['design'], 'kallistors_source_sha256': spec['kallistors_source_sha256'],
        'tools': {name: {k: record[k] for k in ('bytes', 'sha256')} for name, record in spec['tools'].items()},
        'performance': performance, 'agreement': agreement, 'thread_agreement': thread_agreement,
        'agreement_all_replicates': repeated_agreement, 'replicate_variation': replicate_variation,
        'fragment_models': models, 'largest_transcript_differences': largest_changes,
        'accuracy': accuracy, 'single_cell_counts': cell_counts, 'single_cell_mapping': mapping, 'single_cell_parity': parity,
        'bulk_phase_estimates': phases}
    write_json(workspace / 'summary.json', result)
    docs = Path(__file__).resolve().parents[2] / 'docs'
    public = docs / 'kallisto-comparison-tables'
    public.mkdir(exist_ok=True)
    for source in out.iterdir():
        if source.suffix in ('.tsv', '.tex'):
            shutil.copy2(source, public / source.name)
    write_json(docs / 'kallisto-comparison-results.json', result)
    from publish_modern_benchmarks import markdown
    rows = []
    keys = sorted({(r['reference'], r['case'], r['layout'], r['threads']) for r in performance})
    for reference, case, layout, threads in keys:
        pair = {r['tool']: r for r in performance if (r['reference'], r['case'], r['layout'], r['threads']) == (reference, case, layout, threads)}
        a, b = pair['kallisto'], pair['kallistors']
        rows.append([reference, case, layout, threads, f"{a['wall_s']:.3f}", f"{b['wall_s']:.3f}",
                     f"{a['wall_s'] / b['wall_s']:.3f}×", f"{a['rss_gib']:.3f}", f"{b['rss_gib']:.3f}"])
    text = ['# Fresh comparison with original kallisto', '',
        'Apple Silicon; original kallisto 0.52.0 and the same frozen kallistors executable used in the modern-tool panels. '
        'All primary values are medians of three runs after one warmup. Builds use eight workers only; '
        'paired-end, R1-only and single-cell processing use one and eight workers.', '',
        markdown(['Reference', 'Dataset', 'Mode', 'Workers', 'Original, s', 'kallistors, s', 'Original/current time',
                  'Original RSS, GiB', 'kallistors RSS, GiB'], rows), '',
        'GENCODE v50 ALL is used for cDNA and PBMC; TT-seq uses only the official Pachter lab mature-plus-nascent '
        'FASTA assets. Strand settings follow metadata, and both tools receive identical empirical single-end '
        'fragment models and fastp-prepared inputs. The independent simulation is adapter-free and untrimmed.', '',
        'Full EC maps are checked before bulk timing. The primary read comparison uses identical original-built '
        'index bytes; both programs additionally quantify matched prefixes using the freshly kallistors-built '
        'indexes to check index compatibility. Prefix compatibility timings are not primary performance observations.', '',
        'The official release executable supplies original timing controls. The separate instrumented original is '
        'used only for full EC correctness exports. Every full single-cell run has exact weighted BUS and canonical '
        'gene/UMI matrix parity after the same bustools pipeline.', '',
        'Native HDF5 abundance comparisons include maximum absolute differences, p99, p99.9 and p99.99 at transcript '
        'and gene level for every measured replicate. Results include the complete feature universe and the '
        'nonzero union, so zeros cannot conceal differences among expressed features. The largest transcript '
        'differences are labelled in a separate table. Native fragment sampling and EM can produce abundance '
        'differences even when EC counts are exactly equal; realized fragment histograms, within-tool replicate '
        'variation and one-versus-eight-worker comparisons are retained separately. Real-data agreement does '
        'not establish which estimator is more accurate; the independent simulation supplies stated-model truth.', '',
        'Single-end inputs use the same supplied Gaussian mean and SD. Kallistors does not serialize that '
        'model into `aux/fld`, while original kallisto stores an integerized Gaussian histogram. Its absence '
        'is reported as unavailable; no zero histogram is substituted and no histogram equality or L1 '
        'difference is asserted. Native count, effective-length and TPM comparisons remain complete. '
        'Empirical paired-end histograms are compared on their common one-base bin coordinates.', '',
        'See the [complete numeric tables](kallisto-comparison-tables/README.md), '
        '[machine-readable results](kallisto-comparison-results.json), and [methods](modern-benchmarks.md).', '']
    (docs / 'kallisto-comparison-results.md').write_text('\n'.join(text))
    (public / 'README.md').write_text('# Original-kallisto comparison tables\n\n'
        'Warmup plus three measured runs. Wall/CPU seconds, peak RSS and index size are in `performance`. '
        'Agreement and accuracy use native HDF5 f64 values and include zeros. Single-cell mapping and complete '
        'pipelines are separate. See [results](../kallisto-comparison-results.md) for scope.\n\n' +
        '\n'.join(f'* [{p.name}]({p.name})' for p in sorted(public.iterdir()) if p.suffix in ('.tsv', '.tex')) + '\n')
