#!/usr/bin/env python3
"""Publish the two completed modern-reference cohorts without exposing host paths."""
import argparse
import csv
import json
from pathlib import Path
import shutil
from benchmark_common import fingerprint


def read(path):
    return json.loads(path.read_text())


def markdown(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join('---' for _ in headers) + ' |',
                      *['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows]])


def validated(workspace):
    summary = read(workspace / 'tables/publication-summary.json')
    spec = read(workspace / 'runs/experiment.json')
    expected = len(spec['design']['index_threads']) * len(spec['design']['builders'])
    for case in spec['inputs']['cases']:
        single_cell = case['kind'] == 'single_cell'
        workers = spec['design'].get('single_cell_threads', [1, 8]) if single_cell else spec['design']['threads']
        methods = spec['design']['single_cell_methods' if single_cell else 'bulk_methods']
        expected += len(workers) * len(methods)
    rows = [read(p) for p in (workspace / 'runs').glob('**/result.json')]
    if summary['complete_groups'] != expected or len(rows) != expected * 4:
        raise ValueError(f'Publication requires a complete warmup-plus-triplicate panel: {workspace}')
    executables = {item['path']: item for item in spec['tools'].values()}
    for item in executables.values():
        if fingerprint(item['path']) != item:
            raise ValueError('Publication executable differs from the frozen experiment')
    for row in rows:
        if row['case'] != 'index':
            case = next(c for c in spec['inputs']['cases'] if c['id'] == row['case'])
            if case['kind'] == 'single_cell' and row['threads'] not in spec['design']['single_cell_threads']:
                raise ValueError('Single-cell workers differ from the experiment design')
        stages = row['stages'] + row.get('mapping_probe', {}).get('stages', [])
        if any(s['command'][0] not in executables for s in stages):
            raise ValueError('A result used an executable outside the frozen experiment')
        if any(s.get('timing_interrupted') for s in stages):
            raise ValueError('Interrupted timing is not a valid publication observation')
    return summary, spec


def publish(gencode, nascent, output, tables):
    cohorts = {name: (root, *validated(root)) for name, root in [('gencode', gencode), ('nascent', nascent)]}
    tables.mkdir(parents=True, exist_ok=True)
    common = ['publication-performance', 'paper-index', 'paper-bulk-paired', 'paper-bulk-single',
              'paper-mapping', 'stages', 'agreement', 'accuracy', 'thread-agreement']
    copied = []
    exclusions = []
    refreshes = []
    for name, (root, _, _) in cohorts.items():
        refresh = root / 'audit/bus-dispatch-refresh.json'
        if refresh.exists():
            refreshes.append({'reference': name, **read(refresh)})
        for path in sorted((root / 'audit/excluded').glob('**/result.json')):
            original = read(path)
            key = path.parent.relative_to(root / 'audit/excluded')
            replacement = read(root / 'runs' / key / 'result.json')
            exclusions.append({'reference': name, 'trial': str(key),
                'original_wall_s': original['wall_seconds'],
                'replacement_wall_s': replacement['wall_seconds'],
                'reason': replacement['replacement_for']['reason']})
    for name, (root, summary, spec) in cohorts.items():
        names = common + (['paper-single-cell', 'single-cell-counts'] if name == 'gencode' else ['nascent-allocation'])
        for stem in names:
            for suffix in ('.tsv', '.tex'):
                src = root / 'tables' / (stem + suffix)
                if not src.exists():
                    if stem == 'accuracy' and name == 'nascent':
                        continue
                    raise ValueError(f'Missing publication table: {src}')
                dst = tables / (name + '-' + src.name)
                shutil.copy2(src, dst); copied.append(dst.name)
        for src in sorted((root / 'tables').glob('*.json')):
            if src.name not in ('repeat-variability.json', 'processed-input-audit.json') and 'matrix-agreement' not in src.name and 'thread-agreement' not in src.name:
                continue
            dst = tables / (name + '-' + src.name)
            shutil.copy2(src, dst); copied.append(dst.name)

    perf = [(name, p) for name, (_, s, _) in cohorts.items() for p in s['performance']]
    labels = {'kallistors': 'kallistors', 'salmon': 'Salmon', 'piscem': 'piscem',
              'salmon-sketch': 'Salmon sketch', 'salmon-sa': 'Salmon selective alignment',
              'piscem-infer': 'piscem + piscem-infer', 'kallistors-bustools': 'kallistors + bustools',
              'piscem-alevin-fry': 'piscem + alevin-fry'}
    titles = {'gencode': 'GENCODE v50 ALL', 'nascent': 'Official mature + nascent'}
    lines = ['# Modern quantifier benchmark results', '',
        'Apple Silicon; pinned native binaries; one warmup and three measured runs per condition. '
        'All times and peak RSS values below are medians. Index construction uses eight workers only. '
        'Bulk and single-cell mapping use one and eight workers; alevin-fry complete processing uses eight only. '
        'Cell Ranger execution is deferred.', '',
        'Read times cover prepared FASTQ to native abundance or count matrix. Downloading, trimming '
        'and reference preparation are outside the timers. See the [methods](modern-benchmarks.md) '
        'for exact settings, algorithm differences and limitations.', '', '## Index construction', '']
    rows = []
    for name, p in perf:
        if p['case'] == 'index':
            rows.append([titles[name], labels[p['method']], f"{p['wall_s']:.3f}",
                f"{p['wall_min_s']:.3f}–{p['wall_max_s']:.3f}", f"{p['cpu_s']:.2f}",
                f"{p['rss_gib']:.3f}", f"{p['index_gib']:.3f}"])
    lines += [markdown(['Reference', 'Builder', 'Wall, s', 'Range, s', 'CPU, s', 'Peak RSS, GiB', 'Index, GiB'], rows), '']
    for layout, heading in [('paired', 'Paired-end quantification'), ('single', 'Single-end quantification')]:
        lines += ['## ' + heading, '']
        rows = []
        for name, (_, summary, spec) in cohorts.items():
            cases = {c['id']: c for c in spec['inputs']['cases']}
            values = {(p['case'], p['threads'], p['method']): p for p in summary['performance']}
            for case in cases.values():
                if case['kind'] != 'bulk' or case.get('layout', 'paired') != layout:
                    continue
                for method in spec['design']['bulk_methods']:
                    a, b = (values[(case['id'], t, method)] for t in (1, 8))
                    rows.append([case['id'], labels[method], f"{a['wall_s']:.3f}", f"{b['wall_s']:.3f}",
                                 f"{a['rss_gib']:.3f}", f"{b['rss_gib']:.3f}"])
        lines += [markdown(['Dataset', 'Method', 'Wall at 1, s', 'Wall at 8, s', 'RSS at 1, GiB', 'RSS at 8, GiB'], rows), '']
    lines += ['Single-end observations use R1 from the paired library. All methods receive the same '
              'frozen paired-end fragment mean and SD for each real library. Simulation fragment lengths '
              'are known by construction. These layouts are related inputs, not biological replicates.', '',
              '## Single-cell processing', '']
    root, summary, spec = cohorts['gencode']
    with (root / 'tables/single-cell-counts.tsv').open() as handle:
        cells = {(r['method'], int(r['threads'])): r for r in csv.DictReader(handle, delimiter='\t')}
    rows = []
    for p in summary['performance']:
        if p['case'] != 'pbmc_1k_v3' or p['scope'] == 'mapping_only':
            continue
        c = cells[(p['method'], p['threads'])]
        rows.append([labels[p['method']], p['threads'], f"{p['wall_s']:.3f}", f"{p['cpu_s']:.2f}",
            f"{p['rss_gib']:.3f}", c['common_gate_cells'], c['median_umis_called'], c['median_genes_called']])
    lines += [markdown(['Pipeline', 'Workers', 'Wall, s', 'CPU, s', 'Peak RSS, GiB',
                        'Cells passing gate', 'Median UMIs', 'Median detected genes'], rows), '',
        'The common descriptive gate is at least 500 UMIs and 200 detected genes. '
        'Native retained barcodes are not a biological cell count. Barcode correction, multimapper '
        'handling and UMI resolution differ between bustools and alevin-fry.', '']
    sc = read(root / 'tables/pbmc_1k_v3-t8-matrix-agreement.json')
    lines += [f"At eight workers, **{sc['common_called_cells']:,} cells** pass the gate in both workflows; "
        f"called-cell Jaccard is **{sc['called_cell_jaccard']:.5f}**. Matrix differences are reported "
        'over the full barcode/gene union and the common-called subset.', '']
    rows = []
    for label, suffix, maxkey in [('Full nonzero union', 'nonzero_union', 'max_abs'),
                                  ('Common-called nonzero union', 'common_called_nonzero_union', 'common_called_matrix_max_abs')]:
        rows.append([label, *[f"{sc[k + '_' + suffix]:.6g}" for k in ('p99', 'p99_9', 'p99_99')], f"{sc[maxkey]:.6g}"])
    lines += [markdown(['Matrix entries', 'p99', 'p99.9', 'p99.99', 'Maximum absolute difference'], rows), '',
        '## Independent simulation accuracy', '',
        'The following summaries use measured replicate 1 at eight workers, selected before observing '
        'performance. Known fragment origins and nominal molecule TPM come from the independent '
        'simulation. Full percentile, correlation and one-worker results are in the supplementary tables.', '']
    acc = {(r['case'], r['method'], r['threads'], r['metric']): r for r in summary['accuracy']}
    rows = []
    for case in ('fragment-truth', 'fragment-truth-se'):
        for method in spec['design']['bulk_methods']:
            for level in ('transcript', 'gene'):
                a = acc[(case, method, 8, level + '_counts')]; b = acc[(case, method, 8, level + '_tpm')]
                rows.append([case, labels[method], level, f"{a['mae']:.6g}", f"{a['max_abs']:.6g}",
                    f"{b['mae']:.6g}", f"{b['pearson_log1p']:.6g}"])
    lines += [markdown(['Input', 'Method', 'Level', 'Count MAE', 'Count max abs.', 'TPM MAE', 'TPM log1p Pearson'], rows), '',
        'This simulation omits indels, polymorphism, contamination, PCR and sequence biases. '
        'It cannot establish a universal biological accuracy ranking. Real-library and single-cell '
        'comparisons measure agreement only.', '', '## Mapping, agreement and reproducibility', '',
        '* Separate native mapping tables include index loading, decompression and intermediate output. '
        'Aggregate EC, RAD and BUS outputs have different workloads; these timings are not pure mapping-kernel timings.',
        '* Full pairwise agreement, replicate variability and within-method one-versus-eight-worker '
        'agreement include maximum absolute differences, p99, p99.9 and p99.99. Gene-by-state and '
        'mature/nascent subsets are included for TT-seq.',
        '* Wall-time and RSS ranges, CPU time, pipeline stages, retained-barcode and UMI summaries, '
        'and processed-input audits are in the [publication tables](modern-benchmark-tables/README.md).',
        '* Native processed-record counters match each frozen input count. These are technical '
        'triplicates on one host and two distinct references; no cross-host or biological replication is implied.',
        '* Salmon selective alignment uses a transcript-only reference without genomic decoys. '
        'The TT-seq reference is the official Ensembl 108 mature-plus-nascent set, not an augmented GENCODE v50 set.', '']
    if refreshes:
        lines += ['### Executable refresh', '',
            'After fixing stranded BUS dispatch, all kallistors measurements were recollected with the '
            'corrected executable. Unchanged competitor controls were retained after verifying identical '
            'commands, executables, inputs and accounting code. Earlier kallistors measurements are excluded; '
            'their original records and workflow snapshots are retained separately.', '']
    if exclusions:
        lines += ['### Audited exclusions', '']
        for item in exclusions:
            lines += [f"* `{item['reference']}/{item['trial']}`: {item['reason']} "
                f"The excluded observation was {item['original_wall_s']:.3f} s; the identical-command "
                f"replacement was {item['replacement_wall_s']:.3f} s. The original record is retained. "
                'This was a setup-overlap exclusion, not a timing-based outlier rule.']
        lines += ['']
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('\n'.join(lines))
    portable = {name: {**s, 'source_sha256': read(root / 'environment.json')['source_sha256']}
                for name, (root, s, _) in cohorts.items()}
    portable['exclusions'] = exclusions
    portable['executable_refreshes'] = refreshes
    output.with_suffix('.json').write_text(json.dumps(portable, indent=2) + '\n')
    (tables / 'README.md').write_text('# Modern benchmark publication tables\n\n'
        'TSV files retain numeric precision; booktabs LaTeX files provide typeset tables. '
        'The `paper-*` tables use longtable and booktabs; bracketed ranges are observed minima and maxima. '
        'GENCODE and nascent reference panels are separate. Peak RSS is the maximum sequential stage peak. '
        'All reported triplicates require one warmup. See [results](../modern-benchmark-results.md) '
        'and [methods](../modern-benchmarks.md) for scope and interpretation.\n\n' +
        '\n'.join(f'* [{name}]({name})' for name in sorted(copied)) + '\n')
    print(f'Published {len(copied)} table and audit artifacts; {sum(s[1]["measured_trials"] for s in cohorts.values())} measured primary trials.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gencode', type=Path, required=True)
    parser.add_argument('--nascent', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('docs/modern-benchmark-results.md'))
    parser.add_argument('--tables', type=Path, default=Path('docs/modern-benchmark-tables'))
    args = parser.parse_args()
    publish(args.gencode.resolve(), args.nascent.resolve(), args.output, args.tables)
