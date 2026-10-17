#!/usr/bin/env python3
"""Export complete, uninterrupted benchmark results as portable JSON and TSV."""

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

from benchmark_common import fingerprint, summarize, write_json
from compare_ec_counts import compare, load_counts


def portable(value):
    if isinstance(value, dict):
        return {key: (Path(item).name if key == 'path' else portable(item))
                for key, item in value.items() if key != 'pid'}
    if isinstance(value, list):
        return [Path(item).name if isinstance(item, str) and Path(item).is_absolute() else portable(item) for item in value]
    if isinstance(value, str) and Path(value).is_absolute():
        return Path(value).name
    return value


def central_stages(rows):
    ordered = sorted(rows, key=lambda r: r['wall_seconds'])
    middle = ordered[(len(ordered) - 1) // 2:len(ordered) // 2 + 1]
    return {key: statistics.mean(row['stages_seconds'][key] for row in middle)
            for key in middle[0]['stages_seconds']}


def target_partition(path, mature_targets):
    result = {key: {'nonzero_ecs': 0, 'counts': 0} for key in ('mature_only', 'nascent_only', 'mixed')}
    for ec, count in load_counts(path).items():
        kind = 'mature_only' if max(ec) < mature_targets else 'nascent_only' if min(ec) >= mature_targets else 'mixed'
        result[kind]['nonzero_ecs'] += 1
        result[kind]['counts'] += count
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark', type=Path, required=True)
    parser.add_argument('--json', type=Path, required=True)
    parser.add_argument('--tsv', type=Path, required=True)
    args = parser.parse_args()
    try:
        protocol = json.loads((args.benchmark / 'protocol.json').read_text())
        summary = json.loads((args.benchmark / 'summary.json').read_text())
        runs = json.loads((args.benchmark / 'runs.json').read_text())
        if protocol['schema_version'] != 1 or protocol['kind'] not in ('index', 'quantification') or not summary['complete'] or summary['kind'] != protocol['kind']:
            raise ValueError('Incomplete or unsupported benchmark')
        if args.json.resolve() == args.tsv.resolve():
            raise ValueError('JSON and TSV report destinations must differ')
        if args.json.exists() or args.tsv.exists():
            raise ValueError('Report destinations must be new')
        config = protocol['configuration']
        records = config['cases'] if protocol['kind'] == 'quantification' else [{'id': 'index', 'mode': 'index', 'threads': config['threads']}]
        if protocol['kind'] == 'quantification' and {c['id'] for c in summary['cases']} != {c['id'] for c in records}:
            raise ValueError('Incomplete case results')
        if {r.get('case', 'index') for r in runs} != {c['id'] for c in records}:
            raise ValueError('Unexpected or missing run configurations')
        rows = []
        for case in records:
            selected = [r for r in runs if r.get('case', 'index') == case['id']]
            if not selected or any(not math.isfinite(r['wall_seconds']) or r['wall_seconds'] <= 0 for r in selected):
                raise ValueError(f'Invalid run durations: {case["id"]}')
            result = summarize(selected)
            stored = next(c for c in summary['cases'] if c['id'] == case['id']) if protocol['kind'] == 'quantification' else summary
            if any(stored[key] != value for key, value in result.items()):
                raise ValueError(f'Summary disagrees with run records: {case["id"]}')
            if not result['all_measured_runs_uninterrupted']:
                raise ValueError(f'Interrupted timing: {case["id"]}; preserve it and run a new series')
            for tool in ('kallisto', 'kallistors'):
                measured = [r for r in selected if r['tool'] == tool and r['phase'] == 'measured']
                if sorted(r['replicate'] for r in measured) != list(range(1, config['runs'] + 1)):
                    raise ValueError(f'Missing or duplicated measured runs: {case["id"]}/{tool}')
                warmups = [r['replicate'] for r in selected if r['tool'] == tool and r['phase'] == 'warmup']
                start = 0 if protocol['kind'] == 'quantification' else 1
                if sorted(warmups) != list(range(start, config['warmups'] + start)):
                    raise ValueError(f'Missing or duplicated warmups: {case["id"]}/{tool}')
                stats = result['tools'][tool]
                row = {'case': case['id'], 'mode': case['mode'], 'threads': case['threads'], 'tool': tool,
                       'runs': len(measured), **{k: stats[k] for k in ('median_seconds', 'min_seconds', 'max_seconds',
                                                                   'sample_sd_seconds', 'median_peak_rss_bytes',
                                                                   'median_peak_memory_footprint_bytes')}}
                if all('stages_seconds' in r for r in measured):
                    stages = central_stages(measured)
                    row.update(startup=stages['startup_index'] + stages['input_setup'], reads=stages['read_processing'],
                               model_em=stages['model_and_em'], output=stages['output_shutdown'])
                rows.append(row)
            if protocol['kind'] == 'quantification':
                if not selected or not all(r['counters_match'] for r in selected):
                    raise ValueError(f'Quant counters differ: {case["id"]}')
                saved = next(c for c in summary['cases'] if c['id'] == case['id'])
                files = saved['ec']['files']
                paths = {'original': args.benchmark / 'correctness' / case['id'] / 'original/diagnostic.ec_counts.tsv',
                         'candidate': args.benchmark / 'correctness' / case['id'] / 'candidate.ec.tsv'}
                for tool, path in paths.items():
                    actual = fingerprint(path)
                    if any(actual[k] != files[tool][k] for k in ('sha256', 'bytes')):
                        raise ValueError(f'Changed EC evidence: {path}')
                ec = compare(load_counts(paths['original']), load_counts(paths['candidate']))
                if not ec['exact'] or ec != saved['ec']['comparison']:
                    raise ValueError(f'Complete EC parity failed: {case["id"]}')
                control = selected[0]['counters']
                if any(r['counters'] != control or r['counters']['n_pseudoaligned'] != ec['reference_total_counts'] for r in selected):
                    raise ValueError(f'Run counters disagree with complete EC evidence: {case["id"]}')
                for tool, audit in saved['validation']['hdf5'].items():
                    for name, meta in audit['files'].items():
                        path = args.benchmark / 'cases' / case['id'] / f'measured-{config["runs"]}-{tool}.output' / name
                        actual = fingerprint(path)
                        if any(actual[k] != meta[k] for k in ('sha256', 'bytes')):
                            raise ValueError(f'Changed validated quant output: {path}')
                reference = case.get('reference')
                if reference and reference.get('nascent_targets', 0):
                    saved['original_ec_target_partition'] = target_partition(paths['original'], reference['mature_targets'])
        fields = list(dict.fromkeys(key for row in rows for key in row))
        args.tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.tsv.open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        write_json(args.json, portable({'protocol': protocol, 'summary': summary, 'runs': runs}))
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f'Report export failed: {error}\n')


if __name__ == '__main__':
    main()
