"""Time complete identical bustools pipelines and audit every full BUS output."""
import hashlib
import json
from pathlib import Path

import numpy as np

from benchmark_bus import canonical_bus
from benchmark_common import fingerprint, write_json
from modern_benchmark.adapters import sc_plan
from modern_benchmark.artifacts import retire_intermediates
from modern_benchmark.experiment import run_one
from modern_benchmark.metrics import load_sc, sc_agreement, sc_summary


def matrix_identity(item):
    """Order-independent matrix/label identity with full integer UMI counts."""
    genes, cells, matrix = item
    gi, ci = np.argsort(genes), np.argsort(cells)
    canonical = matrix[gi, :][:, ci].tocsr()
    canonical.sort_indices()
    canonical.sum_duplicates()
    canonical.eliminate_zeros()
    if np.any(canonical.data != np.rint(canonical.data)):
        raise ValueError('Expected integer gene/UMI counts from bustools')
    digest = hashlib.sha256()
    for names, order in ((genes, gi), (cells, ci)):
        digest.update(json.dumps([names[i] for i in order], separators=(',', ':')).encode())
        digest.update(b'\0')
    digest.update(json.dumps(canonical.shape).encode())
    for values in (canonical.indptr, canonical.indices, canonical.data):
        digest.update(np.asarray(values, dtype='<i8').tobytes())
    return {'sha256': digest.hexdigest(), 'genes': len(genes), 'barcodes': len(cells),
            'nonzero_entries': int(canonical.nnz), 'summary': sc_summary(item)}


def run_single_cell(workspace, tools, case, reference, index, spec, status):
    root = workspace / 'single-cell'
    indexes = root / 'indexes'
    indexes.mkdir(parents=True, exist_ok=True)
    link = indexes / 'kallistors.idx'
    if not link.exists():
        link.symlink_to(index.resolve())
    if fingerprint(link)['sha256'] != fingerprint(index)['sha256']:
        raise ValueError('Single-cell comparison must use identical index bytes')
    # Verify every prepared lane and the common on-list before any timing.
    prep = case['provenance']
    for lane in prep['lanes']:
        for name in ('read1', 'read2'):
            item = lane[name]
            if fingerprint(item['path']) != item:
                raise ValueError('Single-cell prepared input changed')
    if fingerprint(case['onlist']) != prep['onlist']:
        raise ValueError('Single-cell barcode on-list changed')
    classes, baseline_bus, baseline_matrix, baseline_output = {}, None, None, None
    proofs = []
    # Revalidate previously completed outputs when resuming. Canonical intermediate
    # digests remain tied to immutable run identities; raw BUS is removed only after
    # its full weighted comparison and matrix validation have passed.
    for threads in spec['design']['single_cell_workers']:
        for phase, rep in [('warmup', 0), ('measured', 1), ('measured', 2), ('measured', 3)]:
            order = ('kallisto', 'kallistors') if rep % 2 == 0 else ('kallistors', 'kallisto')
            for tool in order:
                key = f'pbmc_1k_v3/t{threads}/{phase}-{rep}/{tool}'
                dest = root / key
                adapter_tools = {'kallistors': tools[tool], 'bustools': tools['bustools']}
                commands = sc_plan(adapter_tools, 'kallistors-bustools', case, reference, indexes, dest, threads)
                status(f'single-cell {phase} {rep}: {tool}, {threads} workers')
                row = run_one(root, key, case['id'], tool, threads, phase, rep, commands, spec)
                item = load_sc(dest, 'kallistors-bustools')
                matrix = matrix_identity(item)
                if baseline_matrix is None:
                    baseline_matrix, baseline_output = matrix['sha256'], dest
                if matrix['sha256'] != baseline_matrix:
                    write_json(root / 'matrix-parity-failure.json', sc_agreement(
                        load_sc(baseline_output, 'kallistors-bustools'), item))
                    raise ValueError('Complete canonical single-cell matrices differ')
                if row.get('bus_parity'):
                    bus = row['bus_parity']
                    if baseline_bus is None:
                        saved = root / 'parity' / f't{threads}-{phase}-{rep}-{tool}.bus'
                        if fingerprint(saved)['sha256'] != bus['sha256']:
                            raise ValueError('Retained canonical BUS evidence changed')
                    # Restore the canonical class numbering used by saved digests.
                    saved_classes = root / 'canonical-classes.json'
                    classes = {tuple(v): i for i, v in enumerate(json.loads(saved_classes.read_text()))}
                else:
                    canonical = root / 'parity' / (f't{threads}-{phase}-{rep}-{tool}.bus')
                    canonical.parent.mkdir(exist_ok=True)
                    transcripts = (dest / 'mapping/transcripts.txt').read_text().splitlines()
                    bus = canonical_bus(dest / 'mapping', canonical, classes, transcripts, tools['bustools'], threads)
                    # Save the exact transcript-set numbering before retiring raw products.
                    write_json(root / 'canonical-classes.json', [list(v) for v in classes])
                    if baseline_bus is not None and bus['sha256'] != baseline_bus:
                        write_json(root / 'bus-parity-failure.json', {'baseline': baseline_bus, 'candidate': bus})
                        raise ValueError('Full weighted barcode/UMI/transcript-set/flag/count parity failed')
                    if baseline_bus is not None:
                        canonical.unlink()
                baseline_bus = baseline_bus or bus['sha256']
                if bus['sha256'] != baseline_bus or bus['run_info']['n_processed'] != case['pairs']:
                    raise ValueError('Saved BUS evidence or native processed count failed')
                row.update(bus_parity=bus, matrix_identity=matrix, matrix_summary=matrix['summary'])
                write_json(dest / 'result.json', row)
                if 'retired_intermediates' not in row:
                    row['retired_intermediates'] = retire_intermediates(dest)
                    write_json(dest / 'result.json', row)
                proofs.append({'threads': threads, 'phase': phase, 'repetition': rep, 'tool': tool,
                               'bus_sha256': bus['sha256'], 'matrix_sha256': matrix['sha256'],
                               'processed': bus['run_info']['n_processed'], 'exact': True})
                del item
    write_json(root / 'parity.json', {'complete': True, 'all_full_libraries_exact': True,
        'matrix_comparison': 'Sorted gene/barcode labels and full sparse integer values, independent of entry order',
        'bus_comparison': 'Complete weighted barcode/UMI/counted-transcript-set/flag tuples; header text and padding normalized',
        'trials': proofs})
