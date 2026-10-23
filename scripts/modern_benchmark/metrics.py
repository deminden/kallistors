"""Canonical labels, full-universe agreement and independent-truth accuracy."""
import csv
from pathlib import Path
import numpy as np
from scipy import io, sparse, stats

QUANTILES = (0.99, 0.999, 0.9999)


def differences(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    delta = np.abs(a - b)
    return {'features': len(a), 'sum_a': float(a.sum()), 'sum_b': float(b.sum()),
            'mae': float(delta.mean()), 'max_abs': float(delta.max(initial=0)),
            **{name: float(np.quantile(delta, q)) for name, q in zip(('p99', 'p99_9', 'p99_99'), QUANTILES)},
            'pearson_log1p': correlation(np.log1p(a), np.log1p(b)),
            'spearman': correlation(stats.rankdata(a), stats.rankdata(b))}


def correlation(a, b):
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def load_bulk(out, method, reference):
    out = Path(out)
    file = out / ('quant/abundance.tsv' if method == 'kallistors' else
                  'quant/quant.sf' if method.startswith('salmon') else 'quant.quant')
    with file.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    names = [r.get('target_id', r.get('Name', r.get('target_name'))) for r in rows]
    if len(set(names)) != len(names) or None in names:
        raise ValueError(f'Duplicate or missing transcript labels: {file}')
    with Path(reference['transcripts']['path']).open() as handle:
        universe = list(csv.DictReader(handle, delimiter='\t'))
    if set(names) != {r['transcript_id'] for r in universe}:
        raise ValueError(f'Target universe mismatch: {file}')
    lookup = {name: i for i, name in enumerate(names)}
    order = [lookup[r['transcript_id']] for r in universe]
    counts = np.array([float(r.get('est_counts', r.get('NumReads', r.get('ecount')))) for r in rows])[order]
    tpm = np.array([float(r.get('tpm', r.get('TPM'))) for r in rows])[order]
    if not np.isfinite(counts).all() or not np.isfinite(tpm).all() or (counts < 0).any() or (tpm < 0).any():
        raise ValueError(f'Invalid abundance: {file}')
    genes = sorted({r['gene_id'] for r in universe})
    gene_index = {name: i for i, name in enumerate(genes)}
    tx_to_gene = np.array([gene_index[r['gene_id']] for r in universe])
    result = {'transcript_counts': counts, 'transcript_tpm': tpm,
              'gene_counts': np.bincount(tx_to_gene, weights=counts, minlength=len(genes)),
              'gene_tpm': np.bincount(tx_to_gene, weights=tpm, minlength=len(genes))}
    if reference['kind'] == 'official-nac':
        states = sorted({(r['gene_id'], r['biotype']) for r in universe})
        state_index = {label: i for i, label in enumerate(states)}
        tx_to_state = np.array([state_index[(r['gene_id'], r['biotype'])] for r in universe])
        result.update(gene_state_counts=np.bincount(tx_to_state, weights=counts, minlength=len(states)),
                      gene_state_tpm=np.bincount(tx_to_state, weights=tpm, minlength=len(states)))
        for state in ('mature', 'nascent'):
            mask = np.array([r['biotype'] == state for r in universe])
            result[state + '_counts'] = counts[mask]
            result[state + '_tpm'] = tpm[mask]
    return result


def labels(path, barcodes=False):
    values = [line.split('\t')[0].strip() for line in Path(path).read_text().splitlines()]
    if barcodes:
        values = [v.removesuffix('-1') for v in values]
    if len(set(values)) != len(values):
        raise ValueError(f'Label collision: {path}')
    return values


def load_sc(out, method):
    out = Path(out)
    if method == 'kallistors-bustools':
        matrix = io.mmread(out / 'quant.mtx').tocsr()
        genes, cells = labels(out / 'quant.genes.txt'), labels(out / 'quant.barcodes.txt', True)
        # bustools is cells x genes, like alevin-fry; assert dimensions below.
        matrix = matrix.T.tocsr()
    else:
        q = out / 'quant/alevin'
        if not q.exists():
            q = out / 'quant'
        matrix = io.mmread(q / 'quants_mat.mtx').T.tocsr()
        genes, cells = labels(q / 'quants_mat_cols.txt'), labels(q / 'quants_mat_rows.txt', True)
    if matrix.shape != (len(genes), len(cells)):
        raise ValueError(f'Matrix/label shape mismatch in {out}: {matrix.shape}')
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError(f'Invalid count matrix: {out}')
    return genes, cells, matrix


def sc_summary(item):
    genes, cells, matrix = item
    umis = np.asarray(matrix.sum(axis=0)).ravel()
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    called = (umis >= 500) & (detected >= 200)
    return {'raw_barcodes': len(cells), 'nonzero_genes': int((matrix.sum(axis=1) > 0).sum()),
            'total_umis': float(umis.sum()), 'common_gate_cells': int(called.sum()),
            'median_umis_called': float(np.median(umis[called])) if called.any() else None,
            'median_genes_called': float(np.median(detected[called])) if called.any() else None}


def align_sc(item, genes, cells):
    oldgenes, oldcells, matrix = item
    gi, ci = {v: i for i, v in enumerate(genes)}, {v: i for i, v in enumerate(cells)}
    coo = matrix.tocoo()
    return sparse.coo_matrix((coo.data, ([gi[oldgenes[i]] for i in coo.row],
                                       [ci[oldcells[i]] for i in coo.col])), shape=(len(genes), len(cells))).tocsr()


def sparse_quantile(values, total, q):
    values = np.sort(values)
    zeros = total - len(values)
    pos = q * (total - 1)
    lo, hi = int(np.floor(pos)), int(np.ceil(pos))
    def at(i):
        return 0.0 if i < zeros else float(values[i - zeros])
    return at(lo) + (pos - lo) * (at(hi) - at(lo))


def sc_agreement(a, b):
    genes, cells = sorted(set(a[0]) | set(b[0])), sorted(set(a[1]) | set(b[1]))
    x, y = align_sc(a, genes, cells), align_sc(b, genes, cells)
    diff = abs(x - y)
    diff.eliminate_zeros()
    union = x.maximum(y).nnz
    total = x.shape[0] * x.shape[1]
    result = {'genes': len(genes), 'barcodes_union': len(cells), 'max_abs': float(diff.max()),
              'mae_all_entries': float(diff.sum() / total) if total else 0,
              **{f'{n}_all_entries': sparse_quantile(diff.data, total, q) for n, q in
                 zip(('p99', 'p99_9', 'p99_99'), QUANTILES)},
              **{f'{n}_nonzero_union': sparse_quantile(diff.data, union, q) if union else 0 for n, q in
                 zip(('p99', 'p99_9', 'p99_99'), QUANTILES)}}
    masks = []
    for m in (x, y):
        u = np.asarray(m.sum(axis=0)).ravel()
        g = np.asarray((m > 0).sum(axis=0)).ravel()
        masks.append((u >= 500) & (g >= 200))
    intersection, cell_union = masks[0] & masks[1], masks[0] | masks[1]
    result.update(common_called_cells=int(intersection.sum()),
                  called_cell_jaccard=float(intersection.sum() / cell_union.sum()) if cell_union.any() else None,
                  per_barcode_umis=differences(np.asarray(x.sum(axis=0)).ravel(), np.asarray(y.sum(axis=0)).ravel()),
                  per_gene_umis=differences(np.asarray(x.sum(axis=1)).ravel(), np.asarray(y.sum(axis=1)).ravel()),
                  per_barcode_detected_genes=differences(np.asarray((x > 0).sum(axis=0)).ravel(),
                                                        np.asarray((y > 0).sum(axis=0)).ravel()))
    if intersection.any():
        a_called, b_called = x[:, intersection], y[:, intersection]
        delta = abs(a_called - b_called)
        delta.eliminate_zeros()
        called_union = a_called.maximum(b_called).nnz
        result.update(common_called_matrix_max_abs=float(delta.max()),
                      **{f'{n}_common_called_nonzero_union': sparse_quantile(delta.data, called_union, q)
                         for n, q in zip(('p99', 'p99_9', 'p99_99'), QUANTILES)},
                      common_called_per_barcode_umis=differences(np.asarray(a_called.sum(axis=0)).ravel(),
                                                                 np.asarray(b_called.sum(axis=0)).ravel()),
                      common_called_per_barcode_detected_genes=differences(
                          np.asarray((a_called > 0).sum(axis=0)).ravel(),
                          np.asarray((b_called > 0).sum(axis=0)).ravel()))
    return result
