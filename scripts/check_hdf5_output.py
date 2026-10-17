#!/usr/bin/env python3
"""Validate quantification with native HDF5 and summarize saved fragment histograms."""

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np

from benchmark_common import fingerprint, write_json


def check_output(folder, significant_digits=17):
    with (folder / 'abundance.tsv').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'target_id', 'length', 'eff_length', 'est_counts', 'tpm'}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError('Missing abundance columns')
        rows = list(reader)
    if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('Empty or malformed abundance table')
    tpm = np.array([float(row['tpm']) for row in rows])
    if not np.isfinite(tpm).all() or (tpm < 0).any():
        raise ValueError('Invalid TSV TPM values')
    names = [row['target_id'] for row in rows]
    if len(set(names)) != len(names) or not all(names):
        raise ValueError('Missing or duplicate transcript IDs')
    info = json.loads((folder / 'run_info.json').read_text())
    datasets, fragment = [], None
    exact = True
    with h5py.File(folder / 'abundance.h5', 'r') as handle:
        for name, column in [('est_counts', 'est_counts'), ('aux/eff_lengths', 'eff_length'), ('aux/lengths', 'length')]:
            actual = handle[name][()]
            expected = np.array([float(row[column]) for row in rows])
            if actual.shape != (len(rows),) or not np.isfinite(expected).all() or (expected < 0).any():
                raise ValueError(f'Invalid TSV values or HDF5 dimensions: {name}')
            if not np.isfinite(actual).all() or (actual < 0).any():
                raise ValueError(f'Invalid HDF5 values: {name}')
            exact &= bool(np.array_equal(actual, expected))
            if name == 'aux/lengths' or significant_digits == 17:
                np.testing.assert_array_equal(actual, expected, err_msg=name)
            else:
                # Original TSV uses six significant digits, while HDF5 keeps f64.
                rounded = np.array([float(format(float(v), f'.{significant_digits}g')) for v in actual])
                np.testing.assert_array_equal(rounded, expected, err_msg=name)
        if list(handle['aux/ids'].asstr()[()]) != names:
            raise ValueError('HDF5 and TSV transcript identities differ')
        for name, field in [('num_processed', 'n_processed'), ('num_bootstrap', 'n_bootstraps')]:
            if int(handle['aux/' + name][0]) != info[field]:
                raise ValueError(f'Metadata mismatch: {field}')
        if 'n_targets' in info and info['n_targets'] != len(rows):
            raise ValueError('Target count metadata mismatch')
        bootstrap_names = []

        def visit(name, dataset):
            if not isinstance(dataset, h5py.Dataset):
                return
            values = dataset[()]
            if name.startswith('bootstrap/bs'):
                bootstrap_names.append(name)
                if values.shape != (len(rows),) or not np.isfinite(values).all() or (values < 0).any():
                    raise ValueError(f'Invalid bootstrap values: {name}')
            datasets.append({'name': name, 'shape': dataset.shape, 'chunks': dataset.chunks,
                             'compression': dataset.compression, 'compression_level': dataset.compression_opts})
        handle.visititems(visit)
        if set(bootstrap_names) != {f'bootstrap/bs{i}' for i in range(info['n_bootstraps'])}:
            raise ValueError('Missing or extra bootstrap datasets')
        if 'aux/fld' in handle:
            histogram = np.asarray(handle['aux/fld'][()], dtype=float)
            if histogram.ndim != 1 or not np.isfinite(histogram).all() or (histogram < 0).any():
                raise ValueError('Invalid fragment histogram')
            total = float(histogram.sum())
            mean = float(np.dot(histogram, np.arange(len(histogram))) / total) if total else None
            sd = float(np.sqrt(np.dot(histogram, (np.arange(len(histogram)) - mean) ** 2) / total)) if total else None
            fragment = {'fragment_observations': total, 'mean': mean, 'sd': sd}
    return {'native_hdf5': h5py.version.hdf5_version, 'h5py': h5py.__version__, 'targets': len(rows),
            'tsv_values_exact': exact, 'comparison_significant_digits': significant_digits,
            'fragment_model': fragment, 'datasets': datasets,
            'files': {name: fingerprint(folder / name) for name in ('abundance.tsv', 'abundance.h5', 'run_info.json')}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--json', type=Path)
    parser.add_argument('--significant-digits', type=int, choices=range(1, 18), default=17)
    args = parser.parse_args()
    try:
        result = check_output(args.output, args.significant_digits)
        if args.json:
            write_json(args.json, result)
        if args.json:
            print(f'Validated {result["targets"]} transcripts with native HDF5 {result["native_hdf5"]}')
        else:
            print(json.dumps(result, indent=2))
    except (OSError, ValueError, KeyError, AssertionError) as error:
        parser.exit(1, f'HDF5 validation failed: {error}\n')
