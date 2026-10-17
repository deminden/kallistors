#!/usr/bin/env python3
"""Prepare full EC/fragment-model replays from diagnostic ECs and retained HDF5.

Requires native h5py and NumPy. Target filtering affects the metadata report only;
every transcript and every supplied nonzero EC is retained in replay inputs.
"""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ecs', type=Path, required=True)
    parser.add_argument('--original', type=Path, required=True, help='original abundance.h5')
    parser.add_argument('--candidate', type=Path, required=True, help='candidate abundance.h5')
    parser.add_argument('--target-contains', default='')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        parser.error('--out must be new or empty')
    args.out.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.original) as original, h5py.File(args.candidate) as candidate:
        names = original['aux/ids'].asstr()[:]
        np.testing.assert_array_equal(names, candidate['aux/ids'].asstr()[:])
        lengths = original['aux/lengths'][()]
        np.testing.assert_array_equal(lengths, candidate['aux/lengths'][()])
        selected = [i for i, name in enumerate(names)
                    if args.target_contains and args.target_contains in name]
        (args.out / 'lengths').write_text(''.join(f'{int(x)}\n' for x in lengths))
        metadata = {'targets': len(names), 'target_contains': args.target_contains,
                    'selected_targets': [{'id': i, 'name': str(names[i]), 'length': int(lengths[i])}
                                         for i in selected], 'models': {}}
        for label, handle, path in [('original', original, args.original),
                                    ('candidate', candidate, args.candidate)]:
            histogram = handle['aux/fld'][()]
            (args.out / f'{label}.histogram').write_text(''.join(f'{int(x)}\n' for x in histogram))
            metadata['models'][label] = {
                'source_file': path.name, 'source_sha256': digest(path),
                'fragment_observations': int(histogram.sum()),
                'fragment_mean': float(np.dot(histogram, np.arange(len(histogram))) / histogram.sum()),
                'selected_counts': [float(handle['est_counts'][i]) for i in selected],
                'selected_eff_lengths': [float(handle['aux/eff_lengths'][i]) for i in selected]}
    count = 0
    with args.ecs.open() as source, (args.out / 'ecs.tsv').open('w') as output:
        for line in source:
            if line.startswith('#') or not line.strip():
                continue
            _, observations, targets = line.rstrip('\n').split('\t')
            output.write(f'{count}\t{observations}\t{targets}\n')
            count += 1
    metadata['ecs'] = count
    metadata['source_ecs_sha256'] = digest(args.ecs)
    metadata['replay_ecs_sha256'] = digest(args.out / 'ecs.tsv')
    (args.out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
