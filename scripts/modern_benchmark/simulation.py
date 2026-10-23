"""Independent, seeded fragment-origin control; not a model of all biological biases."""
import gzip
from pathlib import Path
import numpy as np
from benchmark_common import fingerprint, write_json
from .environment import read_json
from .reference import fasta_records


def simulate(reference, out, pairs=4000000, seed=48321):
    record = out / 'simulation.json'
    if record.exists():
        saved = read_json(record)
        if saved['reference_sha256'] != reference['fasta']['sha256'] or saved['pairs'] != pairs or saved['seed'] != seed:
            raise ValueError('Simulation definition changed')
        for f in saved['files']:
            if fingerprint(out / Path(f['path']).name)['sha256'] != f['sha256']:
                raise ValueError('Simulation input or truth changed')
        return saved['case']
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.Generator(np.random.PCG64(seed))
    seqs = [s for _, s in fasta_records(reference['fasta']['path'])]
    lengths = np.array([len(s) for s in seqs])
    eligible = np.flatnonzero(lengths >= 400)
    selected = np.sort(rng.choice(eligible, min(50000, len(eligible)), replace=False))
    abundance = rng.lognormal(0, 2, len(selected))
    # Sampling fragments in proportion to molecule abundance times effective length.
    mass = abundance * (lengths[selected] - 249)
    mass /= mass.sum()
    truth = np.zeros(len(seqs), dtype=np.int64)
    tpm = np.zeros(len(seqs), dtype=float)
    tpm[selected] = abundance / abundance.sum() * 1e6
    paths = [out / 'reads_R1.fastq.gz', out / 'reads_R2.fastq.gz']
    complement = bytes.maketrans(b'ACGT', b'TGCA')
    bases = np.frombuffer(b'ACGT', dtype=np.uint8)
    with gzip.open(paths[0], 'wb', compresslevel=1) as a, gzip.open(paths[1], 'wb', compresslevel=1) as b:
        for offset in range(0, pairs, 10000):
            count = min(10000, pairs - offset)
            targets = rng.choice(selected, count, p=mass)
            fragment_lengths = np.clip(np.rint(rng.normal(250, 30, count)).astype(int), 200, 400)
            starts = (rng.random(count) * (lengths[targets] - fragment_lengths + 1)).astype(int)
            reverse = rng.integers(0, 2, count)
            np.add.at(truth, targets, 1)
            reads = [[], []]
            for j, (tx, length, start, rev) in enumerate(zip(targets, fragment_lengths, starts, reverse)):
                fragment = seqs[tx][start:start + length]
                if rev:
                    fragment = fragment.translate(complement)[::-1]
                reads[0].append(fragment[:100])
                reads[1].append(fragment[-100:].translate(complement)[::-1])
            for mate, dst in enumerate((a, b)):
                array = np.frombuffer(b''.join(reads[mate]), dtype=np.uint8).copy().reshape(count, 100)
                errors = rng.random(array.shape) < 0.001
                # Uniform among the other three bases; no no-op substitutions.
                shifts = rng.integers(1, 4, int(errors.sum()))
                old = np.searchsorted(bases, array[errors])
                array[errors] = bases[(old + shifts) % 4]
                lines = []
                for j, seq in enumerate(array):
                    lines.append(f'@truth{offset+j}/{mate+1}\n'.encode() + seq.tobytes() + b'\n+\n' + b'?' * 100 + b'\n')
                dst.write(b''.join(lines))
    np.savez_compressed(out / 'truth.npz', fragment_counts=truth, molecule_tpm=tpm)
    case = {'id': 'fragment-truth', 'kind': 'bulk', 'pairs': pairs,
            'read1': [str(paths[0].resolve())], 'read2': [str(paths[1].resolve())],
            'truth': str((out / 'truth.npz').resolve()), 'preprocessing': 'Adapter-free simulated Q30 reads; no trimming required.'}
    write_json(record, {'seed': seed, 'rng': 'PCG64', 'pairs': pairs,
        'reference_sha256': reference['fasta']['sha256'], 'case': case,
        'model': '50,000 eligible transcripts chosen without replacement; length >=400; lognormal(0,2) molecule abundances; fragment probabilities proportional to abundance times (length-249); fragment normal(250,30) clipped 200..400; uniform start and random orientation; PE100; independent substitution probability 0.001, Q30. No indels, contamination, polymorphism, sequence bias or PCR duplication. Truth is independent of all tested tools. This simplified control cannot establish a universal accuracy ranking.',
        'files': [fingerprint(p) for p in [*paths, out / 'truth.npz']]})
    return case
