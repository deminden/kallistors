#!/usr/bin/env python3
"""Prepare GENCODE v50 cDNA or Pachter Lab human NAC references."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_common import fetch, fingerprint, write_json

BASE = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/'
TRANSCRIPTS_MD5 = '4d5bd41b247480e12c0b370490f52ddd'


def read_fasta(path, destination=None, seen=None):
    seen = set() if seen is None else seen
    targets = bases = current_bases = 0
    line = b''
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as source:
        for line in source:
            if line.startswith(b'>'):
                fields = line[1:].split()
                if not fields or fields[0] in seen or (targets and not current_bases):
                    raise ValueError(f'Missing/duplicate FASTA ID or empty sequence: {path}')
                seen.add(fields[0])
                targets += 1
                current_bases = 0
            else:
                length = len(line.strip())
                if length and not targets:
                    raise ValueError(f'Sequence before the first FASTA header: {path}')
                bases += length
                current_bases += length
            if destination is not None:
                destination.write(line)
    if not targets or not current_bases:
        raise ValueError(f'Empty FASTA: {path}')
    if destination is not None and not line.endswith(b'\n'):
        destination.write(b'\n')
    return {'targets': targets, 'bases': bases}


def official_nac(out):
    base = 'https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/'
    assets = {'human_nac.mature.fa.gz': 58045671, 'human_nac.nascent.fa.gz': 558305249}
    sources = [fetch(base + name, out / name, size=size) for name, size in assets.items()]
    fasta = out / 'human_nac.mature-plus-nascent.fa'
    if fasta.exists():
        raise ValueError(f'Output already exists: {fasta}')
    seen, components = set(), []
    with fasta.open('wb') as destination:
        for name in assets:
            counts = read_fasta(out / name, destination, seen)
            destination.write(b'\n')
            components.append({'file': name, **counts})
    return {'label': 'Pachter Lab human NAC (Ensembl 108)', 'annotation': 'Ensembl 108, GRCh38 primary assembly',
            'definition': 'Pachter Lab mature cDNAs followed by Pachter Lab nascent transcripts; IDs and sequences unchanged',
            'sources': sources, 'components': components,
            'mature_targets': components[0]['targets'], 'nascent_targets': components[1]['targets'],
            'total_targets': len(seen), 'fasta': fingerprint(fasta)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', choices=('gencode50', 'official-nac'), required=True)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--transcripts', type=Path, help='existing GENCODE v50 ALL FASTA.gz for gencode50')
    args = parser.parse_args()
    try:
        roots = {'gencode50': 'data/srr-small',
                 'official-nac': 'data/official-human-nac'}
        out = (args.out or Path(roots[args.reference])).resolve()
        out.mkdir(parents=True, exist_ok=True)
        report = out / 'reference.json'
        if report.exists():
            raise ValueError(f'Reference already prepared: {report}')
        if args.reference == 'official-nac':
            if args.transcripts:
                raise ValueError('--transcripts applies only to gencode50')
            result = official_nac(out)
        else:
            transcripts = (args.transcripts or out / 'gencode.v50.transcripts.fa.gz').resolve()
            source = fetch(BASE + 'gencode.v50.transcripts.fa.gz', transcripts, TRANSCRIPTS_MD5)
            counts = read_fasta(transcripts)
            result = {'label': 'GENCODE v50 ALL transcripts', 'mature_targets': counts['targets'],
                      'nascent_targets': 0, 'total_targets': counts['targets'], 'sources': [source],
                      'fasta': fingerprint(transcripts)}
        write_json(report, {'schema_version': 1, 'reference': args.reference, **result})
        print(json.dumps(result, indent=2))
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f'Reference preparation failed: {error}\n')


if __name__ == '__main__':
    main()
