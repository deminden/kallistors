"""One target universe and identifier policy, shared by every transcriptome index."""

import gzip
from itertools import chain
from pathlib import Path

from benchmark_common import digest, fingerprint, write_json
from .environment import read_json


def fasta_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as handle:
        header, parts = None, []
        for line in handle:
            if line.startswith(b'>'):
                if header is not None:
                    yield header, b''.join(parts).upper()
                header, parts = line[1:].strip(), []
            else:
                if header is None and line.strip():
                    raise ValueError('Sequence before a FASTA header')
                parts.append(line.strip())
        if header is not None:
            yield header, b''.join(parts).upper()


def canonical_reference(source, out, kind='gencode50'):
    sources = list(source) if isinstance(source, (list, tuple)) else [source]
    out.mkdir(parents=True, exist_ok=True)
    record = out / 'reference.json'
    if record.exists():
        saved = read_json(record)
        expected = saved['source'].get('components', [saved['source']])
        if [item['sha256'] for item in expected] != [digest(p) for p in sources]:
            raise ValueError('Reference source changed')
        for key in ['fasta', 'transcripts', 'gene_map']:
            if saved.get(key) and digest(out / Path(saved[key]['path']).name) != saved[key]['sha256']:
                raise ValueError(f'Prepared reference changed: {key}')
        for key in ['fasta', 'transcripts', 'gene_map']:
            if saved.get(key):
                saved[key]['path'] = str((out / Path(saved[key]['path']).name).resolve())
        write_json(record, saved)
        return saved
    fasta, t2g, table = out / 'transcripts.fa', out / 'transcripts-to-genes.tsv', out / 'transcripts.tsv'
    seen, gene_versions = set(), {}
    changed, targets, bases, short = 0, 0, 0, 0
    state_counts = {'mature': 0, 'nascent': 0}
    with fasta.open('wb') as dst, t2g.open('w') as genes, table.open('w') as meta:
        meta.write('transcript_id\tgene_id\tlength\tbiotype\n')
        for header, seq in chain.from_iterable(fasta_records(p) for p in sources):
            fields = header.split()[0].decode('ascii').split('|')
            tx = fields[0]
            if not tx or tx in seen or not seq:
                raise ValueError(f'Empty/duplicate target or empty sequence: {tx}')
            seen.add(tx)
            if kind == 'gencode50':
                if len(fields) < 8:
                    raise ValueError('Expected full GENCODE headers')
                original_gene = fields[1]
                gene = original_gene.rsplit('.', 1)[0]
                previous = gene_versions.setdefault(gene, original_gene)
                if previous != original_gene:
                    raise ValueError(f'Multiple versions of one gene: {gene}')
                biotype = fields[7]
                genes.write(f'{tx}\t{gene}\n')
            elif kind == 'official-nac':
                annotations = dict(word.split(':', 1) for word in header.decode('ascii').split()[1:] if ':' in word)
                original_gene = annotations['gene_id']
                gene = original_gene.rsplit('.', 1)[0]
                previous = gene_versions.setdefault(gene, original_gene)
                if previous != original_gene:
                    raise ValueError(f'Multiple versions of one gene: {gene}')
                biotype = 'nascent' if tx == original_gene else 'mature'
                state_counts[biotype] += 1
                genes.write(f'{tx}\t{gene}\n')
            else:
                raise ValueError(f'Unsupported reference kind: {kind}')
            if any(c not in b'ACGT' for c in seq):
                # Feed every indexer the same deterministic replacement instead
                # of letting each tool apply its own ambiguous-base policy.
                seq = bytes(c if c in b'ACGT' else b'ACGT'[i % 4] for i, c in enumerate(seq))
                changed += 1
            dst.write(b'>' + tx.encode() + b'\n')
            for start in range(0, len(seq), 80):
                dst.write(seq[start:start + 80] + b'\n')
            meta.write(f'{tx}\t{gene}\t{len(seq)}\t{biotype}\n')
            targets += 1
            bases += len(seq)
            short += len(seq) < 31
    source_record = fingerprint(sources[0]) if len(sources) == 1 else {'components': [fingerprint(p) for p in sources]}
    result = {'kind': kind, 'source': source_record, 'targets': targets, 'genes': len(gene_versions),
              'bases': bases, 'targets_shorter_than_k': short, 'ambiguous_sequences_normalized': changed,
              'normalization': 'Uppercase; ID before first pipe, retaining transcript version; every non-ACGT base replaced by ACGT[position mod 4]; no sequence deduplication, polyA clipping or biotype filtering.',
              'fasta': fingerprint(fasta), 'transcripts': fingerprint(table),
              'gene_map': fingerprint(t2g),
              'k': 31, 'minimizers': {'kallistors': 'native default', 'salmon': 19, 'piscem': 19}}
    if kind == 'official-nac':
        result.update(mature_targets=state_counts['mature'], nascent_targets=state_counts['nascent'])
        if (result['mature_targets'], result['nascent_targets']) != (227665, 39546):
            raise ValueError('Official NAC target counts differ from the pinned assets')
    write_json(record, result)
    return result


def official_nac_sources(data):
    """Verify publisher FASTAs directly; never reconstruct intronic sequences."""
    from benchmark_common import fetch
    base = 'https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/'
    assets = [
        ('human_nac.mature.fa.gz', 58045671, '271969e349a5d80657dd94a313e4efc8aea02fdb55193f5e71e9460f73d82338'),
        ('human_nac.nascent.fa.gz', 558305249, '8cdd770313fc8592e4fc17ec207824361a30bb943e1a7844e01d58a39c26450f'),
    ]
    directory = data / 'official-human-nac'
    for name, size, sha in assets:
        fetch(base + name, directory / name, size=size, sha256=sha)
    return [directory / name for name, _, _ in assets]
