"""Reuse checksum-verified prepared libraries; prefixes are validation fixtures only."""
import gzip
import math
from contextlib import ExitStack
from pathlib import Path
from benchmark_common import fingerprint, write_json
from .environment import CONFIG, read_json


def verified_file(record, directory):
    path = directory / Path(record['path']).name
    found = fingerprint(path)
    if any(found[k] != record[k] for k in ('bytes', 'sha256')):
        raise ValueError(f'Input checksum mismatch: {path}')
    return str(path.resolve())


def bulk(data, accession):
    directory = data / accession.lower() / 'trimmed-fastp'
    meta = read_json(directory / 'preparation.json')
    if meta['trimmer']['name'] != 'fastp' or meta['trimmer']['version'] != '1.4.0':
        raise ValueError('Inputs require the pinned fastp preparation')
    model_file = CONFIG / ('nascent-library-models.json' if accession == 'SRR10919051' else 'library-models.json')
    library = read_json(model_file)['libraries'][accession]
    if [v['sha256'] for v in meta['outputs']] != [v['sha256'] for v in library['prepared_inputs']]:
        raise ValueError('Fragment calibration requires the identical prepared paired inputs')
    hist = library['calibration']['histogram_counts']
    n = sum(hist)
    if n <= 0 or any(v < 0 for v in hist):
        raise ValueError('Invalid calibration histogram')
    mean = math.fsum(i * count for i, count in enumerate(hist)) / n
    sd = math.sqrt(math.fsum(count * (i - mean)**2 for i, count in enumerate(hist)) / n)
    fragment_model = {'mean': mean, 'sd': sd, 'source': 'frozen paired-end empirical calibration',
                      'samples': n, 'calibration_specification': fingerprint(model_file),
                      'source_hdf5_sha256': library['calibration']['source_hdf5_sha256']}
    return {'id': accession, 'kind': 'bulk', 'pairs': meta['pairs_after'],
            'strand': library['metadata']['strand'], 'library_metadata': library['metadata'],
            'single_end_fragment_model': fragment_model,
            'read1': [verified_file(meta['outputs'][0], directory)],
            'read2': [verified_file(meta['outputs'][1], directory)], 'provenance': meta}


def single_cell(data):
    directory = data / 'tenx-pbmc-1k-v3'
    meta = read_json(directory / 'preparation.json')
    if not meta['retained_barcode_umi_sequences_unchanged']:
        raise ValueError('Barcode and UMI sequences must be retained')
    return {'id': 'pbmc_1k_v3', 'kind': 'single_cell', 'pairs': meta['total_pairs'],
            'read1': [verified_file(lane['read1'], directory / 'trimmed-fastp') for lane in meta['lanes']],
            'read2': [verified_file(lane['read2'], directory / 'trimmed-fastp') for lane in meta['lanes']],
            'onlist': verified_file(meta['onlist'], directory), 'provenance': meta}


def single_end(case):
    """Use R1 from the same prepared library, with a shared declared model."""
    result = {**case, 'id': case['id'] + '-se', 'layout': 'single',
              'source_case': case['id'], 'reads': case['pairs'], 'read2': [],
              'fragment_model': {'mean': 250, 'sd': 30, 'source': 'known simulation model'}
                                if case.get('truth') else case['single_end_fragment_model']}
    result.pop('pairs')
    return result


def prefix_case(case, out, records=100000):
    """Keep matched prefixes from every lane, including derived single-end inputs."""
    out.mkdir(parents=True, exist_ok=True)
    record = out / 'inputs.json'
    if record.exists():
        saved = read_json(record)
        if saved['source'] != case or saved['limit'] != records:
            raise ValueError('Pilot input specification changed')
        for item in saved['files']:
            verified_file(item, out)
        return saved['case']
    result = {**case, 'id': case['id'] + '-pilot', 'read1': [], 'read2': []}
    single = case.get('layout') == 'single'
    lanes = [(r1,) for r1 in case['read1']] if single else list(zip(case['read1'], case['read2'], strict=True))
    total, files = 0, []
    for lane, inputs in enumerate(lanes):
        paths = [out / f'lane{lane}-R{i+1}.fq.gz' for i in range(len(inputs))]
        with ExitStack() as stack:
            sources = [stack.enter_context(gzip.open(path, 'rb')) for path in inputs]
            outputs = [stack.enter_context(gzip.open(path, 'wb', compresslevel=1)) for path in paths]
            for _ in range(max(1, records // len(lanes))):
                reads = [[src.readline() for _ in range(4)] for src in sources]
                if all(not rec[0] for rec in reads):
                    break
                if any(not line for rec in reads for line in rec):
                    raise ValueError('Truncated or unmatched FASTQ')
                names = [rec[0].split()[0].removesuffix(b'/1').removesuffix(b'/2') for rec in reads]
                if len(names) == 2 and names[0] != names[1]:
                    raise ValueError('FASTQ pairs do not match')
                for dst, rec in zip(outputs, reads, strict=True):
                    dst.writelines(rec)
                total += 1
        result['read1'].append(str(paths[0].resolve()))
        if not single:
            result['read2'].append(str(paths[1].resolve()))
        files += [fingerprint(p) for p in paths]
    result['reads' if single else 'pairs'] = total
    write_json(record, {'source': case, 'limit': records, 'case': result, 'files': files})
    return result
