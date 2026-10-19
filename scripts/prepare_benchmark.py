#!/usr/bin/env python3
"""Retrieve checksum-pinned RNA-seq libraries and prepare paired reads with fastp."""

import argparse
import concurrent.futures
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

from benchmark_common import DATASETS, datasets, digest, fetch, fingerprint, write_json

FASTP_VERSION = '1.4.0'
TENX_PAGE = 'https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0'
TENX_ARCHIVE = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar'
TENX_SHA256 = '873f6511a7cb732c0a709d2fb9201ca4f1230a70f3a04e46c0198a9e4be896fd'
TENX_ONLIST = 'https://github.com/BUStools/getting_started/releases/download/species_mixing/10xv3_whitelist.txt'
TENX_ONLIST_SHA256 = '843a6f7038db8cb3c06f3dc21cc69d04139ffa10689780518d7a9e42dc2e819b'


def retain_barcode_reads(read1, read2, output):
    """Select untouched R1 records by the ordered IDs of fastp's passing R2 records."""
    raw_count, count = 0, 0
    with gzip.open(read1, 'rb') as raw, gzip.open(read2, 'rb') as passed, \
            gzip.open(output, 'wb', compresslevel=4) as dest:
        wanted = [passed.readline() for _ in range(4)]
        while True:
            original = [raw.readline() for _ in range(4)]
            if not original[0]:
                if any(original):
                    raise ValueError('Truncated barcode FASTQ')
                break
            raw_count += 1
            header, seq, plus, qual = original
            if not (header.startswith(b'@') and plus.startswith(b'+')
                    and len(seq.rstrip(b'\r\n')) == len(qual.rstrip(b'\r\n')) == 28):
                raise ValueError(f'Invalid 10xv3 R1 record: {raw_count}')
            if not wanted[0] or header.split()[0] != wanted[0].split()[0]:
                continue
            header, seq, plus, qual = wanted
            if not (header.startswith(b'@') and plus.startswith(b'+')
                    and seq.endswith(b'\n') and qual.endswith(b'\n')
                    and len(seq.rstrip(b'\r\n')) == len(qual.rstrip(b'\r\n'))):
                raise ValueError(f'Invalid transcript FASTQ record: {count}')
            dest.writelines(original)
            count += 1
            wanted = [passed.readline() for _ in range(4)]
        if any(wanted):
            raise ValueError('Transcript IDs are reordered or do not match the barcode reads')
    return raw_count, count


def prepare_tenx(data, binary, trimmer, threads):
    root = data / 'tenx-pbmc-1k-v3'
    record = root / 'preparation.json'
    if record.exists():
        saved = json.loads(record.read_text())
        if saved['trimmer']['sha256'] != trimmer['binary_sha256']:
            raise ValueError('Trimmer changed; use a fresh --data directory')
        for lane in saved['lanes']:
            for key in ('read1', 'read2'):
                if fingerprint(lane[key]['path']) != lane[key]:
                    raise ValueError('Prepared paired reads changed')
        return saved
    archive = root / 'pbmc_1k_v3_fastqs.tar'
    source = fetch(TENX_ARCHIVE, archive, size=5549312000, sha256=TENX_SHA256, gzip_check=False)
    onlist = fetch(TENX_ONLIST, root / '10xv3_whitelist.txt',
                   sha256=TENX_ONLIST_SHA256, gzip_check=False)
    raw, trim = root / 'raw', root / 'trimmed-fastp'
    if trim.exists() and any(trim.iterdir()):
        raise ValueError('Unrecorded trimming outputs exist; use a fresh --data directory')
    raw.mkdir(exist_ok=True)
    trim.mkdir(exist_ok=True)
    expected = {f'pbmc_1k_v3_S1_{lane}_R{mate}_001.fastq.gz'
                for lane in ('L001', 'L002') for mate in (1, 2)}
    extracted = set()
    with tarfile.open(archive) as tar:
        for member in tar:
            name = Path(member.name).name
            if name not in expected:
                continue
            if not member.isfile() or name in extracted:
                raise ValueError('Invalid or duplicate archive member')
            extracted.add(name)
            out = raw / name
            if not out.exists():
                with tar.extractfile(member) as src, out.open('wb') as dst:
                    shutil.copyfileobj(src, dst, 4 * 1024 * 1024)
            if out.stat().st_size != member.size:
                raise ValueError(f'Incomplete archive member: {name}')
    if extracted != expected:
        raise ValueError('Archive is missing a paired lane')
    lanes = []
    for lane in ('L001', 'L002'):
        read1, read2 = [raw / f'pbmc_1k_v3_S1_{lane}_R{mate}_001.fastq.gz' for mate in (1, 2)]
        out = trim / read2.name
        report = trim / (lane + '.json')
        command = [binary, '-i', str(read2), '-o', str(out), '-w', str(threads), '-z', '4',
                   '--cut_front', '--cut_tail', '--cut_mean_quality', '20', '--cut_window_size', '4',
                   '--trim_poly_g', '--length_required', '50', '--qualified_quality_phred', '15',
                   '--unqualified_percent_limit', '40', '--n_base_limit', '5', '--dont_eval_duplication',
                   '--json', str(report), '--html', str(trim / (lane + '.html'))]
        print('Trimming transcript reads:', lane, flush=True)
        start = time.perf_counter()
        with (trim / (lane + '.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        stats = json.loads(report.read_text())['summary']
        print('Validating complete paired-record order:', lane, flush=True)
        retained = trim / read1.name
        before, after = retain_barcode_reads(read1, out, retained)
        if (before != stats['before_filtering']['total_reads']
                or after != stats['after_filtering']['total_reads']):
            raise ValueError('fastp counts do not match complete paired-library validation')
        lanes.append({'lane': lane, 'raw_pairs': before, 'pairs': after,
                      'raw_read1': fingerprint(read1), 'raw_read2': fingerprint(read2),
                      'read1': fingerprint(retained), 'read2': fingerprint(out),
                      'report': fingerprint(report), 'command': command,
                      'preparation_seconds': time.perf_counter() - start})
    result = {'dataset': '10x pbmc_1k_v3', 'source_page': TENX_PAGE, 'download_url': TENX_ARCHIVE,
              'archive': source, 'reported_cells': 1222, 'technology': '10xv3',
              'barcode_bases': 16, 'umi_bases': 12, 'trimmer': fingerprint(binary),
              'trimmer_version': 'fastp ' + FASTP_VERSION, 'lanes': lanes,
              'total_pairs': sum(lane['pairs'] for lane in lanes),
              'retained_barcode_umi_sequences_unchanged': True,
              'pair_filtering': 'R2 fastp pass; retain matching original R1',
              'full_pair_order_and_gzip_validation': True, 'onlist': onlist, 'onlist_source': TENX_ONLIST}
    write_json(record, result)
    return result


def prepare(row, data, binary, trimmer, threads):
    run = row['run_accession']
    root = data / row.get('data_directory', run.lower())
    root.mkdir(parents=True, exist_ok=True)
    (root / 'ena-metadata.json').write_text(json.dumps(row, indent=2) + '\n')
    urls, hashes = row['fastq_ftp'].split(';'), row['fastq_md5'].split(';')
    sizes = row['fastq_bytes'].split(';') if row.get('fastq_bytes') else [None, None]
    if len(urls) != 2 or len(hashes) != 2 or len(sizes) != 2:
        raise RuntimeError(f'Manifest must supply two mates: {run}')
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        files = list(pool.map(lambda item: fetch(item[0], root / Path(item[0]).name, item[1], item[2]),
                              zip(urls, hashes, sizes)))
    trimmed = root / 'trimmed-fastp'
    outputs = [trimmed / f'{run}_{mate}.trimmed.fastq.gz' for mate in (1, 2)]
    command = [binary, '-i', files[0]['path'], '-I', files[1]['path'],
               '-o', str(outputs[0]), '-O', str(outputs[1]),
               '--detect_adapter_for_pe', '--cut_front', '--cut_tail', '--cut_mean_quality', '20',
               '--cut_window_size', '4', '--length_required', '50',
               '--qualified_quality_phred', '15', '--unqualified_percent_limit', '40',
               '--n_base_limit', '5', '--compression', '4', '--thread', str(threads),
               '--json', str(trimmed / 'fastp.json'), '--html', str(trimmed / 'fastp.html')]
    record = trimmed / 'preparation.json'
    if record.exists():
        saved = json.loads(record.read_text())
        if saved['command'] != command or saved['trimmer'] != trimmer:
            raise RuntimeError(f'Trimming settings changed; use a fresh --data directory: {record}')
        if any(not Path(item['path']).is_file() or digest(Path(item['path'])) != item['sha256']
               for item in saved['outputs']):
            raise RuntimeError(f'Trimmed output is missing or changed: {record}')
        return saved
    if trimmed.exists() and any(trimmed.iterdir()):
        raise RuntimeError(f'Unrecorded trimming outputs exist: {trimmed}')
    trimmed.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with (trimmed / 'fastp.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)
    report = json.loads((trimmed / 'fastp.json').read_text())
    result = {'files': files, 'trimmer': trimmer, 'command': command,
              'trim_seconds': time.perf_counter() - start,
              'pairs_before': report['summary']['before_filtering']['total_reads'] // 2,
              'pairs_after': report['summary']['after_filtering']['total_reads'] // 2,
              'outputs': [{'path': str(path), 'bytes': path.stat().st_size,
                           'sha256': digest(path)} for path in outputs]}
    record.write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=DATASETS)
    parser.add_argument('--samples', default='all', help='comma-separated ENA accessions, or all')
    parser.add_argument('--tenx-pbmc-v3', action='store_true', help='prepare the complete 10x PBMC 1k v3 library instead of ENA samples')
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--fastp', default='fastp', help='fastp executable or path')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    available = datasets(args.manifest)
    selected = list(available) if args.samples == 'all' else args.samples.split(',')
    if not selected or len(set(selected)) != len(selected) or set(selected) - available.keys():
        parser.error('select unique accessions from the dataset manifest, or all')
    if args.threads < 1:
        parser.error('--threads must be positive')
    binary = shutil.which(args.fastp)
    if not binary and args.fastp == 'fastp':
        tools = Path(__file__).resolve().parents[1] / 'target/tools'
        for name in ('fastp',):
            binary = shutil.which('fastp', path=str(tools / name / 'bin'))
            if binary:
                break
    if not binary:
        parser.error('fastp is missing; see docs/benchmarks.md#retrieve-and-trim-reads for the pinned build')
    binary = str(Path(binary).resolve())
    version = subprocess.check_output([binary, '--version'], text=True, stderr=subprocess.STDOUT).strip()
    if version != 'fastp ' + FASTP_VERSION:
        parser.error(f'expected fastp {FASTP_VERSION}, got {version!r}')
    trimmer = {'name': 'fastp', 'version': FASTP_VERSION, 'binary_sha256': digest(Path(binary))}
    data = args.data.resolve()
    if args.tenx_pbmc_v3:
        if args.samples != 'all':
            parser.error('--tenx-pbmc-v3 cannot be combined with --samples')
        print(json.dumps(prepare_tenx(data, binary, trimmer, args.threads), indent=2), flush=True)
        return
    for run in selected:
        result = prepare(available[run], data, binary, trimmer, args.threads)
        print(json.dumps({'run': run, **result}, indent=2), flush=True)


if __name__ == '__main__':
    main()
