"""Shared download verification, benchmark execution and result records."""

import gzip
import hashlib
import json
import os
import platform
from pathlib import Path
import re
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time

DATASETS = Path(__file__).parent / 'benchmarks/human_datasets.json'
ORIGINAL_MARKERS = {'index_end': b'[index] number of k-mers:',
                    'reads_begin': b'[quant] finding pseudoalignments for the reads ...',
                    'reads_end': b'[quant] processed ',
                    'em_begin': b'[   em] quantifying the abundances ...',
                    'em_end': b'[   em] the Expectation-Maximization algorithm ran for '}
RUST_MARKERS = {key: f'[quant-stage] {key}'.encode() for key in ORIGINAL_MARKERS}


def digest(path, algorithm='sha256'):
    result = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError(f'Not a file: {path}')
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest(path)}


def executable(value):
    found = shutil.which(str(value))
    if not found:
        raise ValueError(f'Executable not found: {value}')
    # Preserve a virtual environment interpreter symlink when invoking Python.
    return Path(found).absolute()


def fetch(url, path, md5=None, size=None, sha256=None, gzip_check=True):
    path = Path(path)
    url = url if '://' in url else 'https://' + url
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file():
        partial = path.with_suffix(path.suffix + '.partial')
        subprocess.run(['curl', '-fL', '--retry', '5', '--retry-delay', '2',
                        '--continue-at', '-', url, '-o', str(partial)], check=True)
        partial.replace(path)
    if size is not None and path.stat().st_size != int(size):
        raise ValueError(f'Download size mismatch: {path}')
    if md5 and digest(path, 'md5') != md5:
        raise ValueError(f'Download MD5 mismatch: {path}')
    if sha256 and digest(path) != sha256:
        raise ValueError(f'Download SHA-256 mismatch: {path}')
    if gzip_check:
        with gzip.open(path, 'rb') as handle:
            while handle.read(4 * 1024 * 1024):
                pass
    return {**fingerprint(path), 'url': url, **({'md5': md5} if md5 else {})}


def datasets(path=DATASETS):
    rows = json.loads(Path(path).read_text())
    result = {}
    for row in rows:
        name = row['run_accession']
        if name in result or not re.fullmatch(r'[A-Z]+[0-9]+', name):
            raise ValueError(f'Duplicate or invalid accession: {name}')
        if row['library_layout'] != 'PAIRED' or any(len(row[key].split(';')) != 2 for key in ('fastq_ftp', 'fastq_md5')):
            raise ValueError(f'Expected paired ENA records: {name}')
        directory = Path(row.get('data_directory', name.lower()))
        if directory.is_absolute() or '..' in directory.parts:
            raise ValueError(f'Dataset directory must stay inside --data: {name}')
        result[name] = row
    if not result:
        raise ValueError('Empty dataset manifest')
    return result


def prepared_reads(row, data, subset=False):
    root = Path(data) / row.get('data_directory', row['run_accession'].lower()) / 'trimmed-fastp'
    names = [f'subset/mate{i}.fastq.gz' if subset else f'{row["run_accession"]}_{i}.trimmed.fastq.gz'
             for i in (1, 2)]
    return [root / name for name in names]


def timing_environment():
    return {key: value for key, value in os.environ.items() if not key.startswith('KALLISTO_DIAGNOSTIC_')}


def cpu_snapshot(pid=None):
    if not shutil.which('ps'):
        return []
    output = subprocess.run(['ps', '-axo', 'pid,ppid,pcpu,comm'], capture_output=True, text=True)
    if output.returncode:
        return []
    rows = [line.split(None, 3) for line in output.stdout.splitlines()[1:]]
    excluded = {os.getpid(), pid}
    while True:
        children = {int(row[0]) for row in rows if len(row) == 4 and int(row[1]) in excluded}
        if children <= excluded:
            break
        excluded |= children
    return [{'pid': int(row[0]), 'cpu_percent': float(row[2]), 'process': Path(row[3]).name}
            for row in rows if len(row) == 4 and int(row[0]) not in excluded and float(row[2]) >= 25]


def system_context():
    commands = ({'swap': ['sysctl', 'vm.swapusage'], 'vm_stat': ['vm_stat'], 'power': ['pmset', '-g', 'batt']}
                if sys.platform == 'darwin' else {})
    result = {'platform': sys.platform, 'architecture': platform.machine(), 'cpu_count': os.cpu_count()}
    if sys.platform == 'darwin' and result['architecture'] in ('arm64', 'aarch64'):
        result['cpu_family'] = 'Apple Silicon'
    for key, command in commands.items():
        output = subprocess.run(command, capture_output=True, text=True)
        result[key] = output.stdout.strip() if output.returncode == 0 else None
    if sys.platform == 'linux':
        result['meminfo'] = Path('/proc/meminfo').read_text()
    return result


def phase_durations(events, elapsed):
    keys = ('index_end', 'reads_begin', 'reads_end', 'em_begin', 'em_end')
    missing = set(keys) - events.keys()
    if missing:
        raise ValueError(f'Missing phase markers: {sorted(missing)}')
    values = [0.0, *(events[key] for key in keys), elapsed]
    if values != sorted(values):
        raise ValueError(f'Out-of-order phase markers: {events}')
    names = ('startup_index', 'input_setup', 'read_processing', 'model_setup', 'em', 'output_shutdown')
    result = {key: end - begin for key, begin, end in zip(names, values, values[1:])}
    result['model_and_em'] = result['model_setup'] + result['em']
    return result


def run_process(command, prefix, memory=True, markers=None, environment=None):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    memory_file = prefix.with_suffix('.memory.log')
    execution = list(map(str, command))
    if memory:
        if sys.platform not in ('darwin', 'linux') or not Path('/usr/bin/time').is_file():
            raise ValueError('RSS measurement requires macOS/Linux /usr/bin/time; use --no-memory')
        flags = ['-l'] if sys.platform == 'darwin' else ['-f', 'wall=%e\nuser=%U\nsystem=%S\npeak_rss_kib=%M']
        execution = ['/usr/bin/time', *flags, '-o', str(memory_file), *execution]
    before, cpu_before = system_context(), cpu_snapshot()
    events, samples = {}, []
    stopped = threading.Event()
    start = time.perf_counter()
    process = subprocess.Popen(execution, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               env=environment or timing_environment(), start_new_session=os.name == 'posix')

    def sample():
        while not stopped.wait(5):
            samples.append({'elapsed_seconds': time.perf_counter() - start,
                            'other_processes': cpu_snapshot(process.pid)})

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    previous_sigterm = None
    if threading.current_thread() is threading.main_thread():
        previous_sigterm = signal.getsignal(signal.SIGTERM)

        def terminate(signum, _frame):
            raise SystemExit(128 + signum)

        signal.signal(signal.SIGTERM, terminate)
    try:
        pending = b''
        with prefix.with_suffix('.log').open('wb') as log:
            while chunk := os.read(process.stdout.fileno(), 65536):
                stamp = time.perf_counter() - start
                log.write(chunk)
                if markers:
                    pending += chunk
                    for key, marker in markers.items():
                        if key not in events and marker in pending:
                            events[key] = stamp
                    pending = pending[-4096:]
        code = process.wait()
        elapsed = time.perf_counter() - start
    finally:
        stopped.set()
        sampler.join()
        process.stdout.close()
        if process.poll() is None:
            if os.name == 'posix':
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
            process.wait()
        if previous_sigterm is not None:
            signal.signal(signal.SIGTERM, previous_sigterm)
    if code:
        raise RuntimeError(f'Command failed ({code}); see {prefix.with_suffix(".log")}')
    result = {'command': list(map(str, command)), 'wall_seconds': elapsed,
              'events_seconds': events, 'other_cpu_before_run': cpu_before,
              'other_cpu_after_run': cpu_snapshot(), 'other_cpu_during_run': samples,
              'memory_before_run': before, 'memory_after_run': system_context(),
              'peak_rss_bytes': None, 'peak_memory_footprint_bytes': None,
              'timing_interrupted': False}
    if markers:
        result['stages_seconds'] = phase_durations(events, elapsed)
    if memory:
        text = memory_file.read_text()
        if sys.platform == 'darwin':
            rss = re.search(r'(\d+)\s+maximum resident set size', text)
            times = re.search(r'^\s*([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys', text, re.M)
            footprint = re.search(r'(\d+)\s+peak memory footprint', text)
            if footprint:
                result['peak_memory_footprint_bytes'] = int(footprint[1])
            if not rss or not times:
                raise ValueError(f'Cannot parse macOS time output: {memory_file}')
            real, user, system = map(float, times.groups())
            result['peak_rss_bytes'] = int(rss[1])
        else:
            fields = dict(re.findall(r'^(wall|user|system|peak_rss_kib)=([\d.]+)$', text, re.M))
            real, user, system = (float(fields[key]) for key in ('wall', 'user', 'system'))
            result['peak_rss_bytes'] = int(fields['peak_rss_kib']) * 1024
        gap = real - elapsed
        result.update(os_real_seconds=real, process_user_seconds=user, process_system_seconds=system,
                      os_real_minus_stage_clock_seconds=gap, timing_interrupted=gap > 2)
    return result


def summarize(rows):
    result = {'tools': {}}
    for tool in ('kallisto', 'kallistors'):
        measured = [row for row in rows if row['tool'] == tool and row['phase'] == 'measured']
        if not measured:
            raise ValueError(f'No measured runs for {tool}')
        times = [row['wall_seconds'] for row in measured]
        item = {'elapsed_seconds': times, 'median_seconds': statistics.median(times),
                'min_seconds': min(times), 'max_seconds': max(times),
                'sample_sd_seconds': statistics.stdev(times) if len(times) > 1 else None}
        for field in ('peak_rss_bytes', 'peak_memory_footprint_bytes', 'index_bytes'):
            values = [row.get(field) for row in measured]
            item['median_' + field] = statistics.median(values) if all(v is not None for v in values) else None
        if all('stages_seconds' in row for row in measured):
            item['median_stages_seconds'] = {key: statistics.median(row['stages_seconds'][key] for row in measured)
                                             for key in measured[0]['stages_seconds']}
        if all('build_stages_seconds' in row for row in measured):
            item['median_build_stages_seconds'] = {key: statistics.median(row['build_stages_seconds'][key] for row in measured)
                                                  for key in measured[0]['build_stages_seconds']}
        result['tools'][tool] = item
    result['total_speedup'] = result['tools']['kallisto']['median_seconds'] / result['tools']['kallistors']['median_seconds']
    result['all_measured_runs_uninterrupted'] = all(not row['timing_interrupted'] for row in rows if row['phase'] == 'measured')
    return result
