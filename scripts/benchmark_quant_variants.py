#!/usr/bin/env python3
"""Alternate full-process timings of two quantification binaries on identical paired reads."""

import argparse
import json
import re
import sys
from pathlib import Path
import shutil
import statistics
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--control', required=True, type=Path)
parser.add_argument('--candidate', type=Path, default=Path('target/release/kallistors'))
parser.add_argument('--out', required=True, type=Path)
parser.add_argument('--runs', type=int, default=3)
parser.add_argument('--warmups', type=int, default=1)
parser.add_argument('--candidate-first', action='store_true')
parser.add_argument('--no-timings', action='store_true', help='omit the Rust --timings flag when including original kallisto')
parser.add_argument('--threads', type=int, default=8)
parser.add_argument('--memory', action='store_true', help='measure whole-process peak RSS using /usr/bin/time (macOS/Linux)')
parser.add_argument('--index', type=Path, default=Path('data/srr-small/gencode.v50.kallisto.idx'))
parser.add_argument('--reads', type=Path, default=Path('data/srr13638690/trimmed/SRR13638690_1.trimmed.fastq.gz'))
parser.add_argument('--reads2', type=Path, default=Path('data/srr13638690/trimmed/SRR13638690_2.trimmed.fastq.gz'))
args = parser.parse_args()
if args.runs < 1 or args.warmups < 0 or args.threads < 1:
    parser.error('runs/threads must be positive and warmups nonnegative')
if args.out.exists() and any(args.out.iterdir()):
    parser.error('--out must name an empty or new directory to preserve previous evidence')
if args.memory and (sys.platform not in ('darwin', 'linux') or not Path('/usr/bin/time').is_file()):
    parser.error('--memory requires macOS/Linux and /usr/bin/time')
args.out.mkdir(parents=True, exist_ok=True)
binaries = {'control': args.control.resolve(), 'candidate': args.candidate.resolve()}
if args.candidate_first:
    binaries = dict(reversed(list(binaries.items())))
rows = []
expected_counters = None

def run(tool, phase, replicate):
    global expected_counters
    cpu_context = []
    for line in subprocess.check_output(['ps', '-axo', 'pcpu,comm'], text=True).splitlines()[1:]:
        fields = line.strip().split(None, 1)
        if len(fields) == 2 and float(fields[0]) >= 25:
            cpu_context.append({'cpu_percent': float(fields[0]), 'process': Path(fields[1]).name})
    out = args.out / f'{phase}-{tool}-{replicate}'
    command = [str(binaries[tool]), 'quant', '-t', str(args.threads),
               '-i', str(args.index.resolve()), '-o', str(out.resolve()),
               str(args.reads.resolve()), str(args.reads2.resolve())]
    if not args.no_timings:
        command.insert(2, '--timings')
    log_path = args.out / f'{phase}-{tool}-{replicate}.log'
    execution_command = command
    memory_path = args.out / f'{phase}-{tool}-{replicate}.memory.log'
    if args.memory:
        time_flags = ['-l'] if sys.platform == 'darwin' else ['-f', 'peak_rss_kib=%M']
        execution_command = ['/usr/bin/time', *time_flags, '-o', str(memory_path.resolve()), *command]
    with log_path.open('w') as log:
        start = time.perf_counter()
        subprocess.run(execution_command, stdout=log, stderr=log, check=True)
        elapsed = time.perf_counter() - start
    info = json.loads((out / 'run_info.json').read_text())
    counters = {field: info[field] for field in ('n_processed', 'n_pseudoaligned', 'n_unique')}
    if expected_counters is None:
        expected_counters = counters
    if counters != expected_counters:
        raise RuntimeError(f'Counter mismatch: {counters} != {expected_counters}')
    stage = {}
    for line in log_path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) == 2 and fields[1].endswith('s'):
            try: stage[fields[0]] = float(fields[1][:-1])
            except ValueError: pass
    peak_rss_bytes = None
    peak_memory_footprint_bytes = None
    if args.memory:
        memory_text = memory_path.read_text()
        pattern = r'(\d+)\s+maximum resident set size' if sys.platform == 'darwin' else r'peak_rss_kib=(\d+)'
        match = re.search(pattern, memory_text)
        if not match:
            raise RuntimeError(f'Cannot parse peak RSS from {memory_path}')
        peak_rss_bytes = int(match.group(1)) * (1 if sys.platform == 'darwin' else 1024)
        footprint = re.search(r'(\d+)\s+peak memory footprint', memory_text)
        if footprint:
            peak_memory_footprint_bytes = int(footprint.group(1))
    row = {'execution_command': execution_command, 'peak_rss_bytes': peak_rss_bytes,
           'peak_memory_footprint_bytes': peak_memory_footprint_bytes, 'command': command, 'tool': tool, 'phase': phase, 'replicate': replicate, 'wall_seconds': elapsed,
           'counters': counters, 'stages': stage, 'other_cpu_before_run': cpu_context}
    rows.append(row)
    (args.out / 'runs.json').write_text(json.dumps(rows, indent=2) + '\n')
    print(f'{phase} {tool} {replicate}: {elapsed:.6f}s; EM={stage.get("em")}; '
          f'mapping={stage.get("pseudoalign")}; peak_RSS_bytes={peak_rss_bytes}', flush=True)
    if phase != 'measured' or replicate != args.runs:
        shutil.rmtree(out)
    return elapsed

for warmup in range(args.warmups):
    for tool in binaries: run(tool, 'warmup', warmup)
for replicate in range(1, args.runs + 1):
    order = list(binaries)
    if replicate % 2 == 0: order.reverse()
    for tool in order: run(tool, 'measured', replicate)
medians = {tool: statistics.median(r['wall_seconds'] for r in rows
           if r['tool'] == tool and r['phase'] == 'measured') for tool in binaries}
speedup = medians['control'] / medians['candidate']
result = {'binaries': {tool: str(binary) for tool, binary in binaries.items()},
          'threads': args.threads, 'index': str(args.index.resolve()),
          'reads': [str(args.reads.resolve()), str(args.reads2.resolve())],
          'expected_counters': expected_counters, 'medians': medians, 'speedup': speedup, 'speed_gain_percent': 100 * (speedup - 1),
          'runs': args.runs, 'warmups_per_binary': args.warmups,
          'timing_qualified': args.runs >= 3 and args.warmups >= 1 and speedup >= 1.02}
if args.memory:
    result['memory_measurement'] = 'maximum resident set size from /usr/bin/time; values normalized to bytes'
    result['median_peak_rss_bytes'] = {
        tool: statistics.median(row['peak_rss_bytes'] for row in rows
                               if row['tool'] == tool and row['phase'] == 'measured')
        for tool in binaries
    }
    result['peak_rss_ratio_candidate_to_control'] = (
        result['median_peak_rss_bytes']['candidate'] / result['median_peak_rss_bytes']['control']
    )
(args.out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2), flush=True)
