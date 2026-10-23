#!/usr/bin/env python3
"""Run the complete original-kallisto refresh after the modern reference panels."""
import argparse
import json
from pathlib import Path
import sys
import time
import traceback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--gencode', type=Path, required=True)
    parser.add_argument('--nascent', type=Path, required=True)
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--kallisto', type=Path, default=Path('target/tools/kallisto-0.52.0'))
    parser.add_argument('--diagnostic-kallisto', type=Path, default=Path('target/tools/kallisto-src/build/src/kallisto'))
    parser.add_argument('--wait-for', type=Path, help='Wait for a preceding pipeline status JSON to say complete')
    args = parser.parse_args()
    args.workspace = args.workspace.resolve()
    args.workspace.mkdir(parents=True, exist_ok=True)
    def status(phase, **extra):
        (args.workspace / 'status.json').write_text(json.dumps({'phase': phase, **extra}, indent=2) + '\n')
        print(phase, flush=True)
    try:
        if args.wait_for:
            status('waiting for alternatives and nascent comparison to finish')
            while True:
                previous = json.loads(args.wait_for.read_text())
                if previous['phase'] == 'failed':
                    raise RuntimeError('Preceding comparison failed; original comparison held for inspection')
                if previous['phase'] == 'complete':
                    break
                time.sleep(5)
        # Import benchmark code only after waiting: preceding cohorts can install
        # a recorded analysis/reference extension before this new cohort starts.
        from original_comparison.workflow import execute
        execute(args, status)
        status('complete')
    except Exception:
        status('failed', error=traceback.format_exc())
        raise


if __name__ == '__main__':
    main()
