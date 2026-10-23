"""Hash temporary native mapping products, then reclaim their disk space."""
from benchmark_common import fingerprint


def retire_intermediates(out):
    removed = []
    for path in sorted(out.rglob('*')):
        if path.is_file() and path.suffix in ('.rad', '.bus'):
            removed.append(fingerprint(path))
            path.unlink()
    # Quantifier outputs, labels, native metadata, logs and aggregate ECs remain.
    return removed
