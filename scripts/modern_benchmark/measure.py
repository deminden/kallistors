"""Native OS accounting plus explicitly labelled sampled process-tree RSS."""
import threading
import time
from pathlib import Path
import psutil
from benchmark_common import run_process


class TreeSampler:
    def __init__(self):
        self.event = threading.Event()
        self.peak = 0
        self.samples = 0

    def start(self, pid):
        self.root = psutil.Process(pid)
        def sample():
            while not self.event.is_set():
                try:
                    family = [self.root, *self.root.children(recursive=True)]
                    rss = 0
                    for process in family:
                        try:
                            rss += process.memory_info().rss
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                    self.peak = max(self.peak, rss)
                    self.samples += 1
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
                self.event.wait(0.05)
        self.thread = threading.Thread(target=sample, daemon=True)
        self.thread.start()

    def stop(self):
        self.event.set()
        self.thread.join()


def disk_usage(root):
    files = [root] if root.is_file() else sorted(p for p in root.rglob('*') if p.is_file())
    return {'logical_bytes': sum(p.stat().st_size for p in files),
            'allocated_bytes': sum(p.stat().st_blocks * 512 for p in files), 'files': len(files)}


def pipeline(plan, out):
    out.mkdir(parents=True, exist_ok=True)
    stages = []
    start = time.perf_counter()
    for name, workers, command in plan:
        observer = TreeSampler()
        result = run_process(command, out / ('stage-' + name), observer=observer)
        result.update(stage=name, requested_workers=workers,
                      sampled_tree_peak_rss_bytes=observer.peak, rss_sample_count=observer.samples)
        stages.append(result)
    elapsed = time.perf_counter() - start
    # Native RSS is max over sequential stages; sampled tree RSS includes all descendants
    # but can double count shared pages. Never sum sequential stage peaks.
    return {'wall_seconds': sum(s['wall_seconds'] for s in stages),
            'orchestrator_wall_seconds': elapsed,
            'cpu_seconds': sum(s['process_user_seconds'] + s['process_system_seconds'] for s in stages),
            'peak_rss_bytes': max(s['peak_rss_bytes'] for s in stages),
            'sampled_tree_peak_rss_bytes': max(s['sampled_tree_peak_rss_bytes'] for s in stages),
            'mapping_wall_seconds': next((s['wall_seconds'] for s in stages if s['stage'] == 'mapping'), None),
            'stages': stages, 'disk': disk_usage(out)}
