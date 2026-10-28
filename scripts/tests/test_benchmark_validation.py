"""Small corrupted-record fixtures protect the gate between comparison panels."""

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_common import fingerprint
from modern_benchmark.validation import validate_completed_panel


class CompletedPanelTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.binary = self.root / 'fixture-tool'
        self.binary.write_bytes(b'unchanged executable fixture')
        self.spec = {
            'design': {'index_threads': [8], 'builders': ['tool'], 'threads': [1, 8],
                       'single_cell_threads': [8], 'bulk_methods': ['tool'],
                       'single_cell_methods': ['tool'], 'warmups': 1, 'measured_runs': 3},
            'inputs': {'cases': [{'id': 'bulk', 'kind': 'bulk'},
                                 {'id': 'cell', 'kind': 'single_cell'}]},
            'tools': {'tool': fingerprint(self.binary)},
        }
        stage = {'command': [str(self.binary)], 'timing_interrupted': False}
        self.rows = [
            {'case': case, 'method': 'tool', 'threads': threads, 'phase': phase,
             'repetition': repetition, 'stages': [copy.deepcopy(stage)]}
            for case, threads in [('index', 8), ('bulk', 1), ('bulk', 8), ('cell', 8)]
            for phase, repetition in [('warmup', 0), ('measured', 1), ('measured', 2), ('measured', 3)]
        ]
        self.summary = {'complete_groups': 4}

    def validate(self):
        files = {'tables/publication-summary.json': self.summary, 'runs/experiment.json': self.spec}
        files.update({f'runs/row-{i}/result.json': row for i, row in enumerate(self.rows)})
        for name, value in files.items():
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(value), encoding='utf-8')
        return validate_completed_panel(self.root)

    def test_complete_panel(self):
        self.assertEqual(self.validate(), (self.summary, self.spec))

    def test_missing_observation(self):
        self.rows.pop()
        with self.assertRaises(ValueError):
            self.validate()

    def test_incomplete_summary(self):
        self.summary['complete_groups'] = 3
        with self.assertRaises(ValueError):
            self.validate()

    def test_duplicate_cannot_replace_missing_warmup(self):
        self.rows[0] = copy.deepcopy(self.rows[1])
        with self.assertRaises(ValueError):
            self.validate()

    def test_undeclared_method(self):
        self.rows[0]['method'] = 'different-tool'
        with self.assertRaises(ValueError):
            self.validate()

    def test_changed_executable(self):
        self.binary.write_bytes(b'changed executable')
        with self.assertRaisesRegex(ValueError, 'executable differs'):
            self.validate()

    def test_undeclared_cell_threads(self):
        self.rows[-1]['threads'] = 1
        with self.assertRaises(ValueError):
            self.validate()

    def test_interrupted_native_timing(self):
        self.rows[0]['stages'][0]['timing_interrupted'] = True
        with self.assertRaisesRegex(ValueError, 'Interrupted timing'):
            self.validate()

    def test_interrupted_mapping_probe(self):
        self.rows[-1]['mapping_probe'] = {'stages': [
            {'command': [str(self.binary)], 'timing_interrupted': True}]}
        with self.assertRaisesRegex(ValueError, 'Interrupted timing'):
            self.validate()

    def test_unfrozen_native_executable(self):
        self.rows[0]['stages'][0]['command'][0] = '/different/tool'
        with self.assertRaisesRegex(ValueError, 'outside the frozen experiment'):
            self.validate()

    def test_unfrozen_mapping_probe(self):
        self.rows[-1]['mapping_probe'] = {'stages': [{'command': ['/different/tool']}]}
        with self.assertRaisesRegex(ValueError, 'outside the frozen experiment'):
            self.validate()


if __name__ == '__main__':
    unittest.main()
