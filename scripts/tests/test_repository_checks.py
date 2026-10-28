"""Protect the checker against accepting broken links or inspecting literal code as Markdown."""

from contextlib import redirect_stdout
import io
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_repository import check


class RepositoryChecks(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name).resolve()
        (self.root / 'docs').mkdir()
        (self.root / 'scripts').mkdir()
        (self.root / 'docs/guide.md').write_text('# Usage\n\n## Examples\n\n## Examples\n')

    def inspect(self, markdown):
        (self.root / 'README.md').write_text(markdown)
        with redirect_stdout(io.StringIO()):
            return check(self.root)

    def test_local_and_packaged_links_with_duplicate_anchors(self):
        self.assertTrue(self.inspect(
            '[Local](docs/guide.md#examples-1)\n'
            '[Packaged](https://github.com/deminden/kallistors/blob/main/docs/guide.md#usage)\n'
            '```text\n[Literal example](missing.md)\n```\n'))

    def test_missing_destination(self):
        self.assertFalse(self.inspect('[Missing](docs/missing.md)\n'))

    def test_missing_heading(self):
        self.assertFalse(self.inspect('[Missing](docs/guide.md#missing)\n'))

    def test_unclosed_fence(self):
        self.assertFalse(self.inspect('```bash\necho unfinished\n'))

    def test_invalid_shell_example(self):
        self.assertFalse(self.inspect('```bash\nif true; then\n```\n'))

    def test_python_syntax_error(self):
        (self.root / 'scripts/broken.py').write_text('def unfinished(\n')
        self.assertFalse(self.inspect('# Valid Markdown\n'))


if __name__ == '__main__':
    unittest.main()
