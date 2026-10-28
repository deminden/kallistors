#!/usr/bin/env python3
"""Check local documentation links, shell examples and Python syntax without builds or network."""

import re
import shutil
import subprocess
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_LINK = 'https://github.com/deminden/kallistors/blob/main/'
FENCE = re.compile(r'^\s*(`{3,}|~{3,})(.*)$')
LINK = re.compile(r'!?\[[^\]]*\]\((<[^>]+>|[^\s)]+)(?:\s+"[^"]*")?\)')


def markdown_parts(text):
    """Keep code examples out of link/heading checks; retain shell blocks for syntax checks."""
    prose, shell, code = [], [], []
    fence = None
    language = ''
    for line in text.splitlines():
        marker = FENCE.match(line)
        if fence:
            if marker and marker[1][0] == fence[0] and len(marker[1]) >= len(fence) and not marker[2].strip():
                if language in ('bash', 'sh'):
                    shell.append('\n'.join(code) + '\n')
                fence, code = None, []
            else:
                code.append(line)
        elif marker:
            fence, language = marker[1], marker[2].strip()
        else:
            prose.append(line)
    if fence:
        raise ValueError('unclosed Markdown code fence')
    return '\n'.join(prose), shell


def heading_ids(text):
    # Match GitHub's punctuation stripping and duplicate-heading suffixes.
    ids = set()
    for title in re.findall(r'^#{1,6}\s+(.+?)\s*#*$', text, re.MULTILINE):
        title = re.sub(r'<[^>]+>', '', title).lower()
        base = re.sub(r'[^\w\- ]', '', title).replace(' ', '-')
        anchor, suffix = base, 0
        while anchor in ids:
            suffix += 1
            anchor = f'{base}-{suffix}'
        ids.add(anchor)
    return ids


def check(root=ROOT):
    errors = []
    pages = sorted({*root.glob('*.md'), *root.glob('docs/**/*.md'),
                    *root.glob('crates/*/README.md')})
    documents = {}
    for path in pages:
        try:
            documents[path.resolve()] = markdown_parts(path.read_text(encoding='utf-8'))
        except ValueError as error:
            errors.append(f'{path.relative_to(root)}: {error}')
    bash = shutil.which('bash')
    if not bash:
        errors.append('bash is required to check shell-example syntax')
    links, examples = 0, 0
    for path, (prose, blocks) in documents.items():
        for match in LINK.finditer(prose):
            target = match[1].strip('<>')
            if target.startswith(REPOSITORY_LINK):
                base, target = root, target[len(REPOSITORY_LINK):]
            elif re.match(r'\w+://|mailto:', target):
                continue
            else:
                base = path.parent
            destination, _, anchor = target.partition('#')
            linked = (base / unquote(destination)).resolve() if destination else path
            links += 1
            if not linked.is_file():
                errors.append(f'{path.relative_to(root)}: missing link {match[1]}')
            elif anchor and linked.suffix == '.md':
                linked_prose = documents.get(linked, ('', []))[0]
                if unquote(anchor) not in heading_ids(linked_prose):
                    errors.append(f'{path.relative_to(root)}: missing anchor {match[1]}')
        for number, block in enumerate(blocks, 1):
            examples += 1
            if bash:
                result = subprocess.run([bash, '-n'], input=block, text=True,
                                        capture_output=True, timeout=10)
                if result.returncode:
                    errors.append(f'{path.relative_to(root)}: shell example {number}: {result.stderr.strip()}')
    python_files = sorted((root / 'scripts').rglob('*.py'))
    for path in python_files:
        try:
            compile(path.read_bytes(), str(path), 'exec')
        except SyntaxError as error:
            errors.append(f'{path.relative_to(root)}: {error}')
    for error in errors:
        print(error)
    print(f'Checked {len(pages)} Markdown files, {links} links/anchors, '
          f'{examples} shell examples and {len(python_files)} Python files.')
    return not errors


if __name__ == '__main__':
    raise SystemExit(0 if check() else 1)
