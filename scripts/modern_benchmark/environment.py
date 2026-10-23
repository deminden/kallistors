"""Checksum-pinned distributions, executable provenance and immutable source snapshots."""

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tarfile

from benchmark_common import digest, fetch, fingerprint, write_json

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'scripts/benchmarks/modern'


def read_json(path):
    return json.loads(Path(path).read_text())


def content_id(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def source_state():
    files = [p for p in (ROOT / 'crates').rglob('*') if p.is_file()]
    files += [ROOT / 'Cargo.toml', ROOT / 'Cargo.lock', ROOT / '.cargo/config.toml']
    return {str(p.relative_to(ROOT)): digest(p) for p in sorted(files)}


def check_analysis():
    lock = read_json(CONFIG / 'analysis.lock.json')
    found = {name: importlib.metadata.version(name) for name in lock['packages']}
    expected = {name: v['version'] for name, v in lock['packages'].items()}
    if found != expected:
        raise ValueError(f'Analysis dependency mismatch: expected {expected}, found {found}')
    if platform.python_version() != lock['python']:
        raise ValueError(f"Use the locked Python {lock['python']}; found {platform.python_version()}")
    return {'python': platform.python_version(), 'packages': found,
            'requirements_sha256': digest(CONFIG / 'analysis-requirements.lock')}


def platform_tag():
    arch = {'arm64': 'aarch64', 'AMD64': 'x86_64'}.get(platform.machine(), platform.machine())
    suffix = {'darwin': 'apple-darwin', 'linux': 'unknown-linux-gnu'}.get(sys.platform)
    if suffix is None:
        raise ValueError('Resource measurements require macOS or Linux')
    return arch + '-' + suffix


def safe_source(archive, destination):
    """Extract regular source files only, without archive paths escaping the destination."""
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as handle:
        for member in handle:
            path = Path(member.name)
            if path.is_absolute() or '..' in path.parts:
                raise ValueError(f'Unsafe archive entry: {member.name}')
            target = destination.joinpath(*path.parts[1:])
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                with handle.extractfile(member) as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
                target.chmod(member.mode & 0o777)
            else:
                raise ValueError(f'Unexpected non-regular source entry: {member.name}')


def install(name, entry, workspace):
    tools = workspace / 'tools'
    tools.mkdir(parents=True, exist_ok=True)
    dest = tools / 'bin' / name
    dest.parent.mkdir(exist_ok=True)
    if 'assets' in entry:
        asset = entry['assets'].get(platform_tag())
        if asset is None:
            raise ValueError(f'No pinned distribution for {name} on {platform_tag()}')
        archive = tools / 'downloads' / Path(asset['url']).name
        fetch(asset['url'], archive, sha256=asset['sha256'], size=asset['bytes'], gzip_check=False)
        with tarfile.open(archive) as handle:
            members = [m for m in handle if m.isfile() and Path(m.name).name == name]
            if len(members) != 1:
                raise ValueError(f'Expected one {name} executable in {archive}')
            with handle.extractfile(members[0]) as src, dest.open('wb') as dst:
                shutil.copyfileobj(src, dst)
        dest.chmod(0o755)
        return dest
    source = entry['source']
    archive = tools / 'downloads' / f"{name}-{entry['version']}.tar.gz"
    fetch(source['url'], archive, sha256=source['sha256'], size=source['bytes'], gzip_check=False)
    src = tools / 'source' / f"{name}-{entry['version']}"
    if not src.exists():
        safe_source(archive, src)
    log = tools / (name + '.build.log')
    with log.open('w') as handle:
        if name == 'bustools':
            build = src / 'build'
            subprocess.run(['cmake', '-S', str(src), '-B', str(build),
                            '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_POLICY_VERSION_MINIMUM=3.5'],
                           stdout=handle, stderr=handle, check=True)
            subprocess.run(['cmake', '--build', str(build), '--parallel', '8'],
                           stdout=handle, stderr=handle, check=True)
            built = build / 'src/bustools'
        else:
            subprocess.run(['make', '-j8'], cwd=src, stdout=handle, stderr=handle, check=True)
            built = src / 'fastp'
    shutil.copy2(built, dest)
    return dest


def bootstrap(workspace, overrides):
    lock = read_json(CONFIG / 'tools.lock.json')
    if (workspace / 'environment.json').exists():
        saved = load_environment(workspace)
        if saved['source_sha256'] != content_id(source_state()):
            raise ValueError('Source changed since bootstrap; use a new workspace')
        return saved['tools']
    records = {}
    suggestions = {
        'kallistors': ROOT / 'target/release/kallistors',
        'fastp': ROOT / 'target/tools/fastp/bin/fastp',
        'bustools': ROOT / 'target/tools/bustools-build/src/bustools',
    }
    for name in ['kallistors', 'salmon', 'piscem', 'piscem-infer', 'alevin-fry', 'bustools', 'fastp']:
        entry = lock['tools'].get(name, {'version': '0.5.1', 'version_args': ['version']})
        supplied = overrides.get(name)
        if supplied:
            executable = Path(supplied).absolute()
        elif suggestions.get(name, Path('/nonexistent')).is_file():
            executable = suggestions[name]
        else:
            if name == 'kallistors':
                raise ValueError('Build kallistors with cargo build --release --locked, then bootstrap')
            executable = install(name, entry, workspace)
        output = subprocess.run([str(executable), *entry['version_args']],
                                capture_output=True, text=True, check=True)
        version_text = (output.stdout + output.stderr).strip()
        if not re.search(r'(?<![\d.])' + re.escape(entry['version']) + r'(?![\d.])', version_text):
            raise ValueError(f"Expected {name} {entry['version']}, found {version_text}")
        frozen = workspace / 'bin' / name
        frozen.parent.mkdir(parents=True, exist_ok=True)
        if frozen.exists() and digest(frozen) != digest(executable):
            raise ValueError(f'{name} changed; use a new benchmark workspace')
        if not frozen.exists():
            shutil.copy2(executable, frozen)
        records[name] = {'version': entry['version'], 'version_output': version_text,
                         **fingerprint(frozen)}
    records['cellranger'] = {'version': lock['tools']['cellranger']['version'],
                             'status': 'deferred_by_user_request'}
    source = source_state()
    write_json(workspace / 'source-state.json', source)
    write_json(workspace / 'environment.json', {'tools': records, 'analysis': check_analysis(),
                                               'lock_sha256': digest(CONFIG / 'tools.lock.json'),
                                               'source_sha256': content_id(source)})
    patch = subprocess.run(['git', 'diff', '--binary'], cwd=ROOT, capture_output=True, check=True)
    (workspace / 'source.patch').write_bytes(patch.stdout)
    # A patch alone omits untracked implementation files; retain every source file too.
    snapshot = workspace / 'source'
    for relative in source:
        dest = snapshot / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, dest)
    return records


def load_environment(workspace):
    saved = read_json(workspace / 'environment.json')
    check_analysis()
    if saved['lock_sha256'] != digest(CONFIG / 'tools.lock.json'):
        raise ValueError('Tool lock changed; use a new benchmark workspace')
    for name, value in saved['tools'].items():
        if name != 'cellranger' and fingerprint(value['path']) != {
                k: value[k] for k in ('path', 'bytes', 'sha256')}:
            raise ValueError(f'Executable changed: {name}')
    return saved


def tool_paths(environment):
    return {name: value['path'] for name, value in environment['tools'].items() if 'path' in value}
