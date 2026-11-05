"""Archive the local sources, model inputs, and environment for a new result.

This is a source overlay, not a replacement for an installed Python environment.
Capture it before constructing/simulating a system, not retrospectively from a
plot of an older in-memory system.
"""

from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata, util
import json
from pathlib import Path
import platform
import subprocess
import sys
from zipfile import ZIP_DEFLATED, ZipFile

__all__ = ('snapshot_sources',)


# %% Repository identity and exact source payloads
def _git_state(directory):
    """Return repository identity when Git is available; tolerate installed wheels."""
    result = {}
    for key, arguments in (
            ('head', ('rev-parse', 'HEAD')),
            ('status', ('status', '--porcelain', '--untracked-files=normal'))):
        try:
            command = subprocess.run(
                ['git', '-C', str(directory), *arguments],
                capture_output=True, text=True, check=True, timeout=20,
            )
        except (OSError, subprocess.SubprocessError):
            result[key] = None
        else:
            result[key] = command.stdout.strip()
    return result


def _source_payloads():
    """Collect ESAP guides/sources and current shared model inputs/engines."""
    esap = Path(__file__).resolve().parent
    exposan = esap.parent
    payloads = {}
    roots = {'exposan': exposan}
    for folder in (esap, exposan / 'werf'):
        for path in sorted(folder.rglob('*.py')):
            payloads[f'EXPOsan/{path.relative_to(exposan.parent)}'] = path.read_bytes()
    for relative in ('__init__.py', 'utils.py',
                     'esap/REVIEW_20261015.md',
                     'werf/data/initial_conditions.xlsx'):
        path = exposan / relative
        if path.is_file():
            payloads[f'EXPOsan/exposan/{relative}'] = path.read_bytes()
    for path in sorted(esap.glob('*README.md')):
        payloads[f'EXPOsan/exposan/esap/{path.name}'] = path.read_bytes()
    tests = exposan.parent / 'tests'
    for path in sorted(tests.glob('test_esap*.py')):
        if path.is_file():
            payloads[f'EXPOsan/tests/{path.name}'] = path.read_bytes()
    for package in ('qsdsan', 'biosteam', 'thermosteam'):
        spec = util.find_spec(package)
        if spec is None or spec.origin is None:
            raise ModuleNotFoundError(f'Cannot snapshot the {package} source package.')
        root = Path(spec.origin).resolve().parent
        roots[package] = root
        for path in sorted(root.rglob('*.py')):
            payloads[f'{package}/{path.relative_to(root.parent)}'] = path.read_bytes()
        if package == 'qsdsan':
            for path in sorted((root / 'data').rglob('*.tsv')):
                payloads[f'{package}/{path.relative_to(root.parent)}'] = path.read_bytes()
    for package, root in roots.items():
        for name in ('pyproject.toml', 'LICENSE.txt', 'LICENSE'):
            path = root.parent / name
            if path.is_file():
                payloads[f'repository_metadata/{package}/{name}'] = path.read_bytes()
    return payloads, roots


# %% Result provenance (no simulation or model mutation)
def snapshot_sources(path, *, settings=None, extra_files=()):
    """Write an exclusive ZIP of sources, SHA256 hashes, inputs, and versions.

    Parameters
    ----------
    path : path-like
        New .zip path. Existing archives are never overwritten.
    settings : JSON-serializable mapping, optional
        Explicit system/simulation arguments and interpretation of the capture.
        Capturing files cannot prove how an older in-memory system was created.
    extra_files : iterable of path-like, optional
        Other inputs such as the GPS-X workbook. They are copied into the ZIP,
        indexed to avoid collisions, and retain their original paths in metadata.

    Returns
    -------
    pathlib.Path
        Absolute archive path. ``manifest.json`` contains the environment and
        a SHA256 for every other member. Archives can contain private input data.
    """
    settings = {} if settings is None else settings
    json.dumps(settings, allow_nan=False)
    payloads, roots = _source_payloads()
    inputs = {}
    for index, item in enumerate(extra_files, 1):
        source = Path(item).expanduser().resolve(strict=True)
        name = f'inputs/{index:02d}_{source.name}'
        payloads[name] = source.read_bytes()
        inputs[name] = str(source)
    manifest = {
        'schema_version': 1,
        'captured_utc': datetime.now(timezone.utc).isoformat(),
        'settings': settings,
        'inputs': inputs,
        'python': sys.version,
        'executable': sys.executable,
        'platform': platform.platform(),
        'distributions': dict(sorted((d.metadata['Name'], d.version)
                                     for d in metadata.distributions()
                                     if d.metadata['Name'])),
        'sources': {name: {'directory': str(root), **_git_state(root)}
                    for name, root in roots.items()},
        'sha256': {name: sha256(value).hexdigest()
                   for name, value in sorted(payloads.items())},
        'limitations': 'Source overlay, not a binary environment. Includes QSDsan '
                      'TSV data and WERF initial_conditions.xlsx, not all package data. '
                      'Capture before a run; do not edit code/inputs during that run.',
    }
    output = Path(path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(output, 'x', compression=ZIP_DEFLATED) as archive:
        for name, value in sorted(payloads.items()):
            archive.writestr(name, value)
        archive.writestr('manifest.json', json.dumps(manifest, indent=2, allow_nan=False))
    return output
