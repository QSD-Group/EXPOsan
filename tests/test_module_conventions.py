# -*- coding: utf-8 -*-
"""Fast structure checks for EXPOsan module conventions."""

import ast
from pathlib import Path


EXPOSAN_ROOT = Path(__file__).resolve().parents[1] / 'exposan'
SKIP_PARTS = {
    '__pycache__',
    'analyses',
    'data',
    'figures',
    'publication_data',
    'readme_figures',
    'results',
}


def iter_source_files():
    for path in EXPOSAN_ROOT.rglob('*.py'):
        if SKIP_PARTS.intersection(path.parts):
            continue
        yield path


def test_public_system_entry_points_do_not_use_mutable_defaults():
    offenders = []
    for path in iter_source_files():
        tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            if node.name != 'load' and not node.name.startswith('create_system'):
                continue
            defaults = (*node.args.defaults, *(d for d in node.args.kw_defaults if d))
            if any(isinstance(d, (ast.Dict, ast.List, ast.Set)) for d in defaults):
                relpath = path.relative_to(EXPOSAN_ROOT.parent)
                offenders.append(f'{relpath}:{node.lineno} {node.name}')

    assert not offenders, 'Mutable defaults found in public entry points:\n' + '\n'.join(offenders)


def test_init_modules_does_not_create_output_dirs_by_default(tmp_path, monkeypatch):
    import exposan.utils as utils

    monkeypatch.setattr(utils, 'es_path', str(tmp_path))
    module = tmp_path / 'demo'
    module.mkdir()

    paths = utils._init_modules('demo', include_data_path=True, include_figures_path=True)

    assert paths == [
        str(module / 'data'),
        str(module / 'results'),
        str(module / 'figures'),
        ]
    assert not (module / 'results').exists()
    assert not (module / 'figures').exists()


def test_init_modules_can_create_dirs_when_requested(tmp_path, monkeypatch):
    import exposan.utils as utils

    monkeypatch.setattr(utils, 'es_path', str(tmp_path))
    module = tmp_path / 'demo'
    module.mkdir()

    paths = utils._init_modules(
        'demo',
        include_data_path=True,
        include_figures_path=True,
        create=True,
        )

    assert all(Path(p).is_dir() for p in paths)


def test_package_init_files_do_not_create_dirs_at_import_time():
    offenders = []
    for path in EXPOSAN_ROOT.rglob('__init__.py'):
        if SKIP_PARTS.intersection(path.parts):
            continue
        text = path.read_text(encoding='utf-8')
        if 'os.mkdir' in text or 'os.makedirs' in text:
            offenders.append(str(path.relative_to(EXPOSAN_ROOT.parent)))

    assert not offenders, 'Directory creation found in package initializers:\n' + '\n'.join(offenders)
