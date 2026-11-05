"""Exact source/input preservation without requiring a new simulation."""

from hashlib import sha256
import json
from zipfile import ZipFile

import pytest

from exposan.esap import provenance


@pytest.fixture
def minimal_sources(monkeypatch, tmp_path):
    payloads = {'EXPOsan/exposan/esap/system.py': b'calibration = 24.5\n'}
    monkeypatch.setattr(provenance, '_source_payloads', lambda: (
        payloads, {'exposan': tmp_path}))
    monkeypatch.setattr(provenance, '_git_state', lambda path: {
        'head': 'base-commit', 'status': '?? exposan/esap/'})
    return payloads


def test_snapshot_preserves_sources_inputs_settings_and_hashes(minimal_sources, tmp_path):
    sample = tmp_path / 'sample.csv'
    sample.write_bytes(b'NH4,1380.06\n')
    path = provenance.snapshot_sources(
        tmp_path / 'sources.zip', settings={'ad_initial_s_ip': 24.5},
        extra_files=[sample],
    )
    sample.write_bytes(b'changed later\n')
    with ZipFile(path) as archive:
        manifest = json.loads(archive.read('manifest.json'))
        assert manifest['settings'] == {'ad_initial_s_ip': 24.5}
        assert manifest['sources']['exposan']['head'] == 'base-commit'
        assert manifest['sources']['exposan']['status'].startswith('??')
        assert manifest['python']
        assert manifest['distributions']['numpy']
        assert archive.read('inputs/01_sample.csv') == b'NH4,1380.06\n'
        for name, digest in manifest['sha256'].items():
            assert sha256(archive.read(name)).hexdigest() == digest


def test_snapshot_never_overwrites_an_existing_archive(minimal_sources, tmp_path):
    path = provenance.snapshot_sources(tmp_path / 'sources.zip')
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        provenance.snapshot_sources(path)
    assert path.read_bytes() == original


def test_snapshot_rejects_nonfinite_settings_before_writing(minimal_sources, tmp_path):
    path = tmp_path / 'invalid.zip'
    with pytest.raises(ValueError):
        provenance.snapshot_sources(path, settings={'atol': float('nan')})
    assert not path.exists()


def test_snapshot_requires_real_additional_inputs(minimal_sources, tmp_path):
    with pytest.raises(FileNotFoundError):
        provenance.snapshot_sources(tmp_path / 'missing.zip',
                                    extra_files=[tmp_path / 'missing.csv'])
    assert not (tmp_path / 'missing.zip').exists()


def test_actual_source_inventory_includes_engines_and_initial_data():
    payloads, roots = provenance._source_payloads()
    assert set(roots) == {'exposan', 'qsdsan', 'biosteam', 'thermosteam'}
    for name in (
            'EXPOsan/exposan/esap/systems/SVCW.py',
            'EXPOsan/exposan/esap/README.md',
            'EXPOsan/exposan/esap/tricklingfilter_README.md',
            'EXPOsan/exposan/esap/SVCW_README.md',
            'EXPOsan/exposan/esap/sandfilter_README.md',
            'EXPOsan/exposan/werf/data/initial_conditions.xlsx',
            'qsdsan/qsdsan/process_models/_adm1_p_extension.py',
            'qsdsan/qsdsan/unit_operations/dynamic/_junction.py',
            'qsdsan/qsdsan/data/process_data/_masm2d.tsv',
            'biosteam/biosteam/__init__.py',
            'thermosteam/thermosteam/__init__.py'):
        assert payloads[name], name
