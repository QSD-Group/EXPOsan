# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

Pytest configuration: isolate qsdsan/biosteam process-global state across tests.
Under ``pytest-xdist --dist loadfile`` multiple test files share a worker process,
so global state leaks between them and silently changes results in whichever test
runs later:

- EXPOsan modules (e.g. exposan.saf, exposan.biobinder) mutate ``HeatUtility``
  agent prices, ``PowerUtility.price``, and ``CE`` in their process-settings;
- several systems define ``ImpactItem`` objects that share an ID but use different
  characterization factors (whichever loads first otherwise wins);
- tests/doctests that read the default flowsheet (via ``create_system`` /
  ``Flowsheet.flowsheet.default``) pick up units/streams left by an earlier test.

``qs.default()`` resets all of these before each test: a fresh 'default' flowsheet
(clearing units/streams/systems and the flowsheet-tied LCA registries), the utility
agents/prices and CEPCI, and the auto-ID ticket counters (native and LCA). Resetting
to defaults *before* each test prevents leakage just as snapshotting/restoring did,
and additionally covers the flowsheet. This keeps EXPOsan's isolation identical to
QSDsan's conftest.

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
import sys
import subprocess

import qsdsan as qs
import pytest


@pytest.fixture(autouse=True)
def _reset_global_state():
    qs.default()
    yield


# Some stiff dynamic-simulation tests are sensitive to tiny process-global state
# left by an earlier test in the same worker process: e.g. running ``test_asm``
# (which builds/simulates ASM2d) before ``test_bsm1`` perturbs the ASM2d rate
# evaluation by ~1e-15, which the BDF integration of bsm1's ASM2d/PFR config
# amplifies until it no longer converges. The ``qs.default()`` reset above does
# not clear it, and the leaked state has not been pinned to a single resettable
# object, so each of these tests is run in a fresh subprocess to guarantee
# isolation. (``pytest-forked`` would be simpler but needs ``os.fork``, which is
# unavailable on Windows.)
_ISOLATED_TEST_FILES = frozenset({
    'test_adm', 'test_asm', 'test_bsm1', 'test_bsm2',
    'test_metab', 'test_pm2', 'test_werf',
})
_SUBPROCESS_FLAG = 'EXPOSAN_ISOLATED_SUBPROCESS'


@pytest.hookimpl(tryfirst=True)
def pytest_pyfunc_call(pyfuncitem):
    # Already inside the isolated subprocess (or not a test we isolate): run it
    # normally.
    if os.environ.get(_SUBPROCESS_FLAG):
        return None
    if pyfuncitem.path.stem not in _ISOLATED_TEST_FILES:
        return None
    # Re-run just this test node in a fresh process so no earlier test's global
    # state can leak into it.
    cmd = [sys.executable, '-m', 'pytest', pyfuncitem.nodeid,
           '-q', '--no-header', '-p', 'no:cacheprovider']
    basetemp = pyfuncitem.config.option.basetemp
    if basetemp:
        safe = (pyfuncitem.nodeid.replace('/', '_').replace('\\', '_')
                .replace('::', '__').replace('.py', ''))
        cmd += ['--basetemp', os.path.join(str(basetemp), 'isolated', safe)]
    result = subprocess.run(cmd, env={**os.environ, _SUBPROCESS_FLAG: '1'})
    if result.returncode != 0:
        pytest.fail(f'{pyfuncitem.nodeid} failed in an isolated subprocess '
                    f'(exit code {result.returncode}); see output above.',
                    pytrace=False)
    return True  # signal that this hook handled the test call
