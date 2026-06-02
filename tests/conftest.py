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

import qsdsan as qs
import pytest


@pytest.fixture(autouse=True)
def _reset_global_state():
    qs.default()
    yield
