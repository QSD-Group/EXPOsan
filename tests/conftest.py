# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

Pytest configuration: isolate biosteam's process-global utility state across
tests. Several EXPOsan modules (e.g. exposan.saf, exposan.biobinder) mutate
``bst.HeatUtility`` agent prices, ``bst.PowerUtility.price``, and ``bst.CE``
inside their ``_load_process_settings`` calls. With ``pytest-xdist
--dist loadfile`` multiple test files share a worker process, so those
mutations leak between tests and silently change TEA results in whichever
test runs later. The autouse fixture below snapshots the affected globals
before each test and restores them after.

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst
import pytest


@pytest.fixture(autouse=True)
def _reset_biosteam_global_state():
    ce = bst.CE
    power_price = bst.PowerUtility.price
    yield
    bst.CE = ce
    bst.PowerUtility.price = power_price
    bst.HeatUtility.default_agents()
