#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pytest

__all__ = (
    'test_require_package_passes_for_installed_package',
    'test_require_package_raises_friendly_error_for_missing_package',
    'test_get_generic_tanker_truck_fee',
    )


def test_require_package_passes_for_installed_package():
    from exposan.utils import require_package
    # numpy is always installed (transitive via qsdsan); must not raise
    require_package('numpy')


def test_require_package_raises_friendly_error_for_missing_package():
    from exposan.utils import require_package
    with pytest.raises(ModuleNotFoundError) as exc_info:
        require_package(
            'a-package-that-does-not-exist-1234',
            feature='a made-up feature',
        )
    msg = str(exc_info.value)
    assert 'a made-up feature' in msg
    assert 'a-package-that-does-not-exist-1234' in msg
    assert 'pip install a-package-that-does-not-exist-1234' in msg
    assert 'pip install exposan[complete]' in msg


def test_get_generic_tanker_truck_fee():
    from exposan.utils import get_generic_tanker_truck_fee
    fee = get_generic_tanker_truck_fee(6)
    assert fee > 0
