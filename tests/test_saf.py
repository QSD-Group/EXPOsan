#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_saf',)

def test_saf():
    from numpy.testing import assert_allclose
    from exposan import saf

    # Because of different CF settings for ImpactItem with the same ID
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()

    saf.load(configuration='baseline')
    rtol = 0.01
    assert_allclose(saf.get_MFSP(saf.sys), 3.95586679600505, rtol=rtol)
    assert_allclose(saf.get_GWP(saf.sys), -5.394022805849971, rtol=rtol)
    
    saf.load(configuration='EC')
    saf.simulate_and_print(saf.sys)
    assert_allclose(saf.get_MFSP(saf.sys), 11.876241988677974, rtol=rtol)
    assert_allclose(saf.get_GWP(saf.sys), 2.8357334832704386, rtol=rtol)
    
    saf.load(configuration='EC-Future')
    saf.simulate_and_print(saf.sys)
    assert_allclose(saf.get_MFSP(saf.sys), 3.821113328378629, rtol=rtol)
    assert_allclose(saf.get_GWP(saf.sys), -8.475883955624251, rtol=rtol)


if __name__ == '__main__':
    test_saf()