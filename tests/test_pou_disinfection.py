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

__all__ = ('test_pou_disinfection',)

def test_pou_disinfection():
    from numpy.testing import assert_allclose
    from exposan import pou_disinfection as pou

    # Because of different CF settings for ImpactItem with the same ID
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()

    pou.load()
    rtol = 0.01
    assert_allclose(pou.teaA.NPV, -272.3341764464383, rtol=rtol)
    assert_allclose(pou.teaB.NPV, -2144.868750994307, rtol=rtol)
    assert_allclose(pou.teaC.NPV, -16026.374451801337, rtol=rtol)
    assert_allclose(pou.teaD.NPV, -41465.286859741536, rtol=rtol)
    assert_allclose(pou.lcaA.total_impacts['GWP'], 419.1941587933792, rtol=rtol)
    assert_allclose(pou.lcaB.total_impacts['GWP'], 235.51873625370752, rtol=rtol)
    assert_allclose(pou.lcaC.total_impacts['GWP'], 9141.66937353872, rtol=rtol)
    assert_allclose(pou.lcaD.total_impacts['GWP'], 3694.7137460627887, rtol=rtol)


if __name__ == '__main__':
    test_pou_disinfection()