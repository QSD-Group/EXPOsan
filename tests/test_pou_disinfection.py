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
    assert_allclose(pou.teaA.NPV, -403.7783398531811, rtol=rtol)
    assert_allclose(pou.teaB.NPV, -3217.303126491459, rtol=rtol)
    assert_allclose(pou.teaC.NPV, -23793.25094762806, rtol=rtol)
    assert_allclose(pou.teaD.NPV, -62156.62977704776, rtol=rtol)
    assert_allclose(pou.lcaA.total_impacts['GWP'], 628.6388164710542, rtol=rtol)
    assert_allclose(pou.lcaB.total_impacts['GWP'], 353.2781043805612, rtol=rtol)
    assert_allclose(pou.lcaC.total_impacts['GWP'], 13713.001478741517, rtol=rtol)
    assert_allclose(pou.lcaD.total_impacts['GWP'], 5810.257890744621, rtol=rtol)


if __name__ == '__main__':
    test_pou_disinfection()