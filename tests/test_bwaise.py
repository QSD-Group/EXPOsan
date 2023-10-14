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

__all__ = ('test_bwaise',)

def test_bwaise():
    from numpy.testing import assert_allclose
    from exposan import bwaise as bw

    # Because of different CF settings for ImpactItem with the same ID
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()

    bw.load()
    rtol = 0.01
    assert_allclose(bw.teaA.NPV, -42012332.26576174, rtol=rtol)
    assert_allclose(bw.teaB.NPV, -3452718.4819053616, rtol=rtol)
    assert_allclose(bw.teaC.NPV, -65123672.010664254, rtol=rtol)
    assert_allclose(bw.lcaA.total_impacts['GlobalWarming'], 214197344.34531045, rtol=rtol)
    assert_allclose(bw.lcaB.total_impacts['GlobalWarming'], 10181296.352069596, rtol=rtol)
    assert_allclose(bw.lcaC.total_impacts['GlobalWarming'], 55327236.9006904, rtol=rtol)


if __name__ == '__main__':
    test_bwaise()