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

__all__ = ('test_eco_san',)

def test_eco_san():
    from numpy.testing import assert_allclose
    from exposan import eco_san as es

    # Because of different CF settings for ImpactItem with the same ID
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()

    es.load()
    assert_allclose(es.teaA.NPV, -99282.30904085489, rtol=1e-2)
    assert_allclose(es.teaB.NPV, -78656.64103836296, rtol=1e-2)
    assert_allclose(es.teaC.NPV, -64039.20450181706, rtol=1e-2)
    assert_allclose(es.lcaA.total_impacts['GlobalWarming'], 284386.4993382638, rtol=1e-2)
    assert_allclose(es.lcaB.total_impacts['GlobalWarming'], 205035.34951036863, rtol=1e-2)
    assert_allclose(es.lcaC.total_impacts['GlobalWarming'], 243717.22113291494, rtol=1e-2)


if __name__ == '__main__':
    test_eco_san()