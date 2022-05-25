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
    from qsdsan import set_thermo
    from exposan import bwaise as bw

    set_thermo(bw.cmps)

    assert_allclose(bw.teaA.NPV, -42012579.5802784, rtol=1e-2)
    assert_allclose(bw.teaB.NPV, -3466006.2170442184, rtol=1e-2)
    assert_allclose(bw.teaC.NPV, -65107482.77677129, rtol=1e-2)

    assert_allclose(bw.lcaA.total_impacts['GlobalWarming'], 214197344.34534717, rtol=1e-2)
    assert_allclose(bw.lcaB.total_impacts['GlobalWarming'], 10349791.100520123, rtol=1e-2)
    assert_allclose(bw.lcaC.total_impacts['GlobalWarming'], 55187798.215826064, rtol=1e-2)


if __name__ == '__main__':
    test_bwaise()