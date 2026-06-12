#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_g2rt',)

def test_g2rt():
    from numpy.testing import assert_allclose
    from exposan import g2rt

    rtol = 0.01

    modelA = g2rt.create_model('A')
    dfA = modelA.metrics_at_baseline()
    valuesA = [534.561, 1.34665, 396.737]
    assert_allclose(dfA.values, valuesA, rtol=rtol)

    modelB = g2rt.create_model('B')
    dfB = modelB.metrics_at_baseline()
    valuesB = [754.523, 1.66349, 470.437]
    assert_allclose(dfB.values, valuesB, rtol=rtol)


if __name__ == '__main__':
    test_g2rt()
