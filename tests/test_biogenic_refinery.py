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

__all__ = ('test_biogenic_refinery',)

def test_biogenic_refinery():
    from numpy.testing import assert_allclose
    from exposan import biogenic_refinery as br

    modelA = br.create_model('A')
    dfA = modelA.metrics_at_baseline()
    valuesA = [37.58, 64.69, 67.36, 10.98, 28.57, -0.02701, 10.21, -0.3066]
    assert_allclose(dfA.values, valuesA, rtol=1e-2)

    modelB = br.create_model('B')
    dfB = modelB.metrics_at_baseline()
    valuesB = [70.23, 47.37, 72.52, 26.93, 1.886, 0.03668, 0.5631, 1.452]
    assert_allclose(dfB.values, valuesB, rtol=1e-2)

    modelC = br.create_model('C')
    dfC = modelC.metrics_at_baseline()
    valuesC = [0.0, 0.0, 0.0, 11.672, 35.428, 0.049, 11.207, 0.23]
    assert_allclose(dfC.values, valuesC, rtol=1e-2)

    modelD = br.create_model('D')
    dfD = modelD.metrics_at_baseline()
    valuesD = [0.0, 0.0, 0.0, 10.017, 74.644, 0.125, 22.774, 0.525]
    assert_allclose(dfD.values, valuesD, rtol=1e-2)


if __name__ == '__main__':
    test_biogenic_refinery()