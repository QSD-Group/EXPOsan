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

__all__ = ('test_reclaimer',)

def test_reclaimer():
    from numpy.testing import assert_allclose
    from exposan import reclaimer as re

    modelA = re.create_model('A')
    dfA = modelA.metrics_at_baseline()
    valuesA = [0.0, 0.0, 0.0, 1.618, 20.527, 0.089, 1.419,  1.168]
    assert_allclose(dfA.values, valuesA, rtol=1e-2)

    modelB = re.create_model('B')
    dfB = modelB.metrics_at_baseline()
    valuesB = [71.45, 91.79, 17.95, 104.7, 185.2, 0.09283, 31.91, 7.459]
    assert_allclose(dfB.values, valuesB, rtol=1e-2)

    modelC = re.create_model('C')
    dfC = modelC.metrics_at_baseline()
    valuesC = [71.445, 91.795, 17.949, 115.2, 145.9, -0.04704, 29.61, 5.901]
    assert_allclose(dfC.values, valuesC, rtol=1e-2)

    modelD = re.create_model('D')
    dfD = modelD.metrics_at_baseline()
    valuesD = [0.0, 0.0, 0.0, 36.539, 91.144, 0.231, 3.817, 2.762]
    assert_allclose(dfD.values, valuesD, rtol=1e-2)


if __name__ == '__main__':
    test_reclaimer()