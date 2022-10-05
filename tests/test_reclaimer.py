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

    # # Without resource recovery
    # re.INCLUDE_RESOURCE_RECOVERY = False
    
    # modelA = re.create_model('A')
    # dfA = modelA.metrics_at_baseline()
    # valuesA = [0.0, 0.0, 0.0, 4.746, 32.468, 0.131, 2.138, 1.66]
    # assert_allclose(dfA.values, valuesA, rtol=1e-2)

    # modelB = re.create_model('B')
    # dfB = modelB.metrics_at_baseline()
    # valuesB =  [71.445, 91.795, 17.949, 78.061, 200.883, 0.369, 37.283, 6.672]
    # assert_allclose(dfB.values, valuesB, rtol=1e-2)

    # modelC = re.create_model('C')
    # dfC = modelC.metrics_at_baseline()
    # valuesC =  [71.445, 91.795, 17.949, 82.457, 150.254, 0.189, 34.307, 4.657]
    # assert_allclose(dfC.values, valuesC, rtol=1e-2)

    # modelD = re.create_model('D')
    # dfD = modelD.metrics_at_baseline()
    # valuesD = [0.0, 0.0, 0.0, 39.809, 91.145, 0.231, 3.817, 2.762]
    # assert_allclose(dfD.values, valuesD, rtol=1e-2)
    
    # With resource recovery
    re.INCLUDE_RESOURCE_RECOVERY = True
    
    modelA2 = re.create_model('A')
    dfA2 = modelA2.metrics_at_baseline()
    # Same results with/without resource recovery
    valuesA2 = [0.0, 0.0, 0.0, 4.746, 32.468, 0.131, 2.138, 1.66]
    assert_allclose(dfA2.values, valuesA2, rtol=1e-2)

    modelB2 = re.create_model('B')
    dfB2 = modelB2.metrics_at_baseline()
    valuesB2 = [71.45, 91.79, 17.95, 74.74, 178.8, 0.0422, 31.49, 3.347]
    assert_allclose(dfB2.values, valuesB2, rtol=1e-2)

    modelC2 = re.create_model('C')
    dfC2 = modelC2.metrics_at_baseline()
    valuesC2 = [71.445, 91.795, 17.949, 79.134, 128.123, -0.138, 28.512, 1.332]
    assert_allclose(dfC2.values, valuesC2, rtol=1e-2)

    modelD2 = re.create_model('D')
    dfD2 = modelD2.metrics_at_baseline()
    # Same results with/without resource recovery
    valuesD2 = [0.0, 0.0, 0.0, 39.809, 91.144, 0.231, 3.817, 2.762]
    assert_allclose(dfD2.values, valuesD2, rtol=1e-2)


if __name__ == '__main__':
    test_reclaimer()