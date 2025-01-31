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

    # Because of different CF settings for ImpactItem with the same ID
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()

    rtol = 0.01

    # # Without resource recovery
    # br.INCLUDE_RESOURCE_RECOVERY = False

    # modelA = br.create_model('A')
    # dfA = modelA.metrics_at_baseline()
    # valuesA = [37.58, 64.69, 67.36, 11.19, 36.84, 0.04949, 11.78, 0.2286]
    # assert_allclose(dfA.values, valuesA, rtol=rtol)

    # modelB = br.create_model('B')
    # dfB = modelB.metrics_at_baseline()
    # valuesB = [70.234, 47.374, 72.52, 27.153, 14.368, 0.143, 2.676, 2.113]
    # assert_allclose(dfB.values, valuesB, rtol=rtol)

    # modelC = br.create_model('C')
    # dfC = modelC.metrics_at_baseline()
    # valuesC = [0.0, 0.0, 0.0, 11.864, 36.857, 0.05, 11.785, 0.23]
    # assert_allclose(dfC.values, valuesC, rtol=rtol)

    # modelD = br.create_model('D')
    # dfD = modelD.metrics_at_baseline()
    # valuesD = [0.0, 0.0, 0.0, 10.017, 74.644, 0.125, 22.774, 0.525]
    # assert_allclose(dfD.values, valuesD, rtol=rtol)

    # With resource recovery
    br.INCLUDE_RESOURCE_RECOVERY = True

    modelA2 = br.create_model('A')
    dfA2 = modelA2.metrics_at_baseline()
    valuesA2 = [37.58, 64.69, 67.36, 10.41, 28.3, -0.02703, 10.11, -0.3066]
    assert_allclose(dfA2.values, valuesA2, rtol=rtol)
    
    modelB2 = br.create_model('B')
    dfB2 = modelB2.metrics_at_baseline()
    valuesB2 = [81.31, 47.37, 72.52, 26.14, 1.886, 0.03668, 0.5631, 1.452]
    assert_allclose(dfB2.values, valuesB2, rtol=rtol)

    modelC2 = br.create_model('C')
    dfC2 = modelC2.metrics_at_baseline()
    valuesC2 = [0.0, 0.0, 0.0, 11.864, 35.428, 0.049, 11.207, 0.23]
    assert_allclose(dfC2.values, valuesC2, rtol=rtol)

    # Same results with/without resource recovery
    modelD2 = br.create_model('D')
    dfD2 = modelD2.metrics_at_baseline()
    valuesD2 = [0.0, 0.0, 0.0, 10.017, 74.644, 0.125, 22.774, 0.525]
    assert_allclose(dfD2.values, valuesD2, rtol=rtol)


if __name__ == '__main__':
    test_biogenic_refinery()