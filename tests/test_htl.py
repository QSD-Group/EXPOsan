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

__all__ = ('test_htl',)

def test_htl():
    from numpy.testing import assert_allclose
    from exposan import htl

    # Because of different CF settings for ImpactItem with the same ID
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()

    # m1 = htl.create_model('baseline', key_metrics_only=True)
    # df1 = m1.metrics_at_baseline()
    # values1 = [5.117, 169.9, 49380.0, 440.1]
    # assert_allclose(df1.values, values1, rtol=1e-2)
    
    # m2 = htl.create_model('no_P', key_metrics_only=True)
    # df2 = m2.metrics_at_baseline()
    # values2 = [5.482, 201.2, 35210.0, 272.3]
    # assert_allclose(df2.values, values2, rtol=1e-2)

    # Test one should be sufficient (one system is about 1 min),
    # QSDsan checks for the baseline configuration
    m3 = htl.create_model('PSA', key_metrics_only=True)
    df3 = m3.metrics_at_baseline()
    values3 = [4.364, 105.9, 64470.0, 621.0]
    assert_allclose(df3.values, values3, rtol=1e-2)


if __name__ == '__main__':
    test_htl()