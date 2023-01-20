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
    from qsdsan.utils import clear_lca_registries
    from exposan import htl

    clear_lca_registries()
    m1 = htl.create_model('baseline', key_metrics_only=True)
    df1 = m1.metrics_at_baseline()
    values1 = [4.545, 181.4, 48872.2, 441.4]
    assert_allclose(df1.values, values1, rtol=1e-2)
    
    m2 = htl.create_model('no_P', key_metrics_only=True)
    df2 = m2.metrics_at_baseline()
    values2 = [4.870, 212.7, 34854.5, 273.6]
    assert_allclose(df2.values, values2, rtol=1e-2)

    m3 = htl.create_model('PSA', key_metrics_only=True)
    df3 = m3.metrics_at_baseline()
    values3 = [3.877, 117.4, 63819.9, 622.3]
    assert_allclose(df3.values, values3, rtol=1e-2)


if __name__ == '__main__':
    test_htl()