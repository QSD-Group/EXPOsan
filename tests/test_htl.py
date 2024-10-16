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
    
    rtol = 5e-2
    kwargs = dict(
        feedstock='sludge',
        plant_size=False,
        ternary=False,
        high_IRR=False,
        exclude_sludge_compositions=False,
        include_HTL_yield_as_metrics=False,
        include_other_metrics=False,
        include_other_CFs_as_metrics=False,
        include_check=False,
        )
    
    m1 = htl.create_model('baseline', **kwargs)
    df1 = m1.metrics_at_baseline()
    values1 = [3.218, 5.355, 60.279, 411.336]
    assert_allclose(df1.values, values1, rtol=rtol)
    
    m2 = htl.create_model('no_P', **kwargs)
    df2 = m2.metrics_at_baseline()
    values2 = [3.763, 38.939, 46.453, 293.038]
    assert_allclose(df2.values, values2, rtol=rtol)

    m3 = htl.create_model('PSA', **kwargs)
    df3 = m3.metrics_at_baseline()
    values3 = [2.577, -34.192, 69.927, 493.899]
    assert_allclose(df3.values, values3, rtol=rtol)


if __name__ == '__main__':
    test_htl()