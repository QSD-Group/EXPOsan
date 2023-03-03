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
    
    kwargs = dict(
        include_HTL_yield_as_metrics=False,
        include_other_metrics=False,
        include_other_CFs_as_metrics=False,
        include_check=False,
        )
    
    # m1 = htl.create_model('baseline', **kwargs)
    # df1 = m1.metrics_at_baseline()
    # values1 = [5.687, 157.57, -25.312, -321.006]
    # assert_allclose(df1.values, values1, rtol=5e-2)
    
    # m2 = htl.create_model('no_P', **kwargs)
    # df2 = m2.metrics_at_baseline()
    # values2 = [6.181, 188.063, -36.281, -414.864]
    # assert_allclose(df2.values, values2, rtol=5e-2)

    m3 = htl.create_model('PSA', **kwargs)
    df3 = m3.metrics_at_baseline()
    values3 = [4.67, 94.854, -4.472, -142.698]
    assert_allclose(df3.values, values3, rtol=5e-2)


if __name__ == '__main__':
    test_htl()