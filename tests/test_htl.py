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

__all__ = ('test_htl', 'test_geospatial_cooling_tower_inlets')


def _assert_allclose_any(actual, options, rtol):
    '''
    Assert ``actual`` matches at least one of the accepted baseline vectors in
    ``options`` (within ``rtol``).

    BioSTEAM corrected the natural-gas heat-utility air composition from mass to
    molar units (``reset_flow(O2=21, N2=79, units='kg/hr')`` -> ``'kmol/hr'``)
    after 2.53.11. Air is 21:79 O2:N2 by mole, so the fix is a correction, and it
    shifts the HTL TEA metrics (most notably the sludge management price). Until
    EXPOsan can raise its BioSTEAM floor past that fix, accept either the pre-fix
    (released <= 2.53.11) or the post-fix value so the test passes on both. Once
    the floor is raised, drop the pre-fix option and keep only the corrected one.
    '''
    from numpy.testing import assert_allclose
    last = None
    for option in options:
        try:
            assert_allclose(actual, option, rtol=rtol)
            return
        except AssertionError as e:
            last = e
    raise AssertionError(
        f'values matched none of the {len(options)} accepted baselines '
        f'(rtol={rtol}). Last comparison:\n{last}'
    )


def test_htl():
    from exposan import htl

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

    # Each entry: [pre-fix baseline (BioSTEAM <= 2.53.11), post-fix baseline
    # (BioSTEAM with the natural-gas air O2-content correction)]. See
    # ``_assert_allclose_any`` for why both are accepted.
    m1 = htl.create_model('baseline', **kwargs)
    _assert_allclose_any(m1.metrics_at_baseline().values,
                         [[2.698, -26.725, 25.087, 110.220],
                          [2.667, -28.606, 25.162, 110.859]], rtol)

    m2 = htl.create_model('no_P', **kwargs)
    _assert_allclose_any(m2.metrics_at_baseline().values,
                         [[3.251, 7.419, 11.932, -2.338],
                          [3.219, 5.408, 11.148, -9.046]], rtol)

    m3 = htl.create_model('PSA', **kwargs)
    _assert_allclose_any(m3.metrics_at_baseline().values,
                         [[2.058, -66.152, 46.083, 289.873],
                          [2.027, -68.043, 46.295, 291.694]], rtol)


def test_geospatial_cooling_tower_inlets():
    from exposan.htl.geospatial_systems import _create_cooling_tower

    _, makeup_water, chemicals = _create_cooling_tower()

    assert makeup_water.ID == 'cooling_tower_makeup_water'
    assert chemicals.ID == 'cooling_tower_chemicals'
    assert makeup_water.price > 0
    assert chemicals.price > makeup_water.price


if __name__ == '__main__':
    test_htl()
