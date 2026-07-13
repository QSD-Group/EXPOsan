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

__all__ = ('test_htl', 'test_htl_reversed_splitter_recycle')


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

    ``thermo`` (a thermosteam dependency, unpinned in EXPOsan's own
    requirements) released 0.6.1 with a change to ``TDependentProperty``'s
    default liquid-volume extrapolation. Because the
    HTL natural-gas/air streams are built from light-gas components (S_O2,
    S_N2, S_CH4, S_H2) whose "liquid volume" at process conditions is deep
    extrapolation territory, this shifts the ``no_P`` model's metrics by more
    than the existing pre/post-fix baselines account for. Until EXPOsan pins
    ``thermo``, accept the thermo-0.6.1 variant of each baseline too.
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
    # thermo 0.6.1 shifts this model enough that the pre/post-fix baselines
    # above no longer match either; the two extra entries are the same
    # pre-fix/post-fix pair re-measured under thermo 0.6.1. See
    # ``_assert_allclose_any`` for details.
    _assert_allclose_any(m2.metrics_at_baseline().values,
                         [[3.251, 7.419, 11.932, -2.338],
                          [3.219, 5.408, 11.148, -9.046],
                          [3.270026, 8.565268, 11.938782, -2.28163],
                          [3.237357, 6.551209, 11.154833, -8.989472]], rtol)

    m3 = htl.create_model('PSA', **kwargs)
    _assert_allclose_any(m3.metrics_at_baseline().values,
                         [[2.058, -66.152, 46.083, 289.873],
                          [2.027, -68.043, 46.295, 291.694]], rtol)


def test_htl_reversed_splitter_recycle():
    '''
    SP1/RSP1 (``ReversedSplitter``, H2SO4/H2 makeup) compute their split from
    AcidEx/MemDis/HT/HC's demand, but those units run *after* them in the
    network path -- so without an explicit recycle, resimulating the same
    ``System`` at a new ``plant_size`` reports the previous evaluation's
    demand, not the current one. Regression test for that: evaluating the
    same plant_size twice, with a very different plant_size evaluated in
    between, must give the same result both times.
    '''
    from numpy.testing import assert_allclose
    from chaospy import distributions as shape
    from exposan.htl import create_model

    model = create_model(
        plant_size=True, feedstock='sludge', include_CFs_as_metrics=False,
        include_other_metrics=False, include_other_CFs_as_metrics=False,
        )
    plant_size = model.parameters[-1]
    MDSP, GWP = [m for m in model.metrics if m.name in ('MDSP', 'GWP diesel')]

    plant_size.baseline = 150
    model.metrics_at_baseline()
    first = (MDSP.get(), GWP.get())

    plant_size.baseline = 800 # a very different scale in between
    model.metrics_at_baseline()

    plant_size.baseline = 150 # back to the original scale
    model.metrics_at_baseline()
    second = (MDSP.get(), GWP.get())

    assert_allclose(second, first, rtol=1e-3)


if __name__ == '__main__':
    test_htl()
    test_htl_reversed_splitter_recycle()
