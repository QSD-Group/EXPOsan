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


def test_htl():
    from numpy.testing import assert_allclose
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

    # Baselines are for BioSTEAM >= 2.54.0 and thermo >= 0.6.1.
    m1 = htl.create_model('baseline', **kwargs)
    assert_allclose(m1.metrics_at_baseline().values,
                    [2.683049, -27.623159, 24.137262, 102.096472], rtol=rtol)

    m2 = htl.create_model('no_P', **kwargs)
    assert_allclose(m2.metrics_at_baseline().values,
                    [3.235468, 6.434865, 10.982505, -10.464207], rtol=rtol)

    m3 = htl.create_model('PSA', **kwargs)
    assert_allclose(m3.metrics_at_baseline().values,
                    [2.044986, -66.961316, 45.230590, 282.584720], rtol=rtol)


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
