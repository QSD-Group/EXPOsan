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

__all__ = ('test_saf',)

def test_saf():
    from numpy.testing import assert_allclose
    from exposan import saf

    # The crude distillation columns use the FUG shortcut method, whose convergence
    # depends on the relative volatility of the keys as computed by the thermo
    # package. A deterministic Lr-ladder specification (see saf/systems.py) keeps the
    # achieved split near the design intent, but the absolute MFSP/GWP can still drift
    # with the thermo stack (and, for heat-duty driven terms, across platforms) so the
    # expected values below are loose targets re-baselined against the current stack
    # (biosteam 2.53.x). Update them if a deliberate change shifts them beyond tolerance.
    rtol = 0.15

    # The EC config GWP (~1.2 kg CO2e/GGE) is a near-zero balance of much larger opposing
    # terms (feedstock credit ~-12 vs natural gas ~+1.1 per GGE), so the small per-platform
    # drift in natural-gas heat demand (~0.4 kg CO2e/GGE between Windows and the Linux CI)
    # swamps a relative tolerance. Pair rtol with an absolute tolerance, sized to a few
    # percent of the largest GWP magnitude in the set, so the near-zero case stays robust
    # while the large-magnitude GWPs keep an effectively relative check.
    gwp_atol = 0.75

    # Re-baselined 2026-07-15 after fixing 7 confirmed bugs, see relevant commits for details
    saf.load(configuration='baseline', simulate=True)
    assert_allclose(saf.get_MFSP(saf.sys), 3.636344, rtol=rtol)
    assert_allclose(saf.get_GWP(saf.sys), -7.415470, rtol=rtol, atol=gwp_atol)

    saf.load(configuration='EC', simulate=True)
    assert_allclose(saf.get_MFSP(saf.sys), 11.529264, rtol=rtol)
    assert_allclose(saf.get_GWP(saf.sys), 0.785847, rtol=rtol, atol=gwp_atol)

    saf.load(configuration='EC-Future', simulate=True)
    assert_allclose(saf.get_MFSP(saf.sys), 3.495850, rtol=rtol)
    assert_allclose(saf.get_GWP(saf.sys), -10.181651, rtol=rtol, atol=gwp_atol)

    #!!! TODO: below were previous metrics, need to figure out the GWP change reason.
    # saf.load(configuration='baseline')
    # assert_allclose(saf.get_MFSP(saf.sys), 3.95586679600505, rtol=rtol)
    # assert_allclose(saf.get_GWP(saf.sys), -5.394022805849971, rtol=rtol)

    # saf.load(configuration='EC')
    # assert_allclose(saf.get_MFSP(saf.sys), 11.876241988677974, rtol=rtol)
    # assert_allclose(saf.get_GWP(saf.sys), 2.8357334832704386, rtol=rtol)

    # saf.load(configuration='EC-Future')
    # assert_allclose(saf.get_MFSP(saf.sys), 3.821113328378629, rtol=rtol)
    # assert_allclose(saf.get_GWP(saf.sys), -8.475883955624251, rtol=rtol)


if __name__ == '__main__':
    test_saf()
    from exposan import saf