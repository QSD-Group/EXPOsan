# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_metab', )

def test_metab():
    """
    See EXPOsan issue 61 and PR 65 on the metab metric discussion.
    https://github.com/QSD-Group/EXPOsan/pull/65
    
    Regression test for metab. Baselines re-measured 2026-05-26 against
    qsdsan ``ee855ac`` / EXPOsan ``cc21751`` (the merge of ``metab_fix``
    into ``main``), which include:

    * ``UASB._design``'s ``min_depth_to_diameter`` constraint that replaces
      the pancake geometry the original UASB+M sizing produced.
    * ``_F_vol`` in ``exposan.metab.utils.misc`` now uses ``cmps.chem_MW``
      instead of ``cmps.MW``, so it stays correct after qsdsan switched
      ``Component.MW`` to the measurement-basis molar mass (qsdsan commit
      ``bccb8329`` "optional adjustment of MW", 2025-01-28).

    FB_H and PB_P numbers match the Feb-13-2024 pinned stack within 0.01 %.
    UASB_M numbers differ from Feb-2024 by ~27 % on NPV and ~6 % on GWP —
    that's the intentional consequence of the new ``min_depth_to_diameter``
    sizing constraint replacing the pancake geometry the Feb-2024 code
    produced for UASB+M's 400x sidestream recirculation (~46 m diameter,
    ~223 t of stainless steel) with a buildable 2.77 m x 2.77 m reactor
    (~800 kg). Re-baseline after any qsdsan or biosteam version bump that
    touches TEA or LCA accounting.
    """
    from exposan.metab import create_system
    import qsdsan as qs, numpy as np

    rtol = 1e-3
    qs.PowerUtility.price = 0.0913
    UASB_M = create_system(n_stages=2, reactor_type='UASB', gas_extraction='M', tot_HRT=4)
    UASB_M.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    fs = UASB_M.flowsheet.stream

    assert np.isclose(1 - fs.eff_dg.COD / fs.inf.COD, 0.9067142448236405, rtol)
    assert np.isclose(UASB_M.TEA.annualized_NPV, -72844.68, rtol=rtol)
    assert np.isclose(UASB_M.LCA.total_impacts['GWP100'], 3720831.94, rtol=rtol)

    FB_H = create_system(n_stages=2, reactor_type='FB', gas_extraction='H', tot_HRT=4)
    # Might fail the first time it runs, re-running will usually fix the problem
    try: FB_H.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    except: FB_H.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    fs = FB_H.flowsheet.stream
    assert np.isclose(1 - fs.eff_dg.COD / fs.inf.COD, 0.8254350623696006, rtol)
    assert np.isclose(FB_H.TEA.annualized_NPV, -167402.83, rtol=rtol)
    assert np.isclose(FB_H.LCA.total_impacts['GWP100'], 1235954.41, rtol=rtol)

    PB_P = create_system(n_stages=2, reactor_type='PB', gas_extraction='P', tot_HRT=4)
    # Might fail the first time it runs, re-running will usually fix the problem
    try: PB_P.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    except: PB_P.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    fs = PB_P.flowsheet.stream
    assert np.isclose(1 - fs.eff_dg.COD / fs.inf.COD, 0.8261060899768736, rtol)
    assert np.isclose(PB_P.TEA.annualized_NPV, -352553.41, rtol=rtol)
    assert np.isclose(PB_P.LCA.total_impacts['GWP100'], 1535648.71, rtol=rtol)

if __name__ == '__main__':
    test_metab()