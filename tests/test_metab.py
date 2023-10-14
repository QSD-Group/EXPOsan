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
    from exposan.metab import create_system
    import qsdsan as qs, numpy as np
    
    # LCA objs between modules might interference each other
    from qsdsan.utils import clear_lca_registries
    clear_lca_registries()
    
    rtol = 1e-3
    qs.PowerUtility.price = 0.0913
    UASB_M = create_system(n_stages=2, reactor_type='UASB', gas_extraction='M', tot_HRT=4)
    UASB_M.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    fs = UASB_M.flowsheet.stream
    assert np.isclose(1 - fs.eff_dg.COD / fs.inf.COD, 0.9067142448236405, rtol)
    assert np.isclose(UASB_M.TEA.annualized_NPV, -123994.44966853957, rtol)
    assert np.isclose(UASB_M.LCA.total_impacts['GWP100'], 613761.7802056486, rtol)
    
    FB_H = create_system(n_stages=2, reactor_type='FB', gas_extraction='H', tot_HRT=4)
    FB_H.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    fs = FB_H.flowsheet.stream
    assert np.isclose(1 - fs.eff_dg.COD / fs.inf.COD, 0.8254341848530495, rtol)
    assert np.isclose(FB_H.TEA.annualized_NPV, -34968.07019357201, rtol)
    assert np.isclose(FB_H.LCA.total_impacts['GWP100'], 1004100.126558344, rtol)
    
    PB_P = create_system(n_stages=2, reactor_type='PB', gas_extraction='P', tot_HRT=4)
    # Might fail the first time it runs, re-running will usually fix the problem
    try: PB_P.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    except: PB_P.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
    fs = PB_P.flowsheet.stream
    assert np.isclose(1 - fs.eff_dg.COD / fs.inf.COD, 0.8261047568019022, rtol)
    assert np.isclose(PB_P.TEA.annualized_NPV, -50509.157739631584, rtol)
    assert np.isclose(PB_P.LCA.total_impacts['GWP100'], 178566.02455877443, rtol)
    
if __name__ == '__main__':
    test_metab()