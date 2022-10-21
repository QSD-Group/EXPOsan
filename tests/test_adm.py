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

__all__ = ('test_adm',)

def test_adm():
    import numpy as np
    from exposan import adm
    adm.load()
    
    rtol = 1e-4
    t = 200
    sys = adm.sys
    sys.simulate(state_reset_hook='reset_cache',
                 t_span=(0, t),
                 method='BDF')

    AD = sys.flowsheet.unit.AD
    AD_state = AD.state
    assert np.isclose(AD_state['S_su'], 0.01195482763310878, rtol=rtol)
    assert np.isclose(AD_state['S_ch4'], 0.05512674690208366, rtol=rtol)
    assert np.isclose(AD_state['X_aa'], 1.1791717409024813, rtol=rtol)
    assert np.isclose(AD_state['X_ac'], 0.7605261727371091, rtol=rtol)

    bg = sys.flowsheet.stream.Biogas
    assert np.isclose(bg.imass['S_h2'], 0.0011943269556170858, rtol=rtol)
    assert np.isclose(bg.imass['S_ch4'], 189.70280013379164, rtol=rtol)
    assert np.isclose(bg.imass['S_IC'],  19.75185297943967, rtol=rtol)
    assert np.isclose(bg.imass['H2O'],  4.6151987, rtol=rtol)

if __name__ == '__main__':
    test_adm()