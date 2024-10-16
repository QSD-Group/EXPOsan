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
    from numpy.testing import assert_allclose
    from exposan import adm
    adm.load()
    rtol = 1e-2
    t = 200
    sys = adm.sys
    AD = sys.flowsheet.unit.AD
    # AD.algebraic_h2 = True
    AD.algebraic_h2 = False
    sys.simulate(state_reset_hook='reset_cache',
                 t_span=(0, t),
                 method='BDF')
    ad_ss = np.array([
        0.0119548297170, 0.0053147401716, 0.0986214009308, 0.0116250064639, 
        0.0132507296663, 0.0157836662845, 0.1976297169375, 0.0000002359451,
        0.0550887764460, 0.1526778706263*12, 0.1302298158037*14, 0.3286976637215,
        0.3086976637215, 0.0279472404350, 0.1025741061067, 0.0294830497073, 
        0.4201659824546, 1.1791717989237, 0.2430353447194, 0.4319211056360, 
        0.1373059089340, 0.7605626583132, 0.3170229533613, 25.6173953274430,
        0.04, 0.02
        ])
    AD_state = AD._state[:26]
    # assert np.isclose(AD_state['S_su'], 0.01195482763310878, rtol=rtol)
    # assert np.isclose(AD_state['S_ch4'], 0.05512674690208366, rtol=rtol)
    # assert np.isclose(AD_state['X_aa'], 1.1791717409024813, rtol=rtol)
    # assert np.isclose(AD_state['X_ac'], 0.7605261727371091, rtol=rtol)
    assert_allclose(AD_state, ad_ss, rtol=rtol)

    bg = sys.flowsheet.stream.Biogas
    assert np.isclose(bg.imass['S_h2'], 0.0011943269556170858, rtol=rtol)
    assert np.isclose(bg.imass['S_ch4'], 189.70280013379164, rtol=rtol)
    assert np.isclose(bg.imass['S_IC'],  19.75185297943967, rtol=rtol)
    assert np.isclose(bg.imass['H2O'],  4.6151987, rtol=rtol)

if __name__ == '__main__':
    test_adm()
    
    