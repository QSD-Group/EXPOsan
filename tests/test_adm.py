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

# TODO: change to the new loading routine
def test_adm():
    import numpy as np
    from exposan import adm
    # adm.load()

    t = 200
    t_step = 5
    sys = adm.sys
    sys.simulate(state_reset_hook='reset_cache',
             t_span=(0, t),
             t_eval=np.arange(0, t+t_step, t_step),
             method='BDF')

    AD = sys.flowsheet.unit.AD
    AD_state = AD.state
    assert np.isclose(AD_state['S_su'], 0.012037175687149777, rtol=0.01)
    assert np.isclose(AD_state['S_ch4'], 0.04139786593982221, rtol=0.01)
    assert np.isclose(AD_state['X_c'], 0.033989158977463685, rtol=0.01)
    assert np.isclose(AD_state['X_ac'], 0.15391526345651502, rtol=0.01)

    bg = sys.flowsheet.stream.Biogas
    assert np.isclose(bg.imass['S_h2'], 0.000362, rtol=0.01)
    assert np.isclose(bg.imass['S_ch4'], 40.3, rtol=0.01)
    assert np.isclose(bg.imass['S_IC'],  6.1, rtol=0.01)
    assert np.isclose(bg.imass['H2O'],  1.19, rtol=0.01)

if __name__ == '__main__':
    test_adm()