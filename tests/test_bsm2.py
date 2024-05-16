# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_bsm2',)


def test_bsm2():
    from numpy.testing import assert_allclose as ac
    import numpy as np
    from exposan import bsm2
    bsm2.load()
    rtol = 5e-2
    # rtol = 1e-2
    sys = bsm2.sys
       
    t_span = (0, 15) # just 15 days to make the test faster
    sys.simulate(method='RK23', t_span=t_span,
                  # state_reset_hook='reset_cache')
                )
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    
    cmps_asm1 = s.inf.components
    eff_ss_concs = dict(
        S_I = 28.0643,
        S_S = 0.67336,
        X_I = 5.9191,
        X_S = 0.12329,
        X_BH = 8.6614,
        X_BA = 0.6484,
        X_P = 3.7485,
        S_O = 1.3748,
        S_NO = 9.1948,
        S_NH = 0.15845,
        S_ND = 0.55943,
        X_ND = 0.0092428,
        S_ALK = 4.5646*12
        )
    ac(s.effluent.state[:13], cmps_asm1.kwarray(eff_ss_concs)[:13], rtol=rtol)
    assert np.isclose(s.effluent.get_TSS(), 14.3255, rtol=rtol)
    assert np.isclose(s.effluent.F_vol*24, 20640.7791, rtol=rtol)

    sludge_disposal = dict(
        S_I = 130.867,
        S_S = 258.5789,
        X_I = 314239.8855,
        X_S = 47666.1788,
        X_P = 11427.269,
        S_NH = 1442.7882,
        S_ND = 0.54323,
        X_ND = 1841.0745,
        S_ALK = 97.8459*12
        )
    ac(s.digested_sludge.state[:13], cmps_asm1.kwarray(sludge_disposal)[:13], rtol=rtol)
    assert np.isclose(s.digested_sludge.get_TSS(), 280000, rtol=rtol)
    assert np.isclose(s.digested_sludge.F_vol*24, 9.5821, rtol=rtol)

    as_ss_concs = np.array([
        [28.0643, 3.0503, 1532.2609, 63.0433, 2245.0634, 166.6699, 964.8992, 0.0093422, 3.935, 6.8924, 0.95797, 3.8453, 5.4213*12],
        [28.0643, 1.3412, 1532.2609, 58.8579, 2245.3852, 166.5512, 965.6805, 0.00010907, 2.2207, 7.2028, 0.68624, 3.7424, 5.5659*12],
        [28.0643, 0.95531, 1532.2609, 46.2983, 2246.7994, 167.3077, 967.2442, 0.46635, 5.5141, 3.4247, 0.65129, 3.1405, 5.0608*12],
        [28.0643, 0.78055, 1532.2609, 37.3881, 2245.6315, 167.8339, 968.8072, 1.4284, 8.4066, 0.69216, 0.60938, 2.6815, 4.659*12],
        [28.0643, 0.67336, 1532.2609, 31.9144, 2242.1274, 167.8482, 970.3678, 1.3748, 9.1948, 0.15845, 0.55943, 2.3926, 4.5646*12]
        ])
    concs = np.vstack((u.A1._state[:13], u.A2._state[:13], u.O1._state[:13], u.O2._state[:13], u.O3._state[:13]))
    ac(concs, as_ss_concs, rtol=rtol)

    clarifier_tss = np.array([14.3255, 20.8756, 34.2948, 81.0276, 
                              423.2035, 423.2035, 423.2035, 423.2035, 
                              3710.5517, 7348.2757])
    ac(u.C2._state[-10:], clarifier_tss, rtol=rtol)
    
    cmps_adm1 = s.biogas.components
    h2 = cmps_adm1.S_h2
    ch4 = cmps_adm1.S_ch4
    co2 = cmps_adm1.S_IC    
    ad_inf_adm = dict(
        S_aa = 0.04388,
        S_IC = 0.0079326*12,
        S_IN = 0.0019721*14,
        S_I = 0.028067,
        X_ch = 3.7236,
        X_pr = 15.9235,
        X_li = 8.047,
        X_I = 17.0106,
        S_an = 0.0052101
        )
    ac(u.AD1.ins[0].state[:26], cmps_adm1.kwarray(ad_inf_adm)[:26]*1e3, rtol=rtol)

    ad_ss = dict(
        S_su = 0.012394,
        S_aa = 0.0055432,
        S_fa = 0.10741,
        S_va = 0.012333,
        S_bu = 0.014003,
        S_pro = 0.017584,
        S_ac = 0.089315,
        S_h2 = 2.5055e-07,
        S_ch4 = 0.05549,
        S_IC = 0.095149*12,
        S_IN = 1.3226,
        S_I = 0.13087,
        X_c = 0.10792,
        X_ch = 0.020517,
        X_pr = 0.08422,
        X_li = 0.043629,
        X_su = 0.31222,
        X_aa = 0.93167,
        X_fa = 0.33839,
        X_c4 = 0.33577,
        X_pro = 0.10112,
        X_ac = 0.67724,
        X_h2 = 0.28484,
        X_I = 17.2162,
        S_cat = 0.,
        S_an = 0.0052101,
        S_h2_gas = 1.1032e-05 * h2.i_mass / h2.chem_MW,
        S_ch4_gas = 1.6535 * ch4.i_mass / ch4.chem_MW,
        S_IC_gas = 0.01354,
        Q = 178.4674,
        )
    ad_ss = np.array([*ad_ss.values()])
    ad_state = np.delete(u.AD1._state, 26)  # remove "H2O"
    ac(ad_state, ad_ss, rtol=rtol)
    assert np.isclose(s.AD_eff.pH, 7.2631, rtol=rtol)
    assert np.isclose(s.biogas.imass['S_h2']*24 * h2.i_mass, 0.0035541, rtol=rtol)
    assert np.isclose(s.biogas.imass['S_ch4']*24 * ch4.i_mass, 1065.3523, rtol=rtol)
    assert np.isclose(s.biogas.imass['S_IC']*24 * co2.i_mass, 1535.4118, rtol=rtol)
    
    ad_eff_asm = dict(
        S_I = 130.867,
        S_S = 258.5789,
        X_I = 17216.2434,
        X_S = 2611.4843,
        X_P = 626.0652,
        S_NH = 1442.7882,
        S_ND = 0.54323,
        X_ND = 100.8668,
        S_ALK = 97.8459*12,
        )
    ac(u.J2._state[:13], cmps_asm1.kwarray(ad_eff_asm)[:13], rtol=rtol)


if __name__ == '__main__':
    test_bsm2()