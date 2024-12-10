# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_bsm1',)

def test_bsm1():
    from numpy.testing import assert_allclose as ac
    from numpy import array, vstack, absolute
    from exposan import bsm1
    bsm1_ss_concs = array([
        [30, 2.81, 1149., 82.1, 2552., 148., 449., 0.0043, 5.37, 7.92, 1.22, 5.28, 4.93*12],
        [30, 1.46, 1149., 76.4, 2553., 148., 450., 6.31e-5, 3.66, 8.34, 0.882, 5.03, 5.08*12],
        [30, 1.15, 1149., 64.9, 2557., 149., 450., 1.72, 6.54, 5.55, 0.829, 4.39, 4.67*12],
        [30, 0.995, 1149., 55.7, 2559., 150., 451., 2.43, 9.3, 2.97, 0.767, 3.88, 4.29*12],
        [30, 0.889, 1149., 49.3, 2559., 150., 452., 0.491, 10.4, 1.73, 0.688, 3.53, 4.13*12]
        ])
    bsm1_ss_c1tss = array([12.5, 18.1, 29.5, 69.0, 356., 356., 356., 356., 356., 6394])
    bsm1.load()
    rtol = 0.01
    sys = bsm1.sys
    sys.simulate(t_span=(0,100), method='BDF', state_reset_hook='reset_cache')
    assert sys.outs[0].isempty() == False
    ac(float(sys.outs[0].iconc['S_S']), 0.895, rtol=rtol)
    ac(float(sys.outs[1].iconc['X_BH']), 4994.3, rtol=rtol)
    ac(sys.outs[0].COD, 47.5, rtol=rtol)
    ac(sys.outs[1].get_TSS(), 6377.9, rtol=rtol)
    u = sys.flowsheet.unit
    concs = vstack((u.A1._state[:13], u.A2._state[:13], u.O1._state[:13], u.O2._state[:13], u.O3._state[:13]))
    ac(concs, bsm1_ss_concs, rtol=rtol)
    ac(u.C1._state[-10:], bsm1_ss_c1tss, rtol=rtol)
    
    asm2d_ss = {'A1': {'S_O2': 0.0041,
      'S_F': 2.9275,
      'S_A': 4.9273,
      'S_I': 30.0,
      'S_NH4': 21.483,
      'S_N2': 20.2931,
      'S_NO3': 0.2331,
      'S_PO4': 10.3835,
      'S_ALK': 71.7624,
      'X_I': 1686.8,
      'X_S': 141.1854,
      'X_H': 1846.2,
      'X_PAO': 210.1226,
      'X_PP': 60.6935,
      'X_PHA': 6.4832,
      'X_AUT': 115.4611,
      'X_MeOH': 0.0,
      'X_MeP': 0.0},
     'A2': {'S_O2': 7.5729e-06,
      'S_F': 1.1289,
      'S_A': 15.3198,
      'S_I': 30.0,
      'S_NH4': 22.7692,
      'S_N2': 20.5221,
      'S_NO3': 0.0043,
      'S_PO4': 15.12,
      'S_ALK': 70.05839999999999,
      'X_I': 1688.3,
      'X_S': 133.9323,
      'X_H': 1833.3,
      'X_PAO': 209.3625,
      'X_PP': 56.3103,
      'X_PHA': 16.8754,
      'X_AUT': 115.1802,
      'X_MeOH': 0.0,
      'X_MeP': 0.0},
     'A3': {'S_O2': 0.0327,
      'S_F': 0.7436,
      'S_A': 0.5607,
      'S_I': 30.0,
      'S_NH4': 10.6314,
      'S_N2': 27.5147,
      'S_NO3': 5.4225,
      'S_PO4': 9.6308,
      'S_ALK': 58.849199999999996,
      'X_I': 1694.4,
      'X_S': 81.1308,
      'X_H': 1852.3,
      'X_PAO': 213.2479,
      'X_PP': 61.7766,
      'X_PHA': 6.159,
      'X_AUT': 117.0121,
      'X_MeOH': 0.0,
      'X_MeP': 0.0},
     'A4': {'S_O2': 0.00061765,
      'S_F': 0.5047,
      'S_A': 0.1455,
      'S_I': 30.0,
      'S_NH4': 10.7277,
      'S_N2': 28.8822,
      'S_NO3': 4.0583,
      'S_PO4': 9.093,
      'S_ALK': 60.2748,
      'X_I': 1695.0,
      'X_S': 78.092,
      'X_H': 1852.3,
      'X_PAO': 213.628,
      'X_PP': 62.3363,
      'X_PHA': 5.0784,
      'X_AUT': 116.8988,
      'X_MeOH': 0.0,
      'X_MeP': 0.0},
     'O1': {'S_O2': 1.9707,
      'S_F': 0.4748,
      'S_A': 0.0336,
      'S_I': 30.0,
      'S_NH4': 8.0209,
      'S_N2': 29.0603,
      'S_NO3': 6.6395,
      'S_PO4': 7.8953,
      'S_ALK': 55.995599999999996,
      'X_I': 1695.8,
      'X_S': 68.2975,
      'X_H': 1855.5,
      'X_PAO': 214.5319,
      'X_PP': 63.5316,
      'X_PHA': 2.7381,
      'X_AUT': 117.4083,
      'X_MeOH': 0.0,
      'X_MeP': 0.0},
     'O2': {'S_O2': 2.6022,
      'S_F': 0.4137,
      'S_A': 0.0133,
      'S_I': 30.0,
      'S_NH4': 5.3737,
      'S_N2': 29.1839,
      'S_NO3': 9.2563,
      'S_PO4': 7.1891,
      'S_ALK': 51.620400000000004,
      'X_I': 1696.6,
      'X_S': 59.943,
      'X_H': 1857.8,
      'X_PAO': 214.9675,
      'X_PP': 64.2491,
      'X_PHA': 1.2117,
      'X_AUT': 117.9127,
      'X_MeOH': 0.0,
      'X_MeP': 0.0},
     'O3': {'S_O2': 3.1517,
      'S_F': 0.3715,
      'S_A': 0.0071,
      'S_I': 30.0,
      'S_NH4': 3.0099,
      'S_N2': 29.2754,
      'S_NO3': 11.6625,
      'S_PO4': 6.9291,
      'S_ALK': 47.5728,
      'X_I': 1697.4,
      'X_S': 52.9409,
      'X_H': 1859.2,
      'X_PAO': 214.9657,
      'X_PP': 64.5339,
      'X_PHA': 0.4642,
      'X_AUT': 118.3582,
      'X_MeOH': 0.0,
      'X_MeP': 0.0}}
    eff_ss = {'S_O2': 3.1517,
              'S_F': 0.3715,
              'S_A': 0.0071,
              'S_I': 30.0,
              'S_NH4': 3.0099,
              'S_N2': 29.2754,
              'S_NO3': 11.6625,
              'S_PO4': 6.9291,
              'S_ALK': 47.5728,
              'X_I': 6.2584,
              'X_S': 0.1952,
              'X_H': 6.8551,
              'X_PAO': 0.7926,
              'X_PP': 0.2379,
              'X_PHA': 0.0017,
              'X_AUT': 0.4364,
              'X_MeOH': 0.0,
              'X_MeP': 0.0,}
    sludge_ss = {'S_O2': 3.1517,
                 'S_F': 0.3715,
                 'S_A': 0.0071,
                 'S_I': 30.0,
                 'S_NH4': 3.0099,
                 'S_N2': 29.2754,
                 'S_NO3': 11.6625,
                 'S_PO4': 6.9291,
                 'S_ALK': 47.5728,
                 'X_I': 3319.4,
                 'X_S': 103.5298,
                 'X_H': 3635.9,
                 'X_PAO': 420.3813,
                 'X_PP': 126.2008,
                 'X_PHA': 0.9077,
                 'X_AUT': 231.4581,
                 'X_MeOH': 0.0,
                 'X_MeP': 0.0,}    
    c1tss_ss = array([12.8852, 18.5343, 30.1775, 70.8981, 373.2459, 373.2459, 373.2459, 373.2459, 1774.9, 6834.2])

    bsm1.load(reload=True, suspended_growth_model='ASM2d', reactor_model='PFR')
    sys = bsm1.sys
    sys.simulate(t_span=(0,200), method='BDF', state_reset_hook='reset_cache')
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    cmps = s.effluent.components
    asm2d_ss = array([cmps.kwarray(v) for k,v in asm2d_ss.items()])
    
    val = u.AS.state.iloc[:,:-2]
    val[val.abs()<1e-6] = 0 # very small values may lead to inf %
    ac(val, asm2d_ss[:,:-1], rtol=rtol)
    
    val = absolute(s.effluent.state[:-2])
    val[val<1e-6] = 0
    ac(val, cmps.kwarray(eff_ss)[:-1], rtol=rtol)
    
    val = absolute(s.WAS.state[:-2])
    val[val<1e-6] = 0
    ac(val, cmps.kwarray(sludge_ss)[:-1], rtol=rtol)
    
    ac(u.C1._state[-10:], c1tss_ss, rtol=rtol)

if __name__ == '__main__':
    test_bsm1()