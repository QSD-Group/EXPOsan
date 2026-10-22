# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_metro',)

def test_metro():
    import os
    from numpy.testing import assert_allclose as ac
    from exposan.metro import metroASM2d as m

    # Baselines are for 5 d of simulation (explicit RK23) from the `ss` initial conditions,
    # i.e., the system has not reached steady state (phosphorus is still changing, the
    # 15 d simulation takes about 6 min). They are the same (within 1e-8) for BioSTEAM 2.54.0
    # and later versions.
    rtol, atol = 0.01, 1e-4
    comps = ('S_O2', 'S_F', 'S_A', 'S_NH4', 'S_NO3', 'S_PO4', 'S_ALK', 'X_I', 'X_S', 'X_H', 'X_PAO', 'X_PP', 'X_PHA', 'X_AUT')
    tank_concs = { # mg/L
        'O1': [2, 0.27932, 0.0051016, 0.11011, 22.142, 3.1809, 56.13, 2020.2, 49.324, 3024.7, 205.05, 60.043, 0.0071292, 169.08],
        'A1': [0.0028281, 0.21198, 0.027618, 5.4619, 9.8802, 5.4219, 70.338, 2035.4, 65.927, 2971, 202.73, 58.988, 0.036783, 167.12],
        'A2': [5.0606e-05, 3.2953, 66.805, 23.529, 0.081009, 5.2521, 82.238, 588.71, 111.32, 853.68, 57.786, 15.124, 4.3539, 47.541],
        'A3': [1.1211e-07, 1.3528, 79.852, 24.831, 0.0017397, 7.3465, 80.475, 589.93, 104.42, 844.8, 57.584, 13.267, 8.775, 47.377],
        'O2': [2, 2.4086, 3.9941, 15.871, 3.1177, 5.5252, 85.12, 591.23, 79.308, 902.01, 58.855, 14.029, 8.409, 48.181],
        'O3': [2, 0.82746, 0.2616, 11.376, 6.7771, 3.8671, 79.212, 592.53, 59.414, 913.78, 60.114, 15.525, 5.9903, 48.971],
        'O4': [2, 0.51184, 0.045045, 7.388, 10.384, 2.4094, 73.05, 593.85, 44.546, 919.54, 61.296, 16.91, 3.4919, 49.725],
        'O5': [2, 0.3951, 0.016318, 3.9048, 13.681, 1.3651, 67.455, 595.17, 33.995, 922.36, 62.23, 17.925, 1.4362, 50.39],
        'O6': [2, 0.32859, 0.0090179, 1.4657, 16.129, 0.92788, 63.345, 596.49, 26.819, 923.04, 62.611, 18.373, 0.37077, 50.841],
        'O7': [2, 0.30012, 0.0067342, 0.43699, 17.328, 0.90693, 61.421, 597.81, 22.09, 922.17, 62.556, 18.437, 0.078495, 50.984],
        'O8': [2, 0.30748, 0.0063898, 0.17135, 17.849, 1.0228, 60.698, 599.13, 19.037, 920.23, 62.371, 18.383, 0.018475, 50.96],
        'O9': [2, 0.31066, 0.0064091, 0.11751, 18.201, 1.1758, 60.291, 600.45, 17.09, 917.62, 62.161, 18.301, 0.0059533, 50.893],
    }
    C1_TSS = [17.047, 26.769, 42.685, 83.828, 267.2, 267.2, 267.21, 267.22, 267.22, 4840.8] # mg/L, layers of the secondary clarifier
    effluent = dict(
        Q=6.171e+05,
        COD=33.698,
        BOD5=7.2033,
        TN=19.847,
        TP=1.593,
        TSS=16.835,
        NH4=0.11684,
        NO3=18.22,
        PO4=1.0976)
    streams = dict(
        WAS_Q=11356,
        RAS_Q=2.4984e+05,
        sludge_PC_Q=2000,
        sludge_PC_TSS=69793,
        sludge_DU1_Q=675.28,
        sludge_DU1_TSS=2.6859e+05,
        effluent_DU1_Q=2581.4,
        effluent_TC2_Q=10099)
    incinerator = dict( # kg/hr and kJ/hr
        CO2=14695,
        H2O=34789,
        N2=62454,
        O2=4986.1,
        ash=446.22,
        heat_combustion=1.7801e+08)

    sys = m.create_system()
    m.batch_init(sys, os.path.join(os.path.dirname(m.__file__), 'initial_conditions_ASM2d.xlsx'), sheet='ss')
    sys.simulate(t_span=(0, 5), method='RK23', state_reset_hook='reset_cache')
    u, fs = sys.flowsheet.unit, sys.flowsheet.stream
    cmps = u.O1.components

    for ID, concs in tank_concs.items():
        ac([u[ID]._state[cmps.index(i)] for i in comps], concs, rtol=rtol, atol=atol, err_msg=ID)
    ac(u.C1._state[-10:], C1_TSS, rtol=rtol, atol=atol)

    def check(actual, expected):
        for k, v in expected.items():
            ac(actual[k], v, rtol=rtol, atol=atol, err_msg=k)

    e = fs.effluent
    check(dict(Q=e.F_vol*24, COD=e.COD, BOD5=e.BOD, TN=e.TN, TP=e.TP, TSS=e.get_TSS(),
               NH4=e.iconc['S_NH4'], NO3=e.iconc['S_NO3'], PO4=e.iconc['S_PO4']), effluent)
    check(dict(WAS_Q=fs.WAS.F_vol*24, RAS_Q=fs.RAS.F_vol*24,
               sludge_PC_Q=fs.sludge_PC.F_vol*24, sludge_PC_TSS=fs.sludge_PC.get_TSS(),
               sludge_DU1_Q=fs.sludge_DU1.F_vol*24, sludge_DU1_TSS=fs.sludge_DU1.get_TSS(),
               effluent_DU1_Q=fs.effluent_DU1.F_vol*24, effluent_TC2_Q=fs.effluent_TC2.F_vol*24), streams)

    flue_gas, ash = u.IC.outs
    check(dict(CO2=flue_gas.imass['S_CO2'], H2O=flue_gas.imass['H2O'], N2=flue_gas.imass['S_N2'],
               O2=flue_gas.imass['S_O2'], ash=ash.imass['ash'], heat_combustion=u.IC.heat_combustion), incinerator)

if __name__ == '__main__':
    test_metro()
