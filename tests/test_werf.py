# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

def load_mdl(ID):
    import os
    from numpy import load
    from exposan.werf import create_system, baseline_underflows, \
        add_performance_metrics, add_OPEX_metrics, data_path
    from exposan.werf.utils import load_state
    from qsdsan import Model
    sys = create_system(ID)
    u = sys.flowsheet.unit
    
    if 'MT' in u: thickener = u.MT
    else: thickener = u.GT
    thickener.sludge_flow_rate, u.DW.sludge_flow_rate = baseline_underflows[ID]
    
    state_arr = load(os.path.join(data_path, f'ss/baseline_unopt/{sys.ID}.npy'))
    load_state(sys, state_arr=state_arr)

    mdl = Model(sys)
    add_performance_metrics(mdl)
    add_OPEX_metrics(mdl)
    
    return mdl
    
# %%

def test_werf():
    import numpy as np
    from numpy.testing import assert_allclose as ac
    from qsdsan import PowerUtility

    rtol = 1e-3
    atol = 1e-6
    sim_kwargs = dict(t_span=(0,300), method='BDF')
    PowerUtility.price = 0.0782

    b1 = load_mdl('B1')
    try: b1.system.simulate(**sim_kwargs)
    except: b1.system.simulate(**sim_kwargs, state_reset_hook='reset_cache')
    s = b1.system.flowsheet.stream
    cmps = s.SE.components
    n_cmps = len(cmps)-1 # ignore H2O
    eff_ss_concs = dict(
        S_O2 = 2.0,
        S_N2 = 14.705720402898551,
        S_NH4 = 4.210823229298282,
        S_NO3 = 23.680363157704974,
        S_PO4 = 5.100695912463844,
        S_F = 0.2929010469348792,
        S_A = 0.008593462912639268,
        S_I = 17.89701276071496,
        S_IC = 8.13746567555504,
        S_K = 27.99532722346303,
        S_Mg = 49.8103924923869,
        X_I = 6.6670127467992035,
        X_S = 0.4083150425105007,
        X_H = 15.133064776319413,
        X_AUT = 0.8758301733140433,
        S_Ca = 139.47267909982796,
        X_struv = 0.014186836918775132,
        X_ACP = 0.010077701517909736,
        S_Na = 86.98548101575692,
        S_Cl = 424.929073927556,
        )
    ac(s.SE.conc[:n_cmps], cmps.kwarray(eff_ss_concs)[:n_cmps], rtol=rtol, atol=atol)
    bg_ss_mass = dict(
        S_h2 = 0.00026464719739745007,
        S_ch4 = 190.57491717665997,
        S_IC = 23.16455952632622,
        H2O = 4.896805020228398,
        )
    cmps_ad = s.biogas.components
    _n = len(cmps_ad)-1 # ignore H2O
    ac(s.biogas.mass[:_n], cmps_ad.kwarray(bg_ss_mass)[:_n], rtol=rtol, atol=atol)    
    b1.system.flowsheet.clear()
    del b1
    
    e2 = load_mdl('E2')
    e2.system.simulate(**sim_kwargs)
    u = e2.system.flowsheet.unit
    # concentration profiles in activated sludge reactor
    nh4_ss = np.array([7.508173609332913,
                       2.570889381334568,
                       0.5159642619233105,
                       0.13710988536998334,
                       0.09047475387085772,
                       0.08854432510529574])
    ac(u.ASR.state.S_NH4.to_numpy(), nh4_ss, rtol=rtol, atol=atol)
    no3_ss = np.array([16.960594609868963,
                       21.48998773419129,
                       23.54371749088413,
                       24.141086380512917,
                       24.511183512675903,
                       24.88616005106289])
    ac(u.ASR.state.S_NO3.to_numpy(), no3_ss, rtol=rtol, atol=atol)
    aut_ss = np.array([42.913022475043704,
                       43.82994395362909,
                       44.124560127221464,
                       44.055433599505065,
                       43.92452808163498,
                       43.79134682732497])
    ac(u.ASR.state.X_AUT.to_numpy(), aut_ss, rtol=rtol, atol=atol)
    # clarifier TSS profile
    tss = np.array([7.669661492339992,
                    11.338648351775584,
                    16.43779477493079,
                    26.957626560395187,
                    60.61517109271088,
                    262.5204729258797,
                    262.52319696855886,
                    262.5231855786824,
                    262.5231854948057,
                    4522.883792033054])
    ac(u.FC._state[-10:], tss, rtol=rtol, atol=atol)
    aed_ss_concs = dict(
        S_O2 = 1.0,
        S_N2 = 21.898950530449103,
        S_NH4 = 0.08811964003129587,
        S_NO3 = 1240.7722368396628,
        S_PO4 = 670.5021597040704,
        S_F = 0.3321961109262077,
        S_A = 0.0018688311252324798,
        S_I = 17.90007374029111,
        S_IC = 57.13922480965932,
        S_K = 27.999938035507668,
        S_Mg = 49.99934996811206,
        X_I = 8109.195245089123,
        X_S = 54.4656799580144,
        X_H = 2833.462443099531,
        X_AUT = 154.28239245105794,
        S_Ca = 38.10349739263795,
        X_CaCO3 = 3.410249675566548e-07,
        X_struv = 0.0003403964132059657,
        X_newb = 0.008509910702394275,
        X_ACP = 288.44106284587235,
        X_MgCO3 = 3.403973084732852e-07,
        X_AlOH = 3.304167934089965e-10,
        X_AlPO4 = 8.720656685385554e-07,
        X_FeOH = 4.258984383336616e-10,
        X_FePO4 = 8.2018172311499e-07,
        S_Na = 86.99963930990731,
        S_Cl = 424.9982967367339,
        )
    cmps = u.AED.outs[0].components
    ac(u.AED._state[:-13], cmps.kwarray(aed_ss_concs)[:-12], rtol=rtol, atol=atol)
    ac(u.AED._state[-13:-2], cmps.kwarray(aed_ss_concs)[-12:-1], rtol=5e-2, atol=atol)
    e2.system.flowsheet.clear()
    del e2
    
    from math import isclose, isnan
    
    h1 = load_mdl('H1')
    h1.system.simulate(**sim_kwargs)
    metrics_ss = [21.964757971468394, 2.78300193174046, 8.007124655120837,
                  7.541028943961043, 0.1264223633457642, 1.7588373988576718,
                  1.3451256754322871, 7.952053512713743, 10.707765002030902,
                  43.51524011275996, 59.083695567070684, 25.36464254464812,
                  350358.4514610429, np.nan, 292.900827831498, 292.900827831498,
                  67.34442096121754, 17.83573101535616, 1.4859344331495183,
                  67.46251823857952, 0, 0.38707533711568903, 0.35882473670880044,
                  0.5056013448620468, 14.339519197200003, 12.4918596,
                  549.7162736741556, 291.6173830665255, 50.35713172658497,
                  1713.816948635574, 1708.6678372150877, 0,
                  0, 1737.4780143083963, 6051.653588626324]
    for m, val in zip(h1.metrics, metrics_ss):
        if not isnan(val):
            assert isclose(m(), val, rel_tol=rtol, abs_tol=atol), \
                f"System H1's {m.name} should equal {val} not {m()}"
    h1.system.flowsheet.clear()
    del h1
                
    i3 = load_mdl('I3')
    i3.system.simulate(**sim_kwargs)
    metrics_ss = [27.134694762854124, 2.3923400803895634, 8.067051949396125,
                  4.407792093009145, 0.2528977058283219, 0.9957580441871827,
                  0.5337430428853317, 9.822474460391645, 12.7116403614637,
                  np.nan, np.nan, 32.8011474266346,
                  398235.2881234985, np.nan, 332.9260221258218, 332.9260221258218,
                  67.34442096121754, 28.768572979845928, 3.2087257952996944,
                  67.46251823857952, 0, 0, 0.571297748361826, 
                  0, 46.901254679999994, 0, 
                  624.8355583257423, 314.0928694454979, 88.02427478342399,
                  0, 0, 175.38995883385218,
                  12.440491184499706, 2246.87859872447, 3461.6617512974863]
    for m, val in zip(i3.metrics, metrics_ss):
        if not isnan(val):
            assert isclose(m(), val, rel_tol=rtol, abs_tol=atol), \
                f"System I3's {m.name} should equal {val} not {m()}"
    i3.system.flowsheet.clear()
    del i3
    
    n2 = load_mdl('N2')
    n2.system.simulate(**sim_kwargs, state_reset_hook='reset_cache')
    metrics_ss = [19.08888104713469, 0.5350746339774913, 0.850918927266344,
                  5.1963865008605366, 0.2786375777769878, 0.3827746543974357,
                  0.21021289924642703, 6.893598008850556, 6.0612436387383,
                  np.nan, np.nan, 17.74962722779321,
                  1274336.5442744389, 429437.4019113503, 1065.3495789737399, 1494.7259104759896,
                  67.34442096121754, 89.950024318106, 3.3745616645195047,
                  67.46251823857952, 22.48429392933781, 0, 0.9113430564608498,
                  0, 14.227092460799996, 0,
                  2805.3015887813376, 472.0661779573176, 26.701407130429438,
                  7971.24162156081, 0, 0,
                  0, 1215.8494651038347, 12491.160260533728]
    for m, val in zip(n2.metrics, metrics_ss):
        if not isnan(val):
            assert isclose(m(), val, rel_tol=rtol, abs_tol=atol), \
                f"System N2's {m.name} should equal {val} not {m()}"

# %%
if __name__ == '__main__':
    test_werf()