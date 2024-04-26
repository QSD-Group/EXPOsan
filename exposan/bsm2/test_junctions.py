# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

Reference:
    Alex, J.; Benedetti, L.; Copp, J. B.; Gernaey, K. V.; Jeppsson, U.;
    Nopens, I.; Pons, M. N.; Rosen, C.; Steyer, J. P.; Vanrolleghem, P. A.
    Benchmark Simulation Model No. 2 (BSM2).
    http://iwa-mia.org/wp-content/uploads/2022/09/TR3_BSM_TG_Tech_Report_no_3_BSM2_General_Description.pdf.

'''
#%%

def test_adm1junctions():
    
    import qsdsan as qs, numpy as np
    from numpy.testing import assert_allclose as ac
    from qsdsan import (
        processes as pc,
        sanunits as su,
        WasteStream,
        )
    from qsdsan.utils import ospath, load_data
    from exposan.bsm2 import data_path
    
    #%%
    matlab_preAD_adm = {
        'S_su': 0.0, # monosacharides (kg COD/m3)
        'S_aa': 0.04388, # amino acids (kg COD/m3)
        'S_fa': 0.0, # long chain fatty acids (LCFA) (kg COD/m3)
        'S_va': 0.0, # total valerate (kg COD/m3)
        'S_bu': 0.0, # total butyrate (kg COD/m3)
        'S_pro': 0.0, # total propionate (kg COD/m3)
        'S_ac': 0.0, # total acetate (kg COD/m3)
        'S_h2': 0.0, # hydrogen gas (kg COD/m3)
        'S_ch4': 0.0, # methane gas (kg COD/m3)
        'S_IC': 0.0079326*12, # inorganic carbon (kmole C/m3 -> kg C/m3) 0.0951912
        'S_IN': 0.0019721*14, # inorganic nitrogen (kmole N/m3 -> kg N/m3) 0.0276094
        'S_I': 0.028067, # soluble inerts (kg COD/m3)
        'X_c': 0.0, # composites (kg COD/m3)
        'X_ch': 3.7236, # carbohydrates (kg COD/m3)
        'X_pr': 15.9235, # proteins (kg COD/m3)
        'X_li': 8.047, # lipids (kg COD/m3)
        'X_su': 0.0, # sugar degraders (kg COD/m3)
        'X_aa': 0.0, # amino acid degraders (kg COD/m3)
        'X_fa': 0.0, # LCFA degraders (kg COD/m3)
        'X_c4': 0.0, # valerate and butyrate degraders (kg COD/m3)
        'X_pro': 0.0, # propionate degraders (kg COD/m3)
        'X_ac': 0.0, # acetate degraders (kg COD/m3)
        'X_h2': 0.0, # hydrogen degraders (kg COD/m3)
        'X_I': 17.0106, # particulate inerts (kg COD/m3)
        'S_cat': 0.0, # cations (base) (kmole/m3)
        'S_an': 0.0052101, # anions (acid) (kmole/m3)
        # 'Q': 178.4674, # Flow rate (m3/d)
        }
    
    matlab_postAD_adm = {
        'S_su': 0.012394,
        'S_aa': 0.0055432,
        'S_fa': 0.10741,
        'S_va': 0.012333,
        'S_bu': 0.014003,
        'S_pro': 0.017584,
        'S_ac': 0.089315,
        'S_h2': 2.5055e-07,
        'S_ch4': 0.05549,
        'S_IC': 0.095149*12,
        'S_IN': 1.3226,
        'S_I': 0.13087,
        'X_c': 0.10792,
        'X_ch': 0.020517,
        'X_pr': 0.08422,
        'X_li': 0.043629,
        'X_su': 0.31222,
        'X_aa': 0.93167,
        'X_fa': 0.33839,
        'X_c4': 0.33577,
        'X_pro': 0.10112,
        'X_ac': 0.67724,
        'X_h2': 0.28484,
        'X_I': 17.2162,
        'S_cat': 0., #-4.0789e-34,
        'S_an': 0.0052101
        }
    
    matlab_postAD_asm = {
        'S_I': 130.867, # soluble inert organic matter, mg COD/l
        'S_S': 258.5789, # readily biodegradable substrate, mg COD/l
        'X_I': 17216.2434, # particulate inert organic matter, mg COD/l
        'X_S': 2611.4843, # slowly biodegradable substrate, mg COD/l
        'X_BH': 0.0, # active heterotrophic biomass, mg COD/l
        'X_BA': 0.0, # active autotrophic biomass, mg COD/l
        'X_P': 626.0652, # particulate products arising from biomass decay, mg COD/l
        'S_O': 0.0, # dissolved O2, mg -COD/l
        'S_NO': 0.0, # nitrate and nitrite nitrogen, mg N/L
        'S_NH': 1442.7882, # ammonium, mg N/L
        'S_ND': 0.54323, # soluble biodegradable organic nitrogen
        'X_ND': 100.8668, # particulate biodegradable organic nitrogen, mg N/l
        'S_ALK': 97.8459*12, # alkalinity, assumed to be HCO3-, 97.8459, mol HCO3/m3 -> g C/m3
        'S_N2': 0.0, # dissolved O2
        # 'Q': 178.4674, # Flow rate, m3/d
        }


    adm1init = load_data(ospath.join(data_path, 'adm1init.csv'), index_col=0).to_dict('index')
    asm1_default_parameters = dict(
        mu_H = 4.0,
        K_S = 10.0,
        K_OH = 0.2,
        K_NO = 0.5,
        b_H = 0.3,
        mu_A = 0.5,
        K_NH = 1.0,
        K_OA = 0.4,
        b_A = 0.05,
        eta_g = 0.8,
        k_a = 0.05,
        k_h = 3.0,
        K_X = 0.1,
        eta_h = 0.8,
        Y_H = 0.67,
        Y_A = 0.24,
        f_P = 0.08,
        i_XB = 0.08,
        i_XP = 0.06,
        fr_SS_COD = 0.75
        )
    
    T = 273.15 + 35
    cmps_asm1 = pc.create_asm1_cmps()
    asm1 = pc.ASM1(components=cmps_asm1, **asm1_default_parameters)
    preAD_asm = WasteStream('preAD_asm', T=T)
    preAD_asm.set_flow_by_concentration(
        flow_tot=178.4674, 
        concentrations=dict(
            S_I = 28.0665,
            S_S = 48.9526,
            X_I = 10361.7101,
            X_S = 20375.0176,
            X_BH = 10210.0698,
            X_BA = 553.2808,
            X_P = 3204.6601,
            S_O = 0.25225,
            S_NO = 1.6871,
            S_NH = 28.9098,
            S_ND = 4.6834,
            X_ND = 906.0933,
            S_ALK = 7.1549*12
            ),
        units=('m3/d', 'mg/L')
        )
    thermo_asm1 = qs.get_thermo()
    cmps_adm1 = pc.create_adm1_cmps()
    adm1 = pc.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N # slight difference
    cmps_adm1.refresh_constants()
    thermo_adm1 = qs.get_thermo()
    
    J1 = su.ASMtoADM('J1', upstream=preAD_asm, downstream='preAD_adm', 
                     thermo=thermo_adm1, isdynamic=True, adm1_model=adm1,#)
                     T=T, pH=7.2631)
    AD1 = su.AnaerobicCSTR('AD1', ins=J1-0, outs=('biogas', 'postAD_adm'), 
                           isdynamic=True, V_liq=3400, V_gas=300, T=T,
                           model=adm1,)
    AD1.set_init_conc(**adm1init['AD1'])
    # Switch back to ASM1 components
    J2 = su.ADMtoASM('J2', upstream=AD1-1, downstream='postAD_asm', 
                     thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    J2.bio_to_xs = 0.79
    qs.set_thermo(thermo_asm1)
    
    sys = qs.System(path=(J1, AD1, J2))
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 200), method='BDF')
    fs = sys.flowsheet.stream
    
    for ws in sys.streams:
        ws.state[ws.state < 2.2e-16] = 0
        
    ac(cmps_adm1.kwarray(matlab_preAD_adm)[:-1]*1e3, fs.preAD_adm.state[:-2], rtol=1e-4)
    ac(cmps_adm1.kwarray(matlab_postAD_adm)[:-1]*1e3, fs.postAD_adm.state[:-2], rtol=1e-2)
    ac(cmps_asm1.kwarray(matlab_postAD_asm)[:-1], fs.postAD_asm.state[:-2], rtol=1e-3)
    
    h2 = cmps_adm1.S_h2
    ch4 = cmps_adm1.S_ch4
    co2 = cmps_adm1.S_IC
    assert np.isclose(AD1.state['S_h2_gas'] * h2.chem_MW / h2.i_mass, 1.1032e-5, rtol=1e-3)
    assert np.isclose(AD1.state['S_ch4_gas'] * ch4.chem_MW / ch4.i_mass, 1.6535, rtol=1e-2)
    assert np.isclose(AD1.state['S_IC_gas'], 0.01354, rtol=1e-2)
    assert np.isclose(AD1.outs[1].pH, 7.2631, rtol=1e-3)

    assert np.isclose(fs.biogas.imass['S_h2']*24 * h2.i_mass, 0.0035541, rtol=1e-2)
    assert np.isclose(fs.biogas.imass['S_ch4']*24 * ch4.i_mass, 1065.3523, rtol=1e-2)
    assert np.isclose(fs.biogas.imass['S_IC']*24 * co2.i_mass, 1535.4118, rtol=1e-2)
    

#%%

if __name__ == '__main__':
    test_adm1junctions()
