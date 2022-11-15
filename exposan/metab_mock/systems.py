# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
import qsdsan as qs
from qsdsan import sanunits as su, processes as pc, WasteStream, System
from qsdsan.utils import time_printer, ospath
from chemicals.elements import molecular_weight as get_mw
from exposan.metab_mock import DegassingMembrane as DM, METAB_AnCSTR as AB

folder = ospath.dirname(__file__)

__all__ = (
    'create_systems', 
    'default_inf_concs',
    'default_R1_init_conds',
    'default_R2_init_conds',
    'R1_ss_conds',
    'R2_ss_conds',
    'yields_bl', 'mus_bl', 'Ks_bl',
    'biomass_IDs',
    'vfa_IDs'
    )

#%% default values
Q = 5           # influent flowrate [m3/d]
T1 = 273.15+35  # temperature [K]
Vl1 = 5         # liquid volume [m^3]
Vg1 = 0.556     # headspace volume [m^3]
split_1 = 0.75  # split ratio to side-stream
tau_1 = 0.021   # degassing membrane retention time [d]

T2 = 273.15+25    
Vl2 = 75
Vg2 = 5
split_2 = 0.75
tau_2 = 0.021

fermenters = ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro')
methanogens = ('X_ac', 'X_h2')
biomass_IDs = (*fermenters, *methanogens)
vfa_IDs = ('S_va', 'S_bu', 'S_pro', 'S_ac')

C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})

default_inf_concs = {
    'S_su':3.0,
    'S_aa':0.6,
    'S_fa':0.4,
    'S_va':0.4,
    'S_bu':0.4,
    'S_pro':0.4,
    'S_ac':0.4,
    'S_h2':5e-9,
    'S_ch4':5e-6,
    'S_IC':0.04*C_mw,
    'S_IN':0.01*N_mw,
    'S_I':0.02,
    'X_c':0.1,
    'X_ch':0.3,
    'X_pr':0.5,
    'X_li':0.25,
    'X_aa':1e-3,
    'X_fa':1e-3,
    'X_c4':1e-3,
    'X_pro':1e-3, 
    'X_ac':1e-3, 
    'X_h2':1e-3, 
    'X_I':0.025, 
    'S_cat':0.04, 
    'S_an':0.02
    }

yields_bl = {
         'Y_su': 0.1,
         'Y_aa': 0.08,
         'Y_fa': 0.06,
         'Y_c4': 0.06,
         'Y_pro': 0.04,
         'Y_ac': 0.05,
         'Y_h2': 0.06
         }

mus_bl = np.array([5.0e-01, 1.0e+01, 1.0e+01, 1.0e+01, 3.0e+01, 5.0e+01, 6.0e+00,
                   2.0e+01, 2.0e+01, 1.3e+01, 8.0e+00, 3.5e+01, 2.0e-02, 2.0e-02,
                   2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02])

Ks_bl = np.array([5.0e-01, 3.0e-01, 4.0e-01, 2.0e-01, 
                  2.0e-01, 1.0e-01, 1.5e-01, 7.0e-06])

default_R1_init_conds = {
    'S_su': 0.0124*1e3,
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,
    'S_pro': 0.0176*1e3,
    'S_ac': 0.0893*1e3,
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ch': 0.0205*1e3,
    'X_pr': 0.0842*1e3,
    'X_li': 0.0436*1e3,
    'X_su': 1.87*1e3,
    'X_aa': 5.58*1e3,
    'X_fa': 2.03*1e3,
    'X_c4': 2.15*1e3,
    'X_pro': 1.00*1e3,
    }

default_R2_init_conds = {
    'S_su': 0.0124*1e3,
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,
    'S_pro': 0.0176*1e3,
    'S_ac': 0.0893*1e3,
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ac': 8.80*1e3,
    'X_h2': 3.70*1e3,
    }

R1_ss_conds = {
    'S_su': 0.0145871088552909*1e3,
    'S_aa': 0.00643308564144693*1e3,
    'S_fa': 0.634823005990967*1e3,
    'S_va': 0.624510322247682*1e3,
    'S_bu': 1.03793927591996*1e3,
    'S_pro': 1.24676871525373*1e3,
    'S_ac': 2.00250371674824*1e3,
    'S_h2': 0.00850364943532684*1e3,
    'S_ch4': 0.0000422133982597226*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.027310256066728*1e3,
    'X_c': 0.146203507058736*1e3,
    'X_ch': 0.0286018513139117*1e3,
    'X_pr': 0.0467836694957302*1e3,
    'X_li': 0.0247209587890493*1e3,
    'X_su': 4.69052782535406*1e3,
    'X_aa': 1.22829926704024*1e3,
    'X_fa': 0.0147446263753011*1e3,
    'X_c4': 0.0149933579422897*1e3,
    'X_pro': 0.0145343147735253*1e3,
    'X_ac': 0.00098041337766024*1e3,
    'X_h2': 0.00110808891184369*1e3,
    'X_I': 0.0396205121367899*1e3
    }

R2_ss_conds = {
    'S_su': 0.00106990968535691*1e3,
    'S_aa': 0.00125571416517827*1e3,
    'S_fa': 0.121097573221394*1e3,
    'S_va': 0.0132519103137696*1e3,
    'S_bu': 0.0172912281196732*1e3,
    'S_pro': 0.020032163173878*1e3,
    'S_ac': 0.00574002755366853*1e3,
    'S_h2': 3.76969944940856e-08*1e3,
    'S_ch4': 0.0499411746585487*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.105601391746794*1e3,
    'X_c': 0.0897520281015078*1e3,
    'X_ch': 0.00108163641708242*1e3,
    'X_pr': 0.00120204580901502*1e3,
    'X_li': 0.00150204523369107*1e3,
    'X_su': 0.195961987850137*1e3,
    'X_aa': 0.059723477130333*1e3,
    'X_fa': 0.0351858744892462*1e3,
    'X_c4': 0.0812315951844566*1e3,
    'X_pro': 0.0503466475437059*1e3,
    'X_ac': 1.1653549028287*1e3,
    'X_h2': 0.4352809013846*1e3,
    'X_I': 0.196117291164614*1e3
    }

# # R1_split = (0.9, 0.1), <0.2% diff if uses these
# R1_ss_conds = {
#     'S_su': 0.01431332653033124*1e3,
#     'S_aa': 0.006319229324446604*1e3,
#     'S_fa': 0.6348272949518201*1e3,
#     'S_va': 0.6245549045932813*1e3,
#     'S_bu': 1.038036275314748*1e3,
#     'S_pro': 1.2468433077261327*1e3,
#     'S_ac': 2.0025986879937308*1e3,
#     'S_h2': 0.011117814005879386*1e3,
#     'S_ch4': 5.952203539352817e-05*1e3,
#     'S_IC': 0.08330274953356481*C_mw*1e3,
#     'S_IN': 0.21612074810295764*N_mw*1e3,
#     'S_I': 0.02730942113623604*1e3,
#     'X_c': 0.14618842266689933*1e3,
#     'X_ch': 0.028601712932738325*1e3,
#     'X_pr': 0.04678353111459632*1e3,
#     'X_li': 0.024720751217894823*1e3,
#     'X_su': 4.6915202551212944*1e3,
#     'X_aa': 1.2274469593168187*1e3,
#     'X_fa': 0.014302724094919613*1e3,
#     'X_c4': 0.014436386133097296*1e3,
#     'X_pro': 0.014311685377540902*1e3,
#     'X_ac': 0.0009804172325111186*1e3,
#     'X_h2': 0.0011332727803773114*1e3,
#     'X_I': 0.039618842271433245*1e3
#     }

# R2_ss_conds = {
#     'S_su': 0.0008282629572271559*1e3,
#     'S_aa':0.001032694581229776*1e3,
#     'S_fa': 0.12736791568063396*1e3,
#     'S_va': 0.013392575245557583*1e3,
#     'S_bu': 0.01747268169433745*1e3,
#     'S_pro': 0.020406491318645616*1e3,
#     'S_ac': 0.00999472311030857*1e3,
#     'S_h2': 7.706131429350044e-08*1e3,
#     'S_ch4': 0.07044765251475965*1e3,
#     'S_IC': 0.49869691792599646*C_mw*1e3,
#     'S_IN': 0.20439321199252136*N_mw*1e3,
#     'S_I': 0.07350288752133147*1e3,
#     'X_c': 0.06137141129938335*1e3,
#     'X_ch': 0.0007989050950664583*1e3,
#     'X_pr': 0.0009192059152983118*1e3,
#     'X_li': 0.0010780568238771538*1e3,
#     'X_su': 0.19130327796202837*1e3,
#     'X_aa': 0.056026151032630665*1e3,
#     'X_fa': 0.03107211586931466*1e3,
#     'X_c4': 0.07959903792107152*1e3,
#     'X_pro': 0.04955864019101192*1e3,
#     'X_ac': 0.6361761447275743*1e3,
#     'X_h2':0.20772282044789275*1e3,
#     'X_I': 0.13200577400466326*1e3
#     }


# %%
# =============================================================================
# Preliminary analyses with mock METAB configuration
# =============================================================================

def create_systems(flowsheet_A=None, flowsheet_B=None, flowsheet_C=None, 
                   flowsheet_D=None, flowsheet_E=None,
                   inf_concs={}, R1_init_conds={}, R2_init_conds={}, which=None):
    which = which or ('A', 'B', 'C', 'D', 'E')
    if isinstance(which, str): which = (which,)
    ############# load components and set thermo #############
    pc.create_adm1_cmps()
    inf_concs = inf_concs or default_inf_concs
    brewery_ww = WasteStream('BreweryWW_A', T=T1)
    brewery_ww.set_flow_by_concentration(Q, concentrations=inf_concs, units=('m3/d', 'kg/m3'))
    ############# load process model ###########################
    adm1 = pc.ADM1()
    R1_init_conds = R1_init_conds or default_R1_init_conds
    R2_init_conds = R2_init_conds or default_R2_init_conds
    systems = []
    
    if 'A' in which:
        flowsheet_A = flowsheet_A or qs.Flowsheet('METAB_sysA')
        qs.main_flowsheet.set_flowsheet(flowsheet_A)
        ############# create WasteStream objects #################
        eff_A = WasteStream('Effluent_A', T=T2)
        bg1_A = WasteStream('biogas_1A', phase='g')
        bg2_A = WasteStream('biogas_2A', phase='g')
        flowsheet_A.stream.register(brewery_ww.ID, brewery_ww)
           
        ############# sysA unit operation ########################   
        H2E = AB('H2E', ins=[brewery_ww, 'return_1'], outs=('sidestream_1', ''), 
                split=(split_1, 1-split_1), V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                retain_cmps=fermenters)
        DM1 = DM('DM1', ins=H2E-0, outs=(bg1_A, 1-H2E), tau=tau_1)
        CH4E = AB('CH4E', ins=[H2E-1, 'return_2'], outs=('sidestream_2', eff_A), 
                split=(split_2, 1-split_2), V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                retain_cmps=methanogens)
        DM2 = DM('DM2', ins=CH4E-0, outs=(bg2_A, 1-CH4E), tau=tau_2)
        H2E.set_init_conc(**R1_ss_conds)
        CH4E.set_init_conc(**R2_ss_conds)
        # H2E.set_init_conc(**R1_init_conds)
        # CH4E.set_init_conc(**R2_init_conds)
        sysA = System('mock_METAB', 
                      path=(H2E, DM1, CH4E, DM2),
                      recycle=(DM1-1, DM2-1))
        sysA.set_dynamic_tracker(H2E, CH4E, bg1_A, bg2_A)
        systems.append(sysA)
    
    if 'B' in which:
        flowsheet_B = flowsheet_B or qs.Flowsheet('METAB_sysB')
        qs.main_flowsheet.set_flowsheet(flowsheet_B)
    
        ############# sysB streams ########################
        inf_b = brewery_ww.copy('BreweryWW_B')
        eff_B = WasteStream('Effluent_B', T=T2)
        bg1_B = WasteStream('biogas_1B', phase='g')
        bg2_B = WasteStream('biogas_2B', phase='g')
        
        ############# sysB unit operation #################
        AnR1 = su.AnaerobicCSTR('AnR1', ins=inf_b, outs=(bg1_B, ''), 
                                V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                                retain_cmps=fermenters)
        AnR2 = su.AnaerobicCSTR('AnR2', ins=AnR1-1, outs=(bg2_B, eff_B), 
                                V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                                retain_cmps=methanogens)
        # AnR1.set_init_conc(**R1_init_conds)
        # AnR2.set_init_conc(**R2_init_conds)
        AnR1.set_init_conc(**R1_ss_conds)
        AnR2.set_init_conc(**R2_ss_conds)
        sysB = System('baseline', path=(AnR1, AnR2))
        sysB.set_dynamic_tracker(AnR1, AnR2, bg1_B, bg2_B)
        systems.append(sysB)
        
    if 'C' in which:
        flowsheet_C = flowsheet_C or qs.Flowsheet('METAB_sysC')
        qs.main_flowsheet.set_flowsheet(flowsheet_C)
        
        ############# sysC streams ########################
        inf_c = brewery_ww.copy('BreweryWW_C')
        eff_c = WasteStream('Effluent_C', T=T2)
        bgm1 = WasteStream('biogas_mem_1', phase='g')
        bgm2 = WasteStream('biogas_mem_2', phase='g')
        bgh1 = WasteStream('biogas_hsp_1', phase='g')
        bgh2 = WasteStream('biogas_hsp_2', phase='g')
        
        ############# sysC unit operation #################
        R1 = su.AnaerobicCSTR('R1', ins=[inf_c, 'return_1'], 
                              outs=(bgh1, 'sidestream_1', ''), 
                              split=(split_1, 1-split_1),
                              V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                              retain_cmps=fermenters)
        DM1c = DM('DM1_c', ins=R1-1, outs=(bgm1, 1-R1), tau=tau_1)
        # DM1c = DM('DM1_c', ins=R1-1, outs=(bgm1, 1-R1), tau=0.1)    
    
        R2 = su.AnaerobicCSTR('R2', ins=[R1-2, 'return_2'], 
                              outs=(bgh2, 'sidestream_2', eff_c), 
                              split=(split_2, 1-split_2),
                              V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                              retain_cmps=methanogens)
        DM2c = DM('DM2_c', ins=R2-1, outs=(bgm2, 1-R2), tau=tau_2)
        # DM2c = DM('DM2_c', ins=R2-1, outs=(bgm2, 1-R2), tau=0.1)
        R1.set_init_conc(**R1_ss_conds)
        R2.set_init_conc(**R2_ss_conds)
        sysC = System('combined_METAB', path=(R1, DM1c, R2, DM2c),
                      recycle=(DM1c-1, DM2c-1))
        sysC.set_dynamic_tracker(R1, R2, bgm1, bgm2, bgh1, bgh2)
        systems.append(sysC)

    if 'D' in which:
        flowsheet_D = flowsheet_D or qs.Flowsheet('METAB_sysD')
        qs.main_flowsheet.set_flowsheet(flowsheet_D)
        
        ############# sysC streams ########################
        inf_d = brewery_ww.copy('BreweryWW_D')
        eff_d = WasteStream('Effluent_D', T=T2)
        bgm1_d = WasteStream('biogas_mem_1d', phase='g')
        bgm2_d = WasteStream('biogas_mem_2d', phase='g')
        bgh1_d = WasteStream('biogas_hsp_1d', phase='g')
        bgh2_d = WasteStream('biogas_hsp_2d', phase='g')
        
        ############# sysC unit operation #################
        R1d = su.AnaerobicCSTR('R1d', ins=[inf_d, 'return_1d'], 
                              outs=(bgh1_d, 'sidestream_1d', ''), 
                              split=(split_1, 1-split_1),
                              V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                              retain_cmps=fermenters)
        DM1d = DM('DM1d', ins=R1d-1, outs=(bgm1_d, 1-R1d), tau=tau_1)
    
        R2d = su.AnaerobicCSTR('R2d', ins=R1d-2, 
                              outs=(bgh2_d, ''), 
                              V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                              retain_cmps=methanogens)
        DM2d = DM('DM2d', ins=R2d-1, outs=(bgm2_d, eff_d), tau=tau_2)
        R1d.set_init_conc(**R1_ss_conds)
        R2d.set_init_conc(**R2_ss_conds)
        sysD = System('hybrid_METAB', path=(R1d, DM1d, R2d, DM2d),
                      recycle=(DM1d-1, ))
        sysD.set_dynamic_tracker(R1d, R2d, bgm1_d, bgm2_d, bgh1_d, bgh2_d)
        systems.append(sysD)
    
    if 'E' in which:
        flowsheet_E = flowsheet_E or qs.Flowsheet('METAB_sysE')
        qs.main_flowsheet.set_flowsheet(flowsheet_E)
        
        ############# sysC streams ########################
        inf_e = brewery_ww.copy('BreweryWW_E')
        eff_e = WasteStream('Effluent_E', T=T2)
        # bgm1_d = WasteStream('biogas_mem_1d', phase='g')
        bgm2_e = WasteStream('biogas_mem_2e', phase='g')
        bgh1_e = WasteStream('biogas_hsp_1e', phase='g')
        bgh2_e = WasteStream('biogas_hsp_2e', phase='g')
        
        ############# sysC unit operation #################
        R1e = su.AnaerobicCSTR('R1e', ins=inf_e, outs=(bgh1_e, ''),
                               fixed_headspace_P=True,
                               headspace_P=0.1013,
                               V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                               retain_cmps=fermenters)    
        R2e = su.AnaerobicCSTR('R2e', ins=R1e-1, outs=(bgh2_e, ''),                               
                              V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                              retain_cmps=methanogens)
        DM2e = DM('DM2e', ins=R2e-1, outs=(bgm2_e, eff_e), tau=tau_2)
        R1e.set_init_conc(**R1_ss_conds)
        R2e.set_init_conc(**R2_ss_conds)
        sysE = System('Best_METAB', path=(R1e, R2e, DM2e))
        sysE.set_dynamic_tracker(R1e, R2e, bgm2_e, bgh1_e, bgh2_e)
        systems.append(sysE)
    
    return systems

#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    global sysA, sysB, sysC, sysD, sysE
    sysA, sysB, sysC, sysD, sysE = create_systems()
    for sys in (sysA, sysB, sysC, sysD):
        print(f'Simulating {sys.ID}...')
        sys.simulate(state_reset_hook='reset_cache',
                    t_span=(0,t),
                    t_eval=np.arange(0, t+t_step, t_step),
                    method=method,
                    # export_state_to=ospath.join(folder, f'results/{method}_{t}d_{sys.ID[-4:]}.xlsx'),
                    **kwargs)


if __name__ == '__main__':
    t = 200
    t_step = 5
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    run(t, t_step, method=method)