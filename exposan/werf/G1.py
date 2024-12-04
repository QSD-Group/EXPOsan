# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''

import qsdsan as qs
from qsdsan import (
    WasteStream,
    processes as pc,
    sanunits as su,
    )
from qsdsan.utils import time_printer, ospath, load_data, get_SRT
# from exposan.bsm2 import data_path

__all__ = ('create_g1_system',)

#%%
MGD2cmd = 3785.412

Q = 10 * MGD2cmd # influent flowrate [MGD], converted to m3/d
V_tot = 4.7 * MGD2cmd # PFR total volume
V_fractions = [0.014, 0.13, 0.148, 0.148, 0.28, 0.28]

Q_pe = 38210 # primary effluent
Q_ps = 0.074 * MGD2cmd # primary sludge flowrate [MGD]
Q_intr = 40 * MGD2cmd # activated sludge process internal recycle [MGD]
# Q_ras = (Q-Q_ps)*0.4 # recycle sludge flowrate
Q_ras = 15280
Q_was = 0.1 * MGD2cmd # sludge wastage flowrate
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

folder = ospath.dirname(__file__)

# Default initial conditions
dfs = load_data(ospath.join(folder, 'data/G1_init.xlsx'), sheet=None)
inf_concs = dfs['asm'].iloc[0].to_dict()
# c1init = dfs['asm'].iloc[1].to_dict()
asinit = dfs['asm'].iloc[1:]
# asinit = dfs['asm_ss']
adinit = dfs['adm'].iloc[0].to_dict()
c2init = dfs['settler'].to_dict('index')
c2init['s'] = {k:v for k,v in c2init['s'].items() if v>0}
c2init['x'] = {k:v for k,v in c2init['x'].items() if v>0}
c2init['tss'] = [v for k,v in c2init['tss'].items() if v>0]
default_inf_kwargs = {
    'flow_tot': Q,
    'concentrations': inf_concs,
    'units': ('m3/d', 'mg/L'),
    }

# based on GPS-X model
asm_params = dict(
    Y_H=0.666, Y_A=0.18, Y_PAO=0.639, Y_PO4=0.4, Y_PHA=0.2,
    f_XI_H=0.08, f_XI_PAO=0.08, f_XI_AUT=0.08,
    mu_H=3.2, K_F=5.0, K_A_H=1.0, K_O2_H=0.2, K_NH4_H=0.05, K_P_H=0.01, K_NO3_H=0.2, 
    b_H=0.62, eta_NO3_H=0.64, eta_NO3_Hl=0.9,
    mu_AUT=0.9, K_NH4_AUT=0.7, K_O2_AUT=0.25, b_AUT=0.17, K_NO3_AUT=20, 
    eta_NO3_AUTl=0.5,
    mu_PAO=1.0, q_PHA=6.0, K_A_PAO=4.0, K_PP=0.01, K_PHA=0.01, K_O2_PAO=0.2,
    q_PP=1.5, K_MAX=0.34, K_IPP=0.02, K_P_S=0.2, b_PAO=0.2, b_PP=0.05, b_PHA=0.2,
    eta_NO3_PAOl=0.9, eta_NO3_PAO=0.24, 
    q_fe=3.0, eta_fe=0.3, K_fe=4.0, 
    k_h=3.0, K_X=0.1, eta_NO3=0.28, K_NO3=0.5, 
    # K_O2=0.2,
    # K_NO3_PAO=0.5, 
    # K_NH4_PAO=0.05, , 
    # K_P_PAO=0.01, K_P_AUT=0.01, 
    )

# based on GPS-X model
adm_params = dict(
    Y_su=0.18, Y_aa=0.18, Y_fa=0.18,
    f_ac_pro=0.04, Y_ac=0.05, Y_h2=0.06, Y_PO4=0.4, 
    q_PHA=6.0, K_A=4e-3, KS_IP=1e-5,
    # q_dis=0.5, q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
    # k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
    # K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6, 
    # K_PP=32e-5,
    # b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
    # b_PAO=0.2, b_PP=0.2, b_PHA=0.2, 
    # KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4, 
    )

ad_ss = {
    'S_su': 0.010064484725155741,
    'S_aa': 0.004342480337182479,
    'S_fa': 0.07403844415293871,
    'S_va': 0.010317043574002958,
    'S_bu': 0.013479943367162763,
    'S_pro': 0.01581654428555062,
    'S_ac': 0.042411409313503845,
    'S_h2': 1.9878056665990386e-07,
    'S_ch4': 0.17219325832127705,
    'S_IC': 0.9051039756054847,
    'S_IN': 0.8190596762571326,
    'S_IP': 0.7124668885680994,
    'S_I': 0.017788198946147093,
    'X_ch': 0.049922353655263416,
    'X_pr': 0.049922353655263416,
    'X_li': 0.051513392172419854,
    'X_su': 0.8893013022239371,
    'X_aa': 0.7020667535997138,
    'X_fa': 0.5411512808233088,
    'X_c4': 0.24495944677748568,
    'X_pro': 0.1316528841275703,
    'X_ac': 0.751859578120529,
    'X_h2': 0.40685460041930865,
    'X_I': 19.411526537242803,
    'X_PHA': 0.5360191866219218,
    'X_PP': 1.785325909966272e-05,
    'X_PAO': 0.6988820449303602,
    'S_K': 0.492020183107498,
    'S_Mg': 0.0002487848213584736,
    'S_Ca': 0.0009250351303063503,
    'X_CaCO3': 4.669524889784688e-05,
    'X_struv': 3.5717057089679782,
    'X_newb': 0.0,
    'X_ACP': 0.3707585018667417,
    'X_MgCO3': 8.863578420427e-05,
    'X_AlOH': 5.5382484477234465e-06,
    'X_AlPO4': 4.417924989999482e-06,
    'X_FeOH': 3.6730147014914605e-06,
    'X_FePO4': 4.423866992208256e-06,
    'S_Na': 0.08652984225675509,
    'S_Cl': 0.4210353282843405,
    }

#%%
def create_g1_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('G1')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm, 
                    electron_acceptor_dependent_decay=True,
                    # **asm_params
                    )
    thermo_asm = qs.get_thermo()
    
    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(**default_inf_kwargs)
    carb = WasteStream('carbon', T=Temp, units='kg/hr',
                        S_A=60, 
                       # S_A=80,
                       )
    
    PC = su.PrimaryClarifier(
        'PC', ins=(inf, 'reject'),
        outs=('PE', 'PS'),
        isdynamic=True, 
        sludge_flow_rate=Q_ps,
        solids_removal_efficiency=0.6,
        depth_clarifier=3.048,
        )
    
    GT = su.IdealClarifier('GT', PC-1, outs=['', 'thickened_PS'],
                             sludge_flow_rate=0.026*MGD2cmd,
                             solids_removal_efficiency=0.9,)
    
    gstrip = True
    anae_kwargs = dict(aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
    anox_kwargs = dict(aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
    ae_kwargs = dict(aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=gstrip)
    fc_kwargs = dict(
        underflow=Q_ras, wastage=Q_was,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01, 
        maximum_nonsettleable_solids=8.0
        )
    
    S1 = su.Splitter('S1', PC-0, split=0.8)
    
    A1 = su.CSTR('A1', ins=[carb, 'RAS'], V_max=V_tot*V_fractions[0], **anae_kwargs)
    A2 = su.CSTR('A2', [A1-0, S1-0], V_max=V_tot*V_fractions[1], **anae_kwargs)
    A3 = su.CSTR('A3', [A2-0, 'RWW', S1-1], V_max=V_tot*V_fractions[2], **anox_kwargs)
    A4 = su.CSTR('A4', A3-0, V_max=V_tot*V_fractions[3], **anox_kwargs)
    O1 = su.CSTR('O1', A4-0, V_max=V_tot*V_fractions[4], **ae_kwargs)
    O2 = su.CSTR('O2', O1-0, [1-A3, 'treated'], split=[Q_intr, Q+Q_ras],
                  V_max=V_tot*V_fractions[5], **ae_kwargs)

    FC = su.FlatBottomCircularClarifier(
        'FC', O2-1, ['effluent', 1-A1, 'WAS'], 
        **fc_kwargs
        )
    # AS = su.PFR('AS', ins=[PC-0, 'RAS', carb], outs='treated', 
    #             N_tanks_in_series=6,
    #             V_tanks=[f*V_tot for f in V_fractions],
    #             influent_fractions=[
    #                 [0, 0.8, 0.2, 0,0,0], # PC-0
    #                 [1,0,0,0,0,0], # RAS
    #                 [1,0,0,0,0,0], # carb
    #                 ], 
    #             internal_recycles=[(5,2,Q_intr)], 
    #             kLa=[0,0,0,0,180,70],
    #             DO_setpoints=[0,0,0,0,2.0,2.0], DO_ID='S_O2',
    #             suspended_growth_model=asm,
    #             gas_stripping=gstrip)
    
    # FC = su.FlatBottomCircularClarifier(
    #     'FC', AS-0, ['effluent', 1-AS, 'WAS'],
    #     **fc_kwargs        
    #     )
    
    MT = su.IdealClarifier('MT', FC-2, outs=['', 'thickened_WAS'],
                           sludge_flow_rate=0.019*MGD2cmd,
                           solids_removal_efficiency=0.95,)
    M1 = su.Mixer('M1', ins=(GT-1, MT-1))
        
    # Switch to ADM1 components for the anaerobic digester
    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0, 
                   # **adm_params
                   )
    
    # breakpoint()
    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR('AD', ins=J1.outs[0], outs=('biogas', 'digestate'), 
                          # isdynamic=True, V_liq=416.395*9, V_gas=416.395, 
                          isdynamic=True, V_liq=3596.14, V_gas=416.395, 
                          fixed_headspace_P=False, fraction_retain=0,
                          T=T_ad, model=adm,
                          pH_ctrl=7.0,
                          )
    AD.algebraic_h2 = False
    # Switch back to ASM1 components
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)
    
    # Dewatering
    DW = su.PrimaryClarifier('DW', J2-0, outs=['', 'cake'],
                             sludge_flow_rate=5.93e-3*MGD2cmd,
                             solids_removal_efficiency=0.9,)

    M2 = su.Mixer('M2', ins=(GT-0, MT-0, DW-0))
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-PC)

    
    if default_init_conds:
        asdct = asinit.to_dict('index')
        for i in (A1, A2, A3, A4, O1, O2):
            i.set_init_conc(**asdct[i.ID])
        # AS.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**c2init['s'])
        FC.set_init_sludge_solids(**c2init['x'])
        FC.set_init_TSS(c2init['tss'])
        # AD.set_init_conc(**adinit)
        AD.set_init_conc(**default_ad_init)
        # AD.set_init_conc(**ad_ss)
    
    sys = qs.System('G1', 
                    path=(PC, GT, S1, A1, A2, A3, A4, O1, O2, FC, 
                          MT, M1, J1, AD, J2, DW, M2, HD),
                    recycle=(O2-0, FC-1, HD-0)
                    # path=(PC, GT, AS, FC, MT, M1, J1, AD, J2, DW, M2, HD),
                    # recycle=(FC-1, HD-0)
                    )
    # sys.set_tolerance(mol=1e-4, rmol=1e-3)
    # sys.maxiter = 5000
    # sys.set_dynamic_tracker(A1, O2, AD, FC-0)
    # sys.set_dynamic_tracker(AS, AD, FC-0)
    sys.set_dynamic_tracker(FC-0, AD-1)
    
    return sys

#%%
default_asm2d_kwargs = dict(iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
            iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
            iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
            f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
            Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
            Y_A=0.24, f_XI_AUT=0.1,
            K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
            mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
            K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
            q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
            b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
            K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
            K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
            mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
            k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, 
            # path=os.path.join(data_path, '_asm2d.tsv'),
            )
# Saumitra's
# domestic_ww = {
#    'S_I': 20,
#    'X_I': 40,
#    'S_F': 45,
#    'S_A': 63,
#    'X_S': 160,
#    'S_NH4': 25,
#    'S_PO4': 4.5,
#    'X_PP': 0,
#    'X_PHA': 10,
#    'X_H': 10,
#    'X_AUT': 10, 
#    'X_PAO': 5, 
#    'X_MeOH': 32, 
#    'S_ALK':7*12,
#     }

domestic_ww = {
   'S_I': 18,
   'X_I': 71.6,
   'S_F': 57.3,
   'S_A': 8.6,
   'X_S': 202.6,
   'S_NH4': 25,
   'S_PO4': 5,
   'S_ALK':75.6,
    }

steady_state_ad_init_conds = {
      'S_su': 1.094e+01,
       # 'S_aa': 4.713e+00, # original ss
      # 'S_aa': 4.713e+01,
      'S_aa': 4.713e+03,
      'S_fa': 9.150e+01,
      'S_va': 9.512e+00,
      'S_bu': 1.255e+01,
      'S_pro': 1.486e+01,
      'S_ac': 3.749e+01,
       'S_h2': 2.193e-04,
      'S_ch4': 7.692e+01,
      'S_IC': 6.708e+02,
      'S_IN': 1.035e+03,
      'S_IP': 1.035e+00, #random
      'S_I': 2.325e+02,
      # 'X_ch': 6.432e+01, # original ss
      # 'X_ch': 6.432e+00, 
      'X_ch': 6.432e-01,
      'X_pr': 6.719e+01,
      # 'X_li': 1.436e+02, # original ss
      # 'X_li': 1.436e+01,
      'X_li': 1.436e+00,
      'X_su': 1.116e+03,
      'X_aa': 8.755e+02,
      'X_fa': 1.274e+03,
      'X_c4': 3.730e+02,
      'X_pro': 1.751e+02,
      'X_ac': 1.405e+03,
      'X_h2': 6.820e+02,
      'X_I': 2.413e+04,
      'X_PHA': 10.000e-01, 
      'X_PP':  10.000e-01,  
      'X_PAO': 10.000e+00, 
      'S_K': 10.000e-001, 
      'S_Mg': 10.000e-001, 
      'X_MeOH': 10.000e-001, 
      'X_MeP': 10.000e-001, 
      'S_cat': 0.000e+00, 
      'S_an': 4.772e+00
      }

def batch_init(sys, path, sheet):
    # isa = isinstance
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    u.DG.set_init_conc(**steady_state_ad_init_conds)
    # for k in sys.units:
    #     if k.ID.startswith('O', 'A'): 
    #         k.set_init_conc(**dct[k.ID])
    u.AS.set_init_conc(concentrations=df.iloc[:6])
    c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    tss = [v for v in dct['C1_tss'].values() if v>0]
    u.FC.set_init_solubles(**c1s)
    u.FC.set_init_sludge_solids(**c1x)
    u.FC.set_init_TSS(tss)
    
def create_oldg1(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Components and stream
    cmps = pc.create_asm2d_cmps()
    qs.set_thermo(cmps)
    thermo_asm2d = qs.get_thermo()
    
    # dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    # ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    # # acetic_acid = qs.WasteStream('organic_carbon', S_A= 59.55, units='kg/hr', T=Temp)
    
    # # inf_kwargs = inf_kwargs or default_inf_kwargs
    # dom_ww.set_flow_by_concentration(Q_domestic, 
    #                                  concentrations=domestic_ww, 
    #                                  units=('m3/d', 'mg/L'))
    # inf = WasteStream('inf', T=Temp)
    # inf.set_flow_by_concentration(**default_inf_kwargs)
    inf = WasteStream('dom_ww', T=Temp)
    inf.set_flow_by_concentration(Q, domestic_ww, units=('m3/d', 'mg/L'))
    carb = WasteStream('carbon', T=Temp, 
                       S_A=107.6, units='kg/hr')
           
    asm2d = pc.ASM2d()

    PC = su.PrimaryClarifier(
        'PC', ins=(inf, 'reject'),
        outs=('PE', 'PS'),
        isdynamic=True, 
        sludge_flow_rate=Q_ps,
        solids_removal_efficiency=0.6,
        depth_clarifier=3.048,
        )
    
    # GT = su.PrimaryClarifier('GT', PC-1, outs=['', 'thickened_PS'],
    #                          sludge_flow_rate=0.026*MGD2cmd,
    #                          solids_removal_efficiency=0.9,)
    
    AS = su.PFR('AS', ins=[PC-0, 'RAS', carb], outs='treated', 
                N_tanks_in_series=6,
                V_tanks=[f*V_tot for f in V_fractions],
                influent_fractions=[
                    [0, 0.8, 0.2, 0,0,0], # PC-0
                    [1,0,0,0,0,0], # RAS
                    [1,0,0,0,0,0], # carb
                    ], 
                internal_recycles=[(5,2,Q_intr)], kLa=[],
                DO_setpoints=[0,0,0,0,2.0,2.0], DO_ID='S_O2',
                suspended_growth_model=asm2d,
                gas_stripping=False)

    FC = su.FlatBottomCircularClarifier('FC', AS-0, ['effluent', 1-AS, 'WAS'],
                                        underflow=Q_ras, wastage=Q_was, 
                                        surface_area=1580, height=3.66, N_layer=10, 
                                        feed_layer=5, thermo = thermo_asm2d)
    
    # thickener_perc and TSS_removal_perc based on WERF report
    GT = su.Thickener('GT', PC-1, ['sludge_GT', 'eff_GT'], 
                       thickener_perc=7, TSS_removal_perc=92)
    
    # thickener_perc and TSS_removal_perc based on WERF report
    MT = su.Thickener('MT', FC-2, ['sludge_MT', 'eff_MT'],
                       thickener_perc=6, TSS_removal_perc=98)
    
    M1 = su.Mixer('M1', [MT-0, GT-0])
    
    cmps_adm1 = qs.processes.create_adm1_p_extension_cmps()
    cmps_adm1.X_PAO.i_N = 0.07
    cmps_adm1.X_PAO.i_P = 0.02
    cmps_adm1.refresh_constants()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1_p_extension()
    
    J1 = su.ASM2dtomADM1('J1', upstream= [M1-0], thermo=thermo_adm1, isdynamic=True, 
                         adm1_model=adm1, asm2d_model=asm2d)
      
    DG = su.AnaerobicCSTR(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
                          model=adm1, thermo = thermo_adm1, V_liq=3596, V_gas=416)
    DG.algebraic_h2 = True
    
    J2 = su.mADM1toASM2d('J2', upstream = DG-1, thermo=thermo_asm2d, isdynamic=True, 
                         adm1_model=adm1, asm2d_model=asm2d)
    qs.set_thermo(thermo_asm2d)

    # thickener_perc and TSS_removal_perc based on WERF report
    DU = su.Centrifuge('DU', ins = J2.outs[0], outs = ['sludge_DU', 'eff_DU'], 
                       thickener_perc=23, TSS_removal_perc=95)
    
    M2 = su.Mixer('M2', ins=[GT-1, MT-1, DU-1], outs=1-PC)

    sys = qs.System('G1_WERF', path=(PC, GT, AS, FC, MT, M1, J1, DG, J2, DU, M2), 
                    recycle = [FC-1, M2-0,])

    sys.set_tolerance(rmol=1e-6)
    sys.maxiter = 500
    sys.set_dynamic_tracker(FC-0)
    
    path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
    batch_init(sys, path, 
               sheet='G1')
    
    return sys

#%%
default_PE_concs = dict(
    S_N2=18, S_NH4=25, S_F=87.0, S_I=21.8, X_S=112.7+39.6, X_I=29.0, S_PO4=8.0,
    S_IC=75.6, S_Ca=140, S_Mg=50, S_K=28, S_Na=3.76*23, S_Cl=12*35.45
    )

# g1_PE_concs = dict(
#     S_O2	=1.64E-08, S_F=5.66E+01, S_A=8.531, S_I=17.9, S_NH4=29.58, S_N2=17.93, 
#     S_NO3=0.05754, S_PO4=8.191, S_ALK=77.23, X_I=32.3789, X_S=100.6, X_H=0.3969, 
#     X_PAO=0.3778, X_PP=0.08498, X_PHA=0.07941, X_AUT=0.08623, 
#     X_MeOH=1.24E-07*2, X_MeP=1.53E-07*2,
#     )

g1_PE_concs = dict(
    S_O2=1.64E-08, S_F=5.66E+01, S_A=8.531, S_I=17.9, S_NH4=29.58, S_N2=17.93, 
    S_NO3=0.05754, S_PO4=8.191, S_IC=77.23, X_I=32.3789, X_S=100.6, X_H=0.3969, 
    X_PAO=0.3778, X_PP=0.08498, X_PHA=0.07941, X_AUT=0.08623, 
    S_K=29.67, S_Mg=49.8, S_Na=87.561, S_Cl=425.4, S_Ca=139.6, 
    X_AlOH=1.24E-07, X_AlPO4=1.53E-07, X_FeOH=1.23E-07, X_FePO4=1.52E-07, 
    X_CaCO3=1.24E-05, X_ACP=5.21E-02, X_MgCO3=7.01E-06, X_newb=3.61E-05, X_struv=5.40E-01
    )

def create_subsys(modified=True):
    flowsheet = qs.Flowsheet('G1_2nd')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    if modified:
        cmps_asm = pc.create_masm2d_cmps()
        # cmps_asm.S_F.i_C = 0.32
        # cmps_asm.S_F.i_N = 0.03
        # cmps_asm.S_F.i_P = 0.01
        # cmps_asm.S_I.i_C = 0.32
        # cmps_asm.S_I.i_N = 0.01
        # cmps_asm.S_I.i_P = 0.0
        # cmps_asm.X_I.i_C = 0.366697674
        # cmps_asm.X_I.i_N = 0.02
        # cmps_asm.X_I.i_P = 0.01
        # cmps_asm.X_S.i_C = 0.32
        # cmps_asm.X_S.i_N = 0.04
        # cmps_asm.X_S.i_P = 0.01
        # cmps_asm.X_PHA.i_C = 1/3
        # for cmp in (cmps_asm.X_H, cmps_asm.X_PAO, cmps_asm.X_AUT):
        #     cmp.i_C = 0.366
        #     cmp.i_N = 0.07
        #     cmp.i_P = 0.02
        # cmps_asm.refresh_constants()
        asm = pc.mASM2d(components=cmps_asm, 
                        electron_acceptor_dependent_decay=False,
                        # **asm_params
                        )
    else:
        cmps_asm = pc.create_asm2d_cmps()
        asm = pc.ASM2d()
        
    # Influent
    PE = WasteStream('PE', T=Temp)
    PE.set_flow_by_concentration(Q, g1_PE_concs, units=('m3/d', 'mg/L'))
    # if modified: PE.set_flow_by_concentration(Q_pe, g1_PE_concs, units=('m3/d', 'mg/L'))
    # else: PE.set_flow_by_concentration(Q_pe, domestic_ww, units=('m3/d', 'mg/L'))
    carb = WasteStream('carbon', T=Temp, units='kg/hr',
                        S_A=60, 
                       # S_A=107.6, 
                       )
    gstrip = True
    fc_kwargs = dict(
        underflow=Q_ras, wastage=Q_was,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01,
        maximum_nonsettleable_solids=8.0
        )

    # anae_kwargs = dict(aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
    # anox_kwargs = dict(aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
    # ae_kwargs = dict(aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=gstrip)
    
    # S1 = su.Splitter('S1', PE, split=0.8)
    
    # A1 = su.CSTR('A1', ins=[carb, 'RAS'], V_max=V_tot*V_fractions[0], **anae_kwargs)
    # A2 = su.CSTR('A2', [A1-0, S1-0], V_max=V_tot*V_fractions[1], **anae_kwargs)
    # A3 = su.CSTR('A3', [A2-0, 'RWW', S1-1], V_max=V_tot*V_fractions[2], **anox_kwargs)
    # A4 = su.CSTR('A4', A3-0, V_max=V_tot*V_fractions[3], **anox_kwargs)
    # O1 = su.CSTR('O1', A4-0, V_max=V_tot*V_fractions[4], **ae_kwargs)
    # O2 = su.CSTR('O2', O1-0, [1-A3, 'treated'], split=[Q_intr, Q+Q_ras],
    #               V_max=V_tot*V_fractions[5], **ae_kwargs)

    # FC = su.FlatBottomCircularClarifier(
    #     'FC', O2-1, ['effluent', 1-A1, 'WAS'], 
    #     **fc_kwargs
    #     )

    AS = su.PFR('AS', ins=[PE, 'RAS', carb], outs='treated', 
                N_tanks_in_series=6,
                V_tanks=[f*V_tot for f in V_fractions],
                influent_fractions=[
                    [0, 0.8, 0.2, 0,0,0], # PE
                    [1,0,0,0,0,0], # RAS
                    [1,0,0,0,0,0], # carb
                    ], 
                internal_recycles=[(5,2,Q_intr)], 
                kLa=[0,0,0,0,180,70],
                DO_setpoints=[0,0,0,0,2.0,2.0], DO_ID='S_O2',
                suspended_growth_model=asm,
                # gas_stripping=modified,
                gas_stripping=gstrip,
                )
    
    FC = su.FlatBottomCircularClarifier(
        'FC', AS-0, ['effluent', 1-AS, 'WAS'],
        **fc_kwargs
        )

    if modified:
        # asdct = asinit.to_dict('index')
        # for i in (A1, A2, A3, A4, O1, O2):
        #     i.set_init_conc(**asdct[i.ID])
        AS.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**c2init['s'])
        FC.set_init_sludge_solids(**c2init['x'])
        FC.set_init_TSS(c2init['tss'])
    else:
        path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
        df = load_data(path, 'G1')
        dct = df.to_dict('index')
        # for i in (A1, A2, A3, A4, O1, O2):
        #     i.set_init_conc(**asdct[i.ID])
        AS.set_init_conc(concentrations=df.iloc[:6])
        c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
        c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
        tss = [v for v in dct['C1_tss'].values() if v>0]
        FC.set_init_solubles(**c1s)
        FC.set_init_sludge_solids(**c1x)
        FC.set_init_TSS(tss)
    
    sub = qs.System('G1_2nd', path=(AS, FC), recycle=(FC-1,))
    sub.set_dynamic_tracker(FC-0)
    # sub = qs.System('G1_2nd', path=(S1, A1, A2, A3, A4, O1, O2, FC), recycle=(O2-0, FC-1))
    # sub.set_dynamic_tracker(A1, O2, FC-0)
    
    return sub

#%%
default_ad_init = dict(
    S_su=0.01, S_aa=0.005, S_fa=0.1, S_va=0.01, S_bu=0.015, S_pro=0.015, S_ac=0.18, 
    S_ch4=0.05, S_IC=0.024, S_IN=3.2, S_IP=0.005, S_I=6.05, 
    X_ch=1.15, X_pr=1.15, X_li=1.73, X_su=0.85, X_aa=0.6, X_fa=0.7, X_c4=0.3, 
    X_pro=0.135, X_ac=0.9, X_h2=0.435, X_I=46.1, 
    S_Ca=0.01, S_Mg=0.05, S_K=0.02, 
    X_CaCO3=1e-8, X_struv=1e-8, X_newb=1e-8, X_ACP=1e-8, 
    X_MgCO3=1e-8, X_AlOH=1e-8, X_AlPO4=1e-8, X_FeOH=1e-8, 
    X_FePO4=1e-8, 
    )

was_mass = {
    'S_O2': 0.03154509999999999,
    'S_N2': 0.23018691858283896,
    'S_NH4': 0.008935722528940239,
    'S_NO3': 0.11657850984181048,
    'S_PO4': 0.0023555001464213928,
    'S_F': 0.006888621098589779,
    'S_A': 0.00010550461955551915,
    'S_I': 0.2823162503117132,
    'S_IC': 0.035558870264853236,
    'S_K': 0.41624574645374185,
    'S_Mg': 0.7532965348442451,
    'X_I': 74.56537256203178,
    'X_S': 1.4695474553573085,
    'X_H': 39.50405700267845,
    'X_PAO': 20.55009363705345,
    'X_PP': 11.58305521245629,
    'X_PHA': 0.04721213121620314,
    'X_AUT': 3.4599323959750103,
    'S_Ca': 2.201751315280157,
    'X_CaCO3': 1.805048459224936e-05,
    'X_struv': 0.7860568382689709,
    'X_newb': 5.254964663109509e-05,
    'X_ACP': 0.0758399285425468,
    'X_MgCO3': 1.0204473458132814e-05,
    'X_AlOH': 1.807923808593028e-07,
    'X_AlPO4': 2.2300674620789448e-07,
    'X_FeOH': 1.7933670882982369e-07,
    'X_FePO4': 2.2155108513492553e-07,
    'S_Na': 1.3809996197510446,
    'S_Cl': 6.709348205731937,
    'H2O': 15634.227546353619
    }

combined_PI_conc = dict(
    S_O2	=1.643E-08, S_F=56.62187102, S_A=8.531, S_I=17.9, S_NH4=29.58, S_N2=17.93,
    S_NO3=0.05756, S_PO4=8.189, S_IC=77.23, X_I=80.943, X_S=206.45, X_H=0.9921, 
    X_PAO=0.9439, X_PP=0.2123, X_PHA=0.1984, X_AUT=0.21554, S_K=29.67, S_Mg=49.8,
    S_Na=87.561, S_Cl=425.4, S_Ca=139.6, X_AlOH=0.000002911, X_AlPO4=0.000002913, 
    X_FeOH=0.000001932, X_FePO4=0.000002779, X_CaCO3=0.00003092, X_ACP=0.1301, 
    X_MgCO3=0.00001751, X_newb=0.00009022, X_struv=1.348,
    )

def create_solidsys():
    flowsheet = qs.Flowsheet('G1_solids')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm, 
                    electron_acceptor_dependent_decay=False,
                    # **asm_params
                    )
    thermo_asm = qs.get_thermo()
    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(**default_inf_kwargs)
    # infreject = WasteStream('PI', T=Temp)
    # infreject.set_flow_by_concentration(38486.31, combined_PI_conc, ('m3/d', 'mg/L')) # PI in GPS-X
    WAS = WasteStream('WAS', **was_mass)

    
    PC = su.PrimaryClarifier(
        'PC', 
        ins=(inf, 'reject'),
        # ins=infreject,
        outs=('PE', 'PS'),
        isdynamic=True, 
        sludge_flow_rate=Q_ps,
        solids_removal_efficiency=0.6,
        depth_clarifier=3.048,
        )
    
    GT = su.IdealClarifier('GT', PC-1, outs=['', 'thickened_PS'],
                           sludge_flow_rate=0.026*MGD2cmd,
                           solids_removal_efficiency=0.9,)
    
    MT = su.IdealClarifier('MT', WAS, outs=['', 'thickened_WAS'],
                           sludge_flow_rate=0.019*MGD2cmd,
                           solids_removal_efficiency=0.95,)
    M1 = su.Mixer('M1', ins=(GT-1, MT-1))
    
    # Switch to ADM1 components for the anaerobic digester
    cmps_adm = pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0, 
                   # **adm_params
                   )
    
    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR('AD', ins=J1.outs[0], outs=('biogas', 'digestate'), 
                          isdynamic=True, 
                           V_liq=3596.14, 
                          # V_liq=416.395*9, 
                          V_gas=416.395, 
                          fixed_headspace_P=False, fraction_retain=0,
                          T=T_ad, model=adm,
                          pH_ctrl=7.0,
                          )
    AD.algebraic_h2 = False
    # Switch back to ASM1 components
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)
    
    # Dewatering
    DW = su.PrimaryClarifier('DW', J2-0, outs=['', 'cake'],
                             sludge_flow_rate=5.93e-3*MGD2cmd,
                             solids_removal_efficiency=0.9,)

    M2 = su.Mixer('M2', ins=(GT-0, MT-0, DW-0), 
                   # outs=1-PC,
                   # outs=['reject']
                  )
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-PC)
    
    # AD.set_init_conc(**adinit)
    AD.set_init_conc(**default_ad_init)
    # AD.set_init_conc(**ad_ss)
    
    sys = qs.System('G1_solids', 
                    # path=(PC, GT, MT, M1, J1, AD, J2, DW, M2),
                    # recycle=(M2-0,)
                    path=(PC, GT, MT, M1, J1, AD, J2, DW, M2, HD),
                    recycle=(HD-0,)
                    )
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    # sys.maxiter = 5000
    # sys.set_dynamic_tracker(A1, O2, AD, FC-0)
    # sys.set_dynamic_tracker(AS, AD, FC-0)
    sys.set_dynamic_tracker(PC-0, AD-1)
    return sys

#%%

def create_ad():
    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm, 
                    # electron_acceptor_dependent_decay=False,
                    # **asm_params
                    )
    thermo_asm = qs.get_thermo()
    
    cmps_adm = pc.create_adm1p_cmps()
    # thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0, 
                   # **adm_params
                   )
    thickened_sludge = WasteStream('thickened_sludge')
    thickened_sludge.set_flow_by_concentration(
        flow_tot=168.5, 
        concentrations=dict(
            S_su=11.305, S_aa=11.305, S_fa=11.305, 
            S_va=4.26E-07, S_pro=0.001426, S_ac=4.984, 
            S_IC=53.82, S_IN=19.65, S_IP=5.269, S_I=17.9, 
            X_ch=9453.5, X_pr=9453.5, X_li=9453.5, 
            X_su=150.53, X_aa=150.53, X_fa=150.54, 
            X_c4=0.2075, X_pro=0.2075, X_ac=66.15, X_h2=79.53, X_I=19145, 
            X_PHA=32.76, X_PP=1124, X_PAO=3682, 
            S_K=28.5, S_Mg=49.1, S_Ca=139.5, 
            X_CaCO3=0.04714, X_struv=173.6, X_newb=0.4121, X_ACP=16.82, 
            X_MgCO3=0.08948, X_AlOH=0.005591, X_AlPO4=0.00446, X_FeOH=0.003708, 
            X_FePO4=0.004466, S_Na=87.354, S_Cl=425.0455,
            ),
        units=('m3/d', 'mg/L')
        )
    
    AD = su.AnaerobicCSTR('AD', ins=thickened_sludge, outs=('biogas', 'digestate'), 
                          isdynamic=True, 
                          V_liq=3596.14, 
                          # V_liq=416.395*9, 
                          V_gas=416.395, 
                          fixed_headspace_P=False, fraction_retain=0,
                          T=T_ad, model=adm,
                          pH_ctrl=7.0,
                          )
    AD.algebraic_h2 = False
    # Switch back to ASM1 components
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    # qs.set_thermo(thermo_asm)
    # AD.set_init_conc(**adinit)
    AD.set_init_conc(**default_ad_init)

    sys = qs.System('G1_ad', path=(AD, J2))
    sys.set_dynamic_tracker(AD, AD-0, AD-1, J2-0)
    
    return sys

# %%

@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        # t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    
    # biomass_IDs = ('X_BH', 'X_BA')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[sys.flowsheet.stream.digested_sludge],
    #               active_unit_IDs=('C3'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

#%%
if __name__ == '__main__':
    sys = create_g1_system()
    # sys = create_subsys(False)
    # sys = create_subsys()
    # sys = create_solidsys()
    # sys = create_oldg1()
    # sys = create_ad()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 300
    # t = 1
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    # sys._setup()
    # sys.converge()
    # sys.diagram()
    # sys.diagram(file=os.path.join(figures_path, 'bsm2_sys'), format='png')