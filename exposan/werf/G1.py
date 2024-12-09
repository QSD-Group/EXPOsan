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
from exposan.werf import data_path#, default_fctss_init

__all__ = ('create_g1_system',)

ID = 'G1'
#%%
dfs = load_data(
    ospath.join(data_path, 'initial_conditions.xlsx'), 
    sheet=None,
    )
asinit = dfs[ID]
fcinit = asinit.iloc[-1].to_dict()
# default_fctss_init = [18, 235, 2702, 7114, 7114, 9889,
#                       11316, 12362, 13409, 14968]
default_fctss_init = [10, 12, 12, 40, 500, 500,
                      500, 5000, 10000, 13000]
adinit = dfs['adm'].loc[ID].to_dict()

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
Q_was = 0.136 * MGD2cmd # sludge wastage flowrate
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

folder = ospath.dirname(__file__)

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
    # inf = WasteStream('inf', T=Temp)
    # inf.set_flow_by_concentration(**default_inf_kwargs)
    inf = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
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

    S1 = su.Splitter('S1', PC-0, split=0.8)
    
    A1 = su.CSTR('A1', ins=[carb, 'RAS'], V_max=V_tot*V_fractions[0], **anae_kwargs)
    A2 = su.CSTR('A2', [A1-0, S1-0], V_max=V_tot*V_fractions[1], **anae_kwargs)
    A3 = su.CSTR('A3', [A2-0, 'intr', S1-1], V_max=V_tot*V_fractions[2], **anox_kwargs)
    A4 = su.CSTR('A4', A3-0, V_max=V_tot*V_fractions[3], **anox_kwargs)
    O1 = su.CSTR('O1', A4-0, V_max=V_tot*V_fractions[4], **ae_kwargs)
    O2 = su.CSTR('O2', O1-0, [1-A3, 'treated'], split=[Q_intr, Q+Q_ras],
                  V_max=V_tot*V_fractions[5], **ae_kwargs)

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
    
    fc_kwargs = dict(
        underflow=Q_ras, wastage=Q_was,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01, 
        maximum_nonsettleable_solids=8.0
        )
    
    FC = su.FlatBottomCircularClarifier(
        # 'FC', AS-0, ['effluent', 1-AS, 'WAS'],
        'FC', O2-1, ['effluent', 1-A1, 'WAS'], 
        **fc_kwargs        
        )
    
    MT = su.IdealClarifier('MT', FC-2, outs=['', 'thickened_WAS'],
                           sludge_flow_rate=0.0335*MGD2cmd,
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
                          isdynamic=True, V_liq=0.95*MGD2cmd, V_gas=0.11*MGD2cmd, 
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
                             sludge_flow_rate=0.0095*MGD2cmd,
                             solids_removal_efficiency=0.9,)

    M2 = su.Mixer('M2', ins=(GT-0, MT-0, DW-0))
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-PC)

    
    if default_init_conds:
        asdct = asinit.to_dict('index')
        for i in (A1, A2, A3, A4, O1, O2):
            i.set_init_conc(**asdct[i.ID])
        # AS.set_init_conc(concentrations=asinit)
        # FC.set_init_solubles(**c2init['s'])
        # FC.set_init_sludge_solids(**c2init['x'])
        # FC.set_init_TSS(c2init['tss'])
        FC.set_init_solubles(**fcinit)
        FC.set_init_sludge_solids(**fcinit)
        FC.set_init_TSS(default_fctss_init)
        AD.set_init_conc(**adinit)
        # AD.set_init_conc(**default_ad_init)
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
    sys.set_dynamic_tracker(FC-0, AD)
    
    return sys

#%%
# default_ad_init = dict(
#     S_su=0.01, S_aa=0.005, S_fa=0.1, S_va=0.01, S_bu=0.015, S_pro=0.015, S_ac=0.18, 
#     S_ch4=0.05, S_IC=0.024, S_IN=3.2, S_IP=0.005, S_I=6.05, 
#     X_ch=1.15, X_pr=1.15, X_li=1.73, X_su=0.85, X_aa=0.6, X_fa=0.7, X_c4=0.3, 
#     X_pro=0.135, X_ac=0.9, X_h2=0.435, X_I=46.1, 
#     S_Ca=0.01, S_Mg=0.05, S_K=0.02, 
#     X_CaCO3=1e-8, X_struv=1e-8, X_newb=1e-8, X_ACP=1e-8, 
#     X_MgCO3=1e-8, X_AlOH=1e-8, X_AlPO4=1e-8, X_FeOH=1e-8, 
#     X_FePO4=1e-8, 
#     )

# default_ad_init = {k:v*1e3 for k,v in default_ad_init.items()}

# %%

@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    
    sys.simulate(
        # state_reset_hook='reset_cache',
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
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    # biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[WAS],
    #               active_unit_IDs=('A1', 'A2', 'A3', 'A4', 'O1', 'O2'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))