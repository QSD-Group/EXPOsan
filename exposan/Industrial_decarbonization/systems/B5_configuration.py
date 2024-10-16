# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Saumitra Rai <raisaumitra9@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, WasteStream, System
from qsdsan.utils import time_printer, load_data, get_SRT
from exposan.bsm1 import data_path

__all__ = ('biomass_IDs','create_system', 'default_asm2d_kwargs', 'B5_asm2d_inf_kwargs', 
           'default_init_conds', 'Q', 'Q_was', 'Q_ras', 'Temp', 'V_ae',)


# %%

# =============================================================================
# Parameters and util functions
# =============================================================================
Q = 38000       # influent flowrate [m3/d]
Temp = 273.15 + 15.6    # temperature [K] (Based on WERF report)

V_ae = 3785.4  # aerated zone tank volume [m3] (Based on WERF report)
Q_was = 950    # sludge wastage flowrate
Q_ras = 1800      # recycle sludge flowrate
# biomass_IDs = ('X_BH', 'X_BA')
biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

# aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa_20=240, SOTE=0.3, V=V_ae,
#                           T_air=Temp, T_water=Temp, d_submergence=4-0.3)
# aer.A = 8.10765
# aer.B = 1750.286
# aer.C = 235.0
# =============================================================================
# 
# default_inf_kwargs = {
#     'concentrations': {
#         'S_S':69.5,
#         'X_BH':28.17,
#         'X_S':202.32,
#         'X_I':51.2,
#         'S_NH':31.56,
#         'S_I':30,
#         'S_ND':6.95,
#         'X_ND':10.59,
#         'S_ALK':7*12,
#         },
#     'units': ('m3/d', 'mg/L'),
#     }
# =============================================================================

B5_asm2d_inf_kwargs = {
    'concentrations': {
        'S_I': 30,
        'X_I': 0,
        'S_F': 0,
        'S_A': 31.504,
        'X_S': 237.054,
        'S_NH4': 25,
        'S_N2': 0,
        'S_NO3': 0, 
        'S_PO4': 5,
        'X_PP': 0,
        'X_PHA': 54.754,
        'X_H': 0,
        'X_AUT': 0, 
        'X_PAO': 0, 
        'S_ALK':7*12,
        },
    'units': ('m3/d', 'mg/L'),
    }

# =============================================================================
# default_asm_kwargs = dict(
#     Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#     mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
#     eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
#     K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
#     path=os.path.join(data_path, '_asm1.tsv'),
#     )
# 
# =============================================================================

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
            k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, path=os.path.join(data_path, '_asm2d.tsv'),
            )

default_init_conds = {
        'S_F':5,
        'S_A':2,
        'X_I':1000,
        'X_S':100,
        'X_H':500,
        'X_AUT':100,
        #'X_P':100,
        'S_O2':2,
        'S_NO3':20,
        'S_NH4':2,
        'S_ALK':7*12,
    }

def batch_init(sys, path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in sys.units:
        if k.ID.startswith('O'):
            k.set_init_conc(**dct['O'])
        elif k.ID.startswith('A'):
            k.set_init_conc(**dct['A'])           
    c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    tss = [v for v in dct['C1_tss'].values() if v>0]
    u.C1.set_init_solubles(**c1s)
    u.C1.set_init_sludge_solids(**c1x)
    u.C1.set_init_TSS(tss)


#%%

# =============================================================================
# Benchmark Simulation Model No. 1
# =============================================================================

def create_components():
     asm2d_cmps = pc.create_asm2d_cmps()
     CO2 = qs.Component.from_chemical('S_CO2', search_ID='CO2', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
     CH4 = qs.Component.from_chemical('S_CH4', search_ID='CH4', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
     H2 = qs.Component.from_chemical('S_H2', search_ID='H2', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
     cmps1 = qs.Components.load_default()
     ash = cmps1.X_Ig_ISS.copy('ash')
     cmps = qs.Components([*asm2d_cmps, CO2, CH4, H2, ash])
     cmps.compile()
     return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    flowsheet = flowsheet or qs.Flowsheet('bsm1')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    wastewater = WasteStream('wastewater', T=Temp)
    #inf_kwargs = inf_kwargs or default_inf_kwargs
    inf_kwargs = B5_asm2d_inf_kwargs 
    wastewater.set_flow_by_concentration(Q, **inf_kwargs)
    effluent = WasteStream('effluent', T=Temp)
    WAS = WasteStream('WAS', T=Temp)
    RAS = WasteStream('RAS', T=Temp)
    fuel = WasteStream('nat_gas', phase='g', S_CH4=1000)
    air = WasteStream('air', phase='g', S_O2=210, S_N2=780, S_H2=10)
    
    # Process models
    # if aeration_processes:
    #     aer1, aer2, aer3 = aeration_processes
    # else:
    #     aer1 = aer2 = pc.DiffusedAeration('aer1', 'S_O', KLa=240, DOsat=8.0, V=V_ae)
    #     aer3 = pc.DiffusedAeration('aer3', 'S_O', KLa=84, DOsat=8.0, V=V_ae)
    # # asm_kwargs = asm_kwargs or default_asm_kwargs
    asm2d = pc.ASM2d(iP_SF=0.005, iP_XS=0.0085, iP_XI=0.005, iN_XS=0.065)
    
    effluent_GT = WasteStream(ID = 'effluent_GT')
    effluent_MT = WasteStream(ID = 'effluent_MT')
    effluent_DU = WasteStream(ID = 'effluent_DU')
    
    PC = su.PrimaryClarifier(ID='PC', ins = [wastewater, effluent_GT, effluent_MT, effluent_DU], outs = ['sludge_PC', 'effluent_PC'], 
                        isdynamic=True, init_with='WasteStream', thickener_perc= 0.09, 
                        TSS_removal_perc=70)
   
    O1 = su.CSTR('O1', ins = [PC-1, RAS], V_max=V_ae, aeration= 2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration= 2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)
    
    O3 = su.CSTR('O3', O2-0, V_max=V_ae, aeration= 2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)
    
    O4 = su.CSTR('O4', O3-0, V_max=V_ae, aeration= 2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)
    
    O5 = su.CSTR('O5', O4-0, V_max=V_ae, aeration= 2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)
    
    O6 = su.CSTR('O6', O5-0, V_max=V_ae, aeration= 2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)

    C1 = su.FlatBottomCircularClarifier('C1', O6-0, [effluent, RAS, WAS],
    underflow=Q_ras, wastage=Q_was, surface_area=1500, height=4, N_layer=10, 
    feed_layer=5, X_threshold=3000, v_max=474, v_max_practical=250, rh=5.76e-4, 
    rp=2.86e-3, fns=2.28e-3)
    
    #Gravity thickener: 'thickener_perc', 'TSS_removal_perc' based on WERF report 
    GT = su.Thickener(ID='GT', ins = [PC-0], outs = ['sludge_GT', effluent_GT],
                      thickener_perc=7, TSS_removal_perc=92)
    
    #Mechanical thickener: 'thickener_perc', 'TSS_removal_perc' based on WERF report 
    MT = su.Thickener(ID='MT', ins = [WAS], outs = ['sludge_MT', effluent_MT],
                      thickener_perc=6, TSS_removal_perc=98)
    
    #Dewatering Unit: 'thickener_perc', 'TSS_removal_perc' based on WERF report 
    DU = su.Centrifuge(ID='DU', ins = [GT-0, MT-0], outs = ['cake', effluent_DU],
                           thickener_perc=23, TSS_removal_perc=95)
    
    # IC = su.Incinerator(ID ='IC', ins = [DU-0, air, fuel], outs=(), thermo=None, isdynamic=False, 
    #               init_with='WasteStream', F_BM_default=None, process_efficiency=0.90, 
    #               calorific_value_sludge= 12000, calorific_value_fuel=50000, carbon_di_oxide_ID = 'CO2',
    #               ash_component_ID = 'ash')

    sys = qs.System('test', path=(PC, O1, O2, O3, O4, O5, O6, C1, GT, MT, DU, 
                                  # IC
                                  ), 
                    recycle = [effluent_GT, effluent_MT, effluent_DU, RAS])

    return sys
#%%
@time_printer

# =============================================================================
# def create_components():
#      asm1_cmps = pc.create_asm1_cmps()
#      CO2 = qs.Component.from_chemical('S_CO2', search_ID='CO2', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
#      CH4 = qs.Component.from_chemical('S_CH4', search_ID='CH4', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
#      H2 = qs.Component.from_chemical('S_H2', search_ID='H2', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
#      cmps1 = qs.Components.load_default()
#      ash = cmps1.X_Ig_ISS.copy('ash')
#      cmps = qs.Components([*asm1_cmps, CO2, CH4, H2, ash])
#      cmps.compile()
#      return cmps
# 
# =============================================================================
def run(t, t_step, method=None, **kwargs):
    sys = create_system()
    # for u in sys.units:
    #     if u.ID in ('O1', 'A1', 'A2', 'A3', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9'):
    #         u.set_init_conc(**default_init_conds)
    batch_init(sys, "/Users/saumitrarai/Desktop/Important files_QSDsan/initial_conditions_ASM2d.xlsx", sheet='t=10')
    # RAS = sys.flowsheet.stream.RAS
    # C1 = sys.flowsheet.unit.C1
    # sys.set_dynamic_tracker(RAS, C1)
    sys.set_dynamic_tracker(*sys.products)
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    sys.diagram()
    return sys
    
if __name__ == '__main__':
    t = 1
    t_step = 1
    # method = 'RK45'
    # method = 'RK23' 
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    system = run(t, t_step, method=method)

    act_units = [u.ID for u in system.units if isinstance(u, (su.CSTR, su.FlatBottomCircularClarifier))]
    fs = system.flowsheet.stream
    srt = get_SRT(system, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

    #cmps = qs.get_components()
    #cmps = create_components()
    f = qs.main_flowsheet
    unit = system.flowsheet.unit
    fs = system.flowsheet.stream
    #fig, axis = fs.RAS.scope.plot_time_series(( 'S_NH')) 
    # fig, axis = fs.RAS.scope.plot_time_series(('S_S','S_NH')) 
    fig, axis = fs.effluent.scope.plot_time_series(('S_F', 'S_A', 'S_NH4', 'S_NO3')) 
    fig
# =============================================================================
#     fig, axis = unit.C1.scope.plot_time_series(('S_S', 'S_NH')) 
#     fig
# =============================================================================