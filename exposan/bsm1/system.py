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

import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, WasteStream, System
from qsdsan.utils import time_printer, load_data, get_SRT
from exposan.bsm1 import data_path

__all__ = (
    'bio_IDs',
    'create_system',
    'default_asm_kwargs', 'default_inf_kwargs', 'default_init_conds',
    'Q', 'Q_ras', 'Q_was', 'Temp', 'V_an', 'V_ae', 
    )


# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q = 18446           # influent flowrate [m3/d]
Temp = 273.15+20    # temperature [K]

V_an = 1000    # anoxic zone tank volume
V_ae = 1333    # aerated zone tank volume
Q_was = 385    # sludge wastage flowrate
Q_ras = 18446    # recycle sludge flowrate
bio_IDs = ('X_BH', 'X_BA')

# aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa_20=240, SOTE=0.3, V=V_ae,
#                           T_air=Temp, T_water=Temp, d_submergence=4-0.3)
# aer.A = 8.10765
# aer.B = 1750.286
# aer.C = 235.0

default_inf_kwargs = {
    'concentrations': {
        'S_S':69.5,
        'X_BH':28.17,
        'X_S':202.32,
        'X_I':51.2,
        'S_NH':31.56,
        'S_I':30,
        'S_ND':6.95,
        'X_ND':10.59,
        'S_ALK':7*12,
        },
    'units': ('m3/d', 'mg/L'),
    }

default_asm_kwargs = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
    eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )

default_init_conds = {
        'S_S':5,
        'X_I':1000,
        'X_S':100,
        'X_BH':500,
        'X_BA':100,
        'X_P':100,
        'S_O':2,
        'S_NO':20,
        'S_NH':2,
        'S_ND':1,
        'X_ND':1,
        'S_ALK':7*12,
    }

def batch_init(sys, path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in [u.A1, u.A2, u.O1, u.O2, u.O3]:
        k.set_init_conc(**dct[k._ID])
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

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    flowsheet = flowsheet or qs.Flowsheet('bsm1')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    pc.create_asm1_cmps()
    
    PE = WasteStream('Wastewater', T=Temp)
    inf_kwargs = inf_kwargs or default_inf_kwargs
    PE.set_flow_by_concentration(Q, **inf_kwargs)
    
    SE = WasteStream('Effluent', T=Temp)
    WAS = WasteStream('WAS', T=Temp)
    RE = WasteStream('RWW', T=Temp)
    RAS = WasteStream('RAS', T=Temp)
    
    # Process models
    if aeration_processes:
        aer1, aer2, aer3 = aeration_processes
    else:
        aer1 = aer2 = pc.DiffusedAeration('aer1', 'S_O', KLa=240, DOsat=8.0, V=V_ae)
        aer3 = pc.DiffusedAeration('aer3', 'S_O', KLa=84, DOsat=8.0, V=V_ae)
    asm_kwargs = asm_kwargs or default_asm_kwargs
    asm1 = pc.ASM1(**asm_kwargs)
    
    # Create unit operations
    A1 = su.CSTR('A1', ins=[PE, RE, RAS], V_max=V_an,
                  aeration=None, suspended_growth_model=asm1)
    
    A2 = su.CSTR('A2', A1-0, V_max=V_an,
                  aeration=None, suspended_growth_model=asm1)
    
    O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                  DO_ID='S_O', suspended_growth_model=asm1)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                  DO_ID='S_O', suspended_growth_model=asm1)
    
    O3 = su.CSTR('O3', O2-0, [RE, 'treated'], split=[0.6, 0.4],
                 V_max=V_ae, aeration=aer3,
                  DO_ID='S_O', suspended_growth_model=asm1)
    
    C1 = su.FlatBottomCircularClarifier('C1', O3-1, [SE, RAS, WAS],
                                        underflow=Q_ras, wastage=Q_was, surface_area=1500,
                                        height=4, N_layer=10, feed_layer=5,
                                        X_threshold=3000, v_max=474, v_max_practical=250,
                                        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)

    # # Legacy codes for debugging the `Sampler`    
    # C1 = su.FlatBottomCircularClarifier('C1', O3-1, ['', RAS, WAS],
    #                                     underflow=Q_ras, wastage=Q_was, surface_area=1500,
    #                                     height=4, N_layer=10, feed_layer=5,
    #                                     X_threshold=3000, v_max=474, v_max_practical=250,
    #                                     rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)

    # S1 = su.Sampler('S1', C1-0, SE)
    
    # System setup
    sys = System('sys', path=(A1, A2, O1, O2, O3, C1), recycle=(RE, RAS))
    # sys = System('sys', path=(A1, A2, O1, O2, O3, C1, S1), recycle=(RE, RAS))
    # bio = System('bio', path=(A1, A2, O1, O2, O3), recycle=(RE,))
    # sys = System('sys', path=(bio, C1, S1), recycle=(RAS,))

    if init_conds:
        for i in [A1, A2, O1, O2, O3]: i.set_init_conc(**init_conds)
    else: batch_init(sys, os.path.join(data_path, 'initial_conditions.xlsx'), 'default')
    sys.set_dynamic_tracker(A1, SE)
    sys.set_tolerance(rmol=1e-6)
    
    return sys


#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    sys = create_system()
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    srt = get_SRT(sys, bio_IDs)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

if __name__ == '__main__':
    t = 50
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
    run(t, t_step, method=method)