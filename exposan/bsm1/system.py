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
    'biomass_IDs',
    'create_system',
    'default_asm1_kwargs', 'default_asm2d_kwargs',
    'default_asm1_inf_kwargs', 'default_asm2d_inf_kwargs',
    'default_asm1_init_conds', 'default_asm2d_init_conds',
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
biomass_IDs = ('X_BH', 'X_BA')

# aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa_20=240, SOTE=0.3, V=V_ae,
#                           T_air=Temp, T_water=Temp, d_submergence=4-0.3)
# aer.A = 8.10765
# aer.B = 1750.286
# aer.C = 235.0

default_asm1_inf_kwargs = {
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

default_asm2d_inf_kwargs = {
    'concentrations': {
      'S_I': 14,
      'X_I': 26.5,
      'S_F': 20.1,
      'S_A': 94.3,
      'X_S': 409.75,
      'S_NH4': 31,
      'S_N2': 0,
      'S_NO3': 0.266, 
      'S_PO4': 2.8,
      'X_PP': 0.05,
      'X_PHA': 0.5,
      'X_H': 0.15,
      'X_AUT': 0, 
      'X_PAO': 0, 
      'S_ALK':7*12,
        },
    'units': ('m3/d', 'mg/L'),
    }

default_asm1_kwargs = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
    eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )

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
            )

default_asm1_init_conds = {
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

default_asm2d_init_conds = {
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


def batch_init(sys, df):
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

def create_system(
        flowsheet=None, 
        suspended_growth_model='ASM1',
        inf_kwargs={},
        asm_kwargs={},
        init_conds=None,
        aeration_processes=(),
        ):
    '''
    Create the system Benchmark Simulation Model No.1.
    
    Parameters
    ----------
    flowsheet : obj
        Flowsheet where this system will be created on.
    suspended_growth_model : str
        Either "ASM1" using Activated Sludge Model No. 1,
        or "ASM2d" using Activated Sludge Model No. 2d.
    inf_kwargs : dict
        Keyword arguments for influent.
    asm_kwargs : dict
        Keyword arguments for the ASM model (ASM1 or ASM2d).
    init_conds : dict or DataFrame
        For a dict, keyword arguments for initial conditions for all bioreactors in the system
        (the same initial conditions will be used),
        or a pandas.DataFrame that contains initial conditions for each unit.
        Default initial conditions will be used if not given.
    '''
    flowsheet = flowsheet or qs.Flowsheet('bsm1')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    kind = suspended_growth_model.lower().replace('-', '').replace('_', '')
    if kind == 'asm1':
        cmps_f = pc.create_asm1_cmps
        pc_f = pc.ASM1
        default_inf_kwargs = default_asm1_inf_kwargs
        default_asm_kwargs = default_asm1_kwargs
        DO_ID = 'S_O'
    elif kind == 'asm2d':
        cmps_f = pc.create_asm2d_cmps
        pc_f = pc.ASM2d
        default_inf_kwargs = default_asm2d_inf_kwargs
        default_asm_kwargs = default_asm2d_kwargs
        DO_ID = 'S_O2'
    else: raise ValueError('`suspended_growth_model` can only be "ASM1" or "ASM2d", '
                           f'not {suspended_growth_model}.')
    
    cmps_f()
    
    wastewater = WasteStream('wastewater', T=Temp)
    inf_kwargs = inf_kwargs or default_inf_kwargs
    wastewater.set_flow_by_concentration(Q, **inf_kwargs)
    
    effluent = WasteStream('effluent', T=Temp)
    WAS = WasteStream('WAS', T=Temp)
    RWW = WasteStream('RWW', T=Temp)
    RAS = WasteStream('RAS', T=Temp)
    
    # Process models
    if aeration_processes:
        aer1, aer2, aer3 = aeration_processes
    else:
        aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=240, DOsat=8.0, V=V_ae)
        aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=84, DOsat=8.0, V=V_ae)
    asm_kwargs = asm_kwargs or default_asm_kwargs
    asm = pc_f(**asm_kwargs)
    
    # Create unit operations
    A1 = su.CSTR('A1', ins=[wastewater, RWW, RAS], V_max=V_an,
                  aeration=None, suspended_growth_model=asm)
    
    A2 = su.CSTR('A2', A1-0, V_max=V_an,
                  aeration=None, suspended_growth_model=asm)
    
    O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                  DO_ID=DO_ID, suspended_growth_model=asm)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                  DO_ID=DO_ID, suspended_growth_model=asm)
    
    O3 = su.CSTR('O3', O2-0, [RWW, 'treated'], split=[0.6, 0.4],
                 V_max=V_ae, aeration=aer3,
                  DO_ID=DO_ID, suspended_growth_model=asm)
    
    C1 = su.FlatBottomCircularClarifier('C1', O3-1, [effluent, RAS, WAS],
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

    # S1 = su.Sampler('S1', C1-0, effluent)
    
    # System setup
    sys = System('bsm1_sys', path=(A1, A2, O1, O2, O3, C1), recycle=(RWW, RAS))
    # sys = System('bsm1_sys', path=(A1, A2, O1, O2, O3, C1, S1), recycle=(RWW, RAS))
    # bio = System('bsm1_sys_bio', path=(A1, A2, O1, O2, O3), recycle=(RWW,))
    # sys = System('bsm1_sys', path=(bio, C1, S1), recycle=(RAS,))

    if init_conds:
        if type(init_conds) is dict:
            for i in [A1, A2, O1, O2, O3]: i.set_init_conc(**init_conds)
        else:
            df = init_conds
    else: 
        path = os.path.join(data_path, f'initial_conditions_{kind}.xlsx')
        df = load_data(path, sheet='default')
        batch_init(sys, df)
    sys.set_dynamic_tracker(A1, effluent)
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
    srt = get_SRT(sys, biomass_IDs)
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