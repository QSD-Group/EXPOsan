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
    'default_asm_kwargs',
    'default_inf_kwargs',
    'default_init_conds',
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
Q_ras = 1*Q    # recycle sludge flowrate
Q_intr = 3*Q    # internal recirculation

# aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa_20=240, SOTE=0.3, V=V_ae,
#                           T_air=Temp, T_water=Temp, d_submergence=4-0.3)
# aer.A = 8.10765
# aer.B = 1750.286
# aer.C = 235.0

valid_models = ('asm1', 'asm2d')
biomass_IDs = dict.fromkeys(valid_models)
biomass_IDs['asm1'] = ('X_BH', 'X_BA')
biomass_IDs['asm2d'] = ('X_H', 'X_PAO', 'X_AUT')

default_inf_kwargs = dict.fromkeys(valid_models)
default_inf_kwargs['asm1'] = {
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

default_inf_kwargs['asm2d'] = {
    # 'concentrations': {
    #   'S_I': 14,
    #   'X_I': 26.5,
    #   'S_F': 20.1,
    #   'S_A': 94.3,
    #   'X_S': 409.75,
    #   'S_NH4': 31,
    #   'S_NO3': 0.266, 
    #   'S_PO4': 2.8,
    #   'X_PP': 0.05,
    #   'X_PHA': 0.5,
    #   'X_H': 0.15,
    #   'S_ALK':7*12,
    #   },
    'concentrations': dict(
        S_F=41.7,
        S_A=27.8,
        S_I=30.0,
        S_NH4=40.04,
        S_PO4=9.01,
        S_ALK=7.0*12,
        X_I=51.2,
        X_S=202.32,
        X_H=28.17,
        ),
    'units': ('m3/d', 'mg/L'),
    }


default_asm_kwargs = dict.fromkeys(valid_models)
default_asm_kwargs['asm1'] = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
    eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )

default_asm_kwargs['asm2d'] = dict(
    iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
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

default_c1_kwargs = dict(
    underflow=Q_ras, wastage=Q_was, surface_area=1500,
    height=4, N_layer=10, feed_layer=5,
    X_threshold=3000, v_max=474, v_max_practical=250,
    rh=5.76e-4, rp=2.86e-3, fns=2.28e-3
    )

default_init_conds = dict.fromkeys(valid_models)
default_init_conds['asm1'] = {
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

default_init_conds['asm2d'] = {
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
    for k in u:
        if k.ID == 'C1': continue
        elif k.ID == 'AS':
            k.set_init_conc(concentrations=df.iloc[:-3])
        else:
            k.set_init_conc(**dct[k.ID])
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
        reactor_model='CSTR',
        inf_kwargs={},
        asm_kwargs={},
        settler_kwargs={},
        init_conds=None,
        aeration_processes=(),
        ):
    '''
    Create the system as described in Benchmark Simulation Model No.1.
    
    Parameters
    ----------
    flowsheet : obj
        Flowsheet where this system will be created on.
    suspended_growth_model : str
        Either "ASM1" using Activated Sludge Model No. 1,
        or "ASM2d" using Activated Sludge Model No. 2d.
    reactor_model : str
        "CSTR" to model each zone in the activated sludge reactor as a CSTR and 
        model the internal reciruclation explicitly; "PFR" to model the entire
        activated sludge reactor as a single unit, with implicit internal recirculation.
    inf_kwargs : dict
        Keyword arguments for influent.
    asm_kwargs : dict
        Keyword arguments for the ASM model (ASM1 or ASM2d).
    settler_kwargs : dict
        Keyword arguments for the clarifier.
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
    asm_kwargs = asm_kwargs or default_asm_kwargs[kind]
    if kind == 'asm1':
        pc.create_asm1_cmps()
        asm = pc.ASM1(**asm_kwargs)
        DO_ID = 'S_O'
    elif kind == 'asm2d':
        pc.create_asm2d_cmps()
        asm = pc.ASM2d(**asm_kwargs)
        DO_ID = 'S_O2'
    else: 
        raise ValueError('`suspended_growth_model` can only be "ASM1" or "ASM2d", '
                           f'not {suspended_growth_model}.')
    
    wastewater = WasteStream('wastewater', T=Temp)
    inf_kwargs = inf_kwargs or default_inf_kwargs[kind]
    wastewater.set_flow_by_concentration(Q, **inf_kwargs)
    
    # Process models
    if aeration_processes:
        aer1, aer2, aer3 = aeration_processes
        kLa = [aer.KLa for aer in aeration_processes]
    else:
        aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=240, DOsat=8.0, V=V_ae)
        if kind == 'asm1':
            kLa = [240, 240, 84]
            aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=84, DOsat=8.0, V=V_ae)
        else:
            kLa = [240]*3
            aer3 = aer1
    
    # Create unit operations
    c1_kwargs = settler_kwargs or default_c1_kwargs
    if reactor_model == 'CSTR':
        an_kwargs = dict(V_max=V_an, aeration=None, suspended_growth_model=asm)
        ae_kwargs = dict(V_max=V_ae, DO_ID=DO_ID, suspended_growth_model=asm)
        if kind == 'asm1':
            A1 = su.CSTR('A1', ins=[wastewater, 'RWW', 'RAS'], **an_kwargs)       
            A2 = su.CSTR('A2', A1-0, **an_kwargs)        
            O1 = su.CSTR('O1', A2-0, aeration=aer1, **ae_kwargs)
            O2 = su.CSTR('O2', O1-0, aeration=aer2, **ae_kwargs)
            O3 = su.CSTR('O3', O2-0, [1-A1, 'treated'], split=[Q_intr, Q+Q_ras],
                         aeration=aer3, **ae_kwargs)
            C1 = su.FlatBottomCircularClarifier('C1', O3-1, 
                                                ['effluent', 2-A1, 'WAS'],
                                                **c1_kwargs)
            path=(A1, A2, O1, O2, O3, C1)
        else:
            A1 = su.CSTR('A1', ins=[wastewater, 'RAS'], **an_kwargs)       
            A2 = su.CSTR('A2', A1-0, **an_kwargs)        
            A3 = su.CSTR('A3', [A2-0, 'RWW'], **an_kwargs)        
            A4 = su.CSTR('A4', A3-0, **an_kwargs)        
            O1 = su.CSTR('O1', A4-0, aeration=aer1, **ae_kwargs)
            O2 = su.CSTR('O2', O1-0, aeration=aer2, **ae_kwargs)
            O3 = su.CSTR('O3', O2-0, [1-A3, 'treated'], split=[Q_intr, Q+Q_ras],
                         aeration=aer3, **ae_kwargs)
            C1 = su.FlatBottomCircularClarifier('C1', O3-1, 
                                                ['effluent', 1-A1, 'WAS'],
                                                **c1_kwargs)
            path = (A1, A2, A3, A4, O1, O2, O3, C1)
        sys = System('bsm1_sys', path=path, recycle=(O3-0, C1-1))
        sys.set_dynamic_tracker(A1, C1-0)
    elif reactor_model == 'PFR':
        as_kwargs = dict(
            DO_ID=DO_ID, DO_sat=8.0, suspended_growth_model=asm,
            gas_stripping=False
            )
        if kind == 'asm1':
            AS = su.PFR('AS', ins=[wastewater, 'RAS'], outs='treated', 
                        N_tanks_in_series=5,
                        V_tanks=[V_an]*2+[V_ae]*3,
                        influent_fractions=[[1]+[0]*4]*2,
                        internal_recycles=[(4,0,Q_intr)],
                        kLa=[0]*2+kLa, **as_kwargs)
        else:
            AS = su.PFR('AS', ins=[wastewater, 'RAS'], outs='treated', 
                        N_tanks_in_series=7,
                        V_tanks=[V_an]*4+[V_ae]*3,
                        influent_fractions=[[1]+[0]*6]*2,
                        internal_recycles=[(6,2,Q_intr)],
                        kLa=[0]*4+kLa, **as_kwargs)
        C1 = su.FlatBottomCircularClarifier('C1', AS-0, 
                                            ['effluent', 1-AS, 'WAS'],
                                            **c1_kwargs)
        sys = System('bsm1_sys', path=(AS, C1), recycle=(C1-1,))
        sys.set_dynamic_tracker(AS, C1-0)
    else:
        raise ValueError('`reactor_model` can only be "CSTR" or "PFR", '
                           f'not {reactor_model}.')

    sys.set_tolerance(rmol=1e-6)
    
    if init_conds:
        if isinstance(init_conds, dict):
            for i in sys.units: 
                if i.ID == 'C1': continue
                i.set_init_conc(**init_conds)
        else:
            df = init_conds
    else: 
        path = os.path.join(data_path, f'initial_conditions_{kind}.xlsx')
        df = load_data(path, sheet='default')
    batch_init(sys, df)

    return sys


#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    asm = 'ASM1'
    # asm = 'ASM2d'
    rxt = 'CSTR'
    # rxt = 'PFR'
    sys = create_system(suspended_growth_model=asm, reactor_model=rxt)
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    srt = get_SRT(sys, biomass_IDs[asm.lower()])
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    return sys

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
    sys = run(t, t_step, method=method)