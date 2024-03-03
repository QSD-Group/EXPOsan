# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Martin Xu <martinx3@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, WasteStream, System
from qsdsan.utils import time_printer, load_data, get_SRT, ospath
from exposan.mecha import data_path

__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm_kwargs',
    'default_inf_kwargs',
    'default_init_conds',
    'Q_ras', 'Q_was', 'Temp', 'V_an', 'V_fa', 'V_ae',
    )

# %%
# =============================================================================
# Parameters and util functions
# =============================================================================

# Q = 117706.01           # influent flowrate [m3/d]
                          # average from 'train_val_test_dynamic_influent': 54033.37

Temp = 273.15 + 18.36     # temperature [K]
                          # average from 'train_val_test_online_dataset': 18.36

V_an = 10920    # anoxic tank volume [m3]
V_fa = 10100    # facultative tank volume [m3]
V_ae = 11179    # aeration tank volume [m3]

Q_int = 9800*24    # internal recycle from aerobic to anoxic zone flowrate [m3/d]
Q_was = 30.5*24    # sludge wastage flowrate [m3/d]
Q_ras = 1600*24    # recycle sludge flowrate [m3/d]

biomass_IDs = ('X_BH', 'X_BA')

valid_models = ('asm1', 'asm2d')
default_inf_kwargs = dict.fromkeys(valid_models)
default_inf_kwargs['asm1'] = {
    'concentrations': {
        'X_BH':11.832,
        'X_I':17.728,
        'X_ND':14.28,
        'X_P':0.01,
        'X_S':88.74,
        'S_ALK':30,
        'S_I':73.695,
        'S_ND':9.52,
        'S_NH':44.2,
        'S_NO':0.01,
        'S_O':0.30,
        'S_S':221.085,
        'X_BA':0.01,
        },
    'units': ('m3/d', 'mg/L'),
    }

# default_inf_kwargs ['asm2d'] = {
#     'concentrations': {
#       'S_I': 14,
#       'X_I': 26.5,
#       'S_F': 20.1,
#       'S_A': 94.3,
#       'X_S': 409.75,
#       'S_NH4': 31,
#       'S_N2': 0,
#       'S_NO3': 0.266,
#       'S_PO4': 2.8,
#       'X_PP': 0.05,
#       'X_PHA': 0.5,
#       'X_H': 0.15,
#       'X_AUT': 0,
#       'X_PAO': 0,
#       'S_ALK':7*12,
#       },
#     'units': ('m3/d', 'mg/L'),
#     }

default_asm_kwargs = dict.fromkeys(valid_models)
default_asm_kwargs['asm1'] = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=6.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.62,
    eta_g=0.8, eta_h=0.4, k_h=3.0, K_X=0.03, mu_A=0.8,
    K_NH=1.0, b_A=0.1, K_O_A=0.4, k_a=0.08, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )

# default_asm_kwargs['asm1'] = dict(
#     Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#     mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
#     eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
#     K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
#     path=os.path.join(data_path, '_asm1.tsv'),
#     )
################################# Kaggle provided #################################

# # Stoichiometric parameters at 20 degrees celsius
# Y_A = 0.24  # autotrophic yield (default 0.24 gCOD/gCOD)
# Y_H = 0.67  # heterotrophic yield (default 0.67 gCOD/gCOD)
# f_P = 0.08  # fraction of biomass leading to particulate material (default 0.08 -)
# i_XB = 0.08 # nitrogen fraction in biomass (default 0.086 gN/gCOD)
# i_XP = 0.06 # nitrogen fraction in endogenous mass (default 0.01 gN/gCOD)
# # Kinetic parameters
# mu_H = 6.0  # maximum specific growth rate (default 6.0 1/day)                             4.0
# K_S = 20.0  # substrate saturation constant (default 20.0 gCOD/m3)                         10.0
# K_OH = 0.2  # Oxygen saturation constant (default 0.2 gO2/m3)
# K_NO = 0.5  # nitrate saturation constant (default 0.5 gNO3-N/m3)
# b_H = 0.62  # specific decay rate (default 0.62 1/day)                                     0.3
# eta_g = 0.8 # anoxic growth correction factor (default 0.8 -)
# eta_h = 0.4 # anoxic hydrolysis correction factor (default 0.4 -)                          0.8
# k_h = 3.0   # maximum specific hydrolysis rate (default 3.0 1/day)
# K_X = 0.03  # half-saturation coefficient for hydrolysis of XS (default 0.03 -)            0.1
# mu_A = 0.8  # maximum specific growth rate (default 0.8 1/day)                             0.5
# K_NH = 1.0  # ammonium saturation constant (default 1.0 gO2/m3)
# b_A = 0.1   # specific decay rate (default 0.1 1/day)                                      0.05
# K_OA = 0.4  # oxygen saturation constant (default 0.4 gO2/m3)
# k_a = 0.08  # ammonification rate constant (default 0.08 m3/(gCOD*day))                    0.05

###################################################################################

# default_asm_kwargs['asm2d'] = dict(
#     iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
#     iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
#     iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
#     f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
#     Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
#     Y_A=0.24, f_XI_AUT=0.1,
#     K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
#     mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
#     K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
#     q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
#     b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
#     K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
#     K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
#     mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
#     k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5,
#     )

default_init_conds = dict.fromkeys(valid_models)
# default_init_conds['asm1'] = {
#     'S_I':30,
#     'S_S':69.5,
#     'X_I':51.2,
#     'X_S':202.32,
#     'X_BH':28.17,
#     'X_BA':0,
#     'X_P':0,
#     'S_O':0,
#     'S_NO':0,
#     'S_NH':31.56,
#     'S_ND':6.95,
#     'X_ND':10.59,
#     'S_ALK':7*12,
#     }              # same as initial condition excel file

# default_init_conds['asm1'] = {
#     'S_I':30,
#     'S_S':69.5,
#     'X_I':51.2,
#     'X_S':202.32,
#     'X_BH':28.17,
#     'X_BA':0,
#     'X_P':0,
#     'S_O':0,
#     'S_NO':0,
#     'S_NH':0.23,         # lower initial NH from 'train_val_test_online_dataset'
#     'S_ND':6.95,
#     'X_ND':10.59,
#     'S_ALK':7*12,
#     }

default_init_conds['asm1'] = {
    'S_I':30,
    'S_S':69.5,
    'X_I':1000,
    'X_S':202.32,
    'X_BH':500,
    'X_BA':100,
    'X_P':100,
    'S_O':0,
    'S_NO':1,
    'S_NH':0.23,         # initial condition: New
    'S_ND':6.95,
    'X_ND':10.59,
    'S_ALK':7*12,
    }

# default_init_conds['asm2d'] = {
#     'S_F':5,
#     'S_A':2,
#     'X_I':1000,
#     'X_S':100,
#     'X_H':500,
#     'X_AUT':100,
#     #'X_P':100,
#     'S_O2':2,
#     'S_NO3':20,
#     'S_NH4':2,
#     'S_ALK':7*12,
#     }


def batch_init(sys, df):
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in [u.ANO, u.FAC, u.AER1, u.AER2, u.AER3]:
        k.set_init_conc(**dct[k._ID])
    c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    tss = [v for v in dct['C1_tss'].values() if v>0]
    u.C1.set_init_solubles(**c1s)
    u.C1.set_init_sludge_solids(**c1x)
    u.C1.set_init_TSS(tss)


#%%

# =============================================================================
# Tilburg WWTP
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
    Create the system as described in Benchmark Simulation Model No.1.

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
    flowsheet = flowsheet or qs.Flowsheet('tilburg')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    kind = suspended_growth_model.lower().replace('-', '').replace('_', '')
    if kind == 'asm1':
        pc.create_asm1_cmps()
        asm = pc.ASM1(**default_asm_kwargs[kind])
        DO_ID = 'S_O'
    elif kind == 'asm2d':
        pc.create_asm2d_cmps()
        asm = pc.ASM2d(**default_asm_kwargs[kind])
        DO_ID = 'S_O2'
    else: raise ValueError('`suspended_growth_model` can only be "ASM1" or "ASM2d", '
                           f'not {suspended_growth_model}.')

    # wastewater = WasteStream('wastewater', T=Temp)
    # inf_kwargs = inf_kwargs or default_inf_kwargs[kind]
    # wastewater.set_flow_by_concentration(Q, **inf_kwargs)

    DYINF = WasteStream('Dynamic_influent', T=Temp)
    effluent = WasteStream('effluent', T=Temp)
    WAS = WasteStream('WAS', T=Temp)
    INT = WasteStream('INT', T=Temp)
    RAS = WasteStream('RAS', T=Temp)

    # Process models
    if aeration_processes:
        aer = aeration_processes
    else:
        aer1 = pc.DiffusedAeration('aer1', DO_ID, KLa=20, DOsat=8.0, V=V_ae)                  #  need to think about KLa

        aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=70, DOsat=8.0, V=V_ae)                  #  need to think about KLa
        # aer2 = pc.DiffusedAeration('aer2', DO_ID, KLa=240, DOsat=8.0, V=V_ae)                  #  need to think about KLa
        # # aerFac = pc.DiffusedAeration('aerFac', DO_ID, DOsat=8.0, V=V_fa, Q_air=576.1015981*24)
        # aer1 = pc.DiffusedAeration('aer1', DO_ID, DOsat=8.0, V=V_ae, Q_air=764.3221021*24)                  #  need to think about KLa
        # aer2 = pc.DiffusedAeration('aer2', DO_ID, DOsat=8.0, V=V_ae, Q_air=1617.996689*24)                  #  need to think about KLa
        # aer3 = pc.DiffusedAeration('aer3', DO_ID, DOsat=8.0, V=V_ae, Q_air=578.9891202*24)

    # Create unit operations
    WW = su.DynamicInfluent('Waste_Water', outs=[DYINF],
                            data_file=os.path.join(data_path, 'dynamic_influent_train_val_test.tsv'))

    ANO = su.CSTR('ANO', ins=[DYINF, INT, RAS], V_max=V_an,
                 aeration=None, suspended_growth_model=asm)


    FAC = su.CSTR('FAC', ANO-0, V_max=V_fa, aeration=aer1,
                  DO_ID=DO_ID, suspended_growth_model=asm)
    # FAC = su.CSTR('FAC', ANO-0, V_max=V_fa, aeration=None,
    #               DO_ID=DO_ID, suspended_growth_model=asm)
    # FAC = su.CSTR('FAC', ANO-0, V_max=V_fa, aeration=aer1,
    #               DO_ID=DO_ID, suspended_growth_model=asm)
    # FAC = su.CSTR('FAC', ANO-0, V_max=V_fa, aeration=aerFac,
    #               DO_ID=DO_ID, suspended_growth_model=asm)

    AER1 = su.CSTR('AER1', FAC-0, V_max=V_ae, aeration=aer2,
                   DO_ID=DO_ID, suspended_growth_model=asm)
    # AER1 = su.CSTR('AER1', FAC-0, V_max=V_ae, aeration=None,
    #                DO_ID=DO_ID, suspended_growth_model=asm)
    # AER1 = su.CSTR('AER1', FAC-0, V_max=V_ae, aeration=aer2,
    #               DO_ID=DO_ID, suspended_growth_model=asm)
    # AER1 = su.CSTR('AER1', FAC-0, V_max=V_ae, aeration=aer1,
    #             DO_ID=DO_ID, suspended_growth_model=asm)

    AER2 = su.CSTR('AER2', AER1-0, V_max=V_ae, aeration=aer2,
                   DO_ID=DO_ID, suspended_growth_model=asm)
    # AER2 = su.CSTR('AER2', AER1-0, V_max=V_ae, aeration=None,
    #                DO_ID=DO_ID, suspended_growth_model=asm)
    # AER2 = su.CSTR('AER2', AER1-0, V_max=V_ae, aeration=aer2,
    #               DO_ID=DO_ID, suspended_growth_model=asm)
    # AER2 = su.CSTR('AER2', AER1-0, V_max=V_ae, aeration=aer2,
    #               DO_ID=DO_ID, suspended_growth_model=asm)

    AER3 = su.CSTR('AER3', AER2-0, [INT, 'treated'], split=[0.601,0.399],                             # split ratio
                  V_max=V_ae, aeration=aer2,
                  DO_ID=DO_ID, suspended_growth_model=asm)
    # AER3 = su.CSTR('AER3', AER2-0, [INT, 'treated'], split=[0.601,0.399],                             # split ratio
    #               V_max=V_ae, aeration=None,
    #               DO_ID=DO_ID, suspended_growth_model=asm)
    # AER3 = su.CSTR('AER3', AER2-0, [INT, 'treated'], split=[0.601,0.399],                             # split ratio
    #               V_max=V_ae, aeration=aer2,
    #               DO_ID=DO_ID, suspended_growth_model=asm)
    # AER3 = su.CSTR('AER3', AER2-0, [INT, 'treated'], split=[0.601,0.399],                               # split ratio
    #               V_max=V_ae, aeration=aer3,
    #               DO_ID=DO_ID, suspended_growth_model=asm)

    C1 = su.FlatBottomCircularClarifier('C1', AER3-1, [effluent, RAS, WAS],
                                        underflow=Q_ras, wastage=Q_was, surface_area=1700,
                                        height=2.3, N_layer=10, feed_layer=5,
                                        X_threshold=3000, v_max=474, v_max_practical=250,
                                        rh=0.00019, rp=0.0028, fns=0.001)

    # C1 = su.FlatBottomCircularClarifier('C1', AER3-1, [effluent, RAS, WAS],
    #                                 underflow=Q_ras, wastage=Q_was, surface_area=1500,
    #                                 height=4, N_layer=10, feed_layer=5,
    #                                 X_threshold=3000, v_max=474, v_max_practical=250,
    #                                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)

    # num_layers=10,
    # Xf=(lambda x: 0.75*sum(x[2:7])),  # Formula to calculate TSS from the states according to ASM1
    # feedlayer=5,  # Stating to count at 0, thus this is the 6th layer
    # A=1700, z=2.3, v0_max=250.0, v0=474.0, rh=0.00019, rp=0.0028, fns=0.001, Xt=3000.0,

    # System setup
    sys = System('tilburg', path=(WW, ANO, FAC, AER1, AER2, AER3, C1), recycle=(INT, RAS))

    if init_conds:
        if type(init_conds) is dict:
            for i in [ANO, FAC, AER1, AER2, AER3]: i.set_init_conc(**init_conds)
        else:
            df = init_conds
    else:
        path = os.path.join(data_path, f'initial_conditions_{kind}_new.xlsx')
        # path = os.path.join(data_path, f'initial_conditions_{kind}_lownh.xlsx')
        df = load_data(path, sheet='default')
        batch_init(sys, df)
    # sys.set_dynamic_tracker(effluent)
    sys.set_dynamic_tracker(WW, ANO, FAC, AER1, AER2, AER3, C1, effluent)

    sys.set_tolerance(rmol=1e-6)

    return sys

#%%

# t = 100
# t_step = 1
# # method = 'RK45'
# method = 'RK23'
# # method = 'DOP853'
# # method = 'Radau'
# # method = 'BDF'

# sys = create_system()
# sys.simulate(
#     state_reset_hook='reset_cache',
#     t_span=(0,t),
#     t_eval=np.arange(0, t+t_step, t_step),
#     method=method,
#     # rtol=1e-2,
#     # atol=1e-3,
#     print_t=True,
#     export_state_to=f'results/sol_{t}d_{method}.xlsx',
#     )

# #%%

# answer = load_data(ospath.join(data_path, 'train_val_test_online_dataset_.xlsx'), sheet='train_val_test_online_dataset')

# from datetime import datetime
# day0 = datetime.strptime('5-5-2021 11:57 AM', '%m-%d-%Y %I:%M %p') # set the first time stamp as 'day0'
# calc_day = lambda t_stamp: (t_stamp-day0).total_seconds()/60/60/24

# t_stamp = answer.index
# t_intp = calc_day(t_stamp).to_numpy()

# import matplotlib.pyplot as plt, matplotlib.ticker as ticker

# def plot_raw_vs_smooth(x1, y1, x2, y2):
#     fig, ax = plt.subplots(figsize=(15, 3)) # figure size
#     l1, = ax.plot(x1, y1, label='Model') # first plot = raw data
#     l2, = ax.plot(x2, y2, label='Answer') # second plot = smoothed result
#     ax.xaxis.set_major_locator(ticker.MultipleLocator(30)) # x-axis tick interval
#     ax.set_xlabel('Time [day]') # x-axis title
#     ax.legend(handles=[l1,l2])
#     return fig, ax

# # plot_raw_vs_smooth(t_intp, answer['NH4'],
# #                     sys.units.AER3.scope.time_series, sys.units.AER3.scope.record[:,sys.units[1].components.indices(['S_NH'])])

# # plot_raw_vs_smooth(t_intp, answer['NO3'],
# #                     sys.units.AER3.scope.time_series, sys.units.AER3.scope.record[:,sys.units[1].components.indices(['S_NO'])])

# # plot_raw_vs_smooth(t_intp, answer['DO_1'],
# #                     sys.units.AER3.scope.time_series, sys.units.AER3.scope.record[:,sys.units[1].components.indices(['S_O'])])

# plot_raw_vs_smooth(t_intp, answer['NH4'],
#                     sys.flowsheet.stream.effluent.scope.time_series, sys.flowsheet.stream.effluent.scope.record[:,sys.units[1].components.indices(['S_NH'])])

# plot_raw_vs_smooth(t_intp, answer['NO3'],
#                     sys.flowsheet.stream.effluent.scope.time_series, sys.flowsheet.stream.effluent.scope.record[:,sys.units[1].components.indices(['S_NO'])])

# plot_raw_vs_smooth(t_intp, answer['DO_1'],
#                     sys.flowsheet.stream.effluent.scope.time_series, sys.flowsheet.stream.effluent.scope.record[:,sys.units[1].components.indices(['S_O'])])

# # sys.units.AER3.scope.plot_time_series(('S_I'))
# # sys.units.AER3.scope.plot_time_series(('S_S'))
# # sys.units.AER3.scope.plot_time_series(('X_I'))
# # sys.units.AER3.scope.plot_time_series(('X_S'))
# # sys.units.AER3.scope.plot_time_series(('X_BH'))
# # sys.units.AER3.scope.plot_time_series(('X_BA'))
# # sys.units.AER3.scope.plot_time_series(('X_P'))
# # sys.units.AER3.scope.plot_time_series(('S_O'))
# # sys.units.AER3.scope.plot_time_series(('S_NO'))
# # sys.units.AER3.scope.plot_time_series(('S_NH'))
# # sys.units.AER3.scope.plot_time_series(('S_ND'))
# # sys.units.AER3.scope.plot_time_series(('X_ND'))
# # sys.units.AER3.scope.plot_time_series(('S_ALK'))

# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_I'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_S'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('X_I'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('X_S'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('X_BH'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('X_BA'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('X_P'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_O'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_NO'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_NH'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_ND'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('X_ND'))
# sys.flowsheet.stream.effluent.scope.plot_time_series(('S_ALK'))
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
        export_state_to=f'results/sol_{t}d_{method}.xlsx',
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