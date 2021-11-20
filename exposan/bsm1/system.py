# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import os 
import numpy as np
from qsdsan import sanunits as su
from qsdsan import processes as pc
from qsdsan import set_thermo, WasteStream, System
from qsdsan.utils import time_printer
from exposan.bsm1 import bsm1_path


# =============================================================================
# Benchmark Simulation Model No. 1
# =============================================================================

############# load components and set thermo #############
cmps = pc.load_asm1_cmps()
set_thermo(cmps)

############# create WasteStream objects #################
Q = 18446.0           # influent flowrate [m3/d]
Temp = 273.15+20    # temperature [K]

PE = WasteStream('Wastewater', T=Temp)
PE.set_flow_by_concentration(Q,
                              {'S_S':69.5,
                              'X_BH':28.17,
                              'X_S':202.32,
                              'X_I':51.2,
                              'S_NH':31.56,
                              'S_I':30,
                              'S_ND':6.95,
                              'X_ND':10.59,
                              'S_ALK':7*12},
                              units=('m3/d', 'mg/L'))

SE = WasteStream('Effluent', T=Temp)
WAS = WasteStream('WAS', T=Temp)
RE = WasteStream('RWW', T=Temp)
RAS = WasteStream('RAS', T=Temp)

############# load and tailor process models #############
V_an = 1000    # anoxic zone tank volume
V_ae = 1333    # aerated zone tank volume
Q_was = 385.0    # sludge wastage flowrate
Q_ras = 18446.0    # recycle sludge flowrate

# pc.DiffusedAeration.A = 8.10765
# pc.DiffusedAeration.B = 1750.286
# pc.DiffusedAeration.C = 235.0
# aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa_20=240, SOTE=0.3, V=V_ae,
#                           T_air=Temp, T_water=Temp, d_submergence=4-0.3)
aer1 = pc.DiffusedAeration('aer1', 'S_O', KLa=240, DOsat=8.0, V=V_ae)
aer2 = pc.DiffusedAeration('aer2', 'S_O', KLa=84, DOsat=8.0, V=V_ae)
asm1 = pc.ASM1(Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
                mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
                eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
                K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
                path=os.path.join(bsm1_path, '_asm1.tsv'))

############# create unit operations #####################
A1 = su.CSTR('A1', ins=[PE, RE, RAS], V_max=V_an,
              aeration=None, suspended_growth_model=asm1)

A2 = su.CSTR('A2', A1-0, V_max=V_an,
              aeration=None, suspended_growth_model=asm1)

O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
              DO_ID='S_O', suspended_growth_model=asm1)

O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer1,
              DO_ID='S_O', suspended_growth_model=asm1)

O3 = su.CSTR('O3', O2-0, V_max=V_ae, aeration=aer2,
              DO_ID='S_O', suspended_growth_model=asm1)

S1 = su.Splitter('S1', O3-0, [RE, 'treated'], split=0.6, init_with='WasteStream')

C1 = su.FlatBottomCircularClarifier('C1', S1-1, [SE, 'sludge'],
                                    sludge_flow_rate=Q_ras+Q_was, surface_area=1500,
                                    height=4, N_layer=10, feed_layer=5,
                                    X_threshold=3000, v_max=474, v_max_practical=250,
                                    rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)

S2 = su.Splitter('S2', C1-1, [RAS, WAS], split=Q_ras/(Q_ras+Q_was), init_with='WasteStream')


A1.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
A2.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
O1.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
O2.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
O3.set_init_conc(S_I=30.0, S_S=0.88949, X_I=1149.1252, X_S=49.3056, X_BH=2559.3436, X_BA=149.7971, 
                  X_P=452.2111, S_O=0.49094, S_NH=1.7333, S_ND=0.68828, X_ND=3.5272, S_NO=10.4152, 
                  S_ALK=4.1256*12)

C1.set_init_solubles(S_I=30, S_S=0.88949, S_O=0.49094, S_NO=10.4152, S_NH=1.7333, 
                      S_ND=0.68828, S_ALK=4.1256*12)
C1.set_init_sludge_solids(X_I=1507.892825, X_S=89.3944512, X_BH=5913.554621, X_BA=372.6863109, X_P=641.7843247, X_ND=2.326138037)
C1.set_init_TSS([12.4969, 18.1132, 29.5402, 68.9781, 356.0747, 
                  356.0747, 356.0747, 356.0747, 356.0747, 6393.9844])
      
############# system simulation ############################
bsm1 = System('BSM1', path=(A1, A2, O1, O2, O3, S1, C1, S2), recycle=(RE, RAS))
bsm1.set_tolerance(rmol=1e-6)

__all__ = (
    'cmps', 'bsm1', 'bsm1_path', 'asm1', 'aer1', 'aer2',
    'Q', 'PE', 'SE', 'WAS', 'RE', 'RAS', 
    *(i.ID for i in bsm1.units),
    )


#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    if method:
        bsm1.simulate(t_span=(0,t), 
                      t_eval=np.arange(0, t+t_step, t_step),
                      method=method, 
                      export_state_to=f'sol_{t}d_{method}.tsv',
                      **kwargs)
    else:
        bsm1.simulate(solver='odeint', 
                      t=np.arange(0, t+t_step, t_step),
                      export_state_to=f'sol_{t}d_odeint.tsv',
                      **kwargs)


if __name__ == '__main__':
    bsm1.reset_cache()
    t = 50
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    method = 'LSODA'
    # method = None
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    run(t, t_step, method=method)
    
    # If want to see a quick plot of the state variable of a certain unit
    # fig, ax = A1.plot_state_over_time(system=bsm1, state_var=('S_S', 'X_I'))