# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import os, sys
sys.path.insert(0, "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSDsan")
sys.path.insert(0, "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/biosteam")

from biosteam import System
import numpy as np
from qsdsan import sanunits as su
from qsdsan import processes as pc
from qsdsan import WasteStream, set_thermo
from qsdsan.utils import time_printer

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
                path='_asm1.tsv')

# asm1 = pc.ASM1(Y_A=0.24, Y_H=0.666, f_P=0.08, i_XB=0.086, i_XP=0.06,
#                 mu_H=6.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.62,
#                 eta_g=0.8, eta_h=0.4, k_h=3.0, K_X=0.03, mu_A=0.8,
#                 K_NH=1.0, b_A=0.04, K_O_A=0.4, k_a=0.08, fr_SS_COD=1/1.48, 
#                 path='_asm1.tsv')

############# create unit operations #####################
M1 = su.Mixer('M1', ins=[PE, RE, RAS])
# HD = su.HydraulicDelay('HD', ins=M1-0, t_delay=1e-3)
# A1 = su.CSTR('A1', ins=HD-0, V_max=V_an,
#               aeration=None, suspended_growth_model=asm1)
A1 = su.CSTR('A1', ins=M1-0, V_max=V_an,
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

# HD.set_init_conc(S_I=30.0, S_S=14.61, X_I=1149.1252, X_S=89.33, X_BH=2542.171, X_BA=148.46,
#                   X_P=448.18, S_O=0.393, S_NO=8.33, S_NH=7.70, S_ND=1.94, X_ND=5.61, 
#                   S_ALK=4.70*12)
# HD.set_init_flow(92230.0)

A1.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
A2.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
O1.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
O2.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
                  S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)
# O3.set_init_conc(S_I=30, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, X_P=100,
#                   S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20, S_ALK=7*12)

# A1.set_init_conc(S_I=30.0, S_S=2.8082, X_I=1149.1252, X_S=82.1349, X_BH=2551.7658, X_BA=148.3894, 
#                   X_P=448.8519, S_O=0.0042984, S_NH=7.9179, S_ND=1.21661, X_ND=5.2849, S_NO=5.3699, 
#                   S_ALK=4.9277*12)
# A2.set_init_conc(S_I=30.0, S_S=1.4588, X_I=1149.1252, X_S=76.3862, X_BH=2553.3851, X_BA=148.3091, 
#                   X_P=449.5227, S_O=6.3132e-05, S_NH=8.3444, S_ND=0.88206, X_ND=5.0291, S_NO=3.662, 
#                   S_ALK=5.0802*12)
# O1.set_init_conc(S_I=30.0, S_S=1.1495, X_I=1149.1252, X_S=64.8549, X_BH=2557.1314, X_BA=148.9413, 
#                   X_P=450.4183, S_O=1.7184, S_NH=5.5479, S_ND=0.82889, X_ND=4.3924, S_NO=6.5409, 
#                   S_ALK=4.6748*12)
# O2.set_init_conc(S_I=30, S_S=0.99532, X_I=1149.1252, X_S=55.694, X_BH=2559.1826, X_BA=149.5271, 
#                   X_P=451.3147, S_O=2.4289, S_NO=9.299, S_NH=2.9674, S_ND=0.76679, X_ND=3.879, 
#                   S_ALK=4.2935*12)
O3.set_init_conc(S_I=30.0, S_S=0.88949, X_I=1149.1252, X_S=49.3056, X_BH=2559.3436, X_BA=149.7971, 
                  X_P=452.2111, S_O=0.49094, S_NH=1.7333, S_ND=0.68828, X_ND=3.5272, S_NO=10.4152, 
                  S_ALK=4.1256*12)

C1.set_init_solubles(S_I=30, S_S=0.88949, S_O=0.49094, S_NO=10.4152, S_NH=1.7333, 
                      S_ND=0.68828, S_ALK=4.1256*12)
# C1.set_init_sludge_solids(X_I=1149.1252, X_S=89.33, X_BH=2542.17, X_BA=148.46, X_P=448.18, X_ND=5.61)
# C1.set_init_TSS([10, 20, 40, 70, 200, 300, 350, 350, 2000, 4000])
C1.set_init_sludge_solids(X_I=1507.892825, X_S=89.3944512, X_BH=5913.554621, X_BA=372.6863109, X_P=641.7843247, X_ND=2.326138037)
C1.set_init_TSS([12.4969, 18.1132, 29.5402, 68.9781, 356.0747, 
                  356.0747, 356.0747, 356.0747, 356.0747, 6393.9844])
      
############# system simulation ############################

# bio = System('Biological_treatment', path=(A1, A2, O1, O2, O3, S1), recycle=(RE,))
# bsm1 = System('BSM1', path=(bio, C1, S2), recycle=(RAS,))
# bio.set_tolerance(rmol=1e-6)
# bsm1 = System('BSM1', path=(M1, HD, A1, A2, O1, O2, O3, S1, C1, S2), recycle=(RE, RAS))
bsm1 = System('BSM1', path=(M1, A1, A2, O1, O2, O3, S1, C1, S2), recycle=(RE, RAS))
# bsm1.set_tolerance(rmol=1e-6)

#%%
@time_printer
def run(T, t_step, method=None, **kwargs):
    if method:
        bsm1.simulate(t_span=(0,T), 
                      t_eval=np.arange(0, T+t_step, t_step),
                      method=method, 
                      **kwargs)
    else:
        bsm1.simulate(solver='odeint', 
                      t=np.arange(0, T+t_step, t_step), 
                      **kwargs)


if __name__ == '__main__':
    T = 55
    t_step = 0.5
    method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    # method = None
    print(f'\nMethod {method}\n------------')
    print(f'Time span 0-{T}d \n')
    run(T, t_step, method=method)
    os.rename('sol.txt', f'sol_{T}_{method}.txt')
    bsm1.reset_cache()