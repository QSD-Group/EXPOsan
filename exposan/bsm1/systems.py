# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
# import sys
# sys.path.insert(0, "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSDsan")
# sys.path.insert(0, "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/biosteam")


import thermosteam as tmo
from biosteam import System
# import numpy as np
from qsdsan import sanunits as su
from qsdsan import processes as pc
from qsdsan import WasteStream

# =============================================================================
# Benchmark Simulation Model No. 1
# =============================================================================

############# load components and set thermo #############
cmps = pc.load_asm1_cmps()
tmo.settings.set_thermo(cmps)

############# create WasteStream objects #################
Q = 18446           # influent flowrate [m3/d]
Temp = 275.15+15    # temperature [K]
PE = WasteStream('primary_effluent', T=Temp, units='kg/d',
                 H2O=Q*1e3,
                 S_S=69.5*Q*1e-3,
                 X_BH=28.17*Q*1e-3,
                 X_S=202.32*Q*1e-3,
                 X_I=51.2*Q*1e-3,
                 S_NH=31.56*Q*1e-3,
                 S_I=30*Q*1e-3,
                 S_ND=6.95*Q*1e-3,
                 X_ND=10.59*Q*1e-3,
                 S_ALK=7*12*Q*1e-3)

SE = WasteStream('secondary_effluent', T=Temp)
WAS = WasteStream('waste_activated_sludge', T=Temp)
RE = WasteStream('recycled_effluent', T=Temp)
RAS = WasteStream('recycled_activated_sludge', T=Temp)

############# load and tailor process models #############
V_an = 1000    # anoxic zone tank volume
V_ae = 1333    # aerated zone tank volume
Q_was = 385    # sludge wastage flowrate
Q_ras = 18446    # total sludge flowrate

aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa=240, DOsat=8.0)
asm1 = pc.ASM1(components=cmps,
                Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
                mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
                eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
                K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05)

# asm1 = 1
############# create unit operations #####################
A1 = su.CSTR('Anoxic_1', ins=[PE, RE, RAS], V_max=V_an,
             aeration=None, suspended_growth_model=asm1)
# A1.init_state = np.array([30, 2.81, 1149, 82.1, 2552, 148, 449,
#                           4.3e-3, 5.37, 7.92, 1.22, 5.28, 4.93*12, 0, 1e6])

A2 = su.CSTR('Anoxic_2', A1-0, V_max=V_an,
             aeration=None, suspended_growth_model=asm1)
# A2.init_state = np.array([30, 1.46, 1149, 76.4, 2553, 148, 450,
#                           6.31e-5, 3.66, 8.34, 0.882, 5.03, 5.08*12, 0, 1e6])

O1 = su.CSTR('Aerobic_1', A2-0, V_max=V_ae, aeration=aer,
             DO_ID='S_O', suspended_growth_model=asm1)
# O1.init_state = np.array([30, 1.15, 1149, 64.9, 2557, 149, 450,
#                           1.72, 6.54, 5.55, 0.829, 4.39, 4.67*12, 0, 1e6])

O2 = su.CSTR('Aerobic_2', O1-0, V_max=V_ae, aeration=aer,
             DO_ID='S_O', suspended_growth_model=asm1)
# O2.init_state = np.array([30, 0.995, 1149, 55.7, 2559, 150, 451,
#                           2.43, 9.3, 2.97, 0.767, 3.88, 4.29*12, 0, 1e6])

O3 = su.CSTR('Aerobic_3', O2-0, V_max=V_ae, aeration=2.0,
             DO_ID='S_O', suspended_growth_model=asm1)
# O3.init_state = np.array([30, 0.889, 1149, 49.3, 2559, 150, 452,
#                           0.491, 10.4, 1.73, 0.688, 3.53, 4.13*12, 0, 1e6])

S1 = su.Splitter('S1', O3-0, [RE, 'treated'], split=0.6, init_with='WasteStream')

C1 = su.FlatBottomCircularClarifier('C1', S1-1, [SE, 'sludge'],
                                    sludge_flow_rate=Q_ras+Q_was, surface_area=1500,
                                    height=4, N_layer=10, feed_layer=4,
                                    X_threshold=3000, v_max=474, v_max_practical=250,
                                    rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)

S2 = su.Splitter('S2', C1-1, [RAS, WAS], split=Q_ras/(Q_ras+Q_was), init_with='WasteStream')

############# system simulation ############################

bio = System('Biological_treatment', path=(A1, A2, O1, O2, O3))
bsm1 = System('BSM1', path=(bio, S1, C1, S2), recycle=(RE, RAS))
# bsm1 = System('BSM1', path=(A1, A2, O1, O2, O3, S1, C1, S2), recycle=(RE, RAS))

# O3 = su.CSTR('Aerobic_3', A1-0, V_max=V_ae, aeration=2.0,
#              DO_ID='S_O', suspended_growth_model=asm1)
# sys = System('sys', path=(A1, O3, S1, C1, S2), recycle=(RE, RAS))
if __name__ == '__main__':
    bsm1.simulate(t_span = (0, 0.1))