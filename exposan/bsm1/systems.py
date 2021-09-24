# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import sys
sys.path.insert(0, "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSDsan")
sys.path.insert(0, "C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/biosteam")


# import thermosteam as tmo
# import biosteam as bst
# bst.RAISE_GRAPHVIZ_EXCEPTION = True
from biosteam import System
# import numpy as np
from qsdsan import sanunits as su
from qsdsan import processes as pc
from qsdsan import WasteStream, set_thermo

# =============================================================================
# Benchmark Simulation Model No. 1
# =============================================================================

############# load components and set thermo #############
cmps = pc.load_asm1_cmps()
set_thermo(cmps)

############# create WasteStream objects #################
Q = 18446           # influent flowrate [m3/d]
Temp = 275.15+15    # temperature [K]

PE = WasteStream('primary_effluent', T=Temp)
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

SE = WasteStream('secondary_effluent', T=Temp)
WAS = WasteStream('waste_activated_sludge', T=Temp)
RE = WasteStream('recycled_effluent', T=Temp)
RAS = WasteStream('recycled_activated_sludge', T=Temp)

############# load and tailor process models #############
V_an = 1000    # anoxic zone tank volume
V_ae = 1333    # aerated zone tank volume
Q_was = 385    # sludge wastage flowrate
Q_ras = 18446    # recycle sludge flowrate

aer = pc.DiffusedAeration('Fixed_Aeration', 'S_O', KLa=240, DOsat=8.0)
asm1 = pc.ASM1(components=cmps,
                Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
                mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
                eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
                K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05)


############# create unit operations #####################
A1 = su.CSTR('Anoxic_1', ins=[PE, RE, RAS], V_max=V_an,
             aeration=None, suspended_growth_model=asm1)
A1.set_init_conc(S_I=30.0, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
                 S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20.0, S_ALK=7*12)

A2 = su.CSTR('Anoxic_2', A1-0, V_max=V_an,
             aeration=None, suspended_growth_model=asm1)
A2.set_init_conc(S_I=30.0, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
                 S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20.0, S_ALK=7*12)

O1 = su.CSTR('Aerobic_1', A2-0, V_max=V_ae, aeration=aer,
             DO_ID='S_O', suspended_growth_model=asm1)
O1.set_init_conc(S_I=30.0, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
                 S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20.0, S_ALK=7*12)

O2 = su.CSTR('Aerobic_2', O1-0, V_max=V_ae, aeration=aer, 
             DO_ID='S_O', suspended_growth_model=asm1)                 
O2.set_init_conc(S_I=30.0, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
                 S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20.0, S_ALK=7*12)

O3 = su.CSTR('Aerobic_3', O2-0, V_max=V_ae, aeration=2.0,
             DO_ID='S_O', suspended_growth_model=asm1)
O3.set_init_conc(S_I=30.0, S_S=5.0, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
                 S_O=2.0, S_NH=2.0, S_ND=1.0, X_ND=1.0, S_NO=20.0, S_ALK=7*12)

S1 = su.Splitter('S1', O3-0, [RE, 'treated'], split=0.6, init_with='WasteStream')

C1 = su.FlatBottomCircularClarifier('C1', S1-1, [SE, 'sludge'], 
                                    sludge_flow_rate=Q_ras+Q_was, surface_area=1500,  
                                    height=4, N_layer=10, feed_layer=4, 
                                    # height=12, N_layer=3, feed_layer=2,
                                    X_threshold=3000, v_max=474, v_max_practical=250, 
                                    rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)
C1.set_init_solubles(S_I=30, S_S=5.0, S_O=2.0, S_NO=20, S_NH=2.0, S_ALK=7*12)
C1.set_init_TSS([10, 20, 40, 70, 200, 300, 350, 350, 2000, 4000])

S2 = su.Splitter('S2', C1-1, [RAS, WAS], split=Q_ras/(Q_ras+Q_was), init_with='WasteStream')

############# system simulation ############################

bio = System('Biological_treatment', path=(A1, A2, O1, O2, O3))
bsm1 = System('BSM1', path=(bio, S1, C1, S2), recycle=(RE, RAS))

if __name__ == '__main__':
    bsm1.simulate(t_span = (0, 0.05))
