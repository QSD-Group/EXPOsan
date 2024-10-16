# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 07:52:15 2023

@author: joy_c
"""
C_mw = 12.0107
N_mw = 14.0067

ww1 = {
    'S_su':0.23,
    'S_aa':0.23,
    'S_fa':0.23,
    'S_va':0.23,
    'S_bu':0.23,
    'S_pro':0.24,
    'S_ac':0.24,
    # 'S_IC':0.04*C_mw,
    # 'S_IN':0.01*N_mw,
    'S_IC':0.001,
    'S_I':0.02,
    'X_c':0.06,
    'X_ch':0.045,
    'X_pr':0.3,
    'X_li':0.15,
    'X_I':0.015, 
    'S_cat':0.04, 
    'S_an':0.02
    }

ww2 = {
    'S_su':0.41,
    'S_aa':0.36,
    'S_fa':0.41,
    'S_va':0.41,
    'S_bu':0.41,
    'S_pro':0.48,
    'S_ac':0.48,
    # 'S_IC':0.04*C_mw,
    # 'S_IN':0.01*N_mw,
    'S_IC':0.001,
    'S_I':0.024,
    'X_c':0.12,
    'X_ch':0.36,
    'X_pr':0.6,
    'X_li':0.3,
    'X_I':0.03, 
    'S_cat':0.04, 
    'S_an':0.02
    }

ww3 = {
    'S_su':2.4,
    'S_aa':0.36,
    'S_fa':0.48,
    'S_va':0.48,
    'S_bu':0.48,
    'S_pro':0.48,
    'S_ac':0.48,
    # 'S_IC':0.04*C_mw,
    # 'S_IN':0.01*N_mw,
    'S_IC':0.001,
    'S_I':0.024,
    'X_c':0.12,
    'X_ch':0.36,
    'X_pr':0.6,
    'X_li':0.3,
    'X_I':0.03, 
    'S_cat':0.04, 
    'S_an':0.02
    }

#%%

from exposan.metab import create_system

sys = create_system(tot_HRT=1, reactor_type='PB', inf_concs=ww1, Q=5)
sys.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')

