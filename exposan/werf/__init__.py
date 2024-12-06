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

folder = os.path.dirname(__file__)
data_path = os.path.join(folder, 'data')
results_path = os.path.join(folder, 'results')
figures_path = os.path.join(folder, 'figures')
# To save simulation results and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)

# from qsdsan import (
#     ImpactIndicator as IInd, 
#     ImpactItem as IItm
#     )

# _impact_item_loaded = False
# def _load_lca_data():
#     global _impact_item_loaded
#     if _impact_item_loaded:
#         IInd.clear_registry()
#         IItm.clear_registry()
#     ind_path = os.path.join(data_path, 'TRACI_indicators.xlsx')
#     itm_path = os.path.join(data_path, '_impact_items.xlsx')
#     IInd.load_from_file(ind_path, sheet=0)
#     IItm.load_from_file(itm_path)
#     IItm('Stainless_steel', source='stainless_steel')
#     _impact_item_loaded = True


default_as_init = dict(
    S_O2=1, S_F=0.1, S_A=1.2, S_I=30, S_NH4=2.2, S_N2=10, S_NO3=0.4, 
    S_PO4=5, S_IC=84, X_I=1500, X_S=101, X_H=500, X_PAO=100, X_PP=100,
    X_PHA=10, X_AUT=105, S_K=28, S_Mg=50, S_Na=86, S_Cl=425, S_Ca=140,
    X_CaCO3=1e-5, X_struv=1e-5, X_newb=1e-5, X_ACP=1e-5, X_MgCO3=1e-5, 
    X_AlOH=1e-5, X_AlPO4=1e-5, X_FeOH=1e-5, X_FePO4=1e-5
    )

default_aed_init = dict(
    S_O2=1, S_F=0.4, S_A=0.1, S_I=30, S_NH4=0.2, S_N2=20, S_NO3=30, 
    S_PO4=300, S_IC=84, X_I=8000, X_S=100, X_H=4000, X_PAO=100, X_PP=10,
    X_PHA=1, X_AUT=200, S_K=28, S_Mg=10, S_Na=86, S_Cl=425, S_Ca=10,
    X_CaCO3=1e-2, X_struv=10, X_newb=250, X_ACP=5000, X_MgCO3=1e-2, 
    X_AlOH=1e-2, X_AlPO4=1e-2, X_FeOH=1e-2, X_FePO4=1e-2
    )

default_ad_init = dict(     # in kg/m3
    S_su=0.01, S_aa=0.005, S_fa=0.1, S_va=0.01, S_bu=0.015, S_pro=0.015, S_ac=0.18, 
    S_ch4=0.05, S_IC=0.024, S_IN=3.2, S_IP=0.005, S_I=6.05, 
    X_ch=1.15, X_pr=1.15, X_li=1.73, X_su=0.85, X_aa=0.6, X_fa=0.7, X_c4=0.3, 
    X_pro=0.135, X_ac=0.9, X_h2=0.435, X_I=46.1, 
    S_Ca=0.01, S_Mg=0.05, S_K=0.02, 
    X_CaCO3=1e-8, X_struv=1e-8, X_newb=1e-8, X_ACP=1e-8, X_MgCO3=1e-8, 
    X_AlOH=1e-8, X_AlPO4=1e-8, X_FeOH=1e-8, X_FePO4=1e-8, 
    )
default_ad_init = {k:v*1e3 for k,v in default_ad_init.items()} # convert to mg/L

# default_fctss_init = [18, 28, 45, 90, 305, 305, 305, 305, 305, 5800]
# default_fctss_init = [21, 34, 58, 127, 493, 493, 493, 493, 6243, 11329]
default_fctss_init = [  22.02 ,   36.136,   62.289,  142.932,  610.813,  610.813,
        610.813,  610.813, 3608.754, 9295.076]


from . import G1
from .G1 import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *G1.__all__,
	)


