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
    S_O2=1, S_F=0.4, S_A=0.1, S_I=30, S_NH4=0.2, S_N2=20, S_NO3=1200, 
    S_PO4=550, S_IC=84, X_I=8000, X_S=100, X_H=4000, X_PAO=100, X_PP=10,
    X_PHA=1, X_AUT=200, S_K=28, S_Mg=10, S_Na=86, S_Cl=425, S_Ca=10,
    X_CaCO3=1e-2, X_struv=10, X_newb=250, X_ACP=7000, X_MgCO3=1e-2, 
    X_AlOH=1e-2, X_AlPO4=1e-2, X_FeOH=1e-2, X_FePO4=1e-2, H2O=9.97e5
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
default_fctss_init = [22.02, 36.136, 62.289, 142.932, 610.813, 610.813,
                      610.813, 610.813, 3608.754, 9295.076]

baseline_underflows = dict(
    B1=(59.66, 18.54),
    B2=(57.10, 6.99),
    B3=(60.21, 37.05),
    C1=(136.64, 19.67),
    C2=(131.76, 6.69),
    C3=(138.19, 37.70),
    E2=(107.95, 5.26),
    E2P=(53.92, 19.93),
    F1=(54.64, 19.42),
    G1=(91.23, 27.07),
    G2=(89.71, 15.37),
    G3=(68.98, 39.20),
    H1=(66.13, 23.66),
    I1=(105.83, 20.76),
    I2=(111.29, 9.54),
    I3=(116.21, 31.65),
    N1=(191.54, 26.61),
    N2=(156.34, 13.74),
    )

opt_underflows = dict(
    B1=(57.26, 18.28),
    B2=(55.69, 6.88),
    B3=(58.51, 36.64),
    C1=(139.93, 19.56),
    C2=(135.45, 6.83),
    C3=(143.14, 39.04),
    E2=(109.02, 5.33),
    E2P=(54.68, 19.92),
    F1=(55.27, 19.44),
    G1=(87.40, 26.48),
    G2=(85.35, 14.67),
    G3=(67.57, 38.85),
    H1=(63.90, 23.86),
    I1=(133.73, 24.22),
    I2=(129.08, 12.70),
    I3=(95.53, 25.51),
    N1=(191.04, 26.52),
    N2=(155.39, 13.66),
    )

from . import (
    _units,
    B1, B2, B3, 
    C1, C2, C3, 
    E2, E2P, 
    F1, 
    G1, G2, G3, 
    H1, 
    I1, I2, I3, 
    N1, N2
    )

from ._units import *
from .B1 import *
from .B2 import *
from .B3 import *
from .C1 import *
from .C2 import *
from .C3 import *
from .E2 import *
from .E2P import *
from .F1 import *
from .G1 import *
from .G2 import *
from .G3 import *
from .H1 import *
from .I1 import *
from .I2 import *
from .I3 import *
from .N1 import *
from .N2 import *

def create_system(ID):
    f = globals()[f'create_{ID.lower()}_system']
    return f()

from . import models
from .models import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *_units.__all__,
    *models.__all__,
	)


