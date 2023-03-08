# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import (
    processes as pc, 
    WasteStream, System, TEA, LCA, PowerUtility,
    ImpactItem as IItm, 
    StreamImpactItem as SIItm,
    )
from qsdsan.utils import ExogenousDynamicVariable as EDV
from exposan.metab_mock import default_inf_concs
from exposan.metab import (
    _impact_item_loaded,
    load_lca_data,
    flex_rhos_adm1
    )
from exposan.metab.units import *

__all__ = ('get_fug_ch4', 'get_NG_eq', 'add_strm_iitm',
           'kWh', 'MJ', 'add_TEA_LCA',
           )

#%%
C0_bulk = np.array([
    1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
    4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 3.212e-02,
    2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 4.708e+00, 1.239e+00,
    4.838e-01, 1.423e+00, 8.978e-01, 2.959e+00, 1.467e+00, 4.924e-02,
    4.000e-02, 2.000e-02, 9.900e+02
    ])

if not _impact_item_loaded: load_lca_data()
bg_offset_CFs = IItm.get_item('biogas_offset').CFs
NaOCl_item = IItm.get_item('NaOCl')
citric_acid_item = IItm.get_item('citric_acid')

def get_fug_ch4(ws):
    '''Returns kg/hr fugitive CH4.'''
    cmps = ws.components
    return ws.imass['S_ch4']*cmps.S_ch4.i_mass

def get_NG_eq(bg):
    '''Returns m3/hr natural gas equivalent.'''
    cmps = bg.components
    KJ_per_kg = cmps.i_mass/cmps.chem_MW*cmps.LHV
    return -sum(bg.mass*KJ_per_kg)*1e-3/39

def add_strm_iitm(sys):
    for ws in sys.products:
        if ws.phase == 'l':
            SIItm(ID=f'{ws.ID}_fugitive_ch4', linked_stream=ws, 
                  flow_getter=get_fug_ch4,
                  GWP100=28, MIR=0.0143794871794872)
        elif ws.phase == 'g':
            SIItm(ID=f'{ws.ID}_NG_offset', linked_stream=ws, functional_unit='m3',
                  flow_getter=get_NG_eq, # m3/hr natural-gas-equivalent
                  **bg_offset_CFs)
    for u in sys.units:
        if isinstance(u, DegassingMembrane):
            SIItm(linked_stream=u.NaOCl, **NaOCl_item.CFs)
            SIItm(linked_stream=u.citric_acid, **citric_acid_item.CFs)

# Operation items
kWh = lambda lca: lca.system.power_utility.rate*lca.lifetime_hr
def MJ(lca):
    sys = lca.system
    duties = [hu.duty for hu in sys.heat_utilities]
    MJ_per_hr = sum(d for d in duties if d > 0)*1e-3
    return MJ_per_hr*lca.lifetime_hr
    
def add_TEA_LCA(sys, irr, lt):
    TEA(sys, discount_rate=irr, lifetime=lt, simulate_system=False, CEPCI=708)   
    LCA(
        sys, lifetime=lt, simulate_system=False,
        electricity = kWh,
        heat_onsite = MJ
        )

#%% Systems
electricity_price=0.0913
def create_system(n_stages=1, reactor='UASB', gas_extraction='passive', 
                  lifetime=30, discount_rate=0.1, flowsheet=None):
    