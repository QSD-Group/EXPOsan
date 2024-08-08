#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>
    
References:
    
(1) Lang, C.; Lee, B. Heat Transfer Fluid Life Time Analysis of Diphenyl
    Oxide/Biphenyl Grades for Concentrated Solar Power Plants. Energy Procedia
    2015, 69, 672–680. https://doi.org/10.1016/j.egypro.2015.03.077.
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs

__all__ = ('_load_process_settings',)

def _load_process_settings():
# =============================================================================
#     add a heating agent
# =============================================================================
    # use DOWTHERM(TM) A Heat Transfer Fluid (HTF) as the heating agent
    # DOWTHERM(TM) A HTF = 73.5% diphenyl oxide (DPO) + 26.5% Biphenyl (BIP)
    # critical temperature for HTF: 497 C
    # critical pressure for HTF: 313.4 kPa
    # https://www.dow.com/en-us/pdp.dowtherm-a-heat-transfer-fluid.238000z.\
    # html#tech-content (accessed 2022-11-16)
    
    DPO_chem = qs.Chemical('DPO_chem', search_ID='101-84-8')
    BIP_chem = qs.Chemical('BIP_chem', search_ID='92-52-4')
    
    DPO = qs.Component.from_chemical('DPO', chemical=DPO_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    
    BIP = qs.Component.from_chemical('BIP', chemical=BIP_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    
    HTF_thermo = bst.Thermo((DPO, BIP,))
    
    HTF = bst.UtilityAgent('HTF', DPO=0.735, BIP=0.265, T=673.15, P=951477, phase='g',
                           # 400 C (673.15 K) and 138 psig (951477 pa) are max temp and pressure for HTF
                           thermo=HTF_thermo,
                           # T_limit = 495 F (530.372 K) is the highest temp that vapor can exist
                           regeneration_price=1) # Lang
                           # use default heat transfer efficiency (1)
    # Temperature and pressure: https://www.dow.com/content/dam/dcc/documents/\
    # en-us/app-tech-guide/176/176-01334-01-dowtherm-heat-transfer-fluids-\
    # engineering-manual.pdf?iframe=true (accessed 2022-11-16)
    bst.HeatUtility.heating_agents.append(HTF)

    bst.CE = qs.CEPCI_by_year[2020] # use 2020$ to match up with latest PNNL report
    
# =============================================================================
#     set utility prices
# =============================================================================
    bst.PowerUtility.price = 0.06879 # !!! the electricity price can be adjusted here
    
    # # These utilities are provided by CHP thus cost already considered
    # # setting the regeneration price to 0 or not will not affect the final results
    # # as the utility cost will be positive for the unit that consumes it
    # # but negative for HXN/CHP as they produce it
    # for adj in ('low', 'medium', 'high'):
    #     steam = bst.HeatUtility.get_agent(f'{adj}_pressure_steam')
    #     steam.heat_transfer_price = steam.regeneration_price = 0.
