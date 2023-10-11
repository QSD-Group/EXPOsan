#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Bright Elijah <be05055@georgiasouthern.edu & brightcarlelijah@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

#!!! Consider adding option to change water source and ppl (get_number_of_households)

from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    ImpactItem,
    System, TEA, LCA,
    )
from qsdsan.utils import clear_lca_registries
from exposan.pou_disinfection import (
    _load_components,
    _load_lca_data,
    _units as u,
    discount_rate,
    household_size,
    get_number_of_households,
    lifetime,
    start_year,
    )

__all__ = ('create_system',)


# %%

# =============================================================================
# System A: POU Chlorination
# =============================================================================

#!!! Seems like only sysA needs NaClO and PE
def create_systemA(flowsheet=None):
    item = ImpactItem.get_item('NaClO_item')
    A_naclo = WasteStream('A_naclo', phase='g', stream_impact_item=item)

    item = ImpactItem.get_item('Polyethylene_item')
    A_polyethylene = WasteStream('A_polyethylene', phase='g', stream_impact_item=item)

    
    A1 = u.RawWater('A1', outs=('raw_water'), household_size=household_size, 
                    number_of_households=get_number_of_households())
    
    A2 = u.POUChlorination('A2', ins=(A1-0, A_naclo, A_polyethylene), 
                            outs='treated_water', 
                            number_of_households=get_number_of_households())
    
    sysA = System('sysA', path=(A1, A2))    
    
    TEA(system=sysA, discount_rate=discount_rate, start_year=start_year,
               lifetime=lifetime, uptime_ratio=1, lang_factor=None,
               annual_maintenance=0, annual_labor=0)
    
    LCA(system=sysA, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True)
    
    return sysA


# =============================================================================
# System B: AgNP CWF
# =============================================================================

def create_systemB(flowsheet=None):
    B1 = u.RawWater('B1', outs=('raw_water'), household_size=household_size, 
                    number_of_households=get_number_of_households())
    
    B2 = u.AgNP_CWF('B2', ins=B1-0, outs='treated_water', 
                    number_of_households=get_number_of_households())
    
    ############### Simulation, TEA, and LCA ###############
    sysB = System('sysV', path=(B1, B2))
    
    TEA(system=sysB, discount_rate=discount_rate, start_year=start_year,
        lifetime=lifetime, uptime_ratio=1, lang_factor=None,
        annual_maintenance=0, annual_labor=0)
    
    LCA(system=sysB, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True)
    
    return sysB


# =============================================================================
# System C: POU UV
# =============================================================================

def create_systemC(flowsheet=None):
    C1 = u.RawWater('C1', outs=('raw_water'), household_size=household_size, 
                      number_of_households=get_number_of_households())
    
    C2 = u.POU_UV('C2', ins=(C1-0), outs='treated_water', 
                  number_of_households=get_number_of_households())
    
    sysC = System('sysC', path=(C1, C2))
    
    TEA(system=sysC, discount_rate=discount_rate, start_year=start_year,
        lifetime=lifetime, uptime_ratio=1, lang_factor=None,
        annual_maintenance=0, annual_labor=0)
    
    LCA(system=sysC, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True, E_item=lambda: C2.power_utility.rate*(365*12*5))
    
    return sysC

# =============================================================================
#  system D: UV LED
# =============================================================================

def create_systemD(flowsheet=None):
    D1 = u.RawWater('D1', outs=('raw_water'), household_size=household_size, 
                      number_of_households=get_number_of_households())
    
    D2 = u.UV_LED('D2', ins=(D1-0), outs='treated_water', 
                  number_of_households=get_number_of_households())
    
    sysD = System('sysD', path=(D1, D2))
    
    TEA(system=sysD, discount_rate=discount_rate, start_year=start_year,
        lifetime=lifetime, uptime_ratio=1, lang_factor=None,
        annual_maintenance=0, annual_labor=0)
    
    LCA(system=sysD, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True, E_item=lambda: D2.power_utility.rate*(365*24*5))
    
    return sysD


# %%

# =============================================================================
# Wrapper function
# =============================================================================

def create_system(system_ID='A', flowsheet=None):
    ID = system_ID.lower().lstrip('sys').upper()  # so that it'll work for "sysA"/"A"
    reload_lca = False

    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f'pou{ID}'
        if hasattr(main_flowsheet.flowsheet, flowsheet_ID): # clear flowsheet
            getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
            clear_lca_registries()
            reload_lca = True
        flowsheet = Flowsheet(flowsheet_ID)
        main_flowsheet.set_flowsheet(flowsheet)

    _load_components()
    _load_lca_data(reload_lca)

    if system_ID == 'A': f = create_systemA
    elif system_ID == 'B': f = create_systemB
    elif system_ID == 'C': f = create_systemC
    elif system_ID == 'D': f = create_systemC
    else: raise ValueError(f'`system_ID` can only be "A" , "B", "C", or "D", not "{ID}".')

    try: system = f(flowsheet)
    except:
        _load_components(reload=True)
        system = f(flowsheet)

    return system