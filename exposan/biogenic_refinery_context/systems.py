#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Lewis Rowles <stetsonsc@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Lane To <lane20@illinois.edu>
    
    Aaron Marszewski <aaronpm3@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    sanunits as su,
    ImpactItem,
    System, TEA, LCA,
    )
from qsdsan.utils import clear_lca_registries
from exposan.utils import add_fugitive_items
from exposan.biogenic_refinery_context import (
    _load_components,
    _load_lca_data,
    _units as u,
    const_daily_wage,
    const_person_days,
    discount_rate,
    operator_daily_wage
    )

__all__ = ('create_system',)


# %%

# =============================================================================
# Universal units and functions
# =============================================================================

def batch_create_streams(prefix, phases=('liq', 'sol')):
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    WasteStream('CH4', phase='g', stream_impact_item=item)

    item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item', set_as_source=True)
    WasteStream('N2O', phase='g', stream_impact_item=item)

# %%

# =============================================================================
# Scenario A (sysA): pit latrine with 12,000 users
# =============================================================================

##TODO biosolids - BiogenicRefineryControls - BiogenicRefineryHousing 
#- BiogenicRefineryCarbonizerBase - BiogenicRefineryPollutionControl 
#- BiogenicRefineryOHX - BiogenicRefineryHHX - BiogenicRefineryHHXdryer
# Mixer for N2O - Mixer for CH4 - Trucking of biochar (set distance to 0 for now)


# TODO: Scale CAPex by a factor according to the number of biogenic refineries a POTW would need. Equations for factor below.
    #3. sizing of the system, i.e., multiple systems needed for more than solids loading to 1 BR, \
    # Number of units needed = flow_rate_db/ (550 / 20 * .65) # assuming one unit is capable of treating 550 kg/day (35% moisture based on 20 hr of run time)
    # units: (flow_rate_db kg-db/hr) / (550 kg biosolids/d * 1d/20h * .65 kg dry solids/kg solids@35% MC)

def create_systemA(flowsheet=None):
    # Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    streamA = flowsheet.stream
    batch_create_streams('A')

    ##### Biosolids Inputs #####
    A1 = u.Biosolids('A1', outs=('biosolids'))


    ##### Treatment #####
    A2 = u.BiogenicRefineryControls('A2', ins=A1-0, outs='A2_out')

    A3 = u.BiogenicRefineryHousing('A3', ins=A2-0, outs='A3_out',
                                    const_wage=const_daily_wage,
                                    const_person_days=const_person_days)

    A4 = u.BiogenicRefineryCarbonizerBase('A4', outs=(streamA.biochar, 'A4_hot_gas', 'A4_N2O'))

    A5 = su.BiogenicRefineryPollutionControl('A5', ins=(A4-1, A4-2), outs=('A5_hot_gas_pcd', 'A5_N2O'))

    # Update uptime_ratio in all units to follow carbonizer base
    A5_old_cost = A5._cost
    def update_A5_uptime_ratio():
        #
        A8.uptime_ratio = A7.uptime_ratio = A6.uptime_ratio = A5.uptime_ratio = A4.uptime_ratio
        A5_old_cost()
    A5._cost = update_A5_uptime_ratio

    A6 = su.BiogenicRefineryOHX('A6', ins=A5-0, outs='A6_hot_gas')
    A7 = su.BiogenicRefineryHHX('A7', ins=A6-0, outs='A7_hot_gas')
    A8 = su.BiogenicRefineryHHXdryer('A8', ins=(A1-0, A7-0), outs=('waste_out', 'A8_N2O', 'A8_CH4')) 
    A8-0-A4

    A9 = su.Mixer('A9', ins=(A8-2), outs=streamA.CH4)
    A9.add_specification(lambda: add_fugitive_items(A9, 'CH4_item'))
    A9.line = 'fugitive CH4 mixer'

    A10 = su.Mixer('A10', ins=(A4-2, A5-1, A8-1), outs=streamA.N2O)
    A10.add_specification(lambda: add_fugitive_items(A10, 'N2O_item'))
    A10.line = 'fugitive N2O mixer'
    
    ##### Conveyance of Biochar #####
    A11 = su.Trucking('A11', ins=A4-0, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=A2.emptying_period, interval_unit='yr',
                      loss_ratio=0.02)

    ##### Simulation, TEA, and LCA ##### TODO: scale annual_labor to scale_factor?
    sysA = System('sysA', path=(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11))
    teaA = TEA(system=sysA, discount_rate=discount_rate,
               start_year=2020, lifetime=20, uptime_ratio=1,
               lang_factor=None, annual_maintenance=0,
               annual_labor=(operator_daily_wage*3*365))

    # 12 is assuming the device is running 12 hr per day (50% of the time)
    # this isn't adjusted through `uptime_ratio` because other OPEX calculation
    # in this unit needs `uptime_ratio` to be 1
    get_powerA = lambda: sum([(u.power_utility.rate*u.uptime_ratio)
                              for u in sysA.units])*(365*teaA.lifetime)*12
    LCA(system=sysA, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)

    return sysA


# %%

# =============================================================================
# Wrapper function
# =============================================================================

def create_system(system_ID='A', flowsheet=None):
    ID = system_ID.lower().lstrip('sys').upper() # so that it'll work for "sysA"/"A"
    reload_lca = False

    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f'br{ID}'
        if hasattr(main_flowsheet.flowsheet, flowsheet_ID): # clear flowsheet
            getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
            clear_lca_registries()
            reload_lca = True
        flowsheet = Flowsheet(flowsheet_ID)
        main_flowsheet.set_flowsheet(flowsheet)

    _load_components()
    _load_lca_data(reload_lca)

    if system_ID == 'A': f = create_systemA
    else: raise ValueError(f'`system_ID` can only be "A", not "{ID}".')

    try: system = f(flowsheet)
    except:
        _load_components(reload=True)
        system = f(flowsheet)

    return system