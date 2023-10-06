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

import biosteam as bst, qsdsan as qs
from qsdsan import WasteStream, ImpactIndicator, ImpactItem, StreamImpactItem, SimpleTEA, LCA
from . import _units as u
from exposan.POU_dis import results_path

from exposan.bwaise._cmps import cmps



# =============================================================================
# Unit parameters
# =============================================================================

currency = qs.currency = 'USD'
qs.CEPCI = qs.CEPCI_by_year[2018]

discount_rate = 0.05


qs.set_thermo(cmps)


household_size = 6
household_per_container = 1
get_pou_user = lambda: household_size * household_per_container

ppl_1k = 1000
ppl_500 = 500
get_ppl = lambda kind: ppl_1k if kind=='1k' else ppl_500


# =============================================================================
# Prices and GWP CFs
# =============================================================================



price_dct = {
    'Electricity': 0.17,
    'Concrete': 194,
    'Steel': 2.665,
    'NaClO': 1.96, #updated later 
    'Polyethylene': 0, #update later
    }

GWP_dct = {
    'Electricity': 0.1135,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'NaClO': 1.0, #update later
    'Polyethylene': 1.0, #update later
    }



GWP = ImpactIndicator('GWP', unit='kg CO2')

bst.PowerUtility.price = price_dct['Electricity']

#ImpactItem.get_item('NaClO').price = price_dct['NaClO']

# =============================================================================
# Universal units and functions
# =============================================================================

def batch_create_stream_items(kind):
    if kind == 'original':
        for k, v in GWP_dct.items():
            if k == 'Electricity':
                ImpactItem(ID='E_item', functional_unit='kWh', GWP=v)
            else:
                StreamImpactItem(ID=f'{k}_item', GWP=v)
    elif kind == 'new':
        E_factor = {'GW2ECO': 0.000000000532, # global warming to ecosystem
                    'GW2HH': 0.0000000812, # global warming to human health
                    'OD2HH': 0.00134} # stratospheric ozone depletion to human health
        H_factor = {'GW2ECO': 0.0000000028, 'GW2HH': 0.000000928, 'OD2HH': 0.000531}
        I_factor = {'GW2ECO': 0.000000025, 'GW2HH': 0.0000125, 'OD2HH': 0.000237}

        StreamImpactItem(ID='CH4_item',
                         E_EcosystemQuality_Total=E_factor['GW2ECO']*4.8,
                         E_HumanHealth_Total=E_factor['GW2HH']*4.8,
                         H_EcosystemQuality_Total=H_factor['GW2ECO']*34,
                         H_HumanHealth_Total=H_factor['GW2HH']*34,
                         I_EcosystemQuality_Total=I_factor['GW2ECO']*84,
                         I_HumanHealth_Total=I_factor['GW2HH']*84
                         )
        StreamImpactItem(ID='N2O_item',
                         E_EcosystemQuality_Total=E_factor['GW2ECO']*78.8,
                         # From climate change + ozone depletion
                         E_HumanHealth_Total=\
                             E_factor['GW2HH']*78.8+E_factor['OD2HH']*0.017,
                         H_EcosystemQuality_Total=H_factor['GW2ECO']*298,
                         H_HumanHealth_Total=\
                             H_factor['GW2HH']*298+H_factor['OD2HH']*0.011,
                         I_EcosystemQuality_Total=I_factor['GW2ECO']*264,
                         I_HumanHealth_Total=\
                             I_factor['GW2HH']*264+I_factor['OD2HH']*0.007
                         )
    else:
        raise ValueError(f'`kind` can only be "original" or "new", not "{kind}".')

    global _ImpactItem_LOADED
    _ImpactItem_LOADED = True

NaClO_item = StreamImpactItem(ID='NaClO_item', GWP=GWP_dct['NaClO'])
Polyethylene_item = StreamImpactItem(ID='Polyethylene_item', GWP=GWP_dct['Polyethylene'])

def batch_create_streams(prefix):


    stream_dct = {}
    # item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    # stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g', stream_impact_item=item)
    # # CH4.stream_impact_item = ImpactItem.get_item('CH4_item').copy(stream=CH4, set_as_source=True)

    # item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item', set_as_source=True)
    # stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g', stream_impact_item=item)
    # # N2O.stream_impact_item = ImpactItem.get_item('N2O_item').copy(stream=N2O, set_as_source=True)

    item = ImpactItem.get_item('NaClO_item').copy(f'{prefix}_NaClO_item', set_as_source=True)
    stream_dct['NaClO'] = WasteStream(f'{prefix}_NaClO', phase='l', price=price_dct['NaClO'], 
                                      stream_impact_item=item)

    item = ImpactItem.get_item('Polyethylene_item').copy(f'{prefix}_Polyethylene_item', set_as_source=True)
    stream_dct['Polyethylene'] = WasteStream(f'{prefix}_Polyethylene', phase='s', price=price_dct['Polyethylene'], 
                                      stream_impact_item=item)
    
    return stream_dct



# %%

# =============================================================================
# POU Chlorination (sysD): test system
# =============================================================================



flowsheetD = bst.Flowsheet('sysD')
bst.main_flowsheet.set_flowsheet(flowsheetD)
streamsD = batch_create_streams('D')

#################### Human Inputs ####################
D1 = u.RawWater('D1', outs=('raw_water'), household_size=household_size, 
                 number_of_households=(get_ppl('1k')/household_size))

D2 = u.POUChlorination('D2', ins=(D1-0, streamsD['NaClO'], streamsD['Polyethylene']), 
                        outs='treated_water', 
                        number_of_households=(get_ppl('1k')/household_size))

############### Simulation, TEA, and LCA ###############
sysD = bst.System('sysD', path=(D1, D2))


sysD.simulate()


teaD = SimpleTEA(system=sysD, discount_rate=discount_rate, start_year=2018,
                 lifetime=5, uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=0) #12*3e6*12

lcaD = LCA(system=sysD, lifetime=5, lifetime_unit='yr', uptime_ratio=1,
           annualize_construction=True,
           E_item=0)


# =============================================================================
# AgNP CWF (sysE): test system
# =============================================================================



flowsheetE = bst.Flowsheet('sysE')
bst.main_flowsheet.set_flowsheet(flowsheetE)
streamsE = batch_create_streams('E')

#################### Human Inputs ####################
E1 = u.RawWater('E1', outs=('raw_water'), household_size=household_size, 
                  number_of_households=(get_ppl('1k')/household_size))

E2 = u.AgNP_CWF('E2', ins=(E1-0), outs='treated_water', 
                        number_of_households=(get_ppl('1k')/household_size))

############### Simulation, TEA, and LCA ###############
sysE = bst.System('sysE', path=(E1, E2))


sysE.simulate()


teaE = SimpleTEA(system=sysE, discount_rate=discount_rate, start_year=2018,
                  lifetime=5, uptime_ratio=1, lang_factor=None,
                  annual_maintenance=0, annual_labor=0) #12*3e6*12

lcaE = LCA(system=sysE, lifetime=5, lifetime_unit='yr', uptime_ratio=1,
            annualize_construction=True,
            # Assuming all additional WWTP OPEX from electricity
            E_item=0)

# =============================================================================
# POU UV (sysF): test system
# =============================================================================



flowsheetF = bst.Flowsheet('sysF')
bst.main_flowsheet.set_flowsheet(flowsheetF)
streamsF = batch_create_streams('F')

#################### Human Inputs ####################
F1 = u.RawWater('F1', outs=('raw_water'), household_size=household_size, 
                  number_of_households=(get_ppl('1k')/household_size))

F2 = u.POU_UV('F2', ins=(F1-0), outs='treated_water', 
                        number_of_households=(get_ppl('1k')/household_size))

############### Simulation, TEA, and LCA ###############
sysF = bst.System('sysF', path=(F1, F2))


sysF.simulate()
 

teaF = SimpleTEA(system=sysF, discount_rate=discount_rate, start_year=2018,
                  lifetime=5, uptime_ratio=1, lang_factor=None,
                  annual_maintenance=0, annual_labor=0) #12*3e6*12

lcaF = LCA(system=sysF, lifetime=5, lifetime_unit='yr', uptime_ratio=1,
            annualize_construction=True,
            # Assuming all additional WWTP OPEX from electricity
            E_item=0)

# =============================================================================
#  UV LED (sysG): test system
# =============================================================================



flowsheetG = bst.Flowsheet('sysG')
bst.main_flowsheet.set_flowsheet(flowsheetG)
streamsG = batch_create_streams('G')

#################### Human Inputs ####################
G1 = u.RawWater('G1', outs=('raw_water'), household_size=household_size, 
                  number_of_households=(get_ppl('1k')/household_size))

G2 = u.UV_LED('G2', ins=(G1-0), outs='treated_water', 
                        number_of_households=(get_ppl('1k')/household_size))

############### Simulation, TEA, and LCA ###############
sysG = bst.System('sysG', path=(G1, G2))


sysG.simulate()
 

teaG = SimpleTEA(system=sysG, discount_rate=discount_rate, start_year=2018,
                  lifetime=5, uptime_ratio=1, lang_factor=None,
                  annual_maintenance=0, annual_labor=0) #12*3e6*12

lcaG = LCA(system=sysG, lifetime=5, lifetime_unit='yr', uptime_ratio=1,
            annualize_construction=True,
            # Assuming all additional WWTP OPEX from electricity
            E_item=0)

# %%
 
# =============================================================================
# Util functions

def get_total_inputs(unit, multiplier=1):
    if len(unit.ins) == 0: # 
        ins = unit.outs
    else:
        ins = unit.ins
    inputs = {}
    inputs['NaClO'] = sum(i.imass['NaClO'] for i in ins)
    inputs['Polyethylene'] = sum(i.imass['Polyethylene'] for i in ins)
    inputs['Ecoli'] = sum(i.Ecoli*i.F_vol/1e3 for i in ins)
    for i, j in inputs.items():
        inputs[i] = j * multiplier
    return inputs


#breakpoint()
sys_dct = {
    'ppl': dict(sysD=get_ppl('1k'), sysE=get_ppl('1k')),
    'input_unit': dict(sysD=D1,sysE=E1, sysF=F1, sysG=G1),
    'liq_unit': dict(sysD=D1, sysE=E1, sysF=F1, sysG=G1),
    'sol_unit': dict(sysD=None, sysE=E2, sysF=None, sysG=None),
    'gas_unit': dict(sysD=None, sysE=None, sysF=None, sysG=None),
    'stream_dct': dict(sysD=streamsD, sysE=streamsE, sysF=streamsF, sysG=streamsG),
    'TEA': dict(sysD=teaD, sysE=teaE, sysF=teaF, sysG=teaG),
    'LCA': dict(sysD=lcaD, sysE=lcaE, sysF=lcaF, sysG=lcaG),
    'cache': dict(sysD={}, sysE={}, sysF={}, sysG={}),
    }
                                                                                                            





def get_summarizing_functions(system):
    func_dct = {}
    func_dct['get_annual_net_cost'] = lambda tea, ppl: (tea.EAC)/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl: tea.annualized_CAPEX/ppl
    func_dct['get_annual_OPEX'] = lambda tea, ppl: tea.AOC/ppl
    ind = 'GlobalWarming'
    func_dct['get_annual_GWP'] = \
        lambda lca, ppl: lca.total_impacts[ind]/lca.lifetime/ppl
    func_dct['get_constr_GWP'] = \
        lambda lca, ppl: lca.total_construction_impacts[ind]/lca.lifetime/ppl
    func_dct['get_trans_GWP'] = \
        lambda lca, ppl: lca.total_transportation_impacts[ind]/lca.lifetime/ppl  
    func_dct['get_stream_items_emission_GWP'] = \
        lambda sys, lca, ppl: (lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind]) \
            /lca.lifetime/ppl
    func_dct['get_other_GWP'] = \
        lambda lca, ppl: lca.total_other_impacts[ind]/lca.lifetime/ppl
        
    #func_dct['NaClO'] = lambda sys: carbon_dict['NaClO'][sys.ID]
    
   

    return func_dct


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )

    for sys in systems:
        func = get_summarizing_functions(sys)
        sys.simulate()
        ppl = sys_dct['ppl'][sys.ID]
        print(f'\n---------- Summary for {sys} ----------\n')


        tea = sys_dct['TEA'][sys.ID]
        tea.show()

        unit = f'{currency}/cap/yr'
        print(f'\nNet cost: {func["get_annual_net_cost"](tea, ppl):.1f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl):.1f} {unit}.')
        print(f'Operating: {func["get_annual_OPEX"](tea, ppl):.1f} {unit}.')
        print('\n')
                                                                                                              
        lca = sys_dct['LCA'][sys.ID]
        lca.show()
        print('\n')
        
        unit = f'{GWP.unit}/cap/yr'
        print(f'\nNet emission: {func["get_annual_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Construction: {func["get_constr_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Transportation: {func["get_trans_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Stream items emission: {func["get_stream_items_emission_GWP"](sys, lca, ppl):.2f} {unit}.')
        print(f'Other: {func["get_other_GWP"](lca, ppl):.2} {unit}.\n')



def save_all_reports():
    import os
    if not os.path.isdir(results_path):
        os.path.mkdir(results_path)
    for i in (sysD, sysE, lcaD, lcaE):
        if isinstance(i, bst.System):
            i.simulate()
            i.save_report(os.path.join(results_path, f'{i.ID}_report.xlsx'))
        else:
            i.save_report(os.path.join(results_path, f'{i.system.ID}_lca.xlsx'))

__all__ = ('sysD','sysE','teaD','teaE','lcaD','lcaE',
           'print_summaries', 'save_all_reports',
           *(i.ID for i in sysD.units),
           *(i.ID for i in sysE.units)
           )
