#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 14:46:31 2025

@author: blues
"""

import qsdsan as qs
from qsdsan import WasteStream
from qsdsan.utils import ospath, data_path, misc
from exposan.VFA._components import create_components
from exposan.VFA import _sanunits_copy as su
import biosteam as bst
from exposan.htl import _tea

bst.PowerUtility.price = 0.0855 
#$/kWh data US total industrial electricity price from EIA
# Average retail price of electricity by sector, 2025 Jan to Aug average for the industry sector
#!!! check electricity price
#%%
# load components and create wastestreams

create_components()

AD_effluent = WasteStream('AD_effluent', Na=13.67, Cl=21.10, K=0, H2O=7743.28, 
                  Propionate=14.68, Butyrate=17.46, Hexanoate=23.02, 
                  Digestate_solids = 412.274, units='kg/hr') # !!! change units to mass per time
# The influent stream composition is estimated from the mass flowrate of the AD effluent in the AD unit in Joy's agile benchmarking paper
# assuming 5% TSS content
# Propionate MW = 74.08 g/mol
# Butyrate MW = 88.11 g/mol
# Hexanoate MW = 116.1583 g/mol
# =============================================================================
# ws0 = WasteStream('ws1', NaCl=0.0375, KCl=0, H2O=1662.28, 
#                      Propionate=0.0125, Butyrate=0.0125, Hexanoate=0.0125, units='mmol/hr')
# # !!! Double check concentrations taken from raw data: C3C4C6 = 25/25/25 mM,
# # volumetric flowrate is 0.5 ml/min. Assume water density at 25 C is 0.9982 g/ml,
# # and water molar weight is 18.015 g/mol.
# =============================================================================

AC_influent = WasteStream('AC_influent', Na=13.67, Cl=21.10, K=0, H2O=7743.28,
                     Propionate=0, Butyrate=0, Hexanoate=0, units='kg/hr')
# !!! to be changed, accumulating channel influent should be scaled based on the effluent flowrate from the pretreatment
# concentrations taken from raw data: NaCl is 75 mM/min, volumetric flowrate
# is 0.5 ml/min. Assume water density at 25 C. Assume water density at 25 C is 0.9982 g/ml,
# and water molar weight is 18.015 g/mol.



#%%
flowsheet_ID = 'vfa_separation'
# clear flowsheet and registry for reloading
if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
    getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
    misc.clear_lca_registries()
flowsheet = qs.Flowsheet(flowsheet_ID)
stream = flowsheet.stream
qs.main_flowsheet.set_flowsheet(flowsheet)


### load construction impact
qs.ImpactIndicator.load_from_file(ospath.join(data_path, 'sanunit_data/VFA/impact_indicators.csv'))
qs.ImpactItem.load_from_file(ospath.join(data_path, 'sanunit_data/VFA/impact_items.xlsx'))


### Creating SanUnit instances
Centrifuge = su.SolidsSeparation('Centrifuge', ins = (AD_effluent, 'polymer'), outs = ('liquid_stream', 'solids_stream'))
# the second influent of the SolidsSeparation unit is polymer, which is calculated based on the first influent, we would not need to define a separate stream for polymer. Polymer requirement will be calculated from ws0
Centrifuge.ins[1].price = 1.3 / 0.453592 * qs.CEPCI_by_year[2023] / qs.CEPCI_by_year[2014] #unit price of polymer in $/lb from CapdetWorks 2014 database, convert from lb to kg, not including Ferric Chloride cost

Microfiltration = su.SludgeThickening('Microfiltration', ins = Centrifuge-0, outs = ('mf_effluent', 'sludge'), sludge_moisture = 0.5)
BasePump = su.BaseDosing('BasePump', ins = (Microfiltration-0, AD_effluent.copy('AD_ref'), 'base'), outs = ('bp_effluent',))
RedoxED = su.RedoxED('RedoxED', ins = (BasePump-0, AC_influent), outs = ('feeding_channel_out', 'accumulating_channel_out'), voltage = 0.8)
### Match accumulating-channel influent volume to feeding-channel flow
@RedoxED.add_specification(run=True)
def match_ac_to_fc_flow():
    fc = BasePump-0
    ac = AC_influent
    print("Spec called: FC flow =", fc.F_vol, "AC flow before =", ac.F_vol)
    if ac.F_vol > 0:
        ac.scale(fc.F_vol / ac.F_vol)
    print("AC flow after =", ac.F_vol)


### Creating system
sys = qs.System.from_units(ID = 'sys', units = [Centrifuge, Microfiltration, BasePump, RedoxED])


### add stream impact items
# add impact for polymer
qs.StreamImpactItem(ID='polymer_item',
                    linked_stream=stream.polymer,
                    GlobalWarming=5.19)

# add impact for waste sludge
qs.StreamImpactItem(ID='waste_sludge1_item',
                    linked_stream=stream.solids_stream,
                    GlobalWarming=0.292)

qs.StreamImpactItem(ID='waste_sludge2_item',
                    linked_stream=stream.sludge,
                    GlobalWarming=0.292)

# add impact of pH adjustment
qs.StreamImpactItem(ID='NaOH_item',
                    linked_stream=stream.base,
                    GlobalWarming=1.2514)
    

tea = _tea.create_tea(sys, IRR_value=0.03, income_tax_value=0.21, finance_interest_value=0.03, labor_cost_value = 29.32, duration = (2025, 2035))

sys.simulate()

Electricity = qs.ImpactItem('Electricity', 'kWh', GWP=1.1) #!!! include uncertainty range for electricity GWP in the technical assumption data sheet
get_power = sum([u.power_utility.rate for u in sys.units]) * (24 * 365 * 10)

lca = qs.LCA(system=sys, lifetime=10, lifetime_unit='yr',
       Electricity = get_power)
    
    # for income tax, 0.35 is the old federal income tax rate
    # in the future, use 0.21 (new federal income tax rate) as the income tax rate
    # if it is necessary to add a state income tax, see exposan.htl.income_tax

