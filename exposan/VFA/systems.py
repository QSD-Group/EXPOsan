#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 14:46:31 2025

@author: blues
"""

import qsdsan as qs
from qsdsan import WasteStream
from exposan.VFA._components import create_components
from exposan.VFA import _sanunits_copy as su
from qsdsan import sanunits as qsu
import biosteam as bst
from exposan.htl import _tea
# !!! import _tea from VFA folder

bst.PowerUtility.price = 0.0855 
#$/kWh data US total industrial electricity price from EIA
# Average retail price of electricity by sector, 2025 Jan to Aug average for the industry sector
#!!! check electricity price
#%%
# load components and create wastestreams

create_components()

AD_effluent = WasteStream('AD_effluent', Na=13.67, Cl=21.08, K=0, H2O=7743.28, 
                  Propionate=14.68, Butyrate=17.46, Hexanoate=23.02, 
                  Digestate_solids = 233.5, units='kg/hr') # !!! change units to mass per time
# The influent stream composition is estimated from the mass flowrate of the AD effluent in the AD unit in Joy's agile benchmarking paper
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

AC_influent = WasteStream('AC_influent', Na=13.67, Cl=21.08, K=0, H2O=7743.28,
                     Propionate=0, Butyrate=0, Hexanoate=0, units='kg/hr')
# !!! to be changed, accumulating channel influent should be scaled based on the effluent flowrate from the pretreatment
# concentrations taken from raw data: NaCl is 75 mM/min, volumetric flowrate
# is 0.5 ml/min. Assume water density at 25 C. Assume water density at 25 C is 0.9982 g/ml,
# and water molar weight is 18.015 g/mol.



#%%
# Creating SanUnit instances
Centrifuge = su.SolidsSeparation('Centrifuge', ins = (AD_effluent, 'polymer'), outs = ('liquid_stream', 'solids_stream'))
# the second influent of the SolidsSeparation unit is polymer, which is calculated based on the first influent, we would not need to define a separate stream for polymer. Polymer requirement will be calculated from ws0
Centrifuge.ins[1].price = 1.3 / 0.453592 * qs.CEPCI_by_year[2023] / qs.CEPCI_by_year[2014] #unit price of polymer in $/lb from CapdetWorks 2014 database, convert from lb to kg, not including Ferric Chloride cost

Microfiltration = qsu.SludgeThickening('Microfiltration', ins = Centrifuge-0, outs = ('effluent', 'sludge'), sludge_moisture = 0.5)
BasePump = su.BaseDosing('BasePump', ins = (Microfiltration-0, 'base'), outs = ('effluent',))
RedoxED = su.RedoxED('RedoxED', ins = (BasePump-0, AC_influent), outs = ('feeding_channel_out', 'accumulating_channel_out'), voltage = 1)

# Creating system
sys = qs.System.from_units(ID = 'sys', units = [Centrifuge, Microfiltration, BasePump, RedoxED])






    
flowsheet_ID = 'htl_geospatial'

# clear flowsheet and registry for reloading
if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
    getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
    clear_lca_registries()
flowsheet = qs.Flowsheet(flowsheet_ID)
stream = flowsheet.stream
qs.main_flowsheet.set_flowsheet(flowsheet)

_load_components()
_load_process_settings(location=state)

folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))

raw_wastewater = qs.WasteStream('sludge_assumed_in_wastewater', H2O=size, units='MGD', T=25+273.15)
# set H2O equal to the total raw wastewater into the WWTP






# load construction impact


# add stream impact items
# add impact for waste sludge
qs.StreamImpactItem(ID='waste_sludge1_item',
                    linked_stream=stream.sludge,
                    Acidification=0,
                    Ecotoxicity=0,
                    Eutrophication=0,
                    GlobalWarming=-WWTP.ww_2_dry_sludge*waste_GWP/3.79/(10**6),
                    OzoneDepletion=0,
                    PhotochemicalOxidation=0,
                    Carcinogenics=0,
                    NonCarcinogenics=0,
                    RespiratoryEffects=0)

qs.StreamImpactItem(ID='waste_sludge2_item',
                    linked_stream=stream.sludge,
                    Acidification=0,
                    Ecotoxicity=0,
                    Eutrophication=0,
                    GlobalWarming=-WWTP.ww_2_dry_sludge*waste_GWP/3.79/(10**6),
                    OzoneDepletion=0,
                    PhotochemicalOxidation=0,
                    Carcinogenics=0,
                    NonCarcinogenics=0,
                    RespiratoryEffects=0)

# add impact of pH adjustment
qs.StreamImpactItem(ID='NaOH_item',
                    linked_stream=stream.NaOH,
                    Acidification=0.33656,
                    Ecotoxicity=0.77272,
                    Eutrophication=0.00032908,
                    GlobalWarming=1.2514,
                    OzoneDepletion=7.89E-07,
                    PhotochemicalOxidation=0.0033971,
                    Carcinogenics=0.0070044,
                    NonCarcinogenics=13.228,
                    RespiratoryEffects=0.0024543)
    




    

    

    

    

    

_tea.create_tea(sys, IRR_value=0.03, income_tax_value=0.21, finance_interest_value=0.03, labor_cost_value=29.32)

qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
       Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30,
       Cooling=lambda:sys.get_cooling_duty()/1000*30)
    
    # for income tax, 0.35 is the old federal income tax rate
    # in the future, use 0.21 (new federal income tax rate) as the income tax rate
    # if it is necessary to add a state income tax, see exposan.htl.income_tax
sys.simulate()
