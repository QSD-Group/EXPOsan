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

bst.PowerUtility.price = 0.2 #$/kWh
#!!! check electricity price
#%%
# load components and create wastestreams

create_components()

ws0 = WasteStream('ws0', Na=0.0375, Cl=0.0375, K=0, H2O=1662.28, 
                  Propionate=0.0125, Butyrate=0.0125, Hexanoate=0.0125, 
                  Digestate_solids = 5, units='mmol/hr')

# =============================================================================
# ws0 = WasteStream('ws1', NaCl=0.0375, KCl=0, H2O=1662.28, 
#                      Propionate=0.0125, Butyrate=0.0125, Hexanoate=0.0125, units='mmol/hr')
# # !!! Double check concentrations taken from raw data: C3C4C6 = 25/25/25 mM,
# # volumetric flowrate is 0.5 ml/min. Assume water density at 25 C is 0.9982 g/ml,
# # and water molar weight is 18.015 g/mol.
# =============================================================================

ws1 = WasteStream('ws1', Na=0.0375, Cl=0.0375, K=0, H2O=1662.28,
                     Propionate=0, Butyrate=0, Hexanoate=0, units='mmol/hr')
# concentrations taken from raw data: NaCl is 75 mM/min, volumetric flowrate
# is 0.5 ml/min. Assume water density at 25 C. Assume water density at 25 C is 0.9982 g/ml,
# and water molar weight is 18.015 g/mol.



#%%
# Creating SanUnit instances
Centrifuge = su.SolidsSeparation('Centrifuge', ins = (ws0, 'polymer'), outs = ('liquid_stream', 'solids_stream'))
Centrifuge.ins[1].price = 1.3 / 0.453592 * qs.CEPCI_by_year[2023] / qs.CEPCI_by_year[2014] #unit price of polymer in $/lb from CapdetWorks 2014 database, convert from lb to kg, not including Ferric Chloride cost
Microfiltration = qsu.SludgeThickening('Microfiltration', ins = Centrifuge-0, outs = ('eff', 'sludge'), sludge_moisture = 0.5)
BasePump = su.BaseDosing('BasePump', ins = (Microfiltration-0, 'base'), outs = ('eff',))
RedoxED = su.RedoxED('RedoxED', ins = (BasePump-0, ws1), outs = ('fc_out', 'ac_out'), voltage = 1)

# Creating system
sys = qs.System.from_units(ID = 'sys', units = [Centrifuge, Microfiltration, BasePump, RedoxED])
sys.simulate()
