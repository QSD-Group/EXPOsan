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

#%%
# load components and create wastestreams

create_components()

ws0 = WasteStream('ws0', NaCl=0.0375, KCl=0, H2O=1662.28, 
                  Propionate=0.0125, Butyrate=0.0125, Hexanoate=0.0125, 
                  Digestate_solids = 5, units='mmol/hr')

# =============================================================================
# ws1 = WasteStream('ws1', NaCl=0.0375, KCl=0, H2O=1662.28, 
#                      Propionate=0.0125, Butyrate=0.0125, Hexanoate=0.0125, units='mmol/hr')
# # !!! Double check concentrations taken from raw data: C3C4C6 = 25/25/25 mM,
# # volumetric flowrate is 0.5 ml/min. Assume water density at 25 C is 0.9982 g/ml,
# # and water molar weight is 18.015 g/mol.
# =============================================================================

ws2 = WasteStream('ws2', NaCl=0.0375, KCl=0, H2O=1662.28,
                     Propionate=0, Butyrate=0, Hexanoate=0, units='mmol/hr')
# concentrations taken from raw data: NaCl is 75 mM/min, volumetric flowrate
# is 0.5 ml/min. Assume water density at 25 C. Assume water density at 25 C is 0.9982 g/ml,
# and water molar weight is 18.015 g/mol.



#%%
Centrifuge = su.centrifuge('Centrifuge', ins = ws0, outs = ('liquid_stream', 'solids_stream'))
RedoxED = su.redox_ED('RedoxED', ins = (Centrifuge-0, ws2), outs = ('fc_out', 'ac_out'), voltage = 1)

# RedoxED.register_alias('')

sys = qs.System.from_units(ID = 'sys', units = [Centrifuge, RedoxED])
sys.simulate()
