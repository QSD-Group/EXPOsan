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

import thermosteam as tmo
from qsdsan import Component, Components

__all__ = ('components', )

Mg = Component.from_chemical('Mg', search_ID='Mg',
                              phase='l', particle_size='Soluble',
                              degradability='Undegradable', organic=False)

Ca = Component.from_chemical('Ca', search_ID='Ca',
                              phase='l', particle_size='Soluble',
                              degradability='Undegradable', organic=False)

H2O = Component('H2O', phase='l', particle_size='Soluble',
                degradability='Undegradable', organic=False)

def add_V_from_rho(cmp, rho):
    V_model = tmo.functional.rho_to_V(rho, cmp.MW)
    try: cmp.V.add_model(V_model)
    except:
        handle = getattr(cmp.V, cmp.locked_state)
        handle.add_model(V_model)

NaClO = Component('NaClO', phase='l', particle_size = 'Soluble', 
               degradability = 'Undegradable', organic=False)

Ecoli = Component('Ecoli', MW = 1, phase = 's', particle_size = 'Particluate', 
                  degradability = 'Undegradable', organic = False)

add_V_from_rho(Ecoli, 500)
Ecoli.copy_models_from(tmo.Chemical('Glucose'), ('Cn', 'mu'))

Polyethylene = Component('Polyethylene', formula='C2H4',
                         phase='s', particle_size='Particulate',
                         degradability='Slowly', organic=True)

# PVC = Component('PVC', formula='H2C–CHCl',
#                 phase='s', particle_size='Particulate',
#                 degradability='Undegradable', organic=True)

# Quartz = Component('Quartz', formula='SiO2',
#                   phase='s', particle_size='Particulate', 
#                   degradability='Undegradable', organic=False)

AgNp = Component('AgNp', search_ID='Ag',
                 phase='l', particle_size='Particulate',
                 degradability='Undegradable', organic=False)

for cmp in (Mg, Ca,):
    cmp.default()
    cmp.copy_models_from(H2O, ('sigma', 'epsilon', 'kappa', 'Cn', 'mu'))
    add_V_from_rho(cmp, 1e3) # assume the same density as water

components = Components((Mg, Ca, H2O, NaClO, Ecoli, Polyethylene, AgNp))

for i in components:
    if i.HHV is None: i.HHV = 0
    if i.LHV is None: i.LHV = 0
    if i.Hf is None: i.Hf = 0

components.compile()
components.set_alias('H2O', 'Water')
