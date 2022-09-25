# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from . import Junction
from qsdsan import Process, System
from qsdsan.sanunits import DynamicInfluent


cmps_asm1 = qs.processes.create_asm1_cmps()
thermo_asm1 = qs.get_thermo()
s1 = qs.WasteStream('s1')
for ID in cmps_asm1.IDs: s1.imass[ID] = 1


cmps_adm1 = qs.processes.create_adm1_cmps()
s2 = qs.WasteStream('s2')
for ID in cmps_adm1.IDs: s2.imass[ID] = 2


reactions = []
for n, ID in enumerate(cmps_asm1.IDs):
    if n < 13:
        rxn = {ID: -1,
               cmps_adm1.IDs[2*n]: 0.5,
               cmps_adm1.IDs[2*n+1]: 0.5,
               }
        reactions.append(rxn)
reactions.append({'S_N2': -1})



DI = DynamicInfluent('DI', thermo=thermo_asm1)

J1 = Junction('J1', upstream=DI.outs[0], downstream=s2, reactions=reactions, isdynamic=True)

sys = System('sys', path=(DI, J1,))
sys.set_dynamic_tracker(DI, s2, J1)
sys.simulate(
    state_reset_hook='reset_cache',
    t_span=(0, 10),
    t_eval=np.arange(0, 10.5, 0.5),
    )


# p1 = Process('biomass_convert', 
#              reaction='X_BH + X_ND -> [?]X_pr + [?]X_li + [?]X_ch + [?]X_I',
#              ref_component='X_BH',
#              conserved_for=('COD', 'N', 'P', 'mass'))

p1 = Process('biomass_convert', 
              reaction='X_BH + [?]X_ND -> X_pr + [0.32]X_I',
              ref_component='X_BH',
              conserved_for=('N',))


# %%

# # Probably shouldn't use reactions at all,
# # should directly use the array

# from thermosteam.reaction import Reaction as Rxn, ParallelReaction as PRxn
# cmps_compiled = qs.Components((*cmps_asm1, *cmps_adm1))
# cmps_compiled.compile()
# qs.set_thermo(cmps_compiled)

# lst = []
# for n, ID in enumerate(cmps_asm.IDs):
#     if n < 13:
#         lst.append(Rxn(
#             {
#                 ID: -1,
#                 cmps_adm1.IDs[2*n]: 0.5,
#                 cmps_adm1.IDs[2*n+1]: 0.5,
#             },
#             reactant=ID,
#             X=1,
#             ))
# asm2adm1 = PRxn(lst)