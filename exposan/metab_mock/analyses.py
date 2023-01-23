# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan import ImpactIndicator as IInd, ImpactItem as IItm, TEA, LCA
from qsdsan.utils import ospath
from exposan.metab_mock import create_systems, data_path

IInd.load_from_file(ospath.join(data_path, 'TRACI_indicators.xlsx'), sheet=0)
IItm.load_from_file(ospath.join(data_path, '_impact_items.xlsx'))
fugitive_ch4 = IItm(ID='fugitive_ch4', functional_unit='kg', 
                    GWP100=25, MIR=0.0143794871794872)

#%%
lt = 20
irr = 0.1

sysA, = create_systems(which='A')
sysA.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
# teaA = TEA(sysA, discount_rate=irr, lifetime=lt)
# #!!! add operation electricity, fugitive CH4, biogas offset
# lcaA = LCA(sysA, lifetime=lt)
