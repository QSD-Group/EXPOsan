# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import qsdsan as qs
from exposan.hap import create_hap_cmps, SBoulardiiFermenter, HApFermenter
from qsdsan.utils import auom

cmps = create_hap_cmps()

p_sta = 68500  # stadium capacity
f_sta = 1/5

n_sch = 114
p_sch = 58000
f_sch = 1/3

q_sch = p_sch * f_sch * 1.4 / n_sch

#%%
population = 8.5e6
N = 100
Q = 1.4/24 * population * 0.01 # L/hr

urine = qs.WasteStream('urine', T=273.15+30)
urine.set_flow_by_concentration(q_sch, 
                                concentrations=dict(
    Urea=16300,
    Cl=5135,
    Na=2780,
    K=1680,
    Creatinine=1410,
    IS=981.5,
    Hippuric_acid=860,
    IP=770,
    Citric_acid=510,
    Glucuronic_acid=475,
    NH3=465,
    Uric_acid=355,
    Other_COD=3497
    ))

CaCl2 = qs.WasteStream('CaCl2')
inocu = qs.WasteStream('inocu')

HF = HApFermenter('HF', tau=60, precipitate_moisture=90,
                  ins=(urine, inocu, CaCl2), outs=['', 'eff', 'precip'])
HF.simulate()
vent, eff, pre = HF.outs

SBF = SBoulardiiFermenter('SBF', design_production_rate=inocu.imass['Yeast']*n_sch)
SBF.simulate()

# sys = qs.System('sys', path=(HF, SBF))
# sys.simulate()