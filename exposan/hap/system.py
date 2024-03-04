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
from exposan.hap import create_hap_cmps
from exposan.hap.units import *
from qsdsan.utils import auom


__all__ = ('default_urine_concs',
           'create_system',)
# p_sta = 68500  # stadium capacity
# f_sta = 1/5

# n_sch = 114
# p_sch = 58000
# f_sch = 1/3

# q_sch = p_sch * f_sch * 1.4 / n_sch

default_urine_concs = dict(
    Urea=16300,         # mg/L
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
    )



#%%
# population = 8.5e6
# N = 100
# Q = 1.4/24 * population * 0.01 # L/hr
def create_system(total_pe_served=50000, N_locations=90, urination_rate=1.4,
                  urine_concentrations={}, lifetime=10, 
                  income_tax=0.28, flowsheet=None):
    sys_ID = 'sys'
    
    flowsheet = flowsheet or qs.Flowsheet(sys_ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    cmps = create_hap_cmps()
    
    urine = qs.WasteStream('urine', T=273.15+30)
    concs = urine_concentrations or default_urine_concs.copy()
    Q = total_pe_served * urination_rate / N_locations / 24 # L/hr
    urine.set_flow_by_concentration(Q, concentrations=concs)
    
    CaCl2 = qs.WasteStream('CaCl2')
    inocu = qs.WasteStream('inoculum')

    HF = HApFermenter('HF', N_parallel_HApFermenter=N_locations, 
                      tau=60, precipitate_moisture=90,
                      ins=(urine, inocu, CaCl2), 
                      outs=['', 'eff', 'precipitates'])
    
    YP = YeastProduction('YP', N_parallel_HApFermenter=N_locations)
    PP = PrecipitateProcessing('PP', N_parallel_HApFermenter=N_locations)
    CD = CollectionDistribution('CD', N_parallel_HApFermenter=N_locations)
    sys = qs.System(sys_ID, path=(HF,), facilities=(YP, PP, CD))
    
    qs.TEA(sys, lifetime=lifetime, income_tax=income_tax, 
           system_add_OPEX={'Facility rent':10.96},
           simulate_system=False, CEPCI=801)   

    return sys

#%%
# sys = create_system()
# sys.simulate()
