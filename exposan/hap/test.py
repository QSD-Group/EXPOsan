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

cmps = create_hap_cmps()

# Q = 10 # L/hr

# urine = qs.WasteStream('urine')
# urine.set_flow_by_concentration(Q, concentrations=dict(
#     Urea=16300,
#     Cl=5135,
#     Na=2780,
#     K=1680,
#     Creatinine=1410,
#     IS=981.5,
#     Hippuric_acid=860,
#     IP=770,
#     Citric_acid=510,
#     Glucuronic_acid=475,
#     NH3=465,
#     Uric_acid=355,
#     Other_COD=3497
#     ))

# CaCl2 = qs.WasteStream('CaCl2')
# inocu = qs.WasteStream('inocu')

# HF = HApFermenter('HF', ins=(urine, inocu, CaCl2), outs=['', 'eff', 'precip'])

SBF = SBoulardiiFermenter('SBF')
