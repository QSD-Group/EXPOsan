#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 21:19:00 2022

@author: yalinli_cabbi
"""
import numpy as np, pandas as pd, qsdsan as qs
from exposan.metab_mock import (
    create_modelB,
    create_modelC,
    create_systems,
    R1_ss_conds,
    R2_ss_conds,
    )

mdlB = create_modelB()
sysB = mdlB.system
sysB.simulate(state_reset_hook='reset_cache', t_span=(0,200), method='BDF')


mdlC = create_modelC()
sysC = mdlC.system
fC = sysC.flowsheet
inf = fC.stream.BreweryWW_C
eff = fC.stream.Effluent_C
bgs = [s for s in sysC.products if 'biogas' in s.ID]

R1 = fC.unit.R1
R2 = fC.unit.R2
R1.set_init_conc(**R1_ss_conds)
R2.set_init_conc(**R2_ss_conds)

DM1_c = fC.unit.DM1_c
DM2_c = fC.unit.DM2_c

for u in (DM1_c, DM2_c):
    for gas in ('H2', 'CH4', 'CO2'):
        setattr(u, f'{gas}_degas_efficiency', 1)

cmps = qs.get_components()
def get_COD(streams):
    try: iter(streams)
    except: streams = [streams]
    return sum((s.mass*cmps.i_COD).sum() for s in streams)

# # Parameter names
# H2E sidestream split
# CH4E sidestream split
# H2 removal efficiency
# CH4 removal efficiency
# CO2 removal efficiency

def get_p(mdl, name):
    for p in mdlC.parameters:
        if p.name == name: return p
p_R1_split = get_p(mdlC, 'H2E sidestream split')
p_R2_split = get_p(mdlC, 'CH4E sidestream split')

# # Metric names
# SRT
# R1 VFAs
# R1 H2
# R2 CH4
# Effluent COD
# Total COD removal
# H2 production
# CH4 production

def get_m(mdl, name):
    for m in mdl.metrics:
        if m.name == name: return m
mB = get_m(mdlB, 'Total COD removal')
mC = get_m(mdlC, 'Total COD removal')


dcts = cod_rm_dct, cod_eff_dct, cod_bgs_dct = [{}, {}, {}]
splits = [round(i, 1) for i in np.arange(0.1, 1, 0.1)]
for R1_split in splits:
    p_R1_split.setter(R1_split)
    cod_rm_dct[R1_split] = []
    cod_eff_dct[R1_split] = []
    cod_bgs_dct[R1_split] = []
    for R2_split in splits:
        p_R2_split.setter(R2_split)
        R1.set_init_conc(**R1_ss_conds)
        R2.set_init_conc(**R2_ss_conds)
        try: 
            sysC.simulate(state_reset_hook='reset_cache', t_span=(0,200), method='BDF')
            cod_rm_dct[R1_split].append(mC())
            cod_eff_dct[R1_split].append(get_COD(eff))
            cod_bgs_dct[R1_split].append(get_COD(bgs))
        except: 
            for dct in dcts: dct[R1_split].append(np.nan)
  
dfs = [pd.DataFrame.from_dict(dct) for dct in dcts]
for df in dfs:    
    df.columns.name = 'R1 split'
    df.index = splits
    df.index.name = 'R2 split'
    

# R1_split = 0.1
# p_R1_split.setter(R1_split)
# lst = []
# splits = [round(i, 1) for i in np.arange(0.1, 1, 0.1)]
# for R2_split in splits:
#     p_R2_split.setter(R2_split)
#     R1.set_init_conc(**R1_ss_conds)
#     R2.set_init_conc(**R2_ss_conds)
#     try: 
#         sysC.simulate(state_reset_hook='reset_cache', t_span=(0,200), method='BDF')
#         removal = m()
#     except: removal = np.nan
#     print(removal)
#     lst.append(removal)
    
    
# %%

import numpy as np, pandas as pd, qsdsan as qs
from exposan.metab_mock import create_systems, create_modelC

sysA, sysB, sysC = create_systems()

fB = sysB.flowsheet
fC = sysC.flowsheet
inf = fC.stream.BreweryWW_C
cmps = inf.components
eff = fC.stream.Effluent_C
bgs = [s for s in sysC.products if 'biogas' in s.ID]

def get_COD(streams):
    try: iter(streams)
    except: streams = [streams]
    return sum((s.mass*cmps.i_COD).sum() for s in streams)

R1 = fC.unit.R1
R2 = fC.unit.R2

R2.split = (0.1, 0.9)

splits = np.linspace(0.1, 0.9, 9)

inf_CODs = []
eff_CODs = []
bgs_CODs = []
for split in splits:
    R1.split = (split, 1-split)
    # sysC.set_tolerance(rmol=1e-6)
    try:
        sysC.simulate(state_reset_hook='reset_cache', t_span=(0,200), method='BDF')
        inf_CODs.append(get_COD(inf))
        eff_CODs.append(get_COD(eff))
        bgs_CODs.append(get_COD(bgs))
    except:
        inf_CODs.append(np.nan)
        eff_CODs.append(np.nan)
        bgs_CODs.append(np.nan)

rms =[1-i/j for i, j in zip(eff_CODs, inf_CODs)]
outs_CODs = [i+j for i, j in zip(eff_CODs, bgs_CODs)]
balances = [i/j for i, j in zip(outs_CODs, inf_CODs)]


# %%
from exposan.metab_mock import *

sysA, sysB, sysC = create_systems()
dct = globals()
for sys in (sysA, sysB, sysC): dct.update(sys.flowsheet.to_dict())

cmps = R1.components
def get_COD(streams):
    try: iter(streams)
    except: streams = [streams]
    return sum((s.mass*cmps.i_COD).sum() for s in streams)

t = 200
# t = 20 # shorter t doesn't make sense
# t = 2000 # long t doesn't make a difference (sysB is about to close COD MB, sysC a little off)

split = 0.9
R1.split = R2.split = (split, 1-split)

degassing = 1
for u in (DM1_c, DM2_c):
    for gas in ('H2', 'CH4', 'CO2'):
        setattr(u, f'{gas}_degas_efficiency', degassing)
# DM1_c.CO2_degas_efficiency = DM2_c.CO2_degas_efficiency = 0
        

sysB.simulate(state_reset_hook='reset_cache', t_span=(0,t), method='BDF')
sysC.simulate(state_reset_hook='reset_cache', t_span=(0,t), method='BDF')

get_bgs = lambda sys: [s for s in sys.products if 'biogas' in s.ID]

print(f'\nt is {t} d, split is {split}, degassing is {degassing}:')
COD_effB = get_COD(Effluent_B)
COD_effC = get_COD(Effluent_C)
print('sysB eff COD: ', COD_effB)
print('sysC eff COD: ', COD_effC)

print('sysB biogas COD: ', get_COD(get_bgs(sysB)))
print('sysC biogas COD: ', get_COD(get_bgs(sysC)))

COD_in = get_COD(sysB.feeds)
assert COD_in == get_COD(sysC.feeds)
print('tot COD in: ', COD_in)
COD_outB = get_COD(sysB.products)
COD_outC = get_COD(sysC.products)
print('sysB tot COD out: ', COD_outB)
print('sysC tot COD out: ', COD_outC)
print('sysB tot COD recovery: ', COD_outB/COD_in)
print('sysC tot COD recovery: ', COD_outC/COD_in)
print('sysB tot COD removal: ', 1-COD_effB/COD_in)
print('sysC tot COD removal: ', 1-COD_effC/COD_in)

print(f'R1 MB: {get_COD(R1.outs)/get_COD(R1.ins)}')
print(f'R2 MB: {get_COD(R2.outs)/get_COD(R2.ins)}')
print(f'DM1_c MB: {get_COD(DM1_c.outs)/get_COD(DM1_c.ins)}')
print(f'DM2_c MB: {get_COD(DM2_c.outs)/get_COD(DM2_c.ins)}')


# %%

from qsdsan import sanunits as su, WasteStream, processes as pc, System
from exposan.metab_mock import *

Q = 5           # influent flowrate [m3/d]
T1 = 273.15+35  # temperature [K]
Vl1 = 5         # liquid volume [m^3]
Vg1 = 0.556     # headspace volume [m^3]
split_1 = 0.75  # split ratio to side-stream
tau_1 = 0.021   # degassing membrane retention time [d]

T2 = 273.15+25    
Vl2 = 75
Vg2 = 5
split_2 = 0.75
tau_2 = 0.021

fermenters = ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro')
methanogens = ('X_ac', 'X_h2')
biomass_IDs = (*fermenters, *methanogens)

adm1 = pc.ADM1()

inf_b = WasteStream('BreweryWW_B', T=T1)
inf_b.set_flow_by_concentration(Q, concentrations=default_inf_concs, units=('m3/d', 'kg/m3'))
eff_B = WasteStream('Effluent_B', T=T2)
bg1_B = WasteStream('biogas_1B', phase='g')
bg2_B = WasteStream('biogas_2B', phase='g')
AnR1 = su.AnaerobicCSTR('AnR1', ins=inf_b, outs=(bg1_B, ''), 
                        V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                        retain_cmps=fermenters)
# AnR2 = su.AnaerobicCSTR('AnR2', ins=AnR1-1, outs=(bg2_B, eff_B), 
#                         V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
#                         retain_cmps=methanogens)
AnR1.set_init_conc(**R1_ss_conds)
# AnR2.set_init_conc(**R2_ss_conds)
sysB = System('baseline', path=(AnR1,))
sysB.set_dynamic_tracker(AnR1, bg1_B)
# sysB = System('baseline', path=(AnR1, AnR2))
# sysB.set_dynamic_tracker(AnR1, AnR2, bg1_B, bg2_B)

inf_c = inf_b.copy('BreweryWW_C')
eff_c = WasteStream('Effluent_C', T=T2)
bgm1 = WasteStream('biogas_mem_1', phase='g')
bgm2 = WasteStream('biogas_mem_2', phase='g')
bgh1 = WasteStream('biogas_hsp_1', phase='g')
bgh2 = WasteStream('biogas_hsp_2', phase='g')

############# sysC unit operation #################
sc1 = 0.9
sc2 = 0.9
R1 = su.AnaerobicCSTR('R1', ins=[inf_c, 'return_1'], 
                      outs=(bgh1, 'sidestream_1', ''), 
                      split=(sc1, 1-sc1),
                      V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                      retain_cmps=fermenters)
DM1c = DM1_c = DegassingMembrane('DM1_c', ins=R1-1, outs=(bgm1, 1-R1), tau=tau_1)
# DM1c = DM('DM1_c', ins=R1-1, outs=(bgm1, 1-R1), tau=0.1)    

# R2 = su.AnaerobicCSTR('R2', ins=[R1-2, 'return_2'], 
#                       outs=(bgh2, 'sidestream_2', eff_c), 
#                       split=(sc2, 1-sc2),
#                       V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
#                       retain_cmps=methanogens)
# DM2c = DM2_c = DegassingMembrane('DM2_c', ins=R2-1, outs=(bgm2, 1-R2), tau=tau_2)
# DM2c = DM('DM2_c', ins=R2-1, outs=(bgm2, 1-R2), tau=0.1)

R1.set_init_conc(**R1_ss_conds)
# R2.set_init_conc(**R2_ss_conds)

degassing = 1
for u in (DM1_c, ):
    for gas in ('H2', 'CH4', 'CO2'):
        setattr(u, f'{gas}_degas_efficiency', degassing)
# DM1_c.H2_degas_efficiency = DM2_c.H2_degas_efficiency = 1
# DM1_c.CH4_degas_efficiency = DM2_c.CH4_degas_efficiency = 1

sysC = System('combined_METAB', path=(R1, DM1c,),
              recycle=(DM1c-1,))
sysC.set_dynamic_tracker(R1, bgm1, bgh1,)

# sysC = System('combined_METAB', path=(R1, DM1c, R2, DM2c),
#               recycle=(DM1c-1, DM2c-1))
# sysC.set_dynamic_tracker(R1, R2, bgm1, bgm2, bgh1, bgh2)



cmps = AnR1.components
def get_COD(streams):
    try: iter(streams)
    except: streams = [streams]
    return sum((s.mass*cmps.i_COD).sum() for s in streams)

t = 200
# sysB.simulate(state_reset_hook='reset_cache', t_span=(0,t), method='BDF')
sysC.simulate(state_reset_hook='reset_cache', t_span=(0,t), method='BDF')

#%%
import exposan.metab_mock as mm
sysA, sysB, sysC = mm.create_systems()
# sysC._setup()
# sysC.converge()
cmps = sysC.units[0].components
sysC.set_tolerance(rmol=1e-6)
sysC.simulate(t_span=(0,200), method='BDF', state_reset_hook='reset_cache')
inf_cod = sum([(ws.mass * cmps.i_COD).sum() for ws in sysC.feeds])
out_cod = sum([(ws.mass * cmps.i_COD).sum() for ws in sysC.products])
out_cod/inf_cod
