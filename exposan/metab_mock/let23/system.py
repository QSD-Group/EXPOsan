# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import sanunits as su, processes as pc, WasteStream, System
from qsdsan.utils import ExogenousDynamicVariable as EDV, auom
# from chemicals.elements import molecular_weight as get_mw
from exposan.metab_mock import DegassingMembrane as DM, rhos_adm1_ph_ctrl

__all__ = ('create_systems',)

#%% Global Variables

from exposan.metab_mock.systems import (
    Q, T1, Vl1, Vg1, T2, Vl2, Vg2,
    fermenters, methanogens, biomass_IDs, #vfa_IDs,
    default_inf_concs, R1_ss_conds, R2_ss_conds,
    # yields_bl, mus_bl, Ks_bl,
    )

# reactor water depths [m]

h1 = 4.5
h2 = 6

q = Q*auom('m3/d').conversion_factor('ft3/s')
L = 5 # pipe length [ft]
D = 3*auom('inch').conversion_factor('ft')
C = 110 # the Hazen-Williams coefficient
hf = 3.02*L*(4*q/C/np.pi)*D**(-4.87)*auom('ft').conversion_factor('m') # friction head 
dP1 = 1000 * 9.8 * (h1+hf*2)
dP2 = 1000 * 9.8 * (h2-h1+hf*2)
dPs = 1000 * 9.8 * hf*2
dPm = 6*auom('psi').conversion_factor('Pa')  # Hofi membrane water pressure drop at median flowrate through lumen

ph1 = 5.8
ph2 = 7.2

Temp1 = EDV('T1', function=lambda t: T1)
Temp2 = EDV('T2', function=lambda t: T2)
pH1 = EDV('pH1', function=lambda t: ph1)
pH2 = EDV('pH2', function=lambda t: ph2)

#%% Systems

def create_systems(flowsheet_A=None, flowsheet_B=None, flowsheet_C=None, flowsheet_D=None, 
                   inf_concs={}, R1_init_conds={}, R2_init_conds={}, which=None,
                   selective=True):
    
    which = which or ('A', 'B', 'C', 'D')
    if isinstance(which, str): which = (which,)
    
    if selective: 
        R1_retain = fermenters
        R2_retain = methanogens
    else: R1_retain = R2_retain = biomass_IDs
    pc.create_adm1_cmps()
    adm1 = pc.ADM1()
    dyn_params = adm1.rate_function.params.copy()
    adm1.set_rate_function(rhos_adm1_ph_ctrl)
    adm1.rate_function._params = dyn_params
    inf_concs = inf_concs or default_inf_concs.copy()
    brewery_ww = WasteStream('BreweryWW_A')
    brewery_ww.set_flow_by_concentration(Q, concentrations=inf_concs, units=('m3/d', 'kg/m3'))
    R1_init_conds = R1_init_conds or R1_ss_conds
    R2_init_conds = R2_init_conds or R2_ss_conds
    systems = []
    
    if 'A' in which:
        flowsheet_A = flowsheet_A or qs.Flowsheet('A')
        qs.main_flowsheet.set_flowsheet(flowsheet_A)
        flowsheet_A.stream.register(brewery_ww.ID, brewery_ww)
        
        effA = WasteStream('Effluent_A', T=T2)
        bg1A = WasteStream('Biogas_1A', phase='g')
        bg2A = WasteStream('Biogas_2A', phase='g')
        
        # adm1A = pc.ADM1()
        # dyn_params = adm1A.rate_function.params.copy()
        # adm1A.set_rate_function(rhos_adm1_ph_ctrl)
        # adm1A.rate_function._params = dyn_params
        P1A = su.Pump('P1A', ins=brewery_ww, dP_design=dP1, 
                      init_with='WasteStream', isdynamic=True)
        R1A = su.AnaerobicCSTR('R1A', ins=P1A-0, outs=(bg1A, ''), 
                               V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                               retain_cmps=R1_retain, 
                               exogenous_vars=(Temp1, pH1))
        P2A = su.Pump('P2A', ins=R1A-1, dP_design=dP2, 
                      init_with='WasteStream', isdynamic=True)
        R2A = su.AnaerobicCSTR('R2A', ins=P2A-0, outs=(bg2A, effA), 
                               V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                               retain_cmps=R2_retain, 
                               exogenous_vars=(Temp2, pH2))
    
        R1A.set_init_conc(**R1_init_conds)
        R2A.set_init_conc(**R2_init_conds)
    
        sysA = System('sysA', path=(P1A, R1A, P2A, R2A))
        sysA.set_dynamic_tracker(R1A, R2A, bg1A, bg2A, effA)
        
        systems.append(sysA)
    
    if 'B' in which:
        flowsheet_B = flowsheet_B or qs.Flowsheet('B')
        qs.main_flowsheet.set_flowsheet(flowsheet_B)
        
        infB = brewery_ww.copy('BreweryWW_B')
        effB = WasteStream('Effluent_B', T=T2)
        bg1B = WasteStream('Biogas_1B', phase='g')
        bg2B = WasteStream('Biogas_2B', phase='g')

        # adm1B = pc.ADM1()
        # dyn_params = adm1B.rate_function.params.copy()
        # adm1B.set_rate_function(rhos_adm1_ph_ctrl)
        # adm1B.rate_function._params = dyn_params
        P1B = su.Pump('P1B', ins=infB, dP_design=dP1, 
                      init_with='WasteStream', isdynamic=True)
        R1B = su.AnaerobicCSTR('R1B', ins=P1B-0, outs=('', ''), 
                               V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                               fixed_headspace_P=True,
                               headspace_P=0.1013,
                               retain_cmps=R1_retain, 
                               exogenous_vars=(Temp1, pH1))
        VP1B = su.Pump('VP1B', ins=R1B-0, outs=bg1B, P=101325)
        P2B = su.Pump('P2B', ins=R1B-1, dP_design=dP2,
                      init_with='WasteStream', isdynamic=True)
        R2B = su.AnaerobicCSTR('R2B', ins=P2B-0, outs=(bg2B, effB), 
                               V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                               retain_cmps=R2_retain, 
                               exogenous_vars=(Temp2, pH2))
    
        R1B.set_init_conc(**R1_init_conds)
        R2B.set_init_conc(**R2_init_conds)
        
        sysB = System('sysB', path=(P1B, R1B, VP1B, P2B, R2B))
        sysB.set_dynamic_tracker(R1B, R2B, bg1B, bg2B, effB)
        
        systems.append(sysB)
        
    if 'C' in which:
        flowsheet_C = flowsheet_C or qs.Flowsheet('C')
        qs.main_flowsheet.set_flowsheet(flowsheet_C)
        
        infC = brewery_ww.copy('BreweryWW_C')
        effC = WasteStream('Effluent_C', T=T2)
        bg1C = WasteStream('Biogas_1C', phase='g')
        bgh2C = WasteStream('Biogas_hsp_2C', phase='g')
        bgm2C = WasteStream('Biogas_mem_2C', phase='g')
        
        # adm1C = pc.ADM1()
        # dyn_params = adm1C.rate_function.params.copy()
        # adm1C.set_rate_function(rhos_adm1_ph_ctrl)
        # adm1C.rate_function._params = dyn_params
        P1C = su.Pump('P1C', ins=infC, dP_design=dP1, 
                      init_with='WasteStream', isdynamic=True)
        R1C = su.AnaerobicCSTR('R1C', ins=P1C-0, outs=(bg1C, ''), 
                                V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1,
                                retain_cmps=R1_retain, 
                                exogenous_vars=(Temp1, pH1))
        P2C = su.Pump('P2C', ins=R1C-1, dP_design=dP2,
                      init_with='WasteStream', isdynamic=True)
        R2C = su.AnaerobicCSTR('R2C', ins=P2C-0, outs=(bgh2C, ''), 
                                V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                                retain_cmps=R2_retain, 
                                exogenous_vars=(Temp2, pH2))
        P2sC = su.Pump('P2sC', ins=R2C-1, dP_design=dPm+dPs, 
                      init_with='WasteStream', isdynamic=True)
        DM2C = DM('DM2C', ins=P2sC-0, outs=('', effC))
        VP2C = su.Pump('VP2C', ins=DM2C-0, outs=bgm2C, dP_design=101325)
    
        R1C.set_init_conc(**R1_init_conds)
        R2C.set_init_conc(**R2_init_conds)
        
        sysC = System('sysC', path=(P1C, R1C, P2C, R2C, P2sC, DM2C, VP2C))
        sysC.set_dynamic_tracker(R1C, R2C, bg1C, bgh2C, bgm2C, effC)
        
        systems.append(sysC)   
    
    if 'D' in which:
        flowsheet_D = flowsheet_D or qs.Flowsheet('D')
        qs.main_flowsheet.set_flowsheet(flowsheet_D)
        
        infD = brewery_ww.copy('BreweryWW_D')
        effD = WasteStream('Effluent_D', T=T2)
        bgh1D = WasteStream('Biogas_hsp_1D', phase='g')
        bgm1D = WasteStream('Biogas_mem_1D', phase='g')
        bgh2D = WasteStream('Biogas_hsp_2D', phase='g')
        bgm2D = WasteStream('Biogas_mem_2D', phase='g')
        
        # adm1D = pc.ADM1()
        # dyn_params = adm1D.rate_function.params.copy()
        # adm1D.set_rate_function(rhos_adm1_ph_ctrl)
        # adm1D.rate_function._params = dyn_params
        P1D = su.Pump('P1D', ins=infD, dP_design=dP1, 
                      init_with='WasteStream', isdynamic=True)
        R1D = su.AnaerobicCSTR('R1D', ins=[P1D-0, 'Return_1D'], 
                                outs=(bgh1D, 'Sidestream_1D', ''), 
                                split=[400, 1],
                                V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                                retain_cmps=R1_retain, 
                                exogenous_vars=(Temp1, pH1))
        P1sD = su.Pump('P1sD', ins=R1D-1, dP_design=dPm+dPs, 
                      init_with='WasteStream', isdynamic=True)
        DM1D = DM('DM1D', ins=P1sD-0, outs=('', 1-R1D))
        VP1D = su.Pump('VP1D', ins=DM1D-0, outs=bgm1D, dP_design=101325)
        P2D = su.Pump('P2D', ins=R1D-2, dP_design=dP2, 
                      init_with='WasteStream', isdynamic=True)
        R2D = su.AnaerobicCSTR('R2D', ins=P2D-0, outs=(bgh2D, ''), 
                                V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                                retain_cmps=R2_retain, 
                                exogenous_vars=(Temp2, pH2))
        P2sD = su.Pump('P2sD', ins=R2D-1, dP_design=dPm+dPs, 
                      init_with='WasteStream', isdynamic=True)        
        DM2D = DM('DM2D', ins=P2sD-0, outs=('', effD))
        VP2D = su.Pump('VP2D', ins=DM2D-0, outs=bgm2D, dP_design=101325)
    
        R1D.set_init_conc(**R1_init_conds)
        R2D.set_init_conc(**R2_init_conds)
    
        sysD = System('sysD', path=(P1D, R1D, P1sD, DM1D, VP1D, 
                                    P2D, R2D, P2sD, DM2D, VP2D), 
                      recycle=(DM1D-1,))
        sysD.set_dynamic_tracker(R1D, R2D, bgh1D, bgm1D, bgh2D, bgm2D, effD)
        
        systems.append(sysD)
    
    return systems


#%%
if __name__ == '__main__':
    sysA, sysB, sysC, sysD = systems = create_systems()
    for sys in systems:
        sys.simulate(
            t_span=(0,400),
            state_reset_hook='reset_cache',
            method='BDF'
            )
        sys.diagram()
    au = sysA.flowsheet.unit
    bu = sysB.flowsheet.unit
    cu = sysC.flowsheet.unit
    du = sysD.flowsheet.unit

#%%
# cmps = qs.get_thermo().chemicals
# #!!! should we include dissolved gas in COD removal calculation?
# def tCOD(streams):
#     # return sum((s.mass*cmps.i_COD).sum() for s in streams)
#     return sum(s.composite('COD', flow=True) for s in streams)

# def lcod(streams):
#     # return sum((s.mass*cmps.i_COD).sum() for s in streams if s.phase == 'l')
#     return sum(s.composite('COD', flow=True) for s in streams if s.phase == 'l')

# inf_cod = [tCOD(sys.feeds) for sys in systems]
# out_cod = [tCOD(sys.products) for sys in systems]
# balance = [out/inf for inf, out in zip(inf_cod, out_cod)]

# eff_cod = [lcod(sys.products) for sys in systems]
# rcod = [1-eff/inf for inf, eff in zip(inf_cod, eff_cod)]

# rcod = [1-eff/inf for inf, eff in zip(out_cod, eff_cod)]

#%%
# adm1 = du.R1D.model
# sto_h2 = adm1.stoichiometry.S_h2.to_numpy()
# sto_ch4 = adm1.stoichiometry.S_ch4.to_numpy()
# sto_ac = adm1.stoichiometry.S_ac.to_numpy()

# def r(unit, sto):
#     react = sum(unit._tempstate['rhos'] * sto[4:12])
#     transfer = sum(unit._tempstate['gas_transfer'] * sto[-3:])
#     return react, transfer, react+transfer

# r1h2 = []
# r1ac = []
# r1ch4 = []

# r2h2 = []
# r2ac = []
# r2ch4 = []
# for sys in systems:
#     R1, R2 = [u for u in sys.units if isinstance(u, su.AnaerobicCSTR)]
#     r1h2.append(r(R1, sto_h2))
#     r1ac.append(r(R1, sto_ac))
#     r1ch4.append(r(R1, sto_ch4))
#     r2h2.append(r(R2, sto_h2))
#     r2ac.append(r(R2, sto_ac))
#     r2ch4.append(r(R2, sto_ch4))
