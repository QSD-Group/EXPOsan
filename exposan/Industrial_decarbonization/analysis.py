# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 11:47:30 2023

@author: joy_c
"""
from qsdsan.utils import ospath, time_printer, get_SRT, \
    get_oxygen_heterotrophs, get_oxygen_autotrophs, get_airflow, get_P_blower
import qsdsan.sanunits as su
from qsdsan import WasteStream
from metro import data_path
from metro.systems import create_system, batch_init, Q_brewery, biomass_IDs
from metro.utils import industrial_wws_adm1, adm1_2_asm2d

#%%

@time_printer
def evaluate(system, ind_concs, t=30, method='RK45', inf=None):
    fs = system.flowsheet.stream
    fu = system.flowsheet.unit
    fs.industrial_wastewater.set_flow_by_concentration(
        Q_brewery, 
        concentrations=ind_concs, 
        units=('m3/d', 'mg/L'))
    
    system.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        )
    
    act_units = [u.ID for u in system.units if isinstance(u, (su.CSTR, su.FlatBottomCircularClarifier))]

    srt = get_SRT(system, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)    
    
    if inf is None: inf = WasteStream('mix')
    inf.mix_from(fu.A2.ins[1:])
    eff_sCOD = fs.effluent.composite('COD', particle_size='s')
    o2_hetero = get_oxygen_heterotrophs(inf.F_vol*24, inf.COD, eff_sCOD, 
                                        SRT=srt, f_d=0.1, Y_H=0.625, b_H=0.4)
    
    inf_TKN = inf.TN - inf.composite('N', subgroup=('S_NO3',))
    o2_auto = get_oxygen_autotrophs(inf.F_vol*24, inf.COD, eff_sCOD, 
                                    inf_TKN, SRT=srt, Y_H=0.625, U_AUT=1.0, 
                                    b_H=0.4, b_AUT=0.15, f_d=0.1, SF_DO=1.25)
    
    q_air = get_airflow(o2_hetero, o2_auto, oxygen_transfer_efficiency=10) # in m3/min
    P_blower = get_P_blower(q_air=q_air) # kW

    E_aeration = P_blower / sum([s.F_vol for s in system.feeds]) # kWh/m3
    sludge_prod = fs.sludge_DU1.composite('solids', True, particle_size='x', unit='ton/d')
    
    return {
        'o2_hetero':o2_hetero,  # kg/d
        'o2_auto':o2_auto,
        'Q_air':q_air*60,       # m3/hr
        'P_blower': P_blower,
        'E_aeration': E_aeration,
        'sludge_production': sludge_prod # ton TSS/day
        }

#%%
if __name__ == '__main__':
    sys = create_system()
    inf = WasteStream('mix')
    path = ospath.join(data_path, "initial_conditions_ASM2d.xlsx")
    batch_init(sys, path, 
               sheet='ss')
    
    industrial_wws_asm2d = {}
    metrics = {}
    for k, vals in industrial_wws_adm1.items():
        print(f'{k}\n{len(k)*"="}\n')
        industrial_wws_asm2d[k] = concs = adm1_2_asm2d(vals)
        metrics[k] = evaluate(sys, concs,
                              t=15,
                              # t=30,
                              method='RK23',
                              # method = 'RK45',
                              # method = 'BDF',
                              inf=inf
                              )


