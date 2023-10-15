# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    # Saumitra Rai <raisaumitra9@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan import (
    processes as pc,
    sanunits as su,
    WasteStream,
    )
#!!! Need to verify system settings
from exposan.bsm1 import (
    default_asm1_kwargs,
    default_asm1_inf_kwargs,
    default_asm1_init_conds,
    Q, Q_ras, Q_was, Temp, V_an, V_ae, 
    )

from exposan.bsm2 import results_path
from exposan.adm import default_init_conds as default_adm1_init_conds

__all__ = ('create_system',)


# %%


# def create_components():
#      asm2d_cmps = pc.create_asm2d_cmps()
#      asm2d_cmps.X_S.f_BOD5_COD = 0.54
#      S_CO2 = qs.Component('S_CO2', search_ID='CO2', particle_size='Dissolved gas',
#                           degradability='Undegradable', organic=False)
#      cmps1 = qs.Components.load_default()
#      # CH4 = qs.Component('S_CH4', search_ID='CH4', particle_size='Dissolved gas',
#      #                    degradability='Readily', organic=True)
#      # H2 = qs.Component('S_H2', search_ID='H2', particle_size='Dissolved gas',
#      #                   degradability='Readily', organic=False)
#      S_CH4 = cmps1.S_CH4.copy('S_CH4')
#      S_H2 = cmps1.S_H2.copy('S_H2')
#      Ash = cmps1.X_Ig_ISS.copy('Ash')
#      cmps = qs.Components([*asm2d_cmps, S_CO2, S_CH4, S_H2, Ash])
#      cmps.compile()
#      return cmps

def create_system(flowsheet=None):
    flowsheet = flowsheet or qs.Flowsheet('bsm2')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    unit = flowsheet.unit
    
    # ASM1 components and process model
    cmps_asm1 = pc.create_asm1_cmps()
    thermo_asm1 = qs.get_thermo()
    DO_ID = 'S_O'
    asm1 = pc.ASM1(**default_asm1_kwargs)
    aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=240, DOsat=8.0, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=84, DOsat=8.0, V=V_ae)
    
    # Influent
    wastewater = WasteStream('wastewater', T=Temp)
    wastewater.set_flow_by_concentration(Q, **default_asm1_inf_kwargs)
    
    # Primary clarifier
    C1 = su.PrimaryClarifierBSM2(
        'C1',
        ins=(wastewater, 'thickener_recycle', 'storage_recycle'),
        outs=('C1_eff', 'C1_underflow'),
        isdynamic=True,
        )

    # Unit operations in BSM1
    A1 = su.CSTR('A1', ins=[C1-0, 'RAS'], V_max=V_an,
                 aeration=None, suspended_growth_model=asm1)
    
    A2 = su.CSTR('A2', A1-0, V_max=V_an,
                 aeration=None, suspended_growth_model=asm1)
    
    O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    O3 = su.CSTR('O3', O2-0, # [RWW, 'treated'], split=[0.6, 0.4],
                 V_max=V_ae, aeration=aer3,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    C2 = su.FlatBottomCircularClarifier(
        'C2', O3-0, ['effluent', 1-A1, 'WAS'],
        underflow=Q_ras, wastage=Q_was, surface_area=1500,
        height=4, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=474, v_max_practical=250,
        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
        )
    TC1 = su.Thickener('TC1', C2-2, outs=['thickened_sludge', 1-C1], 
                       thickener_perc=7, TSS_removal_perc=97.14)
    M1 = su.Mixer('M1', ins=(C1-1, TC1-0))
        
    # Switch to ADM1 components for the anaerobic digester
    cmps_adm1 = pc.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N    
    
    J1 = su.ASMtoADM('J1', upstream=M1-0, thermo=thermo_adm1, isdynamic=True, adm1_model=adm1) # WAS is C1.outs[2]
    AD1 = su.AnaerobicCSTR('AD1', ins=J1.outs[0], outs=('biogas', 'AD_eff'), isdynamic=True ,model=adm1,                                    
                           retain_cmps=[i for i in cmps_adm1.IDs if i.startswith('X_')])
    AD1.set_init_conc(**default_adm1_init_conds)
    # Switch back to ASM1 components
    J2 = su.ADMtoASM('J2', upstream=AD1-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    qs.set_thermo(thermo_asm1)
    
    C3 = su.Centrifuge(ID='C3', ins=J2-0, outs = ['digested_sludge', 'C3_eff'],
                       thickener_perc=27, TSS_removal_perc=96.29)

    #!!! Should be a stroage tank with HRT = 1
    T1 = su.HydraulicDelay('T1', C3-1, t_delay=1)
    T1-0-2-C1
    
    sys = flowsheet.create_system('bsm2_sys')
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    sys.maxiter = 5000
    sys.set_dynamic_tracker(unit.A1, unit.C1, J1, AD1, J2)
    
    return sys


# %%

if __name__ == '__main__':
    bsm2_sys = create_system()