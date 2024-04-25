# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

Reference:
    Alex, J.; Benedetti, L.; Copp, J. B.; Gernaey, K. V.; Jeppsson, U.;
    Nopens, I.; Pons, M. N.; Rosen, C.; Steyer, J. P.; Vanrolleghem, P. A.
    Benchmark Simulation Model No. 2 (BSM2).
    http://iwa-mia.org/wp-content/uploads/2022/09/TR3_BSM_TG_Tech_Report_no_3_BSM2_General_Description.pdf.

'''

import os, numpy as np, pandas as pd, qsdsan as qs
from qsdsan import (
    processes as pc,
    sanunits as su,
    WasteStream,
    )
from qsdsan.utils import time_printer, load_data, get_SRT
from exposan.bsm2 import data_path, figures_path, results_path
from exposan.bsm1 import data_path as bsm1_path
# from exposan.adm import default_init_conds as default_adm1_init_conds

__all__ = ('create_system',)


# %%

Q = 20648.361 # influent flowrate [m3/d]
Q_intr = 3 * Q # activated sludge process internal recycle [m3/d]
Q_ras = Q # recycle sludge flowrate
Q_was = 300 # sludge wastage flowrate
Temp = 273.15+14.85808 # temperature [K]
V_an = 1500 # anoxic zone tank volume
V_ae = 3000 # aerated zone tank volume

# O2 saturation concentration at 15 degC
SOSAT1 = 8

# Default initial conditions
asm1init = pd.read_csv(os.path.join(data_path, 'asm1init.csv'), index_col=0).to_dict('index')
settler1dinit = pd.read_csv(os.path.join(data_path, 'settler1dinit.csv'), index_col=0).to_dict('index')
adm1init = pd.read_csv(os.path.join(data_path, 'adm1init.csv'), index_col=0).to_dict('index')

default_inf_kwargs = {
    'concentrations': asm1init['inf'],
    'units': ('m3/d', 'mg/L'),
    }

def create_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('bsm2')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    unit = flowsheet.unit
    
    # ASM1 components and process model
    cmps_asm1 = pc.create_asm1_cmps()
    asm1 = pc.ASM1(components=cmps_asm1, 
                   path=os.path.join(bsm1_path, '_asm1.tsv'))
    thermo_asm1 = qs.get_thermo()
    DO_ID = 'S_O'
    #!!! Not sure where KLa are from
    # the BSM STR Table 3.1 suggests kLa=120 for tanks 3&4, kLa=60 for tank 5
    aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=120, DOsat=SOSAT1, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=60, DOsat=SOSAT1, V=V_ae)
    
    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(Q, **default_inf_kwargs)
    carb = WasteStream('carbon', T=Temp)
    carb.set_flow_by_concentration(2, {'S_S':400}, units=('m3/d', 'kg/m3'))
    
    # Primary clarifier using the Otterpohl-Freund model
    # Where are other parameters used?
    C1 = su.PrimaryClarifierBSM2(
        'C1',
        # ins=(inf, 'thickener_recycle', 'reject_water'),
        ins=(inf, 'reject'),
        outs=('C1_eff', 'C1_underflow'),
        isdynamic=True, 
        volume=900,
        # HRT=1/24, #!!! should set V (900 m3), not HRT
        f_corr=0.65,
        ratio_uf=0.007, # f_PS
        )

    # Unit operations in BSM1
    A1 = su.CSTR('A1', ins=[C1-0, 'RWW', 'RAS', carb], V_max=V_an,
    # A1 = su.CSTR('A1', ins=[C1-0, 'RWW', 'RAS'], V_max=V_an,
                 aeration=None, suspended_growth_model=asm1)
    
    A2 = su.CSTR('A2', A1-0, V_max=V_an,
                 aeration=None, suspended_growth_model=asm1)

    O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    O3 = su.CSTR('O3', O2-0, [1-A1, 'treated'], 
                 split=[0.6, 0.4],
                 V_max=V_ae, aeration=aer3,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    # 10-layer one-dimensional settler model, Table 4
    C2 = su.FlatBottomCircularClarifier(
        'C2', O3-1, ['effluent', 2-A1, 'WAS'],
        underflow=Q_ras, wastage=Q_was,
        # Table 4 and the corresponding section
        surface_area=1500, height=4, N_layer=10,
        feed_layer=5, # from top to bottom, 6th if from bottom to top
        X_threshold=3000, v_max=474, v_max_practical=250,
        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
        )
    
    # TC1 = su.Thickener('TC1', C2-2, outs=['thickened_sludge', 1-C1],
    TC1 = su.Thickener('TC1', C2-2, outs=['thickened_sludge', ''],
                       thickening_perc=7, TSS_removal_perc=98)
    M1 = su.Mixer('M1', ins=(C1-1, TC1-0))
        
    # Switch to ADM1 components for the anaerobic digester
    cmps_adm1 = pc.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = pc.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N # slight difference
    
    J1 = su.ASMtoADM('J1', upstream=M1-0, thermo=thermo_adm1, isdynamic=True, adm1_model=adm1) # WAS is C1.outs[2]
    AD1 = su.AnaerobicCSTR('AD1', ins=J1.outs[0], outs=('biogas', 'AD_eff'), isdynamic=True,
                           # Tables 7-8
                           V_liq=3400, V_gas=300, T=308.15,
                           model=adm1,)
    # Switch back to ASM1 components
    J2 = su.ADMtoASM('J2', upstream=AD1-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    J2.bio_to_xs = 0.79
    qs.set_thermo(thermo_asm1)
    
    # Dewatering
    # C3 = su.Centrifuge(ID='C3', ins=J2-0, outs=['digested_sludge', 2-C1],
    C3 = su.Centrifuge(ID='C3', ins=J2-0, outs=['digested_sludge', ''],
                       thickening_perc=28, TSS_removal_perc=96.29)

    M2 = su.Mixer('M2', ins=(TC1-1, C3-1), outs=1-C1)
    
    #!!! Should have a storage tank with HRT = 1,
    # where the outs should have a bypass stream and an out stream.
    # Now equivalent to 100% bypass.
    # T1 = su.HydraulicDelay('T1', C3-1, outs='liquid_recycle', t_delay=0)
    # T1-0-2-C1
    
    if default_init_conds:
        for i in ('C1', 'A1', 'A2', 'O1', 'O2', 'O3'):
            getattr(unit, i).set_init_conc(**asm1init[i])
        C2.set_init_TSS(list(settler1dinit['C2'].values()))
        AD1.set_init_conc(**adm1init['AD1'])

    
    # sys = flowsheet.create_system('bsm2_sys')
    sys = qs.System('bsm2_sys', 
                    path=(C1, A1, A2, O1, O2, O3, C2, 
                    #         TC1, M1, J1, AD1, J2, C3),
                    # recycle=(O3-0, C2-1, TC1-1, C3-1, )
                          TC1, M1, J1, AD1, J2, C3, M2),
                    recycle=(O3-0, C2-1, M2-0)
                    )
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    sys.maxiter = 5000
    sys.set_dynamic_tracker(C1, A1, O3, C2, J1, AD1, J2, C3)
    
    return sys


# %%

@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        # t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    
    # biomass_IDs = ('X_BH', 'X_BA')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[sys.flowsheet.stream.digested_sludge],
    #               active_unit_IDs=('C3'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

#%%
if __name__ == '__main__':
    sys = create_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    # cmps_adm1 = J1.components
    # cmps_asm1 = J2.components
    
    t = 20
    t_step = 1
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    # sys._setup()
    # sys.converge()
    # sys.diagram()
    # sys.diagram(file=os.path.join(figures_path, 'bsm2_sys'), format='png')
    
    # sys.converge()
    # for u in sys.units:
    #     if not hasattr(u, '_state'): u._init_dynamic()
    # y0, idx, nr = sys._load_state()
    

    # J1.state # post ASM2ADM interface
    # J2.state # post ADM2ASM interface