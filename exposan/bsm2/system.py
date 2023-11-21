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

import os, numpy as np, qsdsan as qs
from qsdsan import (
    processes as pc,
    sanunits as su,
    WasteStream,
    )
from qsdsan.utils import time_printer, load_data, get_SRT
from exposan.bsm2 import figures_path, results_path
from exposan.adm import default_init_conds as default_adm1_init_conds

__all__ = ('create_system',)


# %%

Q = 20648.361 # influent flowrate [m3/d]
Q_intr = 3 * Q # activated sludge process internal recycle [m3/d]
Q_ras = Q # recycle sludge flowrate
Q_was = 300 # sludge wastage flowrate
Temp = 273.15+14.85808 # temperature [K]
V_an = 1500 # anoxic zone tank volume
V_ae = 3000 # aerated zone tank volume
biomass_IDs = ('X_BH', 'X_BA')

# Parameters for ASM No. 1 at 15 degC, Tables 2-3
default_asm_kwargs = dict(
    # Table 2 Stoichiometric parameters
    Y_H=0.67,
    Y_A=0.24,
    f_P=0.08,
    i_XB=0.08, #0.086;
    i_XP=0.06,
    fr_SS_COD=0.75, # X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS
    # Table 3 Kinetic parameters
    mu_H=4.0, #6.0;
    K_S=10.0, #20;
    K_O_H=0.2, # K_OH = 0.2;
    K_NO=0.5,
    b_H=0.3, #0.62;
    mu_A=0.5, #0.8;
    K_NH=1.0,
    K_O_A=0.4, # K_OA = 0.4;
    b_A=0.05, #0.2;
    eta_g=0.8, # ny_g
    k_a=0.05, #0.08;
    k_h=3.0,
    K_X=0.1, #0.03;
    eta_h=0.8, # ny_h #0.4;
    # path=os.path.join(data_path, '_asm1.tsv'),
    )

# O2 saturation concentration at 15 degC
SOSAT1 = 8

# Parameters for ADM No. 1, Table 5, same as the default values
default_adm_kwargs = dict(
    # Table 5 Stoichiometric parameters
    # f_sI_sc = 0.1 # 1-f_xI_xc-f_ch_xc-f_pr_xc-f_li_xc
    f_xI_xc=0.2,
    f_ch_xc=0.2,
    f_pr_xc=0.2,
    f_li_xc=0.3,
    N_xc=2.686e-3, # 0.0376/14
    N_I=4.286e-3, # 0.06/14
    N_aa=7e-3,
    #!!! No C_xc, C_SI, C_ch, C_pr, C_li, C_XI, C_su, C_aa
    f_fa_li=0.95,
    # No C_fa, f_h2_su
    f_bu_su=0.13,
    f_pro_su=0.27,
    f_ac_su=0.41,
    # No N_bac, C_bu, C_pro, C_ac, C_bac
    Y_su=0.1,
    # No f_h2_aa
    f_va_aa=0.23,
    f_bu_aa=0.26,
    f_pro_aa=0.05,
    f_ac_aa=0.4,
    # No C_va
    Y_aa=0.08,
    Y_fa=0.06,
    Y_c4=0.06,
    Y_pro=0.04,
    # No C_ch4
    Y_ac=0.05,
    Y_h2=0.06,
    
    #!!! PAUSED
    f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
    q_dis=0.5, q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
    k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
    K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6,
    b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
    KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4,
    pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
    T_base=298.15, pKa_base=[14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86],
    Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],
    kLa=200, K_H_base=[7.8e-4, 1.4e-3, 3.5e-2],
    K_H_dH=[-4180, -14240, -19410],

    )

 
                



# Constant influent
default_inf_kwargs = {
    'concentrations': {
        'S_I': 27.226191,
        'S_S': 58.176186,
        'X_I': 92.499001,
        'X_S': 363.94347,
        'X_BH': 50.683288,
        'S_NH': 23.859466,
        'S_ND': 5.651606,
        'X_ND': 16.129816,
        'S_ALK': 7*12,
        },
    'units': ('m3/d', 'mg/L'),
    }


def create_system(flowsheet=None):
    flowsheet = flowsheet or qs.Flowsheet('bsm2')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    unit = flowsheet.unit
    
    # ASM1 components and process model
    cmps_asm1 = pc.create_asm1_cmps()
    thermo_asm1 = qs.get_thermo()
    DO_ID = 'S_O'
    asm1 = pc.ASM1(**default_asm_kwargs)
    #!!! Not sure where KLa are from
    aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=240, DOsat=SOSAT1, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=84, DOsat=SOSAT1, V=V_ae)
    
    # Influent
    wastewater = WasteStream('wastewater', T=Temp)
    wastewater.set_flow_by_concentration(Q, **default_inf_kwargs)
    
    # Primary clarifier using the Otterpohl-Freund model
    #!!! Should set V, not HRT
    # Where are other parameters used?
    C1 = su.PrimaryClarifierBSM2(
        'C1',
        ins=(wastewater, 'thickener_recycle', 'reject_water'),
        outs=('C1_eff', 'C1_underflow'),
        isdynamic=True,
        f_corr=0.65,
        ratio_uf=0.007, # f_PS
        )

    # Unit operations in BSM1
    A1 = su.CSTR('A1', ins=[C1-0, 'RWW', 'RAS'], V_max=V_an,
                 aeration=None, suspended_growth_model=asm1)
    
    A2 = su.CSTR('A2', A1-0, V_max=V_an,
                 aeration=None, suspended_growth_model=asm1)
    
    O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                 DO_ID=DO_ID, suspended_growth_model=asm1)
    
    O3 = su.CSTR('O3', O2-0, [1-A1, 'treated'], split=[0.6, 0.4],
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
    TC1 = su.Thickener('TC1', C2-2, outs=['thickened_sludge', 1-C1],
                       thickening_perc=7, TSS_removal_perc=98)
    M1 = su.Mixer('M1', ins=(C1-1, TC1-0))
        
    # Switch to ADM1 components for the anaerobic digester
    cmps_adm1 = pc.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = pc.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N    
    
    J1 = su.ASMtoADM('J1', upstream=M1-0, thermo=thermo_adm1, isdynamic=True, adm1_model=adm1) # WAS is C1.outs[2]
    AD1 = su.AnaerobicCSTR('AD1', ins=J1.outs[0], outs=('biogas', 'AD_eff'), isdynamic=True,
                           # Tables 7-8
                           V_liq=3400, V_gas=300, T=308.15,
                           model=adm1,
                           retain_cmps=[i for i in cmps_adm1.IDs if i.startswith('X_')])
    AD1.set_init_conc(**default_adm1_init_conds)
    # Switch back to ASM1 components
    J2 = su.ADMtoASM('J2', upstream=AD1-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    qs.set_thermo(thermo_asm1)
    
    # Dewater
    C3 = su.Centrifuge(ID='C3', ins=J2-0, outs=['digested_sludge', 2-C1],
                       thickening_perc=28, TSS_removal_perc=96.29)

    #!!! Should have a storage tank with HRT = 1,
    # where the outs should have a bypass stream and an out stream.
    # Now equivalent to 100% bypass.
    # T1 = su.HydraulicDelay('T1', C3-1, outs='liquid_recycle', t_delay=0)
    # T1-0-2-C1
    
    sys = flowsheet.create_system('bsm2_sys')
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    sys.maxiter = 5000
    sys.set_dynamic_tracker(unit.A1, unit.C1, J1, AD1, J2)
    
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
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    srt = get_SRT(sys, biomass_IDs)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

if __name__ == '__main__':
    sys = create_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 1
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    
    sys.diagram()
    # sys.diagram(file=os.path.join(figures_path, 'bsm2_sys'), format='png')
    
    sys.converge()
    for u in sys.units:
        if not hasattr(u, '_state'): u._init_dynamic()
    y0, idx, nr = sys._load_state()
    run(sys, t, t_step, method=method)
