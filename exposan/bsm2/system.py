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
    default_init_conds,
    )

from exposan.bsm2 import figures_path, results_path
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

#!!! Not sure if this is already provided as the steady state
# Qin0 = 20648;
Q = 20648.361 # influent flowrate [m3/d]
Q_intr = 3 * Q #!!! what is this?
Q_ras = Q # recycle sludge flowrate
Q_was = 300 # sludge wastage flowrate
Temp = 273.15+14.85808 # temperature [K]
V_an = 1500 # anoxic zone tank volume
V_ae = 3000 # aerated zone tank volume

# Parameters for AS system at 15 degC, based on BSM1
default_asm_kwargs = dict(
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
    Y_H=0.67,
    Y_A=0.24,
    f_P=0.08,
    i_XB=0.08, #0.086;
    i_XP=0.06,
    fr_SS_COD=0.75, # X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS
    # path=os.path.join(data_path, '_asm1.tsv'),
    )


# O2 saturation concentration at 15 degC, based on BSM1
SOSAT1 = 8;

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

# default_adm_kwargs = dict(
#     Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#     mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
#     eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
#     K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
#     )

# def __new__(cls, components=None, path=None, N_xc=2.686e-3, N_I=4.286e-3, N_aa=7e-3,
#             f_ch_xc=0.2, f_pr_xc=0.2, f_li_xc=0.3, f_xI_xc=0.2,
#             f_fa_li=0.95, f_bu_su=0.13, f_pro_su=0.27, f_ac_su=0.41,
#             f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.4,
#             f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
#             Y_su=0.1, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06,
#             q_dis=0.5, q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
#             k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
#             K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6,
#             b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
#             KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4,
#             pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
#             T_base=298.15, pKa_base=[14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86],
#             Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],
#             kLa=200, K_H_base=[7.8e-4, 1.4e-3, 3.5e-2],
#             K_H_dH=[-4180, -14240, -19410],
#             **kwargs):

def create_system(flowsheet=None):
    flowsheet = flowsheet or qs.Flowsheet('bsm2')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    unit = flowsheet.unit
    
    # ASM1 components and process model
    cmps_asm1 = pc.create_asm1_cmps()
    thermo_asm1 = qs.get_thermo()
    DO_ID = 'S_O'
    asm1 = pc.ASM1(**default_asm_kwargs['asm1'])
    # The O2 saturation concentration of 8 is at 15 degC, based on BSM1
    aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=240, DOsat=8.0, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=84, DOsat=8.0, V=V_ae)
    
    # Influent
    wastewater = WasteStream('wastewater', T=Temp)
    wastewater.set_flow_by_concentration(Q, **default_inf_kwargs)
    
    # Primary clarifier using the Otterpohl-Freund model
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
    
    # 10-layer one-dimensional settler model
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
    
    C3 = su.Centrifuge(ID='C3', ins=J2-0, outs=['digested_sludge', 'C3_eff'],
                       thickener_perc=27, TSS_removal_perc=96.29)

    #!!! Should be a stroage tank with HRT = 1
    T1 = su.HydraulicDelay('T1', C3-1, outs='liquid_recycle', t_delay=1)
    T1-0-2-C1
    
    sys = flowsheet.create_system('bsm2_sys')
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    sys.maxiter = 5000
    sys.set_dynamic_tracker(unit.A1, unit.C1, J1, AD1, J2)
    
    return sys


# %%

if __name__ == '__main__':
    bsm2_sys = create_system()
    bsm2_sys.diagram(file=os.path.join(figures_path, 'bsm2_sys'), format='png')