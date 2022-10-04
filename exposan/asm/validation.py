# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan.utils import time_printer
from exposan.asm import data_path, results_path, create_system


# %%

# =============================================================================
# Universal settings
# =============================================================================

# Streams
Temp = 273.15+20    # temperature [K]
Q = 18446           # influent flowrate [m3/d]

# Tanks
V_an = 1000    # anoxic zone tank volume [m3]
V_aer = 1333    # aerated zone tank volume [m3]


# %%

# =============================================================================
# Functions to create the needed configuration
# =============================================================================

def create_asm_validation_system(process_model='ASM1', aerated=False):
    # Load process model (including associated components)
    pc_lower = process_model.lower()
    if pc_lower == 'asm1':
        # Process model
        asm_kwargs = dict(
            Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
            mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
            eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
            K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
            path=os.path.join(data_path, '_asm1.tsv'),
            )
        # Influent
        inf_kwargs = {
            'concentrations': {
                'S_S': 69.5,
                'X_BH': 28.17,
                'X_S': 202.32,
                'X_I': 51.2,
                'S_NH': 31.56,
                'S_I': 30,
                'S_ND': 6.95,
                'X_ND': 10.59,
                'S_ALK': 7*12,
                },
            'units': ('m3/d', 'mg/L'),
            }
        # Initial conditions in the CSTR
        init_conds = {
                'S_I': 30,
                'S_S': 5,
                'X_I': 1000,
                'X_S': 100,
                'X_BH': 500,
                'X_BA': 100,
                'X_P': 100,
                'S_O': 2,
                'S_NO': 20,
                'S_NH': 2,
                'S_ND': 1,
                'X_ND': 1,
                'S_ALK': 7*12,
            }
    elif pc_lower == 'asm2d':
        asm_kwargs = dict(
                iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
                iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
                iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
                f_SI=0.0, Y_H=0.625, f_XI_H=0.1,  # Y_H=0.67 is GPS-X default
                Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
                Y_A=0.24, f_XI_AUT=0.1,
                K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
                mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0, # mu_H=4.0, b_H=0.3 are GPS-X defaults
                K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
                q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
                b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
                K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
                K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
                mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
                k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5,
            )
        inf_kwargs = {
            'concentrations': { # Henze et al., Activated Sludge Models ASM1, ASM2, ASM2d and ASM3, P91
                'S_F': 30,
                'S_A': 20, # 'S_A':0 is GPS-X default
                'S_NH4': 16,
                'S_PO4': 3.6,
                'S_I': 30,
                'S_ALK': 5*12, # mmol/K to mg C/L
                'X_I': 25,
                'X_S': 125,
                'X_H': 30,
                },
            'units': ('m3/d', 'mg/L'),
            }
        init_conds = {
                'S_F': 2, # 'S_F':5, 'S_O2':2 are GPS-X defaults
                'S_A': 5,
                'S_I': 30,
                'S_NH4': 22, # 'S_NH4':2, 'S_N2':0, 'S_NO3':20 are GPS-X defaults
                'S_N2': 21,
                'S_NO3': 0.001,
                'S_PO4': 11, # 'S_PO4':5 is GPS-X defaults
                'X_I': 1800, # 'X_I':1000, 'X_S':100, 'X_H':500 are GPS-X defaults
                'X_S': 150,
                'X_H': 1900,
                'X_PAO': 250, # 'X_PAO':200, 'X_PP':100, 'X_PHA':100 are GPS-X defaults
                'X_PP': 70,
                'X_PHA': 7,
                'X_AUT': 125, # 'X_AUT':100 is GPS-X defaults
                'S_ALK': 6*12, # 'S_ALK':7*12 is GPS-X defaults
            }
    else:
        raise ValueError(f'`process_model` can only be "ASM1" or "ASM2d", not {process_model}.')

    sys = create_system(
        process_model=process_model,
        aerated=aerated,
        asm_kwargs=asm_kwargs,
        inf_kwargs=inf_kwargs,
        init_conds=init_conds,
        )
    return sys


# %%

@time_printer
def run(process_model, aerated, t, t_step, method, simulate=True,
        save_stoichiometry=False, save_states=False, **kwargs):
    global sys
    sys = create_asm_validation_system(process_model, aerated=aerated)
    suffix = 'aer' if aerated else 'an'
    folder = path = os.path.join(results_path, 'validation')
    if not os.path.isdir(folder): os.mkdir(path)
    if save_stoichiometry:
        asm = sys.flowsheet.unit.CSTR.suspended_growth_model
        path = os.path.join(folder, f'{process_model}_{suffix}_stoichiometry.csv')
        asm.stoichiometry.to_csv(path)
    export_state_to = os.path.join(folder, f'{process_model}_{suffix}.xlsx') if save_states else ''
    if simulate:
        sys.simulate(state_reset_hook='reset_cache',
                     t_span=(0, t),
                     t_eval=np.arange(0, t+t_step, t_step),
                     method=method,
                     export_state_to=export_state_to,
                     **kwargs)

if __name__ == '__main__':
    t = 10
    t_step = 0.1
    for process_model in ('ASM1', 'ASM2d'):
        for aerated in (False, True):
            suffix = 'Aerated' if aerated else 'Anoxic'
            msg = f'{process_model}-{suffix}'
            print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
            print(f'Time span 0-{t}d \n')
            # run(process_model, aerated, t, t_step, method='BDF', simulate=True)
            run(process_model, aerated, t, t_step, method='BDF', simulate=True, save_stoichiometry=True, save_states=True)