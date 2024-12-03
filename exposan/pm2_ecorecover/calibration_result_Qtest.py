# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np
from qsdsan import sanunits as su
from qsdsan import processes as pc
from qsdsan import WasteStream, System
from qsdsan.utils import ospath, time_printer, load_data, get_SRT, \
    ExogenousDynamicVariable as EDV

from exposan.pm2_ecorecover import (
    data_path,
    results_path,
    )

def batch_init(path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')

    for k in [MIX, PBR1, PBR2, PBR3, PBR4, PBR5, PBR6, PBR7, PBR8, PBR9, PBR10,
              PBR11, PBR12, PBR13, PBR14, PBR15, PBR16, PBR17, PBR18, PBR19, PBR20, MEV, RET]:
        k.set_init_conc(**dct[k._ID[:3]])

#%%
# =============================================================================
# Phototrophic-Mixotrophic Process Model (PM2)
# =============================================================================

############# load components #############
cmps = pc.create_pm2_cmps()

############# create WasteStream objects #################
# Q = 449.06           # influent flowrate [m3/d]
Temp = 286.08          # temperature [K]

T, I = EDV.batch_init(os.path.join(data_path, 'exo_vars_dynamic_influent_cali.xlsx'), 'linear')

T_mix = T
I_mix = EDV('light_I_mix', function=lambda t: 0)

DYINF = WasteStream('Dynamic_influent', T=Temp)
PHO = WasteStream('To_PBR', T=Temp)

ME = WasteStream('To_membrane', T=Temp)
INT = WasteStream('Internal_recycle', T=Temp)

TE = WasteStream('Effluent', T=Temp)
RETEN = WasteStream('Retentate', T=Temp)

RE = WasteStream('To_return_tank', T=Temp)
CE = WasteStream('To_centrifuge', T=Temp)

CEN = WasteStream('Centrate', T=Temp)
ALG = WasteStream('Harvested_biomass', T=Temp)

RAA = WasteStream('Return_activated_algae', T=Temp)

#%%
############# load and tailor process models #############
V_mix = 66.23
V_pbr = 77.49
V_mem = 7.03
V_ret = 6.21

pm2 = pc.PM2(arr_e=6663.36141724313, K_P=6.06569854392092, f_CH_max=9.60813888591872, exponent=7.56541058257826, q_CH=1.92792246509906, q_LI=26.1535941900048, V_NH=0.150722549179019, V_P=0.540050768528713,
              a_c=0.049, I_n=1500, arr_a=1.8e10, beta_1=2.90,
              beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
              K_N=0.1, K_A=6.3, K_G=6.3, rho=1.186, K_STO=1.566,
              f_LI_max=3.249, m_ATP=10,
              mu_max=1.969, Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
              V_NO=0.003, n_dark=0.7,
              Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
              Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
              Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
              Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
              Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317)  # ecorecover_cali (optuna results) seed777

#%%
############# create unit operations #####################

# Dynamic influent
SE = su.DynamicInfluent('SE', outs=[DYINF],
                        data_file=ospath.join(data_path, 'dynamic_influent_cali_Qtest.tsv'))

MIX = su.CSTR('MIX', ins=[DYINF, RAA], outs=[PHO], V_max=V_mix,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_mix, I_mix))

PBR1 = su.CSTR('PBR1', ins=MIX-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR2 = su.CSTR('PBR2', ins=PBR1-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR3 = su.CSTR('PBR3', ins=PBR2-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR4 = su.CSTR('PBR4', ins=PBR3-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR5 = su.CSTR('PBR5', ins=PBR4-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR6 = su.CSTR('PBR6', ins=PBR5-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR7 = su.CSTR('PBR7', ins=PBR6-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR8 = su.CSTR('PBR8', ins=PBR7-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR9 = su.CSTR('PBR9', ins=PBR8-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR10 = su.CSTR('PBR10', ins=PBR9-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR11 = su.CSTR('PBR11', ins=PBR10-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR12 = su.CSTR('PBR12', ins=PBR11-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR13 = su.CSTR('PBR13', ins=PBR12-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR14 = su.CSTR('PBR14', ins=PBR13-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR15 = su.CSTR('PBR15', ins=PBR14-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR16 = su.CSTR('PBR16', ins=PBR15-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR17 = su.CSTR('PBR17', ins=PBR16-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR18 = su.CSTR('PBR18', ins=PBR17-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR19 = su.CSTR('PBR19', ins=PBR18-0, V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))
PBR20 = su.CSTR('PBR20', ins=PBR19-0, outs=[ME, INT], split=[0.52, 0.48], V_max=V_pbr/20,
              aeration=None, suspended_growth_model=pm2, exogenous_vars=(T, I))

MEM = su.Splitter('MEM', PBR20-0, outs=[TE, RETEN], split=0.39*(1-cmps.x))

MEV = su.CSTR('MEV', ins=MEM-1, V_max=V_mem, aeration=None, suspended_growth_model=None)

POST_MEM = su.Splitter('POST_MEM', MEV-0, outs=[RE, CE], split=0.97)                    # changed compared to previous ver.

CENT = su.Splitter('CENT', POST_MEM-1, outs=[CEN, ALG], split={'X_CHL':0.33,
                                                                'X_ALG':0.33,
                                                                'X_PG':0.33,
                                                                'X_TAG':0.33,
                                                                'S_CO2':0.975,
                                                                'S_A':0.975,
                                                                'S_G':0.975,
                                                                'S_O2':0.975,
                                                                'S_NH':0.975,
                                                                'S_NO':0.975,
                                                                'S_P':0.975,
                                                                'X_N_ALG':0.33,
                                                                'X_P_ALG':0.33,
                                                                'H2O':0.975})

RET = su.CSTR('RET', ins=[PBR20-1, POST_MEM-0, CENT-0], outs=[RAA], V_max=V_ret,
              aeration=None, suspended_growth_model=None)

#%%

_init_conds = {
        'X_CHL':2.53,
        'X_ALG':505.97,
        'X_PG':22.99,
        'X_TAG':93.78,
        'S_CO2':30.0,
        'S_A':5.0,
        'S_G':5.0,
        'S_O2':5.0,
        'S_NH':35.80,
        'S_NO':0.7,
        'S_P':0.36,
        'X_N_ALG':3.23,
        'X_P_ALG':0.19,
    }

batch_init(ospath.join(data_path, 'initial_conditions_pm2_dynamic_influent_cali.xlsx'), 'default')

#%%

eco = System('EcoRecovery',
             path=(SE, MIX,
                   PBR1, PBR2, PBR3, PBR4, PBR5, PBR6, PBR7, PBR8, PBR9, PBR10,
                   PBR11, PBR12, PBR13, PBR14, PBR15, PBR16, PBR17, PBR18, PBR19, PBR20,
                   MEM, MEV, POST_MEM, CENT, RET),
             recycle=(RAA,))

eco.set_dynamic_tracker(SE, MIX, PBR1, PBR20, RET, TE, CEN, ALG)
eco.set_tolerance(rmol=1e-6)
bio_IDs = ('X_ALG',)

__all__ = (
    'cmps', 'pm2', 'DYINF', 'PHO', 'ME', 'INT',
    'TE', 'RETEN', 'RE', 'CE', 'CEN', 'ALG', 'RAA', 'eco',
    *(i.ID for i in eco.units),
    '_init_conds'
    )

#%%

@time_printer
def run(t, t_step, method=None, print_t=False, **kwargs):
    if method:
        eco.simulate(state_reset_hook='reset_cache',
                      t_span=(0,t),
                      t_eval = np.arange(0, t+t_step, t_step),
                      method=method,
                      # rtol=1e-2,
                      # atol=1e-3,
                      export_state_to=ospath.join(results_path, f'sol_{t}d_{method}_calibration_result_111824_Qtest.xlsx'),
                      print_t=print_t,
                      **kwargs)
    else:
        eco.simulate(state_reset_hook='reset_cache',
                      solver='odeint',
                      t=np.arange(0, t+t_step/30, t_step/30),
                      # export_state_to=f'results/sol_{t}d_odeint.xlsx',
                      print_msg=True,
                      print_t=print_t,
                      **kwargs)

        eco.simulate()

    unit_IDs = [u.ID for u in eco.units if isinstance(u, su.CSTR)]
    srt = get_SRT(eco, bio_IDs, active_unit_IDs=unit_IDs)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

if __name__ == '__main__':
    t = 25
    t_step = 0.1
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    # method = None
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    run(t, t_step, method=method,
        print_t = True,
        )

# ALG.get_mass_concentration('mg/l','X_ALG')
# RET.scope.plot_time_series(('X_ALG'))

# simulation ends with "Estimated SRT assuming at steady state is 6.98 days"