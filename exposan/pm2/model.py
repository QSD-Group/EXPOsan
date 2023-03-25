# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from scipy.interpolate import interp1d
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, ospath, load_data

# from qsdsan.utils import DictAttrSetter, AttrGetter, FuncGetter, \
#     get_SRT as srt, ospath, time_printer

from exposan.pm2 import (
    create_system,
    data_path,
    )

__all__ = (
    'sensitive_params',
    'import_scada_data',
    'create_model',
    )

sensitive_params = {
    'arr_e': (6842, 'K', (1000, 10000)),
    'K_P': (1.0, 'g P/m^3', (0.01, 100)),
    'f_CH_max': (0.819, 'g COD/g COD', (0.1, 10)),
    'exponent': (4, '', (1, 10)),
    'q_CH': (1, 'g COD/g COD/d', (0.1, 10)),
    'q_LI': (15, 'g COD/g COD/d', (1.5, 50)),
    'V_NH': (0.1, 'g N/g COD/d', (0.01, 1)),
    'V_P': (0.2, 'g P/g COD/d', (0.01, 1))
    }

# Parameters used for UA & SA
# baseline_values = {
#     'a_c': (0.049, 'm^2/g TSS'),                # (0.0245, 0.0735)
#     'arr_a': (1.8e10, ''),                      # (0.9e10, 2.7e10)    v
#     'arr_e': (6842, 'K'),                       # (3421, 10263)
#     'beta_1': (2.9, ''),                        # (1.45, 4.35)
#     'beta_2': (3.5, ''),                        # (1.75, 5.25)
#     'b_reactor': (0.03, 'm'),                   # (0.015, 0.045)
#     'k_gamma': (1e-05, ''),                     # (0.5e-05, 1.5e-05)
#     'K_N': (0.1, 'g N/m^3'),                    # (0.05, 0.15)
#     'K_P': (1.0, 'g P/m^3'),                    # (0.5, 1.5)
#     'K_A': (6.3, 'g COD/m^3'),                  # (3.15, 9.45)
#     'K_F': (6.3, 'g COD/m^3'),                  # (3.15, 9.45)
#     'rho': (1.186, ''),                         # (0.593, 1.779)
#     'K_STO': (1.566, 'g COD/g COD'),            # (0.783, 2.349)
#     'f_CH_max': (0.819, 'g COD/g COD'),         # (0.4095, 1.2285)
#     'f_LI_max': (3.249, 'g COD/g COD'),         # (1.6245, 4.8735)
#     'Q_N_max': (0.417, 'g N/g COD'),            # (0.2085, 0.6255)
#     'Q_P_max': (0.092, 'g P/g COD'),            # (0.046, 0.138)
#     'exponent': (4, ''),                        # (2, 6)
#     'n_dark': (0.7, ''),                        # (0.35, 1.05)
#     'I_n': (1500, 'uE/m^2/s'),                  # (750, 2250)            # Increased baseline, differ from default_pm2_kwargs
#     'I_opt': (2000, 'uE/m^2/s'),                # (1000, 3000)           # Increased baseline, differ from default_pm2_kwargs
#     'm_ATP': (10, 'g ATP/g COD/d'),             # (5, 15) v
#     'mu_max': (1.969, 'd^(-1)'),                # (0.9845, 2.9535)
#     'q_CH': (1, 'g COD/g COD/d'),               # (0.5, 1.5)          v
#     'q_LI': (15, 'g COD/g COD/d'),              # (7.5, 22.5)         v
#     'V_NH': (0.1, 'g N/g COD/d'),               # (0.05, 0.15)
#     'V_NO': (0.003, 'g N/g COD/d'),             # (0.0015, 0.0045)    v
#     'V_P': (0.2, 'g P/g COD/d')                 # (0.1, 0.3)          v
#     }

# pm2 = pc.PM2(mu_max=2, I_n=1000, I_opt=2000, V_NH=0.12, V_NO=0.0035, V_P=0.24, q_CH=1.5, q_LI=24,
#              arr_a=4*10**10, f_CH_max=1.2, f_LI_max=3.9)    # GY gutfeeling optimized - batch

# pm2 = pc.PM2(mu_max=1, I_n=1200, I_opt=2000, V_NH=0.09, V_NO=0.003, V_P=0.17, q_CH=0.6, q_LI=13,
#              Q_N_max=0.6, arr_e=6650, m_ATP=5)    # GY gutfeeling optimized - conti


# sensitive_params = {
#     'arr_e': (6842, 'K', (1000, 10000)),
#     'K_P': (1.0, 'g P/m^3', (0.0001, 10)),
#     'K_STO': (1.566, 'g COD/g COD', (0.0001, 20)),
#     'exponent': (4, '', (1, 10)),
#     'm_ATP': (15.835, 'g ATP/g COD/d', (1, 50)),
#     'q_CH': (0.594, 'g COD/g COD/d', (0.001, 10)),
#     'V_NH': (0.254, 'g N/g COD/d', (0.0001, 5)),
#     'V_P': (0.016, 'g P/g COD/d', (0.0001, 5))
#     }



#%%
def import_scada_data():

    file = ospath.join(data_path, 'conti_scada_result.xlsx')
    result_scada = load_data(file, sheet=0, index_col=None)

    t_scada = result_scada.t_stamp.to_numpy()

    vss_scada = result_scada.VSS.to_numpy()
    snh_scada = result_scada.S_NH.to_numpy()
    sp_scada = result_scada.S_P.to_numpy()
    # sp_scada = np.maximum(result_scada.S_P.to_numpy(), 0.015)

    return t_scada, vss_scada, snh_scada, sp_scada

#%%

def create_model(system=None, flowsheet=None):
    # sys = create_system(flowsheet)
    sys = system or create_system(flowsheet)

    stream = sys.flowsheet.stream
    DYINF, PHO, ME, INT, TE, RETEN, RE, CE, CEN, ALG, RAA = stream.Dynamic_influent, stream.To_PBR, stream.To_membrane, stream.Internal_recycle, stream.Effluent, stream.Retentate, stream.To_return_tank, stream.To_centrifuge, stream.Centrate, stream.Harvested_biomass, stream.Return_activated_algae

    unit = sys.flowsheet.unit
    SE, MIX, PBR1, PBR2, PBR3, PBR4, PBR5, PBR6, PBR7, PBR8, PBR9, PBR10, PBR11, PBR12, PBR13, PBR14, PBR15, PBR16, PBR17, PBR18, PBR19, PBR20, MEM, MEV, POST_MEM, CENT, RET = unit.SE, unit.MIX, unit.PBR1, unit.PBR2, unit.PBR3, unit.PBR4, unit.PBR5, unit.PBR6, unit.PBR7, unit.PBR8, unit.PBR9, unit.PBR10, unit.PBR11, unit.PBR12, unit.PBR13, unit.PBR14, unit.PBR15, unit.PBR16, unit.PBR17, unit.PBR18, unit.PBR19, unit.PBR20, unit.MEM, unit.MEV, unit.POST_MEM, unit.CENT, unit.RET

    model = qs.Model(system=sys, exception_hook='warn')
    param = model.parameter
    metric = model.metric

    pm2 = MIX.suspended_growth_model

    cmps = MIX.components

    ##### General model with all uncertain variables #####

    # Add Uncertainty Parameters

    for k, v in sensitive_params.items():
        b, units, bounds = v
        D = shape.Uniform(lower=bounds[0], upper=bounds[1])
        param(setter=DictAttrSetter(pm2.rate_function, '_params', k),
              name=k, element=RET, kind='coupled', units=units, distribution=D, bounds=bounds)   # distribution=D (optional)

    ##### Add universal evaluation metrics #####

    idx_vss = cmps.indices(['X_ALG', 'X_CH', 'X_LI', 'X_N_ALG', 'X_P_ALG'])
    idx_snh = cmps.indices(['S_NH'])
    idx_sp = cmps.indices(['S_P'])

    imass = cmps.i_mass[idx_vss]

    # import_scada_data
    scada_data = import_scada_data()

    t_scada, vss_scada, snh_scada, sp_scada = scada_data

    vss_range = max(vss_scada) - min(vss_scada)
    snh_range = max(snh_scada) - min(snh_scada)
    sp_range = max(sp_scada) - min(sp_scada)

    @metric(element='Error')
    def NRMSE_VSS():

        t_simul = RET.scope.time_series
        result_simul_vss = RET.scope.record[:,idx_vss]

        vss = np.sum(result_simul_vss * imass, axis = 1)
        # vss = result_simul[:,0] * cmps.X_ALG.i_mass + result_simul[0:,1] * cmps.X_CH.i_mass + result_simul[0:,2] * cmps.X_LI.i_mass + result_simul[0:,3] + result_simul[0:,4]

        f = interp1d(t_simul, vss)
        vss_simul = f(t_scada)

        rmse_vss = (sum((vss_scada - vss_simul)**2)/len(t_scada))**0.5
        nrmse_vss = rmse_vss/vss_range
        return nrmse_vss

    @metric(element='Error')
    def NRMSE_SNH():
        t_simul = TE.scope.time_series
        result_simul_snh = TE.scope.record[:,idx_snh].flatten()

        f = interp1d(t_simul, result_simul_snh)
        snh_simul = f(t_scada)

        rmse_snh = (sum((snh_scada - snh_simul)**2)/len(t_scada))**0.5
        nrmse_snh = rmse_snh/snh_range
        return nrmse_snh

    @metric(element='Error')
    def NRMSE_SP():
        t_simul = TE.scope.time_series
        result_simul_sp = TE.scope.record[:,idx_sp].flatten()

        f = interp1d(t_simul, result_simul_sp)
        sp_simul = f(t_scada)

        rmse_sp = (sum((sp_scada - sp_simul)**2)/len(t_scada))**0.5
        nrmse_sp = rmse_sp/sp_range
        return nrmse_sp

    return model