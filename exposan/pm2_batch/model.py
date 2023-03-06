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

import os, numpy as np, qsdsan as qs
from scipy.interpolate import interp1d
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, ospath, time_printer, load_data

from exposan.pm2_batch import (
    create_system,
    data_path,
    results_path,
    )

__all__ = ('create_model', 'run_uncertainty',)

baseline_values = {
    'a_c': (0.049, 'm^2/g TSS'),
    'arr_a': (1.8e10, ''),
    'arr_e': (6842, 'K'),
    'beta_1': (2.9, ''),
    'beta_2': (3.5, ''),
    'b_reactor': (0.03, 'm'),
    'k_gamma': (1e-05, ''),
    'K_N': (0.1, 'g N/m^3'),
    'K_P': (1.0, 'g P/m^3'),
    'K_A': (6.3, 'g COD/m^3'),
    'K_F': (6.3, 'g COD/m^3'),
    'rho': (1.186, ''),
    'K_STO': (1.566, 'g COD/g COD'),
    'f_CH_max': (0.819, 'g COD/g COD'),
    'f_LI_max': (3.249, 'g COD/g COD'),
    'Q_N_max': (0.417, 'g N/g COD'),
    # 'Q_N_min': (0.082, 'g N/g COD'),   # Theoretically grounded
    'Q_P_max': (0.092, 'g P/g COD'),
    # 'Q_P_min': (0.0163, 'g P/g COD'),  # Theoretically grounded
    'exponent': (4, ''),
    'n_dark': (0.7, '')
    }

baseline_values_hifrac = {
    'I_n': (1000, 'uE/m^2/s'),              # Increased baseline
    'I_opt': (1500, 'uE/m^2/s'),            # Increased baseline
    'm_ATP': (15.835, 'g ATP/g COD/d'),
    'mu_max': (1.969, 'd^(-1)'),
    'q_CH': (0.594, 'g COD/g COD/d'),
    'q_LI': (0.91, 'g COD/g COD/d'),
    'V_NH': (0.254, 'g N/g COD/d'),
    'V_NO': (0.254, 'g N/g COD/d'),
    'V_P': (0.016, 'g P/g COD/d'),
    }

file = ospath.join(data_path, 'batch_exp_result.xlsx')
result_exp = load_data(file, sheet=0, index_col=None)

t_exp = result_exp.t_stamp.to_numpy()

vss_exp = result_exp.VSS.to_numpy()
snh_exp = result_exp.S_NH.to_numpy()
sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

vss_range = max(vss_exp) - min(vss_exp)
snh_range = max(snh_exp) - min(snh_exp)
sp_range = max(sp_exp) - min(sp_exp)

#%%

def create_model(system=None):
    sys = system or create_system()

    model = qs.Model(system=sys, exception_hook='raise')
    param = model.parameter
    metric = model.metric

    PBR = sys.units[0]
    pm2 = PBR.model

    cmps = PBR.components

    ##### General model with all uncertain variables #####

    # Add Uncertainty Parameters
    get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

    for k, v in baseline_values.items():
        b, units = v
        D = get_uniform_w_frac(b, 0.25)
        param(setter=DictAttrSetter(pm2.rate_function, '_params', k),
              name=k, element=PBR, kind='coupled', units=units, distribution=D)

    for k, v in baseline_values_hifrac.items():
        b, units = v
        D = get_uniform_w_frac(b, 0.50)
        param(setter=DictAttrSetter(pm2.rate_function, '_params', k),
              name=k, element=PBR, kind='coupled', units=units, distribution=D)

    ##### Add universal evaluation metrics #####

    idx_vss = cmps.indices(['X_ALG', 'X_CH', 'X_LI', 'X_N_ALG', 'X_P_ALG'])
    idx_snh = cmps.indices(['S_NH'])
    idx_sp = cmps.indices(['S_P'])

    imass = cmps.i_mass[idx_vss]

    # @metric (name='MAPE_VSS', units='%', element='')
    # def get_MAPE_VSS():
    @metric(element='Error')
    def NRMSE_VSS():
        t_simul = PBR.scope.time_series
        result_simul_vss = PBR.scope.record[:,idx_vss]

        vss = np.sum(result_simul_vss * imass, axis = 1)
        # vss = result_simul[:,0] * cmps.X_ALG.i_mass + result_simul[0:,1] * cmps.X_CH.i_mass + result_simul[0:,2] * cmps.X_LI.i_mass + result_simul[0:,3] + result_simul[0:,4]

        f = interp1d(t_simul, vss)
        vss_simul = f(t_exp)

        # mape_vss = sum(np.abs(vss_exp - vss_simul)/np.abs(vss_exp))/len(t_exp) * 100
        # return mape_vss
        rmse_vss = (sum((vss_exp - vss_simul)**2)/len(t_exp))**0.5
        return rmse_vss/vss_range


    # @metric (name='MAPE_SNH', units='%', element='')
    # def get_MAPE_SNH():
    @metric(element='Error')
    def NRMSE_SNH():
        t_simul = PBR.scope.time_series
        result_simul_snh = PBR.scope.record[:,idx_snh].flatten()

        f = interp1d(t_simul, result_simul_snh)
        snh_simul = f(t_exp)

        # mape_snh = sum(np.abs(snh_exp - snh_simul)/np.abs(snh_exp))/len(t_exp) * 100
        # return mape_snh
        rmse_snh = (sum((snh_exp - snh_simul)**2)/len(t_exp))**0.5
        return rmse_snh/snh_range


    # @metric (name='MAPE_SP', units='%', element='')
    # def get_MAPE_SP():
    @metric(element='Error')
    def NRMSE_SP():
        t_simul = PBR.scope.time_series
        result_simul_sp = PBR.scope.record[:,idx_sp].flatten()

        f = interp1d(t_simul, result_simul_sp)
        sp_simul = f(t_exp)
        # mape_sp = sum(np.abs(sp_exp - sp_simul)/np.abs(sp_exp))/len(t_exp) * 100
        # return mape_sp
        rmse_sp = (sum((sp_exp - sp_simul)**2)/len(t_exp))**0.5
        return rmse_sp/sp_range

    # @metric(name='Biomass', units='mg/L', element='')
    # def get_X_ALG():
    #     return PBR.state['X_ALG']

    # @metric(name='Carbohydrate', units='mg/L', element='')
    # def get_X_CH():
    #     return PBR.state['X_CH']

    # @metric(name='Lipid', units='mg/L', element='')
    # def get_X_LI():
    #     return PBR.state['X_LI']

    # @metric(name='Stored N', units='mg/L', element='')
    # def get_X_N_ALG():
    #     return PBR.state['X_N_ALG']

    # @metric(name='Stored P', units='mg/L', element='')
    # def get_X_P_ALG():
    #     return PBR.state['X_P_ALG']

    # @metric(name='Nitrogen', units='mg/L', element='')
    # def get_S_NH():
    #     return PBR.state['S_NH']

    # @metric(name='Phosphorus', units='mg/L', element='')
    # def get_S_P():
    #     return PBR.state['S_P']

    # @metric(name='VSS', units='mg/L', element='')
    # def get_VSS():
    #     return PBR.state['X_ALG'] * cmps.X_ALG.i_mass + \
    #            PBR.state['X_CH'] * cmps.X_CH.i_mass + \
    #            PBR.state['X_LI'] * cmps.X_LI.i_mass + \
    #            PBR.state['X_N_ALG'] + PBR.state['X_P_ALG']

    return model

#%%

# =============================================================================
# Run uncertainty analysis
# =============================================================================

@time_printer
def run_uncertainty(model, N, T, t_step, method='BDF',
                    metrics_path='', timeseries_path='',
                    rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    model = create_model()
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = metrics_path or ospath.join(results_path, f'table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'time_series_data_{seed}')
        os.mkdir(folder)
        tpath = ospath.join(folder, 'state.npy')
    model.evaluate(
        state_reset_hook='reset_cache',
        t_span=t_span,
        t_eval=t_eval,
        method=method,
        export_state_to=tpath
        )
    model.table.to_excel(mpath)
    return model

#%%
if __name__ == '__main__':
    seed = 120
    t = 7
    t_step = 0.01
    # t_step = 0.25/24
    n = 10
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    model = create_model()
    model = run_uncertainty(model, n, t, t_step, method=method, seed=seed)