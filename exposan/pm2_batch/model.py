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

def import_exp_data(kind=''):
    if kind == 'include':
        file = ospath.join(data_path, 'batch_exp_result_chli.xlsx')
        result_exp = load_data(file, sheet=0, index_col=None)

        t_exp_x = result_exp.t_stamp_x.to_numpy()
        t_exp_nh = result_exp.t_stamp_nh.to_numpy()
        t_exp_p = result_exp.t_stamp_p.to_numpy()

        t_exp_x=t_exp_x[~np.isnan(t_exp_x)]
        t_exp_p=t_exp_p[~np.isnan(t_exp_p)]

        vss_exp = result_exp.VSS.to_numpy()
        ch_exp = result_exp.CH.to_numpy()
        li_exp = result_exp.LI.to_numpy()
        snh_exp = result_exp.S_NH.to_numpy()
        sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

        vss_exp=vss_exp[~np.isnan(vss_exp)]
        ch_exp=ch_exp[~np.isnan(ch_exp)]
        li_exp=li_exp[~np.isnan(li_exp)]
        sp_exp=sp_exp[~np.isnan(sp_exp)]

        # vss_range = max(vss_exp) - min(vss_exp)
        # ch_range = max(ch_exp) - min(ch_exp)
        # li_range = max(li_exp) - min(li_exp)
        # snh_range = max(snh_exp) - min(snh_exp)
        # sp_range = max(sp_exp) - min(sp_exp)

        return (t_exp_x, vss_exp), (t_exp_nh, snh_exp), (t_exp_p, sp_exp), (t_exp_x, ch_exp), (t_exp_x, li_exp)

    else:
        file = ospath.join(data_path, 'batch_exp_result.xlsx')
        result_exp = load_data(file, sheet=0, index_col=None)

        t_exp = result_exp.t_stamp.to_numpy()

        vss_exp = result_exp.VSS.to_numpy()
        snh_exp = result_exp.S_NH.to_numpy()
        sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

        # vss_range = max(vss_exp) - min(vss_exp)
        # snh_range = max(snh_exp) - min(snh_exp)
        # sp_range = max(sp_exp) - min(sp_exp)

        return (t_exp, vss_exp), (t_exp, snh_exp), (t_exp, sp_exp), (), ()

#%%
#### 'kind=include': using metrics of VSS, CH, LI, SNH, SP with unit process experiment data (May 2022) #########
# file = ospath.join(data_path, 'batch_exp_result_chli.xlsx')
# result_exp = load_data(file, sheet=0, index_col=None)

# t_exp_x = result_exp.t_stamp_x.to_numpy()
# t_exp_nh = result_exp.t_stamp_nh.to_numpy()
# t_exp_p = result_exp.t_stamp_p.to_numpy()

# vss_exp = result_exp.VSS.to_numpy()
# ch_exp = result_exp.CH.to_numpy()
# li_exp = result_exp.LI.to_numpy()
# snh_exp = result_exp.S_NH.to_numpy()
# sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

# vss_range = max(vss_exp) - min(vss_exp)
# ch_range = max(ch_exp) - min(ch_exp)
# li_range = max(li_exp) - min(li_exp)
# snh_range = max(snh_exp) - min(snh_exp)
# sp_range = max(sp_exp) - min(sp_exp)

#%%
#### 'kind=exclude': using metrics of VSS, SNH, SP with kinetic assay data (May 2022) #########
# file = ospath.join(data_path, 'batch_exp_result.xlsx')
# result_exp = load_data(file, sheet=0, index_col=None)

# t_exp = result_exp.t_stamp.to_numpy()

# vss_exp = result_exp.VSS.to_numpy()
# snh_exp = result_exp.S_NH.to_numpy()
# sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

# vss_range = max(vss_exp) - min(vss_exp)
# snh_range = max(snh_exp) - min(snh_exp)
# sp_range = max(sp_exp) - min(sp_exp)

#%%
def create_model(system=None, kind=''):
    sys = create_system()
    # sys = system or create_system()

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
        D = get_uniform_w_frac(b, 0.50)
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
    idx_ch = cmps.indices(['X_CH'])
    idx_li = cmps.indices(['X_LI'])

    imass = cmps.i_mass[idx_vss]

    # import_exp_data
    exp_data = import_exp_data(kind)

    vss_pair, snh_pair, sp_pair, ch_pair, li_pair = exp_data

    t_exp_vss, vss_exp = vss_pair
    t_exp_nh, snh_exp = snh_pair
    t_exp_p, sp_exp = sp_pair

    vss_range = max(vss_exp) - min(vss_exp)
    snh_range = max(snh_exp) - min(snh_exp)
    sp_range = max(sp_exp) - min(sp_exp)

    @metric(element='Error')
    def NRMSE_VSS():

        t_simul = PBR.scope.time_series
        result_simul_vss = PBR.scope.record[:,idx_vss]

        vss = np.sum(result_simul_vss * imass, axis = 1)
        # vss = result_simul[:,0] * cmps.X_ALG.i_mass + result_simul[0:,1] * cmps.X_CH.i_mass + result_simul[0:,2] * cmps.X_LI.i_mass + result_simul[0:,3] + result_simul[0:,4]

        f = interp1d(t_simul, vss)
        vss_simul = f(t_exp_vss)

        rmse_vss = (sum((vss_exp - vss_simul)**2)/len(t_exp_vss))**0.5
        nrmse_vss = rmse_vss/vss_range
        return nrmse_vss

    @metric(element='Error')
    def NRMSE_SNH():
        t_simul = PBR.scope.time_series
        result_simul_snh = PBR.scope.record[:,idx_snh].flatten()

        f = interp1d(t_simul, result_simul_snh)
        snh_simul = f(t_exp_nh)

        rmse_snh = (sum((snh_exp - snh_simul)**2)/len(t_exp_nh))**0.5
        nrmse_snh = rmse_snh/snh_range
        return nrmse_snh

    @metric(element='Error')
    def NRMSE_SP():
        t_simul = PBR.scope.time_series
        result_simul_sp = PBR.scope.record[:,idx_sp].flatten()

        f = interp1d(t_simul, result_simul_sp)
        sp_simul = f(t_exp_p)

        rmse_sp = (sum((sp_exp - sp_simul)**2)/len(t_exp_p))**0.5
        nrmse_sp = rmse_sp/sp_range
        return nrmse_sp

    if kind == 'include':
        t_exp_ch, ch_exp = ch_pair
        t_exp_li, li_exp = li_pair

        ch_range = max(ch_exp) - min(ch_exp)
        li_range = max(li_exp) - min(li_exp)

        @metric(element='Error')
        def NRMSE_CH():
            t_simul = PBR.scope.time_series
            result_simul_ch = PBR.scope.record[:,idx_ch].flatten()

            f = interp1d(t_simul, result_simul_ch)
            ch_simul = f(t_exp_ch)
            rmse_ch = (sum((ch_exp - ch_simul)**2)/len(t_exp_ch))**0.5
            nrmse_ch = rmse_ch/ch_range
            return nrmse_ch

        @metric(element='Error')
        def NRMSE_LI():
            t_simul = PBR.scope.time_series
            result_simul_li = PBR.scope.record[:,idx_li].flatten()

            f = interp1d(t_simul, result_simul_li)
            li_simul = f(t_exp_li)
            rmse_li = (sum((li_exp - li_simul)**2)/len(t_exp_li))**0.5
            nrmse_li = rmse_li/li_range
            return nrmse_li


        # vss_pair, snh_pair, sp_pair = exp_data

        # t_exp, vss_exp = vss_pair
        # t_exp, snh_exp = snh_pair
        # t_exp, sp_exp = sp_pair

        # vss_range = max(vss_exp) - min(vss_exp)
        # snh_range = max(snh_exp) - min(snh_exp)
        # sp_range = max(sp_exp) - min(sp_exp)

        # @metric(element='Error')
        # def NRMSE_VSS_t():
        #     t_simul = PBR.scope.time_series
        #     result_simul_vss = PBR.scope.record[:,idx_vss]

        #     vss = np.sum(result_simul_vss * imass, axis = 1)
        #     # vss = result_simul[:,0] * cmps.X_ALG.i_mass + result_simul[0:,1] * cmps.X_CH.i_mass + result_simul[0:,2] * cmps.X_LI.i_mass + result_simul[0:,3] + result_simul[0:,4]

        #     f = interp1d(t_simul, vss)
        #     vss_simul = f(t_exp)

        #     rmse_vss = (sum((vss_exp - vss_simul)**2)/len(t_exp))**0.5
        #     nrmse_vss = rmse_vss/vss_range
        #     return nrmse_vss

        # @metric(element='Error')
        # def NRMSE_SNH_t():
        #     t_simul = PBR.scope.time_series
        #     result_simul_snh = PBR.scope.record[:,idx_snh].flatten()

        #     f = interp1d(t_simul, result_simul_snh)
        #     snh_simul = f(t_exp)

        #     rmse_snh = (sum((snh_exp - snh_simul)**2)/len(t_exp))**0.5
        #     nrmse_shh = rmse_snh/snh_range
        #     return nrmse_shh

        # @metric(element='Error')
        # def NRMSE_SP_t():
        #     t_simul = PBR.scope.time_series
        #     result_simul_sp = PBR.scope.record[:,idx_sp].flatten()

        #     f = interp1d(t_simul, result_simul_sp)
        #     sp_simul = f(t_exp)

        #     rmse_sp = (sum((sp_exp - sp_simul)**2)/len(t_exp))**0.5
        #     nrmse_sp = rmse_sp/sp_range
        #     return nrmse_sp

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
    # model = create_model()
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