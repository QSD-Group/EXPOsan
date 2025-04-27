# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from scipy.interpolate import interp1d
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, ospath, time_printer, load_data

from exposan.pm2abaco2 import (
    create_system,
    data_path,
    results_path,
    )

__all__ = (
    'baseline_values',
    'sensitive_params',
    'import_exp_data',
    'create_model',
    'run_uncertainty',
    )

modified_pm2abaco2_kwargs = dict(
    a_c=0.049, I_n=1500, arr_a=1.8e10, arr_e=6842, beta_1=2.90,
    beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
    K_N=0.1, K_P=1.0, K_A=6.3, K_G=6.3, rho=1.186, K_STO=1.566,
    f_CH_max=0.819, f_LI_max=3.249, m_ATP=15.835,
    mu_max=1.969, q_CH=0.594, q_LI=0.910,
    Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
    V_NH=0.254, V_NO=0.254, V_P=0.016, exponent=4,
    Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
    Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
    Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
    Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
    Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,  # modified baseline for pm2 (I_n, I_opt)
    
    Y_NH_NIT=0.18, Y_NO_NIT=0.19, Y_NH_HET=9.09, Y_G=0.45, Y_O2_NIT=0.09, Y_O2_HET=2.78, f_BAC=0.1,
    mu_max_NIT=0.75, mu_max_HET=3.4, temp_min_NIT=265, temp_max_NIT=311, temp_opt_NIT=293,
    temp_min_HET=270, temp_max_HET=315, temp_opt_HET=298, ph_min_NIT=2, ph_max_NIT=13.4, ph_opt_NIT=9,
    ph_min_HET=6, ph_max_HET=12, ph_opt_HET=9, K_S_O2_NIT=1.08, K_I_O2_NIT=104.9, K_S_O2_HET=1.98,
    K_S_NH_NIT=1.0, K_S_NH_HET=0.5, K_S_G_HET=0.32, theta_NIT=274.1, theta_HET=274.07,  # original baseline for abaco2
    
    # f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
    # K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
    # mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F_H=4.0,
    # K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
    # mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,  # original baseline for asm2d
    path=None,
    )  

# Parameters used for UA & SA
baseline_values = {
    'a_c': (0.049, 'm^2/g TSS'),                # (0.0245, 0.0735)
    'I_n': (1500, 'uE/m^2/s'),                  # (750, 2250)            # Increased baseline, differ from default_pm2_kwargs
    'arr_a': (1.8e10, ''),                      # (0.9e10, 2.7e10)
    'arr_e': (6842, 'K'),                       # (3421, 10263)
    'beta_1': (2.9, ''),                        # (1.45, 4.35)
    'beta_2': (3.5, ''),                        # (1.75, 5.25)
    'b_reactor': (0.03, 'm'),                   # (0.015, 0.045)
    'I_opt': (2000, 'uE/m^2/s'),                # (1000, 3000)           # Increased baseline, differ from default_pm2_kwargs
    'k_gamma': (1e-05, ''),                     # (0.5e-05, 1.5e-05)
    'K_N': (0.1, 'g N/m^3'),                    # (0.05, 0.15)
    'K_P': (1.0, 'g P/m^3'),                    # (0.5, 1.5)
    'K_A': (6.3, 'g COD/m^3'),                  # (3.15, 9.45)
    'K_G': (6.3, 'g COD/m^3'),                  # (3.15, 9.45)
    'rho': (1.186, ''),                         # (0.593, 1.779)
    'K_STO': (1.566, 'g COD/g COD'),            # (0.783, 2.349)
    'f_CH_max': (0.819, 'g COD/g COD'),         # (0.4095, 1.2285)
    'f_LI_max': (3.249, 'g COD/g COD'),         # (1.6245, 4.8735)
    'm_ATP': (15.835, 'g ATP/g COD/d'),            
    'mu_max': (1.969, 'd^(-1)'),                # (0.9845, 2.9535)
    'q_CH': (0.594, 'g COD/g COD/d'),            
    'q_LI': (0.910, 'g COD/g COD/d'),               
    'Q_N_max': (0.417, 'g N/g COD'),            # (0.2085, 0.6255)
    'Q_P_max': (0.092, 'g P/g COD'),            # (0.046, 0.138)  
    'V_NH': (0.254, 'g N/g COD/d'),              
    'V_NO': (0.254, 'g N/g COD/d'),            
    'V_P': (0.016, 'g P/g COD/d'),              
    'exponent': (4, ''),                        # (2, 6)
    'n_dark': (0.7, ''),                        # (0.35, 1.05)
    
    'Y_NH_NIT': (0.18, 'g COD/g N'),
    'Y_NO_NIT': (0.19, 'g COD/g N'),
    'Y_NH_HET': (9.09, 'g COD/g N'),
    'Y_G': (0.45, 'g COD/g COD'),
    'Y_O2_NIT': (0.09, 'g COD/g O2'),
    'Y_O2_HET': (2.78, 'g COD/g O2'), 
    'f_BAC': (0.1,'g COD/g COD'),
    'mu_max_NIT': (0.75, 'd^(-1)'), 
    'mu_max_HET': (3.4, 'd^(-1)'), 
    'temp_min_NIT': (265, 'K'), 
    'temp_max_NIT': (311, 'K'),
    'temp_opt_NIT': (293, 'K'),
    'temp_min_HET': (270, 'K'), 
    'temp_max_HET': (315, 'K'), 
    'temp_opt_HET': (298, 'K'), 
    'ph_min_NIT': (2, 'unitless'),
    'ph_max_NIT': (13.4, 'unitless'),
    'ph_opt_NIT': (9, 'unitless'),
    'ph_min_HET': (6, 'unitless'), 
    'ph_max_HET': (12, 'unitless'), 
    'ph_opt_HET': (9, 'unitless'), 
    'K_S_O2_NIT': (1.08, 'g O2/m^3'), 
    'K_I_O2_NIT': (104.9, 'g O2/m^3'), 
    'K_S_O2_HET': (1.98, 'g O2/m^3'),
    'K_S_NH_NIT': (1.0, 'g N/m^3'), 
    'K_S_NH_HET': (0.5, 'g N/m^3'), 
    'K_S_G_HET': (0.32, 'g COD/m^3'), 
    'theta_NIT': (274.1, 'K'), 
    'theta_HET': (274.07, 'K'),
    }

# Parameters used for calibration - sensitive parameters (baseline_value, units, bounds)
sensitive_params = {
    'arr_e': (6842, 'K', (4000, 8000)),
    'K_P': (1.0, 'g P/m^3', (0.1, 10)),
    'K_A': (6.3, 'g COD/m^3', (1, 10)),
    'rho': (1.186, '', (1, 5)),
    'K_STO': (1.566, 'g COD/g COD', (0.1, 10)),
    'm_ATP': (15.835, 'g ATP/g COD/d', (1, 20)),            
    'mu_max': (1.969, 'd^(-1)', (0.1, 5)),                
    'q_LI': (0.910, 'g COD/g COD/d', (0.1, 20)),               
    'V_NH': (0.254, 'g N/g COD/d', (0.01, 1)),
    'V_NO': (0.254, 'g N/g COD/d', (0.0007, 0.01)),
    'V_P': (0.016, 'g P/g COD/d', (0.001, 0.1)),
    'Y_G': (0.45, 'g COD/g COD', (0.1, 1)),
    'f_BAC': (0.1,'g COD/g COD', (0.01, 0.5)),
    'mu_max_HET': (3.4, 'd^(-1)', (1, 5)), 
    'temp_min_NIT': (265, 'K', (260, 273)), 
    'temp_max_NIT': (311, 'K', (300, 320)),
    'temp_opt_NIT': (293, 'K', (273, 310)),
    'temp_max_HET': (315, 'K', (300, 320)), 
    'temp_opt_HET': (298, 'K', (273, 310)), 
    'ph_max_NIT': (13.4, 'unitless', (10, 14)),
    'ph_min_HET': (6, 'unitless', (4, 7)), 
    }           

#%%
def import_exp_data(kind=''):
    file = ospath.join(data_path, 'batch_exp_result_november.xlsx')
    result_exp = load_data(file, sheet=0, index_col=None)

    t_exp = result_exp.t_stamp.to_numpy()

    vss_exp = result_exp.VSS.to_numpy()
    snh_exp = result_exp.S_NH.to_numpy()
    sno_exp = result_exp.S_NO.to_numpy()
    sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

    return (t_exp, vss_exp), (t_exp, snh_exp), (t_exp, sno_exp), (t_exp, sp_exp)
    
    # if kind == 'include':
    #     file = ospath.join(data_path, 'batch_exp_result_chli.xlsx')
    #     result_exp = load_data(file, sheet=0, index_col=None)

    #     t_exp_x = result_exp.t_stamp_x.to_numpy()
    #     t_exp_nh = result_exp.t_stamp_nh.to_numpy()
    #     t_exp_p = result_exp.t_stamp_p.to_numpy()

    #     t_exp_x=t_exp_x[~np.isnan(t_exp_x)]
    #     t_exp_p=t_exp_p[~np.isnan(t_exp_p)]

    #     vss_exp = result_exp.VSS.to_numpy()
    #     ch_exp = result_exp.CH.to_numpy()
    #     li_exp = result_exp.LI.to_numpy()
    #     snh_exp = result_exp.S_NH.to_numpy()
    #     sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

    #     vss_exp=vss_exp[~np.isnan(vss_exp)]
    #     ch_exp=ch_exp[~np.isnan(ch_exp)]
    #     li_exp=li_exp[~np.isnan(li_exp)]
    #     sp_exp=sp_exp[~np.isnan(sp_exp)]

    #     return (t_exp_x, vss_exp), (t_exp_nh, snh_exp), (t_exp_p, sp_exp), (t_exp_x, ch_exp), (t_exp_x, li_exp)

    # else:
    #     file = ospath.join(data_path, 'batch_exp_result.xlsx')
    #     result_exp = load_data(file, sheet=0, index_col=None)

    #     t_exp = result_exp.t_stamp.to_numpy()

    #     vss_exp = result_exp.VSS.to_numpy()
    #     snh_exp = result_exp.S_NH.to_numpy()
    #     sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

    #     return (t_exp, vss_exp), (t_exp, snh_exp), (t_exp, sp_exp), (), ()

#%%
def create_model(system=None, kind='', analysis=''):
    sys = create_system(kind=kind, pm2abaco2_kwargs=modified_pm2abaco2_kwargs)
    # sys = system or create_system()

    model = qs.Model(system=sys, exception_hook='warn')
    param = model.parameter
    metric = model.metric

    PBR = sys.units[0]
    pm2abaco2 = PBR.model

    cmps = PBR.components

    ##### General model with all uncertain variables #####

    # Add uncertainty parameters
    if analysis == 'uasa':
        get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

        for k, v in baseline_values.items():
            b, units = v
            D = get_uniform_w_frac(b, 0.2)
            param(setter=DictAttrSetter(pm2abaco2.rate_function, '_params', k),
                  name=k, element=PBR, kind='coupled', units=units, distribution=D)

    # Add sensitive parameters for calibration
    else:
        for k, v in sensitive_params.items():
            b, units, bounds = v
            D = shape.Uniform(lower=bounds[0], upper=bounds[1])
            param(setter=DictAttrSetter(pm2abaco2.rate_function, '_params', k),
                  name=k, element=PBR, kind='coupled', units=units, distribution=D, bounds=bounds)   # distribution=D (optional)

    ##### Add universal evaluation metrics #####

    idx_vss = cmps.indices(['X_ALG', 'X_PG', 'X_TAG', 'X_N_ALG', 'X_P_ALG'])
    idx_snh = cmps.indices(['S_NH'])
    idx_sno = cmps.indices(['S_NO'])
    idx_sp = cmps.indices(['S_P'])
    # idx_pg = cmps.indices(['X_PG'])
    # idx_tag = cmps.indices(['X_TAG'])

    imass = cmps.i_mass[idx_vss]

    # Import_exp_data
    exp_data = import_exp_data(kind)

    vss_pair, snh_pair, sno_pair, sp_pair = exp_data

    t_exp_vss, vss_exp = vss_pair
    t_exp_nh, snh_exp = snh_pair
    t_exp_no, sno_exp = sno_pair
    t_exp_p, sp_exp = sp_pair

    vss_range = max(vss_exp) - min(vss_exp)
    snh_range = max(snh_exp) - min(snh_exp)
    sno_range = max(sno_exp) - min(sno_exp)
    sp_range = max(sp_exp) - min(sp_exp)

    @metric(element='Error')
    def NRMSE_VSS():

        t_simul = PBR.scope.time_series
        result_simul_vss = PBR.scope.record[:,idx_vss]

        vss = np.sum(result_simul_vss * imass, axis = 1)
        # vss = result_simul[:,0] * cmps.X_ALG.i_mass + result_simul[0:,1] * cmps.X_PG.i_mass + result_simul[0:,2] * cmps.X_TAG.i_mass + result_simul[0:,3] + result_simul[0:,4]

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
    def NRMSE_SNO():
        t_simul = PBR.scope.time_series
        result_simul_sno = PBR.scope.record[:,idx_sno].flatten()

        f = interp1d(t_simul, result_simul_sno)
        sno_simul = f(t_exp_no)

        rmse_sno = (sum((sno_exp - sno_simul)**2)/len(t_exp_no))**0.5
        nrmse_sno = rmse_sno/sno_range
        return nrmse_sno

    @metric(element='Error')
    def NRMSE_SP():
        t_simul = PBR.scope.time_series
        result_simul_sp = PBR.scope.record[:,idx_sp].flatten()

        f = interp1d(t_simul, result_simul_sp)
        sp_simul = f(t_exp_p)

        rmse_sp = (sum((sp_exp - sp_simul)**2)/len(t_exp_p))**0.5
        nrmse_sp = rmse_sp/sp_range
        return nrmse_sp

    # if kind == 'include':
    #     t_exp_ch, ch_exp = ch_pair
    #     t_exp_li, li_exp = li_pair

    #     ch_range = max(ch_exp) - min(ch_exp)
    #     li_range = max(li_exp) - min(li_exp)

    #     @metric(element='Error')
    #     def NRMSE_CH():
    #         t_simul = PBR.scope.time_series
    #         result_simul_ch = PBR.scope.record[:,idx_pg].flatten()

    #         f = interp1d(t_simul, result_simul_ch)
    #         ch_simul = f(t_exp_ch)
    #         rmse_ch = (sum((ch_exp - ch_simul)**2)/len(t_exp_ch))**0.5
    #         nrmse_ch = rmse_ch/ch_range
    #         return nrmse_ch

    #     @metric(element='Error')
    #     def NRMSE_LI():
    #         t_simul = PBR.scope.time_series
    #         result_simul_li = PBR.scope.record[:,idx_tag].flatten()

    #         f = interp1d(t_simul, result_simul_li)
    #         li_simul = f(t_exp_li)
    #         rmse_li = (sum((li_exp - li_simul)**2)/len(t_exp_li))**0.5
    #         nrmse_li = rmse_li/li_range
    #         return nrmse_li

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
    t = 9
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