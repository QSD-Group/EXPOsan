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
from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, AttrGetter, FuncGetter, \
    get_SRT as srt, ospath, time_printer
      
from exposan.pm2 import (
    biomass_IDs,
    create_system, 
    default_init_conds as _ic, 
    results_path,
    T, I,
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
    'I_n': (1000, 'uE/m^2/s'),  # Increased baseline
    'I_opt': (1500, 'uE/m^2/s'),   # Increased baseline
    'm_ATP': (15.835, 'g ATP/g COD/d'),
    'mu_max': (1.969, 'd^(-1)'),
    'q_CH': (0.594, 'g COD/g COD/d'),   
    'q_LI': (0.91, 'g COD/g COD/d'),    
    'V_NH': (0.254, 'g N/g COD/d'),   
    'V_NO': (0.254, 'g N/g COD/d'),  
    'V_P': (0.016, 'g P/g COD/d'),   
    }

#%%

def create_model(system=None):
    sys = system or create_system()

    model = qs.Model(system=sys, exception_hook='raise')
    param = model.parameter
    metric = model.metric
    
    PBR = sys.units[0]
    pm2 = PBR.model

    ##### General model with all uncertain variables #####

    # Add Uncertainty Parameters    
    get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))
        
    for k, v in baseline_values.items():
        b, units = v
        D = get_uniform_w_frac(b, 0.25)
        param(setter=DictAttrSetter(pm2.rate_function, '_params', k),
              name=k, element=PBR, kind='coupled', units=units, distribution=D)     # element = MIX? INF?

    for k, v in baseline_values_hifrac.items():
        b, units = v
        D = get_uniform_w_frac(b, 0.50)
        param(setter=DictAttrSetter(pm2.rate_function, '_params', k),
              name=k, element=PBR, kind='coupled', units=units, distribution=D)
        
    ##### Add universal evaluation metrics #####   
    @metric(name='Biomass', units='mg/L', element='')
    def get_X_ALG():
        return PBR.state['X_ALG']
    
    @metric(name='Carbohydrate', units='mg/L', element='')
    def get_X_CH():
        return PBR.state['X_CH']
    
    @metric(name='Lipid', units='mg/L', element='')
    def get_X_LI():
        return PBR.state['X_LI']
    
    @metric(name='Nitrogen', units='mg/L', element='')
    def get_S_NH():
        return PBR.state['S_NH']  
    
    @metric(name='Phosphorus', units='mg/L', element='')
    def get_S_P():
        return PBR.state['S_P']  
    
    return model

#%%

# =============================================================================
# Run uncertainty analysis
# =============================================================================

@time_printer
def run_uncertainty(model, N, T, t_step, method='LSODA', 
                    metrics_path='', timeseries_path='', 
                    rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    model = create_model()
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = metrics_path or os.path.join(results_path, f'table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = os.path.join(results_path, f'time_series_data_{seed}')
        os.mkdir(folder)
        tpath = os.path.join(folder, 'state.npy')    
    model.evaluate(
        state_reset_hook='reset_cache',
        t_span=t_span,
        t_eval=t_eval,
        method=method,
        export_state_to=tpath
        )
    model.table.to_excel(mpath)

#%% 
# =============================================================================
# UA with random initial conditions to test steady state
# =============================================================================

'''
def create_model(flowsheet=None):
    sys = create_system(flowsheet)

    unit = sys.flowsheet.unit
    MIX = unit.MIX

    model_ss = qs.Model(system=sys, exception_hook='raise')
    param_ss = model_ss.parameter
    # metric_ss = model_ss.metric
    
    get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))
    
    # Set initial conditions of all bioreactors
    for k, v in _ic.items():
        b = v
        D = get_uniform_w_frac(b, 0.5)
        @param_ss(name='initial '+k, element=MIX, kind='coupled', units='mg/L',
                  baseline=b, distribution=D)
        def ic_setter(conc): pass
    
    return model_ss

@time_printer
def run_wdiff_init(model, N, T, t_step, method='RK23', 
                   metrics_path='', timeseries_path='', 
                   rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    # model = create_model(kind='steady state')
    samples = model.sample(N=N, rule=rule)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'time_series_data_{seed}')
        os.mkdir(folder)
        tpath = os.path.join(folder, 'state.npy')
        
    cmps = model.system.units[0].components
      
    for i, smp in enumerate(samples):
        concs = dict(zip(cmps.IDs[0:-1], smp))
        for u in model._system.units:
            if u.ID in ('MIX', 'PBR1', 'PBR2', 'PBR3', 'PBR4', 'PBR5', 'PBR6', 'PBR7', 'PBR8', 'PBR9', 'PBR10', 
            'PBR11', 'PBR12', 'PBR13', 'PBR14', 'PBR15', 'PBR16', 'PBR17', 'PBR18', 'PBR19', 'PBR20', 'MEV', 'RET'):
                u.set_init_conc(**concs)
            
        model._system.simulate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method,
            export_state_to=tpath,
            sample_id=i,
            )
        
    for u in model._system.units:
        if u.ID in ('MIX', 'PBR1', 'PBR2', 'PBR3', 'PBR4', 'PBR5', 'PBR6', 'PBR7', 'PBR8', 'PBR9', 'PBR10', 
        'PBR11', 'PBR12', 'PBR13', 'PBR14', 'PBR15', 'PBR16', 'PBR17', 'PBR18', 'PBR19', 'PBR20', 'MEV', 'RET'):
            u.set_init_conc(**concs)
'''
 
#%%
if __name__ == '__main__':
    seed = 119
    t = 50
    t_step = 5
    n = 2
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    model = create_model()
    run_uncertainty(model, n, t, t_step, method=method, seed=seed)