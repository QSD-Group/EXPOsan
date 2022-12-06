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
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import ospath, time_printer
from exposan.pm2 import (
    biomass_IDs,
    create_system, 
    default_init_conds as _ic, 
    results_path,
    Q, V_mix, V_pbr, V_mem, V_ret,
    )

__all__ = ('create_model', 'run_wdiff_init',)

'''
## bsm1 ###
from qsdsan.utils import DictAttrSetter, AttrGetter, FuncGetter, \
    get_SRT as srt, time_printer
from exposan.bsm1 import (
    biomass_IDs,
    create_system,
    default_init_conds as _ic,
    results_path,
    Q, Q_was, V_an, V_ae,
    )

__all__ = ('create_model', 'run_uncertainty', 'run_wdiff_init',)
'''

#%%

# =============================================================================
# UA with random initial conditions to test steady state
# =============================================================================

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

'''
    sys = system or create_system(flowsheet)
    stream = sys.flowsheet.stream
    wastewater, effluent, WAS = stream.wastewater, stream.effluent, stream.WAS
    cmps = wastewater.components
    
    unit = sys.flowsheet.unit
    A1, A2, O1, O2, O3, C1 = unit.A1, unit.A2, unit.O1, unit.O2, unit.O3, unit.C1
    aer1, aer3 = O1.aeration, O3.aeration
    
    model = qs.Model(system=sys, exception_hook='raise')
    param = model.parameter
    metric = model.metric
    
    # Add Uncertainty Parameters
    get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))


    ##### UA with random initial conditions to test steady state #####
    elif kind in ('steady state', 'ss'):
        
        # Set initial conditions of all bioreactors
        for k, v in _ic.items():
            # i = cmps.index(k)
            b = v
            D = get_uniform_w_frac(b, 0.5)
            @param(name='initial '+k, element=A1, kind='coupled', units='mg/L',
                      baseline=b, distribution=D)
            def ic_setter(conc): pass
'''


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
 
#%%
if __name__ == '__main__':
    seed = 119
    t = 50
    t_step = 5
    n = 2
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    model_ss = create_model()
    run_wdiff_init(model_ss, n, t, t_step, method=method, seed=seed)