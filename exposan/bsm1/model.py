# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, AttrGetter, FuncGetter, \
    get_SRT as srt, time_printer
from exposan import bsm1 as bm
import pandas as pd
import numpy as np
import os

results_path = bm.results_path

__all__ = ('model_bsm1', 'model_2dv', 'model_ss',
           'run_uncertainty', 'analyze_timeseries', 'run_wdiff_init')

bsm1 = bm.bsm1
s = bm.system

#%%
# =============================================================================
# model with all uncertain variables
# =============================================================================

model_bsm1 = qs.Model(system=bsm1, exception_hook='raise')

########## Add Uncertainty Parameters ##########
param = model_bsm1.parameter
get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

cmps = s.cmps
PE = s.PE

b = 0.08
D = shape.Triangle(lower=0.04, midpoint=b, upper=0.12)
@param(name='Biomass N content i_XB', element=PE, kind='coupled',
       units='g N/g COD', baseline=b, distribution=D)
def set_i_XB(i):
    cmps.X_BH.i_N = cmps.X_BA.i_N = i
    cmps.refresh_constants()

b = 0.6
D = shape.Triangle(lower=0.57, midpoint=b, upper=0.63)
@param(name='Biomass products N content i_XP', element=PE, kind='coupled',
       units='g N/g COD', baseline=b, distribution=D)
def set_i_XP(i):
    cmps.X_P.i_N = cmps.X_I.i_N = i
    cmps.refresh_constants()
 
b = 0.75
D = shape.Triangle(lower=0.7, midpoint=b, upper=0.95)
@param(name='Organic particulates ash content fr_SS_COD', element=PE, kind='coupled',
       units='g SS/g COD', baseline=b, distribution=D)
def set_fr_SS_COD(i):
    cmps.X_I.i_mass = cmps.X_S.i_mass = cmps.X_P.i_mass = cmps.X_BH.i_mass = cmps.X_BA.i_mass = i
    cmps.refresh_constants()

Q, V_an, V_ae = s.Q, s.V_an, s.V_ae
A1, A2, O1, O2, O3, C1 = s.A1, s.A2, s.O1, s.O2, s.O3, s.C1

b = V_an * 2 / Q * 24
D = shape.Uniform(lower=2.34, upper=b)
@param(name='Anoxic zone hydraulic retention time', element=A1, 
       kind='coupled', units='hr', baseline=b, distribution=D)
def set_A1_A2_HRT(i):
    A1._V_max = A2._V_max = i / 24 * Q / 2

b = V_ae * 3 / Q * 24
D = shape.Uniform(lower=4.68, upper=b)
@param(name='Aerobic zone hydraulic retention time', element=O1, 
       kind='coupled', units='hr', baseline=b, distribution=D)
def set_O1_O2_O3_HRT(i):
    O1._V_max = O2._V_max = O3._V_max = i / 24 * Q / 3
    
b = 3
D = shape.Uniform(lower=2.25, upper=3.75)
@param(name='Internal recirculation rate as a fraction of influent', element=O3, 
       kind='coupled', units='', baseline=b, distribution=D)
def set_Q_intr(i):
    intr = i/(1+i+C1._Qras/Q)
    O3.split = [intr, 1-intr]

b = 1
D = shape.Uniform(lower=0.75, upper=1)
@param(name='Sludge recycling as a fraction of influent', element=C1, 
       kind='coupled', units='', baseline=b, distribution=D)
def set_Q_ras(i):
    C1.underflow = Q*i

b = s.Q_was
D = get_uniform_w_frac(b, 0.1)
@param(name='Waste sludge flowrate', element=C1, kind='coupled', units='m3/d',
       baseline=b, distribution=D)
def set_Q_was(i):
    C1.wastage = i


# Kinetic and stoichiometric parameters based on
# Sin, G.; Gernaey, K. V.; Neumann, M. B.; van Loosdrecht, M. C. M.; Gujer, W.
# Uncertainty Analysis in WWTP Model Applications:
# A Critical Discussion Using an Example from Design.
# Water Research 2009, 43 (11), 2894–2906.
# https://doi.org/10.1016/j.watres.2009.03.048.
asm1 = s.asm1
param_ranges = {
    'mu_H': (4, 3, 5, '/d'), # default, min, max, unit
    'K_S': (10, 5, 15, 'g COD/m3'),
    'K_O_H': (0.2, 0.1, 0.3, 'g O2/m3'),
    'K_NO': (0.5, 0.25, 0.75, 'g N/m3'),
    'b_H': (0.3, 0.285, 0.315, '/d'),
    'mu_A': (0.5, 0.475, 0.525, '/d'),
    'K_NH': (1, 0.5, 1.5, 'g N/m3'),
    'K_O_A': (0.4, 0.3, 0.5, 'g COD/m3'),
    'b_A': (0.05, 0.04, 0.06, '/d'),
    'eta_g': (0.8, 0.6, 1, ''),
    'k_a': (0.05, 0.03, 0.08, 'm3/g COD/d'),
    'k_h': (3, 2.25, 3.75, 'g X_S/g X_BH COD/d'),
    'K_X': (0.1, 0.075, 0.125, 'g X_S/g X_BH COD'),
    'eta_h': (0.8, 0.6, 1, ''),
    'Y_H': (0.67, 0.64, 0.7, 'g COD/g COD'),
    'Y_A': (0.24, 0.23, 0.25, 'g COD/g N'),
    }

asm1_baseline_param = asm1.parameters
for name, vals in param_ranges.items():
    b, lb, ub, unit = vals
    set_val = asm1_baseline_param.get(name)
    if set_val != b:
        warn(f'The value for parameter `{name}` set in the system script is {set_val}, '
             f'different from the provided baseline of {b}, is this intentional?')
    D = shape.Triangle(lower=lb, midpoint=b, upper=ub)
    setter = DictAttrSetter(asm1, '_parameters', keys=(name,))
    param(setter=setter, name=f'ASM1 {name}', element=A1, kind='coupled', 
          units=unit, baseline=b, distribution=D)

b = 0.2
D = shape.Triangle(lower=0.15, midpoint=b, upper=0.25)
@param(name='ASM1 f_Pobs', element=A1, kind='coupled', units='',
       baseline=b, distribution=D)
def set_f_Pobs(i):
    params = asm1._parameters
    Y_H = params['Y_H']
    params['f_P'] = i*(1-Y_H) / (1-i*Y_H)

# Aeration
aer1 = s.aer1
b = aer1.KLa
D = shape.Uniform(lower=180, upper=360)
@param(name='O1 and O2 KLa', element=O1, kind='coupled', units='',
       baseline=b, distribution=D)
def set_O1_O2_KLa(i):
    aer1.KLa = i

aer2 = s.aer2
b = aer2.KLa
D = get_uniform_w_frac(b, 0.1)
@param(name='O3 KLa', element=O3, kind='coupled', units='',
       baseline=b, distribution=D)
def set_O3_KLa(i):
    aer2.KLa = i
    
b = 8.0
D = get_uniform_w_frac(b, 0.1)
@param(name='Saturation DO', element=O1, kind='coupled', units='mg/L',
       baseline=b, distribution=D)
def set_DOsat(i):
    aer1.DOsat = i
    aer2.DOsat = i

########## Add Evaluation Metrics ##########
metric = model_bsm1.metric
SE, RAS, WAS, RE = s.SE, s.RAS, s.WAS, s.RE

# Effluent composite variables and daily sludge production
for i in ('COD', 'BOD5', 'TN'):
    metric(getter=AttrGetter(SE, attr=i), name='Effluent '+i, 
           units='mg/L', element='Effluent')

@metric(name='Effluent TKN', units='mg/L', element='Effluent')
def get_TKN():
    return SE.composite('N', subgroup=('S_NH', 'S_ND', 'X_ND'))

metric(getter=SE.get_TSS, name='Effluent TSS', units='mg/L', element='Effluent')

@metric(name='Daily sludge production', units='kg TSS/d', element='WAS')
def get_daily_sludge_production():
    return WAS.get_TSS() * 1e-3 * WAS.F_vol * 24

bio_IDs = ('X_BH', 'X_BA')
@metric(name='SRT', units='d', element='System')
def get_SRT():
    return srt(bsm1, bio_IDs)

#%%
########### Functions for UA #################
@time_printer
def run_uncertainty(model, N, T, t_step, method='LSODA', 
                    metrics_path='', timeseries_path='', 
                    rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
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

def analyze_timeseries(variable_getter, N, folder='', todf=True, **kwargs):
    outputs = {}
    for sample_id in range(N):       
        arr = np.load(os.path.join(folder, f'state_{sample_id}.npy'))
        outputs[sample_id] = variable_getter(arr, **kwargs)
    if todf: outputs = pd.DataFrame.from_dict(outputs)
    return outputs

#%%

# =============================================================================
# A new model for UA with 2 DVs as parameters
# =============================================================================

model_2dv = qs.Model(system=bsm1, exception_hook='raise')
param_2dv = model_2dv.parameter
metric_2dv = model_2dv.metric

b = aer1.Q_air
D = shape.Uniform(lower=2.4e3, upper=b*1.25)
@param_2dv(name='O1 O2 air flowrate', element=O1, kind='coupled', units='m3/d',
           baseline=b, distribution=D)
def set_Q_air(i):
    aer1.Q_air = i
    
b = s.Q_was
D = shape.Uniform(lower=300, upper=900)
param_2dv(setter=set_Q_was.setter, name='Waste sludge flowrate', element=C1, 
          kind='coupled', units='m3/d', baseline=b, distribution=D)

# DO of aerated reactors
DO_getter = lambda unit: float(unit.outs[0].iconc['S_O'])
for u in (O1, O2, O3):
    metric_2dv(getter=FuncGetter(DO_getter, [u]),
               name=u.ID+' DO', units='mg/L', element=u.ID)

# Effluent composite variables and daily sludge production
for i in ('COD', 'BOD5', 'TN'):
    metric_2dv(getter=AttrGetter(SE, attr=i), name='Effluent '+i, 
               units='mg/L', element='Effluent')

metric_2dv(getter=get_TKN.getter, name='Effluent TKN', units='mg/L', element='Effluent')
metric_2dv(getter=SE.get_TSS, name='Effluent TSS', units='mg/L', element='Effluent')
metric_2dv(getter=get_daily_sludge_production.getter,
           name='Daily sludge production', units='kg TSS/d', element='WAS')
metric_2dv(getter=get_SRT.getter, name='SRT', units='d', element='System')

#%%

# =============================================================================
# A new model for UA with random initial conditions to test steady state
# =============================================================================

model_ss = qs.Model(system=bsm1, exception_hook='raise')
param_ss = model_ss.parameter
metric_ss = model_ss.metric

# Set initial conditions of all bioreactors
_ic = s._init_conds
for k, v in _ic.items():
    i = cmps.index(k)
    b = v
    D = get_uniform_w_frac(b, 0.5)
    @param_ss(name='initial '+k, element=A1, kind='coupled', units='mg/L',
              baseline=b, distribution=D)
    def ic_setter(conc): pass


# Effluent composite variables and daily sludge production
for i in ('COD', 'BOD5', 'TN'):
    metric_ss(getter=AttrGetter(SE, attr=i), name='Effluent '+i, 
                units='mg/L', element='Effluent')

metric_ss(getter=get_TKN.getter, name='Effluent TKN', units='mg/L', element='Effluent')
metric_ss(getter=SE.get_TSS, name='Effluent TSS', units='mg/L', element='Effluent')
metric_ss(getter=get_daily_sludge_production.getter,
            name='Daily sludge production', units='kg TSS/d', element='WAS')
metric_ss(getter=get_SRT.getter, name='SRT', units='d', element='System')

@time_printer
def run_wdiff_init(model, N, T, t_step, method='LSODA', 
                   metrics_path='', timeseries_path='', 
                   rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    samples = model.sample(N=N, rule=rule)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    # mpath = metrics_path or os.path.join(results_path, f'table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = os.path.join(results_path, f'time_series_data_{seed}')
        os.mkdir(folder)
        tpath = os.path.join(folder, 'state.npy')
    for i, smp in enumerate(samples):
        concs = dict(zip(cmps.IDs[1:-2], smp))
        for u in model._system.units:
            if u.ID in ('A1', 'A2', 'O1', 'O2', 'O3'):
                u.set_init_conc(**concs)                
        model._system.simulate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method,
            export_state_to=tpath,
            sample_id=i,
            )