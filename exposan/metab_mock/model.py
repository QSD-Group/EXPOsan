# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:08:42 2022

@author: joy_c
"""
import qsdsan as qs
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import ospath, time_printer
from exposan.metab_mock import system as s, results_path
from biosteam.evaluation._utils import var_indices
# import pandas as pd
import numpy as np
import os

__all__ = ('model_ua', 'run_uncertainty', 'cmps', 'H2E')


#%%
# =============================================================================
# model with uncertain ADM1 parameters
# =============================================================================
sysA, sysB = s.create_systems()
model_ua = qs.Model(system=sysB, exception_hook='raise')

########## Add Uncertainty Parameters ##########
param = model_ua.parameter
get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

ws_reg = sysB.flowsheet.stream
inf, eff, bg1, bg2 = ws_reg.BreweryWW_B, ws_reg.Effluent_B, ws_reg.biogas_1B, ws_reg.biogas_2B
cmps = inf.components

u_reg = sysB.flowsheet.unit
H2E, CH4E = u_reg.AnR1, u_reg.AnR2
adm1 = H2E.model
Ys_bl, mus_bl, Ks_bl = s.yields_bl, s.mus_bl, s.Ks_bl
n_Ys = len(Ys_bl)
n_mus = len(mus_bl)

for k,v in Ys_bl.items():
    b = v
    D = get_uniform_w_frac(b, 0.5)
    @param(name=k, element=H2E, kind='coupled', units='',
           baseline=b, distribution=D)
    def Y_setter(Y):
        pass

for i in range(len(mus_bl)):
    b = mus_bl[i]
    D = get_uniform_w_frac(b, 0.5)
    @param(name=f'mu_{adm1.IDs[i]}', element=H2E, kind='coupled', units='d^(-1)',
           baseline=b, distribution=D)
    def mu_setter(mu):
        pass

for j in range(len(Ks_bl)):
    b = Ks_bl[j]
    D = get_uniform_w_frac(b, 0.5)
    @param(name=f'K_{adm1.IDs[j+4]}', element=H2E, kind='coupled', units='kg/m3',
           baseline=b, distribution=D)
    def K_setter(K):
        pass
    
########## Add Evaluation Metrics ##########
metric = model_ua.metric

# Effluent composite variables and daily sludge production
@metric(name='H2E COD removal', units='%', element='H2E')
def get_H2E_rCOD():
    return (1 - H2E.outs[1].COD/inf.COD)*100

@metric(name='CH4E COD removal', units='%', element='CH4E')
def get_CH4E_rCOD():
    return (1 - eff.COD/H2E.outs[1].COD)*100

@metric(name='Total COD removal', units='%', element='System')
def get_rCOD():
    return (1 - eff.COD/inf.COD)*100

# @metric(name='H2E VFAs', units='g/L', element='H2E')
# def get_H2E_VFAs():
#     return H2E.outs[1].composite('COD', subgroup=('S_va', 'S_bu', 'S_pro', 'S_ac'))/1000

# @metric(name='H2E sCOD', units='g/L', element='H2E')
# def get_H2E_sCOD():
#     return H2E.outs[1].composite('COD', particle_size='s')/1000

@metric(name='Effluent COD', units='g/L', element='System')
def get_eff_COD():
    return eff.COD/1000

@metric(name='H2 production', units='kg/d', element='Biogas_1')
def get_QH2():
    return bg1.get_flow('kg/d', 'S_h2')*cmps.S_h2.i_mass

@metric(name='CH4 production', units='kg/d', element='Biogas_2')
def get_QCH4():
    return bg2.get_flow('kg/d', 'S_ch4')*cmps.S_ch4.i_mass

#%%
########### Functions for UA #################
@time_printer
def run_uncertainty(model, N, T, t_step, method='BDF', 
                    metrics_path='', timeseries_path='', 
                    rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
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
    index = model._index
    values = [None] * len(index)
    for i in index:
        smp = samples[i]
        Ys = dict(zip(Ys_bl.keys(), smp[:n_Ys]))
        adm1.set_parameters(**Ys)
        adm1.rate_function._params['rate_constants'][:] = smp[n_Ys: (n_Ys+n_mus)]
        adm1.rate_function._params['half_sat_coeffs'][:] = smp[(n_Ys+n_mus):]            
        model._system.simulate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method,
            export_state_to=tpath,
            sample_id=i,
            )
        values[i] = [i() for i in model.metrics]
    model.table[var_indices(model._metrics)] = values
    model.table.to_excel(mpath)