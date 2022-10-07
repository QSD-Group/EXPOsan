# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import ospath, time_printer
from exposan.metab_mock import systems as s, results_path
from biosteam.evaluation._utils import var_indices
# import pandas as pd
import numpy as np
import os

__all__ = ('create_modelA', 
           'create_modelB',
           'run_modelA',
           'run_modelB'
           )

#%%
sysA, sysB = s.create_systems()
Ys_bl, mus_bl, Ks_bl = s.yields_bl, s.mus_bl, s.Ks_bl
n_Ys = len(Ys_bl)
n_mus = len(mus_bl)

# =============================================================================
# model with uncertain parameters regarding gas extraction
# =============================================================================

def add_degas_params(model, bioreactors, membranes):
    param = model.parameter
    H2E, CH4E = bioreactors
    DM1, DM2 = membranes
    
    b = 0.5
    D = shape.Uniform(0.1, 0.9)
    @param(name='H2E_sidestream_split', element=H2E, kind='coupled', units='',
           baseline=b, distribution=D)
    def H2E_split_setter(s):
        H2E.split = [s, 1-s]
    
    @param(name='CH4E_sidestream_split', element=CH4E, kind='coupled', units='',
           baseline=b, distribution=D)
    def CH4E_split_setter(s):
        CH4E.split = [s, 1-s]
    
    b = 0.75
    D = shape.Uniform(0.25, 0.85)
    @param(name='H2_removal_efficiency', element=DM1, kind='coupled', units='',
           baseline=b, distribution=D)
    def H2_e_rmv(e):
        DM1.H2_degas_efficiency = e
        DM2.H2_degas_efficiency = e
        
    @param(name='CH4_removal_efficiency', element=DM1, kind='coupled', units='',
           baseline=b, distribution=D)
    def CH4_e_rmv(e):
        DM1.CH4_degas_efficiency = e
        DM2.CH4_degas_efficiency = e
        
    @param(name='CO2_removal_efficiency', element=DM1, kind='coupled', units='',
           baseline=b, distribution=D)
    def CO2_e_rmv(e):
        DM1.CO2_degas_efficiency = e
        DM2.CO2_degas_efficiency = e

def add_metrics(model, biogas, wastewater):
    metric = model.metric
    bg1, bg2 = biogas
    inf, eff = wastewater
    S_h2_i_mass = eff.components.S_h2.i_mass
    S_ch4_i_mass = eff.components.S_ch4.i_mass
    
    @metric(name='Effluent COD', units='g/L', element='System')
    def get_eff_COD():
        return eff.COD/1000
    
    @metric(name='Total COD removal', units='%', element='System')
    def get_rCOD():
        return (1 - eff.COD/inf.COD)*100
    
    @metric(name='H2 production', units='kg/d', element='Biogas')
    def get_QH2():
        return (bg1.imass['S_h2'] + bg2.imass['S_h2'])*S_h2_i_mass*24

    @metric(name='CH4 production', units='kg/d', element='Biogas')
    def get_QCH4():
        return (bg1.imass['S_ch4'] + bg2.imass['S_ch4'])*S_ch4_i_mass*24


def create_modelA():
    model = qs.Model(system=sysA, exception_hook='raise')
    ws_reg = sysA.flowsheet.stream
    inf, eff, bg1, bg2 = ws_reg.BreweryWW_A, ws_reg.Effluent_A, ws_reg.biogas_1A, ws_reg.biogas_2A
    u_reg = sysA.flowsheet.unit
    H2E, DM1, CH4E, DM2 = u_reg.H2E, u_reg.DM1, u_reg.CH4E, u_reg.DM2
    add_degas_params(model, (H2E, CH4E), (DM1, DM2))
    add_metrics(model, (bg1, bg2), (inf, eff))   
    return model
    

#%%
# =============================================================================
# model with uncertain ADM1 parameters
# =============================================================================

get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

def add_adm_params(model, adm1, units):
    param = model.parameter
    AnR1, AnR2 = units    
    for k,v in Ys_bl.items():
        b = v
        D = get_uniform_w_frac(b, 0.5)
        @param(name=k, element=AnR1, kind='coupled', units='',
               baseline=b, distribution=D)
        def Y_setter(Y):
            pass
    
    for i in range(len(mus_bl)):
        b = mus_bl[i]
        D = get_uniform_w_frac(b, 0.5)
        @param(name=f'mu_{adm1.IDs[i]}', element=AnR1, kind='coupled', units='d^(-1)',
               baseline=b, distribution=D)
        def mu_setter(mu):
            pass
    
    for j in range(len(Ks_bl)):
        b = Ks_bl[j]
        D = get_uniform_w_frac(b, 0.5)
        @param(name=f'K_{adm1.IDs[j+4]}', element=AnR1, kind='coupled', units='kg/m3',
               baseline=b, distribution=D)
        def K_setter(K):
            pass

def create_modelB():
    model = qs.Model(system=sysB, exception_hook='raise')
    ws_reg = sysB.flowsheet.stream
    inf, eff, bg1, bg2 = ws_reg.BreweryWW_B, ws_reg.Effluent_B, ws_reg.biogas_1B, ws_reg.biogas_2B
    u_reg = sysB.flowsheet.unit
    AnR1, AnR2 = u_reg.AnR1, u_reg.AnR2
    adm1 = AnR1.model
    add_adm_params(model, adm1, (AnR1, AnR2))
    add_metrics(model, (bg1, bg2), (inf, eff))
    return model

#%%
@time_printer
def run_modelA(model, N, T, t_step, method='BDF',
               metrics_path='', timeseries_path='',
               rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = metrics_path or ospath.join(results_path, f'sysA_table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'sysA_time_series_data_{seed}')
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


@time_printer
def run_modelB(model, N, T, t_step, method='BDF', 
               metrics_path='', timeseries_path='', 
               rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = metrics_path or ospath.join(results_path, f'sysB_table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'sysB_time_series_data_{seed}')
        os.mkdir(folder)
        tpath = ospath.join(folder, 'state.npy')
    index = model._index
    values = [None] * len(index)
    adm1 = model._system.flowsheet.unit.AnR1.model
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