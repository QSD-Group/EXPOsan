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
           'create_modelC',
           'create_ss_model',
           'run_model',
           'run_modelB',
           'run_ss_model'
           )
#%%
sysA, sysB, sysC = s.create_systems()
Ys_bl, mus_bl, Ks_bl = s.yields_bl, s.mus_bl, s.Ks_bl
n_Ys = len(Ys_bl)
n_mus = len(mus_bl)

# =============================================================================
# model with uncertain parameters regarding gas extraction
# =============================================================================
get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

def add_degas_params(model, bioreactors, membranes, 
                     b_split=0.5, var_split=0.8, 
                     b_ermv=0.75, bounds_ermv=(0.25, 0.85)):
    param = model.parameter
    H2E, CH4E = bioreactors
    DM1, DM2 = membranes
    
    b = b_split
    D = get_uniform_w_frac(b, var_split)
    @param(name='H2E_sidestream_split', element=H2E, kind='coupled', units='',
           baseline=b, distribution=D)
    def H2E_split_setter(s):
        try: H2E.split = [s, 1-s]
        except: H2E.split = s

    D = get_uniform_w_frac(b, var_split)
    @param(name='CH4E_sidestream_split', element=CH4E, kind='coupled', units='',
           baseline=b, distribution=D)
    def CH4E_split_setter(s):
        try: CH4E.split = [s, 1-s]
        except: CH4E.split = s
    
    b = b_ermv
    D = shape.Uniform(*bounds_ermv)
    @param(name='H2_removal_efficiency', element=DM1, kind='coupled', units='',
           baseline=b, distribution=D)
    def H2_e_rmv(e):
        DM1.H2_degas_efficiency = e
        DM2.H2_degas_efficiency = e

    D = shape.Uniform(*bounds_ermv)        
    @param(name='CH4_removal_efficiency', element=DM1, kind='coupled', units='',
           baseline=b, distribution=D)
    def CH4_e_rmv(e):
        DM1.CH4_degas_efficiency = e
        DM2.CH4_degas_efficiency = e
        
    D = shape.Uniform(*bounds_ermv)
    @param(name='CO2_removal_efficiency', element=DM1, kind='coupled', units='',
           baseline=b, distribution=D)
    def CO2_e_rmv(e):
        DM1.CO2_degas_efficiency = e
        DM2.CO2_degas_efficiency = e

vfa_IDs = ('S_va', 'S_bu', 'S_pro', 'S_ac')

def add_metrics(model, biogas, wastewater, units):
    metric = model.metric
    inf, eff = wastewater
    R1, R2 = units
    S_h2_i_mass = eff.components.S_h2.i_mass
    S_ch4_i_mass = eff.components.S_ch4.i_mass
    
    @metric(name='R1 VFAs', units='g/L', element='Stage_1')
    def get_stage1_VFAs():
        return R1.outs[1].composite('COD', subgroup=vfa_IDs)/1000
    
    @metric(name='R1 H2', units='mg/L', element='Stage_1')
    def get_H2():
        return R1.outs[1].iconc['S_h2']*S_h2_i_mass
    
    @metric(name='R2 CH4', units='mg/L', element='Stage_2')
    def get_CH4():
        return eff.iconc['S_ch4']*S_ch4_i_mass
    
    @metric(name='Effluent COD', units='g/L', element='System')
    def get_eff_COD():
        return eff.COD/1000
    
    @metric(name='Total COD removal', units='%', element='System')
    def get_rCOD():
        return (1 - eff.COD/inf.COD)*100
    
    @metric(name='H2 production', units='kg/d', element='Biogas')
    def get_QH2():
        return sum([bg.imass['S_h2'] for bg in biogas])*S_h2_i_mass*24

    @metric(name='CH4 production', units='kg/d', element='Biogas')
    def get_QCH4():
        return sum([bg.imass['S_ch4'] for bg in biogas])*S_ch4_i_mass*24

def create_modelA():
    model = qs.Model(system=sysA, exception_hook='raise')
    ws_reg = sysA.flowsheet.stream
    inf, eff, bg1, bg2 = ws_reg.BreweryWW_A, ws_reg.Effluent_A, ws_reg.biogas_1A, ws_reg.biogas_2A
    u_reg = sysA.flowsheet.unit
    H2E, DM1, CH4E, DM2 = u_reg.H2E, u_reg.DM1, u_reg.CH4E, u_reg.DM2
    add_degas_params(model, (H2E, CH4E), (DM1, DM2))
    add_metrics(model, (bg1, bg2), (inf, eff), (H2E, CH4E))   
    return model
    
def create_modelC():
    model = qs.Model(system=sysC, exception_hook='raise')
    ws_reg = sysC.flowsheet.stream
    inf, eff, bgm1, bgm2, bgh1, bgh2 = (
        ws_reg.BreweryWW_C, 
        ws_reg.Effluent_C, 
        ws_reg.biogas_mem_1, 
        ws_reg.biogas_mem_2,
        ws_reg.biogas_hsp_1, 
        ws_reg.biogas_hsp_2
        )
    u_reg = sysC.flowsheet.unit
    R1, DM1, R2, DM2 = u_reg.R1, u_reg.DM1_c, u_reg.R2, u_reg.DM2_c
    add_degas_params(model, (R1, R2), (DM1, DM2))
    add_metrics(model, (bgm1, bgm2, bgh1, bgh2), (inf, eff), (R1, R2))
    return model

#%%
# =============================================================================
# model with uncertain ADM1 parameters
# =============================================================================

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
    add_metrics(model, (bg1, bg2), (inf, eff), (AnR1, AnR2))
    return model

#%%
# =============================================================================
# model with random initial conditions
# =============================================================================
def create_ss_model(system, R1_ID, R2_ID, 
                    R1_baseline_init_conds=None, 
                    R2_baseline_init_conds=None, 
                    frac_var=0.5):
    model = qs.Model(system, exception_hook='raise')
    _ic1 = R1_baseline_init_conds or s.R1_ss_conds
    _ic2 = R2_baseline_init_conds or s.R2_ss_conds
    u_reg = system.flowsheet.unit
    R1 = getattr(u_reg, R1_ID)
    R2 = getattr(u_reg, R2_ID)
    param = model.parameter
    for ic, u in zip((_ic1,_ic2), (R1,R2)):
        for k, v in ic.items():
            b = v
            D = get_uniform_w_frac(b, frac_var)
            @param(name=f'{k}_0', element=u, kind='coupled', units='mg/L',
                      baseline=b, distribution=D)
            def ic_setter(conc): pass
    return model

@time_printer
def run_ss_model(model, N, T, t_step, method='BDF', sys_ID=None,
                 R1_baseline_init_conds=None, 
                 R2_baseline_init_conds=None,
                 R1_ID='R1', R2_ID='R2',
                 metrics_path='', timeseries_path='', 
                 rule='L', seed=None, pickle=False):
    _ic1 = R1_baseline_init_conds or s.R1_ss_conds
    _ic2 = R2_baseline_init_conds or s.R2_ss_conds
    k1 = _ic1.keys()
    k2 = _ic2.keys()
    n_ic1 = len(k1)
    u_reg = model._system.flowsheet.unit
    R1 = getattr(u_reg, R1_ID)
    R2 = getattr(u_reg, R2_ID)
    if seed: np.random.seed(seed)
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'ss{sys_ID}_time_series_data_{seed}')
        os.mkdir(folder)
        tpath = ospath.join(folder, 'state.npy')
    for i, smp in enumerate(samples):
        ic1 = dict(zip(k1, smp[:n_ic1]))
        ic2 = dict(zip(k2, smp[n_ic1:]))
        R1.set_init_conc(**ic1)
        R2.set_init_conc(**ic2)
        model._system.simulate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method,
            export_state_to=tpath,
            sample_id=i,
            )
    R1.set_init_conc(**_ic1)
    R2.set_init_conc(**_ic2)
            
#%%
@time_printer
def run_model(model, N, T, t_step, method='BDF', sys_ID=None,
              metrics_path='', timeseries_path='',
              rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = metrics_path or ospath.join(results_path, f'sys{sys_ID}_table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'sys{sys_ID}_time_series_data_{seed}')
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