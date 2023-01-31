# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs, qsdsan.sanunits as su, numpy as np, os
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import ospath, time_printer, MethodSetter, auom
from exposan.metab_mock.let23 import create_systems, results_path
# import pandas as pd

__all__ = ('create_models', 
           'run_model',
           'run_uncertainty')
#%%
ks = {
    'disintegration': (0.5, 3), # baseline, variation factor
    'hydrolysis_carbs': (10, 1),
    'hydrolysis_proteins': (10, 1),
    'hydrolysis_lipids': (10, 3),
    'uptake_LCFA': (6, 3),
    'uptake_propionate': (13, 1),
    'uptake_acetate': (8, 1),
    }
# bs = (0.02, 1)
# fb = (1, 1, 2) # baseline, lower bound, upper bound

Ks = {
      'uptake_LCFA': (0.4, 3),
      'uptake_c4': (0.2, 3),
      'uptake_propionate': (0.1, 1),
      'uptake_acetate': (0.15, 1),
      'uptake_h2': (7e-6, 1)
      }

# KIs_h2 = {
#     'uptake_propionate': (3.5e-6, 0.3),
#     }

# KI_nh3 = (0.0018, 0.3)
# pH_ULs = {
#     'uptake_acetate': (7, 0.3),
#     'uptake_h2': (6, 0.5)
#     }
# pH_LLs = {
#     'uptake_acetate': (6, 0.3),
#     'uptake_h2': (5, 0.3)
#     }

def get_uniform_from_fvar(b, f):
    u = b * (1+f)
    if f >= 1: l =  b*0.01
    else: l = b * (1-f)
    return shape.Uniform(l, u)

def add_params(model):
    param = model.parameter
    R1, R2 = units = [u for u in model.system.units if isinstance(u, su.AnaerobicCSTR)]
    
    #****** technological parameters *******
    b = 0.1
    D = shape.Uniform(0.1, 0.5)
    @param(name='Headspace_P', element=R1, kind='coupled', units='atm',
           baseline=b, distribution=D)
    def set_hsp_P(P):
        if R1.fixed_headspace_P: R1.headspace_P = P*1.01325
    
    b = 400
    D = shape.Uniform(1, 400)
    @param(name='Recirculation_rate', element=R1, kind='coupled', units='',
           baseline=b, distribution=D)
    def set_recir(r):
        if R1.split is not None: R1.split = [r, 1]
    
    b = 0.95
    D = shape.Uniform(0.9, 1.0)
    @param(name='Biomass_retention', element=R1, kind='coupled', units='',
           baseline=b, distribution=D)
    def set_f_retain(f):
        for u in units:
            u._f_retain = (u._f_retain > 0) * f
    
    #****** ADM1 parameters *********
    adm1 = R1.model
    for p, v in ks.items():
        b, f = v
        D = get_uniform_from_fvar(b, f)
        setter = MethodSetter(adm1, 'set_rate_constant', process=p)
        name = f'mu_{p}'
        param(setter=setter, name=name, element=adm1, kind='coupled', units='d^(-1)',
              baseline=b, distribution=D)
        
    b = 0.02
    D = get_uniform_from_fvar(b, 1)
    @param(name='Slow_decay_b', element=adm1, kind='coupled', 
           units='d^(-1)', baseline=b, distribution=D)
    def set_slow_b(k):
        adm1.rate_function._params['rate_constants'][12:19] = k
    
    b = 1
    D = shape.Uniform(1, 2)
    idx = adm1.indices(['decay_Xsu', 'decay_Xaa', 'decay_Xac'])
    @param(name='Fast_to_slow_b_ratio', element=adm1, kind='coupled', units='',
           baseline=b, distribution=D)
    def set_b_ratio(r):
        adm1.rate_function._params['rate_constants'][idx] *= r
    
    for p, v in Ks.items():
        b, f = v
        D = get_uniform_from_fvar(b, f)
        if p.endswith('_c4'):
            @param(name='K_uptake_c4', element=adm1, kind='coupled', units='kgCOD/m3',
                   baseline=b, distribution=D)
            def set_Ks_c4(K):
                adm1.set_half_sat_K(K, process='uptake_butyrate')
                adm1.set_half_sat_K(K, process='uptake_valerate')
        else:
            setter = MethodSetter(adm1, 'set_half_sat_K', process=p)
            param(setter=setter, name=f'K_{p}', element=adm1, kind='coupled',
                  units='kgCOD/m3', baseline=b, distribution=D)
    
    b = 3.5e-6
    D = get_uniform_from_fvar(b, 0.3)
    @param(name='KI_h2_uptake_propionate', element=adm1, kind='coupled', 
           units='kgCOD/m3', baseline=b, distribution=D)
    def set_KIh2_pro(KI):
        adm1.set_h2_inhibit_K(KI, process='uptake_propionate')
            
    b = 0.0018
    D = get_uniform_from_fvar(b, 0.3)
    @param(name='KI_nh3', element=adm1, kind='coupled', units='M', 
           baseline=b, distribution=D)
    def set_KI_nh3(KI):
        adm1.set_KI_nh3(KI)
    
    b = 7
    D = get_uniform_from_fvar(b, 0.3)
    @param(name='pH_UL_uptake_acetate',element=adm1, kind='coupled', units='', 
           baseline=b, distribution=D)
    def set_pH_UL_ac(ul):
        adm1.rate_function._params['pH_ULs'][-2] = ul
    
    b = 6
    D = get_uniform_from_fvar(b, 0.5)
    @param(name='pH_UL_uptake_h2',element=adm1, kind='coupled', units='', 
           baseline=b, distribution=D)
    def set_pH_UL_h2(ul):
        adm1.rate_function._params['pH_ULs'][-1] = ul    
    
    b = 6
    D = get_uniform_from_fvar(b, 0.3)
    @param(name='pH_LL_uptake_acetate',element=adm1, kind='coupled', units='', 
           baseline=b, distribution=D)
    def set_pH_LL_ac(ll):
        ULs = adm1.rate_function._params['pH_ULs']
        LLs = adm1.rate_function._params['pH_LLs']
        ul = ULs[-2]
        if ll < ul: LLs[-2] = ll
        elif ll == ul:
            ULs[-2] = ul + 0.5
            LLs[-2] = ll - 0.5
        else: 
            ULs[-2] = ll
            LLs[-2] = ul
        
    b = 5
    D = get_uniform_from_fvar(b, 0.3)
    @param(name='pH_LL_uptake_h2',element=adm1, kind='coupled', units='', 
           baseline=b, distribution=D)
    def set_pH_LL_h2(ll):
        ULs = adm1.rate_function._params['pH_ULs']
        LLs = adm1.rate_function._params['pH_LLs']
        ul = ULs[-1]
        if ll < ul: LLs[-1] = ll
        elif ll == ul:
            ULs[-1] = ul + 0.5
            LLs[-1] = ll - 0.5
        else: 
            ULs[-1] = ll
            LLs[-1] = ul
        
def add_metrics(model):
    metric = model.metric
    sys = model.system
    inf, = sys.feeds
    eff, = [ws for ws in sys.products if ws.ID.startswith('Effluent')]
    bgs = [ws for ws in sys.products if ws.ID.startswith('Biogas')]
    pumps = [u for u in sys.units if isinstance(u, su.Pump)]
    reactors = [u for u in sys.units if isinstance(u, su.AnaerobicCSTR)]
    
    cmps = eff.components
    S_h2_i_mass = cmps.S_h2.i_mass
    S_ch4_i_mass = cmps.S_ch4.i_mass
    cmps_i_COD = cmps.i_COD
    exclude_gas = cmps.s + cmps.x
    substrate_IDs = ('S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac',
                     'X_c', 'X_ch', 'X_pr', 'X_li')
    
    @metric(name='Soluble COD', units='mg/L', element='Effluent')
    def get_sCOD():
        return eff.composite('COD', particle_size='s')
    
    @metric(name='Particulate COD', units='mg/L', element='Effluent')
    def get_xCOD():
        return eff.composite('COD', particle_size='x')
    
    @metric(name='Substrate COD', units='mg/L', element='Effluent')
    def get_subCOD():
        return eff.composite('COD', subgroup=substrate_IDs)
    
    @metric(name='Total COD removal', units='%', element='System')
    def get_rCOD():
        return (1 - sum(eff.mass*cmps_i_COD*exclude_gas)/sum(inf.mass*cmps_i_COD*exclude_gas))*100
    
    @metric(name='Hourly COD removal', units='kgCOD/hr', element='System')
    def get_rCOD_hr():
        return sum(inf.mass*cmps_i_COD*exclude_gas) - sum(eff.mass*cmps_i_COD*exclude_gas) # kg/hr removal
    
    @metric(name='H2 production', units='kg/d', element='Biogas')
    def get_QH2():
        return sum([bg.imass['S_h2'] for bg in bgs])*S_h2_i_mass*24

    @metric(name='CH4 production', units='kg/d', element='Biogas')
    def get_QCH4():
        return sum([bg.imass['S_ch4'] for bg in bgs])*S_ch4_i_mass*24
    
    @metric(name='H2 yield', units='kg/kgCOD-removed', element='Biogas')
    def get_yH2():
        rCOD = sum(inf.mass*cmps_i_COD*exclude_gas) - sum(eff.mass*cmps_i_COD*exclude_gas) # kg/hr removal
        return sum([bg.imass['S_h2'] for bg in bgs])*S_h2_i_mass/rCOD

    @metric(name='CH4 yield', units='kg/kgCOD-removed', element='Biogas')
    def get_yCH4():
        rCOD = sum(inf.mass*cmps_i_COD*exclude_gas) - sum(eff.mass*cmps_i_COD*exclude_gas) # kg/hr removal
        return sum([bg.imass['S_ch4'] for bg in bgs])*S_ch4_i_mass/rCOD
    
    @metric(name='Pumping electricity', units='kW', element='Energy')
    def get_Epump():
        return sum(p.parallel['self']*p.power_utility.consumption for p in pumps)
    
    @metric(name='Heat', units='MJ/hr', element='Energy')
    def get_Eheat():
        return reactors[0].heat_utilities[0].duty/1e3 # kJ/hr to MJ/hr
    
    @metric(name='Biogas energy', units='MJ/hr', element='Energy')
    def get_E_biogas():
        # kmol/hr * J/mol = kJ/hr
        return sum(sum(bg.mass*cmps.i_mass/cmps.chem_MW*cmps.LHV) for bg in bgs)/1e3

    @metric(name='H2 energy', units='MJ/hr', element='Energy')    
    def get_E_h2():
        return sum(bg.imass['S_h2'] for bg in bgs)*S_h2_i_mass/cmps.S_h2.chem_MW*cmps.S_h2.LHV/1e3
    
    @metric(name='CH4 energy', units='MJ/hr', element='Energy')    
    def get_E_ch4():
        return sum(bg.imass['S_ch4'] for bg in bgs)*S_ch4_i_mass/cmps.S_ch4.chem_MW*cmps.S_ch4.LHV/1e3

    @metric(name='H2O energy', units='MJ/hr', element='Energy')    
    def get_E_h2o():
        return sum(bg.imass['H2O'] for bg in bgs)/cmps.H2O.MW*cmps.H2O.LHV/1e3

    conv = auom('kJ/hr').conversion_factor('kW')
    @metric(name='Net operation energy', units='kW', element='Energy')
    def get_net_E():
        consume = sum(p.parallel['self']*p.power_utility.consumption for p in pumps) \
            + reactors[0].heat_utilities[0].duty * conv
        produce = sum(sum(bg.mass*cmps.i_mass/cmps.chem_MW*cmps.LHV) \
                      for bg in bgs) * conv * 0.6 # assume 60% CHP efficiency
        return produce - consume

#%%
def create_models(systems=None, **kwargs):
    systems = systems or create_systems(**kwargs)
    models = []
    for sys in systems:
        mdl = qs.Model(sys, exception_hook='warn')
        add_params(mdl)
        add_metrics(mdl)
        models.append(mdl)
    return models

@time_printer
def run_model(model, N, T, t_step, method='BDF', 
              metrics_path='', timeseries_path='', prefix='',
              rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    sys_ID = model._system.ID
    samples = model.sample(N=N, rule=rule)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = metrics_path or ospath.join(results_path, f'{sys_ID}_{prefix}_table_{seed}.xlsx')
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'{sys_ID}_{prefix}_ytdata_{seed}')
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

def run_uncertainty(seed, N, models=None, T=400, t_step=5, method='BDF', selective=True):
    mdls = models or create_models(selective=selective)
    if selective: prefix = 'selective'
    else: prefix = 'nonselect'
    for mdl in mdls:
        msg = f'Model {mdl.system.ID[-1]} ({prefix})'
        print(f'\n{msg}\n{"-"*len(msg)}')
        print(f'Seed = {seed} \nN = {N} \nTime span 0-{T}d \n')
        run_model(mdl, N=N, T=T, t_step=t_step, prefix=prefix, seed=seed)
    # return mdls

#%%
if __name__ == '__main__':
    # mdls = create_models(which='D')
    # run_uncertainty(123, 10, models=mdls)
    run_uncertainty(123, 1000, selective=False)