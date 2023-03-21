# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from exposan.metab import (
    Beads,
    create_system,
    # C0_bulk, C1_bulk, C2_bulk,
    fermenters, methanogens,
    results_path
    )
import qsdsan as qs, os, numpy as np
from chaospy import distributions as shape
from qsdsan.utils import FuncGetter, AttrSetter, AttrFuncSetter, MethodSetter, \
    SanUnitScope, categorize_cashflow, categorize_all_impacts
from qsdsan.sanunits import AnaerobicCSTR

__all__ = ('add_discrete_dv', 
           'add_continuous_params',
           'add_metrics',
           'create_model',
           'run_model')

#%%
def reset_init_conc(sys):
    cmps = sys.feeds[0].components
    n = len(cmps)
    u = sys.flowsheet.unit
    if '1' in sys.ID:
        C0 = dict(zip(cmps.IDs, u.R1._concs[:n]))
        u.R1._cached_state = None
        u.R1.set_init_conc(**C0)
    else:
        for unit in (u.R1, u.R2):
            C0 = dict(zip(cmps.IDs, unit._concs[:n]))
            unit._cached_state = None
            unit.set_init_conc(**C0)

def add_discrete_dv(model):
    param = model.parameter
    sys = model.system
    n_stage = 1 if '1' in sys.ID else 2
    reactor_type, gas_xt = sys.ID.rstrip('_edg').split(str(n_stage))
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
       
    if reactor_type in ('PB', 'FB'):
        @param(name='Bead diameter', units='mm', kind='coupled', element='Beads')
        def set_db(d):
            n_dz = 10 if d > 5 else 5
            u.R1.bead_diameter = d
            if n_stage == 1:
                if n_dz != u.R1.n_layer:
                    u.R1.n_layer = n_dz
                    reset_init_conc(sys)
            else: 
                u.R2.bead_diameter = d
                if n_dz != u.R2.n_layer:
                    u.R2.n_layer = n_dz
                    reset_init_conc(sys)
        
        if reactor_type == 'FB':
            @param(name='Voidage', units='', kind='coupled', element='Beads')
            def set_fvoid(f):
                u.R1.voidage = f
                if n_stage == 2: u.R2.voidage = f
    
    if gas_xt == 'H':
        param(setter=AttrSetter(u.R1, 'headspace_P'),
              name='Headspace pressure', units='bar', kind='coupled', element='Gas')
    elif gas_xt == 'M':
        param(setter=AttrFuncSetter(u.R1, 'split', lambda r: [r,1]),
              name='Recirculation ratio', units='', kind='coupled', element='Gas')

    @param(name='Total HRT', units='d', kind='coupled', element='Reactors')
    def set_tau(tau):
        V = s.inf.F_vol * 24 * tau
        if n_stage == 1:
            u.R1.V_liq = V
            u.R1.V_gas = V*0.1
        else:
            u.R1.V_liq = V/12
            u.R1.V_gas = V/12*0.1
            u.R2.V_liq = V*11/12
            u.R2.V_gas = V*11/12*0.1

    
    @param(name='Temperature', units='C', kind='coupled', element='Reactors')
    def set_temp(T):
        u.R1.T = 273.15 + T
        if reactor_type in ('FB','PB'): u.R1._prep_model()
        u.R1._compile_ODE()
        u.R1.scope = SanUnitScope(u.R1)
        if n_stage == 2: 
            u.R2.T = (273.15 + T + u.R2.T_air)/2
            if reactor_type in ('FB','PB'): u.R2._prep_model()
            u.R2._compile_ODE()
            u.R2.scope = SanUnitScope(u.R2)

#%%
def add_continuous_params(model):
    param = model.parameter
    sys = model.system
    n_stage = 1 if '1' in sys.ID else 2
    reactor_type, gas_xt = sys.ID.rstrip('_edg').split(str(n_stage))
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    adm1 = u.R1.model
    
    # config-specific uncertainty/DVs
    b = 0.963
    D = shape.Triangle(0.87, b, 0.995)
    @param(name='UASB solid retention efficacy', units='', kind='coupled', 
           element='UASB', baseline=b, distribution=D)
    def set_f_retain(f):
        if reactor_type == 'UASB':
            u.R1._f_retain = (u.R1._f_retain > 0) * f
            if n_stage == 2:
                u.R2._f_retain = (u.R2._f_retain > 0) * f
    
    b = 16
    D = shape.LogUniform(11, 22)
    @param(name='Max encapsulation density', units='gTSS/L', kind='coupled',
           element='Encapsulation', baseline=b, distribution=D)
    def set_max_tss(tss):
        if reactor_type in ('PB', 'FB'):
            u.R1.max_encapsulation_tss = tss
            if n_stage == 2: u.R2.max_encapsulation_tss = tss
    
    b = 1420
    D = shape.Uniform(990, 1860)
    @param(name='Bead density', units='kg/m3', kind='coupled', 
           element='Encapsulation', baseline=b, distribution=D)
    def set_rho_b(rho):
        Beads._bead_density = rho
    
    b = 0.55
    D = shape.LogUniform(0.2, 1.1)
    @param(name='Bead-to-water diffusivity fraction', units='', kind='coupled',
           element='Encapsulation', baseline=b, distribution=D)
    def set_f_diff(f):
        if reactor_type in ('PB', 'FB'):
            u.R1.f_diff = f
            if n_stage == 2: u.R2.f_diff = f
            
    b = 10
    D = shape.Uniform(1, 10)
    @param(name='Bead lifetime', units='yr', kind='coupled', 
           element='Encapsulation', baseline=b, distribution=D)
    def set_blt(lt):
        if reactor_type in ('PB', 'FB'):
            u.R1.bead_lifetime = lt
            if n_stage == 2: u.R2.bead_lifetime = lt
    
    b = 5
    D = shape.Uniform(1, 5)
    @param(name='Bead diameter', units='mm', kind='coupled', 
           element='Encapsulation', baseline=b, distribution=D)
    def set_db(d):
        if reactor_type in ('PB', 'FB'):
            n_dz = 10 if d > 5 else 5
            u.R1.bead_diameter = d
            if n_stage == 2: u.R2.bead_diameter = d
            if n_dz != u.R1.n_layer:
                u.R1.n_layer = n_dz
                if n_stage == 2: u.R2.n_layer = n_dz
                reset_init_conc(sys)
    b = 0.75
    D = shape.Uniform(0.6, 0.9)
    @param(name='FB voidage', units='', kind='coupled', element='FB',
           baseline=b, distribution=D)
    def set_FB_void(f):
        if reactor_type == 'FB':
            u.R1.voidage = f
            if n_stage == 2: u.R2.voidage = f
    
    b = 1.5
    D = shape.Uniform(0.5, 2)
    @param(name='FB height-to-diameter', units='', kind='coupled', element='FB',
           baseline=b, distribution=D)
    def set_FB_h2d(r):
        if reactor_type == 'FB':
            u.R1.reactor_height_to_diameter = r
            if n_stage == 2: u.R2.reactor_height_to_diameter = r
            
    b = 0.39
    D = shape.Triangle(0.35, b, 0.45)
    @param(name='PB voidage', units='', kind='coupled', element='PB',
           baseline=b, distribution=D)
    def set_PB_void(f):
        if reactor_type == 'PB':
            u.R1.voidage = f
            if n_stage == 2: u.R2.voidage = f

    b = 1.5
    D = shape.Uniform(1, 4)
    @param(name='PB height-to-diameter', units='', kind='coupled', element='PB',
           baseline=b, distribution=D)
    def set_PB_h2d(r):
        if reactor_type == 'PB':
            u.R1.reactor_height_to_diameter = r
            if n_stage == 2: u.R2.reactor_height_to_diameter = r
    
    # common uncertainty/DVs 
    ks = (
        ('k_fa', 'uptake_LCFA', 6, 0.6, 24),
        ('k_pro', 'uptake_propionate', 13, 1.3, 26),
        ('k_ac', 'uptake_acetate', 8, 0.8, 16)
        )
    
    for name, pid, b, lb, ub in ks:
        D = shape.Triangle(lb, b, ub)
        param(setter=MethodSetter(adm1, 'set_rate_constant', key='k', process=pid),
              name=name, units='COD/COD/d', kind='coupled', element='ADM1',
              baseline=b, distribution=D)
    
    Ks = (
        ('K_pro', 'uptake_propionate', 0.1, 0.01, 0.2),
        ('K_ac', 'uptake_acetate', 0.15, 0.015, 0.3)
        )
    
    for name, pid, b, lb, ub in Ks:
        D = shape.Triangle(lb, b, ub)
        param(setter=MethodSetter(adm1, 'set_half_sat_K', key='K', process=pid),
              name=name, units='COD/COD/d', kind='coupled', element='ADM1',
              baseline=b, distribution=D)
    
    start = adm1._find_index('decay_Xsu')
    b = 0.02
    D = shape.Triangle(0.002, b, 0.04)    
    @param(name='k_dec', units='d^(-1)', kind='coupled', element='ADM1',
           baseline=b, distribution=D)
    def set_k_dec(k):
        adm1.rate_function._params['rate_constants'][start:] = k

    b = 0.5
    D = shape.Uniform(0, 1)
    @param(name='Degassing factor', units='-', kind='coupled', element='Membrane',
           baseline=b, distribution=D)
    def set_f_degas(f):
        eh2 = 0.55 + f*(0.7-0.55)
        ech4 = 0.36 + f*(0.55-0.36)
        eco2 = 0.06 + f*(0.2-0.06)
        u.DMe._h2_ermv = eh2
        u.DMe._ch4_ermv = ech4
        u.DMe._co2_ermv = eco2
        if gas_xt == 'M':
            u.DMs._h2_ermv = eh2
            u.DMs._ch4_ermv = ech4
            u.DMs._co2_ermv = eco2
    
    b = 1
    D = shape.Uniform(1/6, 2)
    @param(name='Total HRT', units='d', kind='coupled', element='System',
           baseline=b, distribution=D)
    def set_tau(tau):
        V = s.inf.F_vol * 24 * tau
        if n_stage == 1:
            u.R1.V_liq = V
            u.R1.V_gas = V*0.1
        else:
            u.R1.V_liq = V/12
            u.R1.V_gas = V/12*0.1
            u.R2.V_liq = V*11/12
            u.R2.V_gas = V*11/12*0.1
        if reactor_type in ('FB','PB'): u.R1._prep_model()
        u.R1._compile_ODE()
        u.R1.scope = SanUnitScope(u.R1)
        if n_stage == 2: 
            if reactor_type in ('FB','PB'): u.R2._prep_model()
            u.R2._compile_ODE()
            u.R2.scope = SanUnitScope(u.R2)

#%%
def add_metrics(model, kind='DV'):
    metric = model.metric
    sys = model.system
    sub, = sys.subsystems
    n_stage = 1 if '1' in sys.ID else 2
    reactor_type, gas_xt = sys.ID.rstrip('_edg').split(str(n_stage))
    
    op_hr = sys.operating_hours
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    cmps = s.inf.components
    h2_i_mass = cmps.S_h2.i_mass
    ch4_i_mass = cmps.S_ch4.i_mass
    
    get_encap_tss = lambda R, IDs: R.biomass_tss(IDs)[1]
    def get_overall_tss(R, IDs):
        return R.get_retained_mass(IDs)/(R.V_liq+R.V_beads)
    
    isa = isinstance
    def get_cost(system, bead_lt=None):
        s = system.flowsheet.stream
        qrcod = (s.inf.COD - s.eff_dg.COD) * s.eff_dg.F_vol * 1e-6 # ton/hr
        if bead_lt: 
            for unit in system.units:
                if isa(unit, AnaerobicCSTR): unit.bead_lifetime = bead_lt
        return -system.TEA.annualized_NPV/(qrcod * op_hr)
    
    def get_gwp(system, bead_lt=None):
        s = system.flowsheet.stream
        qrcod = (s.inf.COD - s.eff_dg.COD) * s.eff_dg.F_vol * 1e-6  # ton/hr
        if bead_lt: 
            for unit in system.units:
                if isa(unit, AnaerobicCSTR): unit.bead_lifetime = bead_lt
        gwp_per_hr = system.LCA.get_total_impacts()['GWP100']/system.LCA.lifetime_hr       
        return gwp_per_hr/qrcod
    
    @metric(name='COD removal', units='%', element='Process')
    def get_rcod():
        rcod = 1 - s.eff_dg.COD/s.inf.COD
        if rcod > 0.8:
            u.R1._cache_state()
            if n_stage == 2: u.R2._cache_state()
        else:
            u.R1._cached_state = None
            if n_stage == 2: u.R2._cached_state = None
        return rcod*100
    
    if kind == 'DV':
        @metric(name='H2 yield', units='kg H2/kg rCOD', element='Process')
        def get_yh2():
            qh2 = sum([bg.imass['S_h2'] for bg in sys.products if bg.phase == 'g'])*h2_i_mass # kg H2/h
            qrcod = (s.inf.COD - s.eff_dg.COD) * s.eff_dg.F_vol * 1e-3 # kg COD/h
            return qh2/qrcod
    
        @metric(name='CH4 yield', units='kg CH4/kg rCOD', element='Process')
        def get_ych4():
            qch4 = sum([bg.imass['S_ch4'] for bg in sys.products if bg.phase == 'g'])*ch4_i_mass
            qrcod = (s.inf.COD - s.eff_dg.COD) * s.eff_dg.F_vol * 1e-3 # kg COD/h
            return qch4/qrcod
        
        kwargs = dict(units='g/L', element='Biomass')   
        if reactor_type == 'UASB':
            if n_stage == 1:
                metric(getter=FuncGetter(u.R1.biomass_tss, (fermenters,)),
                       name='Fermenters TSS', **kwargs)
                metric(getter=FuncGetter(u.R1.biomass_tss, (methanogens,)),
                       name='Methanogens TSS', **kwargs)
            else:
                for Ri in (u.R1, u.R2):
                    metric(getter=FuncGetter(Ri.biomass_tss, (fermenters,)),
                           name=f'{Ri.ID} fermenters TSS', **kwargs)
                    metric(getter=FuncGetter(Ri.biomass_tss, (methanogens,)),
                           name=f'{Ri.ID} methanogens TSS', **kwargs)
            metric(getter=FuncGetter(get_cost, (sys,)),
                   name='Levelized cost (w/ degas)', units='$/ton rCOD', element='TEA')
            metric(getter=FuncGetter(get_cost, (sub,)),
                   name='Levelized cost (w/o degas)', units='$/ton rCOD', element='TEA')
            metric(getter=FuncGetter(get_gwp, (sys,)),
                   name='GWP100 (w/ degas)', units='kg CO2eq/ton rCOD', element='LCA')
            metric(getter=FuncGetter(get_gwp, (sub,)),
                   name='GWP100 (w/o degas)', units='kg CO2eq/ton rCOD', element='LCA')
        else:
            iters = zip((fermenters, methanogens), ('fermenters', 'methanogens'))
            if n_stage == 1:
                for group, name in iters:
                    metric(getter=FuncGetter(get_encap_tss, (u.R1, group)),
                           name=f'Encapsulated {name} TSS', **kwargs)
                    metric(getter=FuncGetter(get_overall_tss, (u.R1, group)),
                           name=f'Overall {name} TSS', **kwargs)
            else:
                for group, name in iters:
                    metric(getter=FuncGetter(get_encap_tss, (u.R1, group)),
                           name=f'R1 encapsulated {name} TSS', **kwargs)
                    metric(getter=FuncGetter(get_overall_tss, (u.R1, group)),
                           name=f'R1 overall {name} TSS', **kwargs)
                    metric(getter=FuncGetter(get_encap_tss, (u.R2, group)),
                           name=f'R2 encapsulated {name} TSS', **kwargs)
                    metric(getter=FuncGetter(get_overall_tss, (u.R2, group)),
                           name=f'R2 overall {name} TSS', **kwargs)
            for blt in (1, 10, 30):
                metric(getter=FuncGetter(get_cost, (sys, blt)),
                       name=f'Levelized cost (w/ degas, {blt}yr)', units='$/ton rCOD', element='TEA')
                metric(getter=FuncGetter(get_cost, (sub, blt)),
                       name=f'Levelized cost (w/o degas, {blt}yr)', units='$/ton rCOD', element='TEA')
                metric(getter=FuncGetter(get_gwp, (sys, blt)),
                       name=f'GWP100 (w/ degas, {blt}yr)', units='kg CO2eq/ton rCOD', element='LCA')
                metric(getter=FuncGetter(get_gwp, (sub, blt)),
                       name=f'GWP100 (w/o degas, {blt}yr)', units='kg CO2eq/ton rCOD', element='LCA')
    
    else:
        _cached_metrics = {}
        @metric(name='H2 production', units='kg/d', element='Biogas')
        def get_QH2():
            return sum([bg.imass['S_h2'] for bg in sys.products if bg.phase == 'g'])*h2_i_mass*24

        @metric(name='CH4 production', units='kg/d', element='Biogas')
        def get_QCH4():
            return sum([bg.imass['S_ch4'] for bg in sys.products if bg.phase == 'g'])*ch4_i_mass*24
        
        biomass = (*fermenters, *methanogens)
        kwargs = dict(units='g/L', element='Biomass')
        metric(getter=FuncGetter(get_overall_tss, (u.R1, biomass)),
               name='R1 Overall TSS', **kwargs)
        if reactor_type in ('FB', 'PB'):
            metric(getter=FuncGetter(get_encap_tss, (u.R1, biomass)),
                   name='R1 Encapsulated TSS', **kwargs)
        if n_stage == 2:
            metric(getter=FuncGetter(get_overall_tss, (u.R2, biomass)),
                   name='R2 Overall TSS', **kwargs)
            if reactor_type in ('FB', 'PB'):
                metric(getter=FuncGetter(get_encap_tss, (u.R2, biomass)),
                       name='R2 Encapsulated TSS', **kwargs)
        
        def get_vessel_cost(suffix, system):
            llc_bd = categorize_cashflow(system.TEA)
            tot = llc_bd.pop('total')
            _cached_metrics[f'llc_{suffix}'] = share = {k: v/tot*100 for k, v in llc_bd.items()}
            return share['vessel']
        
        def get_other_cost(suffix, item):
            share = _cached_metrics[f'llc_{suffix}']
            return share[item]

        def get_vessel_gwp(suffix, system):
            gwp_bd = categorize_all_impacts(system.LCA)
            tot = abs(gwp_bd.pop('total'))
            _cached_metrics[f'gwp_{suffix}'] = share = {k: v/tot*100 for k, v in gwp_bd.items()}
            return share['vessel']
        
        def get_other_gwp(suffix, item):
            share = _cached_metrics[f'gwp_{suffix}']
            return share[item]       
        
        other_items = ('beads', 'dm', 'others', 'electricity', 'heat_onsite', 'chemicals', 'biogas_offset')
        
        for suffix, system in (('(w/ degas)', sys), ('(w/o degas)', sub)):
            metric(getter=FuncGetter(get_cost, (system,)),
                   name=f'Levelized cost {suffix}', units='$/ton rCOD', element='TEA')

            metric(getter=FuncGetter(get_vessel_cost, params=(suffix, system)),
                   name=f'llc vessel {suffix}', units='%', element='TEA')
            for i in other_items:
                metric(getter=FuncGetter(get_other_cost, params=(suffix, i)),
                       name=f'llc {i} {suffix}', units='%', element='TEA')

            metric(getter=FuncGetter(get_gwp, (system,)),
                   name=f'GWP100 {suffix}', units='kg CO2eq/ton rCOD', element='LCA')
            
            metric(getter=FuncGetter(get_vessel_gwp, params=(suffix, system)),
                   name=f'gwp vessel {suffix}', units='%', element='LCA')
            for i in (*other_items, 'fug_ch4'):
                metric(getter=FuncGetter(get_other_gwp, params=(suffix, i)),
                       name=f'gwp {i} {suffix}', units='%', element='LCA')


#%%
def create_model(sys=None, kind='DV', exception_hook='warn', **kwargs):
    sys = sys or create_system(**kwargs)
    mdl = qs.Model(sys, exception_hook=exception_hook)
    if kind == 'DV': add_discrete_dv(mdl)
    else: add_continuous_params(mdl)
    add_metrics(mdl, kind=kind)
    return mdl

def run_model(model, sample, T=400, t_step=10, method='BDF', 
              mpath='', tpath=''):
    name = model.system.ID.rstrip('_edg')
    model.load_samples(sample)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = mpath or os.path.join(results_path, f'{name}.xlsx')
    if not tpath:
        folder = os.path.join(results_path, f'ty_data_{name}')
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

