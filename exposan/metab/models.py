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
    create_system,
    C0_bulk, C1_bulk, C2_bulk,
    fermenters, methanogens,
    results_path
    )
import qsdsan as qs, os, numpy as np
from qsdsan.utils import FuncGetter, AttrSetter, AttrFuncSetter
from qsdsan.sanunits import AnaerobicCSTR

__all__ = ('add_discrete_dv', 
           'add_metrics',
           'create_model',
           'run_model')

#%%
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
            u.R1.n_layer = n_dz
            if n_stage == 2: 
                u.R2.bead_diameter = d
                u.R2.n_layer = n_dz
        
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
        u.R1._ODE = None
        if n_stage == 2: 
            u.R2.T = (273.15 + T + u.R2.T_air)/2
            u.R2._ODE = None

def add_metrics(model):
    metric = model.metric
    sys = model.system
    sub, = sys.subsystems
    n_stage = 1 if '1' in sys.ID else 2
    reactor_type, gas_xt = sys.ID.rstrip('_edg').split(str(n_stage))
    
    op_hr = 365*24
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    cmps = s.inf.components
    C0 = dict(zip(cmps.IDs, C0_bulk))
    C1 = dict(zip(cmps.IDs, C1_bulk))
    C2 = dict(zip(cmps.IDs, C2_bulk))
    h2_i_mass = cmps.S_h2.i_mass
    ch4_i_mass = cmps.S_ch4.i_mass
    
    @metric(name='COD removal', units='%', element='Process')
    def get_rcod():
        rcod = 1 - s.eff_dg.COD/s.inf.COD
        if rcod >= 0.7:
            u.R1._cache_state()
            if n_stage == 2: u.R2._cache_state()
        else:
            if n_stage == 1: u.R1.set_init_conc(**C0)
            else:
                u.R1.set_init_conc(**C1)
                u.R2.set_init_conc(**C2)
        return rcod*100
    
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
        get_encap_tss = lambda R, IDs: R.biomass_tss(IDs)[1]
        def get_overall_tss(R, IDs):
            return R.get_retained_mass(IDs)/(R.V_liq+R.V_beads)
        
        if n_stage == 1:
            for group, name in iters:
                metric(getter=FuncGetter(get_encap_tss, (u.R1, group)),
                       name=f'Encapsulated {name} TSS', **kwargs)
                metric(getter=FuncGetter(get_overall_tss, (u.R1, group)),
                       name=f'Overall {name} TSS', **kwargs)
        else:
            for Ri in (u.R1, u.R2):
                for group, name in iters:
                    metric(getter=FuncGetter(get_encap_tss, (Ri, group)),
                           name=f'{Ri.ID} encapsulated {name} TSS', **kwargs)
                    metric(getter=FuncGetter(get_overall_tss, (Ri, group)),
                           name=f'{Ri.ID} overall {name} TSS', **kwargs)
        for blt in (1, 10, 30):
            metric(getter=FuncGetter(get_cost, (sys, blt)),
                   name=f'Levelized cost (w/ degas, {blt}yr)', units='$/ton rCOD', element='TEA')
            metric(getter=FuncGetter(get_cost, (sub, blt)),
                   name=f'Levelized cost (w/o degas, {blt}yr)', units='$/ton rCOD', element='TEA')
            metric(getter=FuncGetter(get_gwp, (sys, blt)),
                   name=f'GWP100 (w/ degas, {blt}yr)', units='kg CO2eq/ton rCOD', element='LCA')
            metric(getter=FuncGetter(get_gwp, (sub, blt)),
                   name=f'GWP100 (w/o degas, {blt}yr)', units='kg CO2eq/ton rCOD', element='LCA')


def create_model(sys=None, kind='DV', **kwargs):
    sys = sys or create_system(**kwargs)
    mdl = qs.Model(sys, exception_hook='warn')
    if kind == 'DV':
        add_discrete_dv(mdl)
    add_metrics(mdl)
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

