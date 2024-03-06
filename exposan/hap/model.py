# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import qsdsan as qs, qsdsan.stats as qst
from time import time
from exposan.hap import (
    create_system,
    results_path
    )
from chaospy import distributions as shape
from qsdsan.utils import AttrSetter, time_printer, ospath, load_data

__all__ = ('create_model', 'run_model', 'rerun_failed_samples')

#%%

def create_model(sys=None, exception_hook='warn', **kwargs):
    sys = sys or create_system(**kwargs)
    mdl = qs.Model(sys, exception_hook=exception_hook)
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    param = mdl.parameter
    metric = mdl.metric
    urination_rate = kwargs.pop('urination_rate', 1.4)
    
    # =============================================================================
    # uncertain parameters
    # =============================================================================
    b = 770
    D = shape.Triangle(470, b, 1070)
    @param(name='urine TP', units='mg/L', kind='coupled', 
            element='HAp', baseline=b, distribution=D)
    def set_TP(c):
        s.urine.imass['IP'] = c/1000*s.urine.F_vol
        
    b = 5
    D = shape.Triangle(b, b, b+2)
    N_rxt = u.HF.N
    tau_0 = u.HF.tau_0
    @param(name='HAp reactor design retention time', units='d', kind='coupled', 
           element='HAp', baseline=b, distribution=D)
    def set_HAp_tau(t):
        t *= 24
        u.HF.tau = (1-1/N_rxt) * t - tau_0/N_rxt
    
    b = 66
    D = shape.Triangle(30, b, 80)
    @param(name='HAp yield', units='% theoretical maximum', kind='coupled',
           element='HAp', baseline=b, distribution=D)
    def set_HAp_yield(p):
        u.HF.f_maximum_hap_yield = p/100
    
    b = 0.5
    D = shape.Triangle(b, b, b*2)
    param(AttrSetter(u.HF, 'inoculum_concentration'),
          name='Inoculum concentration', units='g/L', kind='coupled',
          element='HAp', baseline=b, distribution=D)
    
    b = 95
    D = shape.Uniform(90, 99)
    param(AttrSetter(u.HF, 'precipitate_moisture'),
          name='Precipitate moisture content', units='w%', kind='coupled',
          element='HAp', baseline=b, distribution=D)
    
    b = 1.45
    D = shape.Triangle(b*0.5, b, b*1.2)
    param(AttrSetter(u.YP, 'yield_on_sugar'),
          name='Osteoyeast biomass yield', units='kg/kg sugar fed', kind='coupled',
          element='Osteoyeast', baseline=b, distribution=D)
    
    b = 1.0
    D = shape.Uniform(0.9, 3.0)
    param(AttrSetter(u.YP, '_feedstock_cost_factor'),
          name='Yeast production feedstock cost factor', units='-', kind='isolated',
          element='Osteoyeast', baseline=b, distribution=D)
    
    b = 1.0
    D = shape.Uniform(0.75, 1.25)
    param(AttrSetter(u.YP, '_centrifuge_cost_factor'),
          name='Yeast production centrifuge cost factor', units='-', kind='isolated',
          element='Osteoyeast', baseline=b, distribution=D)
    
    b = 30
    D = shape.Uniform(15, 45)
    param(AttrSetter(u.CD, 'vehicle_vol_capacity'),
          name='Truck volume capacity', units='m3', kind='isolated',
          element='Logistics', baseline=b, distribution=D)
    
    b = 20
    D = shape.Uniform(15, 30)
    param(AttrSetter(u.CD, 'duration_per_location'),
          name='Service time per location', units='min', kind='isolated',
          element='Logistics', baseline=b, distribution=D)
    
    b = 200
    D = shape.Uniform(150, 250)
    param(AttrSetter(u.CD, 'vehicle_rental_price'),
          name='Truck rental cost', units='USD/day', kind='isolated',
          element='Logistics', baseline=b, distribution=D)
    
    b = 0.75
    D = shape.Uniform(0.5, 1.0)
    param(AttrSetter(u.CD, 'vehicle_fuel_cost'),
          name='Truck fuel cost', units='USD/mile', kind='isolated',
          element='Logistics', baseline=b, distribution=D)
    
    b = 45
    D = shape.Triangle(b-10, b, b+10)
    @param(name='Dryer cake solid content', units='w%', kind='coupled',
           element='Postprocessing', baseline=b, distribution=D)
    def set_csc(p):
        u.PP.dryer_cake_moisture = 100-p
    
    b = 1.0
    D = shape.Uniform(0.5, 1.5)
    param(AttrSetter(u.PP, '_dryer_cost_factor'),
          name='Dryer cost factor', units='-', kind='isolated',
          element='Postprocessing', baseline=b, distribution=D)
    
    b = 1.0
    D = shape.Uniform(1.0, 1.5)
    param(AttrSetter(u.PP, '_incinerator_cost_factor'),
          name='Incinerator cost factor', units='-', kind='isolated',
          element='Postprocessing', baseline=b, distribution=D)
    
    b = 5.23
    D = shape.Uniform(b, 6.97)
    param(AttrSetter(u.PP, 'auxiliary_fuel_price'),
          name='Incinerator fuel price', units='USD/gal', kind='isolated',
          element='Postprocessing', baseline=b, distribution=D)
    
    b = 90
    D = shape.DiscreteUniform(30, 150)
    @param(name='Number of deployed locations', units='-', kind='coupled',
           element='System', baseline=b, distribution=D)
    def set_N_loc(N):
        for unit in sys.units:
            unit.N_parallel_HApFermenter = N
    
    b = 45000
    D = shape.Uniform(10000, 80000)
    @param(name='Total capita equivalent served', units='pe', kind='coupled',
           element='System', baseline=b, distribution=D)
    def set_pe(pe):
        N_loc = u.HF.N_parallel_HApFermenter
        Q = pe * urination_rate / N_loc / 24 # L/hr
        s.urine.set_total_flow(Q, 'L/hr')
    
    b = 21.68
    D = shape.Uniform(18.50, 27.75)
    @param(name='Hourly labor wage', units='USD/hr', kind='coupled',
           element='System', baseline=b, distribution=D)
    def set_wage(w):
        for unit in (u.YP, u.CD, u.PP):
            unit.labor_wage = w
    
    b = 0.17
    D = shape.Uniform(0.12, 0.20)
    @param(name='Electricity cost', units='USD/kWh', kind='coupled',
           element='System', baseline=b, distribution=D)
    def set_pu_cost(c):
        qs.PowerUtility.price = c

    
    b = 35000
    D = shape.Uniform(28000, 50000)
    @param(name='Centralized facility rent', units='USD/yr', kind='coupled',
           element='System', baseline=b, distribution=D)
    def set_rent(r):
        sys.TEA.system_add_OPEX['Facility rent'] = r
    
    b = 28
    D = shape.Uniform(24, 32)
    @param(name='Income tax rate', units='%', kind='coupled',
           element='TEA', baseline=b, distribution=D)
    def set_tax(r):
        sys.TEA.income_tax = r/100
    
    # =============================================================================
    # metrics
    # =============================================================================
    
    @metric(name='HAp recovery rate', units='kg/yr', element='System')
    def get_rHAp():
        return s.product.imass['HAP'] * 24 * 365
    
    @metric(name='Average transportation duty', units='mile/location', element='System')
    def get_distance():
        return u.CD.design_results['Total travel distance']/u.CD.N_parallel_HApFermenter
    
    @metric(name='MPSP', units='USD/kg', element='TEA')
    def get_MPSP():
        p = sys.TEA.solve_price([s.product])
        return p * s.product.F_mass / s.product.imass['HAP']
    
    @metric(name='Annualized NPV', units='USD', element='TEA')
    def get_ANPV():
        return sys.TEA.annualized_NPV
    
    @metric(name='CAPEX', units='% ANPV', element='TEA')
    def get_capex():
        return sys.TEA.annualized_CAPEX / (-sys.TEA.annualized_NPV) * 100
    
    @metric(name='OPEX', units='% ANPV', element='TEA')
    def get_opex():
        return sys.TEA.AOC / (-sys.TEA.annualized_NPV) * 100
    
    @metric(name='HAp fermenter CAPEX', units='% total CAPEX', element='CAPEX')
    def get_hf_capex():
        return u.HF.installed_cost / sys.TEA.installed_equipment_cost * 100
    
    @metric(name='Yeast production CAPEX', units='% total CAPEX', element='CAPEX')
    def get_yp_capex():
        return u.YP.installed_cost / sys.TEA.installed_equipment_cost * 100
    
    @metric(name='Post processing CAPEX', units='% total CAPEX', element='CAPEX')
    def get_pp_capex():
        return u.PP.installed_cost / sys.TEA.installed_equipment_cost * 100
    
    def get_unit_opex(unit):
        feed = sum(i.cost for i in unit.ins if i.price)
        util = unit.utility_cost
        others = sum(unit.add_OPEX.values())
        return (feed + util + others)*24*365
    
    @metric(name='HAp fermenter OPEX', units='% total OPEX', element='OPEX')
    def get_hf_opex():
        return get_unit_opex(u.HF)/sys.TEA.AOC * 100
    
    @metric(name='Yeast production OPEX', units='% total OPEX', element='OPEX')
    def get_yp_opex():
        return get_unit_opex(u.YP)/sys.TEA.AOC * 100
    
    @metric(name='Transportation OPEX', units='% total OPEX', element='OPEX')
    def get_cd_opex():
        return get_unit_opex(u.CD)/sys.TEA.AOC * 100
    
    @metric(name='Post processing OPEX', units='% total OPEX', element='OPEX')
    def get_pp_opex():
        return get_unit_opex(u.PP)/sys.TEA.AOC * 100
   
    @metric(name='Other OPEX', units='% total OPEX', element='OPEX')
    def get_other_opex():
        return sum(sys.TEA.system_add_OPEX.values())/sys.TEA.AOC * 100
    
    return mdl

#%%
@time_printer
def run_model(mdl=None, samples=None, N=100, seed=None):
    mdl = mdl or create_model()
    if samples is None:
        if seed is None: seed = int(str(time())[-3:])
        samples = mdl.sample(N=N, rule='L', seed=seed)
    mdl.load_samples(samples)
    print(f'Seed = {seed}\n')
    print('Monte Carlo simulations begin...\n')
    mdl.evaluate()
    mpath = ospath.join(results_path, f'table_{seed}.xlsx')
    mdl.table.to_excel(mpath)
    return seed, mdl

def rerun_failed_samples(mdl=None, seed=None):
    if mdl is None:
        mdl = create_model()
        smp = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                        header=[0,1], skiprows=[2,])
    else:
        smp = mdl.table
    smp = smp[smp.isna().any(axis=1)]
    smp = smp.iloc[:, :len(mdl.parameters)].to_numpy()
    mdl.load_samples(smp)
    u = mdl.system.flowsheet.unit
    u.CD.solve_time = 100
    mdl.evaluate()
    mpath = ospath.join(results_path, f'table_{seed}_rerun.xlsx')
    mdl.table.to_excel(mpath)

def spearman(mdl=None, seed=None):
    if mdl is None:
        mdl = create_model()
        mdl.table = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                              header=[0,1], skiprows=[2,])
    path = ospath.join(results_path, f'spearman_{seed}.xlsx')
    rho, p = qst.get_correlations(mdl, kind='Spearman', file=path)
    return rho, p

#%%
if __name__ == '__main__':
    N = 2000
    seed = None
    mdl = run_model(N=N, seed=seed)
    rerun_failed_samples(mdl)
    r, p = spearman(mdl)
    # rerun_failed_samples(seed=292)
    # r, p = spearman(seed=292)