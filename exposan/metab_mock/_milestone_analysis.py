# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan import ImpactIndicator as IInd, ImpactItem as IItm, TEA, LCA, PowerUtility
from qsdsan.utils import ospath, auom
from exposan.metab_mock import create_systems, data_path, results_path
import pandas as pd

IInd.load_from_file(ospath.join(data_path, 'TRACI_indicators.xlsx'), sheet=0)
IItm.load_from_file(ospath.join(data_path, '_impact_items.xlsx'))
fugitive_ch4 = IItm(ID='fugitive_ch4', functional_unit='kg', 
                    GWP100=28, MIR=0.0143794871794872)

#%%
irr = 0.1
# p_treatment = 0.28 # assume USD/kg COD centralized treatment
p_ng = 0.688 # USD/therm of natural gas
op_hr = 365*24
get = getattr
PowerUtility.price = 0.0913

def run_tea_lca(sys, save_report=True, info='', lt=30):
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
    F = sys.flowsheet
    inf = get(F.stream, f'BreweryWW_{F.ID}')
    eff = get(F.stream, f'Effluent_{F.ID}')
    cmps = eff.components
    
    # Operation items
    fug_ch4 = eff.imass['S_ch4']*cmps.S_ch4.i_mass # kg/h
    kW = sys.power_utility.consumption
    duties = [hu.duty for hu in sys.heat_utilities]
    heat = sum(d for d in duties if d > 0)
    cool = sum(d for d in duties if d < 0)
    MJ_per_hr = (heat + cool * 0.2)*1e-3 # assume 20% efficiency in recoverying heat from R1 eff (i.e., R2 inf)
    offset = sum(sum(bg.mass*cmps.i_mass/cmps.chem_MW*cmps.LHV) \
               for bg in sys.products if bg.phase == 'g') # kJ/hr
    V_offset = offset * 1e-3 / 39      # m3/hr natural-gas-equivalent
    
    rCOD = inf.composite('COD', flow=True) - eff.composite('COD', flow=True)  # kgCOD/hr

    tea = TEA(
        sys, discount_rate=irr, lifetime=lt, simulate_system=False, 
        system_add_OPEX=dict(
            # COD_discharge_reduction = -rCOD*p_treatment*op_hr,
            NG_purchase_reduction = -auom('kJ').convert(offset,'therm')*p_ng*op_hr,
            )
        )
    levelized_cost = -tea.annualized_NPV/(rCOD * op_hr *1e-3)
    print(f'{sys.ID} TEA: ${round(levelized_cost, 1)}/ton rCOD')
    
    lca = LCA(
        sys, lifetime=lt, simulate_system=False,
        electricity = (kW*lt*op_hr, 'kWh'),
        heat_onsite = (MJ_per_hr*lt*op_hr, 'MJ'),
        fugitive_ch4 = (fug_ch4*lt*op_hr, 'kg'),
        biogas_offset = (-V_offset*lt*op_hr, 'm3')
        )
    
    gwp = lca.get_total_impacts()['GWP100']
    gwp_per_rcod = gwp/(rCOD * op_hr * lt *1e-3)
    print(f'{sys.ID} LCA: {round(gwp_per_rcod, 1)} kg CO2eq/ton rCOD')
    
    if save_report:
        tea_path = ospath.join(results_path, f'{F.ID}_{info}_tea.xlsx')
        summary = {
            'lifetime':lt,
            'NPV': tea.NPV,
            'annualized_CAPEX': tea.annualized_CAPEX,
            'annualized_OPEX': tea.AOC,
            'annual utility': tea.utility_cost,
            **tea.system_add_OPEX,
            'levelized_cost': levelized_cost,
            'GWP100': gwp_per_rcod
            }
        with pd.ExcelWriter(tea_path) as writer:
            summary = pd.DataFrame.from_dict({sys.ID: summary})
            summary.to_excel(writer, sheet_name='summary')
            tea.get_cashflow_table().to_excel(writer, sheet_name='cash_flow')
            for u in sys.units:
                rs = u.results()
                rs.index.names = ['l0', 'l1']
                rs.join(pd.DataFrame.from_dict({'lifetime': u.lifetime}),
                        on='l1')
                rs.to_excel(writer, sheet_name=u.ID)
        lca.save_report(ospath.join(results_path, f'{F.ID}_{info}_lca.xlsx'))
    
    # return sys

#%%
if __name__ == '__main__':
    save = True
    # save = False
    info = '35C_20yr_1yr'
    lt = 20
    sysA, sysB = create_systems(which=('A', 'B'))
    print(info)
    run_tea_lca(sysA, save, info, lt)
    run_tea_lca(sysB, save, info, lt)
    # systems = create_systems(which=('C', 'D'))
    # for sys in systems:
    #     u = sys.units[0]
    #     for tau in (1, 2, 4, 8):
    #         u.V_liq = tau*5
    #         u.V_gas = tau*5/10
    #         info = f'HRT{tau}'
    #         run_tea_lca(sys, save, info, lt)
