# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan import ImpactIndicator as IInd, ImpactItem as IItm, \
    StreamImpactItem as SIItm, \
    TEA, LCA, PowerUtility
from qsdsan.utils import ospath, auom
from exposan.metab_mock import create_systems, data_path, results_path
import pandas as pd

IInd.load_from_file(ospath.join(data_path, 'TRACI_indicators.xlsx'), sheet=0)
IItm.load_from_file(ospath.join(data_path, '_impact_items.xlsx'))
bg_offset_CFs = IItm.get_item('biogas_offset').CFs

#%%
irr = 0.1
# p_ng = 0.65165 # USD/therm of natural gas charged to the brewery in 2016 (including tax and delivert etc.)
# p_ng = 0.85*auom('kJ').conversion_factor('therm') # [USD/kJ] 5.47 2021USD/Mcf vs. 4.19 2016USD/Mcf
op_hr = 365*24
get = getattr
PowerUtility.price = 0.0913 # 2021$ for minnesota area

def run_tea_lca(sys, save_report=True, info='', lt=30):
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
    F = sys.flowsheet
    inf = get(F.stream, f'BreweryWW_{F.ID}')
    eff = get(F.stream, f'Effluent_{F.ID}')
    cmps = eff.components
    KJ_per_kg = cmps.i_mass/cmps.chem_MW*cmps.LHV
    for ws in sys.products:
        if ws.phase == 'l':
            SIItm(ID=f'{ws.ID}_fugitive_ch4', linked_stream=ws, 
                  flow_getter=lambda s: s.imass['S_ch4']*cmps.S_ch4.i_mass,
                  GWP100=28, MIR=0.0143794871794872)
        elif ws.phase == 'g':
            fget = lambda s: -sum(s.mass*KJ_per_kg)*1e-3/39
            SIItm(ID=f'{ws.ID}_NG_offset', linked_stream=ws, functional_unit='m3',
                  flow_getter=fget, # m3/hr natural-gas-equivalent
                  **bg_offset_CFs)
    
    # Operation items
    kWh = lambda lca: lca.system.power_utility.consumption*lca.lifetime_hr
    def MJ(lca):
        sys = lca.system
        duties = [hu.duty for hu in sys.heat_utilities]
        heat = sum(d for d in duties if d > 0)
        cool = sum(d for d in duties if d < 0)
        MJ_per_hr = (heat + cool * 0.2)*1e-3 # assume 20% efficiency in recoverying heat from R1 eff (i.e., R2 inf)
        return MJ_per_hr*lca.lifetime_hr
    
    rCOD = inf.composite('COD', flow=True) - eff.composite('COD', flow=True)  # kgCOD/hr

    tea = TEA(sys, discount_rate=irr, lifetime=lt, simulate_system=False, CEPCI=708)
    levelized_cost = -tea.annualized_NPV/(rCOD * op_hr *1e-3)
    print(f'{sys.ID} TEA: ${round(levelized_cost, 1)}/ton rCOD')
    
    lca = LCA(
        sys, lifetime=lt, simulate_system=False,
        electricity = kWh,
        heat_onsite = MJ
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
    
    return sys

#%%
if __name__ == '__main__':
    # save = True
    save = False
    info = '35C_20yr_1yr'
    lt = 20
    sysA, sysB = create_systems(which=('A', 'B'))
    print(info)
    sysA = run_tea_lca(sysA, save, info, lt)
    # run_tea_lca(sysB, save, info, lt)
    # systems = create_systems(which=('C', 'D'))
    # for sys in systems:
    #     u = sys.units[0]
    #     for tau in (1, 2, 4, 8):
    #         u.V_liq = tau*5
    #         u.V_gas = tau*5/10
    #         info = f'HRT{tau}'
    #         run_tea_lca(sys, save, info, lt)
