# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Saumitra Rai <raisaumitra9@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
import qsdsan as qs
from qsdsan import processes as pc, sanunits as su, Model as mod
import numpy as np
from qsdsan.utils import (
    ospath, 
    time_printer, 
    load_data, 
    get_SRT,
    get_P_blower, get_power_utility, 
    get_cost_sludge_disposal,
    get_normalized_energy, 
    get_daily_operational_cost, 
    get_total_operational_cost, 
    get_GHG_emissions_sec_treatment,
    get_GHG_emissions_discharge,
    get_GHG_emissions_electricity,
    get_CO2_eq_WRRF,
    get_total_CO2_eq,
    
    get_aeration_cost,
    get_pumping_cost,
    get_sludge_disposal_costs,
    get_CH4_CO2_eq_treatment, 
    get_N2O_CO2_eq_treatment,
    get_CH4_CO2_eq_discharge,
    get_N2O_CO2_eq_discharge,
    get_CO2_eq_electricity
    )
from chaospy import distributions as shape

# from . import data_path
folder = ospath.dirname(__file__)


__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm2d_kwargs', 'domestic_ww', 'default_init_conds',
    'Q_domestic', 'Q_brewery', 'Q_ras', 'Q_was', 'Temp', 'V_an', 'V_ae', 
    )

# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_domestic = 605000 # influent flowrate [m3/d]
Q_brewery = 12500
Temp = 273.15+20 # temperature [K]
V_an = 25210     # anoxic zone tank volume [m3] 
V_ae = 25210     # aerated zone tank volume [m3] 
Q_was = 11356    # sludge wastage flowrate (0.04% of 69 MGD) [m3/day] 
Q_ras = 249837   # recycle sludge flowrate (0.96% of 69 MGD) [m3/day] 
biomass_IDs = ('X_H', 'X_AUT')
ammonia_ID = 'S_NH4'

domestic_ww = {
   'S_I': 14,
   'X_I': 26.5,
   'S_F': 20.1,
   'S_A': 94.3,
   'X_S': 395,
   'S_NH4': 31,
   'S_N2': 0,
   'S_NO3': 0.266, 
   'S_PO4': 2.8,
   'X_PP': 0.05,
   'X_PHA': 0.5,
   'X_H': 0.15,
   'X_AUT': 0, 
   'X_PAO': 0, 
   'S_ALK':7*12,
    }

bp2_metab = {'S_O2': 0.0, 'S_N2': 0.0, 'S_NH4': 68.13798882904094, 'S_NO3': 0.0, 'S_PO4': 13.988754649480631, 'S_F': 49.65109745857063, 
             'S_A': 45.1984987392647, 'S_I': 91.60971385064741, 'S_ALK': 397.6072879786626, 'X_I': 271.6682810145552, 'X_S': 1436.4718241591547, 
             'X_H': 0.0, 'X_PAO': 0.0, 'X_PP': 0.0, 'X_PHA': 0.0, 'X_AUT': 0.0, 'X_MeOH': 0.0, 'X_MeP': 0.0}

default_asm2d_kwargs = dict(iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
            iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
            iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
            f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
            Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
            Y_A=0.24, f_XI_AUT=0.1,
            K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
            mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
            K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
            q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
            b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
            K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
            K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
            mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
            k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, 
            # path=os.path.join(data_path, '_asm2d.tsv'),
            )

default_init_conds = {
        'S_F':5,
        'S_A':2,
        'X_I':1000,
        'X_S':100,
        'X_H':500,
        'X_AUT':100,
        #'X_P':100,
        'S_O2':2,
        'S_NO3':20,
        'S_NH4':2,
        'S_ALK':7*12,
    }


def batch_init(sys, path, sheet):
    isa = isinstance
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in sys.units:
        if sheet.startswith('ss'):
            if isa(k, su.CSTR): k.set_init_conc(**dct[k.ID])
        else:
            if k.ID.startswith('O'):
                k.set_init_conc(**dct['O'])
            elif k.ID.startswith('A'):
                k.set_init_conc(**dct['A'])           
    c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    tss = [v for v in dct['C1_tss'].values() if v>0]
    u.C1.set_init_solubles(**c1s)
    u.C1.set_init_sludge_solids(**c1x)
    u.C1.set_init_TSS(tss)



#%%
def create_components():
     asm2d_cmps = pc.create_asm2d_cmps(False)
     asm2d_cmps.X_S.f_BOD5_COD = 0.54
     cmps = qs.Components([*asm2d_cmps])
     cmps.compile()
     return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={}, lifetime=100, discount_rate=0.1,
                  aeration_processes=()):
    # flowsheet = flowsheet or qs.Flowsheet('bsm1')
    # qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    #inf_kwargs = inf_kwargs or default_inf_kwargs
    dom_ww.set_flow_by_concentration(Q_domestic, 
                                     concentrations=domestic_ww, 
                                     units=('m3/d', 'mg/L'))
    
    ind_ww.set_flow_by_concentration(Q_brewery, 
                                     concentrations= bp2_metab, 
                                     units=('m3/d', 'mg/L'))
    
    effluent = qs.WasteStream('effluent', T=Temp)
    WAS = qs.WasteStream('WAS', T=Temp)
    RAS = qs.WasteStream('RAS', T=Temp)
    real_RAS = qs.WasteStream('real_RAS', T=Temp)
    only_recycle = qs.WasteStream('only_recycle', T=Temp)
    
  # kLa = {'O2': 317.0091531741395,
  # 'O3': 165.46166852781076,
  # 'O4': 143.2395289908076,
  # 'O5': 130.0786271277892,
  # 'O6': 113.9290436198542,
  # 'O7': 87.23226681670803,
  # 'O8': 55.120812633865086,
  # 'O9': 37.677433171670806,
  # 'O1': 98.96749331543272}
    
    # Process models
    # if aeration_processes:
    #     aer2, aer3, aer9 = aeration_processes
    # else:
    #     aer2 = pc.DiffusedAeration('aer2', 'S_O2', KLa=317, DOsat=8.0, V=V_ae)
    #     aer3 = pc.DiffusedAeration('aer3', 'S_O2', KLa=165.46, DOsat=8.0, V=V_ae)
    #     aer4 = pc.DiffusedAeration('aer4', 'S_O2', KLa=143.24, DOsat=8.0, V=V_ae)
    #     aer5 = pc.DiffusedAeration('aer5', 'S_O2', KLa=130.08, DOsat=8.0, V=V_ae)
    #     aer6 = pc.DiffusedAeration('aer6', 'S_O2', KLa=113.93, DOsat=8.0, V=V_ae)
    #     aer7 = pc.DiffusedAeration('aer7', 'S_O2', KLa=87.23, DOsat=8.0, V=V_ae)
    #     aer8 = pc.DiffusedAeration('aer8', 'S_O2', KLa=55.12, DOsat=8.0, V=V_ae)
    #     aer9 = pc.DiffusedAeration('aer9', 'S_O2', KLa=37.68, DOsat=8.0, V=V_ae)
    #     aer1 = pc.DiffusedAeration('aer1', 'S_O2', KLa=98.97, DOsat=8.0, V=V_ae)
        
    asm_kwargs = asm_kwargs or default_asm2d_kwargs
    
    asm2d = pc.ASM2d(iP_SF=0.005, iP_XS=0.005, iP_XI=0.005, iN_BM=0.1, iTSS_XI=0.72)
    
    # asm2d = pc.ASM2d()
    
    eff_TC1 = qs.WasteStream(ID = 'effluent_TC1')
    
    # PC = su.PrimaryClarifier(ID='PC', ins = [dom_ww, ind_ww, eff_TC1], 
    #                   outs = ['effluent_PC', 'sludge_PC'], isdynamic=True, 
    #                   init_with='WasteStream', thickener_perc= 0.09, 
    #                   TSS_removal_perc=69)
    
    PC = su.PrimaryClarifier(ID='PC', ins = [dom_ww, ind_ww, eff_TC1], 
                      outs = ['sludge_PC', 'effluent_PC'], isdynamic=True, 
                      init_with='WasteStream', thickener_perc= 0.09, 
                      TSS_removal_perc=69)
    
    eff_TC2 = qs.WasteStream(ID = 'effluent_TC2')
    eff_DU = qs.WasteStream(ID = 'effluent_DU')
    
    kwargs_A = dict(V_max=V_an, aeration=None, suspended_growth_model=asm2d)
    kwargs_O = dict(V_max=V_ae, aeration=2, DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O1 = su.CSTR('O1', RAS, **kwargs_O)    
    
    A1 = su.CSTR('A1', O1-0, real_RAS, **kwargs_A)    
    A2 = su.CSTR('A2',  [PC-1, only_recycle],  **kwargs_A)
    A3 = su.CSTR('A3', A2-0, **kwargs_A)    
    
    # O2 = su.CSTR('O2', A3-0, V_max=V_ae, aeration=aer2, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O3 = su.CSTR('O3', O2-0, V_max=V_ae, aeration=aer3, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O4 = su.CSTR('O4', O3-0, V_max=V_ae, aeration=aer4, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O5 = su.CSTR('O5', O4-0, V_max=V_ae, aeration=aer5, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O6 = su.CSTR('O6', O5-0, V_max=V_ae, aeration=aer6, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O7 = su.CSTR('O7', O6-0, V_max=V_ae, aeration=aer7, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O8 = su.CSTR('O8', O7-0, V_max=V_ae, aeration=aer8, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O9 = su.CSTR('O9', O8-0, ['treated'],V_max=V_ae, aeration=aer9, DO_ID='S_O2', suspended_growth_model=asm2d)
    O2 = su.CSTR('O2', A3-0, **kwargs_O)    
    O3 = su.CSTR('O3', O2-0, **kwargs_O)    
    O4 = su.CSTR('O4', O3-0, **kwargs_O)    
    O5 = su.CSTR('O5', O4-0, **kwargs_O)    
    O6 = su.CSTR('O6', O5-0, **kwargs_O)    
    O7 = su.CSTR('O7', O6-0, **kwargs_O)    
    O8 = su.CSTR('O8', O7-0, **kwargs_O)    
    O9 = su.CSTR('O9', O8-0, ['treated'], **kwargs_O)
    
    O1.temperature = 20
    O1.blower_efficiency = 0.7
    O1.inlet_pressure_loss = 1
    O1.diffuser_pressure_loss = 7
    O1.sludge_disposal_cost = 115
    O1.unit_electricity_costs = 0.0577
    
    C1 = su.FlatBottomCircularClarifier('C1', O9-0, [effluent, RAS, WAS],
                                        underflow=Q_ras, wastage=Q_was, 
                                        surface_area=24226, height=6, N_layer=10, 
                                        feed_layer=5)
    
    TC1 = su.Thickener('TC1', PC-0, ['sludge_TC1', eff_TC1], 
                       thickener_perc= 6.1, TSS_removal_perc=97.14)
    
    TC2 = su.Thickener('TC2', WAS, ['sludge_TC2', eff_TC2],
                       thickener_perc=4.3, TSS_removal_perc=95.54)
    
    DU = su.Centrifuge('DU', [TC1-0, TC2-0], ['sludge_DU', eff_DU],
                        thickener_perc=27, TSS_removal_perc=96.29)
    
    DU.DOC_f = 0.45
    DU.MCF = 0.8
    DU.k = 0.06 # temperate and boreal (dry)
    DU.pl = 30
    DU.elec_EF = 0.675 # the value for Minnesota
    DU.CH4_EF_st = 0.0075
    DU.N2O_EF_st = 0.0040
    DU.CH4_EF_dis = 0.009
    DU.N2O_EF_dis = 0.005
    
    Mixer = su.Mixer('Mixer', [eff_TC2, eff_DU, real_RAS], only_recycle)

    sys = qs.System('metro_ASM2d', 
                    path=(PC, A2, A3, O2, O3, O4, O5, O6, O7, O8, O9, 
                          C1, O1, A1, TC1, TC2, DU,  Mixer,), 
                    recycle = [eff_TC1, only_recycle])
    
    return sys
#%%

@time_printer
def run(t, method=None, **kwargs):
    sys = create_system()    
    path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
    batch_init(sys, path, 
                sheet='ss_alt')
                # sheet='t=10')
        
    # RAS = sys.flowsheet.stream.RAS
    # C1 = sys.flowsheet.unit.C1
    # sys.set_dynamic_tracker(RAS, C1)
    sys.set_dynamic_tracker(*sys.products)
    sys.set_dynamic_tracker(*sys.products)
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
    # sys._setup()
    # sys.converge() 
    
    return sys
    
if __name__ == '__main__':
    t = 0.001
    # method = 'RK45'
    method = 'RK23' 
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    sys = run(t, method=method)
    
    sys.diagram()
    
    act_units = [u.ID for u in sys.units if 
                  u.ID.startswith('O') or 
                 u.ID.startswith('A')
                 ]
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit
    srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------

    # # ACCOUNTING FOR AERATION COSTS (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    # aeration_units = [u for u in sys.units if u.ID.startswith('O')]
    
    # # # _dstate of units should be close to zero to imply steady state 
    
    # ACCOUNTING FOR AERATION COSTS (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    aeration_units = [u for u in sys.units if u.ID.startswith('O')]
    cmps = qs.get_thermo().chemicals
    asm2d = aeration_units[0].suspended_growth_model
    idx_DO = cmps.index('S_O2')
    DO_sat = 8
    kLa = {}
    Q_air = {}
    aer = pc.DiffusedAeration('aer', 'S_O2', KLa=240, DOsat=DO_sat, V=V_ae)
    for u in aeration_units:
        conc = u._state[:-1]
        OTR = - asm2d.production_rates_eval(conc)[idx_DO]
        DO = conc[idx_DO]
        kLa[u.ID] = aer.KLa = OTR / (DO_sat - DO)
        Q_air[u.ID] = aer.Q_air
    airflow = sum(Q_air.values())/24/60 # in m3/min
    
    print(f'Airflow = {airflow} m3/min')
    
    power_blower = get_P_blower(q_air=airflow)
    print(f'Required power at steady state is {power_blower:.2f} kW')

    
    # ACCOUNTING FOR PUMPING POWER (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    act_power_units = [u.ID for u in sys.units if \
    isinstance(u, (su.PrimaryClarifier, su.FlatBottomCircularClarifier, su.Thickener, su.Centrifuge))]
        
    required_pumping_power = get_power_utility(sys, active_unit_IDs=act_power_units)
    print(f'Required pumping (and other equipment) power at steady state is {required_pumping_power:.4f} kW\n')
    
    disposed_sludge = (fs.sludge_DU, )
    sludge_disposal_costs = get_cost_sludge_disposal(sludge = disposed_sludge, 
                                                      unit_weight_disposal_cost = 350)
    print(f'Sludge disposal cost = {sludge_disposal_costs} USD/day\n')

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
    
    influent_ST = (fs.effluent_PC,)
    effluent_ST = (fs.effluent, fs.WAS)
    GHG_ST = get_GHG_emissions_sec_treatment(influent = influent_ST, effluent = effluent_ST)
    print(f'CH4 and N2O emissions during secondary treatment equals {GHG_ST[0]} kg CH4/day and\
      {GHG_ST[1]} kg N2O-N/day respectively\n')
     
    effluent = (fs.effluent, )
    GHG_discharge = get_GHG_emissions_discharge(effluent = effluent)
    print(f'CH4 and N2O emissions at discharge equals {GHG_discharge[0]} kg CH4/day and \
      {GHG_discharge[1]} kg N2O-N/day respectively\n')
     
    GHG_electricity = get_GHG_emissions_electricity(system=sys, power_blower=power_blower, 
                                                    power_pump=required_pumping_power)
    print(f'CO2 emissions due to electricity consumption equals {GHG_electricity} kg-CO2-eq/day\n')
    
    # No GHG emission assumed from incineration 
    GHG_sludge_disposal = (0, 0)
    print(f'CH4 emissions due to sludge_disposal equals {GHG_sludge_disposal} kg-CH4/day\n')
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
    
    normalized_energy = get_normalized_energy(system = sys, aeration_power = power_blower, \
                                              pumping_power = required_pumping_power, \
                                              miscellaneous_power = 0)
    print(f'Normalized energy = {normalized_energy} kWh/m3\n')
    
    daily_operational_cost = get_daily_operational_cost(system= sys, aeration_power = power_blower, 
                            pumping_power = required_pumping_power, miscellaneous_power = 0, 
                            sludge_disposal_cost = sludge_disposal_costs)
    print(f'Daily operational costs = {daily_operational_cost} USD/m3 \n')
    
    total_daily_operational_cost = get_total_operational_cost(q_air = airflow, sludge = disposed_sludge, 
                                              system = sys, active_unit_IDs= act_power_units)
    print(f'Total daily operational costs = {total_daily_operational_cost} USD/m3 \n')
    
    
    CO2_eq_WRRF = get_CO2_eq_WRRF (system = sys, GHG_treatment = GHG_ST, GHG_discharge = GHG_discharge, \
                                    GHG_electricity = GHG_electricity, GHG_sludge_disposal = GHG_sludge_disposal)
    print(f'GHG emissions = {CO2_eq_WRRF} kg CO2 eq/m3\n')
    
    total_CO2_eq_WRRF = get_total_CO2_eq(system = sys, q_air = airflow, 
                                    influent_sc = influent_ST, effluent_sc = effluent_ST, effluent_sys = effluent, 
                                    active_unit_IDs=act_power_units, sludge= disposed_sludge)
    print(f'Total GHG emissions = {total_CO2_eq_WRRF} kg CO2 eq/m3\n')
    
    
    print(f'Airflow = {airflow} m3/min\n')
    
    sludge_prod = np.array([sludge.composite('solids', True, particle_size='x', unit='ton/d') \
                            for sludge in disposed_sludge]) # in ton/day
    sludge_prod = np.sum(sludge_prod)
    print(f'Total sludge produced = {sludge_prod} ton/day\n')
    
    industrial_COD = fs.industrial_wastewater.COD
    print(f'Industrial COD = {industrial_COD} mg/L\n')
    
    industrial_TN = fs.industrial_wastewater.TN
    print(f'Industrial TN = {industrial_TN} mg/L\n')
    
    effluent_COD = fs.effluent.COD
    print(f'Effluent COD = {effluent_COD} mg/L\n')
    
    effluent_TN = fs.effluent.TN
    print(f'Effluent TN = {effluent_TN} mg/L\n')
    
    fig, axis = fs.effluent.scope.plot_time_series(('S_F', 'S_A', 'S_NH4', 'S_NO3', 'X_I', 'S_I')) 
    fig
    
    # mass_degradable_sludge = np.array([slg.composite("C", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
    #               degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
    # mass_degradable_sludge = np.sum(mass_degradable_sludge)
    # print(f'Mass degradable sludge = {mass_degradable_sludge} kg/day\n')
    
    mass_N_sludge = np.array([slg.composite("N", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
                  degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
    mass_N_sludge = np.sum(mass_N_sludge)
    print(f'Mass N sludge = {mass_N_sludge} kg/day\n')
    
    # _dstate of units should be close to zero to imply steady state 
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
#     metro_model = mod(sys)
    
#     param = metro_model.parameter
    
#     # T : float
#     #     Air temperature [degree Celsius].
#     dist_T = shape.Uniform(lower = 15, upper = 25)
#     @param(name = 'Air temperature', 
#             element = fu.O1, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'degree celsius', 
#             baseline = fu.O1.temperature, # mean of upper and lower limits
#             distribution = dist_T)
#     def set_air_temperature(i):
#         fu.O1.temperature = i
        
#     # efficiency : float
#     #     Blower efficiency. Default is 0.7. 
#     dist_blower_efficiency = shape.Uniform(lower = 0.6, upper = 0.8)
#     @param(name = 'Blower efficiency', 
#             element = fu.O1, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'unitless', 
#             baseline = fu.O1.blower_efficiency, # mean of upper and lower limits
#             distribution = dist_blower_efficiency)
#     def set_blower_efficiency(i):
#         fu.O1.blower_efficiency = i
        
#     # P_inlet_loss : float
#     #     Head loss at inlet [kPa]. The default is 1 kPa. 
#     dist_P_inlet_loss = shape.Uniform(lower = 0.9, upper = 1.1)
#     @param(name = 'Head loss at inlet', 
#             element = fu.O1, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'kPa', 
#             baseline = fu.O1.inlet_pressure_loss, # mean of upper and lower limits
#             distribution = dist_P_inlet_loss)
#     def set_inlet_pressure_loss(i):
#         fu.O1.inlet_pressure_loss = i
        
#     # P_diffuser_loss : float
#     #     Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
#     dist_P_diffuser_loss = shape.Uniform(lower = 6.3, upper = 7.7)
#     @param(name = 'Head loss due to piping and diffuser', 
#             element = fu.O1, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'kPa', 
#             baseline = fu.O1.diffuser_pressure_loss, # mean of upper and lower limits
#             distribution = dist_P_diffuser_loss)
#     def set_diffuser_pressure_loss(i):
#         fu.O1.diffuser_pressure_loss = i
        
#     dist_disposal_costs = shape.Uniform(lower = 102, upper = 126)
#     @param(name = 'Sludge disposal unit cost', 
#             element = fu.O1, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'USD/ton', 
#             baseline = fu.O1.sludge_disposal_cost, # mean of upper and lower limits
#             distribution = dist_disposal_costs)
#     def set_sludge_disposal_cost(i):
#         fu.O1.sludge_disposal_cost = i
        
#     dist_electricity_cost = shape.Triangle(lower= 0.0577, midpoint= 0.075, upper=0.0878)
#     @param(name = 'Electricity unit cost', 
#             element = fu.O1, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'USD/kWh', 
#             baseline = fu.O1.unit_electricity_costs, # mean of upper and lower limits
#             distribution = dist_electricity_cost)
#     def set_electricity_cost(i):
#         fu.O1.unit_electricity_costs = i
    
# # # -----------------------------------------------------------------------------
# # # -----------------------------------------------------------------------------

#     # DOC_f : float, optional
#     #     fraction of DOC that can decompose (fraction).
#     dist_DOC_f = shape.Uniform(lower = 0.4, upper = 0.5)
#     @param(name = 'Decomposable fraction of DOC', 
#             element = fu.DU, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'unitless', 
#             baseline = fu.DU.DOC_f, # Deafault = 0.38 (Table 2.4 V5 Ch2, IPCC 2019)
#             distribution =  dist_DOC_f)
#     def set_DOC_f(i):
#         fu.DU.DOC_f = i
    
#     # MCF : float, optional
#     #     CH4 correction factor for aerobic decomposition in the year of deposition (fraction). The default is 0.8.
#     dist_MCF = shape.Uniform(lower = 0.4, upper = 1)
#     @param(name = 'MCF', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'unitless', 
#               baseline = fu.DU.MCF, # Table 3.1, V5 Ch3, IPCC 2019
#               distribution =  dist_MCF)
#     def set_MCF(i):
#           fu.DU.MCF = i
    
#     # k : TYPE, optional
#     #     Methane generation rate (k). The default is 0.06. (1/year)
#     dist_k = shape.Triangle(lower= 0.05, midpoint= 0.06, upper=0.2) # Table 3.3, V5 Ch3, IPCC 2019
#     @param(name = 'Methane generation rate', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = '1/year', 
#               baseline = fu.DU.k, 
#               distribution =  dist_k)
#     def set_k(i):
#           fu.DU.k = i
    
#     # pl : float, optional
#     #     The project lifetime over which methane emissions would be calculated. (years)
#     #     The default is 30 years.
#     dist_pl = shape.Uniform(lower = 30, upper = 50) # in years
#     @param(name = 'Project lifetime', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'years', 
#               baseline = fu.DU.pl, 
#               distribution =  dist_pl)
#     def set_pl(i):
#           fu.DU.pl = i
    
#     # CO2_EF : float
#     #     The emission factor used to calculate tier-2 CO2 emissions due to electricity consumption. 
#     #     The default is 0.675 kg-CO2-Eq/kWh.
#     dist_elec_EF = shape.Uniform(lower = 0.244, upper = 0.814) # in years
#     @param(name = 'EF for electrcity', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'kg-CO2-Eq/kWh', 
#               baseline = fu.DU.elec_EF, 
#               distribution =  dist_elec_EF)
#     def set_elec_EF(i):
#           fu.DU.elec_EF = i
         
#     # CH4_EF_sc :  float, optional.
#     #     The emission factor used to calculate methane emissions in secondary treatment. The default is 0.0075 kg CH4/ kg rCOD.
#     dist_CH4_EF_st = shape.Triangle(lower= 0.00075, midpoint=0.0075, upper=0.0225) # in kg CH4/kg CO2 r
#     @param(name = 'EF for CH4 during treatment', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'kg CH4/kg COD removed', 
#               baseline = fu.DU.CH4_EF_st, 
#               distribution =  dist_CH4_EF_st)
#     def set_CH4_EF_st(i):
#           fu.DU.CH4_EF_st = i
          
#     # N2O_EF_sc : float, optional
#     #     The emission factor used to calculate nitrous oxide emissions in secondary treatment. The default is 0.004 kg N2O-N/ kg N.
#     dist_N2O_EF_st = shape.Triangle(lower= 0.0001, midpoint=0.004, upper=0.0221) # in kg N2O/kg TN load
#     @param(name = 'EF for N2O during treatment', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'kg N/kg influent N', 
#               baseline = fu.DU.N2O_EF_st, 
#               distribution =  dist_N2O_EF_st)
#     def set_N2O_EF_st(i):
#           fu.DU.N2O_EF_st = i
        
#     # CH4_EF_discharge :  float, optional.
#     #     The emission factor used to calculate methane emissions in discharge. The default is 0.009 kg CH4/ kg effluent COD. 
#     dist_CH4_EF_dis = shape.Triangle(lower= 0.001, midpoint=0.009, upper=0.015) # in kg CH4/kg CO2 r
#     @param(name = 'EF for CH4 during discharge', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'kg CH4/kg effluent COD', 
#               baseline = fu.DU.CH4_EF_dis, 
#               distribution =  dist_CH4_EF_dis)
#     def set_CH4_EF_dis(i):
#           fu.DU.CH4_EF_dis = i
          
#     # # N2O_EF_discharge : float, optional
#     # #     The emission factor used to calculate nitrous oxide emissions in discharge. The default is 0.005 kg N2O-N/ kg effluent N.
#     # dist_N2O_EF_dis =  shape.Triangle(lower= 0.0005, midpoint=0.005, upper=0.075) # in kg CH4/kg CO2 r
#     # @param(name = 'EF for N2O during discharge', 
#     #           element = fu.DU, 
#     #           kind = 'coupled', # try 'isolated' too
#     #           units = 'kg N/kg effluent N', 
#     #           baseline = fu.DU.N2O_EF_dis, 
#     #           distribution =  dist_N2O_EF_dis)
#     # def set_N2O_EF_dis(i):
#     #       fu.DU.N2O_EF_dis = i
    
#     metric = metro_model.metric
    
#     @metric(name='Total operational cost', units='USD/day', element= fu.O1)
#     def operational_cost():
#         return get_total_operational_cost(q_air = airflow, sludge = disposed_sludge, 
#                                           system = sys, active_unit_IDs= act_power_units, 
                                          
#                                               # uncertain parameters 
                                              
#                                               T= fu.O1.temperature, 
#                                               efficiency= fu.O1.blower_efficiency, 
#                                               P_inlet_loss= fu.O1.inlet_pressure_loss, 
#                                               P_diffuser_loss= fu.O1.diffuser_pressure_loss, 
                                              
#                                               unit_weight_disposal_cost = fu.O1.sludge_disposal_cost,
#                                               unit_electricity_cost = fu.O1.unit_electricity_costs)
    
#     @metric(name='Aeration cost', units='USD/m3', element= fu.O1)
#     def aeration_cost():
#         return get_aeration_cost(q_air = airflow, system = sys, 
                                 
#                               # uncertain parameters 
                                 
#                               T= fu.O1.temperature, 
#                               efficiency= fu.O1.blower_efficiency, 
#                               P_inlet_loss= fu.O1.inlet_pressure_loss, 
#                               P_diffuser_loss= fu.O1.diffuser_pressure_loss, 
#                               unit_electricity_cost = fu.O1.unit_electricity_costs)
    
#     @metric(name='Pumping cost', units='USD/m3', element= fu.O1)
#     def pumping_cost():
#         return get_pumping_cost(system = sys, active_unit_IDs= act_power_units,  
                                
#                                 # uncertain parameter
#                               unit_electricity_cost = fu.O1.unit_electricity_costs)
    
#     @metric(name='Sludge disposal cost', units='USD/m3', element= fu.O1)
#     def sludge_disposal_costs():
#         return get_sludge_disposal_costs(sludge = disposed_sludge,  # sludge disposal costs 
#                                       system= sys, 
#                                       # uncertain parameters 
#                                       unit_weight_disposal_cost = fu.O1.sludge_disposal_cost # sludge disposal costs 
#                                       )
    

#     @metric(name='Total GHG emissions', units='kg CO2 eq./m3', element= fu.DU)
#     def GHG_emissions():
#         return get_total_CO2_eq(system = sys, q_air = airflow, 
#                                 influent_sc = influent_ST, effluent_sc = effluent_ST, effluent_sys = effluent, 
#                                 active_unit_IDs=act_power_units, sludge= None, 
                            
#                                 # uncertain parameters 
#                                 efficiency= fu.O1.blower_efficiency,  
#                                 P_inlet_loss= fu.O1.inlet_pressure_loss, 
#                                 P_diffuser_loss= fu.O1.diffuser_pressure_loss,
                                 
#                                 CH4_EF_sc = fu.DU.CH4_EF_st, 
#                                 N2O_EF_sc = fu.DU.N2O_EF_st, 
#                                 CH4_EF_discharge = fu.DU.CH4_EF_dis,
#                                 N2O_EF_discharge = 0.005,
                                 
#                                 CO2_EF= fu.DU.elec_EF, 
#                                 DOC_f = fu.DU.DOC_f, 
#                                 MCF = fu.DU.MCF, 
#                                 k = fu.DU.k, 
#                                 pl= fu.DU.pl)
    
#     @metric(name='CH4 CO2_eq_treatment', units='kg CO2 eq./m3', element= fu.DU)
#     def CH4_CO2_eq_treatment():
#         return get_CH4_CO2_eq_treatment(system = sys, influent_sc =influent_ST, effluent_sc = effluent_ST,
#                                   # uncertain parameters 
#                                 CH4_EF_sc = fu.DU.CH4_EF_st)
        
#     @metric(name='N2O CO2_eq_treatment', units='kg CO2 eq./m3', element= fu.DU)
#     def N2O_CO2_eq_treatment():
#         return get_N2O_CO2_eq_treatment(system = sys, influent_sc =influent_ST,  
#                                   N2O_EF_sc=fu.DU.N2O_EF_st)
        
#     @metric(name='CH4 CO2_eq_discharge', units='kg CO2 eq./m3', element= fu.DU)
#     def CH4_CO2_eq_discharge():
#         return get_CH4_CO2_eq_discharge(system = sys, effluent_sys = effluent, 
#                                   CH4_EF_discharge=fu.DU.CH4_EF_dis)
    
#     @metric(name='N2O CO2_eq_discharge', units='kg CO2 eq./m3', element= fu.DU)
#     def N2O_CO2_eq_discharge():
#         return get_N2O_CO2_eq_discharge(system = sys, effluent_sys = effluent,
#                                         N2O_EF_discharge  = 0.005)
    
#     @metric(name='CO2 eq_electricity', units='kg CO2 eq./m3', element= fu.DU)
#     def CO2_eq_electricity():
#         return get_CO2_eq_electricity(system = sys, q_air = airflow, active_unit_IDs=act_power_units, 
#                           # uncertain parameters 
#                           efficiency= fu.O1.blower_efficiency,  
#                           P_inlet_loss= fu.O1.inlet_pressure_loss, 
#                           P_diffuser_loss= fu.O1.diffuser_pressure_loss,
#                           CO2_EF= fu.DU.elec_EF)
    						
#     np.random.seed(3221) # setting the seed ensures you getting the same sample
    
#     #seed fixes the sample
#     samples = metro_model.sample(N=10000, rule='L')
#     metro_model.load_samples(samples)
    
#     # sys.isdynamic = False
#     # metro_model.evaluate()
    
#     all_outputs = []
#     for sample in samples:
#         for p, v in zip(metro_model.parameters, sample):
#             p.setter(v)
    
#         out = [m() for m in metro_model.metrics]
#         all_outputs.append(out)
    
#     metro_model.table.iloc[:,-10:] = all_outputs
    
#     metro_model.table

#     import pandas as pd
#     file_path = "/Users/saumitrarai/Desktop/Research/MS_research/metro/results/Metro_comprehensive_uncertainty_analysis.xlsx"  
#     sheet_name = "bp2_metab"  # Change this to your desired sheet name
#     with pd.ExcelWriter(file_path, engine='openpyxl', mode='a') as writer:
#         metro_model.table.to_excel(writer, sheet_name=sheet_name)
    
#     fig, ax = qs.stats.plot_uncertainties(metro_model, x_axis=metro_model.metrics[0], y_axis=metro_model.metrics[1],
#                                       kind='kde-kde', center_kws={'fill': True})
#     fig

#     r_df, p_df = qs.stats.get_correlations(metro_model, kind='Spearman')
#     print(r_df, p_df)
#     fig, ax = qs.stats.plot_correlations(r_df)
#     fig
    
#     sludge_prod = fs.sludge_DU.composite('solids', True, particle_size='x', unit='ton/d') 
#     print(f'sludge_prod = {sludge_prod}')