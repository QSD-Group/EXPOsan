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
import numpy as np
from qsdsan import processes as pc, sanunits as su, Model as mod
from chaospy import distributions as shape

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
    get_GHG_emissions_sludge_disposal,
    get_CO2_eq_WRRF,
    get_total_CO2_eq,
    
    get_aeration_cost,
    get_pumping_cost,
    get_sludge_disposal_costs,
    get_CH4_CO2_eq_treatment, 
    get_N2O_CO2_eq_treatment,
    get_CH4_CO2_eq_discharge,
    get_N2O_CO2_eq_discharge,
    get_CH4_emitted_during_pl, 
    get_CH4_emitted_after_pl,
    get_CO2_eq_electricity
    )

# from exposan.bsm1 import data_path
# from . import data_path
folder = ospath.dirname(__file__)

__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm2d_kwargs', 
    'steady_state_ad_init_conds',
    'domestic_ww', 'Q_domestic', 'Q_brewery',
    'Q_ras', 'Q_was', 'Temp', 'V_ae',
    )
# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_domestic = 38000      # influent flowrate [m3/d] = 10 MGD (WERF report)
Q_brewery = 785         # industrial wastewater [m3/d] flowrate is determined based on a constant 
                        # domestic/industrial flowrate ratio across data from Metropolitan WRRF
                        
Temp = 273.15+20 # temperature [K]

Q_was = 378.54 # [m3/day] 0.1 MGD (WERF report)
Q_ras = 0.40*(39147) # [m3/day] 40% of effluent from PC (WERF report)
V_ae = 4982

biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

domestic_ww = {
   'S_I': 20,
   'X_I': 40,
   'S_F': 45,
   'S_A': 63,
   'X_S': 160,
   'S_NH4': 25,
   'S_PO4': 4.5,
   'X_PP': 0,
   'X_PHA': 10,
   'X_H': 10,
   'X_AUT': 10, 
   'X_PAO': 5, 
   'X_MeOH': 32, 
   'S_ALK':7*12,
    }

bm2 = {'S_O2': 0.0, 'S_N2': 0.0, 'S_NH4': 126.63, 'S_NO3': 0.0, 'S_PO4': 0.0, 'S_F': 0.0, 'S_A': 2404.0000000000027, 'S_I': 52.00000000000006, 
       'S_ALK': 164.79681703664826, 'X_I': 54.00000000000003, 'X_S': 0.0, 'X_H': 0.0, 'X_PAO': 0.0, 'X_PP': 0.0, 'X_PHA': 0.0, 'X_AUT': 0.0, 'X_MeOH': 0.0, 
       'X_MeP': 0.0}

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

# Using effluent DG sludge concentrations from steady state sys consisting of only AD

steady_state_ad_init_conds = {
      'S_su': 1.094e+01,
       # 'S_aa': 4.713e+00, # original ss
      # 'S_aa': 4.713e+01,
       'S_aa': 4.713e+03,
      'S_fa': 9.150e+01,
      'S_va': 9.512e+00,
      'S_bu': 1.255e+01,
      'S_pro': 1.486e+01,
      'S_ac': 3.749e+01,
       'S_h2': 2.193e-04,
      'S_ch4': 7.692e+01,
      'S_IC': 6.708e+02,
      'S_IN': 1.035e+03,
      # --
      'S_IP': 1.035e+00, #random
      # --
      'S_I': 2.325e+02,
      # --
      # 'X_c': 2.171e+02,
      # --
       # 'X_ch': 6.432e+01, # original ss
      # 'X_ch': 6.432e+00, 
       'X_ch': 6.432e-01,
      'X_pr': 6.719e+01,
       # 'X_li': 1.436e+02, # original ss
      # 'X_li': 1.436e+01,
       'X_li': 1.436e+00,
      'X_su': 1.116e+03,
      'X_aa': 8.755e+02,
      'X_fa': 1.274e+03,
      'X_c4': 3.730e+02,
      'X_pro': 1.751e+02,
      'X_ac': 1.405e+03,
      'X_h2': 6.820e+02,
      'X_I': 2.413e+04,
      # --
      'X_PHA': 10.000e-01, 
      'X_PP':  10.000e-01,  
      
      'X_PAO': 10.000e+00, 
      
      'S_K': 10.000e-001, 
      'S_Mg': 10.000e-001, 
      'X_MeOH': 10.000e-001, 
      'X_MeP': 10.000e-001, 
      # --
       'S_cat': 0.000e+00, 
       'S_an': 4.772e+00
      }

def batch_init(sys, path, sheet):
    # isa = isinstance
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    u.DG.set_init_conc(**steady_state_ad_init_conds)
    # for k in sys.units:
    #     if k.ID.startswith('O', 'A'): 
    #         k.set_init_conc(**dct[k.ID])
    u.AS.set_init_conc(concentrations=df.iloc[:6])
    c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    tss = [v for v in dct['C1_tss'].values() if v>0]
    u.C1.set_init_solubles(**c1s)
    u.C1.set_init_sludge_solids(**c1x)
    u.C1.set_init_TSS(tss)

#%%

def create_components():
    cmps = pc.create_asm2d_cmps(False)
    # cmps.X_I.i_N = 0.0600327162
    # cmps.X_I.i_P = 0.01
    # cmps.S_I.i_N = 0.0600327162
    # cmps.S_I.i_P = 0.01
    # cmps.X_PAO.i_N = 0.07
    # cmps.refresh_constants()
    cmps.compile()
    return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_asm2d = qs.get_thermo()
    
    # breakpoint()
    
    dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    # acetic_acid = qs.WasteStream('organic_carbon', S_A= 79.4, units='kg/hr', T=Temp)
    
    # inf_kwargs = inf_kwargs or default_inf_kwargs
    dom_ww.set_flow_by_concentration(Q_domestic, 
                                     concentrations=domestic_ww, 
                                     units=('m3/d', 'mg/L'))
    
    ind_ww.set_flow_by_concentration(Q_brewery, 
                                      concentrations= bm2, 
                                      units=('m3/d', 'mg/L'))
    
    effluent = qs.WasteStream('effluent', T=Temp)
    WAS = qs.WasteStream('WAS', T=Temp)
    RAS = qs.WasteStream('RAS', T=Temp)
           
    asm_kwargs = asm_kwargs or default_asm2d_kwargs
    asm2d = pc.ASM2d()
    
    eff_GT = qs.WasteStream(ID = 'effluent_GT')
    eff_MT = qs.WasteStream(ID = 'effluent_MT')
    eff_DU = qs.WasteStream(ID = 'effluent_DU')
    # split_A2 = qs.WasteStream(ID = 'split_A2')
    # split_A3 = qs.WasteStream(ID = 'split_A3')
    
    # Sludge flow rate and solids removal efficiency based on WERF report
    PC = su.PrimaryClarifier('PC', ins=[dom_ww, 
                                         ind_ww, 
                                        'reject'],
                              outs=('effluent_PC', 'sludge_PC'), isdynamic=True,
                              init_with='WasteStream', thermo=thermo_asm2d,
                              sludge_flow_rate=280, 
                              solids_removal_efficiency=0.6)
    
    # S1 = su.Splitter('S1', ins = PC-0, outs = [split_A2, split_A3], split= 0.8)

# <<<<<<< Updated upstream
#     # A1 = su.CSTR('A1', ins = [acetic_acid, RAS],  V_max=249, aeration=None, suspended_growth_model=asm2d)   
#     # A2 = su.CSTR('A2', ins = [split_A2, A1-0], V_max=2313, aeration=None, suspended_growth_model=asm2d)    
#     # A3 = su.CSTR('A3', [split_A3, A2-0, 'internal_recycle'], V_max=2633, aeration=None, suspended_growth_model=asm2d)    
#     # A4 = su.CSTR('A4', A3-0, V_max=2633, aeration=None, suspended_growth_model=asm2d)    
#     # O1 = su.CSTR('O1', A4-0, V_max=V_ae, aeration=2, DO_ID='S_O2', suspended_growth_model=asm2d)    
#     # O2 = su.CSTR('O2', O1-0, outs=['', 2-A3], split=[0.205, 0.795],V_max=V_ae, aeration=2, DO_ID='S_O2', suspended_growth_model=asm2d)    
# =======
    # A1 = su.CSTR('A1', ins = [acetic_acid, RAS],  V_max=249, aeration=None, suspended_growth_model=asm2d)   
    # A2 = su.CSTR('A2', ins = [split_A2, A1-0], V_max=2313, aeration=None, suspended_growth_model=asm2d)    
    # A3 = su.CSTR('A3', [split_A3, A2-0, 'internal_recycle'], V_max=2633, aeration=None, suspended_growth_model=asm2d)    
    # A4 = su.CSTR('A4', A3-0, V_max=2633, aeration=None, suspended_growth_model=asm2d)    
    # O1 = su.CSTR('O1', A4-0, V_max=V_ae, aeration=2, DO_ID='S_O2', suspended_growth_model=asm2d)    
    # O2 = su.CSTR('O2', O1-0, outs=['', 2-A3], split=[0.232, 0.768],V_max=V_ae, aeration=2, DO_ID='S_O2', suspended_growth_model=asm2d)    
# >>>>>>> Stashed changes
    
    # O1.temperature = 20
    # O1.blower_efficiency = 0.7
    # O1.inlet_pressure_loss = 1
    # O1.diffuser_pressure_loss = 7
    # O1.sludge_disposal_cost = 375
    # O1.unit_electricity_costs = 0.0577
    
    AS = su.PFR('AS', ins=[PC-0, RAS
                            # , acetic_acid
                           ], outs='treated', 
                N_tanks_in_series=6,
                V_tanks=[249,2313,2633,2633,V_ae,V_ae],
                influent_fractions=[[0, .8, .2, 0,0,0], 
                                    [1,0,0,0,0,0]
                                    # ,[1,0,0,0,0,0]
                                    ],
                internal_recycles=[(5,2,151416.47)],
                DO_setpoints=[0]*4+[2.0]*2, DO_ID='S_O2', kLa=[],
                suspended_growth_model=asm2d)
    
    AS.temperature = 20
    AS.blower_efficiency = 0.7
    AS.inlet_pressure_loss = 1
    AS.diffuser_pressure_loss = 7
    AS.sludge_disposal_cost = 375
    AS.unit_electricity_costs = 0.0577
    AS.natural_gas_price = 0.0041

    # C1 = su.FlatBottomCircularClarifier('C1', O2-0, [effluent, RAS, WAS],
    C1 = su.FlatBottomCircularClarifier('C1', AS-0, [effluent, 1-AS, WAS],
                                        underflow=Q_ras, wastage=Q_was, 
                                        surface_area= 1580, height=3.66, N_layer=10, 
                                        feed_layer=5, thermo = thermo_asm2d)
    
    # thickener_perc and TSS_removal_perc based on WERF report
    GT = su.Thickener('GT', PC-1, ['sludge_GT', eff_GT], 
                       thickener_perc= 7, TSS_removal_perc=92)
    
    # thickener_perc and TSS_removal_perc based on WERF report
    MT = su.Thickener('MT', WAS, ['sludge_MT', eff_MT],
                       thickener_perc= 6, TSS_removal_perc=98)
    
    M1 = su.Mixer('M1', [MT-0, GT-0])
    
    cmps_adm1 = qs.processes.create_adm1_p_extension_cmps()
    cmps_adm1.X_PAO.i_N = 0.07
    cmps_adm1.X_PAO.i_P = 0.02
    cmps_adm1.refresh_constants()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1_p_extension()
    
    J1 = su.ASM2dtomADM1('J1', upstream= [M1-0], thermo=thermo_adm1, isdynamic=True, 
                         adm1_model=adm1, asm2d_model=asm2d)
    
    # Volume of AD based on WERF report (G1 configuration)
  
    DG = su.AnaerobicCSTR(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
                          model=adm1, thermo = thermo_adm1, V_liq=3596, V_gas=416)
    DG.algebraic_h2 = True
    
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_asm2d = qs.get_thermo()
    
    J2 = su.mADM1toASM2d('J2', upstream = DG-1, thermo=thermo_asm2d, isdynamic=True, 
                         adm1_model=adm1, asm2d_model=asm2d)

    # thickener_perc and TSS_removal_perc based on WERF report
    DU = su.Centrifuge('DU', ins = J2.outs[0], outs = ['sludge_DU', eff_DU], thermo = thermo_asm2d,
                        thickener_perc=23, TSS_removal_perc=95)
    
    DU.DOC_f = 0.45
    DU.MCF = 0.8
    DU.k = 0.06 # temperate and boreal (dry)
    DU.pl = 30
    DU.elec_EF = 0.675 # the value for Minnesota
    DU.CH4_EF_st = 0.0075
    DU.N2O_EF_st = 0.016
    DU.CH4_EF_dis = 0.009
    DU.N2O_EF_dis = 0.005
    
    M2 = su.Mixer('M2', ins=[eff_GT, eff_MT, eff_DU], outs=2-PC)
    
    # DU.DOC_f = 0.45
    # DU.MCF = 0.8
    # DU.k = 0.06 # temperate and boreal (dry)
    # DU.pl = 30
    # DU.elec_EF = 0.675 # the value for Minnesota
    # DU.CH4_EF_st = 0.0075
    # DU.N2O_EF_st = 0.016
    # DU.CH4_EF_dis = 0.009
    # DU.N2O_EF_dis = 0.005

    sys = qs.System('G1_WERF', path=(PC, GT, AS, C1, MT, M1, J1, DG, J2, DU, M2), 
                    recycle = [RAS, M2-0,])
    # sys = qs.System('G1_WERF', path=(PC, S1, A1, A2, A3, A4, O1, O2, C1, GT, MT, M1, J1, DG, J2, DU, M2), 
    #                 recycle = [M2-0, 1-A3, RAS])
                    # recycle = [eff_GT, eff_MT, eff_DU, RAS])
    sys.set_tolerance(rmol=1e-6)
    sys.maxiter = 500
    return sys
#%%

@time_printer
def run(sys, t, method=None, **kwargs):

    path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
    batch_init(sys, path, 
               sheet='G1')
    
    sys.set_dynamic_tracker(*sys.products)
    
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')    
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
    sys.diagram()
    # return sys
    
if __name__ == '__main__':
    t = 15
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'

    sys = create_system()
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit
    
    # sys._setup()
    # sys.converge()
    run(sys, t, method=method)
    
    sys.diagram()

    act_units = [u.ID for u in sys.units if u.ID.startswith('AS')]
    
    # srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days\n')
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------

    # ACCOUNTING FOR AERATION COSTS (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    aeration_units = [u for u in sys.units if u.ID.startswith('AS')]
    
    # To get ASM2d components
    cmps = create_components()
    qs.set_thermo(cmps)
    cmps = qs.get_thermo().chemicals
    
    asm2d = aeration_units[0].suspended_growth_model
    idx_DO = cmps.index('S_O2')
    DO_sat = 8
    kLa = {}
    Q_air = {}
    aer = pc.DiffusedAeration('aer', 'S_O2', KLa=240, DOsat=DO_sat, V=V_ae)
    
    for u in aeration_units:
        i = 0
        for v in u.DO_setpoints:
            print(f'i = {i}')
            if v > 0:
                conc = u.state.iloc[i].to_numpy()
                print(f'conc[{i}] = {conc}')
                OTR = - asm2d.production_rates_eval(conc)[idx_DO]
                print(f'OTR = {OTR}')
                DO = conc[idx_DO]
                print(f'DO = {DO}')
                kLa[i] = aer.KLa = OTR / (DO_sat - DO)
                Q_air[i] = aer.Q_air
            else:
                kLa[i] = 0
                Q_air[i] = 0
            print(f'kLa = {kLa}')
            print(f'Q_air = {Q_air}')
            i += 1
    
    airflow = sum(Q_air.values())/24/60 # in m3/min
    
    power_blower = get_P_blower(q_air=airflow)
    print(f'Required aeration power at steady state is {power_blower:.2f} kW\n')
    
    # ACCOUNTING FOR PUMPING POWER (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    act_power_units = [u.ID for u in sys.units if \
    isinstance(u, (su.PrimaryClarifier, su.FlatBottomCircularClarifier, su.Thickener, su.Centrifuge))]
        
    required_pumping_power = get_power_utility(sys, active_unit_IDs=act_power_units)
    print(f'Required pumping (and other equipment) power at steady state is {required_pumping_power:.2f} kW\n')
    
    disposed_sludge = (fs.sludge_DU, )
    sludge_disposal_costs = get_cost_sludge_disposal(sludge = disposed_sludge, unit_weight_disposal_cost = 375)
    print(f'Sludge disposal cost = {sludge_disposal_costs} USD/day\n')

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
    
    # influent_ST = (fs.effluent_PC, fs.organic_carbon)
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
    
    
    GHG_sludge_disposal = get_GHG_emissions_sludge_disposal(sludge = disposed_sludge)
    print(f'CH4 emissions due to sludge_disposal equals {GHG_sludge_disposal} kg-CH4/day\n')
    
# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
    
    normalized_energy = get_normalized_energy(system = sys, aeration_power = power_blower, \
                                              pumping_power = required_pumping_power, \
                                              miscellaneous_power = 0)
    print(f'Normalized energy = {normalized_energy} kWh/m3\n')
    
    daily_operational_cost = get_daily_operational_cost(system = sys, aeration_power = power_blower, 
                            pumping_power = required_pumping_power, miscellaneous_power = 0, 
                            sludge_disposal_cost = sludge_disposal_costs)
    print(f'Daily operational costs = {daily_operational_cost} USD/day\n')
    
    total_daily_operational_cost = get_total_operational_cost(q_air = airflow, sludge = disposed_sludge, 
                                              system = sys, active_unit_IDs= act_power_units)
    print(f'Total daily operational costs = {total_daily_operational_cost} USD/day\n')
    
    CO2_eq_WRRF = get_CO2_eq_WRRF(system = sys, GHG_treatment = GHG_ST, GHG_discharge = GHG_discharge, \
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
    
    mass_degradable_sludge = np.array([slg.composite("C", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
                  degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
    mass_degradable_sludge = np.sum(mass_degradable_sludge)
    print(f'Mass degradable sludge = {mass_degradable_sludge} kg/day\n')
    
    mass_N_sludge = np.array([slg.composite("N", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
                  degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
    mass_N_sludge = np.sum(mass_N_sludge)
    print(f'Mass N sludge = {mass_N_sludge} kg/day\n')
    
    fig, axis = fs.effluent.scope.plot_time_series(('S_A', 'S_F', 'X_S', 'S_NH4', 'X_I', 'S_I', 'S_N2')) 
    fig
    
    # _dstate of units should be close to zero to imply steady state 
    
# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------

    G1_model = mod(sys)
    
    param = G1_model.parameter
    
    # T : float
    #     Air temperature [degree Celsius].
    dist_T = shape.Uniform(lower = 15, upper = 25)
    @param(name = 'Air temperature', 
            element = fu.AS, 
            kind = 'coupled', # try 'isolated' too
            units = 'degree celsius', 
            baseline = fu.AS.temperature, # mean of upper and lower limits
            distribution = dist_T)
    def set_air_temperature(i):
        fu.AS.temperature = i
        
    # efficiency : float
    #     Blower efficiency. Default is 0.7. 
    dist_blower_efficiency = shape.Uniform(lower = 0.6, upper = 0.8)
    @param(name = 'Blower efficiency', 
            element = fu.AS, 
            kind = 'coupled', # try 'isolated' too
            units = 'unitless', 
            baseline = fu.AS.blower_efficiency, # mean of upper and lower limits
            distribution = dist_blower_efficiency)
    def set_blower_efficiency(i):
        fu.AS.blower_efficiency = i
        
    # P_inlet_loss : float
    #     Head loss at inlet [kPa]. The default is 1 kPa. 
    dist_P_inlet_loss = shape.Uniform(lower = 0.9, upper = 1.1)
    @param(name = 'Head loss at inlet', 
            element = fu.AS, 
            kind = 'coupled', # try 'isolated' too
            units = 'kPa', 
            baseline = fu.AS.inlet_pressure_loss, # mean of upper and lower limits
            distribution = dist_P_inlet_loss)
    def set_inlet_pressure_loss(i):
        fu.AS.inlet_pressure_loss = i
        
    # P_diffuser_loss : float
    #     Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
    dist_P_diffuser_loss = shape.Uniform(lower = 6.3, upper = 7.7)
    @param(name = 'Head loss due to piping and diffuser', 
            element = fu.AS, 
            kind = 'coupled', # try 'isolated' too
            units = 'kPa', 
            baseline = fu.AS.diffuser_pressure_loss, # mean of upper and lower limits
            distribution = dist_P_diffuser_loss)
    def set_diffuser_pressure_loss(i):
        fu.AS.diffuser_pressure_loss = i
        
    dist_disposal_costs = shape.Uniform(lower = 337.5, upper = 412.5)
    @param(name = 'Sludge disposal unit cost', 
            element = fu.AS, 
            kind = 'coupled', # try 'isolated' too
            units = 'USD/ton', 
            baseline = fu.AS.sludge_disposal_cost, # mean of upper and lower limits
            distribution = dist_disposal_costs)
    def set_sludge_disposal_cost(i):
        fu.AS.sludge_disposal_cost = i
        
    dist_electricity_cost = shape.Triangle(lower= 0.0577, midpoint= 0.075, upper=0.0878)
    @param(name = 'Electricity unit cost', 
            element = fu.AS, 
            kind = 'coupled', # try 'isolated' too
            units = 'USD/kWh', 
            baseline = fu.AS.unit_electricity_costs, # mean of upper and lower limits
            distribution = dist_electricity_cost)
    def set_electricity_cost(i):
        fu.AS.unit_electricity_costs = i
        
    
    dist_natural_gas_price =  shape.Triangle(lower= 0.00359, midpoint=0.0041, upper=0.00598) # Source: EIA for year 2023
    @param(name = 'Price of natural gas', 
              element = fu.AS, 
              kind = 'coupled', # try 'isolated' too
              units = 'USD/m3', 
              baseline = fu.AS.natural_gas_price, 
              distribution =  dist_natural_gas_price)
    def set_dist_natural_gas_price(i):
          fu.AS.natural_gas_price = i
    
# # -----------------------------------------------------------------------------
# # -----------------------------------------------------------------------------

    # DOC_f : float, optional
    #     fraction of DOC that can decompose (fraction).
    dist_DOC_f = shape.Uniform(lower = 0.4, upper = 0.5)
    @param(name = 'Decomposable fraction of DOC', 
            element = fu.DU, 
            kind = 'coupled', # try 'isolated' too
            units = 'unitless', 
            baseline = fu.DU.DOC_f, # Deafault = 0.38 (Table 2.4 V5 Ch2, IPCC 2019)
            distribution =  dist_DOC_f)
    def set_DOC_f(i):
        fu.DU.DOC_f = i
    
    # MCF : float, optional
    #     CH4 correction factor for aerobic decomposition in the year of deposition (fraction). The default is 0.8.
    dist_MCF = shape.Uniform(lower = 0.4, upper = 1)
    @param(name = 'MCF', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = 'unitless', 
              baseline = fu.DU.MCF, # Table 3.1, V5 Ch3, IPCC 2019
              distribution =  dist_MCF)
    def set_MCF(i):
          fu.DU.MCF = i
    
    # k : TYPE, optional
    #     Methane generation rate (k). The default is 0.06. (1/year)
    dist_k = shape.Triangle(lower= 0.05, midpoint= 0.06, upper=0.2) # Table 3.3, V5 Ch3, IPCC 2019
    @param(name = 'Methane generation rate', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = '1/year', 
              baseline = fu.DU.k, 
              distribution =  dist_k)
    def set_k(i):
          fu.DU.k = i
    
    # pl : float, optional
    #     The project lifetime over which methane emissions would be calculated. (years)
    #     The default is 30 years.
    dist_pl = shape.Uniform(lower = 30, upper = 50) # in years
    @param(name = 'Project lifetime', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = 'years', 
              baseline = fu.DU.pl, 
              distribution =  dist_pl)
    def set_pl(i):
          fu.DU.pl = i
    
    # CO2_EF : float
    #     The emission factor used to calculate tier-2 CO2 emissions due to electricity consumption. 
    #     The default is 0.675 kg-CO2-Eq/kWh.
    dist_elec_EF = shape.Uniform(lower = 0.244, upper = 0.814) # in years
    @param(name = 'EF for electrcity', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = 'kg-CO2-Eq/kWh', 
              baseline = fu.DU.elec_EF, 
              distribution =  dist_elec_EF)
    def set_elec_EF(i):
          fu.DU.elec_EF = i
         
    # CH4_EF_sc :  float, optional.
    #     The emission factor used to calculate methane emissions in secondary treatment. The default is 0.0075 kg CH4/ kg rCOD.
    dist_CH4_EF_st = shape.Triangle(lower= 0.00075, midpoint=0.0075, upper=0.0225) # in kg CH4/kg CO2 r
    @param(name = 'EF for CH4 during treatment', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = 'kg CH4/kg COD removed', 
              baseline = fu.DU.CH4_EF_st, 
              distribution =  dist_CH4_EF_st)
    def set_CH4_EF_st(i):
          fu.DU.CH4_EF_st = i
          
    # N2O_EF_sc : float, optional
    #     The emission factor used to calculate nitrous oxide emissions in secondary treatment. The default is 0.016 kg N2O-N/ kg N.
    dist_N2O_EF_st = shape.Triangle(lower= 0.013, midpoint=0.016, upper=0.045) # in kg CH4/kg CO2 r
    @param(name = 'EF for N2O during treatment', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = 'kg N/kg influent N', 
              baseline = fu.DU.N2O_EF_st, 
              distribution =  dist_N2O_EF_st)
    def set_N2O_EF_st(i):
          fu.DU.N2O_EF_st = i
        
    # CH4_EF_discharge :  float, optional.
    #     The emission factor used to calculate methane emissions in discharge. The default is 0.009 kg CH4/ kg effluent COD. 
    dist_CH4_EF_dis = shape.Triangle(lower= 0.001, midpoint=0.009, upper=0.015) # in kg CH4/kg CO2 r
    @param(name = 'EF for CH4 during discharge', 
              element = fu.DU, 
              kind = 'coupled', # try 'isolated' too
              units = 'kg CH4/kg effluent COD', 
              baseline = fu.DU.CH4_EF_dis, 
              distribution =  dist_CH4_EF_dis)
    def set_CH4_EF_dis(i):
          fu.DU.CH4_EF_dis = i
          
    # N2O_EF_discharge : float, optional
    #     The emission factor used to calculate nitrous oxide emissions in discharge. The default is 0.005 kg N2O-N/ kg effluent N.
    # dist_N2O_EF_dis =  shape.Triangle(lower= 0.0005, midpoint=0.005, upper=0.075) # in kg CH4/kg CO2 r
    # @param(name = 'EF for N2O during discharge', 
    #           element = fu.DU, 
    #           kind = 'coupled', # try 'isolated' too
    #           units = 'kg N/kg effluent N', 
    #           baseline = fu.DU.N2O_EF_dis, 
    #           distribution =  dist_N2O_EF_dis)
    # def set_N2O_EF_dis(i):
    #       fu.DU.N2O_EF_dis = i
          
#     # QUESTIONS FOR METRIC 
    
#     # Should the element be 'fu.DU' -> if param and metric should be linked by a common element
#     # Should the element be the function 'get_cost_sludge_disposal', if so I want the 
#     # function 'daily_operational_cost' also to run for every value of get_cost_sludge_disposal
    
#     disposed_sludge = (fs.sludge_DU, )
    
    metric = G1_model.metric
    
    @metric(name='Total operational cost', units='USD/m3', element= fu.AS)
    def operational_cost():
        return get_total_operational_cost(q_air = airflow, sludge = disposed_sludge, 
                                          system = sys, active_unit_IDs= act_power_units, 
                                          
                                              # uncertain parameters 
                                              
                                              T= fu.AS.temperature, 
                                              efficiency= fu.AS.blower_efficiency, 
                                              P_inlet_loss= fu.AS.inlet_pressure_loss, 
                                              P_diffuser_loss= fu.AS.diffuser_pressure_loss, 
                                              
                                              unit_weight_disposal_cost = fu.AS.sludge_disposal_cost,
                                              unit_electricity_cost = fu.AS.unit_electricity_costs)
    
    @metric(name='Aeration cost', units='USD/m3', element= fu.AS)
    def aeration_cost():
        return get_aeration_cost(q_air = airflow, system = sys, 
                                 
                              # uncertain parameters 
                                 
                              T= fu.AS.temperature, 
                              efficiency= fu.AS.blower_efficiency, 
                              P_inlet_loss= fu.AS.inlet_pressure_loss, 
                              P_diffuser_loss= fu.AS.diffuser_pressure_loss, 
                              unit_electricity_cost = fu.AS.unit_electricity_costs)
    
    @metric(name='Pumping cost', units='USD/m3', element= fu.AS)
    def pumping_cost():
        return get_pumping_cost(system = sys, active_unit_IDs= act_power_units,  
                                
                                # uncertain parameter
                              unit_electricity_cost = fu.AS.unit_electricity_costs)
    
    @metric(name='Sludge disposal cost', units='USD/m3', element= fu.AS)
    def sludge_disposal_costs():
        return get_sludge_disposal_costs(sludge = disposed_sludge,  # sludge disposal costs 
                                      system= sys, 
                                      # uncertain parameters 
                                      unit_weight_disposal_cost = fu.AS.sludge_disposal_cost # sludge disposal costs 
                                      )

    @metric(name='Total GHG emissions', units='kg CO2 eq./m3', element= fu.DU)
    def GHG_emissions():
        return get_total_CO2_eq(system = sys, q_air = airflow, 
                                influent_sc = influent_ST, effluent_sc = effluent_ST, effluent_sys = effluent, 
                                active_unit_IDs=act_power_units, sludge= disposed_sludge, 
                            
                                # uncertain parameters 
                                efficiency= fu.AS.blower_efficiency,  
                                P_inlet_loss= fu.AS.inlet_pressure_loss, 
                                P_diffuser_loss= fu.AS.diffuser_pressure_loss,
                                 
                                CH4_EF_sc = fu.DU.CH4_EF_st, 
                                N2O_EF_sc = fu.DU.N2O_EF_st, 
                                CH4_EF_discharge = fu.DU.CH4_EF_dis,
                                # N2O_EF_discharge = fu.DU.N2O_EF_dis,
                                N2O_EF_discharge = 0.005,
                                 
                                CO2_EF= fu.DU.elec_EF, 
                                DOC_f = fu.DU.DOC_f, 
                                MCF = fu.DU.MCF, 
                                k = fu.DU.k, 
                                pl= fu.DU.pl)
    
    @metric(name='CH4 CO2_eq_treatment', units='kg CO2 eq./m3', element= fu.DU)
    def CH4_CO2_eq_treatment():
        return get_CH4_CO2_eq_treatment(system = sys, influent_sc =influent_ST, effluent_sc = effluent_ST,
                                 
                                        # active_unit_IDs=act_power_units, sludge=disposed_sludge, 
                                 
                                  # uncertain parameters 
                                CH4_EF_sc = fu.DU.CH4_EF_st)
        
    @metric(name='N2O CO2_eq_treatment', units='kg CO2 eq./m3', element= fu.DU)
    def N2O_CO2_eq_treatment():
        return get_N2O_CO2_eq_treatment(system = sys, influent_sc =influent_ST,  
                                  N2O_EF_sc=fu.DU.N2O_EF_st)
        
    @metric(name='CH4 CO2_eq_discharge', units='kg CO2 eq./m3', element= fu.DU)
    def CH4_CO2_eq_discharge():
        return get_CH4_CO2_eq_discharge(system = sys, effluent_sys = effluent, 
                                  CH4_EF_discharge=fu.DU.CH4_EF_dis)
    
    @metric(name='N2O CO2_eq_discharge', units='kg CO2 eq./m3', element= fu.DU)
    def N2O_CO2_eq_discharge():
        return get_N2O_CO2_eq_discharge(system = sys, effluent_sys = effluent,
                                        N2O_EF_discharge  = 0.005)
    
    @metric(name='CH4 emitted_during_pl', units='kg CO2 eq./m3', element= fu.DU)
    def CH4_emitted_during_pl():
        return get_CH4_emitted_during_pl(system = sys, sludge=disposed_sludge,
                                          DOC_f = fu.DU.DOC_f, 
                                          MCF = fu.DU.MCF, 
                                          k = fu.DU.k, 
                                          pl= fu.DU.pl)
    
    @metric(name='CH4 emitted_after_pl', units='kg CO2 eq./m3', element= fu.DU)
    def CH4_emitted_after_pl():
        return get_CH4_emitted_after_pl(system = sys, sludge=disposed_sludge,  
                              # uncertain parameters 
                              DOC_f = fu.DU.DOC_f, 
                              MCF = fu.DU.MCF, 
                              k = fu.DU.k, 
                              pl= fu.DU.pl
                              )
    
    @metric(name='CO2 eq_electricity', units='kg CO2 eq./m3', element= fu.DU)
    def CO2_eq_electricity():
        return get_CO2_eq_electricity(system = sys, q_air = airflow, active_unit_IDs=act_power_units, 
                          # uncertain parameters 
                          efficiency= fu.AS.blower_efficiency,  
                          P_inlet_loss= fu.AS.inlet_pressure_loss, 
                          P_diffuser_loss= fu.AS.diffuser_pressure_loss,
                          CO2_EF= fu.DU.elec_EF)
    						
    np.random.seed(3221) # Setting the seed ensures you getting the same sample
    
    #seed fixes the sample
    samples = G1_model.sample(N=10000, rule='L')
    G1_model.load_samples(samples)
    
    # sys.isdynamic = False
    # B1_model.evaluate()
    
    all_outputs = []
    for sample in samples:
        for p, v in zip(G1_model.parameters, sample):
            p.setter(v)
    
        out = [m() for m in G1_model.metrics]
        all_outputs.append(out)
    
    G1_model.table.iloc[:,-12:] = all_outputs
    
    G1_model.table

    import pandas as pd
    file_path = "/Users/saumitrarai/Desktop/Research/MS_research/metro/results/G1_comprehensive_uncertainty_analysis.xlsx"  
    sheet_name = "bm2_new"  # Change this to your desired sheet name
    with pd.ExcelWriter(file_path, engine='openpyxl', mode='a') as writer:
        G1_model.table.to_excel(writer, sheet_name=sheet_name)
    
    fig, ax = qs.stats.plot_uncertainties(G1_model, x_axis=G1_model.metrics[0], y_axis=G1_model.metrics[1],
                                      kind='kde-kde', center_kws={'fill': True})
    fig

    r_df, p_df = qs.stats.get_correlations(G1_model, kind='Spearman')
    print(r_df, p_df)
    fig, ax = qs.stats.plot_correlations(r_df)
    fig