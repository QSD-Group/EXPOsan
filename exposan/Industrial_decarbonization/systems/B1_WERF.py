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

import os, qsdsan as qs
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
    get_GHG_emissions_sec_treatment,
    get_GHG_emissions_discharge,
    get_GHG_emissions_electricity,
    get_GHG_emissions_sludge_disposal,
    get_CO2_eq_WRRF,
    )

from exposan.bsm1 import data_path
from exposan.metab import UASB, flex_rhos_adm1
# from . import data_path
folder = ospath.dirname(__file__)


__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm1_kwargs', 
    # 'default_ad_init_conds', 
    'steady_state_ad_init_conds',
    'domestic_ww', 'default_init_conds',
    'Q_domestic', 'Q_brewery', 'Q_ras', 'Q_was', 'Temp', 'V_ae',
    )



# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_domestic = 38000      # influent flowrate [m3/d] = 10 MGD (WERF report)
Q_brewery = 785         # industrial wastewater [m3/d] flowrate is determined based on a constant 
                        # domestic/industrial flowrate ratio across data from Metropolitan WRRF

Temp = 273.15+20 # temperature [K]

V_ae = 3785/6 # (1 MG) aerated zone tank volume [m3] (WERF report)

Q_was =  0.0175*(39472) # [m3/day] 0.1 MGD (WERF report)
Q_ras =  0.67*(39472)  # [m3/day] 67% of effluent from PC (WERF report)

biomass_IDs = ('X_BH', 'X_BA')

ammonia_ID = 'S_NH'

# Based on WERF report
domestic_ww = {
  'S_I': 10, 
  'S_S': 68,  
  'X_I': 28, 
  'X_S': 242,
  'X_BH': 5, 
  'X_BA': 5, 
  'S_NH': 25,  
  'S_ND': 12,
  'S_ALK':7*12,
    }

baseline_plus_two_after_METAB = {'S_I': 61.30441099087758, 
                          'S_S': 80.87361657821778, 
                          'X_I': 104.0874623234365, 
                          'X_S': 1220.0156256560024, 
                          'X_BH': 0.0, 
                          'X_BA': 0.0, 
                          'X_P': 440.62152373361744, 
                          'S_O': 0.0, 
                          'S_NO': 0.0, 
                          'S_NH': 132.1989432564827, 
                          'S_ND': 0.24140946486420325, 
                          'X_ND': 98.37021792729705, 
                          'S_ALK': 324.4643951852935, 
                          'S_N2': 0.0}

# default_asm2d_kwargs = dict(iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
#             iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
#             iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
#             f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
#             Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
#             Y_A=0.24, f_XI_AUT=0.1,
#             K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
#             mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
#             K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
#             q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
#             b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
#             K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
#             K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
#             mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
#             k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, 
#             # path=os.path.join(data_path, '_asm2d.tsv'),
#             )

default_asm1_kwargs = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
    eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )

default_init_conds = {
        'S_S':7,
        'X_I':1000,
        'X_S':100,
        'X_H':500,
        'X_AUT':100,
        #'X_P':100,
        'S_O':2,
        'S_NO3':20,
        'S_NH4':2,
        'S_ALK':7*1
    }

C_mw = 12
N_mw = 14

default_ad_init_conds = {
      'S_su': 0.0124*1e3,
      'S_aa': 0.0055*1e3,
      'S_fa': 0.1074*1e3,
      'S_va': 0.0123*1e3,
      'S_bu': 0.0140*1e3,
      'S_pro': 0.0176*1e3,
      'S_ac': 0.0893*1e3,
      'S_h2': 2.5055e-7*1e3,
      'S_ch4': 0.0555*1e3,
      'S_IC': 0.0951*C_mw*1e3,
      'S_IN': 0.0945*N_mw*1e3,
      'S_I': 0.1309*1e3,
      'X_ch': 0.0205*1e3,
      'X_pr': 0.0842*1e3,
      'X_li': 0.0436*1e3,
      'X_su': 0.3122*1e3,
      'X_aa': 0.9317*1e3,
      'X_fa': 0.3384*1e3,
      'X_c4': 0.3258*1e3,
      'X_pro': 0.1011*1e3,
      'X_ac': 0.6772*1e3,
      'X_h2': 0.2848*1e3,
      'X_I': 17.2162*1e3
      }

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
      'S_I': 2.325e+02,
      'X_c': 2.171e+02,
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
       'S_cat': 0.000e+00, 
       'S_an': 4.772e+00
      }


def batch_init(sys, path, sheet):
    # isa = isinstance
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    # u.DG.set_init_conc(**default_ad_init_conds)
    u.DG.set_init_conc(**steady_state_ad_init_conds)
    for k in sys.units:
        if sheet.startswith('B1'):
            if k.ID.startswith('O'): k.set_init_conc(**dct[k.ID])
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
     asm1_cmps = pc.create_asm1_cmps(False)
     cmps = qs.Components([*asm1_cmps])
     cmps.compile()
     return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={}, lifetime=100, discount_rate=0.1,
                  aeration_processes=()):
    # flowsheet = flowsheet or qs.Flowsheet('bsm1')
    # qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_asm1 = qs.get_thermo()
    
    dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    
    # inf_kwargs = inf_kwargs or default_inf_kwargs
    dom_ww.set_flow_by_concentration(Q_domestic, 
                                     concentrations=domestic_ww, 
                                     units=('m3/d', 'mg/L'))
    
    ind_ww.set_flow_by_concentration(Q_brewery, 
                                      concentrations= baseline_plus_two_after_METAB, 
                                      units=('m3/d', 'mg/L'))
    
    effluent = qs.WasteStream('effluent', T=Temp)
    WAS = qs.WasteStream('WAS', T=Temp)
    RAS = qs.WasteStream('RAS', T=Temp)

    # fuel = WasteStream('nat_gas', phase='g', S_CH4=1000)
    # air = WasteStream('air', phase='g', S_O2=210, S_N2=780, S_H2=10)
    
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
    #     aer1 = pc.DiffusedAeration('aer1', 'S_O', KLa=98.97, DOsat=8.0, V=V_ae)
    #     aer2 = pc.DiffusedAeration('aer2', 'S_O', KLa=317, DOsat=8.0, V=V_ae)
    #     aer3 = pc.DiffusedAeration('aer3', 'S_O', KLa=165.46, DOsat=8.0, V=V_ae)
    #     aer4 = pc.DiffusedAeration('aer4', 'S_O', KLa=143.24, DOsat=8.0, V=V_ae)
    #     aer5 = pc.DiffusedAeration('aer5', 'S_O', KLa=130.08, DOsat=8.0, V=V_ae)
    #     aer6 = pc.DiffusedAeration('aer6', 'S_O', KLa=113.93, DOsat=8.0, V=V_ae)
       
    # asm_kwargs = asm_kwargs or default_asm2d_kwargs
    # asm2d = pc.ASM2d(iP_SF=0.005, iP_XS=0.005, iP_XI=0.005, iN_BM=0.1, iTSS_XI=0.72)
    
    asm_kwargs = asm_kwargs or default_asm1_kwargs
    asm1 = pc.ASM1(**asm_kwargs)
    
    # asm1 = pc.ASM1(iTSS_XI=3)
    
    eff_GT = qs.WasteStream(ID = 'effluent_GT')
    eff_MT = qs.WasteStream(ID = 'effluent_MT')
    eff_DU = qs.WasteStream(ID = 'effluent_DU')
    
    PC = su.PrimaryClarifier(ID='PC', ins = [dom_ww, ind_ww, 
                                             eff_GT, eff_MT, eff_DU], 
                      outs = ['sludge_PC', 'effluent_PC'], isdynamic=True, 
                      init_with='WasteStream', thickener_perc= 0.09, 
                      TSS_removal_perc=70, thermo = thermo_asm1)
    
    kwargs_O = dict(V_max=V_ae, aeration=2, DO_ID='S_O', suspended_growth_model=asm1)
    
    # O1 = su.CSTR('O1', ins = [RAS, PC-1], V_max=V_ae, aeration=aer1, DO_ID='S_O', suspended_growth_model=asm1)    
    # O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2, DO_ID='S_O', suspended_growth_model=asm1)    
    # O3 = su.CSTR('O3', O2-0, V_max=V_ae, aeration=aer3, DO_ID='S_O', suspended_growth_model=asm1)    
    # O4 = su.CSTR('O4', O3-0, V_max=V_ae, aeration=aer4, DO_ID='S_O', suspended_growth_model=asm1)    
    # O5 = su.CSTR('O5', O4-0, V_max=V_ae, aeration=aer5, DO_ID='S_O', suspended_growth_model=asm1)    
    # O6 = su.CSTR('O6', O5-0, V_max=V_ae, aeration=aer6, DO_ID='S_O', suspended_growth_model=asm1)  
    
    O1 = su.CSTR('O1', [RAS, PC-1],  **kwargs_O)    
    O2 = su.CSTR('O2', O1-0, **kwargs_O)    
    O3 = su.CSTR('O3', O2-0, **kwargs_O)    
    O4 = su.CSTR('O4', O3-0, **kwargs_O)    
    O5 = su.CSTR('O5', O4-0, **kwargs_O)    
    O6 = su.CSTR('O6', O5-0, **kwargs_O)    
    
    C1 = su.FlatBottomCircularClarifier('C1', O6-0, [effluent, RAS, WAS],
                                        underflow=Q_ras, wastage=Q_was, 
                                        surface_area= 1580, height=3.66, N_layer=10, 
                                        feed_layer=5, thermo = thermo_asm1)
    
    # thickener_perc and TSS_removal_perc based on WERF report
    GT = su.Thickener('GT', PC-0, ['sludge_GT', eff_GT], 
                       thickener_perc= 7, TSS_removal_perc=92)
    
    # thickener_perc and TSS_removal_perc based on WERF report
    MT = su.Thickener('MT', WAS, ['sludge_MT', eff_MT],
                       thickener_perc= 6, TSS_removal_perc=98)
    
    Mixer = su.Mixer('Mixer', [MT-0, GT-0])
    
    cmps_adm1 = qs.processes.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    # adm1 = qs.processes.ADM1()
    adm1 = qs.processes.ADM1(flex_rate_function=flex_rhos_adm1)

    J1 = su.ASMtoADM('J1', upstream= [Mixer-0], thermo=thermo_adm1, isdynamic=True, 
                      adm1_model=adm1)
    
    # Volume of AD based on WERF report
    DG = UASB(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
                          model=adm1, thermo = thermo_adm1, V_liq = 3217, V_gas = 321.7, 
                          pH_ctrl = 6.8, fraction_retain=0) 
    
    # DG = su.AnaerobicCSTR(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
    #                       model=adm1, thermo = thermo_adm1, V_liq = 3217, V_gas = 321.7) 
    
    J2 = su.ADMtoASM('J2', upstream = DG-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)

    # thickener_perc and TSS_removal_perc based on WERF report
    DU = su.Centrifuge('DU', ins = J2.outs[0], outs = ['sludge_DU', eff_DU], thermo = thermo_asm1,
                       thickener_perc=23, TSS_removal_perc=95)
    
    DU.sludge_disposal_cost = 440

    sys = qs.System('B1_ASM1', 
                    path=(PC, O1, O2, O3, O4, O5, O6, C1, GT, MT, Mixer, J1, DG, J2, DU), 
                    recycle = [eff_GT, eff_MT, eff_DU, RAS])
        
    
    
    return sys
#%%

@time_printer
def run(t, method=None, **kwargs):
# <<<<<<< Updated upstream
    sys = create_system()
    path = "/Users/saumitrarai/Dropbox/Personal/Gates_WWTP modelling/initial_conditions.xlsx"
    batch_init(sys, path, sheet='t=10')
# =======
#     sys = create_system()    
    
#     batch_init(sys, "/Users/saumitrarai/Dropbox/Personal/Gates_WWTP modelling/initial_conditions.xlsx", sheet='t=10')
# >>>>>>> Stashed changes
    
    
    path = ospath.join(folder, "data/initial_conditions.xlsx")    
    batch_init(sys, path, sheet='B1')
        
    # RAS = sys.flowsheet.stream.RAS
    # C1 = sys.flowsheet.unit.C1
    # sys.set_dynamic_tracker(RAS, C1)
    sys.set_dynamic_tracker(*sys.products)
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
    # sys._setup()
    # sys.converge() 
    
    sys.diagram()
    return sys

#%%
if __name__ == '__main__':
    t = 1
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
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit
    
    

    # act_units = [u.ID for u in sys.units if isinstance(u, su.FlatBottomCircularClarifier) or u.ID.startswith('O')]
    act_units = [u.ID for u in sys.units if u.ID.startswith('O')]
    
    srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days\n')
    
    # Aeration tanks (AT)
    numerator_tank = V_ae*(fs.ws1.get_TSS() + fs.ws3.get_TSS() + fs.ws5.get_TSS() + fs.ws7.get_TSS() + \
                      fs.ws9.get_TSS() + fs.ws11.get_TSS())  # mg/L * m3
    
    # Secondary clarifier (SC)
    N_layer = 10 
    mlss_reactor = fu.C1._state[-N_layer:].mean()
    V_reactor = fu.C1.V_settle
    numerator_sc = mlss_reactor*V_reactor # mg/L * m3
    
    denominator = (fs.WAS.get_TSS()*fs.WAS.F_vol) + (fs.effluent.get_TSS()*fs.effluent.F_vol) # mg/L * m3/hr
    
    SRT_sc = numerator_sc/denominator # in hr
    
    MCRT = (numerator_tank + numerator_sc)/denominator # in hr
    
    SRT = numerator_tank/denominator # in hr
    
    print(f'SRT_sc  = {SRT_sc/24} day\n')
    
    print(f'MCRT  = {MCRT/24} day\n')
    
    print(f'SRT  = {SRT/24} day\n')
    
# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------

    # ACCOUNTING FOR AERATION COSTS (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    aeration_units = [u for u in sys.units if u.ID.startswith('O')]
    
    # To get ASM1 components
    cmps = create_components()
    qs.set_thermo(cmps)
    cmps = qs.get_thermo().chemicals
    
    asm1 = aeration_units[0].suspended_growth_model
    idx_DO = cmps.index('S_O')
    DO_sat = 8
    kLa = {}
    Q_air = {}
    aer = pc.DiffusedAeration('aer', 'S_O', KLa=240, DOsat=DO_sat, V=V_ae)
    for u in aeration_units:
        conc = u._state[:-1]
        OTR = - asm1.production_rates_eval(conc)[idx_DO]
        DO = conc[idx_DO]
        kLa[u.ID] = aer.KLa = OTR / (DO_sat - DO)
        Q_air[u.ID] = aer.Q_air
    airflow = sum(Q_air.values())/24/60 # in m3/min
    
    print(f'Airflow = {airflow} m3/min\n')
    
    power_blower = get_P_blower(q_air=airflow)
    print(f'Required aeration power at steady state is {power_blower:.2f} kW\n')
    
    # ACCOUNTING FOR PUMPING POWER (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    act_power_units = [u.ID for u in sys.units if \
    isinstance(u, (su.PrimaryClarifier, su.FlatBottomCircularClarifier, su.Thickener, su.Centrifuge))]
        
    required_pumping_power = get_power_utility(sys, active_unit_IDs=act_power_units)
    print(f'Required pumping (and other equipment) power at steady state is {required_pumping_power:.2f} kW\n')
    
    disposed_sludge = (fs.sludge_DU, )
    sludge_disposal_costs = get_cost_sludge_disposal(sludge = disposed_sludge, unit_weight_disposal_cost = fu.DU.sludge_disposal_cost)
    print(f'Sludge disposal cost = {sludge_disposal_costs} USD/day\n')

# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
    
#     influent_ST = (fs.effluent_PC,)
#     effluent_ST = (fs.effluent, fs.WAS)
#     GHG_ST = get_GHG_emissions_sec_treatment(influent = influent_ST, effluent = effluent_ST)
#     print(f'CH4 and N2O emissions during secondary treatment equals {GHG_ST[0]} kg CH4/day and\
#       {GHG_ST[1]} kg N2O-N/day respectively\n')
     
#     effluent = (fs.effluent, )
#     GHG_discharge = get_GHG_emissions_discharge(effluent = effluent)
#     print(f'CH4 and N2O emissions at discharge equals {GHG_discharge[0]} kg CH4/day and \
#       {GHG_discharge[1]} kg N2O-N/day respectively\n')
     
#     GHG_electricity = get_GHG_emissions_electricity(system=sys, power_blower=power_blower, 
#                                                     power_pump=required_pumping_power)
#     print(f'CO2 emissions due to electricity consumption equals {GHG_electricity} kg-CO2-eq/day\n')
    
#     sludge = (fs.sludge_DU, )
#     GHG_sludge_disposal = get_GHG_emissions_sludge_disposal(sludge = sludge)
#     print(f'CH4 emissions due to sludge_disposal equals {GHG_sludge_disposal} kg-CH4/day\n')
    
# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
    
#     normalized_energy = get_normalized_energy(system = sys, aeration_power = power_blower, \
#                                               pumping_power = required_pumping_power, \
#                                               miscellaneous_power = 0)
#     print(f'Normalized energy = {normalized_energy} kWh/m3\n')
    
#     daily_operational_cost = get_daily_operational_cost(aeration_power = power_blower, 
#                             pumping_power = required_pumping_power, miscellaneous_power = 0, 
#                             sludge_disposal_cost = sludge_disposal_costs)
#     print(f'Daily operational costs = {daily_operational_cost} USD/day\n')
    
#     CO2_eq_WRRF = get_CO2_eq_WRRF(system = sys, GHG_treatment = GHG_ST, GHG_discharge = GHG_discharge, \
#                                     GHG_electricity = GHG_electricity, GHG_sludge_disposal = GHG_sludge_disposal)
#     print(f'GHG emissions = {CO2_eq_WRRF} kg CO2 eq/m3\n')
    
#     fig, axis = fs.effluent.scope.plot_time_series(('S_S', 'X_S', 'S_NH', 'S_NO', 'X_I', 'S_I', 'S_N2')) 
#     fig
    
#     # _dstate of units should be close to zero to imply steady state 
    
# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
    
#     # B1_model = mod(sys)
    
#     # param = B1_model.parameter
#     # dist_disposal_costs = shape.Uniform(lower = 100, upper = 650)
    
#     # @param(name = 'Landfill disposal costs', 
#     #         element = fu.DU, 
#     #         kind = 'coupled', # try 'isolated' too
#     #         units = 'USD/ton', 
#     #         baseline = 375, 
#     #         distribution = dist_disposal_costs)
#     # def set_sludge_disposal_cost(i):
#     #     fu.DU.sludge_disposal_cost = i
        
#     # # QUESTIONS FOR METRIC 
    
#     # # Should the element be 'fu.DU' -> if param and metric should be linked by a common element
#     # # Should the element be the function 'get_cost_sludge_disposal', if so I want the 
#     # # function 'daily_operational_cost' also to run for every value of get_cost_sludge_disposal
    
#     # disposed_sludge = (fs.sludge_DU, )
    
#     # metric = B1_model.metric
#     # @metric(name='Operational cost', units='USD/day', element= fu.DU)
#     # def get_cost_sludge_dis():
#     #     return get_cost_sludge_disposal(sludge = disposed_sludge, 
#     #                                     unit_weight_disposal_cost = fu.DU.sludge_disposal_cost)
    
#     # import numpy as np
#     # np.random.seed(3221) # setting the seed ensures you getting the same sample

#     # samples = B1_model.sample(N=10, rule='L')
#     # B1_model.load_samples(samples)
        
#     # B1_model.evaluate()
    
#     # B1_model.table
