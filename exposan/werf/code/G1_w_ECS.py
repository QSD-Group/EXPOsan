# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>
    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import qsdsan as qs
import numpy as np
from qsdsan import WasteStream, processes as pc, sanunits as su, Model as mod
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
    get_eq_natural_gas_price, 
    get_GHG_emissions_sec_treatment,
    get_GHG_emissions_discharge,
    get_GHG_emissions_electricity,
    get_GHG_emissions_sludge_disposal,
    get_CO2_eq_WRRF,
    get_total_CO2_eq,
    get_eq_natural_gas_emission, 
    
    get_aeration_cost,
    get_pumping_cost,
    get_sludge_disposal_costs,
    get_carbon_add_cost,
    get_CH4_CO2_eq_treatment, 
    get_N2O_CO2_eq_treatment,
    get_CH4_CO2_eq_discharge,
    get_N2O_CO2_eq_discharge,
    get_CH4_emitted_during_pl, 
    get_CH4_emitted_after_pl,
    get_CO2_eq_electricity,
    get_carbon_add_CO2_eq
    )

# from exposan.werf import data_path

__all__ = ('create_g1_system',)

ID = 'G1'
#%%
dfs = load_data(
    ospath.join('/Users/saumitrarai/Downloads/initial_conditions_werf.xlsx'), 
    sheet=None,
    )
asinit = dfs[ID]
fcinit = asinit.iloc[-1].to_dict()
default_fctss_init = [10, 12, 12, 40, 500, 500, 500, 5e3, 1e4, 1.3e4]
adinit = dfs['adm'].loc[ID].to_dict()

MGD2cmd = 3785.412


n_zones = 6
V_tot = 4.7 * MGD2cmd
fr_V = [0.014, 0.13, 0.148, 0.148, 0.28, 0.28]
V_ae = V_tot*fr_V[4]

Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

ESAP_inc = 'included'
# ESAP_inc = 'not_included'

def create_g1_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    thermo_asm = qs.get_thermo()
    
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,S_Mg=7,
        )
    
    # rww = pc.create_masm2d_inf(
    #     'RWW', 10, 'MGD', T=Temp, 
    #     COD=358, NH4_N=40, PO4_P=5,
    #     fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2, S_Mg=7,
    #     )
    
    carb = WasteStream('carbon', T=Temp, units='kg/hr', S_A=40)

    PC = su.PrimaryClarifier(
        'PC', ins=[rww, 'reject'], 
        outs=('PE', 'PS'),
        isdynamic=True, 
        sludge_flow_rate=0.074*MGD2cmd,
        solids_removal_efficiency=0.6
        )
    
    GT = su.IdealClarifier(
        'GT', PC-1, outs=['', 'thickened_PS'],
        sludge_flow_rate=0.026*MGD2cmd,
        solids_removal_efficiency=0.9,
        )
    
    
    
    gstrip = True
    an_kwargs = dict(aeration=None, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=gstrip)
    ae_kwargs = dict(aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=gstrip)

    S1 = su.Splitter('S1', PC-0, split=0.8)
    
    A1 = su.CSTR('A1', ins=[carb, 'RAS'], V_max=V_tot*fr_V[0], **an_kwargs)
    A2 = su.CSTR('A2', [A1-0, S1-0], V_max=V_tot*fr_V[1], **an_kwargs)
    A3 = su.CSTR('A3', [A2-0, 'intr', S1-1], V_max=V_tot*fr_V[2], **an_kwargs)
    A4 = su.CSTR('A4', A3-0, V_max=V_tot*fr_V[3], **an_kwargs)
    O5 = su.CSTR('O5', A4-0, V_max=V_tot*fr_V[4], **ae_kwargs)
    O6 = su.CSTR('O6', O5-0, [1-A3, 'treated'], split=[40, 14],
                  V_max=V_tot*fr_V[5], **ae_kwargs)
    
    O6.temperature = 20
    O6.blower_efficiency = 0.7
    O6.inlet_pressure_loss = 1
    O6.diffuser_pressure_loss = 7
    O6.sludge_disposal_cost = 375
    O6.unit_electricity_costs = 0.0577
    O6.acetic_acid_price = 0.62 #USD/kg
    O6.nat_gas_price = 0.1617 #USD/m3

    # ASR = su.PFR(
    #     'ASR', ins=[PC-0, 'RAS', carb], outs='treated', 
    #     N_tanks_in_series=n_zones,
    #     V_tanks=[f*V_tot for f in fr_V],
    #     influent_fractions=[
    #         [0, 0.8, 0.2, 0,0,0],   # PC-0
    #         [1,0,0,0,0,0],          # RAS
    #         [1,0,0,0,0,0],          # carb
    #         ], 
    #     internal_recycles=[(5,2,40*MGD2cmd)], 
    #     kLa=[0,0,0,0,180,70],
    #     DO_setpoints=[0,0,0,0,2.0,2.0], DO_ID='S_O2',
    #     suspended_growth_model=asm,
    #     gas_stripping=True)
    
    FC = su.FlatBottomCircularClarifier(
        'FC', O6-1, ['SE', 1-A1, 'WAS'], 
        # 'FC', ASR-0, ['SE', 1-AS, 'WAS'],
        underflow=0.4*10*MGD2cmd, wastage=0.136*MGD2cmd,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01, 
        maximum_nonsettleable_solids=8.0
        )
    
    MT = su.IdealClarifier(
        'MT', FC-2, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.0335*MGD2cmd,
        solids_removal_efficiency=0.95,
        )
    M1 = su.Mixer('M1', ins=(GT-1, MT-1))
        
    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0)
    
    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR(
        'AD', ins=J1-0, outs=('biogas', 'digestate'), 
        V_liq=0.95*MGD2cmd, V_gas=0.11*MGD2cmd, 
        fixed_headspace_P=False, fraction_retain=0,
        T=T_ad, model=adm,
        pH_ctrl=7.0,
        )
    AD.algebraic_h2 = False
    
    if ESAP_inc == 'included':
    
        ESAP = qs.sanunits.ESAP('ESAP', ins = AD-1, component_ID_NH3 = 'S_IN', component_ID_P ='S_IP', component_ID_Mg = "S_Mg", 
                                outs=('recovered_product','loss','ESAP_effluent'),
          recovery={'S_IN':0.8,"S_IP":0.9,"S_Mg":0.8}, loss={'S_IN':0.06,"S_IP":0.05,"S_Mg":0.1}, isdynamic = True)
        
        J2 = su.ADM1ptomASM2d('J2', upstream=ESAP-2, thermo=thermo_asm, isdynamic=True, 
                              adm1_model=adm, asm2d_model=asm)
        qs.set_thermo(thermo_asm)
        
    elif ESAP_inc == 'not_included':
    
        J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                              adm1_model=adm, asm2d_model=asm)
    
    qs.set_thermo(thermo_asm)
    
    DU = su.PrimaryClarifier(
        'DU', J2-0, outs=['', 'cake'],
        sludge_flow_rate=0.0095*MGD2cmd,
        solids_removal_efficiency=0.9,
        )
    
    DU.DOC_f = 0.45
    DU.MCF = 0.8
    DU.k = 0.06 # temperate and boreal (dry)
    DU.pl = 30
    DU.elec_EF = 0.675 # the value for Minnesota
    DU.CH4_EF_st = 0.0075
    DU.N2O_EF_st = 0.004802 # Jason Ren
    DU.CH4_EF_dis = 0.009
    DU.N2O_EF_dis = 0.005
    DU.acetic_acid_EF = 1.669
    DU.nat_gas_EF = 0.454

    M2 = su.Mixer('M2', ins=(GT-0, MT-0, DU-0))
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-PC)

    if default_init_conds:
        asdct = asinit.to_dict('index')
        for i in (A1, A2, A3, A4, O5, O6):
            i.set_init_conc(**asdct[i.ID])
        # ASR.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**fcinit)
        FC.set_init_sludge_solids(**fcinit)
        FC.set_init_TSS(default_fctss_init)
        AD.set_init_conc(**adinit)
    
    if ESAP_inc == 'included':
        sys = qs.System(ID, 
            path=(PC, GT, S1, A1, A2, A3, A4, O5, O6, FC, 
                  MT, M1, J1, AD, ESAP, J2, DU, M2, HD),
            recycle=(O6-0, FC-1, HD-0))
    elif ESAP_inc == 'not_included':
        sys = qs.System(ID, 
            path=(PC, GT, S1, A1, A2, A3, A4, O5, O6, FC, 
                  MT, M1, J1, AD, J2, DU, M2, HD),
            recycle=(O6-0, FC-1, HD-0))
    
    sys.set_dynamic_tracker(FC-0, AD)
    
    return sys

# %%

@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    
    sys.simulate(
        # state_reset_hook='reset_cache',
        t_span=(0,t),
        # t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    
    return sys

#%%
if __name__ == '__main__':
    sys = create_g1_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 300
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit
    
    biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    srt = get_SRT(sys, biomass_IDs,
                  wastage=[fs.WAS],
                  active_unit_IDs=('A1', 'A2', 'A3', 'A4', 'O5', 'O6'))
                  # active_unit_IDs=('ASR'))
    if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))
    
    # ACCOUNTING FOR AERATION COSTS (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    aeration_units = [u for u in sys.units if u.ID.startswith('O')]
    
    # To get ASM2d components
    cmps = pc.create_masm2d_cmps(False)
    cmps.compile()
    qs.set_thermo(cmps)
    cmps = qs.get_thermo().chemicals
    
    mASM2d = aeration_units[0].suspended_growth_model
    idx_DO = cmps.index('S_O2')
    DO_sat = 8
    kLa = {}
    Q_air = {}
    aer = pc.DiffusedAeration('aer', 'S_O2', KLa=240, DOsat=DO_sat, V=V_ae)
    
    for u in aeration_units:
        conc = u._state[:-1]
        OTR = - mASM2d.production_rates_eval(conc)[idx_DO]
        DO = conc[idx_DO]
        kLa[u.ID] = aer.KLa = OTR / (DO_sat - DO)
        Q_air[u.ID] = aer.Q_air
    airflow = sum(Q_air.values())/24/60 # in m3/min
    
    
    power_blower = get_P_blower(q_air=airflow)
    print(f'Required aeration power at steady state is {power_blower:.2f} kW\n')
    
    # ACCOUNTING FOR PUMPING POWER (AND OTHER ELECTRICITY COSTS SUCH AS MOTOR IN CENTRIFUGE)
    act_power_units = [u.ID for u in sys.units if \
    isinstance(u, (su.PrimaryClarifier, su.FlatBottomCircularClarifier, su.IdealClarifier))]
        
    required_pumping_power = get_power_utility(sys, active_unit_IDs=act_power_units)
    print(f'Required pumping (and other equipment) power at steady state is {required_pumping_power:.2f} kW\n')
    
    disposed_sludge = (fs.cake, )
    sludge_disposal_costs = get_cost_sludge_disposal(sludge = disposed_sludge, unit_weight_disposal_cost = 375)
    print(f'Sludge disposal cost = {sludge_disposal_costs} USD/day\n')

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
    
    influent_ST = (fs.PE, fs.carbon)
    effluent_ST = (fs.SE, fs.WAS)
    
    GHG_ST = get_GHG_emissions_sec_treatment(influent = influent_ST, effluent = effluent_ST)
    print(f'CH4 and N2O emissions during secondary treatment equals {GHG_ST[0]} kg CH4/day and\
      {GHG_ST[1]} kg N2O-N/day respectively\n')
    
    effluent = (fs.SE, )
    GHG_discharge = get_GHG_emissions_discharge(effluent = effluent)
    print(f'CH4 and N2O emissions at discharge equals {GHG_discharge[0]} kg CH4/day and \
      {GHG_discharge[1]} kg N2O-N/day respectively\n')
     
    GHG_electricity = get_GHG_emissions_electricity(system=sys, power_blower=power_blower, 
                                                    power_pump=required_pumping_power)
    print(f'CO2 emissions due to electricity consumption equals {GHG_electricity} kg-CO2-eq/day\n')
    
    
    GHG_sludge_disposal = get_GHG_emissions_sludge_disposal(sludge = disposed_sludge)
    print(f'CH4 emissions due to sludge_disposal equals {GHG_sludge_disposal} kg-CH4/day\n')
    
# # # ------------------------------------------------------------------------------------------------------------------
# # # ------------------------------------------------------------------------------------------------------------------
    
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
    
    
    effluent_COD = fs.SE.COD
    print(f'Effluent COD = {effluent_COD} mg/L\n')
    
    effluent_TN = fs.SE.TN
    print(f'Effluent TN = {effluent_TN} mg/L\n')
    
    mass_degradable_sludge = np.array([slg.composite("C", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
                  degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
    mass_degradable_sludge = np.sum(mass_degradable_sludge)
    print(f'Mass degradable sludge = {mass_degradable_sludge} kg/day\n')
    
    mass_N_sludge = np.array([slg.composite("N", flow=True, exclude_gas=True, subgroup=None, particle_size=None,
                  degradability="b", organic=True, volatile=None, specification=None, unit="kg/day") for slg in disposed_sludge])
    mass_N_sludge = np.sum(mass_N_sludge)
    print(f'Mass N sludge = {mass_N_sludge} kg/day\n')
    
    fig, axis = fs.SE.scope.plot_time_series(('S_A', 'S_F', 'X_S', 'S_NH4', 'X_I', 'S_I', 'S_N2')) 
    fig
    
    # _dstate of units should be close to zero to imply steady state 
    
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------

#     G1_model = mod(sys)
    
#     param = G1_model.parameter
    
#     # T : float
#     #     Air temperature [degree Celsius].
#     dist_T = shape.Uniform(lower = 15, upper = 25)
#     @param(name = 'Air temperature', 
#             element = fu.O6, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'degree celsius', 
#             baseline = fu.O6.temperature, # mean of upper and lower limits
#             distribution = dist_T)
#     def set_air_temperature(i):
#         fu.O6.temperature = i
        
#     # efficiency : float
#     #     Blower efficiency. Default is 0.7. 
#     dist_blower_efficiency = shape.Uniform(lower = 0.6, upper = 0.8)
#     @param(name = 'Blower efficiency', 
#             element = fu.O6, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'unitless', 
#             baseline = fu.O6.blower_efficiency, # mean of upper and lower limits
#             distribution = dist_blower_efficiency)
#     def set_blower_efficiency(i):
#         fu.O6.blower_efficiency = i
        
#     # P_inlet_loss : float
#     #     Head loss at inlet [kPa]. The default is 1 kPa. 
#     dist_P_inlet_loss = shape.Uniform(lower = 0.9, upper = 1.1)
#     @param(name = 'Head loss at inlet', 
#             element = fu.O6, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'kPa', 
#             baseline = fu.O6.inlet_pressure_loss, # mean of upper and lower limits
#             distribution = dist_P_inlet_loss)
#     def set_inlet_pressure_loss(i):
#         fu.O6.inlet_pressure_loss = i
        
#     # P_diffuser_loss : float
#     #     Head loss due to piping and diffuser [kPa]. The default is 7 kPa. 
#     dist_P_diffuser_loss = shape.Uniform(lower = 6.3, upper = 7.7)
#     @param(name = 'Head loss due to piping and diffuser', 
#             element = fu.O6, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'kPa', 
#             baseline = fu.O6.diffuser_pressure_loss, # mean of upper and lower limits
#             distribution = dist_P_diffuser_loss)
#     def set_diffuser_pressure_loss(i):
#         fu.O6.diffuser_pressure_loss = i
        
#     dist_disposal_costs = shape.Uniform(lower = 337.5, upper = 412.5)
#     @param(name = 'Sludge disposal unit cost', 
#             element = fu.O6, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'USD/ton', 
#             baseline = fu.O6.sludge_disposal_cost, # mean of upper and lower limits
#             distribution = dist_disposal_costs)
#     def set_sludge_disposal_cost(i):
#         fu.O6.sludge_disposal_cost = i
        
#     dist_electricity_cost = shape.Triangle(lower= 0.0577, midpoint= 0.075, upper=0.0878)
#     @param(name = 'Electricity unit cost', 
#             element = fu.O6, 
#             kind = 'coupled', # try 'isolated' too
#             units = 'USD/kWh', 
#             baseline = fu.O6.unit_electricity_costs, # mean of upper and lower limits
#             distribution = dist_electricity_cost)
#     def set_electricity_cost(i):
#         fu.O6.unit_electricity_costs = i
        
#     # Shaped derived from monthly prices between June 2022 and June 2024. 
#     # 10th, and 90th percentile form the shape
#     dist_external_carbon_price =  shape.Triangle(lower= 0.47, midpoint= 0.62, upper=0.71)
#     @param(name = 'Price of external carbon', 
#               element = fu.O6, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'USD/kg', 
#               baseline = fu.O6.acetic_acid_price, 
#               distribution =  dist_external_carbon_price)
#     def set_dist_external_carbon_price(i):
#           fu.O6.acetic_acid_price = i
    
#     # Shaped derived from eia.gov monthly prices between June 2022 and June 2024. 
#     # 10th, 50th, and 90th percentile form the shape
#     dis_price_nat_gas = shape.Triangle(lower= 0.1269, midpoint= 0.1617, upper=0.3106)
#     @param(name = 'Price of natural gas', 
#               element = fu.O6, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'USD/m3', 
#               baseline = fu.O6.nat_gas_price, 
#               distribution =  dis_price_nat_gas)
#     def set_dist_price_nat_gas(i):
#           fu.O6.nat_gas_price = i
          
# # # # -----------------------------------------------------------------------------
# # # # -----------------------------------------------------------------------------

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
          
#     dist_external_carbon_EF = shape.Uniform(lower = 1.5021, upper = 1.8359)
#     @param(name = 'EF for external carbon addition', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'kg-CO2-Eq/kg acetic acid', 
#               baseline = fu.DU.acetic_acid_EF, 
#               distribution =  dist_external_carbon_EF)
#     def set_external_carbon_EF(i):
#           fu.DU.acetic_acid_EF = i
          
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
#     #     The emission factor used to calculate nitrous oxide emissions in secondary treatment. (used Song et al. 2024 (Jason Ren) SI to find all WRRFs with bioreactor as the filter.)
#     dist_N2O_EF_st = shape.Triangle(lower= 0.000048, midpoint=0.0048, upper=0.0247) # kg N2O-N/ kg N
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
          
#     dist_nat_gas_EF = shape.Triangle(lower = 0.409, midpoint = 0.454, upper = 0.499)
#     @param(name = 'EF for natural gas', 
#               element = fu.DU, 
#               kind = 'coupled', # try 'isolated' too
#               units = 'kg-CO2-Eq/m3 natural gas', 
#               baseline = fu.DU.nat_gas_EF, 
#               distribution =  dist_nat_gas_EF)
#     def set_nat_gas_EF(i):
#           fu.DU.nat_gas_EF = i
          
#     # N2O_EF_discharge : float, optional
#     #     The emission factor used to calculate nitrous oxide emissions in discharge. The default is 0.005 kg N2O-N/ kg effluent N.
#     # dist_N2O_EF_dis =  shape.Triangle(lower= 0.0005, midpoint=0.005, upper=0.075) # in kg CH4/kg CO2 r
#     # @param(name = 'EF for N2O during discharge', 
#     #           element = fu.DU, 
#     #           kind = 'coupled', # try 'isolated' too
#     #           units = 'kg N/kg effluent N', 
#     #           baseline = fu.DU.N2O_EF_dis, 
#     #           distribution =  dist_N2O_EF_dis)
#     # def set_N2O_EF_dis(i):
#     #       fu.DU.N2O_EF_dis = i
    
#     metric = G1_model.metric
    
#     @metric(name='Total operational cost', units='USD/m3', element= fu.O6)
#     def operational_cost():
#         return get_total_operational_cost(q_air = airflow, sludge = disposed_sludge, 
#                                           system = sys, active_unit_IDs= act_power_units, 
                                          
#                                               # uncertain parameters 
                                              
#                                               T= fu.O6.temperature, 
#                                               efficiency= fu.O6.blower_efficiency, 
#                                               P_inlet_loss= fu.O6.inlet_pressure_loss, 
#                                               P_diffuser_loss= fu.O6.diffuser_pressure_loss, 
                                              
#                                               unit_weight_disposal_cost = fu.O6.sludge_disposal_cost,
#                                               unit_electricity_cost = fu.O6.unit_electricity_costs,
#                                               organic_carbon = fs.carbon,
#                                               unit_cost_carbon_source = fu.O6.acetic_acid_price)
    
#     @metric(name='Aeration cost', units='USD/m3', element= fu.O6)
#     def aeration_cost():
#         return get_aeration_cost(q_air = airflow, system = sys, 
                                 
#                               # uncertain parameters 
                                 
#                               T= fu.O6.temperature, 
#                               efficiency= fu.O6.blower_efficiency, 
#                               P_inlet_loss= fu.O6.inlet_pressure_loss, 
#                               P_diffuser_loss= fu.O6.diffuser_pressure_loss, 
#                               unit_electricity_cost = fu.O6.unit_electricity_costs)
    
#     @metric(name='Pumping cost', units='USD/m3', element= fu.O6)
#     def pumping_cost():
#         return get_pumping_cost(system = sys, active_unit_IDs= act_power_units,  
                                
#                                 # uncertain parameter
#                               unit_electricity_cost = fu.O6.unit_electricity_costs)
    
#     @metric(name='Sludge disposal cost', units='USD/m3', element= fu.O6)
#     def sludge_disposal_costs():
#         return get_sludge_disposal_costs(sludge = disposed_sludge,  # sludge disposal costs 
#                                       system= sys, 
#                                       # uncertain parameters 
#                                       unit_weight_disposal_cost = fu.O6.sludge_disposal_cost # sludge disposal costs 
#                                       )
    
#     @metric(name='Carbon addition cost', units='USD/m3', element= fu.O6)
#     def carbon_add_cost():
#         return get_carbon_add_cost(organic_carbon = fs.carbon,
#                                       system= sys, 
#                                       # uncertain parameters 
#                                       unit_cost_carbon_source = fu.O6.acetic_acid_price
#                                       )
    
#     @metric(name='Natural gas cost', units='USD/m3', element= fu.O6)
#     def nat_gas_cost():
#         return get_eq_natural_gas_price(system = sys, gas = fs.biogas, 
#                                         mass_density = 0.68, 
#                                         natural_gas_price = fu.O6.nat_gas_price)


#     @metric(name='Total GHG emissions', units='kg CO2 eq./m3', element= fu.DU)
#     def GHG_emissions():
#         return get_total_CO2_eq(system = sys, q_air = airflow, 
#                                 influent_sc = influent_ST, effluent_sc = effluent_ST, effluent_sys = effluent, 
#                                 active_unit_IDs=act_power_units, sludge= disposed_sludge, 
                            
#                                 # uncertain parameters 
#                                 efficiency= fu.O6.blower_efficiency,  
#                                 P_inlet_loss= fu.O6.inlet_pressure_loss, 
#                                 P_diffuser_loss= fu.O6.diffuser_pressure_loss,
                                 
#                                 CH4_EF_sc = fu.DU.CH4_EF_st, 
#                                 N2O_EF_sc = fu.DU.N2O_EF_st, 
#                                 CH4_EF_discharge = fu.DU.CH4_EF_dis,
#                                 # N2O_EF_discharge = fu.DU.N2O_EF_dis,
#                                 N2O_EF_discharge = 0.005,
                                 
#                                 CO2_EF= fu.DU.elec_EF, 
#                                 DOC_f = fu.DU.DOC_f, 
#                                 MCF = fu.DU.MCF, 
#                                 k = fu.DU.k, 
#                                 pl= fu.DU.pl,
                                
#                                 organic_carbon = fs.carbon, 
#                                 EF_external_carbon = fu.DU.acetic_acid_EF
#                                 )
    
    
#     @metric(name='CH4 CO2_eq_treatment', units='kg CO2 eq./m3', element= fu.DU)
#     def CH4_CO2_eq_treatment():
#         return get_CH4_CO2_eq_treatment(system = sys, influent_sc =influent_ST, effluent_sc = effluent_ST,
                                 
#                                         # active_unit_IDs=act_power_units, sludge=disposed_sludge, 
                                 
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
    
#     @metric(name='CH4 emitted_during_pl', units='kg CO2 eq./m3', element= fu.DU)
#     def CH4_emitted_during_pl():
#         return get_CH4_emitted_during_pl(system = sys, sludge=disposed_sludge,
#                                           DOC_f = fu.DU.DOC_f, 
#                                           MCF = fu.DU.MCF, 
#                                           k = fu.DU.k, 
#                                           pl= fu.DU.pl)
    
#     @metric(name='CH4 emitted_after_pl', units='kg CO2 eq./m3', element= fu.DU)
#     def CH4_emitted_after_pl():
#         return get_CH4_emitted_after_pl(system = sys, sludge=disposed_sludge,  
#                               # uncertain parameters 
#                               DOC_f = fu.DU.DOC_f, 
#                               MCF = fu.DU.MCF, 
#                               k = fu.DU.k, 
#                               pl= fu.DU.pl
#                               )
#     @metric(name='CO2 eq_electricity', units='kg CO2 eq./m3', element= fu.DU)
#     def CO2_eq_electricity():
#         return get_CO2_eq_electricity(system = sys, q_air = airflow, active_unit_IDs=act_power_units, 
#                           # uncertain parameters 
#                           efficiency= fu.O6.blower_efficiency,  
#                           P_inlet_loss= fu.O6.inlet_pressure_loss, 
#                           P_diffuser_loss= fu.O6.diffuser_pressure_loss,
#                           CO2_EF= fu.DU.elec_EF)
    
#     @metric(name='CO2 eq Carbon addition', units='kg CO2 eq./m3', element= fu.DU)
#     def CO2_eq_C_addition():
#         return get_carbon_add_CO2_eq(organic_carbon = fs.carbon, system = sys, 
#                                       EF_external_carbon = fu.DU.acetic_acid_EF)
        
#     @metric(name='Natural gas emissions', units='kg CO2 eq/m3', element= fu.DU)
#     def nat_gas_emissions():
#         return get_eq_natural_gas_emission(system = sys, gas = fs.biogas, 
#                                             mass_density = 0.68, 
#                                             emission_factor = fu.DU.nat_gas_EF)
    						
#     np.random.seed(3221) # Setting the seed ensures you getting the same sample
    
#     #seed fixes the sample
#     samples = G1_model.sample(N=1000, rule='L')
#     G1_model.load_samples(samples)
    
#     # sys.isdynamic = False
#     # B1_model.evaluate()
    
#     all_outputs = []
#     for sample in samples:
#         for p, v in zip(G1_model.parameters, sample):
#             p.setter(v)
    
#         out = [m() for m in G1_model.metrics]
#         all_outputs.append(out)
    
#     G1_model.table.iloc[:,-16:] = all_outputs
    
#     G1_model.table
    
#     import pandas as pd
#     file_path = "/Users/saumitrarai/Desktop/Research/PhD_research/ARPA_E/Results/ECS_deployment_new.xlsx"
    
#     if ESAP_inc == 'included':
#         sheet_name = "G1_w_ESAP"  
#     elif ESAP_inc == 'not_included':
#         sheet_name = "G1" 
        
#     with pd.ExcelWriter(file_path, engine='openpyxl', mode='a') as writer:
#         G1_model.table.to_excel(writer, sheet_name=sheet_name)
    
#     fig, ax = qs.stats.plot_uncertainties(G1_model, x_axis=G1_model.metrics[0], y_axis=G1_model.metrics[1],
#                                       kind='kde-kde', center_kws={'fill': True})
#     fig

#     r_df, p_df = qs.stats.get_correlations(G1_model, kind='Spearman')
#     print(r_df, p_df)
#     fig, ax = qs.stats.plot_correlations(r_df)
#     fig
    
#     disposed_sludge = (fs.cake, )
#     ASP_inf = (fs.PE, )
#     ASP_eff = (fs.SE, fs.WAS, )
    
#     mass_N_sludge = np.array([slg.TN*slg.F_vol for slg in disposed_sludge])
#     mass_N_sludge = np.sum(mass_N_sludge)*24/1000
#     print(f'Mass N sludge = {mass_N_sludge} kg/day\n')
    
#     mass_N_ASP_inf = np.array([slg.TN*slg.F_vol for slg in ASP_inf])
#     mass_N_ASP_inf = np.sum(mass_N_ASP_inf)*24/1000
#     print(f'Mass N in ASP influent = {mass_N_ASP_inf} kg/day\n')
    
#     mass_N_ASP_eff = np.array([slg.TN*slg.F_vol for slg in ASP_eff])
#     mass_N_ASP_eff = np.sum(mass_N_ASP_eff)*24/1000
#     print(f'Mass N in ASP effluent = {mass_N_ASP_eff} kg/day\n')