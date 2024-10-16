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
import numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su
from qsdsan.utils import ospath, time_printer, load_data, get_SRT, \
    get_oxygen_heterotrophs, get_oxygen_autotrophs, get_airflow, get_P_blower
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

Q_domestic = 605000       # influent flowrate [m3/d]
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

# baseline_minus_one_effluent_metab = {
#     'S_O2': 0, 
#     'S_N2': 0,
#     'S_NH4': 101.44,    
#     'S_NO3': 0, 
#     'S_PO4': 0,
#     'S_F': 6.85,
#     'S_A': 59.6,
#     'S_I': 53.92,
#     'S_ALK': 0,
#     'X_I': 88.47,
#     'X_S': 622.611,
#     'X_H': 0,
#     'X_PAO': 0, 
#     'X_PP': 0,
#     'X_PHA': 0,
#     'X_AUT': 0, 
#     'X_MeOH': 0,
#     'X_MeP': 0, 
#     }
bm1_metab = {
    'S_NH4': 101.5638526804,
    'S_F': 34.122,
    'S_A': 32.315,
    'S_I': 53.92,
    'S_ALK': 287.2,
    'X_I': 137.423,
    'X_S': 573.658
    }

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
        if sheet == 'ss':
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
     asm2d_cmps = pc.create_asm2d_cmps()
     # asm2d_cmps.X_S.f_BOD5_COD = 0.54
     CO2 = qs.Component.from_chemical('S_CO2', search_ID='CO2', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
     CH4 = qs.Component.from_chemical('S_CH4', search_ID='CH4', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
     H2 = qs.Component.from_chemical('S_H2', search_ID='H2', particle_size='Dissolved gas', degradability='Undegradable', organic=False)
     cmps1 = qs.Components.load_default()
     ash = cmps1.X_Ig_ISS.copy('ash')
     cmps = qs.Components([*asm2d_cmps, CO2, CH4, H2, ash])
     cmps.compile()
     return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    flowsheet = flowsheet or qs.Flowsheet('bsm1')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    domestic_wastewater = qs.WasteStream('domestic_wastewater', T=Temp)
    industrial_wastewater = qs.WasteStream('industrial_wastewater', T=Temp)
    #inf_kwargs = inf_kwargs or default_inf_kwargs
    domestic_wastewater.set_flow_by_concentration(Q_domestic, 
                                                  concentrations=domestic_ww, 
                                                  units=('m3/d', 'mg/L'))
    # industrial_wastewater.set_flow_by_concentration(Q_brewery, 
    #                                                 concentrations=bm1_metab, 
    #                                                 units=('m3/d', 'mg/L'))
    effluent = qs.WasteStream('effluent', T=Temp)
    WAS = qs.WasteStream('WAS', T=Temp)
    RAS = qs.WasteStream('RAS', T=Temp)
    # fuel = WasteStream('nat_gas', phase='g', S_CH4=1000)
    # air = WasteStream('air', phase='g', S_O2=210, S_N2=780, S_H2=10)
    
    asm2d = pc.ASM2d(iP_SF=0.005, iP_XS=0.005, iP_XI=0.005, iN_BM=0.1, iTSS_XI=0.72)
    
    effluent_TC1 = qs.WasteStream(ID = 'effluent_TC1')
    
    PC = su.Thickener(ID='PC', ins = [domestic_wastewater, industrial_wastewater, effluent_TC1], 
                      outs = ['sludge_PC', 'effluent_PC'], isdynamic=True, 
                      init_with='WasteStream', thickener_perc= 0.09, 
                      TSS_removal_perc=69)
    
    effluent_TC2 = qs.WasteStream(ID = 'effluent_TC2')
    effluent_DU1 = qs.WasteStream(ID = 'effluent_DU1')
    
    O1 = su.CSTR('O1', RAS, V_max=V_ae, aeration=2, DO_ID='S_O2', 
                 suspended_growth_model=asm2d)
    
    A1 = su.CSTR('A1', O1-0, V_max=V_an,
                  aeration=None, suspended_growth_model=asm2d)
    
    A2 = su.CSTR('A2', [A1-0, PC-1, effluent_TC2, effluent_DU1], V_max=V_an,
                  aeration=None, suspended_growth_model=asm2d)
    
    A3 = su.CSTR('A3', A2-0, V_max=V_an,
                  aeration=None, suspended_growth_model=asm2d)
    
    O2 = su.CSTR('O2', A3-0, V_max=V_ae, aeration=2, 
                 DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O3 = su.CSTR('O3', O2-0, V_max=V_ae, aeration=2,
                  DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O4 = su.CSTR('O4', O3-0, V_max=V_ae, aeration=2,
                  DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O5 = su.CSTR('O5', O4-0, V_max=V_ae, aeration=2,
                  DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O6 = su.CSTR('O6', O5-0, V_max=V_ae, aeration=2,
                  DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O7 = su.CSTR('O7', O6-0, V_max=V_ae, aeration=2,
                  DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O8 = su.CSTR('O8', O7-0, V_max=V_ae, aeration=2,
                  DO_ID='S_O2', suspended_growth_model=asm2d)
    
    O9 = su.CSTR('O9', O8-0, ['treated'], V_max=V_ae, aeration=2, 
                  DO_ID='S_O2', suspended_growth_model=asm2d)

    C1 = su.FlatBottomCircularClarifier('C1', O9-0, [effluent, RAS, WAS],
                                        underflow=Q_ras, wastage=Q_was, 
                                        surface_area=24226, height=6, N_layer=10, 
                                        feed_layer=5, X_threshold=3000, v_max=474, 
                                        v_max_practical=250, rh=5.76e-4, rp=2.86e-3, 
                                        fns=2.28e-3)
    
    TC1 = su.Thickener(ID='TC1', ins = PC-0, outs = ['sludge_TC1', effluent_TC1], 
                       thickener_perc= 6.1, TSS_removal_perc=97.14)
    
    TC2 = su.Thickener(ID='TC2', ins = WAS, outs = ['sludge_TC2', effluent_TC2],
                       thickener_perc=4.3, TSS_removal_perc=95.54)
    
    DU1 = su.Centrifuge(ID='DU1', ins = [TC1-0, TC2-0], outs = ['sludge_DU1', effluent_DU1],
                            thickener_perc=27, TSS_removal_perc=96.29)

    # IC = su.Incinerator(ID ='IC', ins = [DU1-0, air, fuel], outs=(), thermo=None, isdynamic=False, 
    #                     init_with='WasteStream', F_BM_default=None, process_efficiency=0.90, 
    #                     calorific_value_sludge= 12000, calorific_value_fuel=50000, carbon_di_oxide_ID = 'CO2',
    #                     ash_component_ID = 'ash')

    sys = qs.System('metro_ASM2d', 
                    # path=(O1, A1, A2, A3, O2, O3, O4, O5, O6, O7, O8, O9, PC, C1, TC1, TC2, DU1, IC), 
                    path=(PC, A2, A3, O2, O3, O4, O5, O6, O7, O8, O9, C1, O1, A1, TC1, TC2, DU1), 
                    recycle = [TC1-1, TC2-1, DU1-1, A1-0])
    
    
    return sys
#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    sys = create_system()
    # for u in sys.units:
    #     if u.ID in ('O1', 'A1', 'A2', 'A3', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9'):
    #         u.set_init_conc(**default_init_conds)
    # batch_init(sys, "/Users/saumitrarai/Desktop/Important files_QSDsan/initial_conditions_ASM2d.xlsx", sheet='t=10')
    path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")
    batch_init(sys, path, 
               sheet='ss')
               # sheet='t=10')
        
    # RAS = sys.flowsheet.stream.RAS
    # C1 = sys.flowsheet.unit.C1
    # sys.set_dynamic_tracker(RAS, C1)
    sys.set_dynamic_tracker(*sys.products)
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    sys.diagram()
    return sys
    
if __name__ == '__main__':
    t = 30
    # t = 15
    t_step = 1
    # method = 'RK45'
    # method = 'RK23' 
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    system = run(t, t_step, method=method)

    # act_units = [u.ID for u in system.units if isinstance(u, (su.CSTR, su.FlatBottomCircularClarifier))]
    # fs = system.flowsheet.stream
    # fu = system.flowsheet.unit
    # srt = get_SRT(system, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # # oxygen_heterotroph = get_oxygen_heterotrophs(system, influent=fu.A2.ins[1:], 
    # #                                              f_d=0.1, b_H=0.4, SRT=srt, Y_H=0.625)
    # inf = qs.WasteStream('mix')
    # inf.mix_from(fu.A2.ins[1:])
    # eff_sCOD = fs.effluent.composite('COD', particle_size='s')
    # oxygen_heterotroph = get_oxygen_heterotrophs(inf.F_vol*24, inf.COD, eff_sCOD, 
    #                                              SRT=srt, f_d=0.1, Y_H=0.625, b_H=0.4)
    # print(f'Required oxygen for heterotrophs at steady state is {oxygen_heterotroph:.2f} kg/day')
    
    # # oxygen_autotroph = get_oxygen_autotrophs(system, influent=fu.A2.ins[1:], 
    # #                                          SRT=srt, ammonia_component_ID=ammonia_ID)
    # inf_TKN = inf.TN - inf.composite('N', subgroup=('S_NO3',))
    # oxygen_autotroph = get_oxygen_autotrophs(inf.F_vol*24, inf.COD, eff_sCOD, 
    #                                          inf_TKN, SRT=srt, Y_H=0.625, U_AUT=1.0, 
    #                                          b_H=0.4, b_AUT=0.15, f_d=0.1, SF_DO=1.25)
    # print(f'Required oxygen for autotrophs at steady state is {oxygen_autotroph:.2f} kg/day')
    
    # airflow = get_airflow(oxygen_heterotrophs = oxygen_heterotroph, oxygen_autotrophs = oxygen_autotroph, oxygen_transfer_efficiency = 15) # in m3/min
    # print(f'Required airflow at steady state is {airflow:.2f} m3/min')
    
    # power_blower = get_P_blower(q_air=airflow)
    # print(f'Required power at steady state is {power_blower:.2f} kW')

    # aeration_energy_demand = power_blower / sum([s.F_vol for s in system.feeds])
    # print(f'Energy demand for aeration is {aeration_energy_demand:.3f} kWh/m3')

    # #fig, axis = fs.RAS.scope.plot_time_series(( 'S_NH')) 
    # # fig, axis = fs.RAS.scope.plot_time_series(('S_S','S_NH')) 
    # fig, axis = fs.effluent.scope.plot_time_series(('S_F', 'S_A', 'S_NH4', 'S_NO3')) 
    # fig
