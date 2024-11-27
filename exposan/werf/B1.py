# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import qsdsan as qs
from qsdsan import (
    # WasteStream,
    processes as pc,
    sanunits as su,
    )
from qsdsan.utils import ospath, time_printer, load_data, get_SRT
from exposan.werf import default_ad_init, default_as_init, default_fctss_init

__all__ = ('create_b1_system',)

#%%
folder = ospath.dirname(__file__)
dfs = load_data(
    ospath.join(folder, 'data/initial_conditions.xlsx'), 
    sheet=None,
    )
asinit = dfs['rBOD']
fcinit = asinit.iloc[-1].to_dict()
adinit = dfs['adm'].iloc[0].to_dict()
# Default initial conditions
# dfs = load_data(ospath.join(folder, 'data/G1_init.xlsx'), sheet=None)
# asinit = dfs['asm'].iloc[1:]
# # asinit = dfs['asm_ss']
# # adinit = dfs['adm'].iloc[0].to_dict()
# c2init = dfs['settler'].to_dict('index')
# c2init['s'] = {k:v for k,v in c2init['s'].items() if v>0}
# c2init['x'] = {k:v for k,v in c2init['x'].items() if v>0}
# c2init['tss'] = [v for k,v in c2init['tss'].items() if v>0]

# adinit = {
#     'S_su': 0.011649409987233706,
#     'S_aa': 0.005224687709048591,
#     'S_fa': 0.09143444853388716,
#     'S_va': 0.010538555495939485,
#     'S_bu': 0.01381502642962188,
#     'S_pro': 0.015627785163003477,
#     'S_ac': 0.055098426143215494,
#     'S_ch4': 0.2233336530302807,
#     'S_IC': 1.0427915551564186,
#     'S_IN': 1.669357973671328,
#     'S_IP': 0.28964898102287245,
#     'S_I': 0.016450418739744545,
#     'X_ch': 0.06047631655400267,
#     'X_pr': 0.06388757232327166,
#     'X_li': 0.08593503930110351,
#     'X_su': 0.9503620480475024,
#     'X_aa': 0.7484248883486445,
#     'X_fa': 0.7586744810131059,
#     'X_c4': 0.31971827836675687,
#     'X_pro': 0.15224043553443883,
#     'X_ac': 1.1263471919494425,
#     'X_h2': 0.4912284536782753,
#     'X_I': 19.166956666296826,
#     'X_PHA': 1.2849572506633961e-06,
#     'X_PP': 5.724052662511763e-15,
#     'X_PAO': 5.645526871794766e-07,
#     'S_K': 0.028198583386767055,
#     'S_Mg': 0.00022947578212156403,
#     'S_Ca': 0.0015405999465588393,
#     'X_CaCO3': 7.572335925505587e-10,
#     'X_struv': 0.5082333428046245,
#     'X_newb': 9.503417736993726e-10,
#     'X_ACP': 3.435265972797314,
#     'X_MgCO3': 7.572335925505587e-10,
#     'X_AlOH': 9.503440160735277e-10,
#     'X_AlPO4': 9.503420444805498e-10,
#     'X_FeOH': 9.503456387448865e-10,
#     'X_FePO4': 1.0453169886394833e-07,
#     'S_Na': 0.07991961153828457,
#     'S_Cl': 0.3903889983759098
#      }


MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

def create_b1_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('B1')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
    thermo_asm = qs.get_thermo()
    
    PC = su.PrimaryClarifier(
        'PC', ins=[rww, 'reject'], 
        outs=('PE', 'PS'),
        sludge_flow_rate=0.066*MGD2cmd,
        solids_removal_efficiency=0.6
        )
    
    GT = su.IdealClarifier(
        'GT', PC-1, outs=['', 'thickened_PS'],
        sludge_flow_rate=0.023*MGD2cmd,
        solids_removal_efficiency=0.9
        )
    
    n_zones = 6
    V_tot = 1.0 * MGD2cmd
    
    # ae_kwargs = dict(V_max=V_tot/n_zones, aeration=2.0, DO_ID='S_O2', 
    #                  suspended_growth_model=asm, gas_stripping=True)
        
    # O1 = su.CSTR('O1', [PC-0, 'RAS'], **ae_kwargs)
    # O2 = su.CSTR('O2', O1-0, **ae_kwargs)
    # O3 = su.CSTR('O3', O2-0, **ae_kwargs)
    # O4 = su.CSTR('O4', O3-0, **ae_kwargs)
    # O5 = su.CSTR('O5', O4-0, **ae_kwargs)
    # O6 = su.CSTR('O6', O5-0, **ae_kwargs)
    
    ASR = su.PFR(
        'ASR', ins=[PC-0, 'RAS'], outs='treated',
        N_tanks_in_series=n_zones,
        V_tanks=[V_tot/n_zones]*n_zones,
        influent_fractions=[[1]+[0]*5]*2,
        internal_recycles=[], DO_ID='S_O2',
        kLa=[50]*n_zones, 
        DO_setpoints=[2.0]*n_zones,
        suspended_growth_model=asm,
        gas_stripping=True
        )
    
    FC = su.FlatBottomCircularClarifier(
        'FC', ins=ASR-0, outs=['SE', 1-ASR, 'WAS'],
        # 'FC', ins=O6-0, outs=['SE', 1-O1, 'WAS'],
        underflow=0.67*10*MGD2cmd, wastage=0.1*MGD2cmd,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=2.5e-3, fns=0.001, 
        maximum_nonsettleable_solids=20.0
        )
    
    # MT = su.Thickener(
    #     'MT', ins=FC-2, outs=['thickened_WAS', ''],
    #     thickener_perc=5, TSS_removal_perc=95,
    #     )
    # M1 = su.Mixer('M1', ins=[GT-1, MT-0])
    MT = su.IdealClarifier(
        'MT', FC-2, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.019*MGD2cmd,
        solids_removal_efficiency=0.95
        )
    M1 = su.Mixer('M1', ins=[GT-1, MT-1])
    
    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0)
    
    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR(
        'AD', ins=J1-0, outs=('biogas', 'digestate'), 
        V_liq=0.85*MGD2cmd, V_gas=0.085*MGD2cmd, 
        fixed_headspace_P=False, fraction_retain=0,
        T=T_ad, model=adm,
        pH_ctrl=7.0,
        )
    AD.algebraic_h2 = False
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)
    
    # DW = su.Centrifuge(
    #     'DW', ins=J2-0, outs=('cake', ''),
    #     thickener_perc=18, TSS_removal_perc=90,
    #     )
    # M2 = su.Mixer('M2', ins=[GT-0, MT-1, DW-1])    
    DW = su.IdealClarifier(
        'DW', J2-0, outs=('', 'cake'),
        sludge_flow_rate=0.0053*MGD2cmd,
        solids_removal_efficiency=0.9
        )
    M2 = su.Mixer('M2', ins=[GT-0, MT-0, DW-0])
    
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-PC)
    
    if default_init_conds:
        # ASR.set_init_conc(**default_as_init)
        # # for unit in (O1, O2, O3, O4, O5, O6):
        # #     unit.set_init_conc(**default_as_init)
        # FC.set_init_solubles(**default_as_init)
        # FC.set_init_sludge_solids(**default_as_init)
        ASR.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**fcinit)
        FC.set_init_sludge_solids(**fcinit)
        FC.set_init_TSS(default_fctss_init)
        # AD.set_init_conc(**default_ad_init)
        AD.set_init_conc(**adinit)
        # FC.set_init_solubles(**c2init['s'])
        # FC.set_init_sludge_solids(**c2init['x'])
        # FC.set_init_TSS(c2init['tss'])
    
    sys = qs.System(
        'B1', 
        path=(PC, GT, ASR, FC, MT, M1, J1, AD, J2, DW, M2, HD),
        # path=(PC, GT, O1, O2, O3, O4, O5, O6, FC, 
        #       MT, M1, J1, AD, J2, DW, M2, HD),
        recycle=(FC-1, HD-0)
        )

    sys.set_dynamic_tracker(FC-0, AD)

    return sys

#%%
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
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    

#%%
if __name__ == '__main__':
    sys = create_b1_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 50
    # t = 1
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    # biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[WAS],
    #               active_unit_IDs=('ASR',))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')