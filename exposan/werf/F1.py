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

__all__ = ('create_f1_system',)

#%%
folder = ospath.dirname(__file__)
# dfs = load_data(
#     ospath.join(folder, 'data/initial_conditions.xlsx'), 
#     sheet=None,
#     )
# asinit = dfs['rBOD']
# fcinit = asinit.iloc[-1].to_dict()
# adinit = dfs['adm'].iloc[0].to_dict()
# Default initial conditions
dfs = load_data(ospath.join(folder, 'data/G1_init.xlsx'), sheet=None)
inf_concs = dfs['asm'].iloc[0].to_dict()
# c1init = dfs['asm'].iloc[1].to_dict()
asinit = dfs['asm'].iloc[1:]
# asinit = dfs['asm_ss']
adinit = dfs['adm'].iloc[0].to_dict()
c2init = dfs['settler'].to_dict('index')
c2init['s'] = {k:v for k,v in c2init['s'].items() if v>0}
c2init['x'] = {k:v for k,v in c2init['x'].items() if v>0}
c2init['tss'] = [v for k,v in c2init['tss'].items() if v>0]

MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

def create_f1_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('F1')
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
        sludge_flow_rate=0.074*MGD2cmd,
        solids_removal_efficiency=0.6
        )
    
    GT = su.IdealClarifier(
        'GT', PC-1, outs=['', 'thickened_PS'],
        sludge_flow_rate=0.026*MGD2cmd,
        solids_removal_efficiency=0.9
        )
    
    n_zones = 6
    V_tot = 2.4 * MGD2cmd
    
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
        rh=4e-4, rp=0.1, fns=0.01, 
        maximum_nonsettleable_solids=8.0
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
        sludge_flow_rate=0.00593*MGD2cmd,
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
        # FC.set_init_solubles(**fcinit)
        # FC.set_init_sludge_solids(**fcinit)
        # FC.set_init_TSS(default_fctss_init)
        # AD.set_init_conc(**default_ad_init)
        AD.set_init_conc(**adinit)
        FC.set_init_solubles(**c2init['s'])
        FC.set_init_sludge_solids(**c2init['x'])
        FC.set_init_TSS(c2init['tss'])
    
    sys = qs.System(
        'F1', 
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
    sys = create_f1_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 50
    # t = 1
    t_step = 1
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    
    run(sys, t, t_step, method=method)
    # biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[WAS],
    #               active_unit_IDs=('ASR',))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')