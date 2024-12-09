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
from exposan.werf import default_aed_init, default_as_init, default_fctss_init

__all__ = ('create_g2_system',)

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

def create_g2_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('G2')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
    carb = qs.WasteStream('carbon', T=Temp, units='kg/hr', S_A=60)
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
    V_tot = 4.7 * MGD2cmd
    fr_V = [0.014, 0.13, 0.148, 0.148, 0.28, 0.28]
    
    # gstrip = True
    # an_kwargs = dict(aeration=None, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=gstrip)
    # ae_kwargs = dict(aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=gstrip)
    
    # S1 = su.Splitter('S1', PC-0, split=0.8)
    
    # A1 = su.CSTR('A1', ins=[carb, 'RAS'], V_max=V_tot*fr_V[0], **an_kwargs)
    # A2 = su.CSTR('A2', [A1-0, S1-0], V_max=V_tot*fr_V[1], **an_kwargs)
    # A3 = su.CSTR('A3', [A2-0, 'intr', S1-1], V_max=V_tot*fr_V[2], **an_kwargs)
    # A4 = su.CSTR('A4', A3-0, V_max=V_tot*fr_V[3], **an_kwargs)
    # O5 = su.CSTR('O5', A4-0, V_max=V_tot*fr_V[4], **ae_kwargs)
    # O6 = su.CSTR('O6', O5-0, [1-A3, 'treated'], split=[40, 14],
    #               V_max=V_tot*fr_V[5], **ae_kwargs)
    
    ASR = su.PFR(
        'ASR', ins=[PC-0, 'RAS', carb], outs='treated',
        N_tanks_in_series=n_zones,
        V_tanks=[f*V_tot for f in fr_V],
        influent_fractions=[
            [0, 0.8, 0.2, 0,0,0],   # PE
            [1,0,0,0,0,0],          # RAS
            [1,0,0,0,0,0],          # carb
            ],
        internal_recycles=[(5,2,40*MGD2cmd)], DO_ID='S_O2',
        kLa=[0, 0, 0, 0, 180, 70], 
        DO_setpoints=[0, 0, 0, 0, 2.0, 2.0],
        suspended_growth_model=asm,
        gas_stripping=True
        )
    
    FC = su.FlatBottomCircularClarifier(
        'FC', ins=ASR-0, outs=['SE', 1-ASR, 'WAS'],
        # 'FC', ins=O6-1, outs=['SE', 1-A3, 'WAS'],
        underflow=0.4*10*MGD2cmd, wastage=0.1*MGD2cmd,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01, 
        maximum_nonsettleable_solids=8.0
        )
    
    # MT = su.Thickener(
    #     'MT', ins=FC-2, outs=['thickened_WAS', ''],
    #     thickener_perc=5, TSS_removal_perc=95,
    #     )
    MT = su.IdealClarifier(
        'MT', FC-2, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.019*MGD2cmd,
        solids_removal_efficiency=0.95
        )
    M1 = su.Mixer('M1', ins=[GT-1, MT-1])
    
    AED = su.AerobicDigester(
        'AED', ins=M1-0, outs='digestate',
        V_max=2.4*MGD2cmd, activated_sludge_model=asm,
        aeration=1.0, DO_ID='S_O2')
    
    # DW = su.Centrifuge(
    #     'DW', ins=AED-0, outs=('cake', ''),
    #     thickener_perc=18, TSS_removal_perc=90,
    #     )
    # M2 = su.Mixer('M2', ins=[GT-0, MT-1, DW-1])    
    DW = su.IdealClarifier(
        'DW', AED-0, outs=('', 'cake'),
        sludge_flow_rate=0.00593*MGD2cmd,
        solids_removal_efficiency=0.9
        )
    M2 = su.Mixer('M2', ins=[MT-0, DW-0])
    
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
        AED.set_init_conc(**default_aed_init)
        FC.set_init_solubles(**c2init['s'])
        FC.set_init_sludge_solids(**c2init['x'])
        FC.set_init_TSS(c2init['tss'])
    
    sys = qs.System(
        'G2', 
        path=(ASR, FC, MT, M1, AED, DW, M2, HD),
        # path=(A1, A2, A3, A4, O5, O6, FC, 
        #       MT, M1, AED, DW, M2, HD),
        recycle=(FC-1, HD-0)
        )

    sys.set_dynamic_tracker(FC-0, AED)

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
    sys = create_g2_system()
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