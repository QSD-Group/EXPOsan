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
    processes as pc,
    sanunits as su,
    )
from qsdsan.utils import ospath, time_printer, load_data, get_SRT
from exposan.werf import data_path

__all__ = ('create_i2_system',)

ID = 'I2'
#%%
dfs = load_data(
    ospath.join(data_path, 'initial_conditions.xlsx'), 
    sheet=None,
    )
asinit = dfs[ID]
fcinit = asinit.iloc[-1].to_dict()
default_fctss_init = [10, 12, 20, 40, 100, 500, 500, 500, 550, 1e4]
aedinit = dfs['AED'].loc[ID].to_dict()

MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]

def create_i2_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
    
    n_zones = 6
    Vs = [0.63, 1.5, 2.0, 2.0, 2.3, 0.21] # MG
    
    # ae_kwargs = dict(aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=True)
    # an_kwargs = dict(aeration=None, DO_ID='S_O2', suspended_growth_model=asm, gas_stripping=True)
        
    # A1 = su.CSTR('A1', [rww, 'RAS', 'reject'], V_max=Vs[0]*MGD2cmd, **an_kwargs)
    # A2 = su.CSTR('A2', [A1-0, 'intr'], V_max=Vs[1]*MGD2cmd, **an_kwargs)
    # O3 = su.CSTR('O3', A2-0, V_max=Vs[2]*MGD2cmd, **ae_kwargs)
    # O4 = su.CSTR('O4', O3-0, ('', 1-A2), split=[30, 17], V_max=Vs[3]*MGD2cmd, **ae_kwargs)
    # A5 = su.CSTR('A5', O4-0, V_max=Vs[4]*MGD2cmd, **an_kwargs)
    # O6 = su.CSTR('A6', A5-0, 'treated', V_max=Vs[5]*MGD2cmd, **ae_kwargs)
    
    ASR = su.PFR(
        'ASR', ins=[rww, 'RAS', 'reject'], outs='treated',
        N_tanks_in_series=n_zones,
        V_tanks=[v*MGD2cmd for v in Vs],
        influent_fractions=[[1]+[0]*5]*3,
        internal_recycles=[(3,1,30*MGD2cmd)], DO_ID='S_O2',
        kLa=[0, 0, 50, 50, 0, 50], 
        DO_setpoints=[0, 0, 2.0, 2.0, 0, 2.0],
        suspended_growth_model=asm,
        gas_stripping=True
        )
    
    FC = su.FlatBottomCircularClarifier(
        'FC', ins=ASR-0, outs=['SE', 1-ASR, 'WAS'],
        # 'FC', ins=O6-0, outs=['SE', 1-O1, 'WAS'],
        underflow=0.67*10*MGD2cmd, wastage=0.233*MGD2cmd,   # 12.3d SRT isn't sufficient for nitrification
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=6,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=3e-4, rp=5.2e-3, fns=0.001, 
        maximum_nonsettleable_solids=20.0
        )

    MT = su.IdealClarifier(
        'MT', FC-2, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.0313*MGD2cmd,
        solids_removal_efficiency=0.95
        )
    
    AED = su.AerobicDigester(
        'AED', ins=MT-1, outs='digestate',
        V_max=2.4*MGD2cmd, activated_sludge_model=asm,
        aeration=1.0, DO_ID='S_O2', gas_stripping=True)
       
    DW = su.IdealClarifier(
        'DW', AED-0, outs=('', 'cake'),
        sludge_flow_rate=3.35e-3*MGD2cmd,    # aim for 17% TS
        solids_removal_efficiency=0.9
        )
    MX = su.Mixer('MX', ins=[MT-0, DW-0])
    
    HD = su.HydraulicDelay('HD', ins=MX-0, outs=2-ASR)
    
    if default_init_conds:
        ASR.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**fcinit)
        FC.set_init_sludge_solids(**fcinit)
        FC.set_init_TSS(default_fctss_init)
        AED.set_init_conc(**aedinit)

    sys = qs.System(
        ID, 
        path=(ASR, FC, MT, AED, DW, MX, HD),
        # path=(A1, A2, O3, O4, A5, O6, FC, 
        #       MT, AED, DW, MX, HD),
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
    sys = create_i2_system()
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
    # biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[WAS],
    #               active_unit_IDs=('A1', 'A2', 'O3', 'O4', 'A5', 'O6'))
    #               # active_unit_IDs=('ASR'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))