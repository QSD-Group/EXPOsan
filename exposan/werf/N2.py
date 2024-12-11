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
    WasteStream,
    processes as pc,
    sanunits as su,
    )
from qsdsan.utils import ospath, time_printer, load_data, get_SRT
from exposan.werf import data_path

__all__ = ('create_n2_system',)

ID = 'N2'
#%%
dfs = load_data(
    ospath.join(data_path, 'initial_conditions.xlsx'), 
    sheet=None,
    )
asinit = dfs['N1']
fcinit = asinit.iloc[-1].to_dict()
aedinit = dfs['AED'].loc['I2'].to_dict()

MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

def create_n2_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
    carb = WasteStream('carbon', T=Temp, units='kg/hr', S_A=50)
    
    n_zones = 5
    V_tot = 2.61 * MGD2cmd
    fr_V = [0.12, 0.18, 0.24, 0.24, 0.18, 0.04]
    Vs = [V_tot*f for f in fr_V]
    
    ASR = su.PFR(
        'ASR', ins=[rww, carb, 'intr', 'reject'], 
        N_tanks_in_series=n_zones,
        V_tanks=Vs[:n_zones],
        influent_fractions=[
            [1,0,0,0,0],          # RWW
            [1,0,0,0,0],          # carb
            [0,0,1,0,0],          # intr from MBR
            [1,0,0,0,0],          # reject
            ],
        internal_recycles=[
            (1,0,10*MGD2cmd), 
            (3,1,20*MGD2cmd)], 
        DO_ID='S_O2',
        kLa=[0, 0, 50, 50, 0], 
        DO_setpoints=[0, 0, 3.0, 3.0, 0],
        suspended_growth_model=asm,
        gas_stripping=True
        )
    
    MBR = su.CompletelyMixedMBR(
        'MBR', ins=ASR-0, outs=('treated', 'WAS'),
        V_max=Vs[-1], pumped_flow=50, solids_capture_rate=0.9999,
        aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm
        )
    
    S1 = su.Splitter('S1', MBR-0, (2-ASR, 'SE'), split=0.8)
    
    MT = su.IdealClarifier(
        'MT', MBR-1, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.019*MGD2cmd,
        solids_removal_efficiency=0.95
        )
        
    AED = su.AerobicDigester(
        'AED', ins=MT-1, outs='digestate',
        V_max=2.4*MGD2cmd, activated_sludge_model=asm,
        aeration=1.0, DO_ID='S_O2', gas_stripping=True)
    # DW = su.Centrifuge(
    #     'DW', ins=J2-0, outs=('cake', ''),
    #     thickener_perc=18, TSS_removal_perc=90,
    #     )
    # M2 = su.Mixer('M2', ins=[GT-0, MT-1, DW-1])    
    DW = su.IdealClarifier(
        'DW', AED-0, outs=('', 'cake'),
        sludge_flow_rate=0.0053*MGD2cmd,    # aim for 17% TS
        solids_removal_efficiency=0.9
        )
    MX = su.Mixer('MX', ins=[MT-0, DW-0], outs=3-ASR)
    
    # HD = su.HydraulicDelay('HD', ins=MX-0)
    
    if default_init_conds:
        ASR.set_init_conc(concentrations=asinit.iloc[:n_zones])
        MBR.set_init_conc(**asinit.iloc[-1].to_dict())
        AED.set_init_conc(**aedinit)
    
    sys = qs.System(
        ID, 
        path=(ASR, MBR, S1, MT, AED, DW, MX),# HD),
        recycle=(S1-0, MX-0)
        )

    sys.set_dynamic_tracker(MBR, S1-1, AED)

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
    sys = create_n2_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    
    t = 10
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
    #               active_unit_IDs=('ASR', 'MBR))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))