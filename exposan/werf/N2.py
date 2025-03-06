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
    carb = WasteStream('carbon', T=Temp, units='kg/hr', S_A=100)
    
    n_zones = 5
    V_tot = 2.61 * MGD2cmd
    fr_V = [0.12, 0.18, 0.24, 0.24, 0.18, 0.04]
    Vs = [V_tot*f for f in fr_V]
    Q_was = 0.17 * MGD2cmd      # SRT = 10.71 d
    Q_intr = 40 * MGD2cmd
    
    ASR = su.PFR(
        'ASR', ins=[rww, carb, 'RAS', 'reject'], 
        N_tanks_in_series=n_zones,
        V_tanks=Vs[:n_zones],
        influent_fractions=[
            [1,0,0,0,0],          # RWW
            # [1,0,0,0,0],          # carb
            [0.5,0,0,0,0.5],          # carb
            [0,0,1,0,0],          # RAS from MBR
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
        'MBR', ins=ASR-0, outs=('SE', ''),
        V_max=Vs[-1], solids_capture_rate=0.9999, pumped_flow=Q_was+Q_intr,
        aeration=2.0, DO_ID='S_O2', gas_stripping=True,
        suspended_growth_model=asm, 
        )
    S1 = su.Splitter('S1', MBR-1, ('', 'WAS'), split=Q_intr/(Q_intr+Q_was))
    HD1 = su.HydraulicDelay('HD1', ins=S1-0, outs=2-ASR)
        
    MT = su.IdealClarifier(
        'MT', S1-1, outs=['', 'thickened_WAS'],
        sludge_flow_rate=0.0359*MGD2cmd,
        solids_removal_efficiency=0.95
        )
        
    AED = su.AerobicDigester(
        'AED', ins=MT-1, outs='digestate',
        V_max=2.4*MGD2cmd, activated_sludge_model=asm,
        aeration=1.0, DO_ID='S_O2', gas_stripping=True)

    DW = su.IdealClarifier(
        'DW', AED-0, outs=('', 'cake'),
        sludge_flow_rate=0.00383*MGD2cmd,    # aim for 17% TS
        solids_removal_efficiency=0.9
        )
    MX = su.Mixer('MX', ins=[MT-0, DW-0])
    HD2 = su.HydraulicDelay('HD2', ins=MX-0, outs=3-ASR)
    
    if default_init_conds:
        ASR.set_init_conc(concentrations=asinit.iloc[:n_zones])
        MBR.set_init_conc(**asinit.iloc[-1].to_dict())
        AED.set_init_conc(**aedinit)
    
    sys = qs.System(
        ID, 
        path=(ASR, MBR, S1, HD1, MT, AED, DW, MX, HD2),
        recycle=(HD1-0, HD2-0)
        )

    sys.set_dynamic_tracker(MBR, MBR-0, AED)

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
    #               active_unit_IDs=('ASR', 'MBR'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')
    
    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))