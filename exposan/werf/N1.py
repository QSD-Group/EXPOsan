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

__all__ = ('create_n1_system',)

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
# c1init = dfs['asm'].iloc[1].to_dict()
asinit = dfs['asm'].iloc[1:]
# asinit = dfs['asm_ss']
adinit = dfs['adm'].iloc[0].to_dict()


MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

def create_n1_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('N1')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    pc.create_masm2d_cmps()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    rww = pc.create_masm2d_inf(
        'RWW', 10, 'MGD', T=Temp, 
        COD=358, NH4_N=25.91, PO4_P=5,
        fr_SI=0.05, fr_SF=0.16, fr_SA=0.024, fr_XI=0.2,
        )
    carb = qs.WasteStream('carbon', T=Temp, units='kg/hr', S_A=60)
    thermo_asm = qs.get_thermo()
    
    PC = su.PrimaryClarifier(
        'PC', ins=[rww, 'reject'], 
        outs=('PE', 'PS'),
        sludge_flow_rate=0.074*MGD2cmd,
        solids_removal_efficiency=0.6
        )
    
    n_zones = 5
    V_tot = 2.61 * MGD2cmd
    fr_V = [0.12, 0.18, 0.24, 0.24, 0.18, 0.04]
    Vs = [V_tot*f for f in fr_V]
    
    ASR = su.PFR(
        'ASR', ins=[PC-0, carb, 'intr'], 
        N_tanks_in_series=n_zones,
        V_tanks=Vs[:n_zones],
        influent_fractions=[
            [1,0,0,0,0,0],          # PE
            [1,0,0,0,0,0],          # carb
            [0,0,1,0,0,0],          # intr
            ],
        internal_recycles=[
            (1,0,10*MGD2cmd), 
            (3,1,40*MGD2cmd)], 
        DO_ID='S_O2',
        kLa=[0, 0, 50, 50, 0], 
        DO_setpoints=[0, 0, 3.0, 3.0, 0],
        suspended_growth_model=asm,
        gas_stripping=True
        )
    
    MBR = su.CompletelyMixedMBR(
        'MBR', ins=ASR-0, outs=('', 'WAS'),
        V_max=Vs[-1], pumped_flow=50, solids_capture_rate=0.9999,
        aeration=2.0, DO_ID='S_O2', suspended_growth_model=asm
        )
    
    S1 = su.Splitter('S1', MBR-0, (2-ASR, 'treated'), split=0.8)
    
    GT = su.IdealClarifier(
        'GT', ins=[PC-1, MBR-1], outs=['', 'thickened_sludge'],
        sludge_flow_rate=0.026*MGD2cmd,
        solids_removal_efficiency=0.85
        )

    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0)
    
    J1 = su.mASM2dtoADM1p('J1', upstream=GT-1, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR(
        'AD', ins=J1-0, outs=('biogas', 'digestate'), 
        V_liq=1.2*MGD2cmd, V_gas=0.12*MGD2cmd, 
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
    MX = su.Mixer('MX', ins=[GT-0, DW-0])
    
    HD = su.HydraulicDelay('HD', ins=MX-0, outs=1-PC)
    
    if default_init_conds:
        # ASR.set_init_conc(**default_as_init)
        ASR.set_init_conc(concentrations=asinit.iloc[:n_zones])
        MBR.set_init_conc(**asinit.iloc[-1].to_dict())
        # AD.set_init_conc(**default_ad_init)
        AD.set_init_conc(**adinit)
    
    sys = qs.System(
        'N1', 
        path=(PC, ASR, MBR, S1, GT, J1, AD, J2, DW, MX, HD),
        recycle=(S1-0, HD-0)
        )

    sys.set_dynamic_tracker(S1-1, AD)

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
    sys = create_n1_system()
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