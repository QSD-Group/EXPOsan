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
    WasteStream,
    )
from qsdsan.utils import time_printer, ospath, load_data, get_SRT
from exposan.bsm2 import data_path

__all__ = ('create_system',)


# %%

Q = 20648.4 # influent flowrate [m3/d]
Q_intr = 3 * Q # activated sludge process internal recycle [m3/d]
Q_ras = Q # recycle sludge flowrate
Q_was = 600 # sludge wastage flowrate
Temp = 273.15+15 # temperature [K]
T_ad = 273.15+35
V_anae = 1000
V_anox = 1500 # anoxic zone tank volume
V_ae = 3000 # aerated zone tank volume

# O2 saturation concentration at 15 degC
SOSAT1 = 8

# Default initial conditions
dfs = load_data(ospath.join(data_path, 'bsm2p_init.xlsx'), sheet=None)
inf_concs = dfs['asm'].iloc[0].to_dict()
c1init = dfs['asm'].iloc[1].to_dict()
asinit = dfs['asm'].iloc[2:]
adinit = dfs['adm'].iloc[0].to_dict()
c2init = dfs['settler'].to_dict('index')
c2init['s'] = {k:v for k,v in c2init['s'].items() if v>0}
c2init['x'] = {k:v for k,v in c2init['x'].items() if v>0}
c2init['tss'] = [v for k,v in c2init['tss'].items() if v>0]
default_inf_kwargs = {
    'flow_tot': Q,
    'concentrations': inf_concs,
    'units': ('m3/d', 'mg/L'),
    }

def create_system(flowsheet=None, default_init_conds=True):
    flowsheet = flowsheet or qs.Flowsheet('bsm2p')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    unit = flowsheet.unit
    
    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm, 
                    electron_acceptor_dependent_decay=True,
                    k_h=2.46, mu_H=4.23, q_fe=2.11, b_H=0.28, mu_PAO=0.82, 
                    q_PP=1.23, q_PHA=2.46, b_PAO=0.14, b_PP=0.14, b_PHA=0.14, 
                    mu_AUT=0.61, b_AUT=0.09)
    thermo_asm = qs.get_thermo()
    
    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(**default_inf_kwargs)
    # carb = WasteStream('carbon', T=Temp)
    # carb.set_flow_by_concentration(2, {'S_A':400}, units=('m3/d', 'kg/m3'))
    
    # Primary clarifier using the Otterpohl-Freund model
    # Where are other parameters used?
    C1 = su.PrimaryClarifierBSM2(
        'C1', ins=(inf, 'reject'),
        outs=('C1_eff', 'C1_underflow'),
        isdynamic=True, 
        volume=900,
        f_corr=0.65,
        ratio_uf=0.007, # f_PS
        )

    DO_ID = 'S_O2'
    # aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=120, DOsat=SOSAT1, V=V_ae)
    # aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=60, DOsat=SOSAT1, V=V_ae)
    # A1 = su.CSTR('A1', ins=[C1-0, 'RWW', 'RAS', carb], V_max=V_an,
    #              aeration=None, suspended_growth_model=asm1)
    
    # A2 = su.CSTR('A2', A1-0, V_max=V_an,
    #              aeration=None, suspended_growth_model=asm1)

    # O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
    #              DO_ID=DO_ID, suspended_growth_model=asm1)
    
    # O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
    #              DO_ID=DO_ID, suspended_growth_model=asm1)
    
    # O3 = su.CSTR('O3', O2-0, [1-A1, 'treated'], 
    #              split=[0.6, 0.4],
    #              V_max=V_ae, aeration=aer3,
    #              DO_ID=DO_ID, suspended_growth_model=asm1)
    
    # AS = su.PFR('AS', ins=[C1-0, 'RAS', carb], outs='treated', 
    AS = su.PFR('AS', ins=[C1-0, 'RAS'], outs='treated', 
                N_tanks_in_series=7,
                V_tanks=[V_anae]*2+[V_anox]*2+[V_ae]*3,
                # influent_fractions=[[1]+[0]*6]*3,
                influent_fractions=[[1]+[0]*6]*2,
                internal_recycles=[(6,2,Q_intr)],
                kLa=[0]*4+[120,120,60], DO_ID=DO_ID, DO_sat=SOSAT1,
                suspended_growth_model=asm,
                gas_stripping=True)
    
    C2 = su.FlatBottomCircularClarifier(
        'C2', AS-0, ['effluent', 1-AS, 'WAS'],
        # 'C2', O3-1, ['effluent', 2-A1, 'WAS'],
        underflow=Q_ras, wastage=Q_was,
        surface_area=1500, height=4, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=474, v_max_practical=250,
        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3
        )
    
    TC1 = su.Thickener('TC1', C2-2, outs=['thickened_sludge', ''],
                       thickening_perc=7, TSS_removal_perc=98)
    M1 = su.Mixer('M1', ins=(C1-1, TC1-0))
        
    # Switch to ADM1 components for the anaerobic digester
    cmps_adm = pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(
        f_bu_su=0.1328, f_pro_su=0.2691, f_ac_su= 0.4076,
        q_ch_hyd=0.3, q_pr_hyd=0.3, q_li_hyd=0.3, 
        )
    
    # breakpoint()
    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR('AD', ins=J1.outs[0], outs=('biogas', 'AD_eff'), isdynamic=True,
                           V_liq=3400, V_gas=300, T=T_ad, model=adm,)
    AD.algebraic_h2 = False
    # Switch back to ASM1 components
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True, 
                          adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)
    
    # Dewatering
    C3 = su.Centrifuge(ID='C3', ins=J2-0, outs=['digested_sludge', ''],
                       thickening_perc=28, TSS_removal_perc=98)

    M2 = su.Mixer('M2', ins=(TC1-1, C3-1), outs=1-C1)
    
    if default_init_conds:
        # for i in ('A1', 'A2', 'O1', 'O2', 'O3'):
        #     getattr(unit, i).set_init_conc(**asinit[i])
        C1.set_init_conc(**c1init)
        AS.set_init_conc(concentrations=asinit)
        C2.set_init_solubles(**c2init['s'])
        C2.set_init_sludge_solids(**c2init['x'])
        C2.set_init_TSS(c2init['tss'])
        AD.set_init_conc(**adinit)
    
    sys = qs.System('bsm2_sys', 
                    # path=(C1, A1, A2, O1, O2, O3, C2, 
                    #       TC1, M1, J1, AD, J2, C3, M2),
                    # recycle=(O3-0, C2-1, M2-0)
                    path=(C1, AS, C2, TC1, M1, J1, AD, J2, C3, M2),
                    recycle=(C2-1, M2-0)
                    )
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    sys.maxiter = 5000
    # sys.set_dynamic_tracker(C1, A1, O3, C2, J1, AD, J2, C3)
    # sys.set_dynamic_tracker(C1, AS, C2, J1, AD, J2, C3)
    sys.set_dynamic_tracker(AS, AD, C2-0)    
    
    return sys

def create_subsys():
    flowsheet = qs.Flowsheet('bsm1p')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm, 
                    electron_acceptor_dependent_decay=True,
                    k_h=2.46, mu_H=4.23, q_fe=2.11, b_H=0.28, mu_PAO=0.82, 
                    q_PP=1.23, q_PHA=2.46, b_PAO=0.14, b_PP=0.14, b_PHA=0.14, 
                    mu_AUT=0.61, b_AUT=0.09)
    
    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(flow_tot=20446, concentrations=c1init, units=('m3/d', 'mg/L'))
    # inf.set_flow_by_concentration(**default_inf_kwargs)

    DO_ID = 'S_O2'

    AS = su.PFR('AS', ins=[inf, 'RAS'], outs='treated', 
                N_tanks_in_series=7,
                V_tanks=[V_anae]*2+[V_anox]*2+[V_ae]*3,
                # influent_fractions=[[1]+[0]*6]*3,
                influent_fractions=[[1]+[0]*6]*2,
                internal_recycles=[(6,2,Q_intr)],
                kLa=[0]*4+[120,120,60], DO_ID=DO_ID, DO_sat=SOSAT1,
                suspended_growth_model=asm,
                gas_stripping=True)
    
    C2 = su.FlatBottomCircularClarifier(
        'C2', AS-0, ['effluent', 1-AS, 'WAS'],
        # 'C2', O3-1, ['effluent', 2-A1, 'WAS'],
        underflow=Q_ras, wastage=Q_was,
        surface_area=1500, height=4, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=474, v_max_practical=250,
        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3
        )

    AS.set_init_conc(concentrations=asinit)
    C2.set_init_solubles(**c2init['s'])
    C2.set_init_sludge_solids(**c2init['x'])
    C2.set_init_TSS(c2init['tss'])
    
    sub = qs.System('bsm1p', path=(AS, C2), recycle=(C2-1, ))
    sub.set_dynamic_tracker(AS, C2-0)
    
    return sub

# %%

@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        # t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    
    # biomass_IDs = ('X_BH', 'X_BA')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[sys.flowsheet.stream.digested_sludge],
    #               active_unit_IDs=('C3'))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

#%%
if __name__ == '__main__':
    # sys = create_system()
    sys = create_subsys()
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
    # sys._setup()
    # sys.converge()
    # sys.diagram()
    # sys.diagram(file=os.path.join(figures_path, 'bsm2_sys'), format='png')