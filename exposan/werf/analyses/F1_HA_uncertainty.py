# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import time as tm, pandas as pd, numpy as np, os
from exposan.werf import (
    create_system, 
    SelectiveRecovery,
    add_performance_metrics, 
    add_OPEX_metrics, 
    add_NH4_recovery_metric,
    add_downstream_uncertainty,
    opt_underflows,
    results_path
    )
from exposan.werf.utils import cache_state, load_state
from qsdsan import get_thermo, WasteStream, Model, System, processes as pc
from qsdsan.utils import get_SRT, ospath, load_data
from biosteam.evaluation._utils import var_columns

# import warnings
# warnings.filterwarnings("ignore")

MGD2cmd = 3785.412
f_rmv = 0.7 # (0.5, 0.8)

# influent strength: low, medium, high
# tech performance: % NH4-N recovery (50-80% reasonable range for ionic strength typical in domestic ww)
# maintain operation as adjusted unless violation?

# contextual uncertainties
# %%
def load_system_with_upstream_uncertainty(ID='F1'):
    sys = create_system(ID)
    cmps = get_thermo().chemicals
    for c in cmps:
        if c.organic:
            if c.i_N: c.i_N *= 0.7
            if c.i_P: c.i_P *= 1.1
    cmps.S_I.f_Vmass_Totmass = cmps.X_I.f_Vmass_Totmass = 0.6
    cmps.X_S.i_mass = cmps.X_I.i_mass = 0.63
    cmps.refresh_constants()
    u = sys.flowsheet.unit
    if 'AD' in u:
        u.J1._compile_reactions()
        u.J2._compile_reactions()
        
    fr_orgs = dict(fr_SI=0.075, fr_SF=0.24, fr_SA=0.05, fr_XI=0.3)
    wws = {}
    wws['low'] = pc.create_masm2d_inf(
        'low', 10, 'MGD', 
        COD=339, NH4_N=14, PO4_P=1.6, S_K=11, S_Cl=39, 
        **fr_orgs
        )
    wws['mid'] = pc.create_masm2d_inf(
        'mid', 10, 'MGD', 
        COD=508, NH4_N=20, PO4_P=2.4, S_K=16, S_Cl=59, 
        **fr_orgs
        )
    wws['high'] = pc.create_masm2d_inf(
        'high', 10, 'MGD', 
        COD=1016, NH4_N=41, PO4_P=4.7, S_K=32, S_Cl=118, 
        **fr_orgs
        )
    return sys, wws

# %%

def create_model(sys):
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    
    mdl = Model(sys)
    add_downstream_uncertainty(mdl)
    add_performance_metrics(mdl)
    add_OPEX_metrics(mdl)
    
    if "PS" in s: 
        if 'AD' in u: cake_tss = 18e4
        elif 'AED' in u: cake_tss = 17e4
        else: cake_tss = 20e4
    else: 
        cake_tss = 17e4
    
    if "thickened_WAS" in s: 
        thickened = s.thickened_WAS
        thickener = u.MT
    else: 
        thickened = s.thickened_sludge
        thickener = u.GT
    
    thickener.sludge_flow_rate, u.DW.sludge_flow_rate = opt_underflows[sys.ID]
    u.ASR.DO_setpoints *= 0
    u.ASR.DO_setpoints += 1
    
    load_state(sys, folder='steady_states/baseline_unopt')
    return mdl, cake_tss, thickened, thickener

def run_model(model=None, n_strength=3, interval=(0,1), samples=None, 
              N=200, rule='L', seed=None, **kwargs):
    if model is None:
        sys, wws = load_system_with_upstream_uncertainty()
        model, cake_tss, thickened, thickener = create_model(sys)
    else:
        wws = kwargs['wws']
        cake_tss = kwargs['cake_tss']
        thickened = kwargs['thickened']
        thickener = kwargs['thickener']
        sys = model.system

    ID = sys.ID
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    cmps = s.RWW.components
    _low = wws['low'].copy('low_strength')
    _high = wws['high'].copy('high_strength')
    
    if samples is None: samples = model.sample(N=N, rule=rule, seed=seed)
    # for strength, qwas in zip(('low', 'mid', 'high'), (0.15, 0.2, 0.3)):
    #     s.RWW.copy_flow(wws[strength])
    #     u.FC.wastage = qwas * MGD2cmd
    metrics = {}
    for high_frac in np.linspace(*interval, n_strength):
        _low.scale(1-high_frac)
        _high.scale(high_frac)
        s.RWW.mix_from([_low, _high])
        s.RWW._init_state()
        _low.copy_like(wws['low'])
        _high.copy_like(wws['high'])
        u.FC.wastage = 0.15 * (1+high_frac) * MGD2cmd
        u.FC._ODE = None
        sys._DAE = None
        # print(f"System {ID}, {strength}-strength ww")
        print(f"System {ID}, influent COD = {s.RWW.COD:.1f} mg/L")
        print("="*45)
        table = []
        start = tm.time()
        print("Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
        try: sys.simulate(t_span=(0,300), method='BDF')
        except: 
            try: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
        end = tm.time()
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
        print('Adjusting TS% ...')
        r_thick = thickened.get_TSS()/5e4
        r_cake = s.cake.get_TSS()/cake_tss
        while abs(r_thick - 1) > 0.01 or abs(r_cake - 1) > 0.01:
            print(f"{r_thick:.3f}  {r_cake:.3f}")
            thickener.sludge_flow_rate *= r_thick
            u.DW.sludge_flow_rate *= r_cake
            try: sys.simulate(t_span=(0,300), method='BDF')
            except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            r_thick = thickened.get_TSS()/5e4
            r_cake = s.cake.get_TSS()/cake_tss
        end2 = tm.time()
        print("Final underflows: ", f"{thickener.sludge_flow_rate:.2f}  {u.DW.sludge_flow_rate:.2f}")
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end2-end)))
        srt = get_SRT(sys, ('X_H', 'X_PAO', 'X_AUT'), wastage=[s.WAS], 
                      active_unit_IDs=('ASR', ))
        print(f'SRT = {srt:.2f} d')
        arr = u.ASR.state.iloc[:,:-1].to_numpy()
        mlss = np.sum(cmps.i_mass * cmps.x * arr, axis=1)
        print(f'MLSS = {np.mean(mlss):.0f} mg/L', '\n')
        # cache_state(sys, f'steady_states/HA_F1/{strength}', 'no_HA')
        cache_state(sys, 'steady_states/HA_F1/', f'{int(high_frac*100)}')

        for smp in samples:
            for p, v in zip(model.parameters, smp): p.setter(v)
            table.append([*smp] + [m() for m in model.metrics])            
        df = pd.DataFrame(table, 
            columns=var_columns([*model.parameters, *model.metrics])
            )
            
        # with pd.ExcelWriter(ospath.join(results_path, f'HA_{ID}_UA-{strength}.xlsx')) as writer:
        #     df.to_excel(writer, sheet_name='no_HA')
        metrics[f'{int(high_frac*100)}'] = df
    
    with pd.ExcelWriter(ospath.join(results_path, f'{ID}_UA.xlsx')) as writer:
        for k, df in metrics.items():
            df.to_excel(writer, sheet_name=k)
            
    return samples
            
#%%
def create_model_with_HA(sys):
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    location_insert = s.PE
    downstream_unit = location_insert.sink
    i_inlet = downstream_unit.ins.index(location_insert)
    
    i_path = sys.path.index(downstream_unit)
    location_insert.disconnect_sink()
    HA = SelectiveRecovery(
        'HydrogelAbsorbent', ins=location_insert, 
        outs=('Recovered_NH4', 'HA_eff'), 
        split={'S_NH4': f_rmv}, init_with='WasteStream'
        )
    HA-1-i_inlet-downstream_unit
    sys_ha = System(sys.ID+'ha', 
                    path=(*sys.path[:i_path], HA, *sys.path[i_path:]),
                    recycle=sys.recycle)
    sys_ha.set_dynamic_tracker(*sys.scope.subjects)
    mdl = Model(sys_ha)
    add_downstream_uncertainty(mdl)
    add_performance_metrics(mdl)
    add_OPEX_metrics(mdl)
    add_NH4_recovery_metric(mdl)
    
    if "PS" in s: 
        if 'AD' in u: cake_tss = 18e4
        elif 'AED' in u: cake_tss = 17e4
        else: cake_tss = 20e4
    else: 
        cake_tss = 17e4
    
    if "thickened_WAS" in s: 
        thickened = s.thickened_WAS
        thickener = u.MT
    else: 
        thickened = s.thickened_sludge
        thickener = u.GT
    
    thickener.sludge_flow_rate, u.DW.sludge_flow_rate = opt_underflows[sys.ID]
    u.ASR.DO_setpoints *= 0
    u.ASR.DO_setpoints += 1
    
    load_state(sys_ha, folder='steady_states/HA_opt')
    return mdl, HA, cake_tss, thickened, thickener

#%%

def run_model_with_HA(model=None, n_strength=3, interval=(0,1), 
                      samples=None, N=200, rule='L', seed=None,
                      removal_efficiencies=np.arange(0.5, 0.8, 0.05), 
                      **kwargs):
    if model is None:
        _sys, wws = load_system_with_upstream_uncertainty()
        model, HA, cake_tss, thickened, thickener = create_model_with_HA(_sys)
        ID = _sys.ID
    else:
        wws = kwargs['wws']
        HA = kwargs['HA']
        cake_tss = kwargs['cake_tss']
        thickened = kwargs['thickened']
        thickener = kwargs['thickener']
        ID = model.system.ID[:-2]
    sys_ha = model.system
    s = sys_ha.flowsheet.stream
    u = sys_ha.flowsheet.unit
    cmps = s.RWW.components
    i_nh4 = cmps.index('S_NH4')
    _low = wws['low'].copy('low_strength')
    _high = wws['high'].copy('high_strength')
    
    if samples is None: samples = model.sample(N=N, rule=rule, seed=seed)
    # for strength, qwas in zip(('low', 'mid', 'high'), (0.15, 0.2, 0.3)):
    #     s.RWW.copy_flow(wws[strength])
    #     u.FC.wastage = qwas * MGD2cmd
    for high_frac in np.linspace(*interval, n_strength):
        _low.scale(1-high_frac)
        _high.scale(high_frac)
        s.RWW.mix_from([_low, _high])
        s.RWW._init_state()
        _low.copy_like(wws['low'])
        _high.copy_like(wws['high'])
        u.FC.wastage = 0.15 * (1+high_frac) * MGD2cmd
        s.RWW._init_state()
        u.FC._ODE = None
        sys_ha._DAE = None
        out = {}
        # print(f"System {ID} + HA, {strength}-strength ww")
        print(f"System {ID} + HA, influent COD = {s.RWW.COD:.1f} mg/L")
        print("="*45)
        for f_rmv in removal_efficiencies:
            table = []
            HA.split[i_nh4] = f_rmv
            start = tm.time()
            print(f"{f_rmv:.0%} Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
            try: sys_ha.simulate(t_span=(0,300), method='BDF')
            except: sys_ha.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            end = tm.time()
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
            print('Adjusting TS% ...')
            r_thick = thickened.get_TSS()/5e4
            r_cake = s.cake.get_TSS()/cake_tss
            while abs(r_thick - 1) > 0.01 or abs(r_cake - 1) > 0.01:
                print(f"{r_thick:.3f}  {r_cake:.3f}")
                thickener.sludge_flow_rate *= r_thick
                u.DW.sludge_flow_rate *= r_cake
                try: sys_ha.simulate(t_span=(0,300), method='BDF')
                except: sys_ha.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
                r_thick = thickened.get_TSS()/5e4
                r_cake = s.cake.get_TSS()/cake_tss
            end2 = tm.time()
            print("Final underflows: ", f"{thickener.sludge_flow_rate:.2f}  {u.DW.sludge_flow_rate:.2f}")
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end2-end)))
            srt = get_SRT(sys_ha, ('X_H', 'X_PAO', 'X_AUT'), wastage=[s.WAS], 
                          active_unit_IDs=('ASR', ))
            print(f'SRT = {srt:.2f} d')
            arr = u.ASR.state.iloc[:,:-1].to_numpy()
            mlss = np.sum(cmps.i_mass * cmps.x * arr, axis=1)
            print(f'MLSS = {np.mean(mlss):.0f} mg/L', '\n')
            # cache_state(sys_ha, f'steady_states/HA_F1/{strength}', f'{f_rmv*100:.0f}')
            try: cache_state(sys_ha, f'steady_states/HA_F1/{int(high_frac*100)}', f'{f_rmv*100:.0f}')
            except:
                os.mkdir(ospath.join(results_path, f'steady_states/HA_F1/{int(high_frac*100)}'))
                cache_state(sys_ha, f'steady_states/HA_F1/{int(high_frac*100)}', f'{f_rmv*100:.0f}')
            
            for smp in samples:
                for p, v in zip(model.parameters, smp): p.setter(v)
                table.append([*smp] + [m() for m in model.metrics])            
            out[f"{f_rmv:.2f}"] = pd.DataFrame(table, 
                columns=var_columns([*model.parameters, *model.metrics])
                )
            
        # with pd.ExcelWriter(ospath.join(results_path, f'HA_{ID}_UA-{strength}.xlsx'), mode='a') as writer:
        with pd.ExcelWriter(ospath.join(results_path, f'HA_{ID}_UA-strength{int(high_frac*100)}.xlsx')) as writer:
            for k, df in out.items():
                df.to_excel(writer, sheet_name=k)
    
# %%
def compile_stats(ID='F1'):
    q = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]
    quantiles = {}
    for strength in ('low', 'mid', 'high'):
        dfs = load_data(ospath.join(results_path, f'HA_{ID}_UA-{strength}.xlsx'), 
                        header=[0,1], skiprows=[2,], sheet=None)
        opex = []
        for k, df in dfs.items():
            opex.append(df[('OPEX', 'Total OPEX [USD/d]')])
        opex = pd.DataFrame(opex, index=dfs.keys()).T
        quantiles[strength] = opex.quantile(q=q).T
        
    save_as = ospath.join(results_path, f'HA_{ID}_UA_opex_stats.xlsx')
    with pd.ExcelWriter(save_as) as writer:
        for k, v in quantiles.items():
            v.to_excel(writer, sheet_name=k)
    
# %%
def compile_stats_by_strength(ID='F1'):
    q = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]
    sys, wws = load_system_with_upstream_uncertainty()
    _low = wws['low'].copy('low_strength')
    _high = wws['high'].copy('high_strength')
    s = sys.flowsheet.stream
    infs = []
    for f in np.linspace(0,1,11):
        _low.scale(1-f)
        _high.scale(f)
        s.RWW.mix_from([_low, _high])
        infs.append([
            f*100, s.RWW.COD, s.RWW.BOD, s.RWW.get_TSS(), 
            s.RWW.TN, s.RWW.iconc['S_NH4'], 
            s.RWW.TP, s.RWW.iconc['S_PO4']
            ])
        _low.copy_like(wws['low'])
        _high.copy_like(wws['high'])
    infs = pd.DataFrame(infs, columns=[
        'Strength', 'COD', 'BOD', 'TSS', 'TN', 'NH4-N', 'TP', 'OP'
        ])
    infs.index = infs.Strength
    infs.drop(columns='Strength', inplace=True)
    
    bl_dfs = load_data(ospath.join(results_path, f'{ID}_UA.xlsx'), 
                       header=[0,1], skiprows=[2,], sheet=None)
    bl_opex = []
    for f in np.linspace(0,1,11):
        df = bl_dfs[f'{int(f*100)}']
        bl_opex.append(df[('OPEX', 'Total OPEX [USD/d]')])
    bl_opex = pd.DataFrame(bl_opex, index=infs.index).T
    bl_qs = bl_opex.quantile(q=q).T
    
    qs = []
    ms = []
    dqs = []
    dms = []
    for f in np.linspace(0,1,11):
        strength = f*100
        dfs = load_data(ospath.join(results_path, f'HA_{ID}_UA-strength{int(strength)}.xlsx'), 
                        header=[0,1], skiprows=[2,], sheet=None)
        noha = bl_opex[strength].to_numpy()
        means = []
        dmeans = []
        opex = []
        dopex = []
        for k, df in dfs.items():
            col = df[('OPEX', 'Total OPEX [USD/d]')]
            opex.append(col.to_numpy())
            dopex.append((noha-col).to_numpy())
            means.append(col.mean())
            dmeans.append((noha-col).mean())
        qs.append(np.quantile(opex, q))
        dqs.append(np.quantile(dopex, q))
        ms.append(means)
        dms.append(dmeans)
    
    qs = pd.DataFrame(qs, columns=q, index=infs.index)
    dqs = pd.DataFrame(dqs, columns=q, index=infs.index)
    ms = pd.DataFrame(ms, columns=[float(k) for k in dfs.keys()], index=infs.index)
    dms = pd.DataFrame(dms, columns=[float(k) for k in dfs.keys()], index=infs.index)

    out = pd.concat(
        [infs, bl_qs, qs, ms, dqs, dms], axis=1, keys=[
            'Influent', 'OPEX_quantiles_noHA', 
            'OPEX_quantiles', 'OPEX_mean_by_frmv', 
            'dOPEX_quantiles', 'dOPEX_mean_by_frmv']
        )
    out.to_excel(ospath.join(results_path, f'HA_{ID}_UA_opex_stats_by_strength.xlsx'))
    return out

# %%

if __name__ == '__main__':
    # smp = load_data(ospath.join(results_path, 'HA_F1_UA-low.xlsx'), header=[0,1], skiprows=[2,])
    # smp = smp.iloc[:,:5].to_numpy()
    # smp = run_model(seed=104)
    # smp = run_model(samples=smp, interval=(0.1, 0.9), n_strength=9)
    # run_model_with_HA(samples=smp, interval=(0.1, 0.9), n_strength=9)
    # compile_stats()
    out = compile_stats_by_strength()
