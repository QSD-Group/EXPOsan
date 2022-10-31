# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 10:39:19 2022

@author: joy_c
"""
# from qsdsan import Components, set_thermo, Process, Processes, WasteStream, System
import qsdsan as qs, qsdsan.processes as pc, qsdsan.sanunits as su, numpy as np, \
    matplotlib.pyplot as plt, matplotlib as mpl#, matplotlib.ticker as tk
from qsdsan.utils import ospath, load_data
from exposan.metab_mock import DegassingMembrane as DM, results_path, figures_path
from chaospy import distributions as shape

def create_cmps():
    adm_cmps = pc.create_adm1_cmps(False)
    cmps = qs.Components([adm_cmps.S_h2, adm_cmps.S_ch4, adm_cmps.S_IC, adm_cmps.H2O])
    cmps.default_compile()
    qs.set_thermo(cmps)
    return cmps

R = 8.3145e-2 # Universal gas constant, [bar/M/K]
T_op = 298.15
rhos = np.zeros(6)
rhos[:3] = [1e-3,1e-2,2e-2]

def f_rhos(state_arr, params):
    # S_h2, S_ch4, S_IC, H2O, g_h2, g_ch4, g_co2, Q = state_arr
    KH = params['K_H']
    kLa = params['kLa']
    biogas_S = state_arr[:3].copy()
    biogas_p = R * T_op * state_arr[4:7]
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    rhos[0] = params['rh2']
    return rhos

def create_processes():
    p1 = qs.Process('H2_production',
                    reaction={'S_h2':1},
                    ref_component='S_h2',
                    conserved_for=())
    p2 = qs.Process('CH4_production', 
                    reaction={'S_ch4':1},
                    ref_component='S_ch4',
                    conserved_for=())
    p3 = qs.Process('CO2_production', 
                    reaction={'S_IC':1},
                    ref_component='S_IC',
                    conserved_for=())
    p4 = qs.Process('H2_transfer',
                    reaction={'S_h2':-1},
                    ref_component='S_h2',
                    conserved_for=())
    p5 = qs.Process('CH4_transfer',
                    reaction={'S_ch4':-1},
                    ref_component='S_ch4',
                    conserved_for=())
    p6 = qs.Process('CO2_transfer',
                    reaction={'S_IC':-1},
                    ref_component='S_IC',
                    conserved_for=())
    
    ps = qs.Processes([p1,p2,p3,p4,p5,p6])
    ps.compile()
    ps.set_rate_function(f_rhos)
    ps.rate_function._params = {
        'K_H': np.array([7.8e-4, 1.4e-3, 3.5e-2]),
        'kLa':200,
        'rh2':1e-3}
    
    ps.__dict__['_biogas_IDs'] = ('S_h2', 'S_ch4', 'S_IC')
    return ps

Q = 5           # influent flowrate [m3/d]
T1 = 273.15+35  # temperature [K]
Vl1 = 5         # liquid volume [m^3]
Vg1 = 0.556     # headspace volume [m^3]
tau_1 = 0.021   # degassing membrane retention time [d]
f_Q = 1    # recirculation rate

def create_system(fixed_hsp_P=False):
    cmps = create_cmps()
    ps = create_processes()
    inf = qs.WasteStream('inf')
    eff = qs.WasteStream('eff')
    bgh = qs.WasteStream('bgh', phase='g')
    bgm = qs.WasteStream('bgm', phase='g')
    inf.set_flow_by_concentration(Q, {}, units=('m3/d', 'kg/m3'))
    R1 = su.AnaerobicCSTR('R1', ins=[inf, 'return_1'], 
                          outs=(bgh, 'sidestream_1', eff), 
                          T=298.15, V_liq=Vl1, V_gas=Vg1, model=ps,
                          split=(f_Q, 1),
                          fixed_headspace_P=fixed_hsp_P)
    R1.set_init_conc(S_h2=8.503, S_ch4=0.0422, S_IC=1141.2)
    DM1 = DM('DM1', ins=R1-1, outs=(bgm, 1-R1), tau=tau_1)
    sys = qs.System('degas', path=(R1, DM1), recycle=(DM1-1,))
    return sys

#%%
def add_params(model, fix_H2_rprod=False):
    param = model.parameter
    R1 = model._system.flowsheet.unit.R1
    ps = R1.model
    
    b = 1
    D = shape.Uniform(1, 9)
    @param(name='side-to-main_flowrate_ratio', element=R1, kind='coupled',
           units='', baseline=b, distribution=D)
    def set_fQ(fq):
        R1.split = [fq, 1]
    
    if not fix_H2_rprod:
        b = 1e-3
        D = shape.LogUniform(-6*np.log(10), -0.1*np.log(10))
        @param(name='H2_production_rate', element=ps, kind='coupled', units='kgCOD/m3/d',
               baseline=b, distribution=D)
        def set_H2_rprod(r):
            ps.rate_function._params['rh2'] = r
    
    return model

KI_h2_fa=5e-6
KI_h2_c4=1e-5
KI_h2_pro=3.5e-6
def get_I(R1, Ki):
    Si = R1._state[0] 
    return pc.non_compet_inhibit(Si, Ki)
    
def add_metrics(model):
    metric = model.metric
    sreg = model._system.flowsheet.stream
    bgm, bgh, eff = sreg.bgm, sreg.bgh, sreg.eff
    R1 = model._system.flowsheet.unit.R1
    h2_i_mass = eff.components.S_h2.i_mass
    
    @metric(name='Sidestream_H2_extraction', units='kg/d', element='Biogas')
    def get_ss_H2():
        return bgm.imass['S_h2'] * h2_i_mass * 24
    
    @metric(name='Headspace_H2_extraction', units='kg/d', element='Biogas')
    def get_hs_H2():
        return bgh.imass['S_h2'] * h2_i_mass * 24
    
    @metric(name='Dissolved_H2', units='kgCOD/m3', element='Reactor')
    def get_S_h2():
        return eff.iconc['S_h2']/1e3
    
    @metric(name='H2_inhibition_fa', units='', element='Reactor')
    def get_Ifa():
        return get_I(R1, KI_h2_fa)
    
    @metric(name='H2_inhibition_c4', units='', element='Reactor')
    def get_Ic4():
        return get_I(R1, KI_h2_c4)
    
    @metric(name='H2_inhibition_pro', units='', element='Reactor')
    def get_Ipro():
        return get_I(R1, KI_h2_pro)
    
    return model



#%%
mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = False
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

def meshgrid_sample(p1, p2, n):
    xl, xu = p1.distribution.lower[0], p1.distribution.upper[0]
    yl, yu = np.log10(p2.distribution.lower[0]), np.log10(p2.distribution.upper[0])
    x = np.linspace(xl, xu, n)
    y = 10**np.linspace(yl, yu, n)
    xx, yy = np.meshgrid(x, y)
    samples = np.array([xx.reshape(n**2), yy.reshape(n**2)])
    return samples.T, xx, yy

def run_mapping(model, n, T, t_step, run=True, method='BDF', mpath=''):
    x = model.parameters[0]
    y = model.parameters[1]
    samples, xx, yy = meshgrid_sample(x, y, n)
    if run:
        model.load_samples(samples)
        t_span = (0, T)
        t_eval = np.arange(0, T+t_step, t_step)
        model.evaluate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method
            )
        if mpath: model.table.to_excel(mpath)
        
    return xx, yy

def fmt(x):
    s = f"{x:.1e}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"

irs = [0, 0, 1, 1, 2, 2]
ics = [0, 1, 0, 1, 0, 1]

def plot_heatmaps(xx, yy, model=None, path='', wide=False):
    if model: data = model.table
    else:
        path = path or ospath.join(results_path, 'table_2dv.xlsx')
        data = load_data(path, header=[0,1])
    zs = data.iloc[:,2:].to_numpy(copy=True)
    n = int(zs.shape[0] ** 0.5)
    zs = zs.T
    zz = zs.reshape((zs.shape[0], n, n))
    if wide:
        plts = zip(zz, ics, irs)
        fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(14,9))
    else:
        plts = zip(zz, irs, ics)
        fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(9,12))

    fig.set_tight_layout({'h_pad':4, 'w_pad': 0})
    for z, ir, ic in plts:
        ax = axes[ir, ic]
        ax.set_yscale('log')
        if ic > 0: 
            norm = 'linear'
            cfmt = lambda x: f"{x:.2g}"
        else: 
            norm = 'log'
            cfmt = fmt
        pos = ax.pcolormesh(xx, yy, z, norm=norm, shading='gouraud')
        cbar = fig.colorbar(pos, ax=ax)
        cbar.ax.tick_params(labelsize=14)
        ax.tick_params(axis='both', which='both', direction='inout', labelsize=14)
        ax2x = ax.secondary_xaxis('top')
        ax2x.tick_params(direction='in', which='both')
        ax2x.xaxis.set_ticklabels([])
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(direction='in', which='both')
        ax2y.yaxis.set_ticklabels([])
        cs = ax.contour(xx, yy, z, 
                        colors='white', norm=norm, origin='lower', 
                        linestyles='dashed', linewidths=1, 
                        extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]))
        ax.clabel(cs, cs.levels, inline=True, fmt=cfmt, fontsize=11)
    fig.savefig(ospath.join(figures_path, 'heatmaps.png'), dpi=300)


def dv_analysis(n=20, T=100, t_step=10, run=True, save_to='table_2dv.xlsx',
                plot=True, wide=True):
    path = ospath.join(results_path, save_to)
    sys = create_system()
    mdl = qs.Model(sys)
    mdl = add_params(mdl)
    mdl = add_metrics(mdl)
    xx, yy = run_mapping(mdl, n, T, t_step, run, mpath=path)
    if plot:
        if run:
            plot_heatmaps(xx, yy, mdl, wide=wide)
            for p in mdl.parameters:
                p.setter(p.baseline)
        else: 
            plot_heatmaps(xx, yy, path=path, wide=wide)
    return xx, yy

#%%
import pandas as pd
mpl.rcParams["figure.autolayout"] = True

def run_scenarios(model, N=50, 
                  T=100, t_step=10, method='BDF', mpath=''):
    X, = model.parameters
    xl, = X.distribution.lower
    xu, = X.distribution.upper
    samples = np.array([np.linspace(xl, xu, N),]).T
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    mpath = mpath or ospath.join(results_path, 'table_1dv.xlsx')
    ps = model._system.flowsheet.unit.R1.model
    dct = {}
    for rh2 in 10**np.linspace(-6,0,7):
        ps.rate_function._params['rh2'] = rh2
        model.evaluate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method
            )
        plot_area(model.table, rh2)
        dct[f'{rh2:.0e}'] = model.table.copy()
    writer=pd.ExcelWriter(mpath)
    for k, df in dct.items():
        df.to_excel(writer, sheet_name=k)
    writer.save()
        
def plot_area(data, rh2):
    name = f'{rh2:.0e}'
    x = data.iloc[:,0].to_numpy(copy=True)
    y_mem = data.iloc[:,1].to_numpy(copy=True)
    y_tot = y_mem + data.iloc[:,2].to_numpy()
    fig, ax = plt.subplots(figsize=(4,3))
    ax.plot(x, y_mem, '-', c='#D81B60', linewidth=0.4)
    area1 = ax.fill_between(x, 0., y_mem, color='#D81B60', alpha=0.7, label='Sidestream')
    ax.plot(x, y_tot, '-', c='#1E88E5', linewidth=0.4)
    area2 = ax.fill_between(x, y_mem, y_tot, color='#1E88E5', alpha=0.7, label='Headspace')
    ax.set_ylim(bottom=0., top=rh2*0.7)
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_xlabel('Recirculation rate [-]', fontsize=11)
    ax.set_ylabel('H2 extraction rate [kg/d]', fontsize=11)
    ax.legend(handles=[area1, area2], fontsize=10)
    ax.set_title(f'H2 production rate = {name} kgCOD/m3/d', fontsize=11)
    ax.tick_params(axis='both', which='both', direction='inout', labelsize=10)
    ax2x = ax.secondary_xaxis('top')
    ax2x.set_xlim(np.min(x), np.max(x))
    ax2x.tick_params(direction='in', which='both')
    ax2x.xaxis.set_ticklabels([])
    ax2y = ax.secondary_yaxis('right')
    ax2y.set_ylim(bottom=0., top=rh2*0.7)
    ax2y.tick_params(direction='in', which='both')
    ax2y.yaxis.set_ticklabels([])
    fig.savefig(ospath.join(figures_path, f'area_rH2_{name}.png'), dpi=300)
    del fig, ax

def dv1_analysis(N=50, T=100, t_step=10, run=True, fixed_hsp_P=False,
                 save_to='table_1dv_fixedP.xlsx', plot=True):
    path = ospath.join(results_path, save_to)
    if run:
        sys = create_system(fixed_hsp_P)
        mdl2 = qs.Model(sys)
        mdl2 = add_params(mdl2, True)
        mdl2 = add_metrics(mdl2)
        run_scenarios(mdl2, N=N, T=T, t_step=t_step, mpath=path)
    elif plot:
        dct = load_data(path, sheet=None, header=[0,1], 
                        skiprows=[2,], index_col=0)
        for k, v in dct.items():
            rh2 = float(k)
            plot_area(v, rh2)
    
    
#%%
if __name__ == '__main__':
    # xx, yy = dv_analysis(n=4, run=False)
    # xx, yy = dv_analysis(run=False)
    # xx, yy = dv_analysis()
    dv1_analysis(fixed_hsp_P=True)