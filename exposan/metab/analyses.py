# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan.utils import ospath, load_data, SanUnitScope
from exposan.metab import create_system, create_model, run_model, meshgrid_sample,\
    data_path, results_path, figures_path
from exposan.metab.utils import categorize_cashflow, categorize_all_impacts
from time import time
import qsdsan as qs, pandas as pd, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns
from scipy.stats import kstest

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = False
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

#%% run simulations
def run_discrete_DVs(samples_path):
    qs.PowerUtility.price = 0.0913
    dct = load_data(samples_path, sheet=None)
    for i in (
            'UASB',
            'FB',
            'PB',
            ):
        for n, j in (
                (1,'P'), 
                (2,'P'), 
                (2,'M'), (2,'H'),
                ):
            sys = create_system(n_stages=n, reactor_type=i, gas_extraction=j)
            print(sys.ID)
            print('='*len(sys.ID))
            mdl = create_model(sys, kind='DV', exception_hook='warn')
            sample = dct[i+'_'+j].to_numpy()
            run_model(mdl, sample)

def run_UA_SA(seed=None, N=1000, rule='L'):
    seed = seed or int(str(time())[-3:])
    qs.PowerUtility.price = 0.0913
    sample = None
    for rt in (
            'UASB', 
            'FB', 
            'PB'
            ):
        sys = create_system(n_stages=1, reactor_type=rt, gas_extraction='P', T=22)
        print(sys.ID)
        print('='*len(sys.ID))
        mdl = create_model(sys, kind='uasa', exception_hook='warn')
        if sample is None: sample = mdl.sample(N=N, rule=rule, seed=seed)
        run_model(mdl, sample, seed=seed)
    return sample

def _rerun_failed_samples(seed, rt='PB'):
    qs.PowerUtility.price = 0.0913  
    sample = load_data(ospath.join(results_path, f'{rt}1P_{seed}.xlsx'),
                       header=[0,1], skiprows=[2,])
    sample = sample[sample.isna().any(axis=1)]
    sample = sample.iloc[:, :18].to_numpy()

    sys = create_system(n_stages=1, reactor_type=rt, gas_extraction='P', T=22)
    cmps = sys.feeds[0].components
    C_bulk = np.array([
        1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
        # 0,0,0,0,0,0,
        4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 3.212e-02,
        # 0,0,0,0,0,0,
        # 2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 
        0,0,0,0,
        9.416, 2.478, 0.968, 2.846, 1.796, 
        1.48 , 0.734, 
        # 1, 0.6,
        # 4.708e+00, 1.239e+00, 4.838e-01, 1.423e+00, 8.978e-01, 
        # 2.959e+00, 1.467e+00,
        # 4.924e-02, 4.000e-02, 2.000e-02, 9.900e+02
        0, 4.000e-02, 2.000e-02, 9.900e+02
        ])
    C = dict(zip(cmps.IDs, C_bulk))
    sys.units[0].set_init_conc(**C)
    mdl = create_model(sys, kind='uasa', exception_hook='raise')
    run_model(mdl, sample, seed='temp')

from exposan.metab import ks_22
def run_PB_over_diffusivity_HRT(n=10, run=True):
    sys = create_system(reactor_type='PB')
    mdl = qs.Model(sys, exception_hook='warn')
    param = mdl.parameter
    metric = mdl.metric
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    u.R1.bead_diameter = 1.0
    u.R1.voidage = 0.45
    u.R1.model.rate_function.params['rate_constants'] = ks_22

    @param(name='Bead-to-water diffusivity fraction', units='', kind='coupled',
           element='Encapsulation', baseline=0.55, bounds=(0.2, 1.1))
    def set_f_diff(f):
        u.R1.f_diff = f
    
    @param(name='HRT', units='d', kind='coupled', element='System',
           baseline=1/3, bounds=(0.2, 1))
    def set_tau(tau):
        V = s.inf.F_vol * 24 * tau
        u.R1.V_liq = V
        u.R1.V_gas = V*0.1
        u.R1._prep_model()
        u.R1._compile_ODE()
        u.R1.scope = SanUnitScope(u.R1)   
    
    @metric(name='COD removal', units='%', element='Process')
    def get_rcod():
        rcod = 1 - s.eff_dg.COD/s.inf.COD
        return rcod*100
    
    @metric(name='CH4 collection', units='kg/d', element='Biogas')
    def get_QCH4():
        return (s.bg.imass['S_ch4'] + s.bge.imass['S_ch4'])*0.2506728377314149*24

    @metric(name='fugitive CH4 emission', units='kg/d', element='Biogas')
    def get_fug_CH4():
        return s.eff_dg.imass['S_ch4']*0.2506728377314149*24
    
    C0_bulk = np.array([
        1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
        4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 3.212e-02,
        2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 
        4.708e+00, 1.239e+00, 4.838e-01, 1.423e+00, 8.978e-01, 
        2.959e+00, 1.467e+00, 
        4.924e-02, 4.000e-02, 2.000e-02, 9.900e+02
        ])
    C = dict(zip(s.inf.components.IDs, C0_bulk))
    samples, gridx, gridy = meshgrid_sample(*mdl.parameters, n)
    mdl.load_samples(samples)
    if run:
        ys = []
        i = 0
        for smp in samples:
            i += 1
            print(f'\n{i}  {"="*20}')
            for p, v in zip(mdl.parameters, smp): p.setter(v)
            u.R1.set_init_conc(**C)
            sys.simulate(state_reset_hook='reset_cache', method='BDF', t_span=(0, 400))
            ys.append([m() for m in mdl.metrics])
        
        mdl.table.iloc[:,-3:] = ys
        mpath = ospath.join(results_path, 'PB_diffusivity_HRT.xlsx')
        mdl.table.to_excel(mpath)
    return mdl, samples


#%% discrete-DV scatter plot
def data_compile(save=True):
    dfs = []
    for i in ('UASB','FB','PB',):
        for n, j in ((1,'P'), (2,'P'), (2,'M'), (2,'H'),):
            df = load_data(ospath.join(results_path, f'{i}{n}{j}.xlsx'), 
                           header=[0,1], skiprows=[2,])
            # df.drop(columns=['Biomass','Process'], level='Element', inplace=True)
            df.drop(columns=['Biomass'], level='Element', inplace=True)            
            df = df.droplevel('Element', axis=1)
            df['Reactor type'] = i
            df['Number of stages'] = n
            df['Gas extraction'] = j
            dfs.append(df)
    out = pd.concat(dfs, ignore_index=True)
    out['id'] = out.index
    out.columns = [i.split(' [')[0] for i in out.columns]
    non_tea = [i for i in out.columns if not i.startswith('Levelized cost')]
    out = pd.melt(out, non_tea, var_name='Effluent degassing', value_name='Levelized cost')
    out.dropna(axis=0, subset=['Levelized cost', ], inplace=True)
    out['Effluent degassing'] = out.apply(lambda row: row['Effluent degassing'].split(' (')[1].rstrip(')'), axis=1)
    non_lca = [i for i in out.columns if not i.startswith('GWP100')]
    out = pd.melt(out, non_lca, var_name='Bead lifetime', value_name='GWP100')
    out.dropna(axis=0, subset=['GWP100', ], inplace=True)
    out['Bead lifetime'] = out.apply(lambda row: row['Bead lifetime'].split(' (')[1].rstrip(')'), axis=1)
    out.drop(out[out['Effluent degassing'] != out['Bead lifetime']].index, inplace=True)
    out['Effluent degassing'] = out.apply(lambda row: 0 if '/o' in row['Effluent degassing'] else 1, axis=1)
    out['Bead lifetime'] = out.apply(lambda row: 1 if '1yr' in row['Bead lifetime'] else \
                                     10 if '10yr' in row['Bead lifetime'] else \
                                      30  if '30yr' in row['Bead lifetime'] else None, axis=1)
    if save:
        out.to_excel(ospath.join(results_path, 'table_compiled.xlsx'))
    return out

def plot_clusters(data=None, save_as='', partial=True):
    if data is None:
        data = load_data(ospath.join(results_path, 'table_compiled.xlsx'))
    if partial: fig, ax = plt.subplots(figsize=(8,8))
    else: fig, ax = plt.subplots(figsize=(18,12))
    # pal = {'UASB':'#60c1cf', 'FB':'#F98F60', 'PB':'#a280b9'}
    # pal = {'UASB': (96,193,207), 'FB': (249,143,96), 'PB':(162,128,185)}
    # edge = ['#000000' if dg else pal[rct] for dg, rct in 
    # edge = [(0,0,0) if dg else pal[rct] for dg, rct in 
            # data.loc[:,['Effluent degassing','Reactor type']].to_numpy()]
    leg = False if partial else 'auto'
    ax = sns.scatterplot(
        data=data,
        x='Levelized cost',
        y='GWP100',
        hue='Reactor type',
        # palette=pal,
        size='Number of stages',
        sizes=(30, 55),
        style='Gas extraction',
        markers={'H':'s', 'M':'^', 'P':'o'},
        # edgecolors=edge,
        legend=leg,
        alpha=0.7,
        linewidths=2,
        ax=ax,
        )
    ax.tick_params(axis='both', which='major', direction='inout', length=10, labelsize=18)
    ax.tick_params(axis='both', which='minor', direction='inout', length=6)
    if partial:
        ax.set_xlim(0, 5000)
        ax.set_ylim(-100, 1500)
    else:
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=-100)
        ax.legend(fontsize=16)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax2x = ax.secondary_xaxis('top')
    ax2x.tick_params(axis='x', which='major', direction='in', length=5)
    ax2x.tick_params(axis='x', which='minor', direction='in', length=3)
    ax2x.xaxis.set_major_formatter(plt.NullFormatter())
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=5)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=3)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    suffix = 'partial' if partial else 'full'
    path = save_as or ospath.join(figures_path, f'clusters_{suffix}.png')
    fig.savefig(path, dpi=300, 
                facecolor='white'
                )
    
    return fig, ax
    
#%% compute diff by discrete-DV
def calc_diff(df, factor, baseline_value, norms):    
    bl = df.xs(key=baseline_value, axis=1, level=factor)
    df.drop(labels=baseline_value, axis=1, level=factor, inplace=True)
    df.index = range(df.shape[0])
    out = []
    # for i in ('Levelized cost', 'GWP100'):
    for i in ('COD removal', 'Levelized cost', 'GWP100'):        
        diff = df.xs(i, axis=1)
        if i == 'COD removal': diff = diff.sub(bl[i].values, axis=0)
        else: diff = diff.sub(bl[i].values, axis=0)/norms[i]
        diff = pd.melt(diff, value_name=i, ignore_index=False)
        diff['id'] = diff.index
        diff = diff.set_index([c for c in diff.columns if c != i])
        out.append(diff)
    # out = out[0].join(out[1])
    out = out[0].join(out[1]).join(out[2])    
    out['pair'] = out.index.get_level_values(factor)   
    out.index = range(out.shape[0])
    return out

groups = (
           # ('Reactor type', ['Voidage'], 'PB'),
           ('Reactor type', ['Bead lifetime','Bead diameter', 'Voidage'], 'UASB'),
           ('Gas extraction', ['Recirculation ratio', 'Headspace pressure'], 'P'),
           ('Headspace pressure', [], 0.4),
           ('Recirculation ratio', [], 1),
           ('Number of stages', [], 1),
           ('Temperature', [], 22),
           ('Effluent degassing', [], 0),
           ('Total HRT', [], 1),
           ('Bead diameter', [], 2),
           ('Voidage', [], 0.9),
           ('Bead lifetime', [], 30),
          )

def compare_DVs(data=None, save_as=''):
    if data is None:
        data = load_data(ospath.join(results_path, 'table_compiled.xlsx'))
    data = data.drop(['id', 'H2 yield', 'CH4 yield'], axis=1)
    # vals = ['Levelized cost', 'GWP100']
    vals = ['COD removal', 'Levelized cost', 'GWP100']
    norms = data.loc[:,vals].max(axis=0) - data.loc[:,vals].min(axis=0)
    dvs = set(data.columns) - set(vals)
    out = {}
    for factor, covar, bl in groups:
        idx = list(dvs - {factor, *covar})
        cols = [factor, *covar]
        df = data.dropna(axis=0, subset=factor)
        df = df.pivot(index=idx, columns=cols, values=vals)
        df.dropna(axis=0, inplace=True)
        out[factor] = calc_diff(df, factor, bl, norms)
    path = save_as or ospath.join(results_path, 'diff_rcod.xlsx')
    with pd.ExcelWriter(path) as writer:
        for k,v in out.items():
            v.to_excel(writer, sheet_name=k)
    return out  
    

#%% contour of diff
boxprops = dict(alpha=1, edgecolor='black')
# flierprops = dict(marker='.', markersize=1, markerfacecolor='#90918e', markeredgecolor='#90918e')
meanprops = dict(marker='x', markersize=6, markerfacecolor='black', markeredgecolor='black', lw=0.6)
medianprops = dict(color='black', lw=0.6)

# palette = ['#60C1CF', '#F98F60', '#79BF82', '#F3C354', '#A280B9', '#ED586F', 
#            '#35767F', '#733763', '#4D7E53', '#AB8937', '#5184EF']

# def plot_joint(df, save_as='', kde=True):
#     g = sns.JointGrid(height=8, ratio=5, space=0, marginal_ticks=True)
#     x, y = df['Levelized cost'], df['GWP100']
#     group = df['group']
#     pal = palette[:group.nunique()]
#     g.refline(x=0, y=0, color='#90918e', lw=1)
#     if kde:
#         sns.kdeplot(x=x, y=y, hue=group, ax=g.ax_joint,
#                     palette=pal,
#                     common_norm=False,
#                     fill=True,
#                     legend=False,
#                     alpha=0.6,
#                     )
#     else:
#         sns.scatterplot(x=x, y=y, hue=group, ax=g.ax_joint,
#                         palette=pal,
#                         legend=False,
#                         size=30,
#                         alpha=0.7,
#                         )
#     # g.ax_joint.set_ylim(top=3000)
#     g.ax_joint.tick_params(axis='both', which='major', direction='inout', length=10, labelsize=18)
#     g.ax_joint.tick_params(axis='both', which='minor', direction='inout', length=6)
#     g.ax_joint.set_xlabel('')
#     g.ax_joint.set_ylabel('')
#     bxp_kwargs = dict(
#         palette=pal,
#         linewidth=0.5,
#         showcaps=True,
#         showmeans=True,
#         showfliers=False,
#         dodge=False,
#         saturation=1,
#         width=0.3,
#         boxprops=boxprops,
#         # flierprops=flierprops,
#         meanprops=meanprops,
#         medianprops=medianprops,
#         )
#     if group.nunique() > 3: bxp_kwargs.pop('width')
#     sns.boxplot(x=x, y=group, hue=group, ax=g.ax_marg_x, **bxp_kwargs)
#     sns.boxplot(y=y, x=group, hue=group, ax=g.ax_marg_y, **bxp_kwargs)
#     for ax in (g.ax_marg_x, g.ax_marg_y): ax.legend_.remove()
#     g.ax_marg_x.tick_params(axis='x', which='major', direction='out', length=5)
#     g.ax_marg_x.tick_params(axis='x', which='minor', direction='out', length=3)    
#     g.ax_marg_x.tick_params(left=False, which='both', labelleft=False)    
#     g.ax_marg_x.spines['left'].set_color('white')
#     # g.ax_marg_y.set_ylim(top=3000)
#     g.ax_marg_y.tick_params(axis='y', which='major', direction='out', length=5)
#     g.ax_marg_y.tick_params(axis='y', which='minor', direction='out', length=3) 
#     g.ax_marg_y.tick_params(bottom=False, which='both', labelbottom=False)    
#     g.ax_marg_y.spines['bottom'].set_color('white')
    
#     path = ospath.join(figures_path, save_as)
#     g.savefig(path, dpi=300, facecolor='white')



# def plot_diff(data=None):
#     if data is None:
#         data = load_data(ospath.join(results_path, 'diff.xlsx'), sheet=None)
#     singles = ('Reactor type', 'Total HRT', 'Bead lifetime')
#     for i in singles:
#         df = data[i]
#         df['group'] = df['pair'].astype(str)
#         plot_joint(df, f'{i}.png', kde=(i != 'Bead lifetime'))
#     others = []
#     for k, v in data.items():
#         if k not in singles:
#             v['group'] = v.apply(lambda row: f"{k}_{row['pair']}", axis=1)
#             others.append(v)
#     df = pd.concat(others)
#     plot_joint(df, 'other_DVs.png')

palette = ['#a280b9', '#f3c354', '#90918E']

def plot_joint(df, save_as='', kde=True):
    g = sns.JointGrid(height=8, ratio=5, space=0, marginal_ticks=True)
    x, y = df['Levelized cost'], df['GWP100']
    group = df['group']
    pal = palette[:group.nunique()]
    g.refline(x=0, y=0, color='#4e4e4e', lw=1)
    if kde:
        sns.kdeplot(x=x, y=y, hue=group, ax=g.ax_joint,
                    palette=pal,
                    common_norm=False,
                    fill=True,
                    # linewidths=1.5,
                    legend=False,
                    # alpha=0.5,
                    thresh=0.05,
                    )
    else:
        sns.scatterplot(x=x, y=y, hue=group, style=group, size=group, ax=g.ax_joint,
                        palette=pal,
                        legend=False,
                        sizes=[100, 200],
                        markers=['P', 'v'],
                        alpha=0.5,
                        )
    g.ax_joint.tick_params(axis='both', which='major', direction='inout', length=14, labelsize=24)
    g.ax_joint.tick_params(axis='both', which='minor', direction='inout', length=8)
    g.ax_joint.set_xlim(-0.2, 1.0)
    g.ax_joint.set_ylim(-0.3, 1.0)
    g.ax_joint.set_xlabel('')
    g.ax_joint.set_ylabel('')
    bxp_kwargs = dict(
        palette=pal,
        width=0.5,
        showcaps=True,
        showmeans=True,
        showfliers=False,
        dodge=False,
        whis=(5,95),
        saturation=1,
        boxprops=boxprops,
        meanprops=meanprops,
        medianprops=medianprops,
        )
    sns.boxplot(x=x, y=group, hue=group, ax=g.ax_marg_x, **bxp_kwargs)
    sns.boxplot(y=y, x=group, hue=group, ax=g.ax_marg_y, **bxp_kwargs)
    # for ax in (g.ax_marg_x, g.ax_marg_y): ax.legend_.remove()
    g.ax_marg_x.tick_params(axis='x', which='major', direction='out', length=7)
    g.ax_marg_x.tick_params(axis='x', which='minor', direction='out', length=4)    
    g.ax_marg_x.tick_params(left=False, which='both', labelleft=False)    
    g.ax_marg_x.spines['left'].set_color('white')
    g.ax_marg_y.tick_params(axis='y', which='major', direction='out', length=7)
    g.ax_marg_y.tick_params(axis='y', which='minor', direction='out', length=4) 
    g.ax_marg_y.tick_params(bottom=False, which='both', labelbottom=False)    
    g.ax_marg_y.spines['bottom'].set_color('white')
    
    path = ospath.join(figures_path, save_as)
    g.savefig(path, dpi=300, facecolor='white')

def plot_diff(data=None):
    if data is None:
        data = load_data(ospath.join(results_path, 'diff_rcod.xlsx'), sheet=None)
    for i, df in data.items():
        df['group'] = df['pair'].astype(str)
        if i in ('Total HRT', 'Bead diameter', 'Reactor type'):
            df = df.sort_values('pair', ascending=False)
        elif i == 'Gas extraction':
            df = df.sort_values('pair', ascending=True)            
        plot_joint(df, f'{i}.png', kde=(i != 'Bead lifetime'))
    
    
#%% cost & impact breakdowns
def best_breakdown():
    kwargs = dict(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
    llc_bd = {}
    imp_bd = {}
    
    # UASB
    mdl1 = create_model(n_stages=2, reactor_type='UASB')
    sys1 = mdl1.system
    sub1, = sys1.subsystems
    
    # best at TEA
    smp = np.array([1, 22])
    for p, v in zip(mdl1.parameters, smp): p.setter(v)
    sys1.simulate(**kwargs)
    llc_bd['uasb_t1'] = categorize_cashflow(sub1.TEA)
    imp_bd['uasb_t1'] = categorize_all_impacts(sub1.LCA)
    
    # best at LCA  
    mdl1 = create_model(n_stages=1, reactor_type='UASB')
    sys1 = mdl1.system
    sub1, = sys1.subsystems
    smp = np.array([4, 22])
    for p, v in zip(mdl1.parameters, smp): p.setter(v)
    sys1.simulate(**kwargs)
    llc_bd['uasb_l1'] = categorize_cashflow(sys1.TEA)
    imp_bd['uasb_l1'] = categorize_all_impacts(sys1.LCA)
    
    # FB
    # best at TEA
    mdl2 = create_model(n_stages=2, reactor_type='FB')
    sys2 = mdl2.system
    u2 = sys2.flowsheet.unit
    u2.R1.bead_lifetime = u2.R2.bead_lifetime = 30
    sub2, = sys2.subsystems
    
    smp = np.array([2, 0.9, 2, 22])
    for p, v in zip(mdl2.parameters, smp): p.setter(v)
    sys2.simulate(**kwargs)
    llc_bd['fb_t1'] = categorize_cashflow(sub2.TEA)
    imp_bd['fb_t1'] = categorize_all_impacts(sub2.LCA)
    
    # best at LCA
    mdl2 = create_model(n_stages=1, reactor_type='FB')
    sys2 = mdl2.system
    u2 = sys2.flowsheet.unit
    u2.R1.bead_lifetime = 30
    sub2, = sys2.subsystems
    smp = np.array([2, 0.9, 4, 22])
    for p, v in zip(mdl2.parameters, smp): p.setter(v)
    sys2.simulate(**kwargs)
    llc_bd['fb_l1'] = categorize_cashflow(sys2.TEA)
    imp_bd['fb_l1'] = categorize_all_impacts(sys2.LCA)
    
    # PB
    # best at TEA & LCA
    mdl3 = create_model(n_stages=2, reactor_type='PB')
    sys3 = mdl3.system
    u3 = sys3.flowsheet.unit
    u3.R1.bead_lifetime = u3.R2.bead_lifetime = 30
    sub3, = sys3.subsystems
    
    smp = np.array([2, 1, 22])
    for p, v in zip(mdl3.parameters, smp): p.setter(v)
    sys3.simulate(**kwargs)
    llc_bd['pb_t1'] = categorize_cashflow(sub3.TEA)
    imp_bd['pb_t1'] = categorize_all_impacts(sub3.LCA)
    
    # best at LCA
    mdl3 = create_model(n_stages=1, reactor_type='PB')
    sys3 = mdl3.system
    u3 = sys3.flowsheet.unit
    u3.R1.bead_lifetime = 30
    sub3, = sys3.subsystems
    smp = np.array([2, 1, 22])
    for p, v in zip(mdl3.parameters, smp): p.setter(v)
    sys3.simulate(**kwargs)
    llc_bd['pb_l1'] = categorize_cashflow(sys3.TEA)
    imp_bd['pb_l1'] = categorize_all_impacts(sys3.LCA)
    
    llc_bd = pd.DataFrame.from_dict(llc_bd).transpose()
    llc_bd['fug_ch4'] = 0
    imp_bd = pd.concat(imp_bd).swaplevel(0, 1)
    with pd.ExcelWriter(ospath.join(results_path, 'breakdown.xlsx')) as writer:
        llc_bd.to_excel(writer, sheet_name='cost')
        for i in sys3.LCA.indicators:
            imp_bd.loc[i.ID].to_excel(writer, sheet_name=i.ID)
    
    return llc_bd, imp_bd

#%% stacked bar of breakdowns
patch_dct = {
    'vessel': ('#60C1CF', r'\\\\\\'),
    'beads': ('#79BF82', ''),
    'dm': ('#35767F', ''),
    'others': ('#ED586F', ''),
    'electricity': ('#A280B9', '////'),
    'heat_onsite': ('#90918e', ''),
    'chemicals': ('#733763', ''),
    'biogas_offset': ('#F98F60', ''),
    'fug_ch4': ('#F3C354', '')
    }

def stacked_bar(data, save_as=''):
    fig, ax = plt.subplots(figsize=(4, 4))
    fig.subplots_adjust(left=0.2, right=0.95)
    x = range(data.shape[0])
    ax.axhline(y=0, color='black', linewidth=1)
    yp = np.zeros(data.shape[0])
    yn = yp.copy()
    for k, v in patch_dct.items():
        c, hat = v
        y = data.loc[:,k]
        y_offset = (y>=0)*yp + (y<0)*yn
        ax.bar(x, y, width=0.55, bottom=y_offset, color=c, hatch=hat)
        yp += (y>=0) * y
        yn += (y<0) * y
    ax.scatter(x, data.loc[:,'total'], marker='x', c='white', zorder=10)

    ax.tick_params(axis='y', which='major', direction='inout', length=10, labelsize=12)
    ax.tick_params(axis='y', which='minor', direction='inout', length=6)
    # ax.ticklabel_format(axis='y', scilimits=[-2,3], useMathText=True)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=5)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=3)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    ax.tick_params(bottom=False, which='both', labelbottom=False)  

    if save_as:
        folder = ospath.join(figures_path, 'breakdown')
        fig.savefig(ospath.join(folder, save_as), dpi=300, transparent=True)
    else: return fig, ax
    
def plot_breakdown(data=None):
    if data is None:
        data = load_data(ospath.join(results_path, 'breakdown.xlsx'), sheet=None)
    for k, df in data.items():
        stacked_bar(df, save_as=f'{k}.png')
    
#%% Monte-Carlo Filtering

def MCF_encap_to_susp(seed, save=True):
    data = {}
    for i in ('UASB', 'FB', 'PB'):
        data[i] = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                            header=[0,1], skiprows=[2,])
    uasb = data['UASB']
    mdl = create_model(kind='uasa')
    nx = len(mdl.parameters)
    x = uasb.iloc[:, :nx]
    pair_cols = [
            ('Process', 'COD removal [%]'),
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    outs = {}
    for rt in ('FB', 'PB'):
        ys = data[rt]
        dys = ys.loc[:,pair_cols] - uasb.loc[:,pair_cols]
        dct = {'D': pd.DataFrame(), 'p': pd.DataFrame(), f'{rt}-UASB': dys}
        for y, col in dys.items():
            lower_xs = x.loc[col <= 0]
            higher_xs = x.loc[col > 0]
            Ds = []
            ps = []
            for i in range(nx):
                D, p = kstest(lower_xs.iloc[:,i], higher_xs.iloc[:,i])
                Ds.append(D)
                ps.append(p)
            y_name = y[1].split(' [')[0]
            dct['D'][y_name] = Ds
            dct['p'][y_name] = ps
        dct['D'].index = dct['p'].index = x.columns
        if rt == 'FB': dct['D'].loc['PB'] = dct['p'].loc['PB'] = None
        else: dct['D'].loc['FB'] = dct['p'].loc['FB'] = None
        dct['D'].loc['Membrane', ['Levelized cost (w/o degas)', 'GWP100 (w/o degas)']] = None
        dct['p'].loc['Membrane', ['Levelized cost (w/o degas)', 'GWP100 (w/o degas)']] = None
        if save:
            with pd.ExcelWriter(ospath.join(results_path, f'KStest_{rt}-v-UASB.xlsx')) as writer:
                for name, df in dct.items(): df.to_excel(writer, sheet_name=name)
        
        out = []
        for name, df in dct.items(): 
            if name in 'Dp': 
                df['Parameter'] = df.index.get_level_values('Feature')
                out.append(df.melt(id_vars='Parameter', var_name='Metric', value_name=name))
        outs[rt] = out[0].merge(out[1])
    return outs

def MCF_pb_to_fb(seed, save=True):
    dys = load_data(ospath.join(results_path, f'diff_3ways_{seed}.xlsx'),
                   sheet='pf', header=[0,1], skiprows=[2,], nrows=1000)
    mdl = create_model(kind='uasa')
    nx = len(mdl.parameters)
    pb = load_data(ospath.join(results_path, f'PB1P_{seed}.xlsx'),
                   header=[0,1], skiprows=[2,], nrows=1000)
    x = pb.iloc[:, :nx]
    pair_cols = [
            ('Process', 'COD removal [%]'),
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    dys = dys.loc[:,pair_cols]
    dct = {'D': pd.DataFrame(), 'p': pd.DataFrame()}
    for y, col in dys.items():
        lower_xs = x.loc[col <= 0]
        higher_xs = x.loc[col > 0]
        if 0 < len(lower_xs) < 1000:
            Ds = []
            ps = []
            for i in range(nx):
                D, p = kstest(lower_xs.iloc[:,i], higher_xs.iloc[:,i])
                Ds.append(D)
                ps.append(p)
            y_name = y[1].split(' [')[0]
            dct['D'][y_name] = Ds
            dct['p'][y_name] = ps
    dct['D'].index = dct['p'].index = x.columns
    dct['D'].loc['UASB'] = dct['p'].loc['UASB'] = None
    dct['D'].loc['Membrane', ['Levelized cost (w/o degas)', 'GWP100 (w/o degas)']] = None
    dct['p'].loc['Membrane', ['Levelized cost (w/o degas)', 'GWP100 (w/o degas)']] = None
    if save:
        with pd.ExcelWriter(ospath.join(results_path, 'KStest_pb_v_fb.xlsx')) as writer:
            for name, df in dct.items(): df.to_excel(writer, sheet_name=name)
    
    out = []
    for name, df in dct.items(): 
        if name in 'Dp': 
            df['Parameter'] = df.index.get_level_values('Feature')
            out.append(df.melt(id_vars='Parameter', var_name='Metric', value_name=name))
    outs = out[0].merge(out[1])
    return outs

def _plot_bubble(df, ax, pal):
    # pal = ('#f3c354', '#a280b9')
    pal = ('#79bf82', '#60c1cf')
    dv = np.array(([0]*7+[1]+[0]*5+[1]*3+[0,1])*df.Metric.nunique())
    sns.scatterplot(df, x='Metric', y='Parameter', 
                    # hue=df.p<0.05, palette=['black', pal], 
                    hue=(df.p<0.05)*(1+dv), palette=['black', *pal],                     
                    size='D', sizes=(0, 350), size_norm=(0,1),
                    legend=False,
                    ax=ax,
                    )
    ax.set_xlim(-0.5, len(set(df.Metric))-0.5)
    ax.set_axisbelow(True)
    ax.grid(True, which='major', color='k', linestyle='--', 
            linewidth=0.7, alpha=0.3)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='both', which='major', direction='inout', length=8)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

def MCF_bubble_plot(data=None, seed=None):
    if data is None: 
        data = MCF_encap_to_susp(seed, False)
    mpl.rcParams['xtick.minor.visible'] = False
    mpl.rcParams['ytick.minor.visible'] = False
    if isinstance(data, dict):
        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(6,9))
        pal = {'FB':'#f3c354', 'PB':'#a280b9'}
        for kv, ax in zip(data.items(), axes):
            k, df = kv
            clr = pal[k]
            _plot_bubble(df, ax, clr)
        save_as = 'MCF_2575.png'
    else:
        fig, ax = plt.subplots(figsize=(2.4,9))
        clr = '#60c1cf'
        _plot_bubble(data, ax, clr)
        save_as = 'MCF_pb_v_fb.png'

    fig.subplots_adjust(wspace=0)
    fig.savefig(ospath.join(figures_path, save_as), dpi=300, transparent=True)

#%% MCF intra-system
def MCF_25_vs_75(seed, save=True):
    data = {}
    for i in ('FB', 'PB'):
        data[i] = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                            header=[0,1], skiprows=[2,], nrows=1000)
    pb = data['PB']
    mdl = create_model(kind='uasa')
    nx = len(mdl.parameters)
    x = pb.iloc[:, :nx]
    pair_cols = [
            ('Process', 'COD removal [%]'),
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    outs = {}
    for rt in ('FB', 'PB'):
        ys = data[rt].loc[:,pair_cols]
        dct = {'D': pd.DataFrame(), 'p': pd.DataFrame(), f'{rt}': ys}
        for y, col in ys.items():
            y_name = y[1].split(' [')[0]
            if y_name == 'COD removal': thres = np.percentile(col, 75)
            else: thres = np.percentile(col, 25)
            lower_xs = x.loc[col <= thres]
            higher_xs = x.loc[col > thres]
            Ds = []
            ps = []
            for i in range(nx):
                D, p = kstest(lower_xs.iloc[:,i], higher_xs.iloc[:,i])
                Ds.append(D)
                ps.append(p)
            dct['D'][y_name] = Ds
            dct['p'][y_name] = ps
        dct['D'].index = dct['p'].index = x.columns
        if rt == 'FB': dct['D'].loc['PB'] = dct['p'].loc['PB'] = None
        else: dct['D'].loc['FB'] = dct['p'].loc['FB'] = None
        dct['D'].loc['UASB'] = dct['p'].loc['UASB'] = None
        dct['D'].loc['Membrane', ['Levelized cost (w/o degas)', 'GWP100 (w/o degas)']] = None
        dct['p'].loc['Membrane', ['Levelized cost (w/o degas)', 'GWP100 (w/o degas)']] = None
        if save:
            with pd.ExcelWriter(ospath.join(results_path, f'KStest_{rt}.xlsx')) as writer:
                for name, df in dct.items(): df.to_excel(writer, sheet_name=name)
        
        out = []
        for name, df in dct.items(): 
            if name in 'Dp': 
                df['Parameter'] = df.index.get_level_values('Feature')
                out.append(df.melt(id_vars='Parameter', var_name='Metric', value_name=name))
        outs[rt] = out[0].merge(out[1])
    return outs

def plot_1dkde(x, groups, x_bounds, y_bounds, prefix=''):
    n = groups.shape[1]
    i = 0
    for m, group in groups.items():
        i += 1
        fig, ax = plt.subplots(figsize=(5, 2))
        sns.kdeplot(x=x, hue=group, ax=ax,
                    palette=['#1A85FF','#D41159'],
                    common_norm=False,
                    fill=True,
                    legend=False,
                    alpha=0.3,
                    cut=0,
                    )
        ax.set_xlim(*x_bounds)
        ax.set_ylim(*y_bounds)
        ax.tick_params(axis='both', which='major', direction='inout', length=8, labelsize=11)
        ax.tick_params(axis='both', which='minor', direction='inout', length=5)
        ax.ticklabel_format(axis='y', scilimits=[-2,3], useMathText=True)
        if i < n: ax.xaxis.set_ticklabels([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax2x = ax.secondary_xaxis('top')
        ax2x.tick_params(direction='in', which='major', length=4)
        ax2x.tick_params(direction='in', which='minor', length=2.5)
        ax2x.xaxis.set_ticklabels([])
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(direction='in', which='major', length=4, labelsize=11)
        ax2y.tick_params(direction='in', which='minor', length=2.5)
        ax2y.yaxis.set_ticklabels([])
        m_name = m[0].replace('/', '')
        fig.savefig(ospath.join(figures_path, f'kde_{prefix}_{m_name}.png'),
                    dpi=300, transparent=True)
        del fig, ax

def plot_univariate_kdes(seed):
    pair_cols = [
            ('Process', 'COD removal [%]'),
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    for i in (
            'FB', 
            # 'PB',
              ):
        df = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                       header=[0,1], skiprows=[2,], nrows=1000)
        ys = df.loc[:,pair_cols]
        # x = df.loc[:,('Encapsulation', 'Bead diameter [mm]')]
        # x = df.loc[:,('System', 'Total HRT [d]')]
        x = 1-df.loc[:,('FB', 'FB voidage')]
        # x = df.loc[:,('PB', 'PB voidage')]
        # x = df.loc[:,('ADM1', 'Uptake k ac [COD/COD/d]')]
        thres = ys.quantile(0.25)
        thres[0] = np.percentile(ys.iloc[:,0], 75)
        groups = ys > thres
        groups.iloc[:,0] = ys.iloc[:,0] < thres[0]
        plot_1dkde(x, groups, [0.03,0.25], [0, 8.5], i+'-void')
        # plot_1dkde(x, groups, [0.35,0.45], [0, 22], i+'-void')
        # plot_1dkde(x, groups, [1/6, 3], [0, 1], i+'-hrt')
        # plot_1dkde(x, groups, [1,5], [0, 0.55], i+'-beaddia')
        # plot_1dkde(x, groups, [3.9,16], [0, 0.16], i+'-kac')
        print(thres)
    
#%% pair-wise comparisons, three-way
def calc_3way_diff(seed, save=True):
    data = {}
    for i in ('UASB', 'FB', 'PB'):
        data[i] = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                            header=[0,1], skiprows=[2,], nrows=1000)
    pair_cols = [
            ('Process', 'COD removal [%]'),
            ('Biogas', 'H2 production [kg/d]'),
            ('Biogas', 'CH4 production [kg/d]'),
            # ('Biomass', 'R1 Overall TSS [g/L]'),
            # ('Biomass', 'R1 Encapsulated TSS [g/L]'),
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    # _col = ('Biomass', 'R1 Encapsulated TSS [g/L]')
    out = {}
    out['pu'] = (data['PB'].loc[:, pair_cols] - data['UASB'].loc[:, pair_cols])
    out['fu'] = (data['FB'].loc[:, pair_cols] - data['UASB'].loc[:, pair_cols])
    # pair_cols.insert(4, _col)
    # out['pu'][_col] = out['fu'][_col] = None
    out['pu'] = out['pu'].loc[:, pair_cols]
    out['fu'] = out['fu'].loc[:, pair_cols]
    out['pf'] = (data['PB'].loc[:, pair_cols] - data['FB'].loc[:, pair_cols])
    
    if save:
        with pd.ExcelWriter(ospath.join(results_path, f'diff_3ways_{seed}.xlsx')) as writer:
            for k, df in out.items(): df.to_excel(writer, sheet_name=k)
    return out

#%%

# pb = load_data(ospath.join(results_path, 'PB1P_364.xlsx'),
#                header=[0,1], skiprows=[2,])

# x = pb.loc[:, ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]')]
# y = pb.loc[:, ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]')]
# # zs = pb.iloc[:, :18]
# sens = [
#         ('ADM1', 'Uptake k pro [COD/COD/d]'),
#         ('ADM1', 'Uptake k ac [COD/COD/d]'),
#         ('System',  'Total HRT [d]'),
#         ('Encapsulation', 'Max encapsulation density [gTSS/L]'),
#         ('Encapsulation', 'Bead lifetime [yr]'),
#         ('Encapsulation', 'Bead diameter [mm]')
#         ]

# for i, col in enumerate(sens):
#     fig, ax = plt.subplots(figsize=(4,3))
#     pos = ax.scatter(x=x, y=y, c=pb.loc[:,col], cmap='viridis', s=1, alpha=0.7)
#     # ax.set_xlim(100)
#     # ax.set_ylim(10)
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.legend('right')
#     fig.colorbar(pos, ax=ax)
#     fig.savefig(ospath.join(figures_path, f'scatter_{i}.png'), dpi=300, facecolor='white')
    
#%% breakdown all cost & gwp

def breakdown_and_sort(data):
    tea = data.loc[:,'TEA (w/ degas)'].copy()
    tea.sort_values(by='Levelized cost (w/ degas) [$/ton rCOD]', inplace=True)
    llc = tea.pop('Levelized cost (w/ degas) [$/ton rCOD]')
    tea = (tea / 100).mul(llc.abs(), axis='rows')
    tea.columns = [col[4:-15].replace(' ', '_') for col in tea.columns]
    tea['fug_ch4'] = 0
    tea['total'] = llc
    
    lca = data.loc[:,'LCA (w/ degas)'].copy()
    lca.sort_values(by='GWP100 (w/ degas) [kg CO2eq/ton rCOD]', inplace=True)
    gwp = lca.pop('GWP100 (w/ degas) [kg CO2eq/ton rCOD]')
    # absolute = (gwp<=0).any()
    # if absolute: 
    lca = (lca / 100).mul(gwp.abs(), axis='rows')
    lca.columns = [col[4:-15].replace(' ', '_') for col in lca.columns]
    lca['total'] = gwp
    return tea, lca#, absolute

def plot_area(df, absolute=False, ylims=None):
    fig, ax = plt.subplots(figsize=(4,4))
    fig.subplots_adjust(left=0.15, right=0.95)
    x = range(df.shape[0])
    ax.plot(x, df.total, color='white', linewidth=1.5, linestyle='dashed')
    ax.axhline(y=0, color='black', linewidth=0.5)
    yp = np.zeros(df.shape[0])
    yn = yp.copy()
    for k, v in patch_dct.items():
        c, hat = v
        y = df.loc[:,k]
        y_offset = (y>=0)*yp + (y<0)*yn
        ax.fill_between(x, y+y_offset, y_offset, facecolor=c, 
                        hatch=hat, linewidth=0.5, zorder=0)
        yp += (y>=0) * y
        yn += (y<0) * y
    ax.set_xlim(0, df.shape[0])
    if ylims: ax.set_ylim(*ylims)
    ax.tick_params(labelsize=12)
    ax.tick_params(axis='y', which='major', direction='inout', length=8)
    ax.tick_params(axis='y', which='minor', direction='inout', length=4)
    ax.set_xlabel('')
    ax.set_ylabel('')
    if absolute: 
        # ax.ticklabel_format(axis='y', scilimits=[-2,3], useMathText=True)
        # ax.yaxis.get_offset_text().set_fontsize(12)
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=4)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=2)
        ax2y.yaxis.set_ticklabels([])
    else:
        ax2y = ax.twinx()
        ax2y.plot(x, df['total'], color='black', linewidth=0.5)
        ax2y.tick_params(axis='y', which='major', direction='inout', length=8)
        ax2y.tick_params(axis='y', which='minor', direction='inout', length=4)
    return fig, ax

def breakdown_uasa(seed):
    for i in (
            # 'UASB', 'FB', 
            'PB',):
        data = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                         header=[0,1], skiprows=[2,], nrows=1000)
        # tea, lca, absolute = breakdown_and_sort(data)
        tea, lca = breakdown_and_sort(data)
        for suffix, df in (('TEA', tea), ('LCA', lca)):
            # fig, ax = plot_area(df, absolute and suffix=='LCA')
            if i in ('FB', 'PB'):
                if suffix == 'TEA': ylims = (-500, 20000)
                else: ylims = (-300, 3500)
            else:
                ylims = None
            fig, ax = plot_area(df, True, ylims)
            # fig, ax = plot_area(df, True)
            fig.savefig(ospath.join(figures_path, f'breakdown/{i}_{suffix}_abs.png'),
                        dpi=300, 
                        # facecolor='white',
                        transparent=True,
                        )

#%% Spearman's rho between cost & GWP

from scipy.stats import spearmanr
def Spearman_corr(seed, save=True):
    cols = [
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    stats = {}
    for i in ('FB', 'PB'):
        data = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                         header=[0,1], skiprows=[2,], nrows=1000)
        m = data.loc[:, cols]
        m.droplevel(1, axis=1)
        stats[i] = spearmanr(m)
    if save:
        with pd.ExcelWriter(ospath.join(results_path, 'spearman.xlsx')) as writer:
            for k,v in stats.items():
                r, p = v
                r = pd.DataFrame(r, index=m.columns, columns=m.columns)
                r.to_excel(writer, sheet_name=f'{k}_rho')
                p = pd.DataFrame(p, index=m.columns, columns=m.columns)
                p.to_excel(writer, sheet_name=f'{k}_p')
    return stats

#%% Spearman's rho between K_TSS and steady-state TSS
# fb = load_data(ospath.join(results_path, 'FB1P_187.xlsx'),
#                header=[0,1], skiprows=[2,], nrows=1000)
# pb = load_data(ospath.join(results_path, 'PB1P_187.xlsx'),
#                header=[0,1], skiprows=[2,], nrows=1000)
# mtss = fb.loc[:,('Encapsulation', 'Max encapsulation density [gTSS/L]')]
# ftss = fb.loc[:,('Biomass', 'R1 encapsulated TSS [g/L]')]
# ptss = pb.loc[:,('Biomass', 'R1 encapsulated TSS [g/L]')]

# sum(ftss / mtss)/1000
# sum(ptss / mtss)/1000
# f_sprho = spearmanr(ftss, mtss)
# p_sprho = spearmanr(ptss, mtss)

#%% heatmaps
# from scipy.interpolate import griddata

def togrid(df, mdl, n):
    zs = df.iloc[:,-len(mdl.metrics):].to_numpy().T
    us = df.iloc[:, :-(len(mdl.metrics)+len(mdl.parameters))].to_numpy().T
    samples, xx, yy = meshgrid_sample(*mdl.parameters, n)
    zzs = np.vstack((us, zs))
    zzs = zzs.reshape((zzs.shape[0], *xx.shape))
    return xx, yy, zzs

def plot_heatmap(xx, yy, z, z2=None, baseline=[], save_as='', specific=False, hrt=True):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    if specific:
        nm = mpl.colors.TwoSlopeNorm(vmin=z.min(), vcenter=np.median(z), vmax=z.max())
    else: nm = 'linear'
    zmin = z.min()
    zmax = z.max()
    if z2 is not None: 
        zmin = min(zmin, z2.min())
        zmax = max(zmax, z2.max())
    pos = ax.pcolormesh(xx, yy, z, shading='gouraud', norm=nm, vmin=zmin, vmax=zmax)
    cbar = fig.colorbar(pos, ax=ax)
    cbar.ax.tick_params(labelsize=11)
    # ax.ticklabel_format(axis='y', scilimits=[-2,3], useMathText=False)
    ax.tick_params(axis='both', which='major', direction='inout', length=6, labelsize=11)
    ax.tick_params(axis='both', which='minor', direction='inout', length=3)
    ax2x = ax.secondary_xaxis('top', zorder=3)
    ax2x.tick_params(direction='in', which='major', length=3)
    ax2x.tick_params(direction='in', which='minor', length=1.5)
    ax2x.xaxis.set_ticklabels([])
    ax2y = ax.secondary_yaxis('right', zorder=3)
    ax2y.tick_params(direction='in', which='major', length=3)
    ax2y.tick_params(direction='in', which='minor', length=1.5)
    ax2y.yaxis.set_ticklabels([])
    if specific: 
        lct = mpl.ticker.FixedLocator(np.arange(0.2, 0.71, 0.1)) if hrt \
            else mpl.ticker.FixedLocator(np.percentile(z, np.linspace(0,100,7)))
        levels = 7
    else:
        lct = mpl.ticker.MaxNLocator(7)
        levels = lct.tick_values(zmin, zmax)
    if z2 is not None: 
        cs2 = ax.contour(xx, yy, z2, 
                        colors='#90918e', origin='lower', 
                        linestyles='dashed', linewidths=1, 
                        extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]),
                        # locator=lct
                        levels=levels,
                        )
        ax.clabel(cs2, cs2.levels, inline=True, fontsize=10)   
    cs = ax.contour(xx, yy, z, 
                    colors='white', origin='lower', 
                    linestyles='dashed', linewidths=1, 
                    extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]),
                    # locator=lct
                    levels=levels,
                    )
    ax.clabel(cs, cs.levels, inline=True, fontsize=10)
     
    if baseline:
        ax.plot(*baseline, marker='^', ms=7, mfc='white', mec='black', mew=0.5)
    fig.savefig(ospath.join(figures_path, save_as), dpi=300, transparent=True)

def mapping(data=None, data2=None, n=20, reactor_type='PB'):
    if data is None:
        if reactor_type == 'PB':
            data = load_data(ospath.join(results_path, f'optimized_{reactor_type}_45.xlsx'),
                             header=[0,1], skiprows=[2,])
            data2 = load_data(ospath.join(results_path, f'optimized_{reactor_type}_35.xlsx'),
                             header=[0,1], skiprows=[2,])
        else:
            data = load_data(ospath.join(results_path, f'optimized_{reactor_type}.xlsx'),
                             header=[0,1], skiprows=[2,])
    sys = create_system(reactor_type=reactor_type)
    mdl = create_model(sys, kind='mapping')
    # bl = [p.baseline for p in mdl.parameters]
    opt = create_model(sys, kind='optimize')
    xx, yy, zzs = togrid(data, mdl, n)
    if data2 is not None: 
        xx2, yy2, zzs2 = togrid(data2, mdl, n)
    for i, m in enumerate((*opt.parameters, *mdl.metrics)):
        z = zzs[i]
        z2 = None if data2 is None else zzs2[i]
        file = f'heatmaps/{reactor_type}/{m.name}.png'
        plot_heatmap(xx, yy, z, z2, save_as=file)
    if reactor_type == 'FB':  z = zzs[1]*50/(1-zzs[0])*zzs[0]
    else: 
        z = zzs[0]*50*(1-0.45)/0.45
        z2 = zzs2[0]*50*(1-0.35)/0.35 if data2 is not None else None
    file = f'heatmaps/{reactor_type}/Total bead volume.png'    
    plot_heatmap(xx, yy, z, z2, save_as=file)

def _map_diff_hrt(mdl, data=None, n=10):
    if data is None:
        data = load_data(ospath.join(results_path, 'PB_diffusivity_HRT.xlsx'),
                         header=[0,1], skiprows=[2,])
    xx, yy, zzs = togrid(data, mdl, n)
    for z, m in zip(zzs, mdl.metrics):
        file = f'heatmaps/PB_diff_hrt/{m.name}.png'
        plot_heatmap(xx, yy, z, save_as=file)
    
    

#%%
if __name__ == '__main__':
    # path = ospath.join(data_path, 'analysis_framework.xlsx')
    # run_discrete_DVs(path)
    # data = data_compile()
    # plot_clusters(partial=True)
    # plot_clusters(partial=False)
    # out = compare_DVs()
    # plot_diff()
    # llc, imp = best_breakdown()
    # plot_breakdown()
    # dt = load_data(ospath.join(results_path, 'table_compiled.xlsx'), nrows=3553)
    # dt = dt[dt.loc[:,'Reactor type'] != 'UASB']
    # encap_out = compare_DVs(dt, save_as='FB_vs_PB.xlsx')
    # smp = run_UA_SA(N=1000)
    # _rerun_failed_samples(187, 'PB')
    # breakdown_uasa(965)
    plot_univariate_kdes(965)
    # out = calc_3way_diff(965)
    # outs = MCF_pb_to_fb(965)
    # data = MCF_25_vs_75(965)
    # MCF_bubble_plot(outs)
    # MCF_bubble_plot(data)
    # mapping(suffix='specific')
    # mapping(suffix='common')
    # mapping(reactor_type='FB', n=20)
    # mapping(reactor_type='PB', n=11)
    # out = Spearman_corr(965, True)
    # mdl, smps = run_PB_over_diffusivity_HRT()
    # _map_diff_hrt(mdl)

    
