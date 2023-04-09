# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan.utils import ospath, load_data
from exposan.metab import create_system, create_model, run_model, \
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
        sys = create_system(n_stages=1, reactor_type=rt, gas_extraction='P')
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

    sys = create_system(n_stages=1, reactor_type=rt, gas_extraction='P')
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


#%% discrete-DV scatter plot
def data_compile(save=True):
    dfs = []
    for i in ('UASB','FB','PB',):
        for n, j in ((1,'P'), (2,'P'), (2,'M'), (2,'H'),):
            df = load_data(ospath.join(results_path, f'{i}{n}{j}.xlsx'), 
                           header=[0,1], skiprows=[2,])
            df.drop(columns=['Biomass','Process'], level='Element', inplace=True)
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
    pal = {'UASB':'#60c1cf', 'FB':'#F98F60', 'PB':'#a280b9'}
    edge = ['#000000' if dg else pal[rct] for dg, rct in 
            data.loc[:,['Effluent degassing','Reactor type']].to_numpy()]
    leg = False if partial else 'auto'
    ax = sns.scatterplot(
        data=data,
        x='Levelized cost',
        y='GWP100',
        hue='Reactor type',
        palette=pal,
        size='Number of stages',
        sizes=(30, 55),
        style='Gas extraction',
        markers={'H':'s', 'M':'^', 'P':'o'},
        edgecolor=edge,
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
    fig.savefig(path, dpi=300, facecolor='white')
    
    return fig, ax
    
#%% compute diff by discrete-DV
def calc_diff(df, factor, baseline_value, norms):    
    bl = df.xs(key=baseline_value, axis=1, level=factor)
    df.drop(labels=baseline_value, axis=1, level=factor, inplace=True)
    df.index = range(df.shape[0])
    out = []
    for i in ('Levelized cost', 'GWP100'):
        diff = df.xs(i, axis=1)
        diff = diff.sub(bl[i].values, axis=0)/norms[i]
        diff = pd.melt(diff, value_name=i, ignore_index=False)
        diff['id'] = diff.index
        diff = diff.set_index([c for c in diff.columns if c != i])
        out.append(diff)
    out = out[0].join(out[1])
    out['pair'] = out.index.get_level_values(factor)   
    out.index = range(out.shape[0])
    return out

groups = (
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
    data.drop('id', axis=1, inplace=True)
    vals = ['Levelized cost', 'GWP100']
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
    path = save_as or ospath.join(results_path, 'diff.xlsx')
    with pd.ExcelWriter(path) as writer:
        for k,v in out.items():
            v.to_excel(writer, sheet_name=k)
    return out  
    

#%% contour of diff
boxprops = dict(alpha=1, edgecolor='black')
# flierprops = dict(marker='.', markersize=1, markerfacecolor='#90918e', markeredgecolor='#90918e')
meanprops = dict(marker='^', markersize=1.5, markerfacecolor='black', markeredgecolor='black', lw=0.5)
medianprops = dict(color='black', lw=0.5)

palette = ['#60C1CF', '#F98F60', '#79BF82', '#F3C354', '#A280B9', '#ED586F', 
           '#35767F', '#733763', '#4D7E53', '#AB8937', '#5184EF']

def plot_joint(df, save_as='', kde=True):
    g = sns.JointGrid(height=8, ratio=5, space=0, marginal_ticks=True)
    x, y = df['Levelized cost'], df['GWP100']
    group = df['group']
    pal = palette[:group.nunique()]
    g.refline(x=0, y=0, color='#90918e', lw=1)
    if kde:
        sns.kdeplot(x=x, y=y, hue=group, ax=g.ax_joint,
                    palette=pal,
                    common_norm=False,
                    fill=True,
                    legend=False,
                    alpha=0.6,
                    )
    else:
        sns.scatterplot(x=x, y=y, hue=group, ax=g.ax_joint,
                        palette=pal,
                        legend=False,
                        size=30,
                        alpha=0.7,
                        )
    # g.ax_joint.set_ylim(top=3000)
    g.ax_joint.tick_params(axis='both', which='major', direction='inout', length=10, labelsize=18)
    g.ax_joint.tick_params(axis='both', which='minor', direction='inout', length=6)
    g.ax_joint.set_xlabel('')
    g.ax_joint.set_ylabel('')
    bxp_kwargs = dict(
        palette=pal,
        linewidth=0.5,
        showcaps=True,
        showmeans=True,
        showfliers=False,
        dodge=False,
        saturation=1,
        width=0.3,
        boxprops=boxprops,
        # flierprops=flierprops,
        meanprops=meanprops,
        medianprops=medianprops,
        )
    if group.nunique() > 3: bxp_kwargs.pop('width')
    sns.boxplot(x=x, y=group, hue=group, ax=g.ax_marg_x, **bxp_kwargs)
    sns.boxplot(y=y, x=group, hue=group, ax=g.ax_marg_y, **bxp_kwargs)
    for ax in (g.ax_marg_x, g.ax_marg_y): ax.legend_.remove()
    g.ax_marg_x.tick_params(axis='x', which='major', direction='out', length=5)
    g.ax_marg_x.tick_params(axis='x', which='minor', direction='out', length=3)    
    g.ax_marg_x.tick_params(left=False, which='both', labelleft=False)    
    g.ax_marg_x.spines['left'].set_color('white')
    # g.ax_marg_y.set_ylim(top=3000)
    g.ax_marg_y.tick_params(axis='y', which='major', direction='out', length=5)
    g.ax_marg_y.tick_params(axis='y', which='minor', direction='out', length=3) 
    g.ax_marg_y.tick_params(bottom=False, which='both', labelbottom=False)    
    g.ax_marg_y.spines['bottom'].set_color('white')
    
    path = ospath.join(figures_path, save_as)
    g.savefig(path, dpi=300, facecolor='white')

def plot_diff(data=None):
    if data is None:
        data = load_data(ospath.join(results_path, 'diff.xlsx'), sheet=None)
    singles = ('Reactor type', 'Total HRT', 'Bead lifetime')
    for i in singles:
        df = data[i]
        df['group'] = df['pair'].astype(str)
        plot_joint(df, f'{i}.png', kde=(i != 'Bead lifetime'))
    others = []
    for k, v in data.items():
        if k not in singles:
            v['group'] = v.apply(lambda row: f"{k}_{row['pair']}", axis=1)
            others.append(v)
    df = pd.concat(others)
    plot_joint(df, 'other_DVs.png')
    
#%% cost & impact breakdowns
def best_breakdown():
    kwargs = dict(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
    llc_bd = {}
    imp_bd = {}
    
    # UASB
    mdl1 = create_model()
    sys1 = mdl1.system
    sub1, = sys1.subsystems
    
    # best at TEA
    smp = np.array([1, 22])
    for p, v in zip(mdl1.parameters, smp): p.setter(v)
    sys1.simulate(**kwargs)
    llc_bd['uasb_t1'] = categorize_cashflow(sub1.TEA)
    imp_bd['uasb_t1'] = categorize_all_impacts(sub1.LCA)
    # best at LCA
    smp = np.array([2, 22])
    for p, v in zip(mdl1.parameters, smp): p.setter(v)
    sys1.simulate(**kwargs)
    llc_bd['uasb_l1'] = categorize_cashflow(sub1.TEA)
    imp_bd['uasb_l1'] = categorize_all_impacts(sub1.LCA)
    
    # FB
    mdl2 = create_model(reactor_type='FB')
    sys2 = mdl2.system
    sys2.units[0].bead_lifetime = 30
    sub2, = sys2.subsystems
    
    # best at TEA
    smp = np.array([2, 0.9, 2, 22])
    for p, v in zip(mdl2.parameters, smp): p.setter(v)
    sys2.simulate(**kwargs)
    llc_bd['fb_t1'] = categorize_cashflow(sub2.TEA)
    imp_bd['fb_t1'] = categorize_all_impacts(sub2.LCA)
    # best at LCA
    smp = np.array([2, 0.75, 1, 22])
    for p, v in zip(mdl2.parameters, smp): p.setter(v)
    sys2.simulate(**kwargs)
    llc_bd['fb_l1'] = categorize_cashflow(sub2.TEA)
    imp_bd['fb_l1'] = categorize_all_impacts(sub2.LCA)
    
    # PB
    mdl3 = create_model(reactor_type='PB')
    sys3 = mdl3.system
    sys3.units[0].bead_lifetime = 30
    sub3, = sys3.subsystems
    
    smp = np.array([2, 1, 22])
    for p, v in zip(mdl3.parameters, smp): p.setter(v)
    sys3.simulate(**kwargs)
    # best at TEA
    llc_bd['pb_t1'] = categorize_cashflow(sub3.TEA)
    imp_bd['pb_t1'] = categorize_all_impacts(sub3.LCA)
    # best at LCA
    # llc_bd['pb_l1'] = categorize_cashflow(sys3.TEA)
    # imp_bd['pb_l1'] = categorize_all_impacts(sys3.LCA)
    
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

    ax.tick_params(axis='y', which='major', direction='inout', length=10, labelsize=12)
    ax.tick_params(axis='y', which='minor', direction='inout', length=6)
    ax.ticklabel_format(axis='y', scilimits=[-2,3], useMathText=True)
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

def MCF_bubble_plot(data=None, seed=None):
    if data is None: 
        data = MCF_encap_to_susp(seed, False)
    mpl.rcParams['xtick.minor.visible'] = False
    mpl.rcParams['ytick.minor.visible'] = False
    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(6,9))
    pal = {'FB':'#F98F60', 'PB':'#a280b9'}
    for kv, ax in zip(data.items(), axes):
        k, df = kv
        sns.scatterplot(df, x='Metric', y='Parameter', 
                        hue=df.p<0.05, palette=['black', pal[k]], 
                        size='D', sizes=(0, 350), size_norm=(0,1),
                        legend=False,
                        ax=ax,
                        )
        ax.set_xlim(-0.5, 4.5)
        ax.set_axisbelow(True)
        ax.grid(True, which='major', color='k', linestyle='--', 
                linewidth=0.7, alpha=0.3)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='both', which='major', direction='inout', length=8)
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.yaxis.set_major_formatter(plt.NullFormatter())

    fig.subplots_adjust(wspace=0)
    fig.savefig(ospath.join(figures_path, 'MCF.png'), dpi=300, transparent=True)

#%% pair-wise comparisons, three-way
def calc_3way_diff(seed, save=True):
    data = {}
    for i in ('UASB', 'FB', 'PB'):
        data[i] = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                            header=[0,1], skiprows=[2,])
    pair_cols = [
            ('Process', 'COD removal [%]'),
            ('Biogas', 'H2 production [kg/d]'),
            ('Biogas', 'CH4 production [kg/d]'),
            ('Biomass', 'R1 Overall TSS [g/L]'),
            # ('Biomass', 'R1 Encapsulated TSS [g/L]'),
            ('TEA (w/ degas)', 'Levelized cost (w/ degas) [$/ton rCOD]'),
            ('LCA (w/ degas)', 'GWP100 (w/ degas) [kg CO2eq/ton rCOD]'),
            ('TEA (w/o degas)', 'Levelized cost (w/o degas) [$/ton rCOD]'),
            ('LCA (w/o degas)',  'GWP100 (w/o degas) [kg CO2eq/ton rCOD]'),    
        ]
    _col = ('Biomass', 'R1 Encapsulated TSS [g/L]')
    out = {}
    out['pu'] = (data['PB'].loc[:, pair_cols] - data['UASB'].loc[:, pair_cols])
    out['fu'] = (data['FB'].loc[:, pair_cols] - data['UASB'].loc[:, pair_cols])
    pair_cols.insert(4, _col)
    out['pu'][_col] = out['fu'][_col] = None
    out['pu'] = out['pu'].loc[:, pair_cols]
    out['fu'] = out['fu'].loc[:, pair_cols]
    out['pf'] = (data['PB'].loc[:, pair_cols] - data['FB'].loc[:, pair_cols])
    
    if save:
        with pd.ExcelWriter(ospath.join(results_path, 'diff_3ways.xlsx')) as writer:
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
    tea.columns = [col.split(' ')[1] for col in tea.columns]
    tea['fug_ch4'] = 0
    tea['total'] = llc
    
    lca = data.loc[:,'LCA (w/ degas)'].copy()
    lca.sort_values(by='GWP100 (w/ degas) [kg CO2eq/ton rCOD]', inplace=True)
    gwp = lca.pop('GWP100 (w/ degas) [kg CO2eq/ton rCOD]')
    # absolute = (gwp<=0).any()
    # if absolute: 
    lca = (lca / 100).mul(gwp.abs(), axis='rows')
    lca.columns = [col.split(' ')[1] for col in lca.columns]
    lca['total'] = gwp
    return tea, lca#, absolute

def plot_area(df, absolute=False):
    fig, ax = plt.subplots(figsize=(4,4))
    x = range(df.shape[0])
    ax.axhline(y=0, color='black', linewidth=0.5)
    yp = np.zeros(df.shape[0])
    yn = yp.copy()
    for k, v in patch_dct.items():
        c, hat = v
        y = df.loc[:,k]
        y_offset = (y>=0)*yp + (y<0)*yn
        ax.fill_between(x, y+y_offset, y_offset, facecolor=c, hatch=hat, linewidth=0.5)
        yp += (y>=0) * y
        yn += (y<0) * y
    ax.set_xlim(0, df.shape[0])
    ax.tick_params(axis='y', which='major', direction='inout', length=6)
    ax.tick_params(axis='y', which='minor', direction='inout', length=3)
    ax.set_xlabel('')
    ax.set_ylabel('')
    if absolute: 
        # ax.yaxis.set_ticklabels([])
        ax.ticklabel_format(axis='y', scilimits=[-2,3], useMathText=True)
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=3)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=1.5)
        ax2y.yaxis.set_ticklabels([])
    else:
        ax2y = ax.twinx()
        ax2y.plot(x, df['total'], color='black', linewidth=0.5)
        ax2y.tick_params(axis='y', which='major', direction='inout', length=6)
        ax2y.tick_params(axis='y', which='minor', direction='inout', length=3)
    return fig, ax

def breakdown_uasa(seed):
    for i in ('UASB', 'FB', 'PB'):
        data = load_data(ospath.join(results_path, f'{i}1P_{seed}.xlsx'),
                         header=[0,1], skiprows=[2,])
        # tea, lca, absolute = breakdown_and_sort(data)
        tea, lca = breakdown_and_sort(data)
        for suffix, df in (('TEA', tea), ('LCA', lca)):
            # fig, ax = plot_area(df, absolute and suffix=='LCA')
            fig, ax = plot_area(df, True)            
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
                         header=[0,1], skiprows=[2,])
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

out = Spearman_corr(364, True)

#%% heatmaps
# from scipy.interpolate import griddata
from exposan.metab.models import meshgrid_sample

def togrid(df, mdl, n):
    zs = df.iloc[:,-len(mdl.metrics):].to_numpy().T
    us = df.iloc[:, :-(-len(mdl.metrics)+len(mdl.parameters))].to_numpy().T
    samples, xx, yy = meshgrid_sample(*mdl.parameters, n)
    zzs = np.vstack((us, zs))
    zzs = zzs.reshape((zzs.shape[0], *xx.shape))
    return xx, yy, zzs

def plot_heatmap(xx, yy, z, save_as=''):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    pos = ax.pcolormesh(xx, yy, z, shading='gouraud')
    cbar = fig.colorbar(pos, ax=ax)
    cbar.ax.tick_params(labelsize=10)
    ax.tick_params(axis='both', which='major', direction='inout', length=6, labelsize=10)
    ax.tick_params(axis='both', which='minor', direction='inout', length=3)
    ax2x = ax.secondary_xaxis('top', zorder=3)
    ax2x.tick_params(direction='in', which='major', length=3)
    ax2x.tick_params(direction='in', which='minor', length=1.5)
    ax2x.xaxis.set_ticklabels([])
    ax2y = ax.secondary_yaxis('right', zorder=3)
    ax2y.tick_params(direction='in', which='major', length=3)
    ax2y.tick_params(direction='in', which='minor', length=1.5)
    ax2y.yaxis.set_ticklabels([])
    cs = ax.contour(xx, yy, z, levels=8,
                    colors='white', origin='lower', 
                    linestyles='dashed', linewidths=0.6, 
                    extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]),
                    algorithm='serial',
                    )
    ax.clabel(cs, cs.levels, inline=True, 
              # fmt=lambda x: f"{x:.2g}", 
              fontsize=7.5)
    fig.savefig(ospath.join(figures_path, save_as), dpi=300, transparent=True)

def mapping(data=None, n=20, reactor_type='PB', suffix=''):
    if data is None:
        data = load_data(ospath.join(results_path, f'optimized_{reactor_type}.xlsx'),
                         header=[0,1], skiprows=[2,])
    sys = create_system(reactor_type=reactor_type)
    mdl = create_model(sys, kind='mapping')
    opt = create_model(sys, kind='optimize')
    xx, yy, zzs = togrid(data, mdl, n)
    for z, m in zip(zzs, (*opt.parameters, *mdl.metrics)):
        file = f'heatmaps/{reactor_type}/{m.name}_{suffix}.png'
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
    # smp = run_UA_SA(seed=364, N=1000)
    # _rerun_failed_samples(364)
    # data = MCF_encap_to_susp(364, False)
    # MCF_bubble_plot(data)
    # breakdown_uasa(364)
    mapping(suffix='common')
