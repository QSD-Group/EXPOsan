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
import qsdsan as qs, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = False
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

#%%
def run_discrete_DVs(samples_path):
    qs.PowerUtility.price = 0.0913
    dct = load_data(samples_path, sheet=None)
    # for i in ('UASB','FB','PB',):
    for i in ('FB',):
        for n, j in (
                (1,'P'), 
                (2,'P'), (2,'M'), (2,'H'),
                ):
            sys = create_system(n_stages=n, reactor_type=i, gas_extraction=j)
            print(sys.ID)
            mdl = create_model(sys, kind='DV', exception_hook='warn')
            sample = dct[i+'_'+j].to_numpy()
            run_model(mdl, sample)

#%%
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
        sizes=(30, 45),
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
        ax.set_ylim(-100, 500)
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
    
#%%
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
    

#%%
boxprops = dict(alpha=0.7, edgecolor='black')
flierprops = dict(marker='.', markersize=1, markerfacecolor='#90918e', markeredgecolor='#90918e')
meanprops = dict(marker='^', markersize=2.5, markerfacecolor='black', markeredgecolor='black')
medianprops = dict(color='black', lw=1)

pal = ['#60C1CF', '#F98F60', '#79BF82', '#F3C354', '#A280B9', '#ED586F', 
       '#35767F', '#733763', '#4D7E53', '#AB8937', '#5184EF']

def plot_joint(df, save_as='', kde=True):
    g = sns.JointGrid(height=8, ratio=5, space=0, marginal_ticks=True)
    x, y = df['Levelized cost'], df['GWP100']
    group = df['group']
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
    g.ax_joint.tick_params(axis='both', which='major', direction='inout', length=10, labelsize=18)
    g.ax_joint.tick_params(axis='both', which='minor', direction='inout', length=6)
    g.ax_joint.set_xlabel('')
    g.ax_joint.set_ylabel('')
    bxp_kwargs = dict(
        palette=pal,
        linewidth=1,
        showcaps=True,
        showmeans=True,
        showfliers=False,
        dodge=False,
        saturation=1,
        width=0.3,
        boxprops=boxprops,
        flierprops=flierprops,
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
    

#%%
if __name__ == '__main__':
    # path = ospath.join(data_path, 'analysis_framework.xlsx')
    # run_discrete_DVs(path)
    # data = data_compile()
    # plot_clusters(partial=True)
    # plot_clusters(partial=False)
    # out = compare_DVs()
    plot_diff()
