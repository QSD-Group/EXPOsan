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
    edge = ['#000000' if i else '#FFFFFF' for i in data['Effluent degassing'].to_numpy()]
    ax = sns.scatterplot(
        data=data,
        x='Levelized cost',
        y='GWP100',
        hue='Reactor type',
        palette={'UASB':'#60c1cf', 'FB':'#F98F60', 'PB':'#a280b9'},
        size='Number of stages',
        sizes=(30, 45),
        style='Gas extraction',
        markers={'H':'s', 'M':'^', 'P':'o'},
        edgecolor=edge,
        legend=False,
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
label_diff = lambda row, bl: f'{row}-{bl}'

def calc_diff(df, factor, baseline_value):    
    bl = df.xs(key=baseline_value, axis=1, level=factor)
    df.drop(labels=baseline_value, axis=1, level=factor, inplace=True)
    df.index = range(df.shape[0])
    out = []
    for i in ('Levelized cost', 'GWP100'):
        diff = df.xs(i, axis=1)
        diff -= bl[i].values
        diff = pd.melt(diff, value_name=i, ignore_index=False)
        diff['id'] = diff.index
        diff = diff.set_index([c for c in diff.columns if c != i])
        out.append(diff)
    out = out[0].join(out[1])    
    out['factor'] = factor
    out['pair'] = [label_diff(fct, baseline_value) for fct in out.index.get_level_values(factor)]
    out.index = range(out.shape[0])
    return out

groups = (
    ('Reactor type', ['Bead lifetime','Bead diameter', 'Voidage'], 'UASB'),
    ('Number of stages', [], 1),
    ('Gas extraction', ['Recirculation ratio', 'Headspace pressure'], 'P'),
    ('Temperature', [], 22),
    ('Effluent degassing', [], 0),
    )

def compare_DVs(data=None):
    if data is None:
        data = load_data(ospath.join(results_path, 'table_compiled.xlsx'))
    vals = ['Levelized cost', 'GWP100']
    dvs = set(data.columns) - set(vals)
    out = []
    for factor, covar, bl in groups:
        idx = list(dvs - {factor, *covar})
        cols = [factor, *covar]
        df = data.pivot(index=idx, columns=cols, values=vals)
        # df.dropna(axis=0, inplace=True)
        #!!! not right
        try: out.append(calc_diff(df, factor, bl))
        except: breakpoint()
    out = pd.concat(out)
    return out
        
#%%
if __name__ == '__main__':
    # path = ospath.join(data_path, 'analysis_framework.xlsx')
    # run_discrete_DVs(path)
    # data = data_compile()
    # plot_clusters(partial=True)
    # plot_clusters(partial=False)
    compare_DVs()
