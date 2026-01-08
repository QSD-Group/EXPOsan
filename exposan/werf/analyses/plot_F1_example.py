# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt
from qsdsan.utils import load_data, ospath, palettes
from exposan.werf import results_path, figures_path

import warnings
warnings.filterwarnings("ignore")

Guest = palettes['Guest']

b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX
p = Guest.purple.HEX
a = Guest.gray.HEX
da = Guest.dgray.HEX
db = Guest.dblue.HEX

#%%
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["figure.autolayout"] = False
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

categories = ['Aeration', 'Pumping', 'Mechanical mixing', 
              'External carbon', 'Coagulant', 'Lime stabilization', 
              'Sludge disposal']
patch_dct = {
    'Aeration energy': (b, ''),
    'Pumping energy': (p, '-----'),
    'Mixing energy': (y, r'\\\\\\'),
    'External carbon': (r, ''),
    'Coagulant': (o, '/////'),
    'Lime stabilization': (g, '|||||'),
    'Sludge disposal': (db, ''),
    # 'Sludge disposal': (a, ''),    
    }

configs=('F1', )
data_handles = {
    'Baseline': ('baseline_unopt_performance', 0),   # file name, sheet
    'Adjusted': ('baseline_opt_performance', 'combined'),
    'HA': ('HA_opt_performance', 0),
    'ECS': ('ECS_opt_performance', 0)
    }

#%%
def compile_data():
    opex = []
    keys = []
    recv = pd.DataFrame()
    for k, v in data_handles.items():
        file, sheet = v
        df = load_data(
            ospath.join(results_path, f'{file}.xlsx'), sheet=sheet,
            header=[0,1], skiprows=[2,]
            )
        if 'baseline' in file: 
            recv[k] = 0
        else:
            recv[k] = df[('System','NH4 recovery [kg-N/d]')]
        df = df.loc[:,['OPEX']]
        opex.append(df)
        keys.append(k)
    
    opex = pd.concat(opex, keys=keys)
    opex.columns = [i[1].split(' cost')[0] for i in opex.columns]
    opex['Lime stabilization'] = opex.loc[:,['Lime', 'Lime stablization energy']].sum(axis=1)
    opex.drop(columns=['Lime stablization energy', 'Total OPEX [USD/d]'], inplace=True)
    return opex, recv

#%%
def stacked_bar(opex=None, recv=(), ID='F1'):
    plt.rcParams['xtick.minor.visible'] = False
    plt.rcParams['ytick.minor.visible'] = True
    if opex is None:
        opex, recv = compile_data()
    fig, (tax, ax) = plt.subplots(nrows=2, sharex=True, figsize=(4, 8), height_ratios=[3,5])
    handles = []
    df = opex.loc[(slice(None), ID),:]
    df.index = df.index.droplevel(1)
    tot = df.sum(axis=1).to_numpy()
    dopex = tot/tot[0]-1
    dopex = [f'{d:.1%}' if d != 0 else '' for d in dopex]
    for i in data_handles.keys():
        if i not in df.index:
            df.loc[i,:] = 0
    x = np.arange(df.shape[0])
    tax.bar(x, recv.loc[ID], width=0.65, color=da)
    tax.set_ylabel('$\mathbf{Ammonia\ recovery}$\n[kg-N·day${^{-1}}$]', 
                   fontname='Arial', fontsize=13, linespacing=0.8,)
    tax.tick_params(axis='x', which='major', direction='inout', length=8)
    tax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=12)
    tax.tick_params(axis='y', which='minor', direction='inout', length=5)
    tax2y = tax.secondary_yaxis('right')
    tax2y.tick_params(axis='y', which='major', direction='in', length=4)
    tax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
    tax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    y_offset = np.zeros(len(x))
    for cat, style in patch_dct.items():
        c, hat = style
        y = df.loc[:,cat].to_numpy()
        patch = ax.bar(x, y, bottom=y_offset, width=0.65, 
                       color=c, hatch=hat, hatch_linewidth=0.5,
                       tick_label=['Baseline', 'Adjusted', 'With HA', 'With ECS'])
        handles.append(patch)
        y_offset += y
    for i in range(len(x)): 
        ax.text(x[i], tot[i]+30, dopex[i], fontsize=11, color='black', ha='center')
    ax.hlines(y=tot[0], xmin=x[0]-0.65/2, xmax=x[-1]+0.65/2, ls='--', color=a)
    ax.set_xlim((-0.75, df.shape[0]-0.25))
    ax.set_ylabel('$\mathbf{WWTP\ OPEX}$ [USD·day${^{-1}}$]', fontname='Arial', fontsize=13)
    ax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=12)
    ax.tick_params(axis='y', which='minor', direction='inout', length=5)
    ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=12, labelrotation=45)
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    fig.legend(handles=handles, labels=categories, reverse=True,
               bbox_to_anchor=(0.9, 0.6), loc='upper left', 
               framealpha=1, fontsize=12)
    fig.subplots_adjust(hspace=0)
    # fig.tight_layout()
    fig.savefig(ospath.join(figures_path, f'{ID}_example_opex_bar.png'),
                bbox_inches='tight',
                dpi=300, transparent=True)



#%%
if __name__ == '__main__':
    opex, recv = compile_data()
    stacked_bar(opex, recv)
    # horizontal_stacked_bar(opex, save_as='opexh.png')
