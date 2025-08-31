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
db = Guest.dblue.HEX

#%%
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["figure.autolayout"] = False
plt.rcParams['mathtext.fontset'] = 'custom'

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

configs=('B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'E2', 'E2P', 'F1', 
        'G1', 'G2', 'G3', 'H1', 'I1', 'I2', 'I3', 'N1', 'N2')
data_handles = {
    'Baseline': ('baseline_unopt_performance', 0),   # file name, sheet
    'Adjusted': ('baseline_opt_performance', 'combined'),
    'UD10': ('UD_opt_performance', '10'),
    # 'UD30': ('UD_opt_performance', '30'),
    'UD100': ('UD_opt_performance', '100'),
    'HA': ('HA_opt_performance', 0),
    'ECS': ('ECS_opt_performance', 0)
    }

#%%
def compile_opex():
    opex = []
    keys = []
    for k, v in data_handles.items():
        file, sheet = v
        df = load_data(
            ospath.join(results_path, f'{file}.xlsx'), sheet=sheet,
            header=[0,1], skiprows=[2,]
            )
        df = df.loc[:,['OPEX']]
        opex.append(df)
        keys.append(k)
    
    opex = pd.concat(opex, keys=keys)
    opex.columns = [i[1].split(' cost')[0] for i in opex.columns]
    opex['Lime stabilization'] = opex.loc[:,['Lime', 'Lime stablization energy']].sum(axis=1)
    opex.drop(columns=['Lime stablization energy', 'Total OPEX [USD/d]'], inplace=True)
    return opex

# %%
MGD2cmd = 3785.412

def single_scenario_stacked_bar(opex=None, scenario='Baseline', save_as=''):
    plt.rcParams['xtick.minor.visible'] = False
    plt.rcParams['ytick.minor.visible'] = True
    if opex is None: 
        opex = compile_opex()
        opex = opex.loc[scenario] / (10*MGD2cmd)
    
    fig, ax = plt.subplots(figsize=(6.5, 4))
    handles = []

    x = np.arange(opex.shape[0])
    y_offset = np.zeros(len(x))
    for cat, style in patch_dct.items():
        c, hat = style
        y = opex.loc[:,cat].to_numpy()
        patch = ax.bar(x, y, bottom=y_offset, width=0.65, 
                       color=c, hatch=hat, hatch_linewidth=0.5)
        handles.append(patch)
        y_offset += y

    ax.set_xlim((-0.75, opex.shape[0]-0.25))
    ax.set_xticks(x, opex.index.values)
    ax.set_xlabel('WRRF configuration', fontname='Arial', weight='bold', fontsize=12)
    ax.set_ylabel('$\mathbf{OPEX}$ [USD·m${^{-3}}$]', fontname='Arial', fontsize=12)
    ax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=11)
    ax.tick_params(axis='y', which='minor', direction='inout', length=5)
    ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=11)
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    fig.legend(handles=handles, labels=categories, reverse=True,
               bbox_to_anchor=(0.163, 0.93), loc='upper left', 
               framealpha=1, fontsize=11)
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
    if save_as:
        fig.savefig(ospath.join(figures_path, 'baseline_vs_gpsx/'+save_as),
                    dpi=300, transparent=True)
    else:
        return fig, ax
#%%
def stacked_bar(opex=None, save_as=''):
    plt.rcParams['xtick.minor.visible'] = False
    plt.rcParams['ytick.minor.visible'] = True
    if opex is None:
        opex = compile_opex()
    fig, axes = plt.subplots(ncols=len(configs), sharey=True, figsize=(22, 5))
    # fig, axes = plt.subplots(ncols=len(configs), sharey=True, figsize=(12, 5))
    handles = []
    for ID, ax in zip(configs, axes):
        df = opex.loc[(slice(None), ID),:]
        df.index = df.index.droplevel(1)
        for i in data_handles.keys():
            if i not in df.index:
                df.loc[i,:] = 0
        x = np.arange(df.shape[0])
        y_offset = np.zeros(len(x))
        for cat, style in patch_dct.items():
            c, hat = style
            y = df.loc[:,cat].to_numpy()
            patch = ax.bar(x, y, bottom=y_offset, width=0.65, 
                           color=c, hatch=hat, hatch_linewidth=0.5,
                           tick_label=df.index.values)
            if ID == 'B1': handles.append(patch)
            y_offset += y
        xx = x[y_offset == 0]
        if xx: ax.scatter(xx, 175, s=12, c='black', marker='x', linewidths=0.75)
        ax.set_title(ID, fontsize=14, fontweight='bold')
        ax.set_xlim((-0.75, df.shape[0]-0.25))
        if ID == 'B1': 
            ax.set_ylabel('$\mathbf{OPEX}$ [USD·day${^{-1}}$]', fontname='Arial', fontsize=14)
        ax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=12)
        ax.tick_params(axis='y', which='minor', direction='inout', length=5)
        ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=11, labelrotation=90)
        if ID == 'N2':
            ax2y = ax.secondary_yaxis('right')
            ax2y.tick_params(axis='y', which='major', direction='in', length=4)
            ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
            ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    fig.legend(handles=handles, labels=categories, reverse=True,
               bbox_to_anchor=(0.135, 0.85), loc='upper left', 
               framealpha=1, fontsize=12)
    fig.subplots_adjust(wspace=0, bottom=0.15)
    if save_as:
        fig.savefig(ospath.join(figures_path, save_as),
                    dpi=300, transparent=True)
    else:
        return fig, axes

#%%
def horizontal_stacked_bar(opex=None, save_as=''):
    plt.rcParams['ytick.minor.visible'] = False
    plt.rcParams['xtick.minor.visible'] = True
    if opex is None:
        opex = compile_opex()
    fig, axes = plt.subplots(nrows=len(configs), sharex=True, figsize=(6, 22))
    handles = []

    for ID, ax in zip(configs, axes):
        df = opex.loc[(slice(None), ID),:]
        df.index = df.index.droplevel(1)
        for i in data_handles.keys():
            if i not in df.index:
                df.loc[i,:] = 0
        x = np.arange(df.shape[0])[::-1]
        y_offset = np.zeros(len(x))
        for cat, style in patch_dct.items():
            c, hat = style
            y = df.loc[:,cat].to_numpy()
            patch = ax.barh(x, y, left=y_offset, height=0.65, 
                            color=c, hatch=hat, hatch_linewidth=0.5,
                            tick_label=df.index.values)
            if ID == 'B1': handles.append(patch)
            y_offset += y
        xx = x[y_offset == 0]
        if len(xx): ax.scatter(175, xx, s=12, c='black', marker='x', linewidths=0.75)
        ax.set_title(ID+'  ', fontsize=14, fontweight='bold',
                     loc='right', y=0.4, pad=0)
        ax.set_ylim((-0.75, df.shape[0]-0.25))
        if ID == 'N2': 
            ax.set_xlabel('$\mathbf{OPEX}$ [USD·day${^{-1}}$]', fontname='Arial', fontsize=14)
        ax.tick_params(axis='x', which='major', direction='inout', length=8, labelsize=12)
        ax.tick_params(axis='x', which='minor', direction='inout', length=5)
        ax.tick_params(axis='y', which='major', direction='out', length=4, labelsize=11)
        if ID == 'B1':
            ax2x = ax.secondary_xaxis('top')
            ax2x.tick_params(axis='x', which='major', direction='in', length=4)
            ax2x.tick_params(axis='x', which='minor', direction='in', length=2.5)
            ax2x.xaxis.set_major_formatter(plt.NullFormatter())
    
    fig.legend(handles=handles, labels=categories, reverse=True,
               bbox_to_anchor=(0.455, 0.88), 
               loc='upper left', framealpha=1, fontsize=12)
    fig.subplots_adjust(hspace=0, left=0.15)
    if save_as:
        fig.savefig(ospath.join(figures_path, save_as),
                    dpi=300, transparent=True)
    else:
        return fig, axes


#%%
if __name__ == '__main__':
    opex = compile_opex()
    stacked_bar(opex, save_as='opex.png')
    # stacked_bar(opex, save_as='opex_wrrfs.png')
    horizontal_stacked_bar(opex, save_as='opexh.png')
    single_scenario_stacked_bar(opex, save_as='opex.tif')
