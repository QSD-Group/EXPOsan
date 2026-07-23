# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os
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

#%%
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["figure.autolayout"] = False
plt.rcParams['mathtext.fontset'] = 'custom'

patch_dct = {
    'Aeration energy': (b, ''),
    'Pumping energy': (p, '-----'),
    'Mixing energy': (y, r'\\\\\\'),
    'External carbon': (r, ''),
    'Coagulant': (o, '/////'),
    'Lime stabilization': (g, '|||||'),
    'Sludge disposal': (a, ''),
    'AD heating': (o, ''),
    'CHP heat recovery': (g, ''),
    'CHP electricity recovery': (b, '/////'),
    }
categories = [
    'Aeration',
    'Pumping',
    'Mechanical mixing',
    'External carbon',
    'Coagulant',
    'Lime stabilization',
    'Sludge disposal',
    'Heating',
    'CHP heat recovery',
    'CHP electricity recovery',
    ]

configs = (
    'B1', 'B1E', 'B2', 'B3',
    'C1', 'C1E', 'C2', 'C3',
    'E2', 'E2P',
    'F1', 'F1E',
    'G1', 'G1E', 'G2', 'G3',
    'H1', 'H1E',
    'I1', 'I1E', 'I2', 'I3',
    'N1', 'N1E', 'N2',
    )
data_handles = {
    'Baseline': ('_baseline_unopt_performance_w_chp', 0),   # file name, sheet
    'Adjusted': ('_baseline_opt_performance_w_chp', 0),
    }

MGD2cmd = 3785.412

#%%
def save_figure(fig, save_as):
    path = ospath.join(figures_path, save_as)
    os.makedirs(ospath.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=300, transparent=True)

#%%
def clean_opex_column(col):
    feature = col[1] if isinstance(col, tuple) else col
    feature = str(feature)
    if ' cost' in feature:
        feature = feature.split(' cost')[0]
    feature = feature.replace(' [USD/d]', '')
    feature = feature.replace('CHP heat credit', 'CHP heat recovery')
    feature = feature.replace('CHP electricity credit', 'CHP electricity recovery')
    feature = feature.replace('CHP total credit', 'CHP total recovery')
    return feature

def compile_opex(normalize_by_flow=True):
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
    opex.columns = [clean_opex_column(i) for i in opex.columns]
    opex['Lime stabilization'] = opex.loc[:,['Lime', 'Lime stablization energy']].sum(axis=1)
    opex.drop(
        columns=[
            'Lime stablization energy',
            'Total OPEX',
            'CHP total recovery',
            ],
        inplace=True,
        errors='ignore',
        )
    for cat in patch_dct:
        if cat not in opex.columns:
            opex[cat] = 0.
    opex = opex.fillna(0.)
    if normalize_by_flow:
        opex /= (10*MGD2cmd)
    return opex

# %%
def stacked_bar(opex=None, save_as=''):
    plt.rcParams['xtick.minor.visible'] = False
    plt.rcParams['ytick.minor.visible'] = True
    if opex is None:
        opex = compile_opex()
    fig, axes = plt.subplots(ncols=len(configs), sharey=True, figsize=(24, 5))
    handles = []

    global_min = 0.0
    global_max = 0.0

    for ID, ax in zip(configs, axes):
        df = opex.loc[(slice(None), ID),:]
        df.index = df.index.droplevel(1)
        for scenario in data_handles:
            if scenario not in df.index:
                df.loc[scenario,:] = 0
        df = df.loc[list(data_handles)]
        x = np.arange(df.shape[0])
        y_offset_pos = np.zeros(len(x), dtype=float)
        y_offset_neg = np.zeros(len(x), dtype=float)
        for cat, style in patch_dct.items():
            c, hat = style
            yvals = df.loc[:,cat].to_numpy()
            y_pos = np.where(yvals >= 0, yvals, 0.0)
            y_neg = np.where(yvals < 0, yvals, 0.0)

            patch_pos = ax.bar(
                x, y_pos, bottom=y_offset_pos, width=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values,
                )
            patch_neg = ax.bar(
                x, y_neg, bottom=y_offset_neg, width=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values,
                )
            if ID == configs[0]:
                if np.any(y_pos):
                    handles.append(patch_pos)
                else:
                    handles.append(patch_neg)
            y_offset_pos += y_pos
            y_offset_neg += y_neg
        y_total = y_offset_pos + y_offset_neg
        ax.bar(x, y_total, width=0.65, ec='black', fc=(1,1,1,0), lw=0.8)
        ax.axhline(0, color='black', lw=0.8)
        xx = x[y_total == 0]
        if len(xx): ax.scatter(xx, 0.01, s=12, c='black', marker='x', linewidths=0.75)
        ax.set_title(ID, fontsize=13, fontweight='bold')
        ax.set_xlim((-0.75, df.shape[0]-0.25))
        if ID == configs[0]:
            ax.set_ylabel(r'$\mathbf{OPEX}$ [USD·m$^{-3}$]', fontname='Arial', fontsize=14)
        ax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=12)
        ax.tick_params(axis='y', which='minor', direction='inout', length=5)
        ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=11, labelrotation=90)
        if ID == configs[-1]:
            ax2y = ax.secondary_yaxis('right')
            ax2y.tick_params(axis='y', which='major', direction='in', length=4)
            ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
            ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        global_min = min(global_min, y_offset_neg.min(), y_total.min())
        global_max = max(global_max, y_offset_pos.max(), y_total.max())

    span = global_max - global_min
    pad = 0.025 * span if span > 0 else 0.1
    for ax in axes:
        ax.set_ylim(global_min - pad, global_max + pad)

    fig.legend(handles=handles, labels=categories, reverse=True,
               bbox_to_anchor=(0.06, 0.9), loc='upper left',
               framealpha=1, fontsize=12)
    fig.subplots_adjust(wspace=0, bottom=0.16, top=0.94, left=0.045, right=0.99)
    if save_as:
        save_figure(fig, save_as)
    else:
        return fig, axes

#%%
def horizontal_stacked_bar(opex=None, save_as=''):
    plt.rcParams['ytick.minor.visible'] = False
    plt.rcParams['xtick.minor.visible'] = True
    if opex is None:
        opex = compile_opex()
    fig, axes = plt.subplots(nrows=len(configs), sharex=True, figsize=(7, 22))
    handles = []

    global_min = 0.0
    global_max = 0.0

    for ID, ax in zip(configs, axes):
        df = opex.loc[(slice(None), ID),:]
        df.index = df.index.droplevel(1)
        for scenario in data_handles:
            if scenario not in df.index:
                df.loc[scenario,:] = 0
        df = df.loc[list(data_handles)]
        x = np.arange(df.shape[0])[::-1]
        y_offset_pos = np.zeros(len(x), dtype=float)
        y_offset_neg = np.zeros(len(x), dtype=float)
        for cat, style in patch_dct.items():
            c, hat = style
            yvals = df.loc[:,cat].to_numpy()
            y_pos = np.where(yvals >= 0, yvals, 0.0)
            y_neg = np.where(yvals < 0, yvals, 0.0)

            patch_pos = ax.barh(
                x, y_pos, left=y_offset_pos, height=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values,
                )
            patch_neg = ax.barh(
                x, y_neg, left=y_offset_neg, height=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values,
                )
            if ID == configs[0]:
                if np.any(y_pos):
                    handles.append(patch_pos)
                else:
                    handles.append(patch_neg)
            y_offset_pos += y_pos
            y_offset_neg += y_neg
        y_total = y_offset_pos + y_offset_neg
        ax.barh(x, y_total, height=0.65, ec='black', fc=(1,1,1,0), lw=0.8)
        ax.axvline(0, color='black', lw=0.8)
        xx = x[y_total == 0]
        if len(xx): ax.scatter(0.01, xx, s=12, c='black', marker='x', linewidths=0.75)
        title = ax.set_title(ID, fontsize=14, fontweight='bold',
                             loc='right', x=1.08, y=0.4, pad=0)
        title.set_clip_on(False)
        ax.set_ylim((-0.75, df.shape[0]-0.25))
        if ID == configs[-1]:
            ax.set_xlabel(r'$\mathbf{OPEX}$ [USD·m${^{-3}}$]', fontname='Arial', fontsize=14)
        ax.tick_params(axis='x', which='major', direction='inout', length=8, labelsize=12)
        ax.tick_params(axis='x', which='minor', direction='inout', length=5)
        ax.tick_params(axis='y', which='major', direction='out', length=4, labelsize=11)
        if ID == configs[0]:
            ax2x = ax.secondary_xaxis('top')
            ax2x.tick_params(axis='x', which='major', direction='in', length=4)
            ax2x.tick_params(axis='x', which='minor', direction='in', length=2.5)
            ax2x.xaxis.set_major_formatter(plt.NullFormatter())
        global_min = min(global_min, y_offset_neg.min(), y_total.min())
        global_max = max(global_max, y_offset_pos.max(), y_total.max())

    span = global_max - global_min
    pad = 0.025 * span if span > 0 else 0.1
    for ax in axes:
        ax.set_xlim(global_min - pad, global_max + pad)

    fig.legend(handles=handles, labels=categories, reverse=True,
               bbox_to_anchor=(0.5, 0.975),
               loc='upper left', framealpha=1, fontsize=12)
    fig.subplots_adjust(hspace=0, left=0.16, right=0.90, bottom=0.04, top=0.99)
    if save_as:
        save_figure(fig, save_as)
    else:
        return fig, axes

#%%
if __name__ == '__main__':
    opex = compile_opex()
    stacked_bar(opex, save_as='arpae_opex_w_chp.tif')
    horizontal_stacked_bar(opex, save_as='arpae_opexh_w_chp.tif')
