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
              'Sludge disposal','Heating','Heat recovery','Electricity recovery']
patch_dct = {
    'Aeration energy': (b, ''),
    'Pumping energy': (p, '-----'),
    'Mixing energy': (y, r'\\\\\\'),
    'External carbon': (r, ''),
    'Coagulant': (o, '/////'),
    'Lime stabilization': (g, '|||||'),
    # 'Sludge disposal': (db, ''),
    'Sludge disposal': (a, ''),
    'AD heating': (o,''),
    # uncomment the following line for heat recovery confifuration:
    'Heat recovery': (g,''),
    # uncomment the following two lines for CHP confifuration:
    # 'CHP heat recovery': (g,''), 
    # 'CHP electricity recovery': (b,'/////'),
    }

configs=('B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'E2', 'E2P', 'F1', 
        'G1', 'G2', 'G3', 'H1', 'I1', 'I2', 'I3', 'N1', 'N2')
data_handles = {
    'Baseline': ('baseline_unopt_performance', 0),   # file name, sheet
    'Adjusted': ('baseline_opt_performance', 0),
    'UD10': ('UD_opt_performance', '10'),
    # 'UD30': ('UD_opt_performance', '30'),
    'UD100': ('UD_opt_performance', '100'),
    'HA': ('HA_opt_performance', 0),
    'ECS': ('ECS_opt_performance', 0)
    }

MGD2cmd = 3785.412

#%%
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
    opex.columns = [i[1].split(' cost')[0] for i in opex.columns]
    opex['Lime stabilization'] = opex.loc[:,['Lime', 'Lime stablization energy']].sum(axis=1)
    # uncomment the following line for heat recovery confifuration:
    opex.drop(columns=['Lime stablization energy', 'CHP heat recovery','CHP electricity recovery','Total OPEX [USD/d]','Total OPEX CHP [USD/d]'], inplace=True)
    # uncomment the following line for CHP confifuration:
    # opex.drop(columns=['Lime stablization energy', 'Heat recovery','Total OPEX [USD/d]','Total OPEX CHP [USD/d]'], inplace=True)
    if normalize_by_flow:
        opex /= (10*MGD2cmd)
    return opex

# %%

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
                       color=c, hatch=hat, 
                       # hatch_linewidth=0.5
                       )
        handles.append(patch)
        y_offset += y

    ax.bar(x, y_offset, width=0.65, ec='black', fc=(1,1,1,0), lw=0.25)
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
    fig, axes = plt.subplots(ncols=len(configs), sharey=True, figsize=(18, 5))
    # fig, axes = plt.subplots(ncols=len(configs), sharey=True, figsize=(12, 5))
    handles = []
    
    global_min = 0.0
    global_max = 0.0
    
    for ID, ax in zip(configs, axes):
        df = opex.loc[(slice(None), ID),:]
        df.index = df.index.droplevel(1)
        for i in data_handles.keys():
            if i not in df.index:
                df.loc[i,:] = 0
        x = np.arange(df.shape[0])
        # y_offset = np.zeros(len(x))
        y_offset_pos = np.zeros(len(x), dtype=float)
        y_offset_neg = np.zeros(len(x), dtype=float)
        for cat, style in patch_dct.items():
            c, hat = style
            y = df.loc[:,cat].to_numpy()
            y_pos = np.where(y >= 0, y, 0.0)
            y_neg = np.where(y < 0,  y, 0.0)
            
            patch_pos = ax.bar(
                x, y_pos, bottom = y_offset_pos, width=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values
            )
            
            patch_neg = ax.bar(
                x, y_neg, bottom = y_offset_neg, width=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values
            )
            # patch = ax.bar(x, y, bottom=y_offset, width=0.65, 
            #                color=c, hatch=hat, 
            #                # hatch_linewidth=0.5,
            #                tick_label=df.index.values)
            if ID == 'B1': 
                if np.any(y_pos):
                    handles.append(patch_pos)
                else:
                    handles.append(patch_neg)
                # handles.append(patch)
            # y_offset += y
            y_offset_pos += y_pos
            y_offset_neg += y_neg
        y_total = y_offset_pos + y_offset_neg
        ax.bar(x, y_total, width=0.65, ec='black', fc=(1,1,1,0), lw=0.8)
        ax.axhline(0, color='black', lw=0.8)
        # ax.bar(x, y_offset, width=0.65, ec='black', fc=(1,1,1,0), lw=0.25)
        xx = x[y_total == 0]
        # xx = x[y_offset == 0]
        if xx: ax.scatter(xx, 0.01, s=12, c='black', marker='x', linewidths=0.75)
        # uncomment the following line for heat recovery confifuration:
        ax.set_title(ID, fontsize=14, fontweight='bold')
        # uncomment the following three lines for CHP confifuration:
        # ID_title = ID + 'E' if ID.endswith('1') else ID
        # ax.set_title(ID_title, fontsize=14, fontweight='bold')
        ax.set_xlim((-0.75, df.shape[0]-0.25))
        if ID == 'B1': 
            ax.set_ylabel('$\mathbf{OPEX}$ [USD·m$^{-3}$]', fontname='Arial', fontsize=14)
        ax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=12)
        ax.tick_params(axis='y', which='minor', direction='inout', length=5)
        ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=11, labelrotation=90)
        if ID == 'N2':
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
               bbox_to_anchor=(0.08, 0.9), loc='upper left', 
               framealpha=1, fontsize=12)
    fig.subplots_adjust(wspace=0, bottom=0.16, top=0.94, left=0.06, right=0.985)
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
    fig, axes = plt.subplots(nrows=len(configs), sharex=True, figsize=(6, 18))
    handles = []
    
    global_min = 0.0
    global_max = 0.0

    for ID, ax in zip(configs, axes):
        df = opex.loc[(slice(None), ID),:]
        df.index = df.index.droplevel(1)
        for i in data_handles.keys():
            if i not in df.index:
                df.loc[i,:] = 0
        x = np.arange(df.shape[0])[::-1]
        # y_offset = np.zeros(len(x))
        # two separate offsets for stacking
        y_offset_pos = np.zeros(len(x), dtype=float)
        y_offset_neg = np.zeros(len(x), dtype=float)
        for cat, style in patch_dct.items():
            c, hat = style
            y = df.loc[:,cat].to_numpy()
            # split positive and negative parts
            y_pos = np.where(y >= 0, y, 0.0)
            y_neg = np.where(y < 0,  y, 0.0)
            
            patch_pos = ax.barh(
                x, y_pos, left=y_offset_pos, height=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values
            )
            
            patch_neg = ax.barh(
                x, y_neg, left=y_offset_neg, height=0.65,
                color=c, hatch=hat,
                tick_label=df.index.values
            )
            # patch = ax.barh(x, y, left=y_offset, height=0.65, 
            #                 color=c, hatch=hat, 
            #                 # hatch_linewidth=0.5,
            #                 tick_label=df.index.values)
            if ID == 'B1': 
                # handles.append(patch)
                if np.any(y_pos):
                    handles.append(patch_pos)
                else:
                    handles.append(patch_neg)
            y_offset_pos += y_pos
            y_offset_neg += y_neg
            # y_offset += y
        y_total = y_offset_pos + y_offset_neg
        ax.barh(x, y_total, height=0.65, ec='black', fc=(1, 1, 1, 0), lw=0.8)
        ax.axvline(0, color='black', lw=0.8)
        # ax.barh(x, y_offset, height=0.65, ec='black', fc=(1,1,1,0), lw=0.25)
        xx = x[y_total == 0]
        # xx = x[y_offset == 0]
        if len(xx): ax.scatter(0.01, xx, s=12, c='black', marker='x', linewidths=0.75)
        # uncomment the following line for heat recovery confifuration:
        ax.set_title(ID+'  ', fontsize=14, fontweight='bold',
                     loc='right', y=0.4, pad=0)
        # uncomment the following three lines for CHP confifuration:
        # ID_title = ID + 'E' if ID.endswith('1') else ID
        # ax.set_title(ID_title, fontsize=14, fontweight='bold',
        #              loc='right', y=0.4, pad=0)    
        ax.set_ylim((-0.75, df.shape[0]-0.25))
        if ID == 'N2': 
            ax.set_xlabel('$\mathbf{OPEX}$ [USD·m${^{-3}}$]', fontname='Arial', fontsize=14)
        ax.tick_params(axis='x', which='major', direction='inout', length=8, labelsize=12)
        ax.tick_params(axis='x', which='minor', direction='inout', length=5)
        ax.tick_params(axis='y', which='major', direction='out', length=4, labelsize=11)
        if ID == 'B1':
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
               bbox_to_anchor=(0.5, 0.965), 
               loc='upper left', framealpha=1, fontsize=12)
    fig.subplots_adjust(hspace=0, left=0.16, right=0.97, bottom=0.045, top=0.985)
    if save_as:
        fig.savefig(ospath.join(figures_path, save_as),
                    dpi=300, transparent=True)
    else:
        return fig, axes
#%%
def screen_chp_electricity_recovery(
    opex,
    recovery_col="CHP electricity recovery",
    energy_cols=("Aeration energy", "Pumping energy", "Mixing energy"),
    print_details=False,
):
    """
    Screen rows where |CHP electricity recovery| >= sum(energy_cols).
    If none meet the criterion, show near-misses with ratio in [0.95, 1.0].
    """
    import numpy as np
    import pandas as pd

    # --- checks ---
    missing = [c for c in (recovery_col, *energy_cols) if c not in opex.columns]
    if missing:
        raise KeyError(f"Missing columns in opex: {missing}")

    # --- compute criterion ---
    energy_sum = opex.loc[:, energy_cols].sum(axis=1)
    rec_abs = opex.loc[:, recovery_col].abs()

    # Avoid division by zero
    energy_sum_safe = energy_sum.replace(0, np.nan)
    ratio = rec_abs / energy_sum_safe

    mask_main = rec_abs >= energy_sum
    energy_sum_str = " + ".join(energy_cols)

    # --- choose result set ---
    if mask_main.any():
        header = f"Pairs where |{recovery_col}| >= {energy_sum_str}:"
        mask_use = mask_main
    else:
        mask_near = ratio.between(0.8, 1.0, inclusive="both")
        if not mask_near.any():
            print(f"No pairs meet |{recovery_col}| >= {energy_sum_str}, "
                  f"and no near-misses found with 0.8 ≤ ratio ≤ 1.0.")
            return pd.DataFrame()

        header = (f"No pairs meet |{recovery_col}| >= {energy_sum_str}. "
                  f"Showing near-misses with 0.8 ≤ |{recovery_col}|/({energy_sum_str}) ≤ 1.0:")
        mask_use = mask_near

    # --- results table ---
    result = opex.loc[mask_use, [recovery_col, *energy_cols]].copy()
    result["energy_sum"] = energy_sum.loc[mask_use]
    result["recovery_abs"] = rec_abs.loc[mask_use]
    result["ratio_abs_recovery_to_energy"] = ratio.loc[mask_use]

    # --- print scenario and ID ---
    print(header)
    for scenario, ID in result.index:
        if print_details:
            row = result.loc[(scenario, ID)]
            print(
                f"  ID={ID}, scenario={scenario} | "
                f"|recovery|={row['recovery_abs']:.4g}, "
                f"energy_sum={row['energy_sum']:.4g}, "
                f"ratio={row['ratio_abs_recovery_to_energy']:.3f}"
            )
        else:
            print(f"  ID={ID}, scenario={scenario}")

    return result


#%%
if __name__ == '__main__':
    opex = compile_opex()
    stacked_bar(opex, save_as='opex_w_CHP_heat_recovery.tif')
    # stacked_bar(opex, save_as='opex_wrrfs.png')
    horizontal_stacked_bar(opex, save_as='opexh_w_CHP_heat_recovery.tif')
    # single_scenario_stacked_bar(save_as='opex.tif')
    hits = screen_chp_electricity_recovery(opex,energy_cols=("Aeration energy",), print_details=True)