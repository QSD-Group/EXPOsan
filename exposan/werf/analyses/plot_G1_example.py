# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 16:31:44 2025

@author: joy_c
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
# %%
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

configs=('B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'E2', 'E2P', 'F1', 
        'G1', 'G2', 'G3', 'H1', 'I1', 'I2', 'I3', 'N1', 'N2')

# %%
# from scipy.interpolate import make_interp_spline as mis

def plot_eff_vs_do(data=None):
    if data is None: 
        data = load_data(
            ospath.join(results_path, 'G1_example_operation.xlsx'), sheet='1d',
            header=[0,1], index_col=None,
            )
        data.columns = [b.split(' [')[0] for a,b in data.columns]
    
    limits = {
        'BOD': (10, '$\mathbf{BOD}$ [mg/L]'),
        'TN': (10, '$\mathbf{TN}$ [mg/L]'),
        'NH4 N': (2, '$\mathbf{NH_4-N}$ [mg/L]'),
        'TP': (2, '$\mathbf{TP}$ [mg/L]'),
        }
    
    x = data.DO
    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    i = 0
    fig, axes = plt.subplots(nrows=4, sharex=True, figsize=(6, 8))
    for col, (y, lab) in limits.items():
        ax = axes[i]
        ax.axhline(y, linestyle='--', color='red')
        ax.plot(x, data[col], linestyle='-', color='blue', marker='o', markersize=3)
        if col == 'NH4 N': ax.set_ylabel(lab, fontname='Arial', fontsize=13, labelpad=10)
        elif col == 'TP': ax.set_ylabel(lab, fontname='Arial', fontsize=13, labelpad=2)
        else: ax.set_ylabel(lab, fontname='Arial', fontsize=13, labelpad=4)
        ax.set_ylim(ax.get_ylim()[0]*0.9, y*1.25)
        ax.xaxis.set_inverted(True)
        
        ax.tick_params(axis='both', direction='inout', **major)
        ax.tick_params(axis='both', direction='inout', **minor)
        
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=4)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
        ax2y.set_yticks(ax.get_yticks())
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        
        if i == 0: ax.legend(labels=['Discharge limit', 'Simulated effluent'], fontsize=12)
        i += 1
    
    ax.set_xlabel('$\mathbf{Aerobic\ zone\ DO\ setpoint}$ [mg/L]',
                  fontname='Arial', fontsize=13)
    fig.subplots_adjust(hspace=0.0)
    fig.savefig(ospath.join(figures_path, 'G1_example_eff1d.png'),
                dpi=300, transparent=True)

def plot_op_vs_do(data=None):
    if data is None: 
        data = load_data(
            ospath.join(results_path, 'G1_example_operation.xlsx'), sheet='1d',
            header=[0,1], index_col=None,
            )
        data.columns = [b.split(' [')[0] for a,b in data.columns]
    
    data['Liquid aeration flowrate'] /= 1e5
    data['Liquid aeration energy'] *= 24/1e3
    data['CH4 production'] *=24/1e3
    cols = {
        'Liquid aeration flowrate': ((1.8, 3.5), '$\mathbf{Aeration\ flowrate}$\n[10$^5$ m$^3$/d]'),
        'Liquid aeration energy': ((3.8, 6.5),'$\mathbf{Aeration\ energy}$\n[10$^3$ kWh/d]'),
        'Sludge production': ((26.2, 28.8), '$\mathbf{Sludge\ production}$\n[tonne/d]'),
        'CH4 production': ((1.17, 1.23), '$\mathbf{Biogas\ production}$\n[10$^3$ kg-CH$_4$/d]'),
        }
    
    x = data.DO
    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    i = 0
    pads = [7, 17, 14, 0]
    lsp = [0.9, 0.9, 1.3, 0.9]
    fig, axes = plt.subplots(nrows=4, sharex=True, figsize=(6, 8))
    for col, (ylim, lab) in cols.items():
        ax = axes[i]
        ax.plot(x, data[col], linestyle='-', color=p, marker='^', markersize=4)
        ax.set_ylim(ylim)
        ax.set_ylabel(lab, fontname='Arial', fontsize=12, labelpad=pads[i], linespacing=lsp[i])
        ax.xaxis.set_inverted(True)
        ax.ticklabel_format(useMathText=True)
        
        ax.tick_params(axis='both', direction='inout', **major)
        ax.tick_params(axis='both', direction='inout', **minor)
        
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=4)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
        ax2y.set_yticks(ax.get_yticks())
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        
        i += 1
    
    ax.set_xlabel('$\mathbf{Aerobic\ zone\ DO\ setpoint}$ [mg/L]',
                  fontname='Arial', fontsize=13)
    fig.subplots_adjust(hspace=0.0, left=0.15)
    fig.savefig(ospath.join(figures_path, 'G1_example_op1d.png'),
                dpi=300, transparent=True)
    
# %%

def togrid(data=None, zcol='Total OPEX', xcol='DO', ycol='Carbon'):
    if data is None: 
        data = load_data(
            ospath.join(results_path, 'G1_example_operation.xlsx'), sheet='2d',
            header=[0,1], index_col=None,
            )
        data.columns = [b.split(' [')[0] for a,b in data.columns]
    
    df = data.pivot(index=xcol, columns=ycol, values=zcol)
    xx, yy = np.meshgrid(df.columns, df.index)
    zz = df.to_numpy()
    df = data.pivot(index=xcol, columns=ycol, values='TN')
    tn = df.to_numpy()
    df = data.pivot(index=xcol, columns=ycol, values='TP')
    tp = df.to_numpy()
    return xx, yy, zz, tn, tp
    
from matplotlib.lines import Line2D

def plot_heatmap(contour=True):
    xx, yy, z, tn, tp = togrid()
    fig, ax = plt.subplots(figsize=(7, 5))
    pos = ax.pcolormesh(xx, yy, z, shading='gouraud',)
    ax.set_ylabel('$\mathbf{Aerobic\ zone\ DO\ setpoint}$ [mg/L]',
                  fontname='Arial', fontsize=13)
    ax.set_xlabel('$\mathbf{External\ carbon\ dosage}$ [kg-COD/hr]',
                  fontname='Arial', fontsize=13)

    cbar = fig.colorbar(pos, ax=ax)
    cbar.ax.tick_params(labelsize=11)
    cbar.ax.set_ylabel('$\mathbf{OPEX}$ [USD/d]', rotation=-90, va="bottom",
                       fontname='Arial', fontsize=13)
    ax.tick_params(axis='both', which='major', direction='inout', length=8, labelsize=12)
    ax.tick_params(axis='both', which='minor', direction='inout', length=5)
    ax2x = ax.secondary_xaxis('top', zorder=3)
    ax2x.tick_params(direction='in', which='major', length=4)
    ax2x.tick_params(direction='in', which='minor', length=2.5)
    ax2x.set_xticks(ax.get_xticks(), labels=[],)
    ax2y = ax.secondary_yaxis('right', zorder=3)
    ax2y.tick_params(direction='in', which='major', length=4)
    ax2y.tick_params(direction='in', which='minor', length=2.5)
    ax2y.set_yticks(ax.get_yticks(), labels=[],)

    if contour:    
        ctn = ax.contour(xx, yy, tn, 
                        colors='white', origin='lower', 
                        linestyles='dashed', linewidths=1, 
                        extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]),
                        levels=[8.0, 10.0, 12.0],
                        )
        ax.clabel(ctn, ctn.levels, inline=True, fontsize=11)
        
        ctp = ax.contour(xx, yy, tp, 
                        colors='black', origin='lower', 
                        linestyles='dotted', linewidths=1, 
                        extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]),
                        levels=[1.0, 1.5, 2.0, 2.5],
                        )
        ax.clabel(ctp, ctp.levels, inline=True, fontsize=11)
        proxies = [Line2D([],[], color=c, linestyle=s, linewidth=1) \
                   for c, s in [('white', '--'), ('black', ':')]]
        ax.legend(handles=proxies, 
                  labels=['Effluent TN [mg/L]', 'Effluent TP [mg/L]'],
                  labelcolor=['white', 'black'], 
                  loc='upper center', framealpha=0)
        save_as = 'G1_example_opex.png'
    else:
        save_as = 'G1_example_opex_noline.png'
    fig.savefig(ospath.join(figures_path, save_as), 
                dpi=300, transparent=True)

# plot_heatmap()
