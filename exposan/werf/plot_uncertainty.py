# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns
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
dg = Guest.dgreen.HEX
dy = Guest.dyellow.HEX

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["figure.autolayout"] = False
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

# %%
# path = ospath.join(results_path, 'HA_F1_UA_opex_stats.xlsx')
# data = load_data(path, sheet=None)
# bl = []
# for strength in ('low', 'mid', 'high'):
#     bl.append(
#         load_data(
#             ospath.join(results_path, f'HA_F1_UA-{strength}.xlsx'), 
#             header=[0,1], skiprows=[2,], sheet='no_HA'
#             )[('OPEX', 'Total OPEX [USD/d]')]
#         )
# bl = pd.DataFrame(bl, index=['low', 'mid', 'high']).T



# %%

colors = {
    'low': (b, db),
    'mid': (g, dg),
    'high': (y, dy)
    }

def plot_linebox(data=None, bl=None):
    if data is None:
        path = ospath.join(results_path, 'HA_F1_UA_opex_stats.xlsx')
        data = load_data(path, sheet=None)
    fig, (ax, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(8, 6), width_ratios=[3,1])
    pos = 1
    for strength, quantiles in data.items():
        x = [float(i)*100 for i in quantiles.index[1:]]
        color, dcolor = colors[strength]
        ax.fill_between(
            x, y1=quantiles.loc['0.50':, 0.05], y2=quantiles.loc['0.50':, 0.95], 
            alpha=0.3, color=color, lw=0, 
            )
        ax.fill_between(
            x, y1=quantiles.loc['0.50':, 0.25], y2=quantiles.loc['0.50':, 0.75], 
            alpha=0.5, color=color, lw=0,
            )
        ax.plot(x, quantiles.loc['0.50':, 0.5], color=dcolor)
        ax.plot(x, quantiles.loc['0.50':, [0.25, 0.75]], color=color, lw=1.5, ls='--')

        lineprops = dict(color=dcolor)
        bplot = ax1.boxplot(
            bl[strength], whis=(5, 95), 
            positions=[pos,], tick_labels=[strength,], widths=0.5,
            showfliers=False, patch_artist=True,
            whiskerprops=lineprops, capprops=lineprops, medianprops=lineprops
            )
        bplot['boxes'][0].set(facecolor=color, edgecolor=dcolor)
        pos += 1
        
    ax.set_xlim((x[0], x[-1]))
    xmin, xmax = ax1.get_xlim()
    ax1.set_xlim((xmin-0.25, xmax))
    ax.set_ylim((0, 1e4))
    
    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    ax.tick_params(axis='both', direction='inout', **major)
    ax.tick_params(axis='both', direction='inout', **minor)
    ax1.tick_params(axis='both', direction='inout' **major)
    ax1.tick_params(axis='y', direction='inout', **minor)
    ax1.tick_params(axis='x', which='minor', length=0)
    
    ax2y = ax1.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax.set_title('Configuration F1 + HA', fontsize=14, fontweight='bold')
    ax1.set_title('F1', fontsize=14, fontweight='bold')
    ax.set_xlabel('$\mathbf{HA\ NH_4\ recovery\ efficiency}$ [%]', fontname='Arial', fontsize=13)
    ax.set_ylabel('$\mathbf{OPEX}$ [USD·day${^{-1}}$]', fontname='Arial', fontsize=13)
    ax1.set_xlabel('Wastewater \n strength', fontname='Arial', fontweight='bold', fontsize=13)
    
    fig.subplots_adjust(wspace=0, bottom=0.25)
    
    fig.savefig(ospath.join(figures_path, 'HA_F1_UA.png'),
                dpi=300, transparent=True)
    
# %%
def compile_dopex(ID='F1', bl=None, save_as=''):
    if bl is None:
        bl = []
        for strength in ('low', 'mid', 'high'):
            bl.append(
                load_data(
                    ospath.join(results_path, f'HA_F1_UA-{strength}.xlsx'), 
                    header=[0,1], skiprows=[2,], sheet='no_HA'
                    )[('OPEX', 'Total OPEX [USD/d]')]
                )
        bl = pd.DataFrame(bl, index=['low', 'mid', 'high']).T
    zz = {}
    xx = None
    for strength in ('low', 'mid', 'high'):
        dfs = load_data(ospath.join(results_path, f'HA_{ID}_UA-{strength}.xlsx'), 
                        header=[0,1], skiprows=[2,], sheet=None)
        dopex = []
        keys = []
        for k, df in dfs.items():
            if k != 'no_HA':
                dopex.append(bl[strength] - df[('OPEX', 'Total OPEX [USD/d]')])
                keys.append(k)
            else:
                if xx is None: 
                    xx = df[('System', 'Electricity price [cents/kWh]')]
                    yy = df[('System', 'Sludge disposal price [USD/wet tonne]')]
        dopex = pd.DataFrame(dopex, index=keys).T
        zz[strength] = dopex
    if save_as:
        with pd.ExcelWriter(save_as) as writer:
            for k, v in zz.items():
                v.to_excel(writer, sheet_name=k)
    return xx, yy, zz

# %%
def compile_all_opex(ID='F1', save_as=''):
    bl_dfs = load_data(ospath.join(results_path, f'{ID}_UA.xlsx'), 
                       header=[0,1], skiprows=[2,], sheet=None)
    
    bl_opex = []
    for k, df in bl_dfs.items():
        bl_opex.append(df[('OPEX', 'Total OPEX [USD/d]')])
    bl_opex = pd.DataFrame(bl_opex, index=[int(k) for k in bl_dfs.keys()]).T
    bl_opex['sample_ID'] = bl_opex.index
    bl_opex = pd.melt(bl_opex, id_vars=['sample_ID'], value_vars=bl_opex.columns[:-1], 
                      var_name='Strength', value_name='OPEX')
    bl_opex['f_rmv'] = np.nan
    bl_opex['HA'] = 0
    
    hasha = []
    for f_str in np.linspace(0,1,11)*100:
        dfs = load_data(ospath.join(results_path, f'HA_{ID}_UA-strength{int(f_str)}.xlsx'), 
                        header=[0,1], skiprows=[2,], sheet=None)
        opex = pd.DataFrame()
        for f_rmv, df in dfs.items():
            opex[f_rmv] = df[('OPEX', 'Total OPEX [USD/d]')]
        opex['sample_ID'] = df.index
        opex = pd.melt(opex, id_vars=['sample_ID'], value_vars=opex.columns[:-1], 
                       var_name='f_rmv', value_name='OPEX')
        opex['Strength'] = int(f_str)
        opex['HA'] = 1
        hasha.append(opex)
    
    allopex = pd.concat([bl_opex, *hasha])
    allopex['Strength'] = allopex.Strength.astype(int)
    if save_as: allopex.to_excel(ospath.join(results_path, save_as))
    else: return allopex
    
# %%
from scipy.interpolate import make_interp_spline as mis

def plot_opex_by_strength(stats=None, data=None, ID='F1'):

    if data is None: data = compile_all_opex()
    if stats is None:
        stats = load_data(
            ospath.join(results_path, f'HA_{ID}_UA_opex_stats_by_strength.xlsx'),
            header=[0,1], skiprows=[2,]
            )
    cod = stats.Influent.COD
    mid = (cod.max()/2-cod.min())/(cod.max()-cod.min()) * 100
    df = stats.dOPEX_quantiles
    x = [int(i) for i in df.index]
    _x = np.linspace(min(x), max(x), 50)
    def smoo(y, k=3):
        spl = mis(x, y, k=k)
        return spl(_x)
    
    fig, (tax, bax) = plt.subplots(2, 1, figsize=(6,11), sharex=True)
    tax = sns.violinplot(
        data=data, x='Strength', y='OPEX', hue='HA', ax=tax,
        palette={0: o, 1: b}, bw_adjust=0.75, 
        linewidth=0.8, linecolor='black', gap=0., 
        split=True, cut=0, inner='quart', 
        inner_kws=dict(linewidth=1),
        saturation=1,
        native_scale=True
        )
    handles, labels = tax.get_legend_handles_labels()
    tax.legend(handles=handles, labels=['Baseline', 'With NH$_4$ recovery by HA'])
    tax.set_ylabel('$\mathbf{OPEX}$ [USD·day${^{-1}}$]', 
                   fontname='Arial', fontsize=13, labelpad=0, linespacing=0.8)

    bax.fill_between(
        _x, y1=smoo(df[0.05]), y2=smoo(df[0.95]), 
        alpha=0.3, color=g, lw=0, 
        )
    bax.fill_between(
        _x, y1=smoo(df[0.25]), y2=smoo(df[0.75]), 
        alpha=0.5, color=g, lw=0,
        )
    bax.plot(_x, smoo(df[0.25]), color=g, lw=1.5, ls='--')
    bax.plot(_x, smoo(df[0.75]), color=g, lw=1.5, ls='--')
    bax.plot(_x, smoo(df[0.5]), color=dg)
    bax.plot(x, df[0.5], marker='o', lw=0, ms=3, mec=dg, mfc=g)
    bax.set_ylabel('$\mathbf{OPEX\ saving}$ [USD·day${^{-1}}$]', 
                    fontname='Arial', fontsize=13, labelpad=8, linespacing=0.8)
    bax.set_xlabel('$\mathbf{Wastewater\ strength}$ [-]', 
                    fontname='Arial', fontsize=13, labelpad=1, linespacing=0.8)
    bax.xaxis.set_ticks([0, mid, 100], ['low', 'medium', 'high'])    

    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    for ax in (tax, bax):
        ax.tick_params(axis='both', direction='inout', **major)
        ax.tick_params(axis='y', direction='inout', **minor)
        ax.tick_params(axis='x', which='minor', length=0)
        
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=4)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    xmin, xmax = bax.get_xlim()
    def get_xlim(xx):
        xmn = xx.min()
        xmx = xx.max()
        return xmn + (xmx-xmn) * xmin/100, xmn + (xmx-xmn) * xmax/100
    
    loc = -0.25    
    for xvar in ('BOD', 'TSS', 'TN', 'NH4-N', 'TP', 'OP'):
        twnx = bax.twiny()
        twnx.spines['bottom'].set_position(('axes', loc))
        twnx.xaxis.tick_bottom()
        twnx.xaxis.set_label_position('bottom')  # Move label to bottom
        twnx.tick_params(axis='x', direction='out', which='major', length=4.5, labelsize=12)
        twnx.tick_params(axis='x', direction='out', which='minor', length=2.5)
        twnx.plot(stats.Influent[xvar], df[0.5], marker='o', lw=0, ms=3, mec=dg, mfc=g)
        twnx.set_xlim(get_xlim(stats.Influent[xvar]))
        if xvar == 'NH4-N': xvar = 'NH_4-N'
        elif xvar == 'OP': xvar = 'Ortho-P'
        twnx.set_xlabel('$\mathbf{%s}$ [mg·L${^{-1}}$]' % xvar, 
                        fontname='Arial', fontsize=13, labelpad=0, linespacing=0.8)
        loc -= 0.25

    fig.subplots_adjust(hspace=0, top=0.95, bottom=0.45, left=0.2)
    fig.savefig(ospath.join(figures_path, f'HA_{ID}_UA_by_strength.png'),
                dpi=300, transparent=True)
    
# %%
    
if __name__ == '__main__':
    # plot_linebox()
    # dopex = compile_dopex(save_as='HA_F1_dopex.xlsx')
    data = compile_all_opex()
    plot_opex_by_strength(data=data)
    