# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>
    Zixuan Wang <wyatt4428@gmail.com

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

def _rgb_to_hex(rgb_tuple):
    return '#%02x%02x%02x' % rgb_tuple

# Dark variants as specified
dblue   = _rgb_to_hex((53, 118, 127))
dgreen  = _rgb_to_hex((77, 126, 83))
dred    = _rgb_to_hex((156, 75, 80))
dorange = _rgb_to_hex((167, 95, 62))
dyellow = _rgb_to_hex((171, 137, 55))
dgray   = _rgb_to_hex((78, 78, 78))
dpurple = _rgb_to_hex((76, 56, 90))

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

def compile_n_recovery(ID='F1'):
    out = []
    for f_str in np.linspace(0,1,11)*100:
        dfs = load_data(ospath.join(results_path, f'HA_{ID}_UA-strength{int(f_str)}.xlsx'), 
                        header=[0,1], skiprows=[2,], sheet=None)
        rcv = []
        fs_rmv = []
        for f_rmv, df in dfs.items():
            fs_rmv.append(f_rmv)
            rcv.append(df[('System','NH4 recovery [kg-N/d]')][0])
        out.append(rcv)
    
    out = pd.DataFrame(out, index=np.linspace(0,1,11)*100, columns=fs_rmv)
    return out

# %%

def plot_nrcv_by_strength(ID='F1'):
    data = compile_n_recovery(ID)
    stats = data.quantile([0.05, 0.25, 0.5, 0.75, 0.95], axis=1).T
    inf_stats = load_data(
        ospath.join(results_path, f'HA_{ID}_UA_opex_stats_by_strength.xlsx'),
        header=[0,1], skiprows=[2,]
        )
    x = [int(i) for i in stats.index]
    _x = np.linspace(min(x), max(x), 50)
    def smoo(y, k=3):
        spl = mis(x, y, k=k)
        return spl(_x)
    
    fig, ax = plt.subplots(figsize=(6,8))

    ax.fill_between(
        _x, y1=smoo(stats[0.05]), y2=smoo(stats[0.95]), 
        alpha=0.3, color=a, lw=0, 
        )
    ax.fill_between(
        _x, y1=smoo(stats[0.25]), y2=smoo(stats[0.75]), 
        alpha=0.5, color=a, lw=0,
        )
    ax.plot(_x, smoo(stats[0.25]), color=a, lw=1.5, ls='--')
    ax.plot(_x, smoo(stats[0.75]), color=a, lw=1.5, ls='--')
    ax.plot(_x, smoo(stats[0.5]), color=da)
    ax.plot(x, stats[0.5], marker='o', lw=0, ms=3, mec=da, mfc=a)
    ax.set_ylabel('$\mathbf{Ammonia\ recovery}$ [kg-N·day$^{-1}$]', 
                    fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
    ax.set_xlabel('$\mathbf{Wastewater\ strength}$ [-]', 
                    fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
    ax.xaxis.set_ticks([0, 24.9, 100], ['low', 'medium', 'high'])

    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    ax.tick_params(axis='both', direction='inout', **major)
    ax.tick_params(axis='y', direction='inout', **minor)
    ax.tick_params(axis='x', which='minor', length=0)
    
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    xmin, xmax = ax.get_xlim()
    def get_xlim(xx):
        xmn = xx.min()
        xmx = xx.max()
        return xmn + (xmx-xmn) * xmin/100, xmn + (xmx-xmn) * xmax/100
    
    loc = -0.23    
    for xvar in ('BOD', 'TSS', 'TN', 'NH4-N', 'TP', 'OP'):
        twnx = ax.twiny()
        twnx.spines['bottom'].set_position(('axes', loc))
        twnx.xaxis.tick_bottom()
        twnx.xaxis.set_label_position('bottom')  # Move label to bottom
        twnx.tick_params(axis='x', direction='out', which='major', length=4.5, labelsize=12)
        twnx.tick_params(axis='x', direction='out', which='minor', length=2.5)
        twnx.plot(inf_stats.Influent[xvar], stats[0.5], marker='o', lw=0, ms=3, mec=da, mfc=a)
        twnx.set_xlim(get_xlim(inf_stats.Influent[xvar]))
        if xvar == 'NH4-N': xvar = 'NH_4^+-N'
        elif xvar == 'OP': xvar = 'Ortho-P'
        twnx.set_xlabel('$\mathbf{%s}$ [mg·L$^{-1}$]' % xvar, 
                        fontname='Arial', fontsize=13, labelpad=0, linespacing=0.8)
        loc -= 0.23

    fig.subplots_adjust(bottom=0.6, top=0.97, left=0.2, right=0.975)
    fig.savefig(ospath.join(figures_path, f'HA_{ID}_UA_nrecovery_industrail_e.tif'),
                dpi=300, transparent=True)

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
    
    fig, (tax, bax) = plt.subplots(2, 1, figsize=(6,11), sharex=True, sharey=False)
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
    tax.legend(handles=handles, labels=['Without NH$_4^+$ recovery', 'With NH$_4^+$ recovery by HA'])
    tax.set_ylim((0, tax.get_ylim()[1]))
    tax.set_ylabel('$\mathbf{OPEX}$ [USD·day$^{-1}$]', 
                   fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
    

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
    bax.set_ylim((0, bax.get_ylim()[1]))
    bax.set_ylabel('$\mathbf{OPEX\ saving}$ [USD·day$^{-1}$]', 
                    fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
    bax.set_xlabel('$\mathbf{Wastewater\ strength}$ [-]', 
                    fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
    bax.xaxis.set_ticks([0, mid, 100], ['low', 'medium', 'high'])
    
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
        if xvar == 'NH4-N': xvar = 'NH_4^+-N'
        elif xvar == 'OP': xvar = 'Ortho-P'
        twnx.set_xlabel('$\mathbf{%s}$ [mg·L$^{-1}$]' % xvar, 
                        fontname='Arial', fontsize=13, labelpad=0, linespacing=0.8)
        loc -= 0.25

    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    for ax in (tax, bax):
        ax.tick_params(axis='both', direction='inout', **major)
        ax.tick_params(axis='y', direction='inout', **minor)
        ax.tick_params(axis='x', which='minor', length=0)

        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=4)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
        ax2y.set_yticklabels([])
        
    fig.subplots_adjust(hspace=0, top=0.975, bottom=0.45, left=0.2, right=0.975)
    fig.savefig(ospath.join(figures_path, f'HA_{ID}_UA_by_strength_industrail_e.tif'),
                dpi=300, transparent=True)

#%%
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.mathtext import _mathtext as mathtext
import math

color_map = LinearSegmentedColormap.from_list(
    'color_map',
    [r,o,y,g,b,db]
)

def plot_opex_by_strength_and_removal_efficiency(
    ID='F1',
    data=None,
):
    mathtext.FontConstantsBase.sup1 = 0.35

    # this is shared for all 3 sheets
    if data is None:
        from __main__ import compile_all_opex
        data = compile_all_opex()

    # loop over the three quantile sheets
    for sheet_name in ('0.05', '0.50', '0.95'):
        stats = load_data(
            ospath.join(results_path, f'HA_{ID}_UA_opex_stats_by_strength&removal_efficiency.xlsx'),
            header=[0, 1],
            skiprows=[2,],
            sheet=sheet_name
        )

        # ----- data for this sheet -----
        x_strengths = [int(i) for i in stats.index]
        reff_keys = list(stats['dOPEX'].columns)
        reff_sorted = sorted(reff_keys, key=lambda s: float(s))
        y_perc = np.array([float(k) * 100 for k in reff_sorted], dtype=float)

        Z = np.vstack([stats[('dOPEX', rk)].to_numpy() for rk in reff_sorted])

        cod = stats['Influent']['COD']
        mid = (cod.max()/2 - cod.min())/(cod.max() - cod.min()) * 100

        # ----- figure -----
        fig, (tax, bax) = plt.subplots(2, 1, figsize=(6, 12), sharex=True, sharey=False)
        tax.set_zorder(3)
        bax.set_zorder(2)

        # ===== top panel (violin, no legend) =====
        vals = np.sort(data['Strength'].unique())
        inner_min, inner_max = vals[1], vals[-2]
        dplot = data[(data['Strength'] >= inner_min) & (data['Strength'] <= inner_max)].copy()

        tax = sns.violinplot(
            data=dplot, x='Strength', y='OPEX', hue='HA', ax=tax,
            palette={0: a, 1: p}, bw_adjust=0.75,
            linewidth=0.8, linecolor='black', gap=0.,
            split=True, cut=0, inner='quart',
            inner_kws=dict(linewidth=1),
            saturation=1,
            native_scale=True
        )
        # remove the legend
        if tax.legend_ is not None:
            tax.legend_.remove()

        tax.margins(x=0)
        tax.set_xlim(0, 100)
        tax.set_ylim((0, max(0, tax.get_ylim()[1])))
        tax.set_ylabel('$\\mathbf{OPEX}$ [USD·day$^{-1}$]',
                       fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)

        # ===== bottom panel (heatmap) =====
        Xg, Yg = np.meshgrid(x_strengths, y_perc)
        cf_levels = 256
        cf = bax.contourf(Xg, Yg, Z, levels=cf_levels, cmap=color_map, antialiased=True)
        cf.set(edgecolor="face", linewidth=0) 
                # ---- dynamic contour/tick levels rounded to 100 ----
        zmin = float(np.nanmin(Z))
        zmax = float(np.nanmax(Z))
        
        # round to 100s
        zmin_150 = math.ceil(zmin / 150.0) * 150
        zmax_150 = math.ceil(zmax / 150.0) * 150
        span = zmax - zmin
        if 500 < span <= 800:
            # make 5 levels, 100 apart
            contour_levels = [zmin_150 + 150 * i for i in range(5)]
        elif 300 < span <= 500:
            # fall back to 75 spacing
            zmin_75 = math.ceil(zmin / 75.0) * 75
            contour_levels = [zmin_75 + 75 * i for i in range(5)]
        else:
            zmin_50 = math.ceil(zmin / 50.0) * 50
            contour_levels = [zmin_50 + 50 * i for i in range(6)]
        cs = bax.contour(Xg, Yg, Z, levels=contour_levels, colors='k', linewidths=0.7)
        bax.clabel(cs, inline=True, fontsize=10, fmt='%.0f')
        bax.plot(
        40, 70,
        marker='D',
        ms=10,
        mfc='white',
        mec='black',
        mew=1.2,
        zorder=10      # above the contourf
        )

        # get the position of bax in figure coords 
        bbox = bax.get_position() 
        # make a new, narrow axes just to the right of bax # [left, bottom, width, height] 
        cax = fig.add_axes([ bbox.x1 + 0.1, # a little to the right of bax's right edge 
                            bbox.y0+0.4, 
                            0.02, # width of the colorbar 
                            bbox.height-0.14 ])
        cb = plt.colorbar(cf, cax=cax, orientation='vertical')
        cb.set_label('$\\mathbf{OPEX\\ saving}$ [USD·day$^{-1}$]', fontsize=13)
        cb.set_ticks(contour_levels)

        bax.set_ylabel('$\\mathbf{NH_4^+\\ recovery\\ efficiency\\ }$[%]',
                       fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
        bax.set_xlabel('$\\mathbf{Wastewater\\ strength}$ [-]',
                       fontname='Arial', fontsize=13, labelpad=2, linespacing=0.8)
        bax.xaxis.set_ticks([0, mid, 100], ['low', 'medium', 'high'])

        # ----- stacked axes under bax -----
        xmin, xmax = bax.get_xlim()
        def get_xlim(xx):
            xmn = xx.min()
            xmx = xx.max()
            return xmn + (xmx - xmn) * xmin / 100, xmn + (xmx - xmn) * xmax / 100

        loc = -0.25
        for xvar in ('BOD', 'TSS', 'TN', 'NH4-N', 'TP', 'OP'):
            twnx = bax.twiny()
            twnx.spines['bottom'].set_position(('axes', loc))
            twnx.xaxis.tick_bottom()
            twnx.xaxis.set_label_position('bottom')
            twnx.tick_params(axis='x', direction='out', which='major', length=4.5, labelsize=12)
            twnx.tick_params(axis='x', direction='out', which='minor', length=2.5)
            twnx.set_xlim(get_xlim(stats['Influent'][xvar]))
            _xvar_label = xvar
            if _xvar_label == 'NH4-N': _xvar_label = 'NH_4^+-N'
            elif _xvar_label == 'OP':  _xvar_label = 'Ortho-P'
            twnx.set_xlabel(r'$\mathbf{%s}$ [mg·L$^{-1}$]' % _xvar_label,
                            fontname='Arial', fontsize=13, labelpad=0, linespacing=0.8)
            loc -= 0.25

        # ticks & right-side secondary y
        major = dict(which='major', length=8, labelsize=12)
        minor = dict(which='minor', length=5)
        for ax in (tax, bax):
            ax.tick_params(axis='both', direction='inout', **major)
            ax.tick_params(axis='y', direction='inout', **minor)
            ax.tick_params(axis='x', which='minor', length=0)
            ax2y = ax.secondary_yaxis('right')
            ax2y.tick_params(axis='y', which='major', direction='in', length=4)
            ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
            ax2y.set_yticklabels([])

        # y-ticks, hide last
        labels = [f'{int(v)}%' if float(v).is_integer() else f'{v:.0f}%' for v in y_perc]
        labels[-1] = ''
        bax.set_yticks(y_perc)
        bax.set_yticklabels(labels)
        
        #         # remove x ticks on tax to prepare for half plot
        # tax.tick_params(axis='x', bottom=False, labelbottom=False, length=0)
        
        # # remove the FIRST y tick label on tax
        # yticks = tax.get_yticks()
        # ylabels = [lab.get_text() for lab in tax.get_yticklabels()]
        
        # if ylabels:  # make sure there is at least one
        #     ylabels[0] = ''  # blank out the first label
        
        # tax.set_yticks(yticks)
        # tax.set_yticklabels(ylabels)

        fig.subplots_adjust(hspace=0, top=0.975, bottom=0.5, left=0.2, right=0.975)

        out_name = f'HA_{ID}_UA_by_strength&removal_efficiency_{sheet_name}_industrail_e.tif'
        fig.savefig(
            ospath.join(figures_path, out_name),
            dpi=300,
            transparent=True,
            bbox_inches='tight',
            pad_inches=0.15
        )
        plt.close(fig)

# %%
    
if __name__ == '__main__':
    # plot_linebox()
    # dopex = compile_dopex(save_as='HA_F1_dopex.xlsx')
    data = compile_all_opex()
    plot_opex_by_strength(data=data)
    plot_nrcv_by_strength()
    plot_opex_by_strength_and_removal_efficiency()

    