# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import qsdsan as qs, numpy as np
from exposan.hap import (
    results_path,
    figures_path
    )
from qsdsan.utils import load_data, ospath
import matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = False
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

#%%
boxprops = dict(alpha=1, facecolor='#60c1cf', edgecolor='black')
meanprops = dict(marker='x', markersize=7, markeredgecolor='black')
flierprops = dict(marker='.', markersize=2, markerfacecolor='#4e4e4e')
medianprops = dict(color='black')

# _market_prices={            # USD/kg
#     'low range': [136,      # puriss https://www.sigmaaldrich.com/US/en/product/sial/04238
#                   194],     # purum  https://www.sigmaaldrich.com/US/en/product/sial/21223
#     'high range': [4000,    # https://www.sigmaaldrich.com/US/en/product/aldrich/900203
#                    9600]    # https://www.sigmaaldrich.com/US/en/product/aldrich/677418
#     }

def plot_mpsp(data=None, seed=None, save=True, 
              market_ranges=[[136, 194], [4000, 9600]],
              single_values=[50, 170]):
    if data is None:
        data = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                          header=[0,1], skiprows=[2,])
    mpsp = data.loc[:,  ('TEA', 'MPSP [USD/kg]')]
    # mpsp = np.random.rand(100) * 100
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(2,6), height_ratios=[1,2])
    fig.subplots_adjust(left=0.3, hspace=0.05)
    for ll, ul in market_ranges:
        ax1.axhspan(ll, ul, alpha=0.3, fill=True, facecolor='#4e4e4e')
    for v in single_values:
        ax1.axhline(v, color='black', linewidth=0.5)
        ax2.axhline(v, color='black', linewidth=0.5)
    ax1.set_ylim(91, 11000)
    ax1.set_yscale('log')
    ax2.boxplot(mpsp,
                whis=(5,95), 
                widths=0.6,
                patch_artist=True,
                showmeans=True,
                boxprops=boxprops,
                meanprops=meanprops,
                flierprops=flierprops,
                medianprops=medianprops
                )
    ax2.set_ylim(0, 55)
    for ax in (ax1, ax2):
        ax.tick_params(axis='y', which='major', direction='inout', length=10, labelsize=12)
        ax.tick_params(axis='y', which='minor', direction='inout', length=6)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=5)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=3)
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        ax.tick_params(bottom=False, which='both', labelbottom=False)
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.tick_params(labeltop=False)
    
    d = 0.35  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
    if save:
        fig.savefig(ospath.join(figures_path, 'msps.png'), 
                    dpi=300, transparent=True)
    else:
        return fig, ax

#%%

def plot_npv_breakdown(data=None, seed=None, save=True):
    pass

#%%
def plot_sensitivity(seed=None, save=True):
    pass

# plot_mpsp(seed=459, save=True)
