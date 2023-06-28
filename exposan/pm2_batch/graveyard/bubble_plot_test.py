# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 16:04:14 2023

@author: joy_c
"""

from qsdsan.utils import DictAttrSetter, ospath, time_printer, load_data

from exposan.pm2_batch import (
    create_system,
    data_path,
    results_path,
    )

import matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns
mpl.rcParams['xtick.minor.visible'] = False
mpl.rcParams['ytick.minor.visible'] = False

def MCF_encap_to_susp(df, save_as=''):
    df['Parameter'] = df.index
    df = df.melt(id_vars='Parameter', var_name='Metric', value_name='D')
    fig, ax = plt.subplots(figsize=(6,9))

    sns.scatterplot(df, x='Metric', y='Parameter', 
                    hue=df.p<0.05, 
                    # palette=['black', pal[k]], 
                    size='D', sizes=(0, 350), size_norm=(0,1),
                    legend=False,
                    ax=ax,
                    )
    ax.set_xlim(-0.5, 4.5)
    ax.set_axisbelow(True)
    ax.grid(True, which='major', color='k', linestyle='--', 
            linewidth=0.7, alpha=0.3)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='both', which='major', direction='inout', length=8)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    fig.savefig(save_as, dpi=300, transparent=True)
    return df

file = ospath.join(data_path, 'bubble_plot_data.xlsx')

bubble_data = load_data(file, sheet=0, index_col=None)

MCF_encap_to_susp(bubble_data, 'test')