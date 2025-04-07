# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd, numpy as np, matplotlib.pyplot as plt
import matplotlib 
from qsdsan.utils import load_data, ospath
from exposan.werf import results_path, figures_path

plt.rcParams['font.sans-serif'] = 'Arial'

#%%

# # diff = load_data(ospath.join(results_path, 'performance.xlsx'), sheet='diff').T

thresholds = {
    'BOD': [10]*18,
    'NH4 N': [40]*6 + [2]*12,
    'TN': [40]*9 + [10]*9,
    'TP': [7]*9 + [2]*9,
    }
thresholds = pd.DataFrame.from_dict(
    thresholds, orient='index', 
    columns=('B1', 'B2', 'B3', 'C1', 'C2', 'C3', 
            'E2', 'E2P', 'F1', 
            'G1', 'G2', 'G3', 'H1', 'I1', 'I2', 'I3', 'N1', 'N2'))


#%%
def heatmap(data, colorbar=True, cmap='viridis', 
            show_ticklabels=True, row_labels=[], col_labels=[],
            annotate=None, valfmt="{x:.1f}",
            txtcolors=None, annotate_kw={}, 
            save_as=None, **kwargs):
    fig, ax = plt.subplots(figsize=(12,6))
    im = ax.imshow(data, cmap=cmap, **kwargs)
    if colorbar:
        cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05)
        if cmap == 'bwr':
            cbar.ax.set_xticks(ticks=[-1,0,1], labels=['-100%', '0%', '+100%'], fontsize=12)
        else:
            cbar.ax.set_xticks(ticks=[0,1], labels=['min', 'max'], fontsize=12)

    # Show all ticks and label them with the respective list entries.
    if show_ticklabels:
        row_labels = row_labels or data.index
        col_labels = col_labels or data.columns
    else:
        row_labels = col_labels = []
        
    ax.set_xticks(
        range(data.shape[1]), 
        labels=col_labels, fontsize=12,
        # rotation=0, ha="right", rotation_mode="anchor"
        )
    ax.set_yticks(
        range(data.shape[0]), 
        labels=row_labels, fontsize=12,
        )

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)   

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    if annotate is not None:
        kw = dict(horizontalalignment="center",
                  verticalalignment="center",
                  size=11)
        kw.update(annotate_kw)
        
        if isinstance(valfmt, str):
            valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
        
        for i in range(annotate.shape[0]):
            row = annotate.index[i]
            for j in range(annotate.shape[1]):
                col = annotate.columns[j]
                val = annotate.at[row, col]
                if txtcolors is None: color = 'black'
                else: 
                    fill = data.at[row, col]
                    color = txtcolors(row, col, val, fill)
                im.axes.text(j, i, valfmt(val), color=color, **kw)
    
    if save_as:
        fig.savefig(ospath.join(figures_path, save_as), 
                    dpi=300, transparent=True)
    else:
        return fig, ax

#%%
def plot_absolute(data=None, suffix='baseline'):
    if data is None:
        data = load_data(
            ospath.join(results_path, 'UD_performance.xlsx'), sheet='0',
            header=[0,1], skiprows=[2,]
            )
        data.columns = [i[1].split(' [')[0] for i in data.columns]
    
    fill = pd.DataFrame()
    for var, col in data.items():
        fill[var] = (col-col.min())/(col.max()-col.min())
    fill.index = data.index
    fill = fill.T
    
    vals = pd.DataFrame()
    for var, col in data.items():
        if col.min() > 1e5:
            vals[var] = (col*1e-5).round(1)
        else:
            vals[var] = col.round(1)
    vals.index = data.index
    vals = vals.T

    def valfmtwnan(val):
        if str(val) == 'nan': return ''
        return f"{val:.1f}"

    def txtcolors(var, config, val, fill):
        if var in thresholds.index:
            if thresholds.at[var, config] < val:
                return 'red'
        if fill <= 0.6: return 'white'
        return 'black'

    eff_vars = ['COD', 'BOD', 'TN', 'NH4 N', 'TP']
    heatmap(fill.loc[eff_vars,:], annotate=vals.loc[eff_vars,:], valfmt=valfmtwnan, 
            txtcolors=txtcolors, save_as=f'eff_{suffix}.png')
    heatmap(fill.iloc[8:,:], show_ticklabels=False, colorbar=False,
            annotate=vals.iloc[8:,:], valfmt=valfmtwnan, 
            txtcolors=txtcolors, save_as=f'op_{suffix}.png')

plot_absolute()

#%%
# diff = diff.T
diff_color = diff.copy()
diff_color[diff_color > 1] = 1
diff_color[diff_color < -1] = -1

#%%
data = diff_color.loc[['COD', 'BOD', 'TN', 'NH4', 'TP']]
fig, ax = plt.subplots(figsize=(12,6))
im = ax.imshow(data, cmap='bwr', vmin=-1, vmax=1)
cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05)
cbar.ax.set_xticks(ticks=[-1,0,1], labels=['-100%', '0%', '+100%'], fontsize=12)

# Show all ticks and label them with the respective list entries.
ax.set_xticks(range(data.shape[1]), labels=data.columns, fontsize=12,
              # rotation=0, ha="right", rotation_mode="anchor"
              )
ax.set_yticks(range(data.shape[0]), labels=data.index, fontsize=12,)

# Let the horizontal axes labeling appear on top.
ax.tick_params(top=True, bottom=False,
               labeltop=True, labelbottom=False)

# Turn spines off and create white grid.
ax.spines[:].set_visible(False)

ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)


kw = dict(horizontalalignment="center",
          verticalalignment="center",
          size=10)

# Get the formatter in case a string is supplied
# valfmt = matplotlib.ticker.StrMethodFormatter("{x:.1%}")

# Loop over the data and create a `Text` for each "pixel".
# Change the text's color depending on the data.
texts = []
for i in range(data.shape[0]):
    var = data.index[i]
    row = data.loc[var]
    for j in range(data.shape[1]):
        config = data.columns[j]
        txt = diff.at[var, config]
        if abs(txt) >= 0.1: txt = f"{txt:.0%}"
        else: txt = f"{txt:.1%}"
        color = None if str(txt) == 'nan' else 'black'
        text = im.axes.text(j, i, txt, color=color, **kw)
        texts.append(text)

fig.savefig(ospath.join(figures_path, 'eff_diff'), 
            dpi=300, 
            transparent=True
            )

del fig, ax, data
#%%
data = diff_color.loc[['biogas', 'CH4', 'sludge', 'asp_air', 'aed_air', 'total_air']]
fig, ax = plt.subplots(figsize=(12,6))
im = ax.imshow(data, cmap='bwr', vmin=-1, vmax=1)
# cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05)
# cbar.ax.set_xticks(ticks=[-1,0,1], labels=['-100%', '0%', '+100%'], fontsize=12)

# Show all ticks and label them with the respective list entries.
ax.set_xticks(range(data.shape[1]), labels=data.columns, fontsize=12,
              # rotation=0, ha="right", rotation_mode="anchor"
              )
ax.set_yticks(range(data.shape[0]), labels=data.index, fontsize=12,)

# Let the horizontal axes labeling appear on top.
ax.tick_params(top=True, bottom=False,
               labeltop=True, labelbottom=False)

# Turn spines off and create white grid.
ax.spines[:].set_visible(False)

ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)


kw = dict(horizontalalignment="center",
          verticalalignment="center",
          size=10)

# Get the formatter in case a string is supplied
# valfmt = matplotlib.ticker.StrMethodFormatter("{x:.1%}")

# Loop over the data and create a `Text` for each "pixel".
# Change the text's color depending on the data.
texts = []
for i in range(data.shape[0]):
    var = data.index[i]
    row = data.loc[var]
    for j in range(data.shape[1]):
        config = data.columns[j]
        txt = diff.at[var, config]
        if str(txt) == 'nan': txt = ''
        else:
            if abs(txt) >= 0.1: txt = f"{txt:.0%}"
            else: txt = f"{txt:.1%}"
            color = 'black'
        text = im.axes.text(j, i, txt, color=color, **kw)
        texts.append(text)

fig.savefig(ospath.join(figures_path, 'op_diff'), 
            dpi=300, 
            transparent=True
            )