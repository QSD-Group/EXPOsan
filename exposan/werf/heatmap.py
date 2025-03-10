# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 08:55:12 2025

@author: joy_c
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import matplotlib 
from qsdsan.utils import load_data, ospath
from exposan.werf import results_path, figures_path

#%%

# def heatmap(data, row_labels, col_labels, ax=None,
#             cbar_kw=None, cbarlabel="", **kwargs):
#     """
#     Create a heatmap from a numpy array and two lists of labels.

#     Parameters
#     ----------
#     data
#         A 2D numpy array of shape (M, N).
#     row_labels
#         A list or array of length M with the labels for the rows.
#     col_labels
#         A list or array of length N with the labels for the columns.
#     ax
#         A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
#         not provided, use current Axes or create a new one.  Optional.
#     cbar_kw
#         A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
#     cbarlabel
#         The label for the colorbar.  Optional.
#     **kwargs
#         All other arguments are forwarded to `imshow`.
#     """

#     if ax is None:
#         ax = plt.gca()

#     if cbar_kw is None:
#         cbar_kw = {}

#     # Plot the heatmap
#     im = ax.imshow(data, **kwargs)

#     # Create colorbar
#     cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
#     cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

#     # Show all ticks and label them with the respective list entries.
#     ax.set_xticks(range(data.shape[1]), labels=col_labels,
#                   rotation=-30, ha="right", rotation_mode="anchor")
#     ax.set_yticks(range(data.shape[0]), labels=row_labels)

#     # Let the horizontal axes labeling appear on top.
#     ax.tick_params(top=True, bottom=False,
#                    labeltop=True, labelbottom=False)

#     # Turn spines off and create white grid.
#     ax.spines[:].set_visible(False)

#     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
#     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
#     ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
#     ax.tick_params(which="minor", bottom=False, left=False)

#     return im, cbar


# def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
#                      textcolors=("black", "white"),
#                      threshold=None, **textkw):
#     """
#     A function to annotate a heatmap.

#     Parameters
#     ----------
#     im
#         The AxesImage to be labeled.
#     data
#         Data used to annotate.  If None, the image's data is used.  Optional.
#     valfmt
#         The format of the annotations inside the heatmap.  This should either
#         use the string format method, e.g. "$ {x:.2f}", or be a
#         `matplotlib.ticker.Formatter`.  Optional.
#     textcolors
#         A pair of colors.  The first is used for values below a threshold,
#         the second for those above.  Optional.
#     threshold
#         Value in data units according to which the colors from textcolors are
#         applied.  If None (the default) uses the middle of the colormap as
#         separation.  Optional.
#     **kwargs
#         All other arguments are forwarded to each call to `text` used to create
#         the text labels.
#     """

#     if not isinstance(data, (list, np.ndarray)):
#         data = im.get_array()

#     # Normalize the threshold to the images color range.
#     if threshold is not None:
#         threshold = im.norm(threshold)
#     else:
#         threshold = im.norm(data.max())/2.

#     # Set default alignment to center, but allow it to be
#     # overwritten by textkw.
#     kw = dict(horizontalalignment="center",
#               verticalalignment="center")
#     kw.update(textkw)

#     # Get the formatter in case a string is supplied
#     if isinstance(valfmt, str):
#         valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

#     # Loop over the data and create a `Text` for each "pixel".
#     # Change the text's color depending on the data.
#     texts = []
#     for i in range(data.shape[0]):
#         for j in range(data.shape[1]):
#             kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
#             text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
#             texts.append(text)

#     return texts

#%%
baseline = load_data(ospath.join(results_path, 'performance.xlsx'), sheet='baseline')
diff = load_data(ospath.join(results_path, 'performance.xlsx'), sheet='diff').T

bl_color = pd.DataFrame()
for var, col in baseline.items():
    bl_color[var] = (col-col.min())/(col.max()-col.min())

bl_color.index = baseline.index
bl_color = bl_color.T

bl_txt = pd.DataFrame()
for var, col in baseline.items():
    if col.min() > 1e5:
        bl_txt[var] = (col*1e-5).round(1)
    else:
        bl_txt[var] = col.round(1)

bl_txt.index = baseline.index
bl_txt = bl_txt.T

thresholds = {
    'BOD': [10]*18,
    'NH4': [40]*6 + [2]*12,
    'TN': [40]*9 + [10]*9,
    'TP': [7]*9 + [2]*9,
    }
thresholds = pd.DataFrame.from_dict(thresholds, orient='index', columns=baseline.index)

bl_highlight = bl_txt.copy()
bl_highlight.loc[:, :] = 0
for var, col in baseline.items():
    if var in thresholds.index:
        bl_highlight.loc[var, col > thresholds.loc[var]] = 1

plt.rcParams['font.sans-serif'] = 'Arial'

#%%

data = bl_color.loc[['COD', 'BOD', 'TN', 'NH4', 'TP']]
fig, ax = plt.subplots(figsize=(12,6))
im = ax.imshow(data)
cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05)
cbar.ax.set_xticks(ticks=[0,1], labels=['min', 'max'], fontsize=12)

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
          size=11)

# Get the formatter in case a string is supplied
valfmt = matplotlib.ticker.StrMethodFormatter("{x:.1f}")

# Loop over the data and create a `Text` for each "pixel".
# Change the text's color depending on the data.
texts = []
for i in range(data.shape[0]):
    var = data.index[i]
    row = bl_txt.loc[var]
    for j in range(data.shape[1]):
        config = data.columns[j]
        txt = bl_txt.at[var, config]
        color = 'red' if bl_highlight.at[var, config] else 'grey'
        text = im.axes.text(j, i, valfmt(txt, None), color=color, **kw)
        texts.append(text)

fig.savefig(ospath.join(figures_path, 'eff_baseline'), 
            dpi=300, 
            transparent=True
            )

del fig, ax, data

#%%
data = bl_color.loc[['biogas', 'CH4', 'sludge', 'asp_air', 'aed_air', 'total_air']]
fig, ax = plt.subplots(figsize=(12,6))
im = ax.imshow(data)
# cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05)
# cbar.ax.set_xticks(ticks=[0,1], labels=['min', 'max'], fontsize=12)

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
          size=11)

# Get the formatter in case a string is supplied
valfmt = matplotlib.ticker.StrMethodFormatter("{x:.1f}")

# Loop over the data and create a `Text` for each "pixel".
# Change the text's color depending on the data.
texts = []
for i in range(data.shape[0]):
    var = data.index[i]
    row = bl_txt.loc[var]
    for j in range(data.shape[1]):
        config = data.columns[j]
        txt = bl_txt.at[var, config]
        if str(txt) == 'nan': txt = ''
        else: 
            txt = f"{txt:.1f}"
            color = 'grey'
        text = im.axes.text(j, i, txt, color=color, **kw)
        texts.append(text)

fig.savefig(ospath.join(figures_path, 'op_baseline'), 
            dpi=300, 
            transparent=True
            )
del fig, ax, data

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