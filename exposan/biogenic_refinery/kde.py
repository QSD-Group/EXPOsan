# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 20:38:24 2023

@author: lanet
"""

# libraries
import seaborn as sns
import pandas as pd
import matplotlib as mpl



df = pd.read_excel('/Users/stetsonrowles/Dropbox/Mac (3)/Documents/GitHub/EXPOsan/exposan/biogenic_refinery/data/kde_combined.xlsx')
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

# set plot theme
sns.set_theme(style="ticks", font='DejaVu Sans')

pal = [(0.376, 0.757, 0.812), 
       (0.475, 0.749, 0.510)]

g = sns.JointGrid(data=df, x=df.mass_bc, y=df.mass_c, hue="old_new", ratio=10, space=0, 
                  xlim=(0, 500), ylim=(0, 80), palette=pal, marginal_ticks=True
                  )
g.plot_joint(sns.kdeplot, fill=True)
g.ax_joint.legend_.set_title(None)

sns.boxplot(data=df, x=g.hue, y=g.y, ax=g.ax_marg_y,  palette=pal, showfliers=False)
sns.boxplot(data=df, y=g.hue, x=g.x, ax=g.ax_marg_x,  palette=pal, showfliers=False)
g.ax_marg_x.set(yticks=[])
g.ax_marg_y.set(xticks=[])

# g = sns.JointGrid(data=df, x=df.mass_bc, y=df.mass_c, ratio=10, space=0, # xlim=(0, 800), ylim=(0, 80), 
#                   marginal_ticks=True)
# g.plot_joint(sns.kdeplot, fill=True, color=(0.376, 0.757, 0.812))
# g.plot_marginals(sns.boxplot, color= (0.376, 0.757, 0.812))

sns.despine(ax=g.ax_marg_x, left=True)
sns.despine(ax=g.ax_marg_y, bottom=True)

# set axes and ticks 
g.ax_joint.tick_params(axis='both', which='both', direction='inout')
# g.ax_joint.tick_params(axis='both', which='minor', direction='inout')
g.ax_marg_x.tick_params(axis='x', which='both', direction='out')
g.ax_marg_y.tick_params(axis='y', which='both', direction='out')
# g.ax_joint.set_yticks(np.arange(0,80,10))
# g.ax_joint.set_yticks(np.arange(5,80,10), minor=True)
# g.ax_joint.set_xticks(np.arange(0,850,100))
# g.ax_joint.set_xticks(np.arange(50,800,100), minor=True)
g.ax_joint.set_xlabel("Mass biochar produced [kg\N{DOT OPERATOR}d\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")
g.ax_joint.set_ylabel("Mass carbon sequestered [kg\N{DOT OPERATOR}d\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")

g.savefig('/Users/lanet/Desktop/GitHub/EXPOsan/exposan/biogenic_refinery/results/kde_combined.svg', format='svg')