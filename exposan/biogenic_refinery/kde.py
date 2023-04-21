# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 20:38:24 2023

@author: lanet
"""

# libraries
import seaborn as sns
import pandas as pd
import numpy as np

df = pd.read_excel('/Users/lanet/Desktop/GitHub/EXPOsan/exposan/biogenic_refinery/data/kde_combined.xlsx')

# set plot theme
sns.set_theme(style='ticks')


g = sns.JointGrid(data=df, x=df.mass_bc, y=df.mass_c, ratio=10, space=0, # xlim=(0, 800), ylim=(0, 80), 
                  marginal_ticks=True)
g.plot_joint(sns.kdeplot, fill=True, color=(0.376, 0.757, 0.812))
g.plot_marginals(sns.boxplot, color=(0.376, 0.757, 0.812))

# UDDT latrine: x=df.pit_biochar, y=df.pit_carbon, color=(0.475, 0.749, 0.510)
# UDDT: x=df.UDDT_biochar, y=df.UDDT_carbon, color=(0.376, 0.757, 0.812)
# septic tank: x=df.septic_biochar, y=df.septic_carbon, color=(0.565, 0.569, 0.557)

# set axes and ticks 
g.ax_joint.tick_params(which='both', direction='inout')
g.ax_marg_x.tick_params(axis='x', which='both', direction='out')
g.ax_marg_y.tick_params(axis='y', which='both', direction='out')
# g.ax_joint.set_yticks(np.arange(0,80,10))
# g.ax_joint.set_yticks(np.arange(5,80,10), minor=True)
# g.ax_joint.set_xticks(np.arange(0,850,100))
# g.ax_joint.set_xticks(np.arange(50,800,100), minor=True)
g.ax_joint.set_xlabel("Mass biochar produced [kg/day]")
g.ax_joint.set_ylabel("Mass carbon sequestered [kg/day]")

g.savefig('/Users/lanet/Desktop/GitHub/EXPOsan/exposan/biogenic_refinery/results/kde_combined.svg', format='svg')