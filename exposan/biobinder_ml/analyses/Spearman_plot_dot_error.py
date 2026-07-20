# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 14:26:23 2025

@author: aliah
"""

import pandas as pd
import matplotlib.pyplot as plt

# Data based on the final selections (parameters passing |ρ| > 0.1 and p < 0.05)
data = {
    'Parameter': ['IRR', 'Electricity_Price', 'Electrode_cost', 'EO_voltage'],
    'CHCU_No_EC_MSP': [0.11, None, None, None],
    'CHCU_EC_MSP': [0.23, 0.10, 0.34, None],
    'DHCU_EC_MSP': [0.41, 0.21, 0.14, None],
    'CHCU_EC_GWP': [None, None, None, 0.47],
    'DHCU_EC_GWP': [None, None, None, 0.28],
}

df = pd.DataFrame(data)
df.set_index('Parameter', inplace=True)

# Plotting
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), sharey=True)

configs_msp = ['CHCU_No_EC_MSP', 'CHCU_EC_MSP', 'DHCU_EC_MSP']
configs_gwp = ['CHCU_EC_GWP', 'DHCU_EC_GWP']
colors_msp = ['#1f77b4', '#ff7f0e', '#2ca02c']
colors_gwp = ['#d62728', '#9467bd']

# Plot for MSP
for config, color in zip(configs_msp, colors_msp):
    axes[0].errorbar(df[config], df.index, xerr=0.02, fmt='o', label=config.replace('_MSP',''), color=color)
axes[0].axvline(0, color='black')
axes[0].axvline(0.1, color='gray', linestyle='--')
axes[0].axvline(-0.1, color='gray', linestyle='--')
axes[0].set_title('MSP(Error bars: ±0.02)')
axes[0].set_xlabel('Spearman ρ')
axes[0].legend(title='Configuration')

# # Plot for GWP
# for config, color in zip(configs_gwp, colors_gwp):
#     axes[1].errorbar(df[config], df.index, xerr=0.02, fmt='o', label=config.replace('_GWP',''), color=color)
# axes[1].axvline(0, color='black')
# axes[1].axvline(0.1, color='gray', linestyle='--')
# axes[1].axvline(-0.1, color='gray', linestyle='--')
# axes[1].set_title('GWP(Error bars: ±0.02)')
# axes[1].set_xlabel('Spearman ρ')
# axes[1].legend(title='Configuration')

plt.tight_layout()
plt.show()

