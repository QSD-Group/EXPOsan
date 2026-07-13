# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

# data path local
data_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\Billion Ton"

# Map the file paths dynamically
files = {
    'max_top_ratio': os.path.join(data_dir, 'county_biobinder_cutoff_sensitivity_max_top_ratio_20260710_1712.xlsx'),
    'min_err': os.path.join(data_dir, 'county_biobinder_cutoff_sensitivity_min_err_20260710_1727.xlsx'),
    'max_Hr_then_max_Lr': os.path.join(data_dir, 'county_biobinder_cutoff_sensitivity_max_Hr_then_max_Lr_20260711_1151.xlsx'),
    'max_Lr_then_max_Hr': os.path.join(data_dir, 'county_biobinder_cutoff_sensitivity_max_Lr_then_max_Hr_20260711_1209.xlsx')
}

compiled_records = []

policy_mapping = {
    'max_top_ratio': 'Max Biofuel\n(max_top_ratio)',
    'min_err': 'Baseline Emulation\n(min_err)',
    'max_Hr_then_max_Lr': 'Max HK Recovery\n(max_Hr_then_max_Lr)',
    'max_Lr_then_max_Hr': 'Max LK Recovery\n(max_Lr_then_max_Hr)'
}

# Gather and combine data from files
for policy, filepath in files.items():
    df_summary = pd.read_excel(filepath, sheet_name='National Cut Summary')
    for idx, row in df_summary.iterrows():
        scenario = row['Availability_case']
        scenario_clean = 'All Waste Available' if 'All Waste' in scenario else 'Only Landfilled Portion'
        interval = int(row['Cutoff_frac_medium'] * 100)
        
        policy_clean = policy_mapping[policy]
        
        for fs, col_name in [('Food Waste', 'Food_binder_dry_ton_per_year'), 
                             ('Green Waste', 'Green_binder_dry_ton_per_year'), 
                             ('Manure Waste', 'Manure_binder_dry_ton_per_year')]:
            compiled_records.append({
                'Policy': policy_clean,
                'Scenario': scenario_clean,
                'Interval_Num': interval,
                'Quality Interval': f"{interval}%",
                'Feedstock': fs,
                'Binder (Million tons/yr)': row.get(col_name, 0) / 1e6
            })

df_stacked = pd.DataFrame(compiled_records)

# Define design parameters
scenarios = ['All Waste Available', 'Only Landfilled Portion']
policies = [
    'Max Biofuel\n(max_top_ratio)', 
    'Max HK Recovery\n(max_Hr_then_max_Lr)', 
    'Max LK Recovery\n(max_Lr_then_max_Hr)', 
    'Baseline Emulation\n(min_err)'
]

intervals_all = [30, 40, 50, 60, 70]
feedstocks = ['Food Waste', 'Green Waste', 'Manure Waste']
colors = ['#E06D43', '#28A745', '#4A2E1B'] 

fig, axes = plt.subplots(len(scenarios), len(policies), figsize=(16, 7.5), sharex=False, sharey='row')

# Base theme setup
sns.set_theme(style="white", context="notebook")

# Plot each grid tile
for s_idx, scenario in enumerate(scenarios):
    for p_idx, policy in enumerate(policies):
        ax = axes[s_idx, p_idx]
        df_sub = df_stacked[(df_stacked['Scenario'] == scenario) & (df_stacked['Policy'] == policy)]
        
        plot_data = {fs: [] for fs in feedstocks}
        valid_intervals = []
        
        for i in intervals_all:
            df_int = df_sub[df_sub['Interval_Num'] == i]
            if not df_int.empty:
                valid_intervals.append(f"{i}%")
                for fs in feedstocks:
                    val = df_int[df_int['Feedstock'] == fs]['Binder (Million tons/yr)'].values[0]
                    plot_data[fs].append(val)
        
        if valid_intervals:
            bottom = [0] * len(valid_intervals)
            for fs_idx, fs in enumerate(feedstocks):
                ax.bar(valid_intervals, plot_data[fs], bottom=bottom, color=colors[fs_idx], 
                       edgecolor=colors[fs_idx], label=fs, width=0.5)
                bottom = [b + v for b, v in zip(bottom, plot_data[fs])]
        
        if s_idx == 0:
            ax.set_title(policy, fontsize=13, fontweight='bold')
        if p_idx == 0:
            ax.set_ylabel(f"{scenario}\nBinder Yield (MM short tons/yr)", fontsize=13, fontweight='bold', labelpad=12)
        else:
            ax.set_ylabel("") 
        
        ax.set_xlabel("Biocrude Medium Cut", fontsize=13, fontweight='bold')
        ax.set_xticklabels(valid_intervals, fontsize=12)
        ax.tick_params(labelsize=11)
        
        # 1. REMOVE ALL GRID LINES Completely
        ax.grid(False)
        
        # 2. RESTORE CRISP BLACK AXIS LINES (BORDERS)
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(1.0)
            
        # 3. RESTORE THE TICK MARKS ON THE AXIS OUTLINES
        ax.tick_params(axis='both', which='both', direction='out', length=4, width=1, left=True, bottom=True, colors='black')

# Extract legend elements and insert 'Feedstock Type:' as an inline label
handles, labels = axes[0, 0].get_legend_handles_labels()
dummy_patch = mpatches.Patch(color='none', label='Feedstock Type:')
handles = [dummy_patch] + handles
labels = ['Feedstock Type:'] + labels

# Build single line inline legend
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=4, 
           frameon=True, handletextpad=0.5, columnspacing=1.5, prop={'size': 12, 'weight': 'bold'})

plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.show()