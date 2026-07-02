# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:27:48 2025

@author: aliah
"""

import os
import pandas as pd
from scipy.stats import spearmanr

# Directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# Configuration and file paths
config_files = {
    'CHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_lhs_1000_CHCU_No_EC_20250521_1430.xlsx'),
    'CHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_lhs_1000_CHCU_EC_20250521_1430.xlsx'),
    'DHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_lhs_1000_DHCU_No_EC_20250521_1430.xlsx'),
    'DHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_lhs_1000_DHCU_EC_20250521_1430.xlsx'),
}

# Spearman correlation calculator
def calculate_spearman(X, Y):
    spearman_values = {}
    p_values = {}
    for param in X.columns:
        corr, p_val = spearmanr(X[param], Y)
        spearman_values[param] = corr
        p_values[param] = p_val
    return spearman_values, p_values

# Loop through each configuration
for config_name, file_path in config_files.items():
    if not os.path.exists(file_path):
        print(f" File not found: {file_path}")
        continue

    print(f" Processing {config_name}...")

    data = pd.read_excel(file_path)

    # Automatically identify parameter columns (excluding MSP and GWP)
    param_cols = [col for col in data.columns if col not in ['MSP ($/kg)', 'GWP (kg CO2e/kg)']]
    X = data[param_cols]
    MSP = data['MSP ($/kg)']
    GWP = data['GWP (kg CO2e/kg)']

    # Compute Spearman correlations
    spearman_msp, spearman_pval_msp = calculate_spearman(X, MSP)
    spearman_gwp, spearman_pval_gwp = calculate_spearman(X, GWP)

    # Combine results
    results_df = pd.DataFrame({
        'Parameter': param_cols,
        'Spearman_MSP': [spearman_msp[p] for p in param_cols],
        'Spearman_MSP_p-value': [spearman_pval_msp[p] for p in param_cols],
        'Spearman_GWP': [spearman_gwp[p] for p in param_cols],
        'Spearman_GWP_p-value': [spearman_pval_gwp[p] for p in param_cols],
    })

    # Save
    out_file = os.path.join(output_dir, f'spearman_correlation_{config_name}.xlsx')
    results_df.to_excel(out_file, index=False)
    print(f"✅ Saved: {out_file}")

    


