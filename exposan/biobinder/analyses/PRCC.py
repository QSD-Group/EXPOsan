# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:27:48 2025

@author: aliah
"""

import os
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

# Create the results directory
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# Configurations and their file paths
config_files = {
    'CHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_results_CHCU_No_EC.xlsx'),
    'CHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_results_CHCU_EC.xlsx'),
    'DHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_results_DHCU_No_EC.xlsx'),
    'DHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_results_DHCU_EC.xlsx'),
}

# Function to calculate PRCC
def calculate_prcc(X, Y):
    residuals_X = {}
    residuals_Y = []

    # Calculate residuals of each parameter regressed on other parameters
    for param in X.columns:
        other_params = X.drop(columns=param)
        model = LinearRegression()
        model.fit(other_params, X[param])
        residuals_X[param] = X[param] - model.predict(other_params)

    # Calculate residuals of Y regressed on all parameters
    model_Y = LinearRegression()
    model_Y.fit(X, Y)
    residuals_Y = Y - model_Y.predict(X)

    # Calculate Spearman correlation between residuals
    prcc_values = {}
    for param in residuals_X:
        prcc, _ = spearmanr(residuals_X[param], residuals_Y)
        prcc_values[param] = prcc

    return prcc_values

# Process each configuration
for config_name, file_path in config_files.items():
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Load the data
    data = pd.read_excel(file_path)

    # Separate parameters and targets
    parameters = data.columns[:-2]  # Exclude MSP and GWP
    MSP = data['MSP ($/kg)']
    GWP = data['GWP (kg CO2e/kg)']

    # Calculate PRCC for MSP and GWP
    prcc_msp = calculate_prcc(data[parameters], MSP)
    prcc_gwp = calculate_prcc(data[parameters], GWP)

    # Combine PRCC results into a dataframe
    prcc_df = pd.DataFrame({
        'Parameter': parameters,
        'PRCC_MSP': [prcc_msp[param] for param in parameters],
        'PRCC_GWP': [prcc_gwp[param] for param in parameters]
    })

    # Save PRCC results to Excel
    prcc_file_path = os.path.join(output_dir, f'prcc_results_{config_name}.xlsx')
    prcc_df.to_excel(prcc_file_path, index=False)
    print(f"PRCC results saved to: {prcc_file_path}")

    # Generate Tornado Chart for MSP
    plt.figure(figsize=(10, 6))
    prcc_df.sort_values(by='PRCC_MSP', key=abs, ascending=False, inplace=True)
    plt.barh(prcc_df['Parameter'], prcc_df['PRCC_MSP'], color='royalblue')
    plt.xlabel('PRCC for MSP')
    plt.ylabel('Parameters')
    plt.title(f'Tornado Chart for {config_name} (MSP)')
    plt.tight_layout()

    # Save the tornado chart as an image
    tornado_chart_path = os.path.join(output_dir, f'tornado_chart_{config_name}_MSP.png')
    plt.savefig(tornado_chart_path)
    plt.close()
    print(f"Tornado chart saved to: {tornado_chart_path}")

