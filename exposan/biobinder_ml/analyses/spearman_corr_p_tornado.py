# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:27:48 2025

@author: aliah
"""

import os
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

# Create the results directory
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# Configurations and their file paths
config_files = {
    'CHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_CHCU_No_EC.xlsx'),
    'CHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_CHCU_EC.xlsx'),
    'DHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_DHCU_No_EC.xlsx'),
    'DHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_DHCU_EC.xlsx'),
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

    # Calculate Spearman correlation and p-values between residuals
    prcc_values = {}
    p_values = {}
    for param in residuals_X:
        prcc, p_val = spearmanr(residuals_X[param], residuals_Y)
        prcc_values[param] = prcc
        p_values[param] = p_val

    return prcc_values, p_values

# Function to calculate Spearman correlation
def calculate_spearman(X, Y):
    spearman_values = {}
    p_values = {}
    for param in X.columns:
        spearman_corr, p_val = spearmanr(X[param], Y)
        spearman_values[param] = spearman_corr
        p_values[param] = p_val
    return spearman_values, p_values

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
    prcc_msp, prcc_pval_msp = calculate_prcc(data[parameters], MSP)
    prcc_gwp, prcc_pval_gwp = calculate_prcc(data[parameters], GWP)

    # Calculate Spearman for MSP and GWP
    spearman_msp, spearman_pval_msp = calculate_spearman(data[parameters], MSP)
    spearman_gwp, spearman_pval_gwp = calculate_spearman(data[parameters], GWP)

    # Combine results into a dataframe
    combined_df = pd.DataFrame({
        'Parameter': parameters,
        'PRCC_MSP': [prcc_msp[param] for param in parameters],
        'PRCC_MSP_p-value': [prcc_pval_msp[param] for param in parameters],
        'PRCC_GWP': [prcc_gwp[param] for param in parameters],
        'PRCC_GWP_p-value': [prcc_pval_gwp[param] for param in parameters],
        'Spearman_MSP': [spearman_msp[param] for param in parameters],
        'Spearman_MSP_p-value': [spearman_pval_msp[param] for param in parameters],
        'Spearman_GWP': [spearman_gwp[param] for param in parameters],
        'Spearman_GWP_p-value': [spearman_pval_gwp[param] for param in parameters]
    })

    # Save results to Excel
    results_file_path = os.path.join(output_dir, f'correlation_results_{config_name}.xlsx')
    combined_df.to_excel(results_file_path, index=False)
    print(f"Correlation results saved to: {results_file_path}")



    # # Create the tornado chart
    # plt.figure(figsize=(10, 6))
    # plt.barh(combined_df['Parameter'], combined_df['Spearman_MSP'], color='royalblue')
    # plt.xlabel('Spearman Correlation for MSP')
    # plt.ylabel('Parameters')
    # plt.title(f'Tornado Chart for {config_name} (Spearman - MSP)')
    # plt.tight_layout()

    # # Save the tornado chart as an image
    # tornado_chart_path = os.path.join(output_dir, f'tornado_chart_{config_name}_Spearman_MSP.png')
    # plt.savefig(tornado_chart_path)
    # plt.close()

    # print(f"Tornado chart saved to: {tornado_chart_path}")
    # Function to calculate Spearman correlation
    def calculate_spearman(X):
        spearman_corr_matrix = X.corr(method='spearman')
        return spearman_corr_matrix

    # Configurations and their file paths
    config_files = {
        'CHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_CHCU_No_EC.xlsx'),
        'CHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_CHCU_EC.xlsx'),
        'DHCU_No_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_DHCU_No_EC.xlsx'),
        'DHCU_EC': os.path.join(output_dir, 'uncertainty_analysis_results_1000_DHCU_EC.xlsx'),
    }

    # Process each configuration
    for config_name, file_path in config_files.items():
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        # Load the data
        data = pd.read_excel(file_path)

        # Separate parameters and MSP
        parameters = data.columns[:-2]  # Exclude MSP and GWP
        MSP = data['MSP ($/kg)']

        # Create a DataFrame with parameters and MSP only
        selected_data = data[list(parameters) + ['MSP ($/kg)']]

        # Compute Spearman correlation matrix
        corr_matrix = calculate_spearman(selected_data)

        # Create a mask for the upper triangle
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool))

        # Set up the matplotlib figure
        plt.figure(figsize=(18, 14))

        # Draw the heatmap with bold labels and increased font size
        ax = sns.heatmap(corr_matrix, mask=mask, annot=True, fmt=".2f", cmap='coolwarm',
                         vmin=-1, vmax=1, cbar_kws={'label': 'Spearman Correlation Coefficient'},
                         annot_kws={"size": 14})  # Normal font for correlation values
                         

        plt.gca().collections[0].colorbar.ax.tick_params(labelsize=14)  # Increase font size
        plt.gca().collections[0].colorbar.set_label('Spearman Correlation Coefficient', size=16, weight='bold')  # Set bold and bigger size
        # Format plot with bold axis labels
        plt.xticks([])
        plt.yticks(np.arange(len(selected_data.columns)) + 0.5, selected_data.columns, rotation=0, fontsize=16, fontweight='bold')

        # Save the heatmap as an image
        heatmap_path = os.path.join(output_dir, f'heatmap_{config_name}_MSP.png')
        plt.savefig(heatmap_path, bbox_inches='tight')
        plt.close()

        print(f"Heatmap saved to: {heatmap_path}")
 


    


