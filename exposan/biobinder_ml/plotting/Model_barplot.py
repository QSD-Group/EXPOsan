# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 17:46:28 2025

@author: aliah
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Load test-set results ===
file_path = "results/Results_10Fold_DNN_Softmax_RF_GBR_20251022_1728.xlsx"
df = pd.read_excel(file_path, sheet_name="TestSet_Results")

# Rename for clarity
df["Model"] = df["Model"].replace({
    "DNN(Softmax)_Test": "DNN",
    "RandomForest_Test": "RF",
    "GradientBoosting_Test": "GBR"
})

df["Yield"] = df["Target"].replace({
    "Biocrude wt%": "BioY%",
    "Aqueous wt%": "AqY%",
    "Gas wt%": "Gas%",
    "Solids wt%": "Solids%"
})

# Remove any macro rows
if "__macro__" in df["Yield"].values:
    df = df[df["Yield"] != "__macro__"]

# === Define metrics and colors ===
metrics = ["R2", "RMSE", "MAE", "MAPE"]
palette = {"DNN": "#1f77b4", "RF": "#ff7f0e", "GBR": "#2ca02c"}

# === Output folder ===
os.makedirs("results/testset_barplots", exist_ok=True)

# === Generate one plot per metric ===
for m in metrics:
    plt.figure(figsize=(6, 4))
    ax = sns.barplot(
        data=df, x="Yield", y=m, hue="Model",
        palette=palette, edgecolor="black", linewidth=0.8
    )
    
    # Add numeric labels
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", label_type="edge", fontsize=9, padding=2)
    
    # Reverse axis for error metrics
    if m in ["RMSE", "MAE", "MAPE"]:
        ax.invert_yaxis()
    
    # Styling: bold axis labels & ticks
    ax.set_xlabel("Yield Type", fontsize=12, fontweight="bold")
    ax.set_ylabel(m, fontsize=12, fontweight="bold")
    ax.tick_params(axis='x', labelsize=11, width=1.2)
    ax.tick_params(axis='y', labelsize=10, width=1.2)
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')


    # Legend and title
    ax.legend(title="Model", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
    # plt.title(f"{m} Across Models for Each Yield", fontsize=14, fontweight="bold")
    
    plt.tight_layout()
    plt.savefig(f"results/testset_barplots/{m}_GroupedBar_TestSet_Bold.png", dpi=300, bbox_inches="tight")
    plt.show()

