# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:40:48 2025

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi
import os

# ------------------------------------------------------------------
# 1️⃣ Load SHAP summary Excel
# ------------------------------------------------------------------
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

file_path = os.path.join(output_dir, "Feature_Importance_SHAP.xlsx")
df = pd.read_excel(file_path, sheet_name="Raw_SHAP")

output_dir = "results_radar_short"
os.makedirs(output_dir, exist_ok=True)

# Short labels map
short_labels = {
    "Temperature (C)": "T(°C)",
    "Lipids wt%": "Lip%",
    "Protein wt%": "Prot%",
    "Carbohydrates wt%": "Carb%",
    "Ash wt%": "Ash%",
    "Solid content (w/w) %": "Sol%",
    "Reactor Volume (mL)": "RV",
    "Residence Time": "t_res",
    "Pre-processing": "PreP",
    "Catalyst": "Cat",
    "Solvent": "Solv",
    "Reactor Type": "RType",
    "HHV Biomass": "HHV",
    "C%": "C%",
    "H%": "H%",
    "O%": "O%",
    "N%": "N%"
}

models_full = ["DNN(MLP)", "RandomForest", "GradientBoosting"]
models_short = {"DNN(MLP)": "DNN", "RandomForest": "RF", "GradientBoosting": "GBR"}
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
yield_list = ["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"]
N_top = 6

# ------------------------------------------------------------
# 3️⃣ Radar plotting function
# ------------------------------------------------------------
def plot_radar_for_yield(yield_name):
    subset = df[df["Yield"] == yield_name]
    pivot = subset.pivot_table(values="Mean|SHAP|", index="Feature", columns="Model", aggfunc="mean").fillna(0)
    pivot["avg"] = pivot.mean(axis=1)
    pivot = pivot.sort_values("avg", ascending=False).head(N_top).drop(columns="avg")
    pivot_norm = pivot / pivot.max()

    features = [short_labels.get(f, f) for f in pivot_norm.index]
    N = len(features)
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    plt.figure(figsize=(7,7))
    ax = plt.subplot(111, polar=True)
    plt.ylim(0, 1)

    # Plot each model polygon
    for model, color in zip(models_full, colors):
        if model not in pivot_norm.columns:
            continue
        vals = pivot_norm[model].tolist() + [pivot_norm[model].tolist()[0]]
        ax.plot(angles, vals, linewidth=2.5, label=models_short[model], color=color)
        ax.fill(angles, vals, color=color, alpha=0.25)

    # --- Feature labels ---
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(features, fontsize=14, fontweight="bold", color="black")
    ax.tick_params(axis='x', pad=12)

    # --- Radial rings ---
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"],
                       fontsize=12, fontweight="bold", color="gray")
    ax.spines["polar"].set_visible(True)
    ax.grid(color="gray", linestyle="dotted", linewidth=0.8)
    ax.xaxis.set_tick_params(labelsize=13)

    # 🧩 Remove degree labels only (keep features)
    ax.set_thetagrids([])  # remove auto-deg lines
    ax.set_xticks(angles[:-1])  # reapply custom ticks
    ax.set_xticklabels(features, fontsize=14, fontweight="bold", color="black")

    # --- Title and legend ---
    # plt.title(f"Feature Importance — {yield_name}", fontsize=17, fontweight='bold', y=1.10)
    legend = plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.05),
                        fontsize=13, frameon=False, labelspacing=0.3)
    for t in legend.get_texts():
        t.set_fontweight('bold')

    plt.tight_layout()
    fname = f"Radar_{yield_name.replace('%','pct').replace(' ','_')}_cleanrings3.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved {fname}")

# ------------------------------------------------------------
# 4️⃣ Generate all yields
# ------------------------------------------------------------

for y in yield_list:
    plot_radar_for_yield(y)

print("\n🎯 All radar plots saved in 'results_radar_short/' (short labels, bold fonts, 300 dpi).")
