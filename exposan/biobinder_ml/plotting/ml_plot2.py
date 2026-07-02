# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 15:30:49 2026

@author: aliah
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ===============================================================
# File path
# ===============================================================
file = r"results/Results_10Fold_DNN_RF_GBR_Extended.xlsx"
out_dir = "results/bonus_defense_figures"
os.makedirs(out_dir, exist_ok=True)

# ===============================================================
# Load CV macro summary
# ===============================================================
cv_macro = pd.read_excel(file, sheet_name="CV_Macro_Summary")

# Keep only CV rows and selected models
model_map = {
    "DNN(MLP)": "DNN",
    "GradientBoosting": "GBR",
    "RandomForest": "RF"
}

plot_df = cv_macro[cv_macro["Model"].isin(model_map.keys())].copy()
plot_df["ModelLabel"] = plot_df["Model"].map(model_map)

# Metrics to compare
metrics = ["R2", "MAPE", "Pct_lt_5wt"]
metric_labels = ["R²", "MAPE (%)", "Coverage <5 wt%"]

# Reorder models
model_order = ["DNN", "GBR", "RF"]
plot_df = plot_df.set_index("ModelLabel").loc[model_order].reset_index()

# Build value matrix
vals = np.array([
    plot_df["R2"].values,
    plot_df["MAPE"].values,
    plot_df["Pct_lt_5wt"].values
]).T  # shape = (models, metrics)

# ===============================================================
# Plot
# ===============================================================
x = np.arange(len(metrics))
width = 0.22

fig, ax = plt.subplots(figsize=(9, 5.5))

for i, model in enumerate(model_order):
    y = [plot_df.loc[plot_df["ModelLabel"] == model, m].values[0] for m in metrics]
    bars = ax.bar(x + (i - 1) * width, y, width=width, label=model)
    for b, v in zip(bars, y):
        ax.text(
            b.get_x() + b.get_width()/2,
            b.get_height() + 0.8 if v > 5 else b.get_height() + 0.01,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=11
        )

ax.set_xticks(x)
ax.set_xticklabels(metric_labels, fontsize=12, fontweight="bold")
ax.set_ylabel("Metric value", fontsize=13, fontweight="bold")
ax.set_title("10-fold CV comparison: RF vs GBR vs DNN", fontsize=15, fontweight="bold")
ax.legend(frameon=False, fontsize=11)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()

save_path = os.path.join(out_dir, "CV_macro_RF_vs_GBR_vs_DNN.png")
plt.savefig(save_path, dpi=300, bbox_inches="tight")
plt.show()

print(f"✅ Saved: {save_path}")