# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:54:54 2026

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# File path
# =========================================================
file_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\summary_outputs\modelB_IRR_shap_direction_summary_20260330_1731.xlsx"

# =========================================================
# Load data
# =========================================================
df = pd.read_excel(file_path, sheet_name="cross_case_direction")

# =========================================================
# Sanity checks
# =========================================================
required_cols = ["feature_shap_corr", "direction_summary"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing column: {col}")

# detect feature column automatically
if "Feature" in df.columns:
    feature_col = "Feature"
elif "feature" in df.columns:
    feature_col = "feature"
else:
    raise ValueError("Feature column not found")

# =========================================================
# Sort (nice visual order)
# =========================================================
df = df.sort_values("feature_shap_corr", ascending=True).reset_index(drop=True)

# =========================================================
# Color mapping
# =========================================================
color_map = {
    "higher values tend to decrease IRR": "#E45756",  # red
    "higher values tend to increase IRR": "#54A24B",  # green
    "mixed / context-dependent": "#9D9D9D",           # gray
}

colors = df["direction_summary"].map(color_map)

rename_map = {
    "Feedstock_price_$/tonne": "Feedstock price",
    "product_ratio_biobinder_over_biofuel": "Binder/Fuel ratio",
    "Y_biocrude": "Biocrude yield",
    "Temperature (C)": "Temperature",
    "Solid content (w/w) %": "Solids content",
    "Electricity_Price": "Electricity price",
    "Natural_Gas_Price": "Natural gas price",
    "Uptime_Ratio": "Uptime",
    "income_tax": "Tax rate",
}

df[feature_col] = df[feature_col].replace(rename_map)

# =========================================================
# Plot
# =========================================================
fig, ax = plt.subplots(figsize=(10, 6.5), dpi=300)

y = np.arange(len(df))
x = df["feature_shap_corr"].values
labels = df[feature_col].astype(str).values

# stems
for yi, xi, ci in zip(y, x, colors):
    ax.hlines(y=yi, xmin=min(0, xi), xmax=max(0, xi),
              color=ci, linewidth=2)

# points
ax.scatter(x, y, s=80, c=colors, zorder=3)

# zero reference
ax.axvline(0, color="black", linewidth=1.2)

# labels
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=11, fontweight="bold")

ax.set_xlabel("Feature–SHAP correlation", fontsize=13, fontweight="bold")
ax.set_ylabel("")
# ax.set_title("Direction of feature effects on IRR", fontsize=14, fontweight="bold")

# axis styling
for tick in ax.get_xticklabels():
    tick.set_fontweight("bold")
    tick.set_fontsize(11)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(False)

# =========================================================
# Legend
# =========================================================
from matplotlib.lines import Line2D

legend_handles = [
    Line2D([0], [0], marker='o', color=color_map["higher values tend to decrease IRR"],
           linestyle='None', markersize=8, label='Higher values decrease IRR'),
    Line2D([0], [0], marker='o', color=color_map["higher values tend to increase IRR"],
           linestyle='None', markersize=8, label='Higher values increase IRR'),
    Line2D([0], [0], marker='o', color=color_map["mixed / context-dependent"],
           linestyle='None', markersize=8, label='Mixed / context-dependent')
]

ax.legend(handles=legend_handles, loc="lower right", frameon=False)

plt.tight_layout()

# save
plt.savefig("direction_lollipop_real_features.png", dpi=300, bbox_inches="tight")
plt.savefig("direction_lollipop_real_features.pdf", bbox_inches="tight")

plt.show()