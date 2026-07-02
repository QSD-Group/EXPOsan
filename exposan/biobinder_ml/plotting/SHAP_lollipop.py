# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:49:52 2026

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# File path
# =========================================================
file_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\summary_outputs\modelB_IRR_shap_summary_20260330_1731.xlsx"

# =========================================================
# Load sheet
# =========================================================
df = pd.read_excel(file_path, sheet_name="cross_case_ranking")

# If ranking column exists, use it; otherwise sort by importance
if "mean_abs_shap_mean" not in df.columns:
    raise ValueError("Column 'mean_abs_shap_mean' not found in cross_case_ranking sheet.")

df = df.sort_values("mean_abs_shap_mean", ascending=True).reset_index(drop=True)

# Choose label column
if "Feature" in df.columns:
    feature_col = "Feature"
elif "feature" in df.columns:
    feature_col = "feature"
else:
    raise ValueError("Could not find feature name column ('Feature' or 'feature').")


# =========================================================
# Plot
# =========================================================
fig, ax = plt.subplots(figsize=(10, 6.5), dpi=300)

y = np.arange(len(df))
x = df["mean_abs_shap_mean"].values
labels = df[feature_col].astype(str).values

blue= "#1f77b4"

# stems (thicker)
ax.hlines(y=y, xmin=0, xmax=x, linewidth=3.5, color=blue)
# points (slightly bigger, filled)
ax.scatter(x, y, s=80, color=blue, edgecolor="none", zorder=3)
ax.axvline(0, color="black", linewidth=1.2, alpha=0.6)
# styling
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=11, fontweight="bold")
ax.set_xlabel("Mean absolute SHAP value", fontsize=13, fontweight="bold")
ax.set_ylabel("")

for tick in ax.get_xticklabels():
    tick.set_fontweight("bold")
    tick.set_fontsize(11)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(False)

plt.tight_layout()

# save
plt.savefig("cross_case_ranking_lollipop.png", dpi=300, bbox_inches="tight")
plt.savefig("cross_case_ranking_lollipop.pdf", bbox_inches="tight")

plt.show()