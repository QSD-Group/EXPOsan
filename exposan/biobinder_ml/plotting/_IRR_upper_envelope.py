# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 12:24:09 2026

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Put the 5 META paths here (CHCU) ---
meta_paths = {
    "sludge": r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_sludge_CHCU_No_EC_prefer-max_top_ratio_20260203_1543.xlsx",
    "manure": r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_manure_CHCU_No_EC_prefer-max_top_ratio_20260203_1543.xlsx",
    "green":  r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_green_CHCU_No_EC_prefer-max_top_ratio_20260203_1543.xlsx",
    "food":   r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_food_CHCU_No_EC_prefer-max_top_ratio_20260204_1409.xlsx",
    "fog":    r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_fog_CHCU_No_EC_prefer-max_top_ratio_20260204_1409.xlsx",
}

def frontier_xy(df, xcol="product_ratio_biobinder_over_biofuel", ycol="IRR_pct", irr_cut=100.0):
    # Keep only OK + finite
    d = df[df["OK"] == True].copy()
    d = d[np.isfinite(d[ycol]) & np.isfinite(d[xcol])].copy()

    # Remove rows with IRR outside [-100, +100]
    d = d[(d[ycol].astype(float) > -irr_cut) & (d[ycol].astype(float) < irr_cut)].copy()

    # Pull x/y (product ratio can be >1; just require non-negative + finite)
    x = d[xcol].astype(float).values
    y = d[ycol].astype(float).values

    m = np.isfinite(x) & np.isfinite(y) & (x >= 0.0)
    x, y = x[m], y[m]

    if len(x) == 0:
        return None, None, None, None

    # Sort by x
    order = np.argsort(x)
    x, y = x[order], y[order]

    # Pareto frontier: best IRR achievable at/under each x (upper envelope)
    y_front = np.maximum.accumulate(y)
    keep = np.r_[True, y_front[1:] > y_front[:-1]]
    return x, y, x[keep], y_front[keep]

plt.figure(figsize=(8.5, 5.6))
plt.figure(dpi=300) 
for feedstock, path in meta_paths.items():
    df = pd.read_excel(path)

    x, y, x_f, y_f = frontier_xy(df)
    if x is None:
        print(f"⚠️ {feedstock}: no valid points after filtering.")
        continue

    # Scatter cloud (optional; comment out if you only want lines)
    plt.scatter(x, y, s=10, alpha=0.18)

    # Frontier line (one per feedstock)
    plt.plot(x_f, y_f, linewidth=2.3, label=feedstock)

plt.xlabel("product_ratio_biobinder_over_biofuel (binder / fuel)")
plt.ylabel("IRR (%)")
# plt.title("IRR vs product ratio frontier (filtered to -100 < IRR < 100)")
plt.grid(alpha=0.25)
plt.legend(title="Feedstock")
plt.tight_layout()
plt.show()
