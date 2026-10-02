# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 12:11:34 2026

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Put the 5 META paths here (CHCU) ---
meta_paths = {
    "sludge": r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_sludge_CHCU_No_EC_prefer-max_top_ratio_20260129_1601.xlsx",
    "manure": r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_manure_CHCU_No_EC_prefer-max_top_ratio_20260129_1601.xlsx",
    "green":  r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_green_CHCU_No_EC_prefer-max_top_ratio_20260129_1601.xlsx",
    "food":   r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_food_CHCU_No_EC_prefer-max_top_ratio_20260129_1601.xlsx",
    "fog":    r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_fog_CHCU_No_EC_prefer-max_top_ratio_20260129_1601.xlsx",
}

# --------------------------
# Helpers
# --------------------------
def clean_xy(df, xcol="product_ratio_biobinder_over_biofuel", ycol="IRR_pct",
             irr_cut=100.0):
    """Filter to OK runs and finite x/y; clip extreme IRR for plotting stability."""
    d = df[df["OK"] == True].copy()
    d = d[np.isfinite(d[xcol]) & np.isfinite(d[ycol])].copy()
    d[xcol] = d[xcol].astype(float)
    d[ycol] = d[ycol].astype(float)

    # Optional: clip IRR to avoid dominating the plot; keep symmetric range
    d = d[(d[ycol] > -irr_cut) & (d[ycol] < irr_cut)].copy()

    # binder/fuel ratio must be >= 0
    d = d[d[xcol] >= 0.0].copy()

    x = d[xcol].values
    y = d[ycol].values
    return x, y

def envelope_by_bins(x, y, nbins=45):
    """
    Empirical 'best-case' envelope: max IRR in each x-bin.
    Returns x_centers, y_max (finite only).
    """
    if len(x) == 0:
        return None, None

    xmin, xmax = float(np.nanmin(x)), float(np.nanmax(x))
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin == xmax:
        return None, None

    edges = np.linspace(xmin, xmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    y_max = np.full(nbins, np.nan)
    for i in range(nbins):
        m = (x >= edges[i]) & (x < edges[i + 1])
        if np.any(m):
            y_max[i] = np.nanmax(y[m])

    ok = np.isfinite(y_max)
    return centers[ok], y_max[ok]

def baseline_irr_from_low_x(x, y, low_quantile=0.05):
    """
    Baseline IRR = best achievable IRR in the lowest-x slice
    (proxy for fuel-max designs).
    """
    if len(x) == 0:
        return np.nan
    x0 = np.quantile(x, low_quantile)
    m = x <= x0
    if not np.any(m):
        # fallback: take best IRR among the smallest 1% points
        k = max(5, int(0.01 * len(x)))
        idx = np.argsort(x)[:k]
        return float(np.nanmax(y[idx]))
    return float(np.nanmax(y[m]))

# --------------------------
# Plot ΔIRR penalty curves
# --------------------------
plt.figure(figsize=(8.5, 5.6), dpi=300)

for feedstock, path in meta_paths.items():
    df = pd.read_excel(path)
    x, y = clean_xy(df)

    if len(x) == 0:
        print(f"⚠️ {feedstock}: no valid points after filtering.")
        continue

    # Baseline: "fuel-max proxy" = best IRR at very low binder ratio
    irr0 = baseline_irr_from_low_x(x, y, low_quantile=0.05)

    # Best-case envelope across binder ratios
    x_env, y_env = envelope_by_bins(x, y, nbins=55)
    if x_env is None:
        print(f"⚠️ {feedstock}: envelope could not be computed.")
        continue

    # ΔIRR penalty (typically <= 0 if binder is always worse)
    dIRR = y_env - irr0

    plt.plot(x_env, dIRR, linewidth=2.4, label=feedstock)

# Reference: no-penalty line
plt.axhline(0.0, linewidth=1.2, alpha=0.6)

plt.xlabel("product_ratio_biobinder_over_biofuel (binder / fuel)")
plt.ylabel("ΔIRR (%) relative to low-binder baseline")
plt.grid(alpha=0.25)
plt.legend(title="Feedstock")
plt.tight_layout()
plt.show()
