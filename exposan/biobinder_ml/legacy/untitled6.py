# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 15:25:59 2026

@author: aliah
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

meta_paths = {
    "fog": r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_fog_CHCU_No_EC_prefer-max_top_ratio_20260224_1626.xlsx",
    "sludge": r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_sludge_CHCU_No_EC_prefer-max_top_ratio_20260224_2014.xlsx",
    "manure":  r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_manure_CHCU_No_EC_prefer-max_top_ratio_20260224_2358.xlsx",
    "green":   r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_green_CHCU_No_EC_prefer-max_top_ratio_20260225_0345.xlsx",
    "food":    r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_food_CHCU_No_EC_prefer-max_top_ratio_20260226_1041.xlsx",
}

XCOL = "product_ratio_biobinder_over_biofuel"
YCOL = "IRR_pct"
ZCOL = "Y_biocrude"

def clean(df, irr_cut=100.0):
    d = df[df["OK"] == True].copy()
    d = d[np.isfinite(d[XCOL]) & np.isfinite(d[YCOL]) & np.isfinite(d[ZCOL])].copy()
    d[XCOL] = d[XCOL].astype(float)
    d[YCOL] = d[YCOL].astype(float)
    d[ZCOL] = d[ZCOL].astype(float)
    d = d[(d[YCOL] > -irr_cut) & (d[YCOL] < irr_cut) & (d[XCOL] >= 0)].copy()
    return d

# gather global yield min/max for consistent sizing
all_z = []
for p in meta_paths.values():
    df = pd.read_excel(p)
    d = clean(df)
    all_z.append(d[ZCOL].values)
zmin = float(np.nanmin(np.concatenate(all_z)))
zmax = float(np.nanmax(np.concatenate(all_z)))

def size_map(z, smin=12, smax=140):
    if zmax == zmin:
        return np.full_like(z, (smin + smax) / 2.0)
    return smin + (z - zmin) * (smax - smin) / (zmax - zmin)

# Consistent feedstock colors (match matplotlib defaults or choose your own)
feedstock_colors = {
    "sludge": "#1f77b4",
    "manure": "#ff7f0e",
    "green":  "#2ca02c",
    "food":   "#d62728",
    "fog":    "#9467bd",
}

fig, ax = plt.subplots(figsize=(8.5, 5.6), dpi=300)

# --- scatter plot ---
for feedstock, path in meta_paths.items():
    df = pd.read_excel(path)
    d = clean(df)
    s = size_map(d[ZCOL].values)
    ax.scatter(
        d[XCOL], d[YCOL],
        s=s, alpha=0.22,
        color=feedstock_colors[feedstock],
        label=feedstock
    )

ax.set_xlabel("binder / fuel")
ax.set_ylabel("IRR (%)")
ax.grid(alpha=0.25)

# Reference line at 10% IRR
ax.axhline(10.0, linewidth=1.2, alpha=0.9)
ax.text(ax.get_xlim()[0], 10.0, "  IRR = 10%", va="bottom", fontsize=9)

# --- Feedstock legend as LINES ---
line_handles = [
    Line2D([0], [0], color=feedstock_colors[k], lw=3)
    for k in meta_paths.keys()
]
line_labels = list(meta_paths.keys())

leg_feed = ax.legend(
    line_handles, line_labels,
    title="Feedstock",
    loc="upper right",
    frameon=True
)
ax.add_artist(leg_feed)

# --- Bubble-size legend (inside lower right) ---
ref_vals = [0.20, 0.35, 0.50, 0.65]  # adjust if needed
size_handles, size_labels = [], []

for rv in ref_vals:
    size_handles.append(
        ax.scatter(
            [], [],
            s=float(size_map(np.array([rv]))[0]),
            facecolors="none",
            edgecolors="gray",
            linewidths=1.0
        )
    )
    size_labels.append(f"Y_biocrude = {rv:.2f}")

leg_size = ax.legend(
    size_handles, size_labels,
    title="Biocrude yield",
    loc="lower right",
    frameon=True
)

fig.tight_layout()
plt.show()

