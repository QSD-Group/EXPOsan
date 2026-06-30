#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 23:12:32 2026

@author: saumitrarai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make paired comparison plots for bead lifetimes 0.5, 0.75, and 1.0 yr
using two Excel files:

1) Side-by-side contour-style heatmaps of USD_per_tonne_COD_removed
   - left: degas_30
   - right: degas_6.5
   - common color scale within each pair

2) Side-by-side stacked breakdown plots vs simulation rank
   - left: degas_30
   - right: degas_6.5
   - common y-limits within each pair

All figures are saved at 300 dpi.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ==========================
# USER INPUT
# ==========================
FILE_30 = "/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/results/METAB_TEA_FB_1stage_passive_degas_30.xlsx"
FILE_65 = "/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/results/METAB_TEA_FB_1stage_passive_degas_6.5.xlsx"
SAVE_DIR = "/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/figures"

BEAD_LIFETIMES = [0.5, 0.75, 1.0]

# Annualized cost columns to use in the breakdown
COST_COLS = [
    "ANNUAL__DMe_Module",
    "ANNUAL__DMe_Vacuum_pump_Rotary-vane_pump_one_stage",
    "ANNUAL__DMe_Water_pump_Motor",
    "ANNUAL__DMe_Water_pump_Pump",
    "ANNUAL__R1_Carbon_steel",
    "ANNUAL__R1_DoubleMembraneGasHolder_GH_Membrane",
    "ANNUAL__R1_DoubleMembraneGasHolder_GH_Slab_concrete",
    "ANNUAL__R1_HDPE_pipes",
    "ANNUAL__R1_IronSpongeTreatment_IST_Compressor",
    "ANNUAL__R1_IronSpongeTreatment_IST_Control_system",
    "ANNUAL__R1_IronSpongeTreatment_IST_Iron_sponge",
    "ANNUAL__R1_IronSpongeTreatment_IST_Vessel",
    "ANNUAL__R1_Pump_ins0_Motor",
    "ANNUAL__R1_Pump_ins0_Pump",
    "ANNUAL__R1_Pump_recirculation_Motor",
    "ANNUAL__R1_Pump_recirculation_Pump",
    "ANNUAL__R1_R1_beads",
    "ANNUAL__R1_Rockwool",
    "ANNUAL__R1_Slab_concrete",
    "ANNUAL__R1_Stainless_steel",
    "ANNUAL__R1_Wall_concrete",
    "ANNUAL__biogas_offset",
    "ANNUAL__chemicals",
    "ANNUAL__electricity",
    "ANNUAL__heat_onsite",
]

BEAD_COL = "ANNUAL__R1_R1_beads"
COD_REMOVED_COL = "COD_removed_tonne_yr"
TOTAL_COST_COL = "USD_per_tonne_COD_removed"

# Colors
COLOR_BEADS = "#f15a80"    # pink-ish
COLOR_OTHERS = "#53c68c"   # green-ish
COLOR_CREDITS = "#4c78a8"  # blue-ish

# ==========================
# LOAD DATA
# ==========================
df_30 = pd.read_excel(FILE_30)
df_65 = pd.read_excel(FILE_65)

# ==========================
# PRINT MEDIAN COSTS (FULL DATASET)
# ==========================
median_30 = df_30[TOTAL_COST_COL].median()
median_65 = df_65[TOTAL_COST_COL].median()

print("\n=== Overall median USD/tonne COD removed ===")
print(f"30  : {median_30:.2f}")
print(f"6.5 : {median_65:.2f}")

# Check required columns exist in both files
required_cols = set(
    COST_COLS + ["bead_lifetime_yr", COD_REMOVED_COL, TOTAL_COST_COL, "HRT_d", "COD_g_L"]
)

missing_30 = [c for c in required_cols if c not in df_30.columns]
missing_65 = [c for c in required_cols if c not in df_65.columns]

if missing_30:
    raise KeyError(f"Missing required columns in FILE_30: {missing_30}")
if missing_65:
    raise KeyError(f"Missing required columns in FILE_65: {missing_65}")

# Make sure save directory exists
os.makedirs(SAVE_DIR, exist_ok=True)


# ==========================
# HELPER FUNCTION FOR BREAKDOWN
# ==========================
def make_breakdown_df(df_sub):
    """
    Convert annualized cost components to USD/tonne COD removed
    and build sorted dataframe for stacked plotting.
    """
    contrib = df_sub[COST_COLS].div(df_sub[COD_REMOVED_COL], axis=0)

    # Beads kept separate
    beads = contrib[BEAD_COL].copy()

    # Everything except beads
    others_cols = [c for c in COST_COLS if c != BEAD_COL]
    others_matrix = contrib[others_cols]

    # Positive contributions above zero
    others_pos = others_matrix.clip(lower=0).sum(axis=1)

    # Negative contributions below zero
    credits = others_matrix.clip(upper=0).sum(axis=1)

    keep_cols = [TOTAL_COST_COL, "COD_g_L", "HRT_d", "bead_lifetime_yr"]
    if "system_ID" in df_sub.columns:
        keep_cols.insert(0, "system_ID")
    if "Q_m3_d" in df_sub.columns:
        keep_cols.append("Q_m3_d")

    plot_df = df_sub[keep_cols].copy()
    plot_df["Beads"] = beads.values
    plot_df["Others"] = others_pos.values
    plot_df["Credits"] = credits.values

    # Total positive stack height
    plot_df["Positive_total"] = plot_df["Beads"] + plot_df["Others"]

    # Sort from lowest to highest total cost
    plot_df = plot_df.sort_values(TOTAL_COST_COL, ascending=True).reset_index(drop=True)

    return plot_df


# ==========================
# MAKE PAIRED PLOTS FOR EACH LIFETIME
# ==========================
for blt in BEAD_LIFETIMES:

    # -----------------------------------------
    # FILTER DATA FOR THIS BEAD LIFETIME
    # -----------------------------------------
    df30_sub = df_30[np.isclose(df_30["bead_lifetime_yr"], blt)].copy()
    df65_sub = df_65[np.isclose(df_65["bead_lifetime_yr"], blt)].copy()

    if df30_sub.empty or df65_sub.empty:
        print(f"Skipping bead lifetime = {blt:.2f} yr because one dataset is empty.")
        continue

    # =========================================
    # 1) PAIRED CONTOUR HEATMAPS
    # =========================================
    heat30 = df30_sub.pivot(index="HRT_d", columns="COD_g_L", values=TOTAL_COST_COL)
    heat65 = df65_sub.pivot(index="HRT_d", columns="COD_g_L", values=TOTAL_COST_COL)

    heat30 = heat30.sort_index().sort_index(axis=1)
    heat65 = heat65.sort_index().sort_index(axis=1)

    # Common color scale across the pair
    pair_min = min(np.nanmin(heat30.values), np.nanmin(heat65.values))
    pair_max = max(np.nanmax(heat30.values), np.nanmax(heat65.values))
    levels = np.linspace(pair_min, pair_max, 20)

    X30, Y30 = np.meshgrid(heat30.columns.values, heat30.index.values)
    Z30 = heat30.values

    X65, Y65 = np.meshgrid(heat65.columns.values, heat65.index.values)
    Z65 = heat65.values

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True, sharey=True)

    cf1 = axes[0].contourf(
        X30, Y30, Z30,
        levels=levels,
        vmin=pair_min,
        vmax=pair_max,
        cmap="viridis"
    )
    axes[0].set_title(f"Degas = 30, bead lifetime = {blt:.2f} yr")
    axes[0].set_xlabel("COD (g/L)")
    axes[0].set_ylabel("HRT (d)")

    cf2 = axes[1].contourf(
        X65, Y65, Z65,
        levels=levels,
        vmin=pair_min,
        vmax=pair_max,
        cmap="viridis"
    )
    axes[1].set_title(f"Degas = 6.5, bead lifetime = {blt:.2f} yr")
    axes[1].set_xlabel("COD (g/L)")

    cbar = fig.colorbar(cf2, ax=axes, shrink=0.95)
    cbar.set_label(r"USD tonne$^{-1}$ COD removed")

    fig.suptitle(f"FB heatmaps, bead lifetime = {blt:.2f} yr", y=1.02)

    save_path = os.path.join(
        SAVE_DIR,
        f"FB_heatmap_pair_bead_lifetime_{blt:.2f}yr.png"
    )
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {save_path}")

    # =========================================
    # 2) PAIRED STACKED BREAKDOWN PLOTS
    # =========================================
    plot30 = make_breakdown_df(df30_sub)
    plot65 = make_breakdown_df(df65_sub)

    # Common y-limits across the pair
    y_max = max(plot30["Positive_total"].max(), plot65["Positive_total"].max())
    y_min = min(plot30["Credits"].min(), plot65["Credits"].min())

    # Small padding so bars do not touch boundaries
    y_pad_upper = 0.05 * y_max if y_max > 0 else 1.0
    y_pad_lower = 0.05 * abs(y_min) if y_min < 0 else 1.0

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

    # -------- Left panel: degas 30 --------
    x30 = np.arange(len(plot30))

    axes[0].bar(
        x30,
        plot30["Beads"].values,
        width=1.0,
        color=COLOR_BEADS,
        edgecolor="none"
    )
    axes[0].bar(
        x30,
        plot30["Others"].values,
        width=1.0,
        bottom=plot30["Beads"].values,
        color=COLOR_OTHERS,
        edgecolor="none"
    )
    axes[0].bar(
        x30,
        plot30["Credits"].values,
        width=1.0,
        color=COLOR_CREDITS,
        edgecolor="none"
    )
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_title(f"Degas = 30, bead lifetime = {blt:.2f} yr")
    axes[0].set_xlabel("Sample")
    axes[0].set_ylabel(r"USD tonne$^{-1}$ COD removed")
    axes[0].set_xlim(-0.5, len(x30) - 0.5)
    axes[0].set_ylim(y_min - y_pad_lower, y_max + y_pad_upper)

    # -------- Right panel: degas 6.5 --------
    x65 = np.arange(len(plot65))

    axes[1].bar(
        x65,
        plot65["Beads"].values,
        width=1.0,
        color=COLOR_BEADS,
        edgecolor="none"
    )
    axes[1].bar(
        x65,
        plot65["Others"].values,
        width=1.0,
        bottom=plot65["Beads"].values,
        color=COLOR_OTHERS,
        edgecolor="none"
    )
    axes[1].bar(
        x65,
        plot65["Credits"].values,
        width=1.0,
        color=COLOR_CREDITS,
        edgecolor="none"
    )
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_title(f"Degas = 6.5, bead lifetime = {blt:.2f} yr")
    axes[1].set_xlabel("Sample")
    axes[1].set_xlim(-0.5, len(x65) - 0.5)
    axes[1].set_ylim(y_min - y_pad_lower, y_max + y_pad_upper)

    legend_handles = [
        Patch(facecolor=COLOR_BEADS, edgecolor="none", label="Beads"),
        Patch(facecolor=COLOR_OTHERS, edgecolor="none", label="Others"),
        Patch(facecolor=COLOR_CREDITS, edgecolor="none", label="Credits (biogas offset)")
    ]
    axes[1].legend(handles=legend_handles, loc="upper left", frameon=False)

    fig.suptitle(f"FB cost breakdowns, bead lifetime = {blt:.2f} yr", y=1.02)

    save_path = os.path.join(
        SAVE_DIR,
        f"FB_breakdown_pair_bead_lifetime_{blt:.2f}yr.png"
    )
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {save_path}")