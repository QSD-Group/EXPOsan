#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make two kinds of plots for bead lifetimes 0.5, 0.75, and 1.0 yr:

1) Separate contour-style heatmaps of USD_per_tonne_COD_removed
2) Separate stacked breakdown plots vs simulation rank
   - x-axis: simulation index after sorting from lowest to highest USD_per_tonne_COD_removed
   - y-axis: cost contribution in USD/tonne COD removed
   - positive groups above x-axis: Beads and Others
   - negative group below x-axis: Credits (e.g., biogas offset)

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
FILE = "/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/results/METAB_TEA_FB_1stage_passive_degas_30.xlsx"
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
df = pd.read_excel(FILE)

# Check required columns exist
required_cols = set(
    COST_COLS + ["bead_lifetime_yr", COD_REMOVED_COL, TOTAL_COST_COL, "HRT_d", "COD_g_L"]
)
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise KeyError(f"Missing required columns: {missing}")

# Make sure save directory exists
os.makedirs(SAVE_DIR, exist_ok=True)

# Fixed color scale across all heatmaps
vmin = df[TOTAL_COST_COL].min()
vmax = df[TOTAL_COST_COL].max()
levels = np.linspace(vmin, vmax, 20)

# ==========================
# MAKE BOTH PLOTS FOR EACH LIFETIME
# ==========================
for blt in BEAD_LIFETIMES:

    # -----------------------------------------
    # FILTER DATA FOR THIS BEAD LIFETIME
    # -----------------------------------------
    df_sub = df[np.isclose(df["bead_lifetime_yr"], blt)].copy()

    if df_sub.empty:
        print(f"No data found for bead lifetime = {blt}")
        continue

    # =========================================
    # 1) CONTOUR HEATMAP
    # =========================================
    heat = df_sub.pivot(
        index="HRT_d",
        columns="COD_g_L",
        values=TOTAL_COST_COL
    )

    heat = heat.sort_index().sort_index(axis=1)

    X, Y = np.meshgrid(heat.columns.values, heat.index.values)
    Z = heat.values

    fig, ax = plt.subplots(figsize=(7, 5.5))

    cf = ax.contourf(
        X, Y, Z,
        levels=levels,
        vmin=vmin,
        vmax=vmax,
        cmap="viridis"
    )

    ax.set_title(f"Bead lifetime = {blt:.2f} yr")
    ax.set_ylabel("HRT (d)")
    ax.set_xlabel("COD (g/L)")

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label(r"USD tonne$^{-1}$ COD removed")

    save_path = os.path.join(SAVE_DIR, f"FB_heatmap_bead_lifetime_{blt:.2f}yr.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {save_path}")

    # =========================================
    # 2) STACKED BREAKDOWN PLOT
    # =========================================

    # Convert annualized cost components to USD/tonne COD removed
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

    # Build plotting dataframe
    plot_df = df_sub[
        ["system_ID", "Q_m3_d", "COD_g_L", "HRT_d", "bead_lifetime_yr", TOTAL_COST_COL]
    ].copy()

    plot_df["Beads"] = beads.values
    plot_df["Others"] = others_pos.values
    plot_df["Credits"] = credits.values

    # Sort from lowest to highest total cost
    plot_df = plot_df.sort_values(TOTAL_COST_COL, ascending=True).reset_index(drop=True)

    # Simulation index
    x = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # Positive stack
    ax.bar(
        x,
        plot_df["Beads"].values,
        width=1.0,
        color=COLOR_BEADS,
        edgecolor="none",
        label="Beads"
    )

    ax.bar(
        x,
        plot_df["Others"].values,
        width=1.0,
        bottom=plot_df["Beads"].values,
        color=COLOR_OTHERS,
        edgecolor="none",
        label="Others"
    )

    # Negative stack below x-axis
    ax.bar(
        x,
        plot_df["Credits"].values,
        width=1.0,
        color=COLOR_CREDITS,
        edgecolor="none",
        label="Credits"
    )

    # Zero line
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_title(f"Bead lifetime = {blt:.2f} yr")
    ax.set_xlabel("Sample")
    ax.set_ylabel(r"USD tonne$^{-1}$ COD removed")
    ax.set_xlim(-0.5, len(x) - 0.5)

    legend_handles = [
        Patch(facecolor=COLOR_BEADS, edgecolor="none", label="Beads"),
        Patch(facecolor=COLOR_OTHERS, edgecolor="none", label="Others"),
        Patch(facecolor=COLOR_CREDITS, edgecolor="none", label="Credits (biogas offset)")
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False)

    save_path = os.path.join(SAVE_DIR, f"30_FB_breakdown_sorted_bead_lifetime_{blt:.2f}yr.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {save_path}")