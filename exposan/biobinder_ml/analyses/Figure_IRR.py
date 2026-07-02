# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 13:36:30 2025

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Config
# -----------------------------
DPI = 300
FIGSIZE = (8, 6)
OUTLIER_MIN_IRR = -30   # remove < this
OUTLIER_MAX_IRR = 90    # remove > this

# -----------------------------
# Helper: filtered dataframe
# -----------------------------
def load_and_filter(path):
    df = pd.read_excel(path, sheet_name="Sheet1")  #load sheet
    df = df.dropna(subset=["IRR (%)", "Y_Bio", "Product Ratio"]).copy()
    df = df[(df["IRR (%)"] > OUTLIER_MIN_IRR) & (df["IRR (%)"] < OUTLIER_MAX_IRR)]
    return df

# -----------------------------
# Plot function with styling
# -----------------------------
def plot_bubble(df, out_path):
    plt.figure(figsize=FIGSIZE)

    # Bubble size (scaled); color encodes Product Ratio
    bubble_size = df["Product Ratio"].values * 300

    sc = plt.scatter(
        df["Y_Bio"].values,
        df["IRR (%)"].values,
        s=bubble_size,
        c=df["Product Ratio"].values,   # color by Product Ratio
        cmap="viridis",                 # same style as earlier
        alpha=0.8,
        marker="o",
        edgecolors="black",
        linewidths=0.5,
    )

    # Axis labels (bold, 18 pt)
    plt.xlabel("Biocrude Yield", fontsize=18, fontweight="bold")
    plt.ylabel("IRR (%)", fontsize=18, fontweight="bold")

    # Colorbar with bold label/ticks
    cbar = plt.colorbar(sc)
    cbar.set_label("Product Ratio", fontsize=16, fontweight="bold")
    cbar.ax.tick_params(labelsize=14, width=1.5)
    for lab in cbar.ax.get_yticklabels():
        lab.set_fontweight("bold")

    # No title
    plt.title("")

    # Thick axes
    ax = plt.gca()
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)

    # Bold tick labels
    plt.xticks(fontsize=14, fontweight="bold")
    plt.yticks(fontsize=14, fontweight="bold")

    plt.tight_layout()
    plt.savefig(out_path, dpi=DPI)
    plt.close()

# -----------------------------
# Paths (update if needed)
# -----------------------------
dhcu_path = "Ybio_sweep_Lr98_dynamicCutoff_FixedLHK_DHCU_No_EC_20250813_1207.xlsx"
chcu_path = "Ybio_sweep_Lr98_dynamicCutoff_FixedLHK_CHCU_No_EC_20250813_1216.xlsx"

# -----------------------------
# Generate both plots with same styling
# -----------------------------
dhcu_df = load_and_filter(dhcu_path)
chcu_df = load_and_filter(chcu_path)

plot_bubble(dhcu_df, "bubble_DHCU_No_EC_300dpi_bold.png")
plot_bubble(chcu_df, "bubble_CHCU_No_EC_300dpi_bold.png")
