# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:53:15 2026

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

# ==============================
# File path
# ==============================
file_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\master_IRR_10000_GWP_summary_cleaned_final.xlsx"

# ==============================
# Load data
# ==============================
df = pd.read_excel(file_path, sheet_name="Long_All_Valid_IRR_GWP")

feedstock_order = ["food", "sludge", "manure", "green"]
combo_order = [
    ("CHCU", "Biobinder"),
    ("CHCU", "SAF"),
    ("DHCU", "Biobinder"),
    ("DHCU", "SAF"),
]

colors = {
    ("CHCU", "Biobinder"): "#4C78A8",
    ("CHCU", "SAF"): "#F58518",
    ("DHCU", "Biobinder"): "#54A24B",
    ("DHCU", "SAF"): "#E45756",
}

df["Feedstock"] = pd.Categorical(df["Feedstock"], categories=feedstock_order, ordered=True)
df = df.sort_values(["Feedstock", "Config", "Product set"])

# ==============================
# Figure setup
# ==============================
fig, ax = plt.subplots(figsize=(13.5, 7), dpi=300)

ax.grid(False)

# keep left and bottom axes, remove rectangle top/right
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    left=True,
    bottom=True,
    length=4,
    width=1
)

base_positions = np.arange(len(feedstock_order)) * 1.45
offsets = [-0.48, -0.16, 0.16, 0.48]

# geometry controls
halfeye_width = 0.14   # slab width
halfpoint_width = 0.06    # horizontal spread of points
halfeye_nudge = 0.065     # halfeye goes to the right
halfpoint_nudge = -0.065  # points go to the left
box_width = 0.12

rng = np.random.default_rng(42)

# ==============================
# Draw grouped rainclouds
# ==============================
for i, feed in enumerate(feedstock_order):
    for off, (cfg, prod) in zip(offsets, combo_order):
        vals = df[
            (df["Feedstock"] == feed) &
            (df["Config"] == cfg) &
            (df["Product set"] == prod)
        ]["IRR_pct"].dropna().to_numpy()

        if len(vals) == 0:
            continue

        x0 = base_positions[i] + off
        c = colors[(cfg, prod)]

        # --------------------------
        # 1) centered boxplot
        # --------------------------
        bp = ax.boxplot(
            [vals],
            positions=[x0],
            widths=box_width,
            vert=True,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(linewidth=2.2, color=c),
            whiskerprops=dict(linewidth=1.8, color=c),
            capprops=dict(linewidth=1.8, color=c),
            boxprops=dict(linewidth=1.8, color=c),
        )

        bp["boxes"][0].set_facecolor("white")
        bp["boxes"][0].set_edgecolor(c)
        bp["boxes"][0].set_linewidth(1.8)

        # --------------------------
        # 2) half-eye style density
        #    nudged to the right
        # --------------------------
        if len(vals) >= 2 and np.std(vals) > 0:
            y_grid = np.linspace(vals.min() - 0.05*(np.ptp(vals) + 1e-9),
                                 vals.max() + 0.05*(np.ptp(vals) + 1e-9), 300)

            kde = gaussian_kde(vals)
            kde.set_bandwidth(bw_method=kde.factor * 0.50)
            dens = kde(y_grid)
            dens = dens / dens.max() * halfeye_width

            ax.fill_betweenx(
                y_grid,
                x0 + halfeye_nudge,
                x0 + halfeye_nudge + dens,
                color=c,
                alpha=0.70,
                linewidth=0,
                zorder=1,
            )
            ax.plot(
                x0 + halfeye_nudge + dens,
                y_grid,
                color=c,
                linewidth=1.6,
                alpha=0.9,
                zorder=2
            )

        # --------------------------
        # 3) half-point raw data
        #    on left side only
        # --------------------------
        x_jitter = rng.uniform(
            x0 + halfpoint_nudge - halfpoint_width/2,
            x0 + halfpoint_nudge + halfpoint_width/2,
            size=len(vals)
        )

        ax.plot(
            x_jitter,
            vals,
            "o",
            markersize=2.8,
            alpha=0.45,
            color=c,
            markeredgewidth=0,
            zorder=3,
        )

# ==============================
# Axis styling
# ==============================
ax.set_xticks(base_positions)
ax.set_xticklabels(
    [f.capitalize() for f in feedstock_order],
    fontsize=14,
    fontweight="bold"
)

ax.set_ylabel("IRR (%)", fontsize=13, fontweight="bold")
ax.set_xlabel("")

for tick in ax.get_yticklabels():
    tick.set_fontweight("bold")
    tick.set_fontsize(14)

ax.axhline(0, color="black", linewidth=1.2, alpha=0.6)
# ax.grid(axis="y", alpha=0.22)

# for x in (base_positions[:-1] + base_positions[1:]) / 2:
#     ax.axvline(x, color="grey", linewidth=0.8, alpha=0.10)

# ==============================
# Legend
# ==============================
legend_handles = [
    Line2D([0], [0], color=colors[("CHCU", "Biobinder")], lw=10, label="c-HTL – Biobinder"),
    Line2D([0], [0], color=colors[("CHCU", "SAF")], lw=10, label="c-HTL – SAF"),
    Line2D([0], [0], color=colors[("DHCU", "Biobinder")], lw=10, label="d-HTL – Biobinder"),
    Line2D([0], [0], color=colors[("DHCU", "SAF")], lw=10, label="d-HTL – SAF"),
]

ax.legend(
    handles=legend_handles,
    loc="upper right",
    frameon=False,
    prop={"size": 14, "weight": "bold"}
)

plt.tight_layout(pad=0.5)
# ax.margins(x=0.02)
# ==============================
# Save
# ==============================
plt.savefig("IRR_vertical_halfeye_halfpoint_2.png", dpi=300, bbox_inches="tight")
plt.savefig("IRR_vertical_halfeye_halfpoint_2.pdf", bbox_inches="tight")
plt.show()