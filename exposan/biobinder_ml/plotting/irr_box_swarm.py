# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:30:03 2026

@author: aliah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ==============================
# File path
# ==============================
file_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\master_IRR_summary_cleaned_final.xlsx"

# ==============================
# Load data
# ==============================
df = pd.read_excel(file_path, sheet_name="Long_All_Valid_IRR")

feedstock_order = ["food", "sludge", "manure", "green"]
config_order = ["CHCU", "DHCU"]
product_order = ["Biobinder", "SAF"]

df["Feedstock"] = pd.Categorical(df["Feedstock"], categories=feedstock_order, ordered=True)
df["Config"] = pd.Categorical(df["Config"], categories=config_order, ordered=True)
df["Product set"] = pd.Categorical(df["Product set"], categories=product_order, ordered=True)
df = df.sort_values(["Feedstock", "Config", "Product set"])

# ==============================
# Plot settings
# ==============================
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

fig, ax = plt.subplots(figsize=(13,7), dpi=300)

base_positions = np.arange(len(feedstock_order)) * 1.6
offsets = [-0.42, -0.14, 0.14, 0.42]
width = 0.22

positions = []
data = []
meta = []

# ==============================
# Prepare boxplot data
# ==============================
for i, feed in enumerate(feedstock_order):
    for off, combo in zip(offsets, combo_order):
        cfg, prod = combo

        vals = df[
            (df["Feedstock"] == feed) &
            (df["Config"] == cfg) &
            (df["Product set"] == prod)
        ]["IRR_pct"].dropna().values

        if len(vals) == 0:
            continue

        positions.append(base_positions[i] + off)
        data.append(vals)
        meta.append((feed, cfg, prod, base_positions[i] + off))

# ==============================
# Draw boxplots
# ==============================
bp = ax.boxplot(
    data,
    positions=positions,
    widths=width,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(linewidth=2.4),
    whiskerprops=dict(linewidth=2.0),
    capprops=dict(linewidth=2.0),
    boxprops=dict(linewidth=2.2),
)

# Apply colors and thicker outlines
for k, (_, cfg, prod, _) in enumerate(meta):
    c = colors[(cfg, prod)]

    bp["boxes"][k].set_facecolor(c)
    bp["boxes"][k].set_alpha(0.28)
    bp["boxes"][k].set_edgecolor(c)
    bp["boxes"][k].set_linewidth(2.2)

    bp["medians"][k].set_color(c)
    bp["medians"][k].set_linewidth(2.4)

    bp["whiskers"][2 * k].set_color(c)
    bp["whiskers"][2 * k + 1].set_color(c)
    bp["whiskers"][2 * k].set_linewidth(2.0)
    bp["whiskers"][2 * k + 1].set_linewidth(2.0)

    bp["caps"][2 * k].set_color(c)
    bp["caps"][2 * k + 1].set_color(c)
    bp["caps"][2 * k].set_linewidth(2.0)
    bp["caps"][2 * k + 1].set_linewidth(2.0)

# ==============================
# Add swarm-style jittered points
# ==============================
rng = np.random.default_rng(42)

for feed, cfg, prod, x0 in meta:
    vals = df[
        (df["Feedstock"] == feed) &
        (df["Config"] == cfg) &
        (df["Product set"] == prod)
    ]["IRR_pct"].dropna().values

    if len(vals) == 0:
        continue

    jitter = rng.uniform(-0.045, 0.045, size=len(vals))

    ax.plot(
        np.full(len(vals), x0) + jitter,
        vals,
        "o",
        markersize=3.5,
        alpha=0.45,
        color=colors[(cfg, prod)],
        markeredgewidth=0,
        zorder=3,
    )

# ==============================
# Axis styling
# ==============================
ax.set_xticks(base_positions)
ax.set_xticklabels(
    [f.capitalize() for f in feedstock_order],
    fontsize=12,
    fontweight="bold"
)

ax.set_ylabel(
    "IRR (%)",
    fontsize=13,
    fontweight="bold"
)

# No x-axis title
ax.set_xlabel("")

ax.axhline(0, color="black", linewidth=1.2, alpha=0.6)
ax.grid(axis="y", alpha=0.25)

# Optional: make y tick labels bold too
for tick in ax.get_yticklabels():
    tick.set_fontweight("bold")
    tick.set_fontsize(11)

# ==============================
# Legend
# ==============================
legend_handles = [
    Line2D([0], [0], color=colors[("CHCU", "Biobinder")], lw=8, alpha=0.5, label="c-HTL – Biobinder"),
    Line2D([0], [0], color=colors[("CHCU", "SAF")], lw=8, alpha=0.5, label="c-HTL – SAF"),
    Line2D([0], [0], color=colors[("DHCU", "Biobinder")], lw=8, alpha=0.5, label="d-HTL – Biobinder"),
    Line2D([0], [0], color=colors[("DHCU", "SAF")], lw=8, alpha=0.5, label="d-HTL – SAF"),
]

ax.legend(
    handles=legend_handles,
    loc="upper right",
    frameon=False,
    fontsize=11
)

plt.tight_layout()

# ==============================
# Save figure
# ==============================
plt.savefig("IRR_boxplots_swarm.png", dpi=300, bbox_inches="tight")
plt.savefig("IRR_boxplots_swarm.pdf", dpi=300, bbox_inches="tight")

plt.show()