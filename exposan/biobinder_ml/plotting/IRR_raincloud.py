# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:03:48 2026

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
# Group order and colors
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

# ==============================
# Figure setup
# ==============================
fig, ax = plt.subplots(figsize=(14, 7), dpi=300)

base_positions = np.arange(len(feedstock_order)) * 1.85
offsets = [-0.52, -0.18, 0.18, 0.52]

violin_width = 0.22
box_width = 0.10
rain_shift = -0.06   # shift points to one side
cloud_shift = 0.04   # shift violin slightly to the other side

positions = []
data = []
meta = []

for i, feed in enumerate(feedstock_order):
    for off, combo in zip(offsets, combo_order):
        cfg, prod = combo

        vals = df[
            (df["Feedstock"] == feed) &
            (df["Config"] == cfg) &
            (df["Product set"] == prod)
        ]["IRR_pct"].dropna().to_numpy()

        if len(vals) == 0:
            continue

        x0 = base_positions[i] + off
        positions.append(x0)
        data.append(vals)
        meta.append((feed, cfg, prod, x0))

# ==============================
# Violin "cloud"
# ==============================
vp = ax.violinplot(
    data,
    positions=[x + cloud_shift for x in positions],
    widths=violin_width,
    showmeans=False,
    showmedians=False,
    showextrema=False,
)

for body, (_, cfg, prod, x0) in zip(vp["bodies"], meta):
    c = colors[(cfg, prod)]
    body.set_facecolor(c)
    body.set_edgecolor(c)
    body.set_alpha(0.80)
    body.set_linewidth(1.4)

    # clip left half to imitate half-violin / raincloud
    verts = body.get_paths()[0].vertices
    center = x0 + cloud_shift
    verts[verts[:, 0] < center, 0] = center

# ==============================
# Boxplot layer
# ==============================
bp = ax.boxplot(
    data,
    positions=positions,
    widths=box_width,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(linewidth=2.3),
    whiskerprops=dict(linewidth=1.9),
    capprops=dict(linewidth=1.9),
    boxprops=dict(linewidth=2.0),
)

for k, (_, cfg, prod, _) in enumerate(meta):
    c = colors[(cfg, prod)]

    bp["boxes"][k].set_facecolor("white")
    bp["boxes"][k].set_edgecolor(c)
    bp["boxes"][k].set_linewidth(2.0)

    bp["medians"][k].set_color(c)
    bp["medians"][k].set_linewidth(2.3)

    bp["whiskers"][2 * k].set_color(c)
    bp["whiskers"][2 * k + 1].set_color(c)
    bp["whiskers"][2 * k].set_linewidth(1.9)
    bp["whiskers"][2 * k + 1].set_linewidth(1.9)

    bp["caps"][2 * k].set_color(c)
    bp["caps"][2 * k + 1].set_color(c)
    bp["caps"][2 * k].set_linewidth(1.9)
    bp["caps"][2 * k + 1].set_linewidth(1.9)

# ==============================
# "Rain" points
# ==============================
rng = np.random.default_rng(42)

for feed, cfg, prod, x0 in meta:
    vals = df[
        (df["Feedstock"] == feed) &
        (df["Config"] == cfg) &
        (df["Product set"] == prod)
    ]["IRR_pct"].dropna().to_numpy()

    if len(vals) == 0:
        continue

    jitter = rng.uniform(-0.035, 0.035, size=len(vals))

    ax.plot(
        np.full(len(vals), x0 + rain_shift) + jitter,
        vals,
        "o",
        markersize=2.7,
        alpha=0.35,
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

ax.set_ylabel("IRR (%)", fontsize=13, fontweight="bold")
ax.set_xlabel("")

for tick in ax.get_yticklabels():
    tick.set_fontweight("bold")
    tick.set_fontsize(11)

ax.axhline(0, color="black", linewidth=1.2, alpha=0.6)
ax.grid(axis="y", alpha=0.22)

for x in (base_positions[:-1] + base_positions[1:]) / 2:
    ax.axvline(x, color="grey", linewidth=0.8, alpha=0.10)

# ==============================
# Legend
# ==============================
legend_handles = [
    Line2D([0], [0], color=colors[("CHCU","Biobinder")], lw=10, label="c-HTL – Biobinder"),
    Line2D([0], [0], color=colors[("CHCU","SAF")], lw=10, label="c-HTL – SAF"),
    Line2D([0], [0], color=colors[("DHCU","Biobinder")], lw=10, label="d-HTL – Biobinder"),
    Line2D([0], [0], color=colors[("DHCU","SAF")], lw=10, label="d-HTL – SAF"),
]

ax.legend(
    handles=legend_handles,
    loc="upper right",
    frameon=False,
    fontsize=11
)

plt.tight_layout()

# ==============================
# Save
# ==============================
plt.savefig("IRR_raincloud.png", dpi=300, bbox_inches="tight")
plt.savefig("IRR_raincloud.pdf", bbox_inches="tight")
plt.show()