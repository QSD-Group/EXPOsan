# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# INPUT
# ============================================================
file_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\MSW_like_By_State_BillionTon2023_near-term_dry_tons_expanded (1).xlsx"
sheet_name = 0

plant_capacity_dry_tpd = 110
uptime_ratio = 0.90
plant_capacity_dry_tpy = plant_capacity_dry_tpd * 365 * uptime_ratio

save_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\HTL_plants_circular_bar_top5.png"

# ============================================================
# LOAD DATA
# ============================================================
df = pd.read_excel(file_path, sheet_name=sheet_name)

required_cols = ["State", "Food", "Paper", "Yard"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df = df[required_cols].copy()

for c in ["Food", "Paper", "Yard"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["State", "Food", "Paper", "Yard"]).copy()

# ============================================================
# CALCULATE PLANTS
# ============================================================
usable_fraction = 0.50

df["target_dry_tons"] = df["Food"] + df["Paper"] + df["Yard"]
df["usable_dry_tons"] = usable_fraction * df["target_dry_tons"]

df["plants_possible"] = df["usable_dry_tons"] / plant_capacity_dry_tpy
df["plants_integer"] = np.floor(df["plants_possible"]).astype(int)

df = df.sort_values("plants_possible", ascending=True).reset_index(drop=True)

df = df.sort_values("plants_possible", ascending=True).reset_index(drop=True)

# ============================================================
# TRANSFORM FOR PLOTTING ONLY
# ============================================================
df["plot_height"] = np.sqrt(df["plants_possible"])

values_real = df["plants_possible"].values
values_plot = df["plot_height"].values
labels = df["State"].values
N = len(df)

angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
inner_radius = 8.5
bar_width = (2 * np.pi / N) * 0.84

# ============================================================
# COLOR PALETTE
# ============================================================
# Use a warm palette closer to your sample figure
# norm = plt.Normalize(values_real.min(), values_real.max() * 0.85)
# colors = plt.cm.viridis(norm(values_real))
from matplotlib.colors import LinearSegmentedColormap

base_cmap = plt.cm.viridis
# colors_trunc = base_cmap(np.linspace(0.08, 0.98, 256))
colors_trunc = base_cmap(np.linspace(0.30, 0.95, 256))
# colors_trunc = base_cmap(np.linspace(0.25, 0.98, 256))
# colors_trunc = base_cmap(np.linspace(0.35, 0.90, 256))
# colors_trunc = base_cmap(np.linspace(0.30, 0.92, 256))
trunc_cmap = LinearSegmentedColormap.from_list("trunc_viridis", colors_trunc)

norm = plt.Normalize(values_real.min(), values_real.max())
colors = trunc_cmap(norm(values_real))

# from matplotlib.colors import LinearSegmentedColormap

# base_cmap = plt.cm.BuGn
# colors_trunc = base_cmap(np.linspace(0.15, 0.95, 256))
# trunc_cmap = LinearSegmentedColormap.from_list("trunc_BuGn", colors_trunc)

# norm = plt.Normalize(values_real.min(), values_real.max())
# colors = trunc_cmap(norm(values_real))

# ============================================================
# PLOT
# ============================================================
fig = plt.figure(figsize=(12, 12))
ax = plt.subplot(111, polar=True)
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

bars = ax.bar(
    angles,
    values_plot,
    width=bar_width,
    bottom=inner_radius,
    color=colors,
    edgecolor="white",
    linewidth=1.4,
    align="edge"
)

ax.set_xticks([])
ax.set_yticks([])
ax.spines["polar"].set_visible(False)
ax.grid(False)
ax.set_ylim(0, inner_radius + values_plot.max() + 3)

# ============================================================
# LABELS (ALL VALUES INSIDE — INTEGER)
# ============================================================
for i, (angle, value_real, value_plot, label) in enumerate(zip(angles, values_real, values_plot, labels)):
    angle_mid = angle + bar_width / 2
    rotation = np.degrees(angle_mid)

    if 90 <= rotation <= 270:
        rotation += 180
        ha = "right"
    else:
        ha = "left"

    # state label outside
    ax.text(
        angle_mid,
        inner_radius + value_plot + 0.9,
        label,
        rotation=rotation,
        rotation_mode="anchor",
        ha=ha,
        va="center",
        fontsize=9,
        fontweight="bold"
    )

    # integer value inside each bar
    value_int = df.loc[i, "plants_integer"]
    y_pos = inner_radius + value_plot * 0.55

    if value_plot < 2:
        fontsize = 11
        weight = "bold"
    elif value_plot < 4:
        fontsize = 11
        weight = "bold"
    else:
        fontsize = 11
        weight = "bold"

    ax.text(
        angle_mid,
        y_pos,
        f"{value_int}",
        rotation=rotation,
        rotation_mode="anchor",
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        color="black"
    )

# ============================================================
# CENTER CIRCLE
# ============================================================
theta = np.linspace(0, 2*np.pi, 600)
ax.fill(theta, np.full_like(theta, inner_radius - 0.25), color="white", zorder=3)
ax.plot(theta, np.full_like(theta, inner_radius), color="black", linewidth=2.0, zorder=4)

center_text = (
    "Potential HTL plants\n"
    "(landfilled part)\n"
    "110 dry tpd, 90% uptime"
)

ax.text(
    0, 0,
    center_text,
    ha="center",
    va="center",
    fontsize=20,
    weight = "bold",
    color='gray'
)

# plt.title(
#     "Number of HTL plants supported by state-level dry feedstock supply",
#     y=1.08,
#     fontsize=17,
#     fontweight="bold"
# )

plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches="tight")
plt.show()

# ============================================================
# OPTIONAL SAVE TABLE
# ============================================================
table_path = save_path.replace(".png", "_table.xlsx")
df.to_excel(table_path, index=False)

print(f"Plant capacity (dry tons/year): {plant_capacity_dry_tpy:,.0f}")
print(f"Figure saved to: {save_path}")
print(f"Table saved to: {table_path}")