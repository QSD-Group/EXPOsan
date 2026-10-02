# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:52:18 2026

@author: aliah
"""

# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

# ============================================================
# INPUT
# ============================================================
htl_file = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\MSW_like_By_State_BillionTon2023_near-term_dry_tons_expanded (1).xlsx"
lf_file = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\EPA_Landfill_State_Summary_for_HTL.xlsx"

save_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\HTL_vs_active_landfills_circular_gap_labels.png"

plant_capacity_dry_tpd = 110
uptime_ratio = 0.90
plant_capacity_dry_tpy = plant_capacity_dry_tpd * 365 * uptime_ratio
usable_fraction = 0.50

# ============================================================
# LOAD DOE HTL DATA
# ============================================================
df = pd.read_excel(htl_file)

required_cols = ["State", "Food", "Paper", "Yard"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required DOE columns: {missing}")

df = df[required_cols].copy()

for c in ["Food", "Paper", "Yard"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["State", "Food", "Paper", "Yard"]).copy()

df["target_dry_tons"] = df["Food"] + df["Paper"] + df["Yard"]
df["usable_dry_tons"] = usable_fraction * df["target_dry_tons"]
df["plants_possible"] = df["usable_dry_tons"] / plant_capacity_dry_tpy
df["plants_integer"] = np.floor(df["plants_possible"]).astype(int)
df["HTL_capacity_dry_tpy"] = df["plants_integer"] * plant_capacity_dry_tpy

# ============================================================
# LOAD EPA LANDFILL SUMMARY
# ============================================================
lf = pd.read_excel(lf_file)

# Fix shifted workbook if needed
if "Unnamed: 0" in lf.columns and not lf["Unnamed: 0"].isna().all():
    if lf["Unnamed: 0"].astype(str).str.len().between(2, 2).any():
        lf = lf.rename(columns={
            "Unnamed: 0": "State_abbrev",
            "State": "Active_Landfills_fixed",
            "Active_Landfills": "Landfill_Tonnage_tpy_fixed",
            "Landfill_Tonnage_tpy": "Waste_in_Place_tons_fixed",
            "Waste_in_Place_tons": "Organic_Tonnage_tpy_51_4pct_fixed",
            "Organic_Tonnage_tpy_51_4pct": "Organic_Waste_in_Place_tons_51_4pct_fixed"
        })

        lf = lf[[
            "State_abbrev",
            "Active_Landfills_fixed",
            "Landfill_Tonnage_tpy_fixed",
            "Waste_in_Place_tons_fixed",
            "Organic_Tonnage_tpy_51_4pct_fixed",
            "Organic_Waste_in_Place_tons_51_4pct_fixed"
        ]].copy()

        lf.columns = [
            "State",
            "Active_Landfills",
            "Landfill_Tonnage_tpy",
            "Waste_in_Place_tons",
            "Organic_Tonnage_tpy_51_4pct",
            "Organic_Waste_in_Place_tons_51_4pct"
        ]

needed_lf_cols = [
    "State",
    "Active_Landfills",
    "Landfill_Tonnage_tpy",
    "Organic_Tonnage_tpy_51_4pct"
]

missing = [c for c in needed_lf_cols if c not in lf.columns]
if missing:
    raise ValueError(f"Missing landfill columns: {missing}")

lf = lf[needed_lf_cols].copy()

for c in ["Active_Landfills", "Landfill_Tonnage_tpy", "Organic_Tonnage_tpy_51_4pct"]:
    lf[c] = pd.to_numeric(lf[c], errors="coerce")

lf = lf.dropna(subset=["State"]).copy()

# ============================================================
# MERGE DOE POTENTIAL HTL WITH EPA LANDFILLS
# ============================================================
merged = pd.merge(df, lf, on="State", how="left")

merged["Active_Landfills"] = merged["Active_Landfills"].fillna(0)
merged["Landfill_Tonnage_tpy"] = merged["Landfill_Tonnage_tpy"].fillna(0)
merged["Organic_Tonnage_tpy_51_4pct"] = merged["Organic_Tonnage_tpy_51_4pct"].fillna(0)

merged["HTL_to_LF_count_ratio"] = (
    merged["plants_integer"] /
    merged["Active_Landfills"].replace(0, np.nan)
)

merged["HTL_capacity_to_organic_LF_ratio"] = (
    merged["HTL_capacity_dry_tpy"] /
    merged["Organic_Tonnage_tpy_51_4pct"].replace(0, np.nan)
)

# Sort by HTL plants for readable circular layout
merged = merged.sort_values("plants_possible", ascending=True).reset_index(drop=True)

# ============================================================
# PLOTTING TRANSFORM
# ============================================================
merged["HTL_plot_height"] = np.sqrt(merged["plants_possible"])
merged["LF_plot_height"] = np.sqrt(merged["Active_Landfills"])

labels = merged["State"].values
N = len(merged)

angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
bar_width = (2 * np.pi / N) * 0.82

htl_values_plot = merged["HTL_plot_height"].values
lf_values_plot = merged["LF_plot_height"].values

organic_tonnage = merged["Organic_Tonnage_tpy_51_4pct"].values / 1e6

# ============================================================
# COLOR BY ORGANIC LANDFILLED WASTE
# ============================================================
base_cmap = plt.cm.viridis
colors_trunc = base_cmap(np.linspace(0.30, 0.95, 256))
trunc_cmap = LinearSegmentedColormap.from_list("trunc_viridis", colors_trunc)

norm = plt.Normalize(np.nanmin(organic_tonnage), np.nanmax(organic_tonnage))
colors = trunc_cmap(norm(organic_tonnage))

lf_colors = colors.copy()
lf_colors[:, :3] = lf_colors[:, :3] * 0.85
lf_colors[:, 3] = 1.0

# ============================================================
# PLOT
# ============================================================
fig = plt.figure(figsize=(15, 15))
ax = plt.subplot(111, polar=True)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

inner_radius = 6.2
ring_gap = 1.25

lf_base = inner_radius
htl_base = inner_radius + lf_values_plot.max() + ring_gap

# Inner ring: active landfills
ax.bar(
    angles,
    lf_values_plot,
    width=bar_width,
    bottom=lf_base,
    color=lf_colors,
    alpha=0.95,
    edgecolor="white",
    linewidth=1.2,
    align="edge",
    # label="Active landfills"
)

# Outer ring: potential HTL plants
ax.bar(
    angles,
    htl_values_plot,
    width=bar_width,
    bottom=htl_base,
    color=colors,
    alpha=1.0,
    edgecolor="white",
    linewidth=1.4,
    align="edge",
    # label="Potential HTL plants"
)

ax.set_xticks([])
ax.set_yticks([])
ax.spines["polar"].set_visible(False)
ax.grid(False)

ax.set_ylim(0, htl_base + htl_values_plot.max() + 2.0)

# ============================================================
# LABELS
# ============================================================
state_font = 9.5
htl_font = 10
lf_font = 9.5
label_offset = 0.70
manual_offset = {
    "PA": 0.30,
    "GA": 0.30,
    "IL": 0.30,
    "OH": 0.45,
    "MI": 0.30,
    "NY": 0.30,
    "IN": 0.30,
    "VA": 0.30,
    "SC": 0.30,
    "MD": 0.30,
    "NV": 0.30,
    "MN": 0.45,
    "AR": 0.30,
    "MA": 0.30,
    "NE": 0.30,
    "WV": 0.30,
}


for i, angle in enumerate(angles):
    angle_mid = angle + bar_width / 2
    rotation = np.degrees(angle_mid)

    if 90 <= rotation <= 270:
        rotation += 180
        ha = "right"
    else:
        ha = "left"

    state = merged.loc[i, "State"]
    htl_int = int(merged.loc[i, "plants_integer"])
    lf_int = int(merged.loc[i, "Active_Landfills"])
    offset = label_offset + manual_offset.get(state, 0)

    # State label outside HTL ring
    ax.text(
        angle_mid,
        htl_base + htl_values_plot[i] + offset,
        state,
        rotation=rotation,
        rotation_mode="anchor",
        ha=ha,
        va="center",
        fontsize=state_font,
        fontweight="bold"
    )

    # HTL plant number inside outer ring
    ax.text(
        angle_mid,
        htl_base + htl_values_plot[i] * 0.55,
        f"{htl_int}",
        rotation=rotation,
        rotation_mode="anchor",
        ha="center",
        va="center",
        fontsize=htl_font,
        fontweight="bold",
        color="black"
    )

    # LF number in white gap outside inner ring
    if lf_int > 0:
        lf_number_radius = lf_base + lf_values_plot[i] + ring_gap * 0.45

        ax.text(
            angle_mid,
            lf_number_radius,
            f"{lf_int}",
            rotation=rotation,
            rotation_mode="anchor",
            ha="center",
            va="center",
            fontsize=lf_font,
            fontweight="heavy",
            color="brown"
        )

# ============================================================
# CENTER CIRCLE
# ============================================================
theta = np.linspace(0, 2 * np.pi, 600)
ax.fill(theta, np.full_like(theta, inner_radius - 0.25), color="white", zorder=3)
ax.plot(theta, np.full_like(theta, inner_radius), color="black", linewidth=2.0, zorder=4)

# ============================================================
# COLORBAR
# ============================================================
sm = ScalarMappable(cmap=trunc_cmap, norm=norm)
sm.set_array([])

cbar = plt.colorbar(sm, ax=ax, pad=0.015, shrink=0.62)
cbar.set_label("Landfilled organic waste, million tons/year", fontsize=12, fontweight= 'bold')
cbar.ax.tick_params(labelsize=11)

# ============================================================
# RING LABELS
# ============================================================
label_angle = np.deg2rad(75)   # move between 60–90 if needed

# Label for outer HTL ring
ax.text(
    label_angle,
    htl_base + htl_values_plot.max() * 0.55,
    "Potential HTL plants",
    rotation=15,
    ha="left",
    va="center",
    fontsize=12,
    fontweight="bold",
    color="black"
)

# Label for inner landfill ring / LF numbers
ax.text(
    label_angle,
    lf_base + lf_values_plot.max() * 0.75,
    "Landfill sites",
    rotation=15,
    ha="left",
    va="center",
    fontsize=12,
    fontweight="bold",
    color="brown"
)
# ============================================================
# LEGEND
# ============================================================
# ax.legend(
#     loc="lower center",
#     bbox_to_anchor=(0.5, -0.08),
#     ncol=2,
#     frameon=False,
#     fontsize=11
# )

plt.tight_layout()
# plt.subplots_adjust(
#     left=0.02,
#     right=0.92,
#     top=0.98,
#     bottom=0.05
# )
plt.savefig(save_path, dpi=600, bbox_inches="tight")
plt.savefig(save_path.replace(".png", ".pdf"), bbox_inches="tight")
plt.savefig(save_path.replace(".png", ".svg"), bbox_inches="tight")

plt.show()

# ============================================================
# SAVE TABLE
# ============================================================
table_path = save_path.replace(".png", "_merged_table.xlsx")
merged.to_excel(table_path, index=False)

print(f"Plant capacity, dry tons/year: {plant_capacity_dry_tpy:,.0f}")
print(f"Figure saved to: {save_path}")
print(f"PDF saved to: {save_path.replace('.png', '.pdf')}")
print(f"SVG saved to: {save_path.replace('.png', '.svg')}")
print(f"Table saved to: {table_path}")