# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 13:52:32 2025

@author: aliah
"""
import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from shapely.geometry import Point

# -------------------------------------------------
# Setup paths
# -------------------------------------------------
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

excel_file = os.path.join(output_dir, "MSW_Type10_Revised_FoodFee_Savings.xlsx")
shapefile = os.path.join(output_dir, "tl_2023_us_state.shp")
out_map = os.path.join(output_dir, "FoodWaste_Savings_Map_USA_Green.png")

# -------------------------------------------------
# Load data
# -------------------------------------------------
df = pd.read_excel(excel_file)
df["State"] = df["State"].astype(str).str.upper().str.strip()

# Compute % savings if not present
if "% Savings on Food Waste Fee" not in df.columns:
    for c in ["Tipping Fee ($/ton)", "Total Wet (composition-based moisture)",
              "Food (dry tons)", "Paper (dry tons)", "Yard (dry tons)"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    denom = df["Food (dry tons)"] + df["Paper (dry tons)"] + df["Yard (dry tons)"]
    df["Food Fraction (composition-based)"] = np.where(denom > 0, df["Food (dry tons)"]/denom, np.nan)
    df["Food Waste (composition-based, wet tons)"] = df["Total Wet (composition-based moisture)"] * df["Food Fraction (composition-based)"]
    df["Revised Food Fee ($/ton)"] = 14.096
    old_cost = df["Tipping Fee ($/ton)"] * df["Food Waste (composition-based, wet tons)"] * 0.5 / 1e6
    new_cost = 14.096 * df["Food Waste (composition-based, wet tons)"] * 0.5 / 1e6
    df["% Savings on Food Waste Fee"] = (old_cost - new_cost) / old_cost * 100

# -------------------------------------------------
# Load shapefile
# -------------------------------------------------
gdf = gpd.read_file(shapefile)
exclude = {"PR", "GU", "VI", "AS", "MP"}
gdf = gdf[~gdf["STUSPS"].isin(exclude)].copy()

df["State"] = df["State"].astype(str).str.strip().str.upper()

# 2. If Excel has full state names, convert them to USPS abbreviations
us_state_abbrev = {
    'ALABAMA': 'AL', 'ALASKA': 'AK', 'ARIZONA': 'AZ', 'ARKANSAS': 'AR',
    'CALIFORNIA': 'CA', 'COLORADO': 'CO', 'CONNECTICUT': 'CT', 'DELAWARE': 'DE',
    'DISTRICT OF COLUMBIA': 'DC', 'FLORIDA': 'FL', 'GEORGIA': 'GA', 'HAWAII': 'HI',
    'IDAHO': 'ID', 'ILLINOIS': 'IL', 'INDIANA': 'IN', 'IOWA': 'IA',
    'KANSAS': 'KS', 'KENTUCKY': 'KY', 'LOUISIANA': 'LA', 'MAINE': 'ME',
    'MARYLAND': 'MD', 'MASSACHUSETTS': 'MA', 'MICHIGAN': 'MI', 'MINNESOTA': 'MN',
    'MISSISSIPPI': 'MS', 'MISSOURI': 'MO', 'MONTANA': 'MT', 'NEBRASKA': 'NE',
    'NEVADA': 'NV', 'NEW HAMPSHIRE': 'NH', 'NEW JERSEY': 'NJ', 'NEW MEXICO': 'NM',
    'NEW YORK': 'NY', 'NORTH CAROLINA': 'NC', 'NORTH DAKOTA': 'ND', 'OHIO': 'OH',
    'OKLAHOMA': 'OK', 'OREGON': 'OR', 'PENNSYLVANIA': 'PA', 'RHODE ISLAND': 'RI',
    'SOUTH CAROLINA': 'SC', 'SOUTH DAKOTA': 'SD', 'TENNESSEE': 'TN',
    'TEXAS': 'TX', 'UTAH': 'UT', 'VERMONT': 'VT', 'VIRGINIA': 'VA',
    'WASHINGTON': 'WA', 'WEST VIRGINIA': 'WV', 'WISCONSIN': 'WI', 'WYOMING': 'WY'
}

df["State"] = df["State"].replace(us_state_abbrev)

# 3. Filter to valid USPS codes
df = df[df["State"].isin(gdf["STUSPS"])].copy()

# Merge data
merged = gdf.merge(df, left_on="STUSPS", right_on="State", how="left")

# -------------------------------------------------
# Bin data and set palette
# -------------------------------------------------
merged["% Savings on Food Waste Fee"] = pd.to_numeric(merged["% Savings on Food Waste Fee"], errors="coerce")
bins = [0, 20, 40, 60, 80, np.inf]
labels = ["0–20%", "20–40%", "40–60%", "60–80%", "80%+"]
merged["Savings_bin"] = pd.cut(merged["% Savings on Food Waste Fee"], bins=bins, labels=labels, include_lowest=True)

color_mapping = {
    "0–20%": "#e5f5e0",
    "20–40%": "#c7e9c0",
    "40–60%": "#a1d99b",
    "60–80%": "#74c476",
    "80%+":  "#238b45",
}

# -------------------------------------------------
# Blog-style subplot layout
# -------------------------------------------------
fig = plt.figure(figsize=(18, 10))
ax_main = plt.subplot2grid((2, 2), (0, 0), colspan=2, fig=fig)
ax_ak = plt.subplot2grid((2, 2), (1, 0), fig=fig)
ax_hi = plt.subplot2grid((2, 2), (1, 1), fig=fig)

contiguous = merged[~merged["STUSPS"].isin(["AK", "HI"])]
alaska = merged[merged["STUSPS"] == "AK"]
hawaii = merged[merged["STUSPS"] == "HI"]

# -------------------------------------------------
# Plot each region
# -------------------------------------------------
contiguous.plot(ax=ax_main, color=contiguous["Savings_bin"].map(color_mapping),
                edgecolor="white", linewidth=0.5)
ax_main.set_xlim(-130, -65)
ax_main.set_ylim(24, 55)
ax_main.set_title(
    "Savings from Revised Food Waste Tipping Fee (Composition-Based)\n"
    "Percent reduction vs. current fees — by U.S. state",
    fontsize=16, fontweight="bold"
)
colors_ak = alaska["Savings_bin"].map(color_mapping).astype("object").fillna("#d9d9d9")
alaska.plot(ax=ax_ak, color=colors_ak, edgecolor="white", linewidth=0.5)
ax_ak.set_xlim(-200, -100)
ax_ak.set_ylim(50, 73)
ax_ak.set_title("Alaska", fontsize=11)
colors_hi = hawaii["Savings_bin"].map(color_mapping).astype("object").fillna("#d9d9d9")
hawaii.plot(ax=ax_hi, color=colors_hi, edgecolor="white", linewidth=0.5)
ax_hi.set_xlim(-162, -152)
ax_hi.set_ylim(18, 24)
ax_hi.set_title("Hawaii", fontsize=11)
from shapely import affinity

# --- Rescale and reposition Alaska and Hawaii ---
alaska.geometry = alaska.geometry.scale(xfact=0.32, yfact=0.32, origin=(0, 0))
alaska.geometry = alaska.geometry.translate(xoff=58, yoff=-52)

hawaii.geometry = hawaii.geometry.scale(xfact=0.55, yfact=0.55, origin=(0, 0))
hawaii.geometry = hawaii.geometry.translate(xoff=32, yoff=-17)

# Replot resized versions
# Replot resized versions safely (avoid Categorical fillna)
ax_ak.clear(); ax_hi.clear()

colors_ak = alaska["Savings_bin"].map(color_mapping).astype("object").fillna("#d9d9d9")
colors_hi = hawaii["Savings_bin"].map(color_mapping).astype("object").fillna("#d9d9d9")

alaska.plot(ax=ax_ak, color=colors_ak, edgecolor="white", linewidth=0.5)
hawaii.plot(ax=ax_hi, color=colors_hi, edgecolor="white", linewidth=0.5)

ax_ak.set_title("Alaska", fontsize=11)
ax_hi.set_title("Hawaii", fontsize=11)
for ax in (ax_ak, ax_hi):
    ax.set_axis_off()
for ax in (ax_ak, ax_hi):
    ax.set_axis_off()
# -------------------------------------------------
# Annotate state labels with arrows for small states
# -------------------------------------------------
from shapely.geometry import Point

# Compute projected centroids for better label placement
proj = merged.to_crs(epsg=5070).copy()
proj["centroid"] = proj.geometry.centroid
merged = merged.copy()
merged["centroid"] = proj["centroid"].to_crs(merged.crs)


# Define small, dense northeastern states
small_states = {"NH", "VT", "MA", "RI", "CT", "NJ", "DE", "MD", "DC"}

for _, row in merged.iterrows():
    st = row["STUSPS"]
    if st in ["AK", "HI"]:
        continue

    centroid = row["centroid"]
    if not isinstance(centroid, Point):
        continue

    x, y = centroid.x, centroid.y
    val = row["% Savings on Food Waste Fee"]

    if pd.notna(val):
        label = f"{st}\n{val:.0f}%"

        if st in small_states:
            # Offset placement for tiny northeastern states
            ax_main.annotate(
                label,
                xy=(x, y),
                xytext=(x + 2.2, y + 1.2),  # gentle diagonal offset
                fontsize=6.5,
                ha="left",
                va="center",
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.6, alpha=0.7),
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7)
            )
        else:
            # Normal in-map labels for all other states
            ax_main.text(
                x, y, label,
                fontsize=6.5,
                ha="center", va="center",
                color="black", fontweight="medium",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6)
            )

# -------------------------------------------------
# Legend and footer
# -------------------------------------------------
legend_handles = [mpatches.Patch(color=color_mapping[lbl], label=lbl) for lbl in labels]
fig.legend(handles=legend_handles, loc="lower center",
           bbox_to_anchor=(0.5, 0.02), ncol=len(labels),
           frameon=False, title="% Savings bins")

overall_mean = merged["% Savings on Food Waste Fee"].mean(skipna=True)
fig.text(0.5, 0.065,
         f"Mean savings across states: {overall_mean:.1f}%   •   Revised food fee = $14.096/ton   •   Data: your Excel + U.S. Census TIGER/Line 2023",
         ha="center", va="center", fontsize=10)

for ax in (ax_main, ax_ak, ax_hi):
    ax.set_axis_off()

plt.tight_layout(rect=[0, 0.09, 1, 1])
plt.savefig(out_map, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved blog-style map to: {out_map}")

