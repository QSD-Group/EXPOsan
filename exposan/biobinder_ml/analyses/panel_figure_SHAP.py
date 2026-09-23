# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# -----------------------------
# User settings
# -----------------------------
DATA_DIR = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results"
OUTPUT_FIG = "htl_shap_multiscale_figure.png"
OUTPUT_SUMMARY = "htl_shap_multiscale_summary.xlsx"

COMMON_FEEDSTOCKS = {"food", "green", "sludge", "manure"}

# -----------------------------
# File discovery
# -----------------------------
import re

def get_latest_rf_files(data_dir, common_feedstocks):
    all_rf = glob.glob(os.path.join(data_dir, "SHAP_RF_*.xlsx"))

    latest = {}
    for f in all_rf:
        name = os.path.basename(f)

        m = re.match(r"SHAP_RF_(?P<feedstock>[^_]+)_(?P<ts>\d{8}_\d{4})\.xlsx$", name)
        if not m:
            continue

        feedstock = m.group("feedstock").lower()
        ts = m.group("ts")

        if feedstock not in common_feedstocks:
            continue

        if feedstock not in latest or ts > latest[feedstock]["ts"]:
            latest[feedstock] = {"file": f, "ts": ts}

    return sorted(v["file"] for v in latest.values())


def get_latest_irr_files(data_dir, common_feedstocks):
    all_irr = glob.glob(os.path.join(data_dir, "SHAP3_SURROGATE_B_*_IRR_pct_*.xlsx"))

    latest = {}
    for f in all_irr:
        name = os.path.basename(f)

        m = re.match(
            r"SHAP3_SURROGATE_B_(?P<n>\d+)_(?P<feedstock>[^_]+)_(?P<scenario>.+?)_IRR_pct_.*?(?P<suffix>MERGED|\d{8}_\d{4})\.xlsx$",
            name
        )
        if not m:
            continue

        feedstock = m.group("feedstock").lower()
        scenario = m.group("scenario")
        suffix = m.group("suffix")

        if feedstock not in common_feedstocks:
            continue

        key = (feedstock, scenario)
        is_merged = (suffix == "MERGED")

        if key not in latest:
            latest[key] = {
                "file": f,
                "is_merged": is_merged,
                "suffix": suffix,
            }
        else:
            prev = latest[key]

            if is_merged and not prev["is_merged"]:
                latest[key] = {
                    "file": f,
                    "is_merged": True,
                    "suffix": suffix,
                }
            elif is_merged == prev["is_merged"]:
                if suffix > prev["suffix"]:
                    latest[key] = {
                        "file": f,
                        "is_merged": is_merged,
                        "suffix": suffix,
                    }

    return sorted(v["file"] for v in latest.values())


rf_files = get_latest_rf_files(DATA_DIR, COMMON_FEEDSTOCKS)
irr_files = get_latest_irr_files(DATA_DIR, COMMON_FEEDSTOCKS)

# -----------------------------
# Harmonization maps
# -----------------------------
process_to_system_map = {
    "Lipids wt%": "FS_factor_Lipid",
    "Protein wt%": "FS_factor_Prot",
    "Proteins wt%": "FS_factor_Prot",
    "Carbohydrates wt%": "FS_factor_Carb",
    "Ash wt%": "FS_factor_Ash",
    "Temperature (C)": "Temperature (C)",
    "Residence Time": "Residence Time",
    "Solid content (w/w) %": "Solid content (w/w) %",
    "Pre-processing": "Pre-processing",
    "Catalyst": "Catalyst",
    "Reactor Type": "Reactor Type",
    "Solvent": "Solvent",
    "Reactor Volume (mL)": "Reactor Volume (mL)",
    "C%": "C%",
    "H%": "H%",
    "O%": "O%",
    "N%": "N%",
    "HHV Biomass": "HHV Biomass",
}

display_label_map = {
    "FS_factor_Lipid": "Lipids",
    "FS_factor_Prot": "Proteins",
    "FS_factor_Carb": "Carbohydrates",
    "FS_factor_Ash": "Ash",
    "FS_factor_Water": "Water",
    "Temperature (C)": "Temperature",
    "Residence Time": "Residence time",
    "Solid content (w/w) %": "Solids loading",
    "Pre-processing": "Pre-processing",
    "Catalyst": "Catalyst",
    "Reactor Type": "Reactor type",
    "Solvent": "Solvent",
    "Reactor Volume (mL)": "Reactor volume",
    "C%": "C%",
    "H%": "H%",
    "O%": "O%",
    "N%": "N%",
    "HHV Biomass": "HHV biomass",
    "Biofuel_price": "Biofuel price",
    "Natural_Gas_Price": "Natural gas price",
    "Electricity_Price": "Electricity price",
    "Uptime_Ratio": "Uptime ratio",
    "income_tax": "Income tax",
    "product_ratio_biobinder_over_biofuel": "Product ratio",
    "Y_biocrude": "Biocrude yield",
}

# -----------------------------
# Helpers
# -----------------------------
def find_first_existing(df_columns, candidates):
    """Return first matching column name from candidates."""
    for c in candidates:
        if c in df_columns:
            return c
    return None

def load_sheet_with_columns(filepath, candidate_sets):
    """
    Search all sheets and return (df, sheet_name, resolved_cols)
    where resolved_cols maps logical names to actual column names.

    candidate_sets example:
    {
        "Feedstock": ["Feedstock", "feedstock"],
        "Feature": ["Feature", "feature"],
        "Biocrude_SHAP": ["Biocrude_SHAP", "biocrude_shap", "SHAP", "mean_abs_shap"],
    }
    """
    xls = pd.ExcelFile(filepath)
    for sheet in xls.sheet_names:
        tmp = pd.read_excel(filepath, sheet_name=sheet, nrows=5)
        resolved = {}
        ok = True
        for logical_name, candidates in candidate_sets.items():
            found = find_first_existing(tmp.columns, candidates)
            if found is None:
                ok = False
                break
            resolved[logical_name] = found
        if ok:
            df = pd.read_excel(filepath, sheet_name=sheet)
            return df, sheet, resolved

    raise ValueError(
        f"Could not find required columns in any sheet of {os.path.basename(filepath)}.\n"
        f"Looked for: {candidate_sets}\n"
        f"Available sheets: {xls.sheet_names}"
    )

# -----------------------------
# 1) Process layer: RF SHAP
# -----------------------------
# -----------------------------
# 1) Process layer: RF SHAP
# -----------------------------
proc_rows = []

rf_candidate_sets = {
    "Feedstock": ["Feedstock", "feedstock"],
    "Feature": ["Feature", "feature"],
    "Biocrude_SHAP": ["Biocrude_SHAP", "biocrude_shap", "SHAP", "mean_abs_shap"],
}

for f in rf_files:
    fname = os.path.basename(f).lower()

    # Skip files not in common feedstocks BEFORE opening them
    matched_feed = None
    for fs in COMMON_FEEDSTOCKS:
        if f"shap_rf_{fs}_" in fname:
            matched_feed = fs
            break

    if matched_feed is None:
        print(f"Skipping RF file not in COMMON_FEEDSTOCKS: {os.path.basename(f)}")
        continue

    try:
        df, sheet_name, cols = load_sheet_with_columns(f, rf_candidate_sets)
    except ValueError as e:
        print(f"Skipping RF file with unsupported format: {os.path.basename(f)}")
        print(e)
        continue

    work = df[[cols["Feedstock"], cols["Feature"], cols["Biocrude_SHAP"]]].copy()
    work.columns = ["Feedstock", "Feature", "Biocrude_SHAP"]

    # Trust filename feedstock more than sheet contents
    work["Feedstock"] = matched_feed

    work["Feature_std"] = work["Feature"].map(process_to_system_map)
    work = work[work["Feature_std"].notna()].copy()
    work["abs_shap"] = pd.to_numeric(work["Biocrude_SHAP"], errors="coerce").abs()
    work = work.dropna(subset=["abs_shap"])

    proc_rows.append(work[["Feedstock", "Feature", "Feature_std", "Biocrude_SHAP", "abs_shap"]])

if not proc_rows:
    raise RuntimeError("No RF files remained after filtering to common feedstocks.")
proc = pd.concat(proc_rows, ignore_index=True)

proc_summary = (
    proc.groupby("Feature_std", as_index=False)["abs_shap"]
    .mean()
    .rename(columns={"abs_shap": "process_abs_shap"})
)

proc_summary["Label"] = proc_summary["Feature_std"].map(display_label_map)\
    .fillna(proc_summary["Feature_std"])
# -----------------------------
# 2) System layer: surrogate IRR SHAP
# -----------------------------
irr_rows = []
dep_rows = []
r2_rows = []

irr_candidate_sets = {
    "Sample_ID": ["Sample_ID", "sample_id"],
    "Feedstock": ["Feedstock", "feedstock"],
    "Scenario": ["Scenario", "scenario", "Config", "config"],
    "Feature": ["Feature", "feature"],
    "Feature value": ["Feature value", "Feature_value", "feature_value"],
    "IRR_pct_SHAP": ["IRR_pct_SHAP", "IRR_SHAP", "shap_value"],
    "product_ratio_biobinder_over_biofuel": [
        "product_ratio_biobinder_over_biofuel"
    ],
    "Target_R2_test": ["Target_R2_test", "R2_test", "test_r2"],
}

for f in irr_files:
    df, sheet_name, cols = load_sheet_with_columns(f, irr_candidate_sets)

    work = df[[
        cols["Sample_ID"], cols["Feedstock"], cols["Scenario"], cols["Feature"],
        cols["Feature value"], cols["IRR_pct_SHAP"], cols["Target_R2_test"]
    ]].copy()

    work.columns = [
        "Sample_ID", "Feedstock", "Scenario", "Feature",
        "Feature value", "IRR_pct_SHAP", "Target_R2_test"
    ]

    feedstock = str(work["Feedstock"].iloc[0]).lower()
    if feedstock not in COMMON_FEEDSTOCKS:
        continue

    r2_rows.append({
        "file": os.path.basename(f),
        "sheet": sheet_name,
        "feedstock": feedstock,
        "scenario": work["Scenario"].iloc[0],
        "r2_test": float(work["Target_R2_test"].iloc[0]),
    })

    work["abs_shap"] = pd.to_numeric(work["IRR_pct_SHAP"], errors="coerce").abs()

    agg = (
        work.groupby(["Feedstock", "Scenario", "Feature"], as_index=False)["abs_shap"]
        .mean()
    )
    irr_rows.append(agg)

    pr = work[work["Feature"] == "product_ratio_biobinder_over_biofuel"][
        ["Sample_ID", "Feedstock", "Scenario", "Feature value", "IRR_pct_SHAP"]
    ].rename(columns={
        "Feature value": "product_ratio",
        "IRR_pct_SHAP": "product_ratio_shap"
    })

    lip = work[work["Feature"] == "FS_factor_Lipid"][
        ["Sample_ID", "Feedstock", "Scenario", "Feature value"]
    ].rename(columns={"Feature value": "lipid_factor"})

    dep_rows.append(pr.merge(lip, on=["Sample_ID", "Feedstock", "Scenario"], how="left"))

irr = pd.concat(irr_rows, ignore_index=True)
dep = pd.concat(dep_rows, ignore_index=True)
r2_df = pd.DataFrame(r2_rows)

irr_summary = (
    irr.groupby("Feature", as_index=False)["abs_shap"]
    .mean()
    .rename(columns={"Feature": "Feature_std", "abs_shap": "irr_abs_shap"})
)
irr_summary["Label"] = irr_summary["Feature_std"].map(display_label_map).fillna(irr_summary["Feature_std"])

# -----------------------------
# 3) Merge layers for driver shift panel
# -----------------------------
driver = pd.merge(proc_summary, irr_summary, on=["Feature_std", "Label"], how="outer").fillna(0.0)
driver["proc_norm"] = driver["process_abs_shap"] / max(driver["process_abs_shap"].max(), 1e-12)
driver["irr_norm"] = driver["irr_abs_shap"] / max(driver["irr_abs_shap"].max(), 1e-12)
driver["max_norm"] = driver[["proc_norm", "irr_norm"]].max(axis=1)

selected_driver_labels = [
    "Ash", "Lipids", "Carbohydrates", "Proteins",
    "Temperature", "Residence time", "Solids loading",
    "Pre-processing", "Reactor type", "Catalyst",
    "Product ratio", "Biocrude yield", "Natural gas price"
]
driver_plot = driver[driver["Label"].isin(selected_driver_labels)].copy()
driver_plot["Label"] = pd.Categorical(driver_plot["Label"], categories=selected_driver_labels[::-1], ordered=True)
driver_plot = driver_plot.sort_values("Label")

top_proc = proc_summary.sort_values("process_abs_shap", ascending=True).tail(7)
top_irr = irr_summary.sort_values("irr_abs_shap", ascending=True).tail(7)

# # -----------------------------
# # 4) Plot
# # -----------------------------
# fig = plt.figure(figsize=(14, 10))
# gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.35], width_ratios=[1.0, 1.15], hspace=0.28, wspace=0.22)

# # Panel A
# axA = fig.add_subplot(gs[0, 0])
# axA.axis("off")
# steps = [
#     ("Feedstock & process inputs", "composition, catalyst,\nreactor type, residence time"),
#     ("RF HTL prediction", "biocrude, aqueous,\ngas, char"),
#     ("Process simulation", "selective upgrading\nand separations"),
#     ("System outcome", "IRR across sampled\ntechno-economic space"),
# ]
# y_positions = [0.83, 0.58, 0.33, 0.08]
# box_colors = ["#dfe8f7", "#e6f4ea", "#fff2cc", "#f4e1f5"]

# for (title, subtitle), y, c in zip(steps, y_positions, box_colors):
#     box = FancyBboxPatch(
#         (0.08, y), 0.78, 0.14,
#         boxstyle="round,pad=0.02,rounding_size=0.03",
#         linewidth=1.2, edgecolor="0.35", facecolor=c,
#         transform=axA.transAxes
#     )
#     axA.add_patch(box)
#     axA.text(0.11, y + 0.09, title, fontsize=12, fontweight="bold",
#              transform=axA.transAxes, va="center")
#     axA.text(0.11, y + 0.04, subtitle, fontsize=10, color="0.25",
#              transform=axA.transAxes, va="center")
# for y1, y2 in zip(y_positions[:-1], y_positions[1:]):
#     axA.annotate(
#         "", xy=(0.47, y2 + 0.14), xytext=(0.47, y1),
#         xycoords=axA.transAxes, textcoords=axA.transAxes,
#         arrowprops=dict(arrowstyle="-|>", lw=1.3, color="0.35")
#     )
# axA.text(0.0, 1.02, "A", transform=axA.transAxes, fontsize=16, fontweight="bold")
# axA.text(0.06, 1.02, "Integrated pathway from process drivers to system IRR",
#          transform=axA.transAxes, fontsize=13, fontweight="bold")

# # Panel B
# axB = fig.add_subplot(gs[0, 1])
# y = np.arange(len(driver_plot))
# axB.scatter(np.zeros_like(y), y, s=800 * driver_plot["proc_norm"] + 20, alpha=0.75)
# axB.scatter(np.ones_like(y), y, s=800 * driver_plot["irr_norm"] + 20, alpha=0.75)
# axB.set_yticks(y)
# axB.set_yticklabels(driver_plot["Label"])
# axB.set_xticks([0, 1])
# axB.set_xticklabels(["Process\n(Biocrude SHAP)", "System\n(IRR SHAP)"])
# axB.set_xlim(-0.45, 1.45)
# axB.grid(axis="x", alpha=0.2)
# for spine in ["top", "right"]:
#     axB.spines[spine].set_visible(False)
# axB.set_title("B  Driver shift across scales", loc="left", fontweight="bold")
# for size, label in zip([0.25, 0.6, 1.0], ["Low", "Medium", "High"]):
#     axB.scatter([], [], s=800 * size + 20, label=label)
# axB.legend(title="Normalized mean |SHAP|", frameon=False, loc="lower right")

# # Panel C
# subgs = gs[1, 0].subgridspec(2, 1, hspace=0.38)
# axC1 = fig.add_subplot(subgs[0, 0])
# axC2 = fig.add_subplot(subgs[1, 0])

# axC1.barh(top_proc["Label"], top_proc["process_abs_shap"])
# axC1.set_title("C  Top process drivers of biocrude yield", loc="left", fontweight="bold")
# axC1.set_xlabel("Mean |Biocrude SHAP| across food, green, sludge")
# for spine in ["top", "right"]:
#     axC1.spines[spine].set_visible(False)

# axC2.barh(top_irr["Label"], top_irr["irr_abs_shap"])
# axC2.set_title("   Top system drivers of IRR", loc="left", fontweight="bold")
# axC2.set_xlabel("Mean |IRR SHAP| across food, green, sludge and CHCU/DHCU")
# for spine in ["top", "right"]:
#     axC2.spines[spine].set_visible(False)

# # Panel D
# axD = fig.add_subplot(gs[1, 1])
# dep_plot = dep.sample(min(5000, len(dep)), random_state=42)
# sc = axD.scatter(
#     dep_plot["product_ratio"], dep_plot["product_ratio_shap"],
#     c=dep_plot["lipid_factor"], s=12, alpha=0.5
# )
# axD.axhline(0, color="0.6", lw=1)
# axD.set_xlabel("Product ratio (biobinder / biofuel)")
# axD.set_ylabel("SHAP contribution to IRR from product ratio")
# axD.set_title("D  Nonlinear system response: product ratio vs IRR contribution", loc="left", fontweight="bold")
# for spine in ["top", "right"]:
#     axD.spines[spine].set_visible(False)
# cbar = fig.colorbar(sc, ax=axD, pad=0.02)
# cbar.set_label("FS_factor_Lipid")

# fig.suptitle("Multi-scale drivers of value creation in HTL-based waste valorization systems",
#              fontsize=16, fontweight="bold", y=0.98)
# fig.savefig(OUTPUT_FIG, dpi=220, bbox_inches="tight")

# # -----------------------------
# # 5) Save summary workbook
# # -----------------------------
# with pd.ExcelWriter(OUTPUT_SUMMARY, engine="openpyxl") as writer:
#     proc.to_excel(writer, sheet_name="RF_raw_common3", index=False)
#     proc_summary.to_excel(writer, sheet_name="RF_summary", index=False)
#     irr.to_excel(writer, sheet_name="IRR_summary_by_file", index=False)
#     irr_summary.to_excel(writer, sheet_name="IRR_summary", index=False)
#     driver.to_excel(writer, sheet_name="Driver_shift", index=False)
#     dep.to_excel(writer, sheet_name="Dependence_data", index=False)
#     r2_df.to_excel(writer, sheet_name="IRR_R2_by_file", index=False)

# print(f"Saved figure: {OUTPUT_FIG}")
# print(f"Saved summary workbook: {OUTPUT_SUMMARY}")


# ==============================
# Save Panel C separately
# ==============================
# figC, (axC1, axC2) = plt.subplots(
#     2, 1,
#     figsize=(6.5, 7),
#     dpi=300,
#     gridspec_kw={"hspace": 0.35}
# )

# # --- Top: Process drivers ---
# axC1.barh(top_proc["Label"], top_proc["process_abs_shap"])
# axC1.set_title("Top process drivers of biocrude yield", loc="left", fontweight="bold")
# axC1.set_xlabel("Mean |Biocrude SHAP|")

# for spine in ["top", "right"]:
#     axC1.spines[spine].set_visible(False)

# # --- Bottom: System drivers ---
# axC2.barh(top_irr["Label"], top_irr["irr_abs_shap"])
# axC2.set_title("Top system drivers of IRR", loc="left", fontweight="bold")
# axC2.set_xlabel("Mean |IRR SHAP|")

# for spine in ["top", "right"]:
#     axC2.spines[spine].set_visible(False)

# plt.tight_layout()

# # Save
# figC.savefig("Panel_C_drivers.png", dpi=300, bbox_inches="tight")
# figC.savefig("Panel_C_drivers.pdf", bbox_inches="tight")

# plt.show()
figC, (axC1, axC2) = plt.subplots(
    2, 1,
    figsize=(6.5, 7),
    dpi=300,
    gridspec_kw={"hspace": 0.45}
)

label_map = {
    "Feedstock_price_$/tonne": "Feedstock price",
    "product_ratio_biobinder_over_biofuel": "Product ratio",
    "Y_biocrude": "Biocrude yield",
    "Solid content (w/w) %": "Solids loading",
    "Natural_Gas_Price": "Natural gas price",
    "Electricity_Price": "Electricity price",
    "Biofuel_price": "Biofuel price",
    "income_tax": "Income tax",
    "Uptime_Ratio": "Uptime ratio",
    "Residence Time": "Residence time",
    "Temperature (C)": "Temperature",
}

top_proc["Label"] = top_proc["Label"].replace(label_map)
top_irr["Label"]  = top_irr["Label"].replace(label_map)
# Optional: reverse order so highest appears on top
# top_proc = top_proc.sort_values("process_abs_shap", ascending=True)
# top_irr  = top_irr.sort_values("irr_abs_shap", ascending=True)

# --- Top: Process drivers ---
axC1.barh(top_proc["Label"], top_proc["process_abs_shap"])
axC1.set_title(
    "Process-level drivers (HTL yield model)",
    loc="center",
    fontweight="bold"
)
axC1.set_xlabel("Mean |SHAP value| (biocrude yield)", fontsize=12, fontweight='bold')
axC1.tick_params(axis='both', labelsize=9)

for spine in ["top", "right"]:
    axC1.spines[spine].set_visible(False)

# --- Bottom: System drivers ---
axC2.barh(top_irr["Label"], top_irr["irr_abs_shap"])
axC2.set_title(
    "System-level drivers (surrogate TEA model)",
    loc="center",
    fontweight="bold"
)
axC2.set_xlabel("Mean |SHAP value| (IRR)", fontsize=12, fontweight='bold')
axC2.tick_params(axis='both', labelsize=9)

for spine in ["top", "right"]:
    axC2.spines[spine].set_visible(False)

for ax in [axC1, axC2]:
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')


# figC.text(
#     0.5, 0.51,
#     "Shift from chemistry-driven → economics-driven importance",
#     ha="center",
#     va="center",
#     fontsize=10,
#     style="italic"
# )

# # --- Optional SHAP comparability note ---
# figC.text(
#     0.01, 0.01,
#     "Note: SHAP magnitudes are not directly comparable across models.",
#     fontsize=8,
#     color="gray"
# )

plt.tight_layout()

# Save
figC.savefig("Panel_C_drivers.png", dpi=300, bbox_inches="tight")
figC.savefig("Panel_C_drivers.pdf", bbox_inches="tight")

plt.show()

# -----------------------------
# Export FULL driver table
# -----------------------------

origin_heatmap_table = driver.copy()

origin_heatmap_table = origin_heatmap_table[[
    "Label",
    "process_abs_shap",
    "irr_abs_shap",
    "proc_norm",
    "irr_norm"
]]

origin_heatmap_table.columns = [
    "Variable",
    "Mean_abs_SHAP_Biocrude_Yield",
    "Mean_abs_SHAP_IRR",
    "Normalized_SHAP_Biocrude_Yield",
    "Normalized_SHAP_IRR"
]

# sort by largest importance in either model
origin_heatmap_table["Max_SHAP"] = (
    origin_heatmap_table[
        ["Normalized_SHAP_Biocrude_Yield",
         "Normalized_SHAP_IRR"]
    ].max(axis=1)
)

origin_heatmap_table = (
    origin_heatmap_table
    .sort_values("Max_SHAP", ascending=False)
    .drop(columns="Max_SHAP")
)

output_excel = os.path.join(
    DATA_DIR,
    "Origin_SHAP_heatmap_table_FULL.xlsx"
)

origin_heatmap_table.to_excel(
    output_excel,
    index=False
)

print(f"Saved to:\n{output_excel}")

# -----------------------------
# Export LONG-FORM table for OriginPro Bubble + Color Mapped plot
# -----------------------------

# Keep normalized columns only
wide_for_origin = origin_heatmap_table[[
    "Variable",
    "Normalized_SHAP_Biocrude_Yield",
    "Normalized_SHAP_IRR"
]].copy()

# Optional: sort variables by combined importance before reshaping
wide_for_origin["Total_Importance"] = (
    wide_for_origin["Normalized_SHAP_Biocrude_Yield"]
    + wide_for_origin["Normalized_SHAP_IRR"]
)

wide_for_origin = wide_for_origin.sort_values(
    "Total_Importance",
    ascending=False
)

# Convert wide format to long format
origin_long = wide_for_origin.melt(
    id_vars=["Variable", "Total_Importance"],
    value_vars=[
        "Normalized_SHAP_Biocrude_Yield",
        "Normalized_SHAP_IRR"
    ],
    var_name="Model",
    value_name="Importance"
)

# Clean model names for figure labels
origin_long["Model"] = origin_long["Model"].replace({
    "Normalized_SHAP_Biocrude_Yield": "Yield Model",
    "Normalized_SHAP_IRR": "IRR Model"
})

# Bubble + color mapped in Origin needs separate size and color columns
origin_long["Bubble_Size"] = origin_long["Importance"]
origin_long["Bubble_Color"] = origin_long["Importance"]

# Preserve sorting order for variables
origin_long["Variable"] = pd.Categorical(
    origin_long["Variable"],
    categories=wide_for_origin["Variable"].tolist(),
    ordered=True
)

origin_long = origin_long.sort_values(
    ["Variable", "Model"]
)

# Save files
output_long_excel = os.path.join(
    DATA_DIR,
    "Origin_SHAP_bubble_long_table.xlsx"
)

output_long_csv = os.path.join(
    DATA_DIR,
    "Origin_SHAP_bubble_long_table.csv"
)

origin_long.to_excel(output_long_excel, index=False)
origin_long.to_csv(output_long_csv, index=False)

print("\nSaved OriginPro long-form bubble table:")
print(output_long_excel)
print(output_long_csv)

print("\nPreview:")
print(origin_long.head(20).to_string(index=False))