# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 14:45:35 2026

@author: aliah
"""
import os
import numpy as np
import pandas as pd

# ============================================================
# CONFIG
# ============================================================
combined_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\COMBINED_YIELD_TEA_SHAP_biocrude_WITH_KNOB_NORMALIZATION.xlsx"
# ^ change filename if yours differs

out_long = os.path.join(os.path.dirname(combined_path), "SUMMARY_KNOBNORM_LONG.xlsx")
out_wide = os.path.join(os.path.dirname(combined_path), "SUMMARY_KNOBNORM_WIDE.xlsx")

TOP_N = 5  # top knobs to list per feedstock

KNOBS = [
    "Temperature (C)",
    "Residence Time",
    "Solid content (w/w) %",
    "Catalyst",
    "Reactor Type",
    "Pre-processing",
    "Solvent",
]
KNOB_SET = set(KNOBS)

# ============================================================
# LOAD
# ============================================================
xl = pd.ExcelFile(combined_path)
# Prefer ALL_FEEDSTOCKS if it exists (best), else build from per-feedstock sheets.
if "ALL_FEEDSTOCKS" in xl.sheet_names:
    df = pd.read_excel(combined_path, sheet_name="ALL_FEEDSTOCKS")
else:
    # Concatenate everything except SUMMARY (if present)
    frames = []
    for sh in xl.sheet_names:
        if sh.upper() in ("SUMMARY",):
            continue
        tmp = pd.read_excel(combined_path, sheet_name=sh)
        # If no Feedstock column, infer from sheet name
        if "Feedstock" not in tmp.columns:
            tmp.insert(0, "Feedstock", sh)
        frames.append(tmp)
    df = pd.concat(frames, ignore_index=True)

# Basic checks
need_cols = {"Feedstock", "Feature", "Yield_KnobNorm_%", "MSP_KnobNorm_%"}
missing = need_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in combined file: {missing}\nColumns found: {list(df.columns)}")

# Keep knobs only
dfk = df[df["Feature"].astype(str).isin(KNOB_SET)].copy()

# Clean numeric
for c in ["Yield_KnobNorm_%", "MSP_KnobNorm_%"]:
    dfk[c] = pd.to_numeric(dfk[c], errors="coerce")

dfk = dfk[np.isfinite(dfk["Yield_KnobNorm_%"]) | np.isfinite(dfk["MSP_KnobNorm_%"])].copy()

# Fill NaNs with 0 for easier pivot/plotting (if a knob is missing)
dfk["Yield_KnobNorm_%"] = dfk["Yield_KnobNorm_%"].fillna(0.0)
dfk["MSP_KnobNorm_%"] = dfk["MSP_KnobNorm_%"].fillna(0.0)

# Add helpful derived columns
dfk["Delta_MSP_minus_Yield"] = dfk["MSP_KnobNorm_%"] - dfk["Yield_KnobNorm_%"]
dfk["Ratio_MSP_div_Yield"] = dfk["MSP_KnobNorm_%"] / dfk["Yield_KnobNorm_%"].replace(0.0, np.nan)

# Optional: force consistent knob order
dfk["Feature"] = pd.Categorical(dfk["Feature"], categories=KNOBS, ordered=True)

# ============================================================
# WIDE PIVOTS
# ============================================================
yield_wide = dfk.pivot_table(index="Feature", columns="Feedstock", values="Yield_KnobNorm_%", aggfunc="mean").reindex(KNOBS)
msp_wide   = dfk.pivot_table(index="Feature", columns="Feedstock", values="MSP_KnobNorm_%", aggfunc="mean").reindex(KNOBS)
delta_wide = dfk.pivot_table(index="Feature", columns="Feedstock", values="Delta_MSP_minus_Yield", aggfunc="mean").reindex(KNOBS)
ratio_wide = dfk.pivot_table(index="Feature", columns="Feedstock", values="Ratio_MSP_div_Yield", aggfunc="mean").reindex(KNOBS)

# ============================================================
# TOP-N TABLES PER FEEDSTOCK
# ============================================================
top_rows = []
for feed, g in dfk.groupby("Feedstock"):
    g2 = g.copy()

    top_y = g2.sort_values("Yield_KnobNorm_%", ascending=False).head(TOP_N)
    for _, r in top_y.iterrows():
        top_rows.append({
            "Feedstock": feed,
            "RankingType": f"TOP{TOP_N}_Yield_KnobNorm",
            "Feature": r["Feature"],
            "Percent": float(r["Yield_KnobNorm_%"]),
        })

    top_m = g2.sort_values("MSP_KnobNorm_%", ascending=False).head(TOP_N)
    for _, r in top_m.iterrows():
        top_rows.append({
            "Feedstock": feed,
            "RankingType": f"TOP{TOP_N}_MSP_KnobNorm",
            "Feature": r["Feature"],
            "Percent": float(r["MSP_KnobNorm_%"]),
        })

df_top = pd.DataFrame(top_rows).sort_values(["Feedstock", "RankingType", "Percent"], ascending=[True, True, False])

# ============================================================
# SAVE
# ============================================================
with pd.ExcelWriter(out_long, engine="openpyxl") as w:
    dfk.sort_values(["Feedstock", "Feature"]).to_excel(w, sheet_name="KNOBNORM_LONG", index=False)

with pd.ExcelWriter(out_wide, engine="openpyxl") as w:
    yield_wide.to_excel(w, sheet_name="YIELD_KNOBNORM")
    msp_wide.to_excel(w, sheet_name="MSP_KNOBNORM")
    delta_wide.to_excel(w, sheet_name="DELTA_MSP_MINUS_YIELD")
    ratio_wide.to_excel(w, sheet_name="RATIO_MSP_DIV_YIELD")
    df_top.to_excel(w, sheet_name="TOPK_PER_FEEDSTOCK", index=False)

print("✅ Saved:")
print(out_long)
print(out_wide)
