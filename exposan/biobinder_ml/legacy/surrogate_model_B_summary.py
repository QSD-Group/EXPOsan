# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 16:10:29 2026

@author: aliah
"""

import os
import re
import glob
import numpy as np
import pandas as pd

# -----------------------------
# Settings
# -----------------------------
INPUT_DIR = r"."   # change to your local folder
OUTPUT_DIR = os.path.join(INPUT_DIR, "summary_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

pattern = os.path.join(INPUT_DIR, "SHAP3_SURROGATE_B_*_IRR_pct_*.xlsx")
files = sorted(glob.glob(pattern))

if not files:
    raise FileNotFoundError(f"No files found matching: {pattern}")

# -----------------------------
# Helper to parse file name
# -----------------------------
def parse_filename(path):
    fname = os.path.basename(path)
    m = re.search(
        r"SHAP3_SURROGATE_B_(.+?)_(CHCU_No_EC|DHCU_No_EC)_IRR_pct",
        fname
    )
    if not m:
        return {"feedstock": None, "config": None, "filename": fname}
    return {
        "feedstock": m.group(1),
        "config": m.group(2),
        "filename": fname,
    }

# -----------------------------
# Per-file summaries
# -----------------------------
per_file_feature_rows = []
case_summary_rows = []

for f in files:
    meta = parse_filename(f)
    df = pd.read_excel(f)

    shap_cols = [c for c in df.columns if c.endswith("_SHAP")]
    if not shap_cols:
        raise ValueError(f"No SHAP column found in {f}")
    shap_col = shap_cols[0]

    r2 = df["Target_R2_test"].iloc[0] if "Target_R2_test" in df.columns else np.nan

    # mean absolute SHAP per feature for this file
    agg = (
        df.groupby("Feature")[shap_col]
        .apply(lambda x: np.mean(np.abs(x)))
        .reset_index(name="mean_abs_shap")
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )

    agg["feedstock"] = meta["feedstock"]
    agg["config"] = meta["config"]
    agg["filename"] = meta["filename"]
    agg["r2_test"] = r2
    agg["rank_within_case"] = np.arange(1, len(agg) + 1)

    per_file_feature_rows.append(agg)

    # compact case summary (top 5)
    top5 = agg.head(5)
    row = {
        "filename": meta["filename"],
        "feedstock": meta["feedstock"],
        "config": meta["config"],
        "r2_test": r2,
    }
    for i, (_, r) in enumerate(top5.iterrows(), start=1):
        row[f"top{i}_feature"] = r["Feature"]
        row[f"top{i}_mean_abs_shap"] = r["mean_abs_shap"]
    case_summary_rows.append(row)

df_case_feature = pd.concat(per_file_feature_rows, ignore_index=True)
df_case_summary = pd.DataFrame(case_summary_rows)

# -----------------------------
# Cross-case feature ranking
# -----------------------------
df_cross = (
    df_case_feature
    .groupby("Feature")["mean_abs_shap"]
    .agg(
        mean_abs_shap_mean="mean",
        mean_abs_shap_median="median",
        mean_abs_shap_max="max",
        mean_abs_shap_min="min",
        n_cases="count",
    )
    .sort_values("mean_abs_shap_mean", ascending=False)
    .reset_index()
)

# -----------------------------
# Feedstock-level summary
# -----------------------------
df_by_feedstock = (
    df_case_feature
    .groupby(["feedstock", "Feature"])["mean_abs_shap"]
    .mean()
    .reset_index()
)

df_by_feedstock["rank_within_feedstock"] = (
    df_by_feedstock.groupby("feedstock")["mean_abs_shap"]
    .rank(ascending=False, method="dense")
)

# -----------------------------
# Config-level summary
# -----------------------------
df_by_config = (
    df_case_feature
    .groupby(["config", "Feature"])["mean_abs_shap"]
    .mean()
    .reset_index()
)

df_by_config["rank_within_config"] = (
    df_by_config.groupby("config")["mean_abs_shap"]
    .rank(ascending=False, method="dense")
)

# -----------------------------
# Wide table for plotting later
# -----------------------------
# Useful for driver-shift or grouped bar charts
df_plot_wide = df_case_feature.pivot_table(
    index="Feature",
    columns=["feedstock", "config"],
    values="mean_abs_shap",
    aggfunc="mean"
)

# -----------------------------
# Save outputs
# -----------------------------
excel_path = os.path.join(OUTPUT_DIR, "modelB_IRR_shap_summary.xlsx")
with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    df_case_feature.to_excel(writer, sheet_name="per_case_feature_summary", index=False)
    df_case_summary.to_excel(writer, sheet_name="per_case_top5", index=False)
    df_cross.to_excel(writer, sheet_name="cross_case_ranking", index=False)
    df_by_feedstock.to_excel(writer, sheet_name="by_feedstock", index=False)
    df_by_config.to_excel(writer, sheet_name="by_config", index=False)
    df_plot_wide.to_excel(writer, sheet_name="plot_wide_matrix")

print(f"Saved summary workbook: {excel_path}")

# -----------------------------
# Print quick console summary
# -----------------------------
print("\n=== Cross-case ranking ===")
print(df_cross.to_string(index=False))

print("\n=== Per-case top 5 ===")
print(df_case_summary.to_string(index=False))