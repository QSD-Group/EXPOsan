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
import re
import glob
import numpy as np
import pandas as pd
from datetime import datetime


# -----------------------------
# Settings
# -----------------------------
INPUT_DIR = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results"
OUTPUT_DIR = os.path.join(INPUT_DIR, "summary_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M")


# -----------------------------
# Find SHAP3 files
# -----------------------------
pattern = os.path.join(INPUT_DIR, "SHAP3_SURROGATE_B_10000_*_IRR_pct_*.xlsx")
all_files = glob.glob(pattern)

if not all_files:
    raise FileNotFoundError(f"No files found matching: {pattern}")


# -----------------------------
# Helper to parse file name
# -----------------------------
def parse_filename(path):
    fname = os.path.basename(path)
    m = re.search(r"SHAP3_SURROGATE_B_(.+?)_(CHCU_No_EC|DHCU_No_EC)_IRR_pct", fname)

    if not m:
        return {"feedstock": None, "config": None, "filename": fname}

    return {
        "feedstock": m.group(1),
        "config": m.group(2),
        "filename": fname,
    }


# -----------------------------
# Keep latest file per feedstock × config
# -----------------------------
latest_files = {}

for f in all_files:
    meta = parse_filename(f)

    if meta["feedstock"] is None or meta["config"] is None:
        print(f"Skipping unrecognized file: {os.path.basename(f)}")
        continue

    key = (meta["feedstock"], meta["config"])

    if key not in latest_files or os.path.getmtime(f) > os.path.getmtime(latest_files[key]):
        latest_files[key] = f

files = sorted(
    latest_files.values(),
    key=lambda f: (parse_filename(f)["feedstock"], parse_filename(f)["config"]),
)

print(f"\nFound {len(all_files)} total SHAP files.")
print(f"Using {len(files)} latest feedstock × configuration files.\n")

for f in files:
    meta = parse_filename(f)
    print(f"{meta['feedstock']:15s} | {meta['config']:12s} | {meta['filename']}")

if len(files) != 8:
    print(f"\nWARNING: Expected 8 unique cases but found {len(files)}.")


# -----------------------------
# Helper to read model metric
# -----------------------------
def get_metric(df, column):
    if column not in df.columns:
        return np.nan

    values = pd.to_numeric(df[column], errors="coerce").dropna()
    return values.iloc[0] if len(values) else np.nan


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

    # Test-set model performance
    r2 = get_metric(df, "Target_R2_test")
    rmse = get_metric(df, "Target_RMSE_test")
    mae = get_metric(df, "Target_MAE_test")

    # Mean absolute SHAP per feature for this file
    df[shap_col] = pd.to_numeric(df[shap_col], errors="coerce")

    agg = (
        df.dropna(subset=["Feature", shap_col])
        .groupby("Feature")[shap_col]
        .apply(lambda x: np.mean(np.abs(x)))
        .reset_index(name="mean_abs_shap")
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )

    agg["feedstock"] = meta["feedstock"]
    agg["config"] = meta["config"]
    agg["filename"] = meta["filename"]
    agg["r2_test"] = r2
    agg["rmse_test"] = rmse
    agg["mae_test"] = mae
    agg["rank_within_case"] = np.arange(1, len(agg) + 1)

    per_file_feature_rows.append(agg)

    # Compact case summary
    top5 = agg.head(5)

    row = {
        "filename": meta["filename"],
        "feedstock": meta["feedstock"],
        "config": meta["config"],
        "r2_test": r2,
        "rmse_test": rmse,
        "mae_test": mae,
    }

    for i, (_, r) in enumerate(top5.iterrows(), start=1):
        row[f"top{i}_feature"] = r["Feature"]
        row[f"top{i}_mean_abs_shap"] = r["mean_abs_shap"]

    case_summary_rows.append(row)


df_case_feature = pd.concat(per_file_feature_rows, ignore_index=True)
df_case_summary = pd.DataFrame(case_summary_rows)


# -----------------------------
# Model performance summary
# -----------------------------
df_performance = (
    df_case_summary[
        ["feedstock", "config", "r2_test", "rmse_test", "mae_test", "filename"]
    ]
    .sort_values(["feedstock", "config"])
    .reset_index(drop=True)
)


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

# Normalize to most influential cross-case feature
max_shap = df_cross["mean_abs_shap_mean"].max()
df_cross["normalized_shap"] = df_cross["mean_abs_shap_mean"] / max_shap
df_cross["rank"] = np.arange(1, len(df_cross) + 1)

df_cross = df_cross[
    [
        "rank",
        "Feature",
        "mean_abs_shap_mean",
        "normalized_shap",
        "mean_abs_shap_median",
        "mean_abs_shap_max",
        "mean_abs_shap_min",
        "n_cases",
    ]
]


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
df_plot_wide = df_case_feature.pivot_table(
    index="Feature",
    columns=["feedstock", "config"],
    values="mean_abs_shap",
    aggfunc="mean",
)


# -----------------------------
# Save outputs
# -----------------------------
excel_path = os.path.join(OUTPUT_DIR, f"modelB_IRR_shap_summary_{timestamp}.xlsx")

with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    df_performance.to_excel(writer, sheet_name="model_performance", index=False)
    df_case_feature.to_excel(writer, sheet_name="per_case_feature_summary", index=False)
    df_case_summary.to_excel(writer, sheet_name="per_case_top5", index=False)
    df_cross.to_excel(writer, sheet_name="cross_case_ranking", index=False)
    df_by_feedstock.to_excel(writer, sheet_name="by_feedstock", index=False)
    df_by_config.to_excel(writer, sheet_name="by_config", index=False)
    df_plot_wide.to_excel(writer, sheet_name="plot_wide_matrix")

print(f"\nSaved summary workbook: {excel_path}")


# -----------------------------
# Print quick console summaries
# -----------------------------
print("\n=== Model performance ===")
print(df_performance[["feedstock", "config", "r2_test", "rmse_test", "mae_test"]].to_string(index=False))

print("\n=== Cross-case ranking ===")
print(df_cross[["rank", "Feature", "mean_abs_shap_mean", "normalized_shap", "n_cases"]].to_string(index=False))

print("\n=== Per-case top 5 ===")
print(df_case_summary.to_string(index=False))


# =============================================================================
# SHAP DIRECTION ANALYSIS
# =============================================================================

all_rows = []

for f in files:
    meta = parse_filename(f)
    df = pd.read_excel(f)

    shap_cols = [c for c in df.columns if c.endswith("_SHAP")]

    if not shap_cols:
        raise ValueError(f"No SHAP column found in {f}")

    shap_col = shap_cols[0]

    r2 = get_metric(df, "Target_R2_test")
    rmse = get_metric(df, "Target_RMSE_test")
    mae = get_metric(df, "Target_MAE_test")

    # One row per feature summary within each case
    for feat, g in df.groupby("Feature"):
        x = pd.to_numeric(g["Feature value"], errors="coerce")
        s = pd.to_numeric(g[shap_col], errors="coerce")

        valid = x.notna() & s.notna()
        x_valid = x[valid]
        s_valid = s[valid]

        corr = np.nan

        if len(x_valid) > 2 and x_valid.nunique() > 1 and s_valid.nunique() > 1:
            corr = x_valid.corr(s_valid)

        s_nonmissing = s.dropna()

        all_rows.append({
            "filename": meta["filename"],
            "feedstock": meta["feedstock"],
            "config": meta["config"],
            "Feature": feat,
            "r2_test": r2,
            "rmse_test": rmse,
            "mae_test": mae,
            "mean_abs_shap": np.mean(np.abs(s_nonmissing)) if len(s_nonmissing) else np.nan,
            "mean_signed_shap": np.mean(s_nonmissing) if len(s_nonmissing) else np.nan,
            "median_signed_shap": np.median(s_nonmissing) if len(s_nonmissing) else np.nan,
            "positive_share": np.mean(s_nonmissing > 0) if len(s_nonmissing) else np.nan,
            "negative_share": np.mean(s_nonmissing < 0) if len(s_nonmissing) else np.nan,
            "feature_shap_corr": corr,
        })


df_case = pd.DataFrame(all_rows)


# -----------------------------
# Cross-case direction aggregation
# -----------------------------
df_cross_direction = (
    df_case
    .groupby("Feature")
    .agg(
        mean_abs_shap=("mean_abs_shap", "mean"),
        mean_signed_shap=("mean_signed_shap", "mean"),
        median_signed_shap=("median_signed_shap", "mean"),
        positive_share=("positive_share", "mean"),
        negative_share=("negative_share", "mean"),
        feature_shap_corr=("feature_shap_corr", "mean"),
        n_cases=("Feature", "count"),
    )
    .reset_index()
)


# -----------------------------
# Direction label
# -----------------------------
def direction_label(row):
    corr = row["feature_shap_corr"]
    pos = row["positive_share"]
    ms = row["mean_signed_shap"]

    if pd.notna(corr):
        if corr > 0.2:
            return "higher values tend to increase IRR"
        elif corr < -0.2:
            return "higher values tend to decrease IRR"

    if pd.notna(ms) and pd.notna(pos):
        if ms > 0 and pos > 0.6:
            return "mostly positive"
        elif ms < 0 and pos < 0.4:
            return "mostly negative"

    return "mixed / context-dependent"


df_cross_direction["direction_summary"] = df_cross_direction.apply(direction_label, axis=1)

# Rank by importance
df_cross_direction = df_cross_direction.sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
df_cross_direction["rank"] = np.arange(1, len(df_cross_direction) + 1)

# Normalized SHAP
max_direction_shap = df_cross_direction["mean_abs_shap"].max()
df_cross_direction["normalized_shap"] = df_cross_direction["mean_abs_shap"] / max_direction_shap

df_cross_direction = df_cross_direction[
    [
        "rank",
        "Feature",
        "mean_abs_shap",
        "normalized_shap",
        "mean_signed_shap",
        "median_signed_shap",
        "feature_shap_corr",
        "positive_share",
        "negative_share",
        "direction_summary",
        "n_cases",
    ]
]


# -----------------------------
# Save direction outputs
# -----------------------------
out_path = os.path.join(OUTPUT_DIR, f"modelB_IRR_shap_direction_summary_{timestamp}.xlsx")

with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    df_case.to_excel(writer, sheet_name="per_case_direction", index=False)
    df_cross_direction.to_excel(writer, sheet_name="cross_case_direction", index=False)

print(f"\nSaved direction summary: {out_path}")


# -----------------------------
# Print direction summary
# -----------------------------
print("\n=== Cross-case SHAP direction ===")

print(
    df_cross_direction[
        [
            "rank",
            "Feature",
            "mean_abs_shap",
            "normalized_shap",
            "mean_signed_shap",
            "feature_shap_corr",
            "positive_share",
            "direction_summary",
        ]
    ].to_string(index=False)
)