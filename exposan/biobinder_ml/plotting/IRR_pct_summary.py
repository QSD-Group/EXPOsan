# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 12:22:56 2026

@author: aliah
"""

import os
import glob
import numpy as np
import pandas as pd

# =========================================================
# USER SETTINGS
# =========================================================
results_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results"

out_excel = os.path.join(results_dir, "master_IRR_10000_GWP_summary_cleaned_final.xlsx")
out_csv   = os.path.join(results_dir, "master_IRR_10000_GWP_summary_cleaned_final.csv")

# Explicit artifact removals only
ARTIFACT_BIOBINDER_MANURE = 97.90195447
ARTIFACT_SAF_EXTREME = 1609372.922

# Explicit feedstock removal only
EXCLUDE_FEEDSTOCKS = {"fog"}

feedstocks = ["food", "sludge", "manure", "green"]
configs = ["CHCU", "DHCU"]
prefer = "max"   # matches prefer-max_top_ratio...

# Change this if your column is named differently
GWP_CANDIDATES = ["GWP", "GWP_total", "GWP_kgCO2eq", "GWP100"]

# =========================================================
# FIND LATEST META FILES
# =========================================================
def find_latest_meta(feedstock_id, config_name, prefer, output_dir="results"):
    patterns = [
        # most general (handles No_EC + extra tags)
        os.path.join(
            output_dir,
            f"META_SHAP3_10000_{feedstock_id}_{config_name}_*_prefer-{prefer}*.xlsx"
        ),
        # fallback (older naming)
        os.path.join(
            output_dir,
            f"META_SHAP3_10000_{feedstock_id}_{config_name}_prefer-{prefer}*.xlsx"
        ),
    ]

    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    files = sorted(set(files))

    if not files:
        raise FileNotFoundError(
            f"No Biobinder META files found for {feedstock_id} | {config_name}\n"
            f"Patterns tried:\n" + "\n".join(patterns)
        )

    # prioritize MERGED if exists
    merged = [f for f in files if "MERGED" in os.path.basename(f)]
    if merged:
        return sorted(merged)[-1]

    return files[-1]

def find_latest_meta_saf(feedstock_id, config_name, prefer, output_dir="results"):
    patterns = [
        os.path.join(
            output_dir,
            f"META_SHAP3_SAF_{feedstock_id}_{config_name}_*_prefer-{prefer}*.xlsx"
        ),
        os.path.join(
            output_dir,
            f"META_SHAP3_SAF_{feedstock_id}_{config_name}_prefer-{prefer}*.xlsx"
        ),
    ]

    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    files = sorted(set(files))

    if not files:
        raise FileNotFoundError(
            f"No SAF META files found for {feedstock_id} | {config_name}\n"
            f"Patterns tried:\n" + "\n".join(patterns)
        )

    return files[-1]

def find_gwp_column(columns):
    for col in GWP_CANDIDATES:
        if col in columns:
            return col
    return None

# =========================================================
# STORAGE
# =========================================================
long_rows = []
summary_rows = []

# =========================================================
# LOOP OVER ALL COMBINATIONS
# =========================================================
for feed in feedstocks:

    if feed in EXCLUDE_FEEDSTOCKS:
        continue

    for cfg in configs:

        # -----------------------------
        # BIOBINDER
        # -----------------------------
        try:
            file_path = find_latest_meta(feed, cfg, prefer, results_dir)
            df = pd.read_excel(file_path)

            total_rows = len(df)
            irr = pd.to_numeric(df["IRR_pct"], errors="coerce")
            gwp_col = find_gwp_column(df.columns)

            df_valid = df.loc[irr.notna()].copy()
            df_valid["IRR_pct"] = pd.to_numeric(df_valid["IRR_pct"], errors="coerce")

            if gwp_col is not None:
                df_valid["GWP"] = pd.to_numeric(df_valid[gwp_col], errors="coerce")
            else:
                df_valid["GWP"] = np.nan

            df_valid["Feedstock"] = feed
            df_valid["Config"] = cfg
            df_valid["Product set"] = "Biobinder"
            df_valid["Source file"] = os.path.basename(file_path)

            long_rows.append(
                df_valid[["Feedstock", "Config", "Product set", "Source file", "IRR_pct", "GWP"]]
            )

            summary_rows.append({
                "Feedstock": feed,
                "Config": cfg,
                "Product set": "Biobinder",
                "Source file": os.path.basename(file_path),
                "Total rows": total_rows,
                "IRR_reason = ok": (df["IRR_reason"] == "ok").sum() if "IRR_reason" in df.columns else np.nan
            })

        except FileNotFoundError as e:
            print(e)

        # -----------------------------
        # SAF
        # -----------------------------
        try:
            file_path = find_latest_meta_saf(feed, cfg, prefer, results_dir)
            df = pd.read_excel(file_path)

            total_rows = len(df)
            irr = pd.to_numeric(df["IRR_pct"], errors="coerce")
            gwp_col = find_gwp_column(df.columns)

            df_valid = df.loc[irr.notna()].copy()
            df_valid["IRR_pct"] = pd.to_numeric(df_valid["IRR_pct"], errors="coerce")

            if gwp_col is not None:
                df_valid["GWP"] = pd.to_numeric(df_valid[gwp_col], errors="coerce")
            else:
                df_valid["GWP"] = np.nan

            df_valid["Feedstock"] = feed
            df_valid["Config"] = cfg
            df_valid["Product set"] = "SAF"
            df_valid["Source file"] = os.path.basename(file_path)

            long_rows.append(
                df_valid[["Feedstock", "Config", "Product set", "Source file", "IRR_pct", "GWP"]]
            )

            summary_rows.append({
                "Feedstock": feed,
                "Config": cfg,
                "Product set": "SAF",
                "Source file": os.path.basename(file_path),
                "Total rows": total_rows,
                "IRR_reason = ok": (df["IRR_reason"] == "ok").sum() if "IRR_reason" in df.columns else np.nan
            })

        except FileNotFoundError as e:
            print(e)

# =========================================================
# COMBINE LONG DATA
# =========================================================
if not long_rows:
    raise ValueError("No valid META files were found.")

long_df = pd.concat(long_rows, ignore_index=True)

# =========================================================
# CLEANING
# =========================================================
# Remove only explicitly requested artifacts
mask_bio_manure_artifact = (
    (long_df["Product set"] == "Biobinder") &
    (long_df["Feedstock"] == "manure") &
    (np.isclose(long_df["IRR_pct"], ARTIFACT_BIOBINDER_MANURE, atol=1e-9))
)
long_df = long_df.loc[~mask_bio_manure_artifact].copy()

mask_saf_artifact = (
    (long_df["Product set"] == "SAF") &
    (np.isclose(long_df["IRR_pct"], ARTIFACT_SAF_EXTREME, atol=1e-6))
)
long_df = long_df.loc[~mask_saf_artifact].copy()

# =========================================================
# SUMMARY STATS
# =========================================================
group_cols = ["Feedstock", "Config", "Product set", "Source file"]

# ---- IRR summary ----
irr_summary = (
    long_df.groupby(group_cols)["IRR_pct"]
    .agg(
        **{
            "Valid IRR rows": "count",
            "IRR min": "min",
            "IRR median": "median",
            "IRR max": "max",
        }
    )
    .reset_index()
)

irr_pct = (
    long_df.groupby(group_cols)["IRR_pct"]
    .quantile([0.10, 0.25, 0.75, 0.90])
    .unstack()
    .reset_index()
)
irr_pct.columns = group_cols + ["IRR p10", "IRR p25", "IRR p75", "IRR p90"]

# ---- GWP summary ----
gwp_valid = long_df.dropna(subset=["GWP"]).copy()

gwp_summary = (
    gwp_valid.groupby(group_cols)["GWP"]
    .agg(
        **{
            "Valid GWP rows": "count",
            "GWP min": "min",
            "GWP median": "median",
            "GWP max": "max",
        }
    )
    .reset_index()
)

gwp_pct = (
    gwp_valid.groupby(group_cols)["GWP"]
    .quantile([0.10, 0.25, 0.75, 0.90])
    .unstack()
    .reset_index()
)
gwp_pct.columns = group_cols + ["GWP p10", "GWP p25", "GWP p75", "GWP p90"]

# ---- Merge all ----
summary = irr_summary.merge(irr_pct, on=group_cols, how="left")
summary = summary.merge(gwp_summary, on=group_cols, how="left")
summary = summary.merge(gwp_pct, on=group_cols, how="left")

meta_df = pd.DataFrame(summary_rows).drop_duplicates(subset=group_cols)
summary = summary.merge(meta_df, on=group_cols, how="left")

summary = summary[
    [
        "Feedstock", "Config", "Product set", "Source file",
        "Total rows", "IRR_reason = ok",
        "Valid IRR rows", "IRR min", "IRR p10", "IRR p25", "IRR median", "IRR p75", "IRR p90", "IRR max",
        "Valid GWP rows", "GWP min", "GWP p10", "GWP p25", "GWP median", "GWP p75", "GWP p90", "GWP max",
    ]
]

# =========================================================
# SORTING
# =========================================================
feedstock_order = ["food", "sludge", "manure", "green"]
config_order = ["CHCU", "DHCU"]
product_order = ["Biobinder", "SAF"]

summary["Feedstock"] = pd.Categorical(summary["Feedstock"], categories=feedstock_order, ordered=True)
summary["Config"] = pd.Categorical(summary["Config"], categories=config_order, ordered=True)
summary["Product set"] = pd.Categorical(summary["Product set"], categories=product_order, ordered=True)

long_df["Feedstock"] = pd.Categorical(long_df["Feedstock"], categories=feedstock_order, ordered=True)
long_df["Config"] = pd.Categorical(long_df["Config"], categories=config_order, ordered=True)
long_df["Product set"] = pd.Categorical(long_df["Product set"], categories=product_order, ordered=True)

summary = summary.sort_values(["Product set", "Feedstock", "Config"]).reset_index(drop=True)
long_df = long_df.sort_values(["Product set", "Feedstock", "Config", "IRR_pct"]).reset_index(drop=True)

# =========================================================
# SAVE
# =========================================================
with pd.ExcelWriter(out_excel, engine="openpyxl") as writer:
    summary.to_excel(writer, sheet_name="Master_Summary", index=False)
    long_df.to_excel(writer, sheet_name="Long_All_Valid_IRR_GWP", index=False)

summary.to_csv(out_csv, index=False)

print("Saved:")
print(out_excel)
print(out_csv)