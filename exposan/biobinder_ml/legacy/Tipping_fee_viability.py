# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 12:17:46 2026

@author: aliah
"""

import os
import glob
import re
import numpy as np
import pandas as pd

# =========================================================
# USER SETTINGS
# =========================================================
results_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results"
out_excel = os.path.join(results_dir, "tipping_fee_viability_dataset.xlsx")
out_csv   = os.path.join(results_dir, "tipping_fee_viability_summary.csv")

feedstocks = ["food", "sludge", "manure", "green"]
configs = ["CHCU", "DHCU"]
prefer = "max"

TARGET_IRR = 10.0

# Explicit removals only
EXCLUDE_FEEDSTOCKS = {"fog"}
ARTIFACT_BIOBINDER_MANURE = 97.90195447
ARTIFACT_SAF_EXTREME = 1609372.922

# =========================================================
# HELPERS
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


def first_existing(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    return None


def estimate_threshold_linear(df_group, target_irr=10.0):
    """
    Fit: IRR = a + b * tipping_fee

    Returns:
        intercept (a)
        slope (b)
        threshold tipping fee for IRR = target_irr
        R²
        RMSE
        N used
    """
    sub = df_group.dropna(subset=["IRR_pct", "Tipping_fee_$/wet_tonne"]).copy()

    if len(sub) < 5:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(sub)

    x = sub["Tipping_fee_$/wet_tonne"].to_numpy(dtype=float)
    y = sub["IRR_pct"].to_numpy(dtype=float)

    # Avoid degenerate case
    if np.allclose(np.std(x), 0):
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(sub)

    # =========================
    # Linear fit
    # =========================
    slope, intercept = np.polyfit(x, y, 1)

    y_pred = intercept + slope * x

    # =========================
    # R²
    # =========================
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    if np.isclose(ss_tot, 0):
        r2 = np.nan
    else:
        r2 = 1 - (ss_res / ss_tot)

    # =========================
    # RMSE
    # =========================
    rmse = np.sqrt(np.mean((y - y_pred) ** 2))

    # =========================
    # Threshold
    # =========================
    if np.isclose(slope, 0):
        threshold = np.nan
    else:
        threshold = (target_irr - intercept) / slope

    return intercept, slope, threshold, r2, rmse, len(sub)

# =========================================================
# COLLECT LONG DATA
# =========================================================
all_rows = []
file_log = []

for feed in feedstocks:
    if feed in EXCLUDE_FEEDSTOCKS:
        continue

    for cfg in configs:

        # -------------------------
        # Biobinder
        # -------------------------
        try:
            f = find_latest_meta(feed, cfg, prefer, results_dir)
            df = pd.read_excel(f)

            irr_col = first_existing(df.columns, ["IRR_pct"])
            price_col = first_existing(df.columns, ["Feedstock_price_$/tonne", "Feedstock_price", "feedstock_price_$/tonne"])

            if irr_col is None:
                raise ValueError(f"'IRR_pct' column not found in {os.path.basename(f)}")

            work = df.copy()
            work["IRR_pct"] = pd.to_numeric(work[irr_col], errors="coerce")

            if price_col is not None:
                work["Feedstock_price_$/tonne"] = pd.to_numeric(work[price_col], errors="coerce")
            else:
                work["Feedstock_price_$/tonne"] = np.nan

            work["Tipping_fee_$/wet_tonne"] = -work["Feedstock_price_$/tonne"]
            work["Feedstock"] = feed
            work["Config"] = cfg
            work["Product set"] = "Biobinder"
            work["Source file"] = os.path.basename(f)

            all_rows.append(
                work[["Feedstock", "Config", "Product set", "Source file",
                      "IRR_pct", "Feedstock_price_$/tonne", "Tipping_fee_$/wet_tonne"]]
            )

            file_log.append({
                "Feedstock": feed,
                "Config": cfg,
                "Product set": "Biobinder",
                "Source file": os.path.basename(f),
                "Total rows": len(df)
            })

        except Exception as e:
            print(f"[Biobinder] {feed} | {cfg} -> {e}")

        # -------------------------
        # SAF
        # -------------------------
        try:
            f = find_latest_meta_saf(feed, cfg, prefer, results_dir)
            df = pd.read_excel(f)

            irr_col = first_existing(df.columns, ["IRR_pct"])
            price_col = first_existing(df.columns, ["Feedstock_price_$/tonne", "Feedstock_price", "feedstock_price_$/tonne"])

            if irr_col is None:
                raise ValueError(f"'IRR_pct' column not found in {os.path.basename(f)}")

            work = df.copy()
            work["IRR_pct"] = pd.to_numeric(work[irr_col], errors="coerce")

            if price_col is not None:
                work["Feedstock_price_$/tonne"] = pd.to_numeric(work[price_col], errors="coerce")
            else:
                work["Feedstock_price_$/tonne"] = np.nan

            work["Tipping_fee_$/wet_tonne"] = -work["Feedstock_price_$/tonne"]
            work["Feedstock"] = feed
            work["Config"] = cfg
            work["Product set"] = "SAF"
            work["Source file"] = os.path.basename(f)

            all_rows.append(
                work[["Feedstock", "Config", "Product set", "Source file",
                      "IRR_pct", "Feedstock_price_$/tonne", "Tipping_fee_$/wet_tonne"]]
            )

            file_log.append({
                "Feedstock": feed,
                "Config": cfg,
                "Product set": "SAF",
                "Source file": os.path.basename(f),
                "Total rows": len(df)
            })

        except Exception as e:
            print(f"[SAF] {feed} | {cfg} -> {e}")

if not all_rows:
    raise ValueError("No valid META files were found.")

long_df = pd.concat(all_rows, ignore_index=True)

# =========================================================
# CLEANING
# =========================================================
# Keep valid numeric IRR only
long_df = long_df.dropna(subset=["IRR_pct"]).copy()

# Explicit artifact removal only
mask_bio_artifact = (
    (long_df["Feedstock"] == "manure") &
    (long_df["Product set"] == "Biobinder") &
    (np.isclose(long_df["IRR_pct"], ARTIFACT_BIOBINDER_MANURE, atol=1e-9))
)
long_df = long_df.loc[~mask_bio_artifact].copy()

mask_saf_artifact = (
    (long_df["Product set"] == "SAF") &
    (np.isclose(long_df["IRR_pct"], ARTIFACT_SAF_EXTREME, atol=1e-6))
)
long_df = long_df.loc[~mask_saf_artifact].copy()

# Viability flag
long_df["IRR_ge_10"] = long_df["IRR_pct"] >= TARGET_IRR

# =========================================================
# SUMMARY TABLE FOR THRESHOLD PLOTTING
# =========================================================
group_cols = ["Feedstock", "Config", "Product set", "Source file"]

summary_rows = []
for keys, grp in long_df.groupby(group_cols):
    feed, cfg, pset, src = keys

    valid_fee = grp.dropna(subset=["Tipping_fee_$/wet_tonne"])
    viable = grp[grp["IRR_ge_10"]]
    viable_fee = viable.dropna(subset=["Tipping_fee_$/wet_tonne"])

    intercept, slope, threshold_fee, r2, rmse, n_fit = estimate_threshold_linear(grp, TARGET_IRR)

    summary_rows.append({
        "Feedstock": feed,
        "Config": cfg,
        "Product set": pset,
        "Source file": src,
        "N all": len(grp),
        "N viable (IRR>=10)": len(viable),
        "Viability fraction": len(viable) / len(grp) if len(grp) > 0 else np.nan,
        "Valid fee rows": len(valid_fee),
        "Valid viable fee rows": len(viable_fee),
        "IRR min": grp["IRR_pct"].min(),
        "IRR p10": grp["IRR_pct"].quantile(0.10),
        "IRR p25": grp["IRR_pct"].quantile(0.25),
        "IRR median": grp["IRR_pct"].median(),
        "IRR p75": grp["IRR_pct"].quantile(0.75),
        "IRR p90": grp["IRR_pct"].quantile(0.90),
        "IRR max": grp["IRR_pct"].max(),
        "Tipping fee min": valid_fee["Tipping_fee_$/wet_tonne"].min() if len(valid_fee) else np.nan,
        "Tipping fee p10": valid_fee["Tipping_fee_$/wet_tonne"].quantile(0.10) if len(valid_fee) else np.nan,
        "Tipping fee p25": valid_fee["Tipping_fee_$/wet_tonne"].quantile(0.25) if len(valid_fee) else np.nan,
        "Tipping fee median": valid_fee["Tipping_fee_$/wet_tonne"].median() if len(valid_fee) else np.nan,
        "Tipping fee p75": valid_fee["Tipping_fee_$/wet_tonne"].quantile(0.75) if len(valid_fee) else np.nan,
        "Tipping fee p90": valid_fee["Tipping_fee_$/wet_tonne"].quantile(0.90) if len(valid_fee) else np.nan,
        "Tipping fee max": valid_fee["Tipping_fee_$/wet_tonne"].max() if len(valid_fee) else np.nan,
        "Linear fit intercept": intercept,
        "Linear fit slope": slope,
        "Linear fit R2": r2,
        "Linear fit RMSE": rmse,
        "N used in fit": n_fit,
        "Estimated tipping fee for IRR=10": threshold_fee,
    })

summary_df = pd.DataFrame(summary_rows)

# =========================================================
# SORTING
# =========================================================
feedstock_order = ["food", "sludge", "manure", "green"]
config_order = ["CHCU", "DHCU"]
product_order = ["Biobinder", "SAF"]

for df_ in [long_df, summary_df]:
    df_["Feedstock"] = pd.Categorical(df_["Feedstock"], categories=feedstock_order, ordered=True)
    df_["Config"] = pd.Categorical(df_["Config"], categories=config_order, ordered=True)
    df_["Product set"] = pd.Categorical(df_["Product set"], categories=product_order, ordered=True)

long_df = long_df.sort_values(["Product set", "Feedstock", "Config", "IRR_pct"]).reset_index(drop=True)
summary_df = summary_df.sort_values(["Product set", "Feedstock", "Config"]).reset_index(drop=True)

viable_df = long_df[long_df["IRR_ge_10"]].copy()

# =========================================================
# SAVE
# =========================================================
with pd.ExcelWriter(out_excel, engine="openpyxl") as writer:
    long_df.to_excel(writer, sheet_name="Long_All_Cases", index=False)
    viable_df.to_excel(writer, sheet_name="Viable_IRR_ge_10", index=False)
    summary_df.to_excel(writer, sheet_name="Threshold_Summary", index=False)
    pd.DataFrame(file_log).to_excel(writer, sheet_name="File_Log", index=False)

summary_df.to_csv(out_csv, index=False)

print("Saved:")
print(out_excel)
print(out_csv)