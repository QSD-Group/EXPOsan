# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 12:29:40 2026

@author: aliah
"""

import os
import re
import glob
import numpy as np
import pandas as pd

# ============================================================
# CONFIG
# ============================================================
# Example paths you gave (used ONLY to infer the directory)
example_yield_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\SHAP_RF_sludge_20260112_1206.xlsx"
example_tea_path   = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\SURROGATE_SHAP_sludge_MSP_20260109_1550.xlsx"

results_dir = os.path.dirname(example_yield_path)  # assume all files live here

# If your yield SHAP file has multiple targets (biocrude/aqueous/gas/char),
# set this to one of them to filter. Otherwise leave as None to aggregate all.
YIELD_TARGET_FILTER = None  # e.g., "biocrude" or "gas" or None

# Output file
out_path = os.path.join(results_dir, "COMBINED_YIELD_TEA_SHAP_PCT.xlsx")


# ============================================================
# HELPERS
# ============================================================
def _pick_latest(files):
    """Pick the latest file by modified time."""
    if not files:
        return None
    return max(files, key=lambda p: os.path.getmtime(p))

def _extract_feedstock_from_yield_filename(path):
    # SHAP_RF_<feedstock>_YYYYMMDD_HHMM.xlsx
    base = os.path.basename(path)
    m = re.match(r"SHAP_RF_(.+?)_\d{8}_\d{4}\.xlsx$", base, flags=re.IGNORECASE)
    return m.group(1) if m else None

def _extract_feedstock_from_tea_filename(path):
    # SURROGATE_SHAP_<feedstock>_MSP_YYYYMMDD_HHMM.xlsx
    base = os.path.basename(path)
    m = re.match(r"SURROGATE_SHAP_(.+?)_MSP_\d{8}_\d{4}\.xlsx$", base, flags=re.IGNORECASE)
    return m.group(1) if m else None

def _find_shap_column(df, preferred=None):
    """
    Find a SHAP column robustly.
    - If preferred is present, use it.
    - Else: choose the first column that ends with '_SHAP' or equals 'SHAP' (case-insensitive),
      preferring MSP_SHAP if present.
    """
    cols = list(df.columns)

    if preferred and preferred in cols:
        return preferred

    # Prefer exact MSP_SHAP if exists
    for c in cols:
        if str(c).strip().lower() == "msp_shap":
            return c

    # Prefer any *_SHAP columns
    shap_like = []
    for c in cols:
        cl = str(c).strip().lower()
        if cl == "shap" or cl.endswith("_shap"):
            shap_like.append(c)

    if not shap_like:
        raise ValueError(
            "Could not find a SHAP column. Expected 'MSP_SHAP', 'SHAP', or something ending with '_SHAP'. "
            f"Columns found: {cols}"
        )

    # If multiple, pick the first deterministically
    shap_like_sorted = sorted(shap_like, key=lambda x: str(x))
    return shap_like_sorted[0]

def _mean_abs_shap_by_feature(df, feature_col="Feature", shap_col=None, target_filter=None, target_col="Target"):
    """
    Return Series indexed by Feature: mean(abs(SHAP)).
    Optionally filter by Target if target_filter is not None and Target column exists.
    """
    if feature_col not in df.columns:
        raise ValueError(f"Missing '{feature_col}' column. Columns found: {list(df.columns)}")

    if shap_col is None:
        shap_col = _find_shap_column(df)

    d = df.copy()

    # Optional target filter
    if target_filter is not None and target_col in d.columns:
        d = d[d[target_col].astype(str).str.lower() == str(target_filter).lower()].copy()

    # Clean SHAP values
    d = d[np.isfinite(pd.to_numeric(d[shap_col], errors="coerce"))].copy()
    d[shap_col] = pd.to_numeric(d[shap_col], errors="coerce")

    # Mean(|SHAP|)
    out = d.groupby(feature_col)[shap_col].apply(lambda x: float(np.mean(np.abs(x.values))))
    out = out.sort_values(ascending=False)
    return out

def _to_percent_share(mean_abs_series):
    """Convert mean(|SHAP|) Series to percent contributions summing to 100."""
    total = float(mean_abs_series.sum())
    if total <= 0 or not np.isfinite(total):
        # return zeros to avoid blowups
        return mean_abs_series * 0.0
    return (mean_abs_series / total) * 100.0


# ============================================================
# DISCOVER FILES
# ============================================================
yield_files = glob.glob(os.path.join(results_dir, "SHAP_RF_*.xlsx"))
tea_files   = glob.glob(os.path.join(results_dir, "SURROGATE_SHAP_*_MSP_*.xlsx"))

# Map feedstock -> latest file
yield_by_feed = {}
for p in yield_files:
    f = _extract_feedstock_from_yield_filename(p)
    if f:
        yield_by_feed.setdefault(f, []).append(p)
yield_by_feed = {f: _pick_latest(paths) for f, paths in yield_by_feed.items()}

tea_by_feed = {}
for p in tea_files:
    f = _extract_feedstock_from_tea_filename(p)
    if f:
        tea_by_feed.setdefault(f, []).append(p)
tea_by_feed = {f: _pick_latest(paths) for f, paths in tea_by_feed.items()}

# Intersection: feedstocks that have BOTH yield + TEA files
feedstocks = sorted(set(yield_by_feed).intersection(set(tea_by_feed)))
if not feedstocks:
    raise RuntimeError(
        f"No matching feedstocks found with BOTH SHAP_RF_*.xlsx and SURROGATE_SHAP_*_MSP_*.xlsx in:\n{results_dir}\n"
        f"Found yield feedstocks: {sorted(yield_by_feed.keys())}\n"
        f"Found TEA feedstocks:   {sorted(tea_by_feed.keys())}"
    )

print("Feedstocks found:", feedstocks)
print("Writing:", out_path)


# ============================================================
# BUILD TABLES + SAVE
# ============================================================
with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    combined_rows = []

    for feed in feedstocks:
        y_path = yield_by_feed[feed]
        t_path = tea_by_feed[feed]

        df_y = pd.read_excel(y_path)
        df_t = pd.read_excel(t_path)

        # Yield: mean abs SHAP by feature -> percent
        y_shap_col = _find_shap_column(df_y)  # yield files may use 'SHAP' or something *_SHAP
        y_mean_abs = _mean_abs_shap_by_feature(
            df_y, feature_col="Feature", shap_col=y_shap_col,
            target_filter=YIELD_TARGET_FILTER, target_col="Target"
        )
        y_pct = _to_percent_share(y_mean_abs).rename("Yield_SHAP_%")

        # TEA: MSP mean abs SHAP by feature -> percent
        t_shap_col = _find_shap_column(df_t, preferred="MSP_SHAP")
        t_mean_abs = _mean_abs_shap_by_feature(df_t, feature_col="Feature", shap_col=t_shap_col)
        t_pct = _to_percent_share(t_mean_abs).rename("MSP_SHAP_%")

        # Align features (outer join)
        df_out = pd.concat([y_pct, t_pct], axis=1)

        # Ratio
        # (avoid division by 0; if Yield_SHAP_% == 0 => NaN)
        df_out["Ratio"] = df_out["MSP_SHAP_%"] / df_out["Yield_SHAP_%"].replace(0.0, np.nan)

        # Final formatting/order
        df_out = df_out.reset_index().rename(columns={"index": "Feature"})
        df_out = df_out[["Feature", "Yield_SHAP_%", "MSP_SHAP_%", "Ratio"]]

        # Sort: most economically important first (or change to Yield_SHAP_% if you prefer)
        df_out = df_out.sort_values("MSP_SHAP_%", ascending=False)

        # Write one sheet per feedstock
        sheet = feed[:31]  # Excel sheet name limit
        df_out.to_excel(writer, sheet_name=sheet, index=False)

        # Build combined sheet rows
        tmp = df_out.copy()
        tmp.insert(0, "Feedstock", feed)
        combined_rows.append(tmp)

        print(f"[OK] {feed}: Yield file='{os.path.basename(y_path)}' | TEA file='{os.path.basename(t_path)}' "
              f"| features={len(df_out)}")

    # Combined long-form sheet
    df_all = pd.concat(combined_rows, ignore_index=True)
    df_all.to_excel(writer, sheet_name="ALL_FEEDSTOCKS", index=False)

print("✅ Done. Output saved to:")
print(out_path)
