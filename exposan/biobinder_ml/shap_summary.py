# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 14:35:52 2026

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
example_yield_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\SHAP_RF_sludge_20260112_1206.xlsx"
example_tea_path   = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\SURROGATE_SHAP_sludge_MSP_20260109_1550.xlsx"

results_dir = os.path.dirname(example_yield_path)

# Choose which yield target you want from the RF yield SHAP files:
# options typically: "biocrude", "aqueous", "gas", "char"
YIELD_TARGET = "biocrude"   # <-- change if desired

# Knobs you vary (fair comparison set)
KNOB_FEATURES = [
    "Temperature (C)",
    "Residence Time",
    "Solid content (w/w) %",
    "Catalyst",
    "Reactor Type",
    "Pre-processing",
    "Solvent",
]
KNOB_SET = set(KNOB_FEATURES)

# Output
out_path = os.path.join(results_dir, f"COMBINED_YIELD_TEA_SHAP_{YIELD_TARGET}_WITH_KNOB_NORMALIZATION.xlsx")


# ============================================================
# HELPERS
# ============================================================
def _pick_latest(files):
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

def _safe_to_numeric(s):
    return pd.to_numeric(s, errors="coerce")

def _percent_share(absvals: pd.Series) -> pd.Series:
    tot = float(absvals.sum())
    if not np.isfinite(tot) or tot <= 0:
        return absvals * 0.0
    return (absvals / tot) * 100.0

def _yield_shap_from_wide(df_yield: pd.DataFrame, target: str) -> pd.Series:
    """
    Yield SHAP files are wide: one row per Feature, columns like Biocrude_SHAP, Aqueous_SHAP, ...
    Return Series indexed by Feature with abs(SHAP) for the selected target.
    """
    if "Feature" not in df_yield.columns:
        raise ValueError(f"Yield file missing 'Feature'. Columns: {list(df_yield.columns)}")

    # Map target -> expected column name
    t = str(target).strip().lower()
    col_map = {
        "biocrude": "Biocrude_SHAP",
        "aqueous": "Aqueous_SHAP",
        "gas": "Gas_SHAP",
        "char": "Char_SHAP",
    }
    if t not in col_map:
        raise ValueError(f"Unknown YIELD_TARGET='{target}'. Use one of {list(col_map.keys())}")

    shap_col = col_map[t]
    if shap_col not in df_yield.columns:
        raise ValueError(
            f"Yield file does not contain '{shap_col}'. Columns: {list(df_yield.columns)}"
        )

    d = df_yield.copy()
    d[shap_col] = _safe_to_numeric(d[shap_col])
    d = d[np.isfinite(d[shap_col])].copy()

    s = d.set_index("Feature")[shap_col].apply(lambda x: float(abs(x)))
    return s

def _msp_shap_from_long(df_tea: pd.DataFrame) -> pd.Series:
    """
    TEA SHAP files are long: many rows per Feature with MSP_SHAP.
    Return Series indexed by Feature with mean(abs(MSP_SHAP)).
    """
    if "Feature" not in df_tea.columns:
        raise ValueError(f"TEA file missing 'Feature'. Columns: {list(df_tea.columns)}")

    # Find MSP shap column robustly
    shap_col = None
    for c in df_tea.columns:
        if str(c).strip().lower() == "msp_shap":
            shap_col = c
            break
    if shap_col is None:
        # fallback: any *_SHAP
        candidates = [c for c in df_tea.columns if str(c).strip().lower().endswith("_shap")]
        if not candidates:
            raise ValueError(f"TEA file missing MSP_SHAP / *_SHAP. Columns: {list(df_tea.columns)}")
        shap_col = sorted(candidates)[0]

    d = df_tea.copy()
    d[shap_col] = _safe_to_numeric(d[shap_col])
    d = d[np.isfinite(d[shap_col])].copy()

    s = d.groupby("Feature")[shap_col].apply(lambda x: float(np.mean(np.abs(x.values))))
    s = s.sort_values(ascending=False)
    return s

def _knob_norm(abs_series: pd.Series, knob_set: set) -> pd.Series:
    """
    Normalize abs SHAP *within knobs only*.
    Returns Series for all features; non-knobs become NaN (so they don't confuse the view).
    """
    knob_sum = float(abs_series[abs_series.index.isin(knob_set)].sum())
    out = pd.Series(index=abs_series.index, dtype=float)
    if not np.isfinite(knob_sum) or knob_sum <= 0:
        out[:] = np.nan
        return out
    out.loc[abs_series.index.isin(knob_set)] = (abs_series.loc[abs_series.index.isin(knob_set)] / knob_sum) * 100.0
    out.loc[~abs_series.index.isin(knob_set)] = np.nan
    return out

def _knob_share(abs_series: pd.Series, knob_set: set) -> float:
    tot = float(abs_series.sum())
    knob = float(abs_series[abs_series.index.isin(knob_set)].sum())
    if not np.isfinite(tot) or tot <= 0:
        return float("nan")
    return (knob / tot) * 100.0


# ============================================================
# DISCOVER FILES (LATEST PER FEEDSTOCK)
# ============================================================
yield_files = glob.glob(os.path.join(results_dir, "SHAP_RF_*.xlsx"))
tea_files   = glob.glob(os.path.join(results_dir, "SURROGATE_SHAP_*_MSP_*.xlsx"))

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
# BUILD OUTPUT
# ============================================================
with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    all_rows = []
    summary_rows = []

    for feed in feedstocks:
        y_path = yield_by_feed[feed]
        t_path = tea_by_feed[feed]

        df_y = pd.read_excel(y_path)
        df_t = pd.read_excel(t_path)

        # --- Yield abs SHAP (selected target)
        y_abs = _yield_shap_from_wide(df_y, YIELD_TARGET)
        y_pct_all = _percent_share(y_abs).rename("Yield_SHAP_%")
        y_pct_knob = _knob_norm(y_abs, KNOB_SET).rename("Yield_KnobNorm_%")
        y_knob_share = _knob_share(y_abs, KNOB_SET)  # % of yield SHAP mass explained by knobs

        # --- MSP abs SHAP
        m_abs = _msp_shap_from_long(df_t)
        m_pct_all = _percent_share(m_abs).rename("MSP_SHAP_%")
        m_pct_knob = _knob_norm(m_abs, KNOB_SET).rename("MSP_KnobNorm_%")
        m_knob_share = _knob_share(m_abs, KNOB_SET)  # THIS is your "MSP knob share"

        # --- Combine
        df_out = pd.concat([y_pct_all, m_pct_all, y_pct_knob, m_pct_knob], axis=1)

        df_out["Ratio_All"] = df_out["MSP_SHAP_%"] / df_out["Yield_SHAP_%"].replace(0.0, np.nan)
        df_out["Ratio_KnobNorm"] = df_out["MSP_KnobNorm_%"] / df_out["Yield_KnobNorm_%"].replace(0.0, np.nan)

        df_out = df_out.reset_index().rename(columns={"index": "Feature"})
        df_out["IsKnob"] = df_out["Feature"].isin(KNOB_SET)

        # attach knob share scalars (repeated column so it’s visible anywhere)
        df_out["Yield_KnobShare_%"] = float(y_knob_share)
        df_out["MSP_KnobShare_%"] = float(m_knob_share)

        # reorder
        df_out = df_out[
            ["Feature", "IsKnob",
             "Yield_SHAP_%", "MSP_SHAP_%", "Ratio_All",
             "Yield_KnobNorm_%", "MSP_KnobNorm_%", "Ratio_KnobNorm",
             "Yield_KnobShare_%", "MSP_KnobShare_%"]
        ]

        # sort: emphasize MSP importance
        df_out = df_out.sort_values(["MSP_SHAP_%"], ascending=False)

        sheet = feed[:31]
        df_out.to_excel(writer, sheet_name=sheet, index=False)

        tmp = df_out.copy()
        tmp.insert(0, "Feedstock", feed)
        all_rows.append(tmp)

        summary_rows.append({
            "Feedstock": feed,
            "Yield_file": os.path.basename(y_path),
            "TEA_file": os.path.basename(t_path),
            "YIELD_TARGET": YIELD_TARGET,
            "Yield_KnobShare_%": float(y_knob_share),
            "MSP_KnobShare_%": float(m_knob_share),
            "Top_MSP_Feature": df_out.iloc[0]["Feature"] if len(df_out) else "",
            "Top_MSP_SHAP_%": float(df_out.iloc[0]["MSP_SHAP_%"]) if len(df_out) else np.nan,
        })

        print(f"[OK] {feed}: {os.path.basename(y_path)} | {os.path.basename(t_path)} | nfeat={len(df_out)}")

    # combined long-form
    df_all = pd.concat(all_rows, ignore_index=True)
    df_all.to_excel(writer, sheet_name="ALL_FEEDSTOCKS", index=False)

    # summary
    df_sum = pd.DataFrame(summary_rows)
    df_sum.to_excel(writer, sheet_name="SUMMARY", index=False)

print("✅ Done. Saved:")
print(out_path)
