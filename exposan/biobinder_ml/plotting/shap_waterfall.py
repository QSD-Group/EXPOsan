# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 15:27:18 2026

@author: aliah
"""

"""
Analyze SURROGATE_SHAP_MSP_20260109_1440.xlsx

What this script does
---------------------
1) Loads your long-form SHAP table (one row per sample-feature).
2) Computes:
   - Global feature importance = mean(|SHAP|)
   - Directionality = corr(feature_value, SHAP) per feature
   - Category importance (composition vs operating vs catalyst/process)
3) Builds per-sample SHAP wide table (samples × features) and:
   - Picks a representative sample (largest total |SHAP|) unless you specify one
   - Saves a "waterfall-like" plot of top contributions for that sample
4) Writes summary tables to an output Excel and saves plots as PNGs.

Notes
-----
- Your file produced columns like:
    "Feature", "Feature value", "MSP_SHAP", "MSP_value", "Sample_ID", ...
- A true SHAP waterfall needs the model expected value (base value). Your Excel
  doesn’t include it, so the plot here is a contribution waterfall centered at 0
  (still excellent for *sensitivity* interpretation: which features push MSP up/down
  and by how much for a given sample).

Run
---
python analyze_surrogate_shap_msp.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# User inputs
# -----------------------------
INPUT_XLSX = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\SURROGATE_SHAP_MSP_20260109_1440.xlsx"

# Optional: set to an int Sample_ID you want a waterfall for; otherwise auto-pick
SAMPLE_ID_FOR_WATERFALL = None  # e.g., 123

# How many top features to show on the sample waterfall plot
TOP_N_WATERFALL = 15

# Output folder (defaults to same folder as input)
OUT_DIR = os.path.dirname(INPUT_XLSX)
OUT_TAG = os.path.splitext(os.path.basename(INPUT_XLSX))[0].replace("SURROGATE_SHAP_", "ANALYSIS_")


# -----------------------------
# Category map (edit if you want)
# -----------------------------
CATEGORY_MAP = {
    # Composition / feedstock
    "Carbohydrates wt%": "Composition",
    "Protein wt%": "Composition",
    "Lipids wt%": "Composition",
    "Ash wt%": "Composition",
    "C%": "Composition",
    "H%": "Composition",
    "O%": "Composition",
    "N%": "Composition",
    "HHV Biomass": "Composition",

    # Operating conditions / scale
    "Temperature (C)": "Operating conditions",
    "Residence Time": "Operating conditions",
    "Solid content (w/w) %": "Operating conditions",
    "Reactor Volume (mL)": "Operating conditions",

    # Categorical / process choice
    "Catalyst": "Process choices",
    "Solvent": "Process choices",
    "Pre-processing": "Process choices",
    "Reactor Type": "Process choices",
}


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Return Pearson corr, safely handling constants/NaNs."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < 3:
        return np.nan
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def ensure_columns(df: pd.DataFrame, required: list[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns in input Excel:\n"
            f"{missing}\n\n"
            f"Available columns:\n{list(df.columns)}"
        )


def sample_waterfall_plot(sample_shap: pd.Series, out_png: str, title: str, top_n: int = 15):
    """
    Create a contribution waterfall-like chart:
    sorted SHAP contributions, cumulative sum, centered at 0.

    This is not the classic SHAP waterfall (no base value), but it correctly shows
    which features push MSP up/down for that sample and by how much.
    """
    s = sample_shap.dropna().copy()
    # sort by abs magnitude
    s = s.reindex(s.abs().sort_values(ascending=False).index)
    s_top = s.iloc[:top_n]
    others = s.iloc[top_n:].sum() if len(s) > top_n else 0.0

    # build display series
    if len(s) > top_n:
        s_plot = pd.concat([s_top, pd.Series({"(all other features)": others})])
    else:
        s_plot = s_top

    # order by signed contribution (keep abs-sorted for interpretation)
    contrib = s_plot.values
    labels = s_plot.index.tolist()

    # cumulative positions for "waterfall"
    cum = np.cumsum(contrib)
    start = np.concatenate([[0.0], cum[:-1]])
    end = cum

    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()

    # bars: green for negative (reduces MSP), red for positive (increases MSP)
    colors = ["tab:red" if v > 0 else "tab:green" for v in contrib]
    ax.barh(range(len(contrib)), contrib, left=start, color=colors)

    ax.axvline(0, linewidth=1)
    ax.set_yticks(range(len(contrib)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("SHAP contribution to MSP (surrogate units)")
    ax.set_title(title)

    # annotate end points
    for i, (st, ed, v) in enumerate(zip(start, end, contrib)):
        ax.text(ed, i, f" {v:+.3g}", va="center", fontsize=9)

    plt.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    # -----------------------------
    # Load
    # -----------------------------
    df = pd.read_excel(INPUT_XLSX)
    ensure_columns(df, ["Sample_ID", "Feature", "Feature value", "MSP_SHAP"])

    # standardize numeric columns
    df["Sample_ID"] = pd.to_numeric(df["Sample_ID"], errors="coerce").astype("Int64")
    df["Feature value"] = pd.to_numeric(df["Feature value"], errors="coerce")
    df["MSP_SHAP"] = pd.to_numeric(df["MSP_SHAP"], errors="coerce")

    # drop rows missing core fields
    df = df.dropna(subset=["Sample_ID", "Feature", "MSP_SHAP"]).copy()

    # Category assignment
    df["Category"] = df["Feature"].map(CATEGORY_MAP).fillna("Other / Unmapped")

    # -----------------------------
    # Global feature importance: mean(|SHAP|)
    # -----------------------------
    feat_importance = (
        df.groupby("Feature")["MSP_SHAP"]
        .apply(lambda s: float(np.mean(np.abs(s.values))))
        .sort_values(ascending=False)
        .rename("mean_abs_shap")
        .reset_index()
    )

    # share of total (nice for reporting)
    total = feat_importance["mean_abs_shap"].sum()
    feat_importance["share"] = feat_importance["mean_abs_shap"] / total if total else np.nan

    # -----------------------------
    # Directionality: corr(feature_value, SHAP)
    # -----------------------------
    # correlation computed on rows where feature values exist
    direction = (
        df.groupby("Feature")
        .apply(lambda g: safe_corr(g["Feature value"].values, g["MSP_SHAP"].values))
        .rename("corr_value_shap")
        .reset_index()
    )

    # Merge global summary
    feat_summary = feat_importance.merge(direction, on="Feature", how="left")

    # Optional: a simple "direction label"
    def dir_label(r):
        c = r["corr_value_shap"]
        if not np.isfinite(c):
            return "non-monotonic/NA"
        if c > 0.15:
            return "↑ feature → ↑ MSP"
        if c < -0.15:
            return "↑ feature → ↓ MSP"
        return "weak/unstable direction"

    feat_summary["direction_interpretation"] = feat_summary.apply(dir_label, axis=1)

    # -----------------------------
    # Category importance: sum mean(|SHAP|) by category
    # -----------------------------
    cat_importance = (
        df.groupby(["Category", "Feature"])["MSP_SHAP"]
        .apply(lambda s: float(np.mean(np.abs(s.values))))
        .rename("mean_abs_shap")
        .reset_index()
    )

    cat_total = (
        cat_importance.groupby("Category")["mean_abs_shap"]
        .sum()
        .sort_values(ascending=False)
        .rename("category_sum_mean_abs_shap")
        .reset_index()
    )

    cat_total_sum = cat_total["category_sum_mean_abs_shap"].sum()
    cat_total["share"] = cat_total["category_sum_mean_abs_shap"] / cat_total_sum if cat_total_sum else np.nan

    # -----------------------------
    # Create wide SHAP table: samples × features
    # -----------------------------
    shap_wide = df.pivot_table(index="Sample_ID", columns="Feature", values="MSP_SHAP", aggfunc="sum")
    # A "total activity" metric to pick a representative sample:
    total_abs = shap_wide.abs().sum(axis=1).sort_values(ascending=False)

    if SAMPLE_ID_FOR_WATERFALL is None:
        sample_id = int(total_abs.index[0])
    else:
        sample_id = int(SAMPLE_ID_FOR_WATERFALL)
        if sample_id not in shap_wide.index:
            raise ValueError(f"Sample_ID {sample_id} not found in file.")

    sample_shap = shap_wide.loc[sample_id].dropna()

    # -----------------------------
    # Plots: global bar + category bar + sample waterfall
    # -----------------------------
    # Global top features bar plot
    top_features = feat_importance.head(20).copy()
    fig = plt.figure(figsize=(9, 6))
    ax = plt.gca()
    ax.barh(top_features["Feature"][::-1], top_features["mean_abs_shap"][::-1])
    ax.set_xlabel("Mean(|SHAP|) for MSP")
    ax.set_title("Global MSP sensitivity (Top 20 features)")
    plt.tight_layout()
    out_png_global = os.path.join(OUT_DIR, f"{OUT_TAG}_GLOBAL_TOP20.png")
    fig.savefig(out_png_global, dpi=200)
    plt.close(fig)

    # Category importance bar plot
    fig = plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.barh(cat_total["Category"][::-1], cat_total["category_sum_mean_abs_shap"][::-1])
    ax.set_xlabel("Sum of mean(|SHAP|) within category")
    ax.set_title("MSP sensitivity by category")
    plt.tight_layout()
    out_png_cat = os.path.join(OUT_DIR, f"{OUT_TAG}_CATEGORY.png")
    fig.savefig(out_png_cat, dpi=200)
    plt.close(fig)

    # Waterfall-like contribution plot for chosen sample
    out_png_wf = os.path.join(OUT_DIR, f"{OUT_TAG}_WATERFALL_SAMPLE_{sample_id}.png")
    sample_waterfall_plot(
        sample_shap=sample_shap,
        out_png=out_png_wf,
        title=f"MSP SHAP contributions (Sample_ID={sample_id})",
        top_n=TOP_N_WATERFALL,
    )

    # -----------------------------
    # Save outputs to Excel
    # -----------------------------
    out_xlsx = os.path.join(OUT_DIR, f"{OUT_TAG}_SUMMARY.xlsx")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="raw_long", index=False)
        feat_summary.to_excel(writer, sheet_name="feature_summary", index=False)
        feat_importance.to_excel(writer, sheet_name="feature_importance", index=False)
        cat_total.to_excel(writer, sheet_name="category_importance", index=False)
        cat_importance.sort_values(["Category", "mean_abs_shap"], ascending=[True, False]).to_excel(
            writer, sheet_name="category_feature_detail", index=False
        )
        shap_wide.to_excel(writer, sheet_name="shap_wide_samples_x_features")
        total_abs.rename("total_abs_shap").reset_index().to_excel(writer, sheet_name="sample_total_abs", index=False)

        # save the selected sample’s contributions
        sample_tbl = sample_shap.rename("MSP_SHAP").reset_index().rename(columns={"index": "Feature"})
        sample_tbl["abs"] = sample_tbl["MSP_SHAP"].abs()
        sample_tbl.sort_values("abs", ascending=False).to_excel(writer, sheet_name=f"sample_{sample_id}_top", index=False)

    print(f"✅ Loaded: {INPUT_XLSX}")
    print(f"✅ Wrote summary Excel: {out_xlsx}")
    print(f"✅ Saved plots:\n  - {out_png_global}\n  - {out_png_cat}\n  - {out_png_wf}")
    print(f"⭐ Waterfall sample used: Sample_ID = {sample_id}")


if __name__ == "__main__":
    main()
