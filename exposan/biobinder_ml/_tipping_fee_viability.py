# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""
Build tipping-fee viability dataset and plots from latest META files.

Outputs:
1) Excel workbook with:
   - Summary
   - Long_All_Cases
   - Viable_IRR_ge_10
2) CSV summary
3) One density plot per feedstock x config x product set:
   - IRR >= 10%
   - IRR < 10%
   - non-parametric intersection threshold (if identifiable)

Notes:
- Uses latest META file per feedstock/config/product
- Biobinder finder prioritizes MERGED files
- SAF finder uses latest timestamped file
- Tipping fee is derived as:
      Tipping_fee_$/wet_tonne = - Feedstock_price_$/tonne
- Only explicit removals are applied:
      97.90195447 for Biobinder + manure
      1609372.922 for SAF
      fog feedstock excluded
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, ks_2samp

# =========================================================
# USER SETTINGS
# =========================================================
output_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results"
prefer = "max"
irr_target = 10.0

feedstocks = ["food", "sludge", "manure", "green"]
configs = ["CHCU", "DHCU"]

# Explicit removals only
EXCLUDE_FEEDSTOCKS = {"fog"}
ARTIFACT_BIOBINDER_MANURE = 97.90195447
ARTIFACT_SAF_EXTREME = 1609372.922

# Output files
out_excel = os.path.join(output_dir, "tipping_fee_thresholds_intersection.xlsx")
out_csv = os.path.join(output_dir, "tipping_fee_thresholds_intersection.csv")
plot_dir = os.path.join(output_dir, "tipping_fee_threshold_plots_all")
os.makedirs(plot_dir, exist_ok=True)

# =========================================================
# META FINDERS
# =========================================================
def find_latest_meta(feedstock_id, config_name, prefer, output_dir="results"):
    patterns = [
        os.path.join(
            output_dir,
            f"META_SHAP3_10000_{feedstock_id}_{config_name}_*_prefer-{prefer}*.xlsx"
        ),
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

# =========================================================
# HELPERS
# =========================================================
def first_existing(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    return None


def load_meta_file(path, feedstock, config, product):
    """
    Read one META file and return a standardized long dataframe.
    Tipping fee is derived from Feedstock_price_$/tonne.
    """
    df = pd.read_excel(path)

    irr_col = first_existing(
        df.columns,
        ["IRR_pct", "IRR", "IRR %"]
    )

    price_col = first_existing(
        df.columns,
        ["Feedstock_price_$/tonne", "Feedstock_price", "Feedstock price"]
    )

    if irr_col is None:
        raise ValueError(
            f"'IRR_pct' column not found in {path}\n"
            f"Available columns:\n{list(df.columns)}"
        )

    if price_col is None:
        raise ValueError(
            f"'Feedstock_price_$/tonne' column not found in {path}\n"
            f"Available columns:\n{list(df.columns)}"
        )

    out = pd.DataFrame({
        "Feedstock": feedstock,
        "Config": config,
        "Product set": product,
        "Source file": os.path.basename(path),
        "IRR_pct": pd.to_numeric(df[irr_col], errors="coerce"),
        "Feedstock_price_$/tonne": pd.to_numeric(df[price_col], errors="coerce"),
    })

    # Derive tipping fee exactly as in your earlier workflow
    out["Tipping_fee_$/wet_tonne"] = -out["Feedstock_price_$/tonne"]

    return out


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

    if np.allclose(np.std(x), 0):
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(sub)

    slope, intercept = np.polyfit(x, y, 1)
    y_pred = intercept + slope * x

    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    if np.isclose(ss_tot, 0):
        r2 = np.nan
    else:
        r2 = 1 - (ss_res / ss_tot)

    rmse = np.sqrt(np.mean((y - y_pred) ** 2))

    if np.isclose(slope, 0):
        threshold = np.nan
    else:
        threshold = (target_irr - intercept) / slope

    return intercept, slope, threshold, r2, rmse, len(sub)


def compute_intersection(df, irr_target=10.0, gridsize=1000, bw_method=None):
    """
    Non-parametric threshold from intersection of KDEs:
      - viable: IRR >= irr_target
      - non-viable: IRR < irr_target

    Returns:
        threshold
        ks_stat
        ks_p
        n_viable
        n_nonviable
        xgrid
        y_viable
        y_nonviable
        all_intersections
    """
    viable = df.loc[df["IRR_pct"] >= irr_target, "Tipping_fee_$/wet_tonne"].dropna().to_numpy(dtype=float)
    nonviable = df.loc[df["IRR_pct"] < irr_target, "Tipping_fee_$/wet_tonne"].dropna().to_numpy(dtype=float)

    if len(viable) < 10 or len(nonviable) < 10:
        return np.nan, np.nan, np.nan, len(viable), len(nonviable), None, None, None, np.array([])

    xmin = min(viable.min(), nonviable.min())
    xmax = max(viable.max(), nonviable.max())
    xgrid = np.linspace(xmin, xmax, gridsize)

    kde_v = gaussian_kde(viable, bw_method=bw_method)
    kde_nv = gaussian_kde(nonviable, bw_method=bw_method)

    y_v = kde_v(xgrid)
    y_nv = kde_nv(xgrid)

    diff = y_v - y_nv
    idx = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]

    intersections = []
    for i in idx:
        x1, x2 = xgrid[i], xgrid[i + 1]
        y1, y2 = diff[i], diff[i + 1]

        if np.isclose(y2 - y1, 0):
            x_star = x1
        else:
            x_star = x1 - y1 * (x2 - x1) / (y2 - y1)

        intersections.append(x_star)

    if len(intersections) == 0:
        chosen = np.nan
    else:
        intersections = np.array(intersections)
        med_mid = 0.5 * (np.median(viable) + np.median(nonviable))
        chosen = intersections[np.argmin(np.abs(intersections - med_mid))]

    ks_stat, ks_p = ks_2samp(viable, nonviable)

    return chosen, ks_stat, ks_p, len(viable), len(nonviable), xgrid, y_v, y_nv, np.array(intersections)


def make_density_plot(
    df_case,
    feedstock,
    config,
    product,
    irr_target,
    threshold_np,
    threshold_linear,
    ks_p,
    out_png,
    out_pdf,
    bw_method=None,
):
    viable = df_case.loc[df_case["IRR_pct"] >= irr_target, "Tipping_fee_$/wet_tonne"].dropna().to_numpy(dtype=float)
    nonviable = df_case.loc[df_case["IRR_pct"] < irr_target, "Tipping_fee_$/wet_tonne"].dropna().to_numpy(dtype=float)

    if len(viable) < 10 or len(nonviable) < 10:
        print(f"Skipping plot for {feedstock} | {config} | {product}: not enough viable/nonviable cases.")
        return

    xmin = min(viable.min(), nonviable.min())
    xmax = max(viable.max(), nonviable.max())
    xgrid = np.linspace(xmin, xmax, 1000)

    kde_v = gaussian_kde(viable, bw_method=bw_method)
    kde_nv = gaussian_kde(nonviable, bw_method=bw_method)

    y_v = kde_v(xgrid)
    y_nv = kde_nv(xgrid)

    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=300)

    # Colors similar to the reference style
    color_viable = "#1f77b4"     # blue
    color_nonviable = "#d62728"  # red

    ax.plot(xgrid, y_v, color=color_viable, linewidth=2.4, label=f"IRR ≥ {irr_target:.0f}%")
    ax.plot(xgrid, y_nv, color=color_nonviable, linewidth=2.4, label=f"IRR < {irr_target:.0f}%")

    ymax = max(y_v.max(), y_nv.max())

    # Non-parametric intersection threshold
    if np.isfinite(threshold_np):
        ax.axvline(threshold_np, linestyle="--", linewidth=1.8, color="black")
        ax.annotate(
            f"target = {threshold_np:.1f}",
            xy=(threshold_np, ymax * 0.93),
            xytext=(threshold_np + 0.05 * (xgrid.max() - xgrid.min()), ymax * 1.02),
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color="black"),
            fontsize=10.5,
            fontweight="bold",
        )

    # Optional linear threshold as a faint gray line
    if np.isfinite(threshold_linear):
        ax.axvline(threshold_linear, linestyle=":", linewidth=1.4, color="gray")

    ax.set_title(f"{feedstock.capitalize()} | {config} | {product}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Tipping fee [$/wet metric tonne]", fontsize=12, fontweight="bold")
    # ax.set_ylabel("Density", fontsize=12, fontweight="bold")
    ax.set_ylabel("")


    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold")
        tick.set_fontsize(10.5)
    # for tick in ax.get_yticklabels():
    #     tick.set_fontweight("bold")
    #     tick.set_fontsize(10.5)
    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.legend(frameon=False, loc="upper left", fontsize=10.5)

    if np.isfinite(ks_p):
        ax.text(
            0.98, 0.96, f"KS p = {ks_p:.3g}",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=10, fontweight="bold"
        )

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# =========================================================
# BUILD LONG DATASET
# =========================================================
all_rows = []
file_log = []

for fs in feedstocks:
    if fs in EXCLUDE_FEEDSTOCKS:
        continue

    for cfg in configs:
        # -------------------------
        # Biobinder
        # -------------------------
        try:
            fpath = find_latest_meta(fs, cfg, prefer, output_dir)
            df_case = load_meta_file(fpath, fs, cfg, "Biobinder")
            all_rows.append(df_case)

            file_log.append({
                "Feedstock": fs,
                "Config": cfg,
                "Product set": "Biobinder",
                "Source file": os.path.basename(fpath),
                "N rows loaded": len(df_case)
            })
        except Exception as e:
            print(f"[Biobinder] {fs} | {cfg} -> {e}")

        # -------------------------
        # SAF
        # -------------------------
        try:
            fpath = find_latest_meta_saf(fs, cfg, prefer, output_dir)
            df_case = load_meta_file(fpath, fs, cfg, "SAF")
            all_rows.append(df_case)

            file_log.append({
                "Feedstock": fs,
                "Config": cfg,
                "Product set": "SAF",
                "Source file": os.path.basename(fpath),
                "N rows loaded": len(df_case)
            })
        except Exception as e:
            print(f"[SAF] {fs} | {cfg} -> {e}")

if not all_rows:
    raise ValueError("No valid META files were found.")

long_df = pd.concat(all_rows, ignore_index=True)

# =========================================================
# CLEANING
# =========================================================
# Keep valid numeric IRR and tipping-fee rows
long_df = long_df.dropna(subset=["IRR_pct", "Tipping_fee_$/wet_tonne"]).copy()

# Explicit artifact removals only
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
long_df["IRR_ge_10"] = long_df["IRR_pct"] >= irr_target

# =========================================================
# SUMMARY TABLE
# =========================================================
summary_rows = []

group_cols = ["Feedstock", "Config", "Product set", "Source file"]

for keys, grp in long_df.groupby(group_cols):
    fs, cfg, prod, src = keys

    intercept, slope, threshold_linear, r2, rmse, n_fit = estimate_threshold_linear(grp, irr_target)

    threshold_np, ks_stat, ks_p, n_viable, n_nonviable, _, _, _, intersections = compute_intersection(
        grp, irr_target=irr_target, gridsize=1000, bw_method=None
    )

    summary_rows.append({
        "Feedstock": fs,
        "Config": cfg,
        "Product set": prod,
        "Source file": src,
        "N all": len(grp),
        "N viable (IRR>=10)": int((grp["IRR_pct"] >= irr_target).sum()),
        "Viability fraction": float((grp["IRR_pct"] >= irr_target).mean()),
        "IRR min": grp["IRR_pct"].min(),
        "IRR p10": grp["IRR_pct"].quantile(0.10),
        "IRR p25": grp["IRR_pct"].quantile(0.25),
        "IRR median": grp["IRR_pct"].median(),
        "IRR p75": grp["IRR_pct"].quantile(0.75),
        "IRR p90": grp["IRR_pct"].quantile(0.90),
        "IRR max": grp["IRR_pct"].max(),
        "Tipping fee min": grp["Tipping_fee_$/wet_tonne"].min(),
        "Tipping fee p10": grp["Tipping_fee_$/wet_tonne"].quantile(0.10),
        "Tipping fee p25": grp["Tipping_fee_$/wet_tonne"].quantile(0.25),
        "Tipping fee median": grp["Tipping_fee_$/wet_tonne"].median(),
        "Tipping fee p75": grp["Tipping_fee_$/wet_tonne"].quantile(0.75),
        "Tipping fee p90": grp["Tipping_fee_$/wet_tonne"].quantile(0.90),
        "Tipping fee max": grp["Tipping_fee_$/wet_tonne"].max(),
        "Linear fit intercept": intercept,
        "Linear fit slope": slope,
        "Linear fit R2": r2,
        "Linear fit RMSE": rmse,
        "N used in fit": n_fit,
        "Estimated tipping fee for IRR=10 (linear)": threshold_linear,
        "Intersection tipping fee for IRR=10 (nonparametric)": threshold_np,
        "KS statistic": ks_stat,
        "KS p-value": ks_p,
        "n viable in intersection": n_viable,
        "n nonviable in intersection": n_nonviable,
        "n intersections found": len(intersections),
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
nonviable_df = long_df[~long_df["IRR_ge_10"]].copy()
file_log_df = pd.DataFrame(file_log)

# ========================================================
# SAVE EXCEL / CSV
# =========================================================
with pd.ExcelWriter(out_excel, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    long_df.to_excel(writer, sheet_name="Long_All_Cases", index=False)
    viable_df.to_excel(writer, sheet_name="Viable_IRR_ge_10", index=False)
    nonviable_df.to_excel(writer, sheet_name="Nonviable_IRR_lt_10", index=False)
    file_log_df.to_excel(writer, sheet_name="File_Log", index=False)

summary_df.to_csv(out_csv, index=False)

print("Saved:")
print(out_excel)
print(out_csv)

# =========================================================
# PLOTTING
# =========================================================
for _, row in summary_df.iterrows():
    fs = row["Feedstock"]
    cfg = row["Config"]
    prod = row["Product set"]

    case_df = long_df[
        (long_df["Feedstock"] == fs) &
        (long_df["Config"] == cfg) &
        (long_df["Product set"] == prod)
    ].copy()

    png_name = f"{fs}_{cfg}_{prod}_intersection_density.png"
    pdf_name = f"{fs}_{cfg}_{prod}_intersection_density.pdf"

    make_density_plot(
        df_case=case_df,
        feedstock=str(fs),
        config=str(cfg),
        product=str(prod),
        irr_target=irr_target,
        threshold_np=row["Intersection tipping fee for IRR=10 (nonparametric)"],
        threshold_linear=row["Estimated tipping fee for IRR=10 (linear)"],
        ks_p=row["KS p-value"],
        out_png=os.path.join(plot_dir, png_name),
        out_pdf=os.path.join(plot_dir, pdf_name),
        bw_method=None,  # try 0.5 for smoother curves if needed
    )

print(f"Plots saved to: {plot_dir}")

# =========================================================
# ORIGIN
# =========================================================
for _, row in summary_df.iterrows():
    fs = row["Feedstock"]
    cfg = row["Config"]
    prod = row["Product set"]

    case_df = long_df[
        (long_df["Feedstock"] == fs) &
        (long_df["Config"] == cfg) &
        (long_df["Product set"] == prod)
    ].copy()

    viable = case_df.loc[
        case_df["IRR_pct"] >= irr_target,
        "Tipping_fee_$/wet_tonne"
    ].reset_index(drop=True)

    nonviable = case_df.loc[
        case_df["IRR_pct"] < irr_target,
        "Tipping_fee_$/wet_tonne"
    ].reset_index(drop=True)

    origin_df = pd.concat(
        [
            viable.rename("IRR_ge_10"),
            nonviable.rename("IRR_lt_10")
        ],
        axis=1
    )

    origin_df.to_excel(
        os.path.join(plot_dir, f"{fs}_{cfg}_{prod}_Origin_KDE.xlsx"),
        index=False
    )
n_v = len(viable)
n_nv = len(nonviable)

print(fs, cfg, prod, "viable:", n_v, "nonviable:", n_nv)

if n_v == 0:
    print(f"Only nonviable cases for {fs}_{cfg}_{prod}; no IRR ≥ {irr_target}% curve.")

if n_nv == 0:
    print(f"Only viable cases for {fs}_{cfg}_{prod}; no IRR < {irr_target}% curve.")
    
# =========================================================
# ORIGIN WIDE FORMAT: 32 COLUMNS
# =========================================================
config_label = {
    "CHCU": "c-HTL",
    "DHCU": "d-HTL",
}

product_label = {
    "Biobinder": "BB",
    "SAF": "SAF",
}

origin_cols = []

for fs in feedstock_order:
    for cfg in config_order:
        for prod in product_order:

            case_df = long_df[
                (long_df["Feedstock"].astype(str) == fs) &
                (long_df["Config"].astype(str) == cfg) &
                (long_df["Product set"].astype(str) == prod)
            ].copy()

            fs_label = fs.capitalize()
            cfg_label = config_label.get(cfg, cfg)
            prod_label = product_label.get(prod, prod)

            viable = case_df.loc[
                case_df["IRR_pct"] >= irr_target,
                "Tipping_fee_$/wet_tonne"
            ].reset_index(drop=True)

            nonviable = case_df.loc[
                case_df["IRR_pct"] < irr_target,
                "Tipping_fee_$/wet_tonne"
            ].reset_index(drop=True)

            origin_cols.append(
                viable.rename(f"{fs_label}/{cfg_label}/{prod_label} (IRR_ge_10)")
            )

            origin_cols.append(
                nonviable.rename(f"{fs_label}/{cfg_label}/{prod_label} (IRR_lt_10)")
            )

origin_wide_df = pd.concat(origin_cols, axis=1)

origin_wide_excel = os.path.join(
    output_dir,
    "Origin_Tipping_Fee_Threshold_32cols.xlsx"
)

origin_wide_df.to_excel(origin_wide_excel, index=False)

print("Saved Origin wide-format file:")
print(origin_wide_excel)