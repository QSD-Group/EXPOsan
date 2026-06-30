#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot pairwise cost difference heatmaps between FB and UASB systems.

Difference plotted:
    FB - UASB

Metric:
    USD_per_tonne_COD_removed

Panels:
    bead_lifetime_yr = 0.50, 0.75, 1.00

Interpretation:
    Positive values  -> FB is more expensive than UASB
    Negative values  -> FB is cheaper than UASB
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import griddata


# =========================================================
# USER INPUTS
# =========================================================
FB_FILE = r"/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/results/METAB_TEA_FB_1stage_passive_degas_6.5.xlsx"
UASB_FILE = r"/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/results/METAB_TEA_UASB_1stage_passive.xlsx"

OUT_DIR = r"/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/figures"
OUT_FILE = os.path.join(OUT_DIR, "6.5_FB_minus_UASB_pairwise_difference_interpolated_300dpi.png")

LIFETIMES_TO_PLOT = [1.00]

XCOL = "COD_g_L"
YCOL = "HRT_d"
VALUE_COL = "USD_per_tonne_COD_removed"
LT_COL = "bead_lifetime_yr"

FIG_W = 15
FIG_H = 4.8
DPI = 300


# =========================================================
# HELPERS
# =========================================================
def require_columns(df, cols, df_name):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{df_name} is missing required columns: {missing}")


def check_unique_grid(df, name, include_lt=False):
    keys = [XCOL, YCOL]
    if include_lt:
        keys = [LT_COL] + keys

    dup = df.duplicated(subset=keys, keep=False)
    if dup.any():
        raise ValueError(
            f"{name} has duplicate rows for grid keys {keys}. "
            f"Fix duplicates before plotting.\n"
            f"Example duplicates:\n{df.loc[dup, keys + [VALUE_COL]].head()}"
        )


def build_delta_grid(fb_sub, uasb_sub):
    merged = fb_sub.merge(
        uasb_sub,
        on=[XCOL, YCOL],
        how="inner",
        suffixes=("_fb", "_uasb"),
    )

    if merged.empty:
        raise ValueError("Merge between FB and UASB returned no rows.")

    merged["delta"] = merged[f"{VALUE_COL}_fb"] - merged[f"{VALUE_COL}_uasb"]

    pivot = merged.pivot(index=YCOL, columns=XCOL, values="delta")
    pivot = pivot.sort_index().sort_index(axis=1)

    x = pivot.columns.to_numpy(dtype=float)
    y = pivot.index.to_numpy(dtype=float)
    z = pivot.to_numpy(dtype=float)

    return x, y, z, merged

# =========================================================
# MAIN
# =========================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    fb = pd.read_excel(FB_FILE)
    uasb = pd.read_excel(UASB_FILE)

    require_columns(fb, [XCOL, YCOL, LT_COL, VALUE_COL], "FB file")
    require_columns(uasb, [XCOL, YCOL, VALUE_COL], "UASB file")

    fb = fb[[XCOL, YCOL, LT_COL, VALUE_COL]].copy()
    uasb = uasb[[XCOL, YCOL, VALUE_COL]].copy()

    fb[XCOL] = fb[XCOL].astype(float)
    fb[YCOL] = fb[YCOL].astype(float)
    fb[LT_COL] = fb[LT_COL].astype(float)
    fb[VALUE_COL] = fb[VALUE_COL].astype(float)

    uasb[XCOL] = uasb[XCOL].astype(float)
    uasb[YCOL] = uasb[YCOL].astype(float)
    uasb[VALUE_COL] = uasb[VALUE_COL].astype(float)

    check_unique_grid(fb, "FB file", include_lt=True)
    check_unique_grid(uasb, "UASB file", include_lt=False)

    # Build all delta grids first so all panels use the same color scale
    grids = {}
    all_abs_vals = []

    for lt in LIFETIMES_TO_PLOT:
        fb_sub = fb[np.isclose(fb[LT_COL], lt)].copy()
        if fb_sub.empty:
            raise ValueError(f"No FB rows found for bead_lifetime_yr = {lt}")

        x, y, z, merged = build_delta_grid(fb_sub, uasb)
        grids[lt] = (x, y, z, merged)

        finite = z[np.isfinite(z)]
        if finite.size > 0:
            all_abs_vals.append(np.abs(finite))

    if not all_abs_vals:
        raise ValueError("All delta grids are empty or NaN.")

    vmax = np.nanmax(np.concatenate(all_abs_vals))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, axes = plt.subplots(1, 3, figsize=(FIG_W, FIG_H), constrained_layout=True)

    mappable = None

    for ax, lt in zip(axes, LIFETIMES_TO_PLOT):
        x, y, z, merged = grids[lt]

        

        # Original points
        X, Y = np.meshgrid(x, y)
        points = np.column_stack([X.flatten(), Y.flatten()])
        values = z.flatten()
        
        # Fine grid
        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        XI, YI = np.meshgrid(xi, yi)
        
        # Interpolate
        ZI = griddata(points, values, (XI, YI), method="linear")
        ZI_nearest = griddata(points, values, (XI, YI), method="nearest")
        ZI = np.where(np.isnan(ZI), ZI_nearest, ZI)
        
        # Plot
        mappable = ax.contourf(
            XI,
            YI,
            ZI,
            levels=50,
            cmap="RdBu_r",
            norm=norm,
        )

        ax.set_title(f"FB ({lt:.2f} yr) - UASB", fontsize=11)
        ax.set_xlabel("COD (g/L)")
        ax.set_ylabel("HRT (d)")

        ax.set_xticks(x)
        ax.set_yticks(y)

    cbar = fig.colorbar(mappable, ax=axes, shrink=0.95)
    cbar.set_label("Δ USD/tonne COD removed (FB - UASB)")

    fig.suptitle("Pairwise cost difference across COD-HRT space", fontsize=13)
    fig.savefig(OUT_FILE, dpi=DPI, bbox_inches="tight")
    plt.show()

    print(f"Saved: {OUT_FILE}")


if __name__ == "__main__":
    main()