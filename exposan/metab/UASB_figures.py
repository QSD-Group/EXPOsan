#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make two contour-style heatmaps for the UASB system:

1) USD_per_tonne_COD_removed
2) COD_removal_pct

Axes:
    x-axis = COD_g_L
    y-axis = HRT_d

Adds constraint line:
    COD / HRT = 160  (only on COD_removal plot)

Both figures are saved at 300 dpi.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ==========================
# USER INPUT
# ==========================
FILE = "/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/results/METAB_TEA_UASB_1stage_passive_small_HRT.xlsx"
SAVE_DIR = "/Users/saumitrarai/Desktop/Research/EXPOsan/exposan/metab/figures"

XCOL = "COD_g_L"
YCOL = "HRT_d"

PLOTS = [
    {
        "value_col": "USD_per_tonne_COD_removed",
        "title": "UASB system",
        "cbar_label": r"USD tonne$^{-1}$ COD removed",
        "outfile": "UASB_heatmap_USD_per_tonne_COD_removed.png",
        "cmap": "viridis",
        "add_constraint": False,
    },
    {
        "value_col": "COD_removal_pct",
        "title": "UASB system",
        "cbar_label": "COD removal (%)",
        "outfile": "UASB_heatmap_COD_removal_pct.png",
        "cmap": "viridis",
        "add_constraint": True,
    },
]

FIGSIZE = (7, 5.5)
DPI = 300
N_LEVELS = 20

SHOW_POINTS = False


# ==========================
# HELPERS
# ==========================
def require_columns(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")


def build_grid(df, xcol, ycol, value_col):
    heat = df.pivot(index=ycol, columns=xcol, values=value_col)
    heat = heat.sort_index().sort_index(axis=1)

    x = heat.columns.to_numpy(dtype=float)
    y = heat.index.to_numpy(dtype=float)
    z = heat.to_numpy(dtype=float)

    if np.all(np.isnan(z)):
        raise ValueError(f"All values are NaN for {value_col}")

    return x, y, z


def make_heatmap(df, xcol, ycol, value_col, title, cbar_label, outfile, cmap, add_constraint):
    x, y, z = build_grid(df, xcol, ycol, value_col)
    X, Y = np.meshgrid(x, y)

    vmin = np.nanmin(z)
    vmax = np.nanmax(z)

    if np.isclose(vmin, vmax):
        levels = np.linspace(vmin - 1e-9, vmax + 1e-9, N_LEVELS)
    else:
        levels = np.linspace(vmin, vmax, N_LEVELS)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    cf = ax.contourf(
        X,
        Y,
        z,
        levels=levels,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )

    if SHOW_POINTS:
        ax.scatter(X.flatten(), Y.flatten(), s=10, c="k")

    # ==========================
    # ADD CONSTRAINT LINE
    # ==========================
    if add_constraint:
        x_line = np.linspace(x.min(), x.max(), 500)
        y_line = x_line / 160.0

        mask = (y_line >= y.min()) & (y_line <= y.max())

        ax.plot(
            x_line[mask],
            y_line[mask],
            linestyle=":",
            color="red",
            linewidth=2,
            label=r"$\mathrm{COD}/\mathrm{HRT} = 160$",
        )

        ax.legend(loc="upper right")

    ax.set_title(title)
    ax.set_xlabel("COD (g/L)")
    ax.set_ylabel("HRT (d)")

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label(cbar_label)

    save_path = os.path.join(SAVE_DIR, outfile)
    plt.savefig(save_path, dpi=DPI, bbox_inches="tight")
    plt.show()

    print(f"Saved: {save_path}")


# ==========================
# MAIN
# ==========================
def main():
    os.makedirs(SAVE_DIR, exist_ok=True)

    df = pd.read_excel(FILE)

    required = {XCOL, YCOL}
    required.update(p["value_col"] for p in PLOTS)
    require_columns(df, list(required))

    for p in PLOTS:
        make_heatmap(
            df=df,
            xcol=XCOL,
            ycol=YCOL,
            value_col=p["value_col"],
            title=p["title"],
            cbar_label=p["cbar_label"],
            outfile=p["outfile"],
            cmap=p["cmap"],
            add_constraint=p["add_constraint"],
        )


if __name__ == "__main__":
    main()