# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 15:25:14 2026

@author: aliah
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===============================================================
# File paths
# ===============================================================
cv_file = r"results/Results_10Fold_DNN_RF_GBR_Extended.xlsx"   # base RF CV workbook
test_file = r"results/Results_10Fold_DNN_Softmax_RF_GBR_20260418_1337.xlsx"  # tuned RF test workbook

out_dir = "results/bonus_defense_figures"
os.makedirs(out_dir, exist_ok=True)

# ===============================================================
# Settings
# ===============================================================
yield_order = ["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"]
yield_labels = ["BioY%", "AqY%", "Gas%", "Solids%"]
metrics = ["R2", "RMSE", "MAE", "MAPE"]

cv_model_map = {
    "DNN(MLP)": "DNN CV",
    "GradientBoosting": "GBR CV",
    "RandomForest": "RF CV",
}

test_model_map = {
    "DNN(Softmax)_Test": "DNN Test",
    "GradientBoosting_Test": "GBR Test",
    "RandomForest_Tuned_Test": "Tuned RF Test",
}

# ===============================================================
# Load CV summary (per-target means from extended workbook)
# ===============================================================
cv_target = pd.read_excel(cv_file, sheet_name="CV_Target_Summary")

# Keep only rows we need
cv_target = cv_target[cv_target["Model"].isin(cv_model_map.keys())].copy()
cv_target["ModelLabel"] = cv_target["Model"].map(cv_model_map)

# ===============================================================
# Load Test results (raw per-target from tuned workbook)
# ===============================================================
test_raw = pd.read_excel(test_file, sheet_name="TestSet_Results")

# Keep only rows we need, excluding macro rows
test_raw = test_raw[
    (test_raw["Model"].isin(test_model_map.keys())) &
    (test_raw["Target"] != "__macro__")
].copy()
test_raw["ModelLabel"] = test_raw["Model"].map(test_model_map)

# ===============================================================
# Helper to plot grouped bars
# ===============================================================
def make_grouped_barplot(df, metric, split_name, model_order, title_suffix, save_name):
    """
    df must contain columns:
    - ModelLabel
    - Target
    - metric
    """
    # Pivot into [target x model]
    piv = df.pivot(index="Target", columns="ModelLabel", values=metric).reindex(yield_order)
    piv = piv[model_order]

    x = np.arange(len(yield_order))
    width = 0.24

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, model in enumerate(model_order):
        vals = piv[model].values
        bars = ax.bar(x + (i - 1) * width, vals, width=width, label=model)
        for b, v in zip(bars, vals):
            if pd.notna(v):
                ax.text(
                    b.get_x() + b.get_width() / 2,
                    b.get_height() + (0.01 if metric == "R2" else max(vals) * 0.01),
                    f"{v:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=10
                )

    ax.set_xticks(x)
    ax.set_xticklabels(yield_labels, fontsize=12, fontweight="bold")
    ax.set_xlabel("Yield Type", fontsize=13, fontweight="bold")
    ax.set_ylabel(metric, fontsize=13, fontweight="bold")
    ax.set_title(f"{metric} by yield type ({title_suffix})", fontsize=16, fontweight="bold")
    ax.legend(frameon=False, fontsize=11)

    # Better axis formatting for R2
    if metric == "R2":
        ymin = max(0, np.nanmin(piv.values) - 0.05)
        ymax = min(1.0, np.nanmax(piv.values) + 0.08)
        ax.set_ylim(ymin, ymax)
    else:
        ymax = np.nanmax(piv.values) * 1.18
        ax.set_ylim(0, ymax)

    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    save_path = os.path.join(out_dir, save_name)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"✅ Saved: {save_path}")

# ===============================================================
# Generate 4 CV plots
# ===============================================================
cv_order = ["DNN CV", "GBR CV", "RF CV"]

for metric in metrics:
    make_grouped_barplot(
        df=cv_target[["ModelLabel", "Target", metric]].copy(),
        metric=metric,
        split_name="CV",
        model_order=cv_order,
        title_suffix="10-fold CV",
        save_name=f"{metric}_CV_barplot.png"
    )

# ===============================================================
# Generate 4 Test plots
# ===============================================================
test_order = ["DNN Test", "GBR Test", "Tuned RF Test"]

for metric in metrics:
    make_grouped_barplot(
        df=test_raw[["ModelLabel", "Target", metric]].copy(),
        metric=metric,
        split_name="Test",
        model_order=test_order,
        title_suffix="holdout test set",
        save_name=f"{metric}_Test_barplot.png"
    )

print("\n🎉 Done. Generated 8 defense-ready figures.")