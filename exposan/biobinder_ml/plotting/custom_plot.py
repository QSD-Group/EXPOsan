# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from shap import TreeExplainer
from joblib import load
from exposan.biobinder_ml._feedstocks import get_feedstock_composition
from exposan.biobinder_ml._ml_features import FEATURE_COLS, make_feature_row, PREPROCESSING_MAP
from exposan.biobinder_ml._ml_features import CATALYST_MAP, SOLVENT_MAP, REACTOR_TYPE_MAP
from datetime import datetime
import os


def decode_feature_value(feature, value):
    try:
        iv = int(value)
        is_int_like = abs(float(value) - iv) < 1e-9
    except Exception:
        is_int_like = False
        iv = value

    if feature == "Pre-processing" and is_int_like:
        return {0: "No", 1: "Yes"}.get(iv, f"Unknown ({iv})")

    if feature == "Catalyst" and is_int_like:
        return CATALYST_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Solvent" and is_int_like:
        return SOLVENT_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Reactor Type" and is_int_like:
        return REACTOR_TYPE_MAP.get(iv, f"Unknown ({iv})")

    return value

def format_display_value(feature, value):
    decoded = decode_feature_value(feature, value)
    if feature == "Reactor Volume (mL)":
        return f"{float(value) / 1000:.2f}"
    if isinstance(decoded, (int, float, np.integer, np.floating)):
        return f"{float(decoded):.2f}"
    return str(decoded)

# -------------------------------------------------
# Custom waterfall plot
# -------------------------------------------------
def custom_waterfall_plot(
    shap_values,
    feature_names,
    feature_display_values,
    base_value,
    final_value,
    save_path,
    max_display=10,
    figsize=(11, 7),
    pos_color="#ff0051",
    neg_color="#1e88e5",
    final_label="final prediction",
    baseline_label="baseline value",
):
    import numpy as np
    import matplotlib.pyplot as plt

    shap_values = np.asarray(shap_values).flatten()
    feature_names = list(feature_names)
    feature_display_values = list(feature_display_values)

    # -------------------------------------------------
    # Sort by absolute SHAP magnitude
    # -------------------------------------------------
    order = np.argsort(np.abs(shap_values))  # ascending

    # Keep top features only for the main waterfall
    # "Remaining features" is shown separately at the top
    if len(order) > max_display - 1:
        top_idx = order[-(max_display - 1):]
        other_idx = order[:-(max_display - 1)]

        top_shap = np.array([shap_values[i] for i in top_idx], dtype=float)
        top_names = [feature_names[i] for i in top_idx]
        top_vals = [feature_display_values[i] for i in top_idx]

        others_shap = float(np.sum(shap_values[other_idx]))
        n_other = len(other_idx)
        show_other = True
    else:
        top_idx = order
        top_shap = np.array([shap_values[i] for i in top_idx], dtype=float)
        top_names = [feature_names[i] for i in top_idx]
        top_vals = [feature_display_values[i] for i in top_idx]

        others_shap = 0.0
        n_other = 0
        show_other = False

    # -------------------------------------------------
    # Reverse so largest abs SHAP is at the bottom
    # -------------------------------------------------
    plot_shap = top_shap[::-1]
    plot_names = top_names[::-1]
    plot_vals = top_vals[::-1]

    n_main = len(plot_shap)
    n_total = n_main + (1 if show_other else 0)

    fig, ax = plt.subplots(figsize=figsize)

    # -------------------------------------------------
    # Build waterfall from baseline using ONLY main features
    # -------------------------------------------------
    running = base_value
    starts = []
    ends = []

    for sv in plot_shap:
        starts.append(running)
        running += sv
        ends.append(running)

    # -------------------------------------------------
    # Draw main bars
    # -------------------------------------------------
    for i, (sv, x0, x1) in enumerate(zip(plot_shap, starts, ends)):
        left = min(x0, x1)
        width = abs(sv)
        color = pos_color if sv >= 0 else neg_color

        ax.barh(
            y=i,
            width=width,
            left=left,
            height=0.82,
            color=color,
            edgecolor="white",
            linewidth=1.2,
            zorder=3,
        )

        txt = f"{sv:+.2f}"
        xm = left + width / 2

        if width > 0.90:
            ax.text(
                xm, i, txt,
                ha="center", va="center",
                fontsize=13, fontweight="bold",
                color="white", zorder=4
            )
        else:
            ax.text(
                x1 + (0.18 if sv >= 0 else -0.18),
                i,
                txt,
                ha="left" if sv >= 0 else "right",
                va="center",
                fontsize=13, fontweight="bold",
                color=color,
                zorder=4
            )

    # -------------------------------------------------
    # Plot "Remaining features" at the TOP
    # -------------------------------------------------
    all_x_for_limits = [base_value, final_value] + starts + ends

    if show_other:
        y_other = n_main
        other_name = "Remaining features"
        other_val = ""
        other_color = pos_color if others_shap >= 0 else neg_color

        other_left = final_value if others_shap >= 0 else final_value + others_shap
        other_width = abs(others_shap)

        ax.barh(
            y=y_other,
            width=other_width,
            left=other_left,
            height=0.82,
            color=other_color,
            edgecolor="white",
            linewidth=1.2,
            zorder=3,
            alpha=0.95,
        )

        txt = f"{others_shap:+.2f}"
        xm = other_left + other_width / 2

        if other_width > 0.90:
            ax.text(
                xm, y_other, txt,
                ha="center", va="center",
                fontsize=13, fontweight="bold",
                color="white", zorder=4
            )
        else:
            label_offset = 0.50
            ax.text(
                final_value + (label_offset if others_shap >= 0 else -label_offset),
                y_other,
                txt,
                ha="left" if others_shap >= 0 else "right",
                va="center",
                fontsize=13, fontweight="bold",
                color=other_color,
                zorder=4
            )

        plot_names = plot_names + [other_name]
        plot_vals = plot_vals + [other_val]

        all_x_for_limits += [other_left, other_left + other_width]

    # -------------------------------------------------
    # Left labels
    # -------------------------------------------------
    ylabels = []
    for val, name in zip(plot_vals, plot_names):
        if val == "":
            ylabels.append(name)
        else:
            ylabels.append(f"{name} ({val})")

    ax.set_yticks(np.arange(n_total))
    ax.set_yticklabels(ylabels, fontsize=15, fontweight="bold")

    # -------------------------------------------------
    # Gridlines
    # -------------------------------------------------
    ax.grid(axis="y", linestyle=(0, (1, 4)), color="0.72", linewidth=0.8, zorder=0)
    ax.grid(False, axis="x")

    # -------------------------------------------------
    # Guide lines
    # -------------------------------------------------
    ax.axvline(base_value, color="0.72", linestyle="--", linewidth=1.2, zorder=1)
    ax.axvline(final_value, color="0.72", linestyle="--", linewidth=1.2, zorder=1)

    # -------------------------------------------------
    # Labels for baseline and final prediction
    # -------------------------------------------------
    plt.subplots_adjust(bottom=0.18, top=0.90)

    ax.text(
        base_value, -0.08,
        f"{baseline_label} = {base_value:.2f}",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top",
        fontsize=16, fontweight="bold",
        color="black",
        clip_on=False,
    )

    ax.text(
        final_value, 1.035,
        f"{final_label} = {final_value:.2f}",
        transform=ax.get_xaxis_transform(),
        ha="center", va="bottom",
        fontsize=17, fontweight="bold",
        color="black",
        clip_on=False,
    )

    # -------------------------------------------------
    # Axis styling
    # -------------------------------------------------
    ax.tick_params(axis="x", labelsize=15, width=1.2)
    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold")

    ax.tick_params(axis="y", length=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.2)

    # -------------------------------------------------
    # Limits / padding
    # -------------------------------------------------
    xmin, xmax = min(all_x_for_limits), max(all_x_for_limits)
    xrng = xmax - xmin if xmax > xmin else 1.0
    pad_left = 0.28 * xrng
    pad_right = 0.22 * xrng
    ax.set_xlim(xmin - pad_left, xmax + pad_right)

    ax.set_ylim(-0.5, n_total - 0.5 + 0.35)

    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    # -------------------------------------------------
    # Load model
    # -------------------------------------------------
    rf_model = load("results/rf_yield_model.joblib")
    explainer = TreeExplainer(rf_model)
    print("RF model loaded")

    # -------------------------------------------------
    # Feedstock → dry-basis composition
    # -------------------------------------------------
    feedstock_id = "food"
    wet_comp = get_feedstock_composition(feedstock_id)

    dry_frac = 1.0 - wet_comp["Water"]

    carb_wt    = wet_comp["Carbohydrates"] / dry_frac * 100
    protein_wt = wet_comp["Proteins"]      / dry_frac * 100
    lipids_wt  = wet_comp["Lipids"]        / dry_frac * 100
    ash_wt     = wet_comp["Ash"]           / dry_frac * 100

    process_conditions = {
        "Temperature (C)": 280,
        "Residence Time": 15,
        "Solid content (w/w) %": 20,
    }

    X = make_feature_row(
        carb_wt=carb_wt,
        protein_wt=protein_wt,
        lipids_wt=lipids_wt,
        ash_wt=ash_wt,
        process=process_conditions,
    )

    y_pred = rf_model.predict(X)[0]

    shap_vals = explainer.shap_values(X)
    feature_values = X.iloc[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    feedstock = feedstock_id

    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    output_names = ["Biocrude", "Aqueous", "Gas", "Char"]

    if isinstance(shap_vals, list):
        shap_arr = np.column_stack([np.array(v)[0, :] for v in shap_vals])
    else:
        shap_vals_arr = np.array(shap_vals)
        if shap_vals_arr.ndim == 3:
            shap_arr = shap_vals_arr[0]
        else:
            raise ValueError(f"Unexpected SHAP shape: {shap_vals_arr.shape}")

    output_idx = 0
    output_name = output_names[output_idx]

    local_shap = shap_arr[:, output_idx]
    base_values = np.array(explainer.expected_value).reshape(-1)
    base_value = float(base_values[output_idx])
    pred_value = float(y_pred[output_idx])

    pretty_feature_names = [
        "Carbohydrates wt%" if f == "Carbohydrates" else
        "Proteins wt%" if f == "Proteins" else
        "Lipids wt%" if f == "Lipids" else
        "Ash wt%" if f == "Ash" else
        "Temperature (°C)" if f == "Temperature (C)" else
        "Residence Time (min)" if f == "Residence Time" else
        "Solid loading(%)" if f == "Solid content (w/w) %" else
        "Pre-processing" if f == "Pre-processing" else
        "Catalyst" if f == "Catalyst" else
        "Solvent" if f == "Solvent" else
        "Reactor Type" if f == "Reactor Type" else
        "Reactor Volume (L)" if f == "Reactor Volume (mL)" else
        "C%" if f == "C%" else
        "H%" if f == "H%" else
        "N%" if f == "N%" else
        "S%" if f == "S%" else
        "O%" if f == "O%" else
        f
        for f in FEATURE_COLS
    ]

    decoded_feature_values = [
        format_display_value(f, v)
        for f, v in zip(FEATURE_COLS, feature_values.values)
    ]

    custom_waterfall_path = os.path.join(
        output_dir, f"WATERFALL_{feedstock}_{output_name}_{timestamp}.png"
    )

    custom_waterfall_plot(
        shap_values=local_shap,
        feature_names=pretty_feature_names,
        feature_display_values=decoded_feature_values,
        base_value=base_value,
        final_value=pred_value,
        save_path=custom_waterfall_path,
        max_display=10,
    )

    print("✅ Custom waterfall plot saved")
    print(f"   → {custom_waterfall_path}")
    print(f"   Base value ({output_name}): {base_value:.4f}")
    print(f"   Prediction ({output_name}): {pred_value:.4f}")
    print(f"   Check base + shap sum: {(base_value + local_shap.sum()):.4f}")