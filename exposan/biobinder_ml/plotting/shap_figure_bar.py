# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 15:52:58 2026

@author: aliah
"""

import os, pandas as pd, numpy as np
import matplotlib.pyplot as plt

in_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\SUMMARY_KNOBNORM_WIDE.xlsx"
out_png = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\KnobNorm_Yield_vs_TEA_StackedBars_300dpi.png"
out_pdf = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\KnobNorm_Yield_vs_TEA_StackedBars.pdf"

# ----------------------------
# Load wide tables
# ----------------------------
yield_wide = pd.read_excel(in_path, sheet_name="YIELD_KNOBNORM", index_col=0)
tea_wide   = pd.read_excel(in_path, sheet_name="MSP_KNOBNORM", index_col=0)

# Align
yield_wide.index = yield_wide.index.astype(str)
tea_wide.index   = tea_wide.index.astype(str)

features = [f for f in yield_wide.index if f in tea_wide.index]
yield_wide = yield_wide.loc[features]
tea_wide   = tea_wide.loc[features]

feedstocks = [c for c in yield_wide.columns if c in tea_wide.columns]
yield_wide = yield_wide[feedstocks]
tea_wide   = tea_wide[feedstocks]

# Rank feedstocks by TEA concentration = max feature share (optional but useful)
feedstocks_sorted = list(tea_wide.max(axis=0).sort_values(ascending=False).index)
yield_wide = yield_wide[feedstocks_sorted]
tea_wide   = tea_wide[feedstocks_sorted]

# Numeric
yield_wide = yield_wide.apply(pd.to_numeric, errors="coerce").fillna(0.0)
tea_wide   = tea_wide.apply(pd.to_numeric, errors="coerce").fillna(0.0)

# ----------------------------
# Plot geometry
# ----------------------------
n = len(feedstocks_sorted)
x = np.arange(n)* 0.85
bar_w = 0.22

# Colors (7 features → 7 colors)
cmap = plt.get_cmap("tab20")
colors = [cmap(i) for i in range(0, 14, 2)]  # 7 distinct colors

# IMPORTANT: make width depend on n so it doesn't squash
fig, ax = plt.subplots(figsize=(max(12, n * 1.8), 7.2))

bottom_y = np.zeros(n)
bottom_t = np.zeros(n)

# ----------------------------
# Stacked bars
# ----------------------------
for i, feat in enumerate(features):
    y_vals = yield_wide.loc[feat].values
    t_vals = tea_wide.loc[feat].values

    # Label only the first set so we can auto-collect handles if needed
    ax.bar(x - bar_w/2, y_vals, bar_w, bottom=bottom_y, color=colors[i], label=feat)
    ax.bar(x + bar_w/2, t_vals, bar_w, bottom=bottom_t, color=colors[i])

    bottom_y += y_vals
    bottom_t += t_vals

# ----------------------------
# Axis formatting
# ----------------------------
ax.set_xticks(x)

# If you want to enforce specific feedstock text styling, you can map here:
# (tick labels are bold anyway)
label_map = {fs: fs for fs in feedstocks_sorted}
# e.g., capitalize if you want: label_map["fog"] = "fog"; label_map["sludge"] = "sludge"
ax.set_xticklabels([label_map[fs] for fs in feedstocks_sorted], fontsize=12, fontweight="bold")

ax.set_ylabel("Normalized SHAP share (%)", fontsize=12, fontweight="bold")

ax.tick_params(axis='y', labelsize=11)
for lab in ax.get_yticklabels():
    lab.set_fontweight('bold')

# Put Yield / MSP under each pair (and make them bold)
for xi in x:
    ax.text(xi - bar_w/2, -2.5, "Yield", ha="center", va="top",
            fontsize=12, fontstyle="italic", fontweight="bold", clip_on=False)
    ax.text(xi + bar_w/2, -2.5, "MSP", ha="center", va="top",
            fontsize=12, fontstyle="italic", fontweight="bold", clip_on=False)

ax.set_ylim(-6, 100)
ax.grid(axis="y", alpha=0.25)

# ----------------------------
# Legend at bottom (FIGURE legend = no axis shrink)
# ----------------------------
handles, labels = ax.get_legend_handles_labels()
# Keep only unique labels in order
seen = set()
uniq_handles, uniq_labels = [], []
for h, l in zip(handles, labels):
    if l not in seen:
        uniq_handles.append(h)
        uniq_labels.append(l)
        seen.add(l)

legend = fig.legend(
    uniq_handles, uniq_labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.08),   # inside the figure bottom margin
    ncol=len(features),           # one row
    frameon=False,
    fontsize=9
)
plt.setp(legend.get_texts(), fontweight="bold")

# Reserve space for bottom legend and bottom "Yield/MSP" texts
# (This is the key to avoid skinny plot)
fig.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.20)

# ----------------------------
# Save
# ----------------------------
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf, bbox_inches="tight")

print("Saved:", out_png)
print("Saved:", out_pdf)
