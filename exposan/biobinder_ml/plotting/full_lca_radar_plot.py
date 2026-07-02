# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:41:48 2025

@author: aliah
"""

import numpy as np 
import matplotlib.pyplot as plt

categories = [
    "Global Warming", "Acidification", "Ecotoxicity", "Eutrophication",
    "Ozone Depletion", "Smog Formation", "Carcinogens",
    "Non-Carcinogens", "PM2.5 Formation"
]


data = {
    "c-HTL no EC": [-0.994, -0.302, -2.482, -3.380, -3.165, -0.307, -0.706, -2.145, -0.564],
    "d-HTL no EC": [-0.997, -0.237, -2.482, -3.383, -3.230, -0.030, -0.782, -2.149, -2.259],
    "c-HTL no LFC": [ 0.171, -0.004,  0.450,  1.609, -2.097, -0.099,  0.576,  0.539, -0.249],
    "AI (2019)":    [0]*9
}

angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
angles += angles[:1]  # close loop

plt.figure(figsize=(9, 9), dpi=300)
ax = plt.subplot(111, polar=True)

colors = ["navy", "royalblue", "darkorange", "red"]

# Plot scenarios
for (name, values), color in zip(data.items(), colors):
    vals = values + values[:1]  # close loop
    ax.plot(angles, vals, linewidth=2, color=color, label=name)
    ax.fill(angles, vals, color=color, alpha=0.15)

# Category labels
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontweight="bold", fontsize=11)

# --------- Custom radial labels ---------

ax.set_ylim(-4, 2)

radial_ticks = [-3, -2, -1, 0, 1, 2]
radial_labels = ["-3", "-2", "-1", "0", "1", "2"]

ax.set_yticks(radial_ticks)
ax.set_yticklabels(radial_labels, fontsize=13, fontweight="bold")
# ----------------------------------------
ax.yaxis.grid(True, linewidth=2)   # increase gridline thickness
# Legend
leg = plt.legend(
    fontsize=14,
    frameon=False,
    loc="upper right",
    bbox_to_anchor=(1.25, 1.10)
)

# Make scenario names bold
for text in leg.get_texts():
    text.set_fontweight("bold")

plt.tight_layout()
plt.show()
