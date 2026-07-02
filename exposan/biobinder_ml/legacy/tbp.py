# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 11:56:41 2025

@author: aliah
"""

#Zhnag Group Data
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

output_dir = "results"
os.makedirs(output_dir, exist_ok=True)
data = [
    {
        "Sample name": "2022 SDW",
        "C_wt%": 74.2, "H_wt%": 11.21, "N_wt%": 0.75, "S_wt%": 0.01, "O_wt%": 13.84,
        "HHV_MJkg": 37.62,
        "Gasoline_<150C_%": 0.9141, "Kerosene_150-230C_%": 6.079, "Diesel_230-350C_%": 46.17,
        "Lube_350-375C_%": 7.392, "HFO_375-550C_%": 37.26, "Residue_>550C_%": 2.185,
        "Visc_cP": 654.33, "Density_gml": 0.9278, "TAN_mgKOHg": 132.04
    },
    {
        "Sample name": "2022 HMFW",
        "C_wt%": 74.12, "H_wt%": 10.65, "N_wt%": 3.16, "S_wt%": 0.12, "O_wt%": 11.7,
        "HHV_MJkg": 37.34,
        "Gasoline_<150C_%": 1.704, "Kerosene_150-230C_%": 7.54, "Diesel_230-350C_%": 52.46,
        "Lube_350-375C_%": 9.135, "HFO_375-550C_%": 22.69, "Residue_>550C_%": 6.471,
        "Visc_cP": 4448.33, "Density_gml": 0.9456, "TAN_mgKOHg": 140.66
    },
    {
        "Sample name": "2023 SDW",
        "C_wt%": 75.57, "H_wt%": 11.31, "N_wt%": 0.64, "S_wt%": 0.04, "O_wt%": 12.45,
        "HHV_MJkg": 38.37,
        "Gasoline_<150C_%": 0.56, "Kerosene_150-230C_%": 8.37, "Diesel_230-350C_%": 66.74,
        "Lube_350-375C_%": 3.9, "HFO_375-550C_%": 17.24, "Residue_>550C_%": 3.16,
        "Visc_cP": 142.67, "Density_gml": 0.926, "TAN_mgKOHg": 185.17
    },
    {
        "Sample name": "2023 HMFW",
        "C_wt%": 75.24, "H_wt%": 10.64, "N_wt%": 2.12, "S_wt%": 0.01, "O_wt%": 12.0,
        "HHV_MJkg": 37.62,
        "Gasoline_<150C_%": 0.73, "Kerosene_150-230C_%": 8.44, "Diesel_230-350C_%": 65.85,
        "Lube_350-375C_%": 4.59, "HFO_375-550C_%": 12.68, "Residue_>550C_%": 7.71,
        "Visc_cP": 3306.33, "Density_gml": 0.977, "TAN_mgKOHg": 242.74
    },
    {
        "Sample name": "2024 SDW",
        "C_wt%": 76.74, "H_wt%": 11.15, "N_wt%": 1.4, "S_wt%": 0.02, "O_wt%": 10.68,
        "HHV_MJkg": 38.85,
        "Gasoline_<150C_%": 0.82, "Kerosene_150-230C_%": 8.3, "Diesel_230-350C_%": 73.88,
        "Lube_350-375C_%": 3.34, "HFO_375-550C_%": 9.05, "Residue_>550C_%": 4.6,
        "Visc_cP": np.nan, "Density_gml": 0.93, "TAN_mgKOHg": np.nan
    },
    {
        "Sample name": "2024 HMFW",
        "C_wt%": 75.56, "H_wt%": 10.27, "N_wt%": 3.3, "S_wt%": 0.1, "O_wt%": 10.78,
        "HHV_MJkg": 37.51,
        "Gasoline_<150C_%": 1.09, "Kerosene_150-230C_%": 6.64, "Diesel_230-350C_%": 60.06,
        "Lube_350-375C_%": 6.11, "HFO_375-550C_%": 16.17, "Residue_>550C_%": 9.94,
        "Visc_cP": 2078.0, "Density_gml": 0.94, "TAN_mgKOHg": 255.97
    },
]

df = pd.DataFrame(data)

# Physics features
df["H_over_C_atomic"] = (12*df["H_wt%"])/df["C_wt%"]
df["O_over_C_atomic"] = (3*df["O_wt%"])/(4*df["C_wt%"])
df["S_over_C_atomic"] = (3*df["S_wt%"])/(8*df["C_wt%"])

# API (approx; density likely near 15-20C; treat as SG)
df["API"] = 141.5/df["Density_gml"] - 131.5

# Simple Kw estimate using MeABP proxy from mid-boiling of each cut (very rough):
# We'll compute a volume-weighted "mean boiling temperature" from the cut shares and use Kw = (Tb_mean[K])^(1/3)/SG
# Define representative midpoints (°C) for each interval
midpoints = {
    "Gasoline_<150C_%": 100.0,
    "Kerosene_150-230C_%": 190.0,
    "Diesel_230-350C_%": 290.0,
    "Lube_350-375C_%": 362.5,
    "HFO_375-550C_%": 462.5,
    "Residue_>550C_%": 600.0,  # placeholder for >550 tail
}
for k in midpoints:
    df[f"{k}_TmidC"] = midpoints[k]

# Compute mean boiling temp (volume-weighted by TGA shares, assuming vol% ≈ wt%)
shares_cols = [c for c in df.columns if c.endswith("_%") and ("Gasoline" in c or "Kerosene" in c or "Diesel" in c or "Lube" in c or "HFO" in c or "Residue" in c)]
wts = df[shares_cols].values/100.0
Tmids = np.vstack([midpoints[c] for c in shares_cols]).T  # shape (M cuts)-> broadcast later
Tb_mean_C = (wts * Tmids).sum(axis=1)
df["Tb_mean_C"] = Tb_mean_C
df["Tb_mean_K"] = df["Tb_mean_C"] + 273.15
df["Kw_est"] = (df["Tb_mean_K"]**(1/3))/df["Density_gml"]


boundaries = [150, 230, 350, 375, 550, 650]  # use 650C as an upper cap to close CDF
cut_order = ["Gasoline_<150C_%","Kerosene_150-230C_%","Diesel_230-350C_%","Lube_350-375C_%","HFO_375-550C_%","Residue_>550C_%"]

cum = np.cumsum(df[cut_order].values, axis=1)  # cumulative wt% ≈ vol%
TBP_points = pd.DataFrame(cum, columns=[f"CumVol%@{T}C" for T in boundaries])
TBP_points.insert(0, "Sample name", df["Sample name"])

# Merge physics features for display
feat_cols = ["Sample name","C_wt%","H_wt%","O_wt%","S_wt%","H_over_C_atomic","O_over_C_atomic","S_over_C_atomic","Density_gml","API","Tb_mean_C","Kw_est"]
features = df[feat_cols]

# Save CSVs
features_path = os.path.join(output_dir, f"biocrude_CHOS_physics_features.csv")
tbp_path = os.path.join(output_dir, f"biocrude_approx_TBP_from_TGA.csv")
features.to_csv(features_path, index=False)
TBP_points.to_csv(tbp_path, index=False)


plt.figure(figsize=(7,5))
for i, row in TBP_points.iterrows():
    temps = boundaries
    cumv = row[[f"CumVol%@{T}C" for T in boundaries]].values
    plt.plot(temps, cumv, marker='o', label=row["Sample name"])
plt.xlabel("Temperature (°C)")
plt.ylabel("Cumulative wt% (approx)")
plt.title("Approximate TBP from TGA cuts")
plt.legend()
plot_path = os.path.join(output_dir, f"approx_TBP_from_TGA.png")
plt.tight_layout()
plt.savefig(plot_path)
plot_path
