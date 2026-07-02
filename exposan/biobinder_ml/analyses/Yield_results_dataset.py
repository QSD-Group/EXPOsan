# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 17:05:41 2025

@author: aliah
"""

import qsdsan as qs
from chaospy import distributions as shape
import numpy as np
import pandas as pd
from scipy.stats import qmc
import os
from datetime import datetime
from exposan.biobinder_ml import (
    # create_system,
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    BiobinderTEA,
    central_dry_flowrate as default_central,
    data_path,
    feedstock_composition,
    HTL_yields,
    pilot_dry_flowrate as default_pilot,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    )
from exposan.biobinder_ml.Dist_flex import (
    create_system)

# Create a directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

yield_table_path = os.path.join(output_dir, "HTL_all_yield_experimental.xlsx")  # adjust if needed
yt = pd.read_excel(yield_table_path, sheet_name="Sheet1")
yt.columns = yt.columns.str.strip()

required_cols = [ "Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]
missing_req = [c for c in required_cols if c not in yt.columns]
if missing_req:
    raise ValueError(f"Missing required columns in HTL_yield_type.xlsx: {missing_req}")

# Drop rows fully empty in yield columns, convert wt% -> fractions
yt = yt.dropna(subset=["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"], how="all").copy()
for col in ["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]:
    yt[col] = pd.to_numeric(yt[col], errors="coerce") / 100.0

# Sum check with ±5% tolerance (no normalization; skip invalid rows)
yt["Sum_fractions"] = yt[["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]].sum(axis=1)
tolerance = 0.05
valid_mask = (yt["Sum_fractions"] >= 1 - tolerance) & (yt["Sum_fractions"] <= 1 + tolerance)
valid_rows = yt.loc[valid_mask].copy()
skipped_rows = yt.loc[~valid_mask].copy()

# Prepare fraction columns for clarity
valid_rows["Solids_f"]   = valid_rows["Solids wt%"]
valid_rows["Biocrude_f"] = valid_rows["Biocrude wt%"]
valid_rows["Aqueous_f"]  = valid_rows["Aqueous wt%"]
valid_rows["Gas_f"]      = valid_rows["Gas wt%"]

# Assign unique index for traceability
valid_rows = valid_rows.reset_index(drop=True)
valid_rows["Sample_ID"] = valid_rows.index + 1  # start from 1 for readability

timestamp = datetime.now().strftime("%Y%m%d_%H%M")

if __name__ == "__main__":        
    # Configuration definitions
     configs = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    # {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
    # {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    # {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
      ]

for config in configs:
    config_name = config["name"]
    print(f"\n🔧 Running for configuration: {config_name}")

    all_results = []  # accumulate ALL feedstock rows for this config

    for idx, row in valid_rows.iterrows():
        sample_id = int(row["Sample_ID"])  # new tracking ID
        Y_char = float(row["Solids_f"])   # solids=char
        Y_bio  = float(row["Biocrude_f"])
        Y_aq   = float(row["Aqueous_f"])
        Y_gas  = float(row["Gas_f"])

        # Update HTL yields
        HTL_yields.update({
            "biocrude": Y_bio,
            "gas":      Y_gas,
            "char":     Y_char,
            "aqueous":  Y_aq,
        })

        try:
            sys = create_system(**config["config_kwargs"])
            tea = sys.TEA
            lca = sys.LCA
            tea.IRR = tea_kwargs["IRR"]
            tea.income_tax = tea_kwargs["income_tax"]

            sys.simulate()
            biobinder = sys.flowsheet.stream.biobinder
            biofuel   = sys.flowsheet.stream.biofuel
            Column    = sys.flowsheet.unit.CrudeHeavyDis
            Cutoff    = sys.flowsheet.unit.BiocrudeSplitter._cutoff_fracs

            MSP  = tea.solve_price(biobinder)
            IRR  = tea.solve_IRR() * 100
            GWP  = lca.get_allocated_impacts(
                streams=(biobinder,), operation_only=True, annual=True
            )["GWP"] / (biobinder.F_mass * lca.system.operating_hours)
            ratio = biobinder.F_mass / biofuel.F_mass if biofuel.F_mass else float('nan')

            all_results.append({
                "Sample_ID": sample_id,
                "Solids (char) [frac]": Y_char,
                "Biocrude [frac]": Y_bio,
                "Aqueous [frac]": Y_aq,
                "Gas [frac]": Y_gas,
                "Solids wt%": Y_char * 100,
                "Biocrude wt%": Y_bio * 100,
                "Aqueous wt%": Y_aq * 100,
                "Gas wt%": Y_gas * 100,
                "Sum wt%": (Y_char + Y_bio + Y_aq + Y_gas) * 100,
                "MSP ($/kg)": MSP,
                "IRR (%)": IRR,
                "GWP (kg CO2e/kg)": GWP,
                "Product Ratio (biobinder/biofuel)": ratio,
                "Cutoff Fracs": Cutoff,
                # "LHK": Column._LHK,  
            })
            # print(f"    ✅ ftype={ftype} | idx={idx} | MSP=${MSP:.2f} | IRR={IRR:.2f}% | GWP={GWP:.4f}")

        except Exception as e:
            # print(f"    ❌ ftype={ftype} | idx={idx} failed: {e}")
            all_results.append({
                "Sample_ID": sample_id,
                "Solids (char) [frac]": Y_char,
                "Biocrude [frac]": Y_bio,
                "Aqueous [frac]": Y_aq,
                "Gas [frac]": Y_gas,
                "MSP ($/kg)": np.nan,
                "IRR (%)": np.nan,
                "GWP (kg CO2e/kg)": np.nan,
                "Error": str(e),
            })

    # Save ONE file per config, with results + skipped rows as separate sheets
    out_name = f"YieldsFromExcel_{config_name}_{timestamp}.xlsx"
    out_path = os.path.join(output_dir, out_name)
    with pd.ExcelWriter(out_path) as w:
        pd.DataFrame(all_results).to_excel(w, sheet_name="results", index=False, float_format="%.15f")
        if not skipped_rows.empty:
            # Save skipped exactly as in input units (wt%), plus the sum column (in fractions)
            skip_out = skipped_rows.copy()
            # convert back to wt% for readability
            for col in ["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]:
                skip_out[col] = skip_out[col] * 100.0
            skip_out.rename(columns={"Sum_fractions": "Sum (fraction)"}, inplace=True)
            skip_out.to_excel(w, sheet_name="skipped_rows", index=False)

    print(f"  📁 Saved: {out_path}")
    if not skipped_rows.empty:
        print(f"  ⚠️ Skipped rows included in sheet: 'skipped_rows'")