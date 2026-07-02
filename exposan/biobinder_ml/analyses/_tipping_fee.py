# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 12:59:43 2025

@author: aliah
"""

"""
scan_tipping_IRR.py
-----------------------------------------------------------
Scan tipping fees ($/ton) for each Biobinder configuration
to find which one yields IRR ≈ 10% at MSP = $0.10/kg.
-----------------------------------------------------------
"""
from datetime import datetime
import numpy as np
import pandas as pd
import os, biosteam as bst, qsdsan as qs
from qsdsan import sanunits as qsu
from scipy.optimize import fsolve
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea
from exposan.biobinder_ml import (
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    BiobinderTEA,
    cutoff_fracs_from_feed,
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
_psi_to_Pa = 6894.76
SnowdenSwan_factor= 1.0646188596312727
configure_utilities()

# ---------------------------------------------------------------------
#  Function: scan tipping fees
# ---------------------------------------------------------------------
def scan_tipping_for_IRR(config_kwargs, factor=1.0,
                         biobinder_price=0.10, target_IRR=0.10,
                         start=-40, stop=0, step=1):
    """Build fresh system each step, record IRR for each tipping fee."""
    results = []

    for tipping in np.arange(start, stop + step, step):
        try:
            # Fresh system for every run
            sys = create_system(**config_kwargs)
            tea, biobinder = sys.TEA, sys.flowsheet.stream.biobinder

            price_dct["tipping"] = tipping / 1e3 * factor
            biobinder.price = biobinder_price

            sys.simulate()
            IRR = tea.solve_IRR()
            results.append((tipping, IRR * 100))
            print(f"Tipping = {tipping:>6.1f} $/ton → IRR = {IRR*100:5.2f}%")

        except Exception as e:
            print(f"Tipping {tipping} failed: {e}")
            results.append((tipping, np.nan))

    df = pd.DataFrame(results, columns=["Tipping ($/ton)", "IRR (%)"])

    # Interpolate tipping that gives target IRR (10%)
    df_valid = df.dropna()
    tipping_interp = np.nan
    if len(df_valid) >= 2:
        tipping_interp = np.interp(
            target_IRR * 100,
            df_valid["IRR (%)"],
            df_valid["Tipping ($/ton)"]
        )

    return tipping_interp, df

# ---------------------------------------------------------------------
#  Configurations
# ---------------------------------------------------------------------
configs = [
{"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
# {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
# {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
# # {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
 ]

# ---------------------------------------------------------------------
#  Main routine
# ---------------------------------------------------------------------
if __name__ == "__main__":
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    results_all = []

    for cfg in configs:
        name = cfg["name"]
        print(f"\n🔧 Running configuration: {name}")

        try:
            tipping_est, df = scan_tipping_for_IRR(
                cfg["config_kwargs"],
                biobinder_price=0.10,
                target_IRR=0.10,
                start=-40,
                stop=0,
                step=1
            )

            # Append the interpolated row after interpolation only
            if not np.isnan(tipping_est):
                df.loc[len(df)] = [tipping_est, 10.0]

            df["Configuration"] = name
            path = os.path.join(output_dir, f"IRR_TippingScan_{name}_{timestamp}.xlsx")
            df.to_excel(path, index=False)
            print(f"📁 Saved full scan to {path}")
            print(f"≈ Tipping fee for 10% IRR: {tipping_est:.2f} $/ton")

            results_all.append({
                "Configuration": name,
                "Tipping Fee for 10% IRR ($/ton)": tipping_est
            })

        except Exception as e:
            print(f"❌ {name} failed: {e}")
            results_all.append({
                "Configuration": name,
                "Error": str(e)
            })

    summary = pd.DataFrame(results_all)
    summary_path = os.path.join(output_dir, f"TippingSummary_{timestamp}.xlsx")
    summary.to_excel(summary_path, index=False)
    print(f"\n✅ Summary saved to: {summary_path}")