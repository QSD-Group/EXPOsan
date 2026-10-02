# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 17:48:19 2025

@author: aliah
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

# Output directory
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# === Configuration list ===
configs = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
        "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
        "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
]

# === Solver Function ===
def solve_tipping_for_IRR_MSP(config, target_IRR=0.5, target_MSP=0.10, factor=1.0):
    sys = create_system(**config["config_kwargs"])
    tea = sys.TEA
    lca = sys.LCA
    biobinder = sys.flowsheet.stream.biobinder

    def equations(x):
        tipping_fee = x[0]
        price_dct["tipping"] = tipping_fee / 1e3 * factor
        sys.simulate()

        IRR = tea.solve_IRR()
        MSP = tea.solve_price(biobinder)
        return [IRR - target_IRR, MSP - target_MSP]

    # Initial guess (start from moderate tipping fee)
    guess = [-40]
    sol = fsolve(equations, guess)
    tipping_opt = float(sol[0])

    # Recalculate impacts at solution
    price_dct["tipping"] = tipping_opt / 1e3 * factor
    sys.simulate()
    IRR = tea.solve_IRR()
    MSP = tea.solve_price(biobinder)
    GWP = lca.get_allocated_impacts(
        streams=(biobinder,), operation_only=True, annual=True
    )["GWP"] / (biobinder.F_mass * lca.system.operating_hours)

    print(f"✅ {config['name']}: Tipping={tipping_opt:.2f} $/ton | IRR={IRR*100:.2f}% | MSP=${MSP:.2f}/kg | GWP={GWP:.3f}")
    return {"Configuration": config["name"], "Tipping Fee ($/ton)": tipping_opt,
            "IRR (%)": IRR*100, "MSP ($/kg)": MSP, "GWP (kg CO2e/kg)": GWP}


# === Run Solver for Each Configuration ===
results = []
for cfg in configs:
    try:
        res = solve_tipping_for_IRR_MSP(cfg, target_IRR=0.10, target_MSP=0.20)
        results.append(res)
    except Exception as e:
        print(f"❌ {cfg['name']} failed: {e}")
        results.append({"Configuration": cfg["name"], "Error": str(e)})

# === Save Results ===
df = pd.DataFrame(results)
file_path = os.path.join(output_dir, f"Tipping_Solver_IRR_MSP_{timestamp}.xlsx")
df.to_excel(file_path, index=False)
print(f"📁 Results saved to: {file_path}")
