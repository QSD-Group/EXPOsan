# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 17:55:30 2025

@author: aliah
"""

"""
Script: solve_tipping_for_fixed_IRR_price.py
Goal: For each Biobinder configuration, find the tipping fee ($/ton)
      that yields IRR = 10% and MSP = $0.10/kg.
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

# -----------------------------
# 1️⃣  Extended TEA class
# -----------------------------
class ExtendedTEA(BiobinderTEA):
    """Subclass that solves for tipping fee given target IRR and biobinder price."""

    def solve_tipping_fee(self, target_IRR, biobinder_price, price_dct, factor=1.0):
        biobinder = self.system.flowsheet.stream.biobinder
        biobinder.price = biobinder_price  # lock price (MSP target)

        def objective(t_fee):
            price_dct["tipping"] = t_fee / 1e3 * factor
            try:
                self.system.simulate()
                IRR = self.solve_IRR()
                return IRR - target_IRR
            except Exception:
                # Return a large residual if simulation fails
                return 1e3

        sol = fsolve(objective, x0=-40)  # initial guess = -40 $/ton
        return float(sol[0])

# -----------------------------
# 2️⃣  Configurations
# -----------------------------
configs = [
    {"name": "CHCU_No_EC", "config_kwargs": {
        "flowsheet": None,
        "central_dry_flowrate": default_central,
        "decentralized_HTL": False,
        "decentralized_upgrading": False,
        "skip_EC": True,
        "generate_H2": False,
        "EC_config": None,
    }},
    {"name": "CHCU_EC", "config_kwargs": {
        "flowsheet": None,
        "central_dry_flowrate": default_central,
        "decentralized_HTL": False,
        "decentralized_upgrading": False,
        "skip_EC": False,
        "generate_H2": False,
        "EC_config": None,
    }},
    {"name": "DHCU_No_EC", "config_kwargs": {
        "flowsheet": None,
        "central_dry_flowrate": default_central,
        "decentralized_HTL": True,
        "decentralized_upgrading": False,
        "skip_EC": True,
        "generate_H2": False,
        "EC_config": None,
    }},
    {"name": "DHCU_EC", "config_kwargs": {
        "flowsheet": None,
        "central_dry_flowrate": default_central,
        "decentralized_HTL": True,
        "decentralized_upgrading": False,
        "skip_EC": False,
        "generate_H2": False,
        "EC_config": None,
    }},
]

# -----------------------------
# 3️⃣  Main routine
# -----------------------------
def run_solver(target_IRR=0.10, biobinder_price=0.10, factor=1.0):
    results = []
    for cfg in configs:
        cfg_name = cfg["name"]
        print(f"\n🔧 Running configuration: {cfg_name}")

        try:
            # Fresh system each time
            sys = create_system(**cfg["config_kwargs"])
            base_tea = sys.TEA
            tea = ExtendedTEA.__new__(ExtendedTEA)
            tea.__dict__.update(base_tea.__dict__)
            tea.system = sys
            if not hasattr(tea, "_IRR"):
               tea.IRR = getattr(base_tea, "IRR", 0.10)
            if not hasattr(tea, "_NPV"):
               tea._NPV = None
            if not hasattr(tea, "_NPV_cache"):
               tea._NPV_cache = {}
            lca = sys.LCA
            biobinder = sys.flowsheet.stream.biobinder
            

            # Solve for tipping fee
            tipping_opt = tea.solve_tipping_fee(target_IRR, biobinder_price, price_dct, factor)

            # Recalculate key outputs
            price_dct["tipping"] = tipping_opt / 1e3 * factor
            sys.simulate()
            IRR = tea.solve_IRR()
            MSP = tea.solve_price(biobinder)
            GWP = lca.get_allocated_impacts(
                streams=(biobinder,), operation_only=True, annual=True
            )["GWP"] / (biobinder.F_mass * lca.system.operating_hours)

            print(f"✅ {cfg_name}: Tipping={tipping_opt:.2f} $/ton | IRR={IRR*100:.2f}% | MSP=${MSP:.2f}/kg | GWP={GWP:.3f}")

            results.append({
                "Configuration": cfg_name,
                "Tipping Fee ($/ton)": tipping_opt,
                "IRR (%)": IRR * 100,
                "MSP ($/kg)": MSP,
                "GWP (kg CO2e/kg)": GWP,
            })

        except Exception as e:
            print(f"❌ {cfg_name} failed: {e}")
            results.append({
                "Configuration": cfg_name,
                "Error": str(e),
            })

    return pd.DataFrame(results)

# -----------------------------
# 4️⃣  Execute & save results
# -----------------------------
if __name__ == "__main__":
    df = run_solver(target_IRR=0.10, biobinder_price=0.10)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, f"TippingFee_FixedIRR_MSP_{timestamp}.xlsx")
    df.to_excel(file_path, index=False, float_format="%.6f")
    print(f"\n📁 Results saved to: {file_path}")
