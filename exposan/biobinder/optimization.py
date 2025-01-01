# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 18:22:30 2024

@author: aliah
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from exposan.biobinder.systems import create_system, simulate_and_print


output_dir = "simulation_results"
os.makedirs(output_dir, exist_ok=True)

plant_sizes = np.arange(1, 151, 0.5)  # 1 to 150 tpd

EC_future_config = {
    'EO_voltage': 2.5,  # Optimistic assumptions
    'ED_voltage': 2.5,
    'electrode_cost': 225,
    'anion_exchange_membrane_cost': 0,
    'cation_exchange_membrane_cost': 0,
    'electricity_price': 0.03,
    'electricity_GHG': 0,
}

configurations = [
    {"name": "CHCU_no_EC", "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None},
    {"name": "CHCU_EC", "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None},
    {"name": "CHCU_EC_H2", "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": True, "EC_config": None},
    {"name": "CHCU_EC_Futures", "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": True, "EC_config": EC_future_config},
    {"name": "DHCU_no_EC", "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None},
    {"name": "DHCU_EC", "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None},
    {"name": "DHCU_EC_H2", "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": True, "EC_config": None},
    {"name": "DHCU_EC_Futures", "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": True, "EC_config": EC_future_config},
]

simulation_results = []
optimization_results = []


def optimize_plant_size(config):
    def objective(plant_size):
        try:
            config_kwargs = dict(
                flowsheet=None,
                central_dry_flowrate=plant_size * 907.185 / (24 * 0.9),
                pilot_dry_flowrate=None,
                **{k: v for k, v in config.items() if k != "name"},
            )
            sys = create_system(**config_kwargs)
            sys.simulate()
            MSP = sys.TEA.solve_price(sys.flowsheet.stream.biobinder)
            return abs(MSP - 0.10)  # Objective: minimize the absolute difference from 0.10 $/kg
        except Exception as e:
            print(f"Optimization failed for plant size {plant_size} tpd: {e}")
            return float('inf')  # Return a large value if simulation fails

    result = minimize_scalar(
        objective, bounds=(10, 150), method='bounded', options={'xatol': 1e-2}
    )
    return result


for config in configurations:
    config_name = config["name"]
    
    for plant_size in plant_sizes:
        try:
            config_kwargs = dict(
                flowsheet=None,
                central_dry_flowrate=plant_size * 907.185 / (24 * 0.9),
                pilot_dry_flowrate=None,
                **{k: v for k, v in config.items() if k != "name"},
            )
            sys = create_system(**config_kwargs)
            sys.simulate()
            MSP = sys.TEA.solve_price(sys.flowsheet.stream.biobinder)
            all_impacts = sys.LCA.get_allocated_impacts(streams=(sys.flowsheet.stream.biobinder,), operation_only=True, annual=True)
            GWP = all_impacts['GWP'] / (sys.flowsheet.stream.biobinder.F_mass * sys.LCA.system.operating_hours)

            simulation_results.append({
                "Configuration": config_name,
                "Plant Size (tpd)": plant_size,
                "MSP ($/kg)": MSP,
                "GWP (kg CO2e/kg)": GWP
            })
        except Exception as e:
            print(f"Simulation failed for plant size {plant_size} tpd, config {config_name}: {e}")

    # Optimization only for "CHCU_no_EC" and "CHCU_EC_Futures"
    if config_name in ["CHCU_no_EC", "CHCU_EC_Futures"]:
        result = optimize_plant_size(config)
        optimal_plant_size = result.x if result.success else None
        optimal_msp = result.fun + 0.10 if result.success else None

        optimization_results.append({
            "Configuration": config_name,
            "Optimal Plant Size (tpd)": optimal_plant_size,
            "Optimal MSP ($/kg)": optimal_msp
        })

simulation_df = pd.DataFrame(simulation_results)
simulation_file = os.path.join(output_dir, "simulation_results.xlsx")
simulation_df.to_excel(simulation_file, index=False)

# Save optimization results to Excel
optimization_df = pd.DataFrame(optimization_results)
optimization_file = os.path.join(output_dir, "optimization_results.xlsx")
optimization_df.to_excel(optimization_file, index=False)

print("Simulation and optimization completed.")
