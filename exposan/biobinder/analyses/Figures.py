# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 12:30:04 2025

@author: aliah
"""

import qsdsan as qs
from chaospy import distributions as shape
import numpy as np
import pandas as pd
import os
from exposan.biobinder import (
    create_system,
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
)

# Create results directory
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# Define parameter values
irr_values = np.arange(0.0, 0.21, 0.01)
electricity_prices = np.arange(0.00, 0.11, 0.01)  # 0.00 to 0.10 for EC configs
electrode_costs = np.arange(0, 50001, 1000)  # 0 to 50,000 for EC configs

# List of configurations
configs = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": {}}},  # EC config enabled
    {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": {}}},  # EC config enabled
]

# Loop over each configuration
for config in configs:
    config_name = config["name"]
    config_kwargs = config["config_kwargs"]
    is_EC = config_kwargs["EC_config"] is not None  # Check if EC is enabled
    
    print(f"Running parameter sensitivity analysis for {config_name}...")

    # Run the system **once** for this configuration
    sys = create_system(**config_kwargs)
    sys.simulate()

    # Access TEA and Biobinder stream
    tea = sys.TEA
    biobinder = sys.flowsheet.stream.biobinder

    
    irr_results = []
    for irr in irr_values:
        tea.IRR = irr  # Set new IRR value
        MSP = tea.solve_price(biobinder)  # Compute MSP
        irr_results.append({'IRR': irr, 'MSP ($/kg)': MSP})
        print(f"[{config_name}] IRR: {irr:.2f}, MSP: ${MSP:.2f}/kg")

    irr_df = pd.DataFrame(irr_results)
    irr_file = os.path.join(output_dir, f"IRR_vs_MSP_{config_name}.xlsx")
    irr_df.to_excel(irr_file, index=False)
    print(f"Saved IRR results to {irr_file}")

 
    if is_EC:
        elec_results = []
        for elec_price in electricity_prices:
            qs.PowerUtility.price = elec_price  # Update electricity price
            MSP = tea.solve_price(biobinder)  # Compute MSP
            elec_results.append({'Electricity Price ($/kWh)': elec_price, 'MSP ($/kg)': MSP})
            print(f"[{config_name}] Electricity Price: {elec_price:.2f}, MSP: ${MSP:.2f}/kg")

        elec_df = pd.DataFrame(elec_results)
        elec_file = os.path.join(output_dir, f"Electricity_vs_MSP_{config_name}.xlsx")
        elec_df.to_excel(elec_file, index=False)
        print(f"Saved Electricity Price results to {elec_file}")

    
    if is_EC:
        electrode_results = []
        for cost in electrode_costs:
            # **For CHCU → Only Upgrading_EC**
            if "CHCU" in config_name:
                sys.flowsheet.unit.Upgrading_EC.electrode_cost = cost
            # **For DHCU → Both HTL_EC & Upgrading_EC**
            elif "DHCU" in config_name:
                sys.flowsheet.unit.HTL_EC.electrode_cost = cost
                sys.flowsheet.unit.Upgrading_EC.electrode_cost = cost

            MSP = tea.solve_price(biobinder)  # Compute MSP
            electrode_results.append({'Electrode Cost ($/m²)': cost, 'MSP ($/kg)': MSP})
            print(f"[{config_name}] Electrode Cost: {cost}, MSP: ${MSP:.2f}/kg")

        electrode_df = pd.DataFrame(electrode_results)
        electrode_file = os.path.join(output_dir, f"Electrode_vs_MSP_{config_name}.xlsx")
        electrode_df.to_excel(electrode_file, index=False)
        print(f"Saved Electrode Cost results to {electrode_file}")

print("Sensitivity analysis completed for all configurations.")