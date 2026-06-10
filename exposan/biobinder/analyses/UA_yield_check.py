# -*- coding: utf-8 -*-
"""
Created on Wed May 21 15:36:09 2025

@author: aliah
"""

import qsdsan as qs
from chaospy import distributions as shape
import numpy as np
import pandas as pd
from scipy.stats import qmc
import os
from datetime import datetime
from exposan.biobinder import (
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
)
from exposan.biobinder.analyses.Dist_flex import (
    create_system)

# Create a directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M")
# Define the function to create the model

class BiobinderModel:
    def __init__(self, config_kwargs, param_distributions=None):
        self.config_kwargs = config_kwargs
        self.param_distributions = param_distributions or {}

    def add_parameter(self, name, element_str, attr, kind, units, baseline, distribution):
        self.param_distributions[name] = {
            "element_str": element_str,  # e.g., "stream.biofuel" or "unit.HTL"
            "attr": attr,
            "kind": kind,
            "units": units,
            "baseline": baseline,
            "distribution": distribution,
        }

    def generate_samples(self, n_samples):
        """
        Generate Latin Hypercube Samples (LHS) for uncertainty analysis.
        """
        num_params = len(self.param_distributions)
        sampler = qmc.LatinHypercube(d=num_params, seed=42)  # d = number of parameters
        lhs_samples = sampler.random(n=n_samples)  # Generates normalized LHS samples (0,1)

        # Convert normalized LHS samples into actual parameter values
        samples = []
        param_names = list(self.param_distributions.keys())

        for i in range(n_samples):
              sample = {}
              for j, param in enumerate(param_names):
                  dist = self.param_distributions[param]["distribution"]
                  scaled_value = dist.ppf(lhs_samples[i, j])
                  sample[param] = float(scaled_value)
              samples.append(sample)

        return samples

    def simulate_with_sample(self, sample):
        try:
            sys = create_system(**self.config_kwargs)
            tea = sys.TEA
            lca = sys.LCA
            biobinder = sys.flowsheet.stream.biobinder
            tea.IRR = tea_kwargs["IRR"]
            tea.income_tax = tea_kwargs["income_tax"]

            for name, value in sample.items():
                element_str = self.param_distributions[name]["element_str"]
                attr = self.param_distributions[name]["attr"]
                element = eval(f"sys.flowsheet.{element_str}")
                setattr(element, attr, value)

            sys.simulate()
            MSP = tea.solve_price(biobinder)
            GWP = lca.get_allocated_impacts(
                streams=(biobinder,), operation_only=True, annual=True
            )["GWP"] / (biobinder.F_mass * lca.system.operating_hours)

            return {**sample, "MSP ($/kg)": MSP, "GWP (kg CO2e/kg)": GWP}

        except Exception as e:
            print(f"Simulation failed: {e}")
            return {**sample, "MSP ($/kg)": np.nan, "GWP (kg CO2e/kg)": np.nan, "Error": str(e)}

    def run_uncertainty_analysis(self, n_samples):
        print(f"Running uncertainty analysis with {n_samples} samples...")
        results = []
        samples = self.generate_samples(n_samples)
        for i, sample in enumerate(samples, 1):
            print(f"▶️ Sample {i} of {n_samples}")
            results.append(self.simulate_with_sample(sample))
        return pd.DataFrame(results)

if __name__ == "__main__":        
    # Configuration definitions
     configs = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
     ]

# Biocrude yield loop
biocrude_yields = np.round(np.arange(0.30, 0.70 + 0.02, 0.02), 2) # 0.6680000000000004, 0.5690000000000002, 0.7000000000000004 works
gas_yield = 0.1756
char_yield = 0.0100  
n_samples = 2

valid_yields = []
Y_gas = 0.1756
Y_char = 0.01

for Y_bio in np.arange(0.30, 0.71, 0.0001):
    Y_aq = 1 - (Y_bio + Y_gas + Y_char)
    if Y_aq < 0 or Y_aq > 0.55:
        continue  # skip invalid yields

    HTL_yields.update({
        'biocrude': Y_bio,
        'gas': Y_gas,
        'char': Y_char,
        'aqueous': Y_aq
    })

    try:
        sys = create_system(**configs[0]["config_kwargs"])  # Pick one config to test
        sys.simulate()
        print(f"✅ Y_bio={Y_bio:.2f} succeeded.")
        valid_yields.append(Y_bio)
    except Exception as e:
        print(f"❌ Y_bio={Y_bio:.2f} failed: {e}")


    

    


