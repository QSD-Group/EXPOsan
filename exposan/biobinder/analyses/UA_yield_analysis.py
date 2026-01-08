# -*- coding: utf-8 -*-
"""
Created on Thu May 22 11:58:13 2025

@author: aliah
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 20 16:19:42 2025

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

# Create a directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

CHCU_yield_df = pd.read_excel(os.path.join(output_dir, "CHCU_Ybio_values.xlsx"))
DHCU_yield_df = pd.read_excel(os.path.join(output_dir, "DHCU_Ybio_values.xlsx"))

timestamp = datetime.now().strftime("%Y%m%d_%H%M")
# Define the function to create the model

class BiobinderModel:
    def __init__(self, config_kwargs, param_distributions=None):
        self.config_kwargs = config_kwargs
        self.param_distributions = param_distributions or {}

    def add_parameter(self, name, element_str, attr, kind, units, baseline, distribution):
        self.param_distributions[name] = {
            "element_str": element_str,  # e.g., "stream.biofuel", "unit.HTL", "PowerUtility"
            "attr": attr,
            "kind": kind,
            "units": units,
            "baseline": baseline,
            "distribution": distribution,
        }

    # def generate_samples(self, n_samples):
    #     samples = []
    #     for _ in range(n_samples):
    #         sample = {
    #             name: float(info["distribution"].sample())
    #             for name, info in self.param_distributions.items()
    #         }
    #         samples.append(sample)
    #     return samples
    
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
                param_info = self.param_distributions[name]
                element_str = param_info["element_str"]
                attr = param_info["attr"]

                if element_str == "PowerUtility":
                    setattr(qs.PowerUtility, attr, value)
                elif element_str == "TEA":
                    setattr(tea, attr, value)
                elif element_str == "LCA":
                    setattr(lca, attr, value)
                else:
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
        samples = self.generate_samples(n_samples)
        results = []
        for i, sample in enumerate(samples, 1):
            print(f"▶️ Sample {i} of {n_samples}")
            results.append(self.simulate_with_sample(sample))
        return pd.DataFrame(results)
 
if __name__ == "__main__":               
     
    configs = [
        {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
        {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
        {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
        {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
    ]
    
    
    for config in configs:
        config_name = config["name"]
        is_chcu = "CHCU" in config_name

        # Select correct Y_Bio source
        yield_df = CHCU_yield_df if is_chcu else DHCU_yield_df
    
    
        for Y_bio in yield_df['Y_Bio']:
            print(repr(Y_bio))
            print(f"{config_name} | Y_Bio = {repr(Y_bio)}")
            Y_gas = 0.1756
            Y_char = 0.0100  #fixed
            Y_aq = 1 - (Y_bio + Y_gas + Y_char)
    
            if Y_aq < 0:  # skip invalid combinations
               continue

            # Update global yield dictionary
            HTL_yields['biocrude'] = Y_bio
            HTL_yields['gas'] = Y_gas
            HTL_yields['char'] = Y_char
            HTL_yields['aqueous'] = Y_aq

            print(f"\Biocrude yield: {Y_bio:.4f}, Aqueous yield: {Y_aq:.4f}")
        
            print(f"\Running uncertainty analysis for {config['name']} configuration...")
            model = BiobinderModel(config_kwargs=config["config_kwargs"])

            # Temporary system for accessing units
            sys = create_system(**config["config_kwargs"])

            model.add_parameter("Biofuel_price", "stream.biofuel", "price", "isolated", "$", 1.13, shape.Uniform(1.02, 1.25))
            model.add_parameter("N_Fertilizer_Price", "stream.recovered_N", "price", "isolated", "$/kg N", 0.9, shape.Triangle(0.85, 0.9, 1.19))
            model.add_parameter("Natural_Gas_Price", "stream.natural_gas", "price", "isolated", "$/kg", 0.14, shape.Normal(0.14, 0.0211))
            model.add_parameter("Electricity_Price", "PowerUtility", "price", "isolated", "$/kWh", 0.074, shape.Triangle(0.03, 0.074, 0.1059))
            model.add_parameter("Uptime_Ratio", "unit.HTL", "Uptime_ratio", "coupled", "ratio", 0.9, shape.Uniform(0.8, 1.0))

            if "CHCU" in config["name"]:
                model.add_parameter("Feedstock_Transportation_Cost", "unit.FeedstockTrans", "cost", "isolated", "$/km/kg", 0.053, shape.Uniform(0.048, 0.059))
            elif "DHCU" in config["name"]:
                model.add_parameter("Biocrude_Transportation_Cost", "unit.BiocrudeTrans", "cost", "isolated", "$/km/kg", 0.099, shape.Uniform(0.089, 0.109))

            model.add_parameter("IRR", "TEA", "IRR", "coupled", "fraction", 0.10, shape.Uniform(0.05, 0.15))
            model.add_parameter("income_tax", "TEA", "income_tax", "coupled", "fraction", 0.21, shape.Uniform(0.15, 0.35))

            if config["name"] in ["CHCU_EC", "DHCU_EC"]:
                for unit_name in ["HTL_EC", "Upgrading_EC"]:
                    try:
                        _ = getattr(sys.flowsheet.unit, unit_name)  # just to confirm existence
                        model.add_parameter(f"{unit_name}_EO_voltage", f"unit.{unit_name}", "EO_voltage", "coupled", "V", 5, shape.Uniform(2.5, 10))
                        model.add_parameter(f"{unit_name}_ED_voltage", f"unit.{unit_name}", "ED_voltage", "coupled", "V", 30, shape.Uniform(2.5, 50))
                        model.add_parameter(f"{unit_name}_electrode_cost", f"unit.{unit_name}", "electrode_cost", "coupled", "$/m2", 40000, shape.Triangle(225, 40000, 80000))
                        model.add_parameter(f"{unit_name}_AEM_cost", f"unit.{unit_name}", "Anion_exchange_membrane_cost", "coupled", "$/m2", 170, shape.Triangle(0, 170, 250))
                        model.add_parameter(f"{unit_name}_CEM_cost", f"unit.{unit_name}", "Cation_exchange_membrane_cost", "coupled", "$/m2", 190, shape.Triangle(0, 190, 250))
                        model.add_parameter(f"{unit_name}_COD_removal", f"unit.{unit_name}", "COD_removal", "coupled", "fraction", 0.95, shape.Uniform(0.9, 1))
                        model.add_parameter(f"{unit_name}_N_recovery", f"unit.{unit_name}", "N_recovery", "coupled", "fraction", 0.8, shape.Uniform(0.72, 0.88))
                    except AttributeError:
                        print(f"{unit_name} not found in the flowsheet. Skipping EC parameters.")

            # Run uncertainty analysis and save to Excel
            df = model.run_uncertainty_analysis(n_samples=100)
            file_path = os.path.join(output_dir, f"UA_100_{config_name}_Ybio_{Y_bio:.5f}_{timestamp}.xlsx")
            df.to_excel(file_path, index=False)
            print(f"Saved results to {file_path}")




    


