#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
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

# Create a directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# Define the function to create the model
class BiobinderModel:
    def __init__(self, config_kwargs, param_distributions=None):
        self.config_kwargs = config_kwargs
        self.param_distributions = param_distributions or {}
        self.results = []
        
    def create_system(self):
        """
        Create the system based on the configuration.
        """
        sys = create_system(**self.config_kwargs)
        self.sys = sys
        self.tea = sys.TEA
        self.lca = sys.LCA
        self.biobinder = sys.flowsheet.stream.biobinder
        self.tea.IRR = tea_kwargs["IRR"]
        self.tea.income_tax = tea_kwargs["income_tax"]
        return sys

    def add_parameter(self, name, element, kind, units, baseline, distribution):
        """
        Add a parameter for uncertainty analysis.
        """
        self.param_distributions[name] = {
            "element": element,
            "kind": kind,
            "units": units,
            "baseline": baseline,
            "distribution": distribution,
        }

    def simulate_and_record(self, sample):
        """
        Simulate the system and record key metrics.
        """
        try:
           self.sys.simulate()
           MSP = self.tea.solve_price(self.biobinder)
           GWP = self.lca.get_allocated_impacts(
            streams=(self.biobinder,), operation_only=True, annual=True
        )["GWP"] / (self.biobinder.F_mass * self.lca.system.operating_hours)

        # Combine parameter values with results
           result = sample.copy()
           result["MSP ($/kg)"] = MSP
           result["GWP (kg CO2e/kg)"] = GWP

           self.results.append(result)
           print(f"Simulation successful: MSP = ${MSP:.2f}/kg, GWP = {GWP:.4f} kg CO2e/kg")
        except Exception as e:
           print(f"Simulation failed: {e}")


    def run_uncertainty_analysis(self, n_samples=1000, random_seed=30):
        """
        Run uncertainty analysis based on parameter distributions.
        """
        np.random.seed(random_seed)
        samples = self.generate_samples(n_samples)
        for sample in samples:
          for param, value in sample.items():
            element = self.param_distributions[param]["element"]

            if isinstance(value, np.ndarray):
                value = value.item()
         
            if hasattr(element, "price"):
                element.price = value
            else:
                setattr(element, param, value)
            
            self.simulate_and_record(sample)


    def generate_samples(self, n_samples):
         """
         Generate random samples based on the defined distributions.
         """
         samples = []
         for _ in range(n_samples):
           sample = {name: float(param_info["distribution"].sample()) for name, param_info in self.param_distributions.items()}
           samples.append(sample)
         return samples
     
    def single_point_sensitivity(self, etol=0.01, **kwargs):
         """Perform single-point sensitivity analysis."""
         
         if not self.results:
             print("Running baseline simulation for sensitivity analysis.")
             baseline_sample = {param: info["baseline"] for param, info in self.param_distributions.items()}
             self.simulate_and_record(baseline_sample)
             
         parameters = list(self.param_distributions.keys())
         bounds = [(param["distribution"].lower, param["distribution"].upper) for param in self.param_distributions.values()]
         metrics = ["MSP ($/kg)", "GWP (kg CO2e/kg)"]
         values_lb = np.zeros((len(parameters), len(metrics)))
         values_ub = np.zeros_like(values_lb)
         baseline_results = self.results[0]
         
         for i, param in enumerate(parameters):
             sample = baseline_results.copy()
             sample[param] = bounds[i][0]
             self.simulate_and_record(sample)
             values_lb[i, 0] = self.results[-1]["MSP ($/kg)"]
             values_lb[i, 1] = self.results[-1]["GWP (kg CO2e/kg)"]

             sample[param] = bounds[i][1]
             self.simulate_and_record(sample)
             values_ub[i, 0] = self.results[-1]["MSP ($/kg)"]
             values_ub[i, 1] = self.results[-1]["GWP (kg CO2e/kg)"]
         
         return pd.DataFrame({
             "Parameter": parameters,
             "Baseline MSP": baseline_results["MSP ($/kg)"],
             "MSP Lower Bound": values_lb[:, 0],
             "MSP Upper Bound": values_ub[:, 0],
             "Baseline GWP": baseline_results["GWP (kg CO2e/kg)"],
             "GWP Lower Bound": values_lb[:, 1],
             "GWP Upper Bound": values_ub[:, 1],
         })
     
    def distribution_sensitivity_analysis(self, n_samples=10, random_seed=42):
     """
     Perform sensitivity analysis by sampling the distribution of one parameter 
     while keeping others at baseline values.
     Saves results for each parameter.
    
     Parameters:
     n_samples (int): Number of samples to generate for each parameter.
     random_seed (int): Seed for reproducibility.
     """
     np.random.seed(random_seed)
    
     
     results = []
     parameters = list(self.param_distributions.keys())

     print(f"Running distribution sensitivity analysis with {n_samples} samples per parameter...")

     for param in parameters:
        print(f"Analyzing sensitivity for parameter: {param}")
        
        # Generate samples for the parameter being analyzed
        param_info = self.param_distributions[param]
        samples = [float(param_info["distribution"].sample()) for _ in range(n_samples)]
        
        for sample_value in samples:
            # Create a baseline sample
            sample = {p: info["baseline"] for p, info in self.param_distributions.items()}
            # Replace the value of the current parameter with the sampled value
            sample[param] = sample_value
            
            # Simulate and record
            self.simulate_and_record(sample)
            
            # Add results to the list
            results.append({
                "Parameter": param,
                "Sample Value": sample_value,
                "MSP ($/kg)": self.results[-1]["MSP ($/kg)"],
                "GWP (kg CO2e/kg)": self.results[-1]["GWP (kg CO2e/kg)"],
            })

    
     results_df = pd.DataFrame(results)
     
     return results_df



if __name__ == "__main__":
    configs = [
        {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
        {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
        {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
        {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
    ]

    for config in configs:
        print(f"Running uncertainty analysis for {config['name']} configuration...")
        model = BiobinderModel(config["config_kwargs"])
        model.create_system()

        
        # General parameters
        general_parameters = [
            ("Biofuel_price", model.sys.flowsheet.stream.biofuel, "isolated", "$", 1.132979577, shape.Uniform(1.019681619, 1.246277535)),
            ("N_Fertilizer_Price", model.sys.flowsheet.stream.recovered_N, "isolated", "$/kg N", 0.9, shape.Triangle(0.85, 0.9, 1.19)),
            ("Natural_Gas_Price", model.sys.flowsheet.stream.natural_gas, "isolated", "$/kg", 0.14, shape.Normal(0.14, 0.0211)),
            ("Electricity_Price", qs.PowerUtility, "isolated", "$/kWh", 0.074, shape.Triangle(0.03, 0.074, 0.1059)),
            ("Uptime_Ratio", model.sys.flowsheet.unit.HTL, "coupled", "ratio", 0.9, shape.Uniform(0.8, 1.0)),
        ]
        for name, element, kind, units, baseline, distribution in general_parameters:
            model.add_parameter(name=name, element=element, kind=kind, units=units, baseline=baseline, distribution=distribution)
        
        if "CHCU" in config['name']:
            model.add_parameter(
                name="Feedstock_Transportation_Cost",
                element=model.sys.flowsheet.unit.FeedstockTrans,
                kind="isolated",
                units="$/km/kg",
                baseline=0.053367347,
                distribution=shape.Uniform(0.048030612, 0.058704082),
            )
        elif "DHCU" in config['name']:
            model.add_parameter(
                name="Biocrude_Transportation_Cost",
                element=model.sys.flowsheet.unit.BiocrudeTrans,
                kind="isolated",
                units="$/km/kg",
                baseline=0.098632868,
                distribution=shape.Uniform(0.088769582, 0.108496155),
            )
        
        model.add_parameter(name="IRR", element=model.tea, kind="coupled", units="fraction", baseline=0.1, distribution=shape.Uniform(0.05, 0.15))
        model.add_parameter(name="income_tax", element=model.tea, kind="coupled", units="fraction", baseline=0.21, distribution=shape.Uniform(0.15, 0.35))
        
        if config['name'] in ["CHCU_EC", "DHCU_EC"]:
            ec_units = {"HTL_EC": None, "Upgrading_EC": None}

            for unit_name in ec_units:
                try:
                    ec_units[unit_name] = getattr(model.sys.flowsheet.unit, unit_name)
                    print(f"{unit_name} found in the flowsheet registry.")
                except AttributeError:
                    print(f"{unit_name} not found in the flowsheet registry. Skipping parameter addition.")

            ec_parameters = [
                ("EO_voltage", "V", 5, shape.Uniform(2.5, 10)),
                ("ED_voltage", "V", 30, shape.Uniform(2.5, 50)),
                ("electrode_cost", "$/m2", 40000, shape.Triangle(225, 40000, 80000)),
                ("Anion_exchange_membrane_cost", "$/m2", 170, shape.Triangle(0, 170, 250)),
                ("Cation_exchange_membrane_cost", "$/m2", 190, shape.Triangle(0, 190, 250)),
                ("COD_removal", "fraction", 0.95, shape.Uniform(0.9, 1)),
                ("N_recovery", "fraction", 0.8, shape.Uniform(0.72, 0.88)),
            ]

            for unit, ec_unit in ec_units.items():
                if ec_unit:
                    for param_name, units, baseline, distribution in ec_parameters:
                        model.add_parameter(name=param_name, element=ec_unit, kind="coupled", units=units, baseline=baseline, distribution=distribution)

        model.run_uncertainty_analysis(n_samples=1000)

        # Save results to Excel
        results_df = pd.DataFrame(model.results)
        results_file = os.path.join(output_dir, f"uncertainty_analysis_results_1000_{config['name']}.xlsx")
        results_df.to_excel(results_file, index=False)
        
        sensitivity_df = model.single_point_sensitivity()
        sensitivity_file = os.path.join(output_dir, f"sensitivity_results_{config['name']}.xlsx")
        sensitivity_df.to_excel(sensitivity_file, index=False)
        
        
        sensitivity_results = model.distribution_sensitivity_analysis(n_samples=1000)
        output_file = os.path.join(output_dir, f"sensitivity_analysis_results_{config['name']}.xlsx")
        sensitivity_results.to_excel(output_file, index=False)
                
        # print(f"Results saved to {results_file}")
        # print(f"Results saved to {sensitivity_file}")
        print(f"Sensitivity analysis results saved to {output_file}")
        
        

