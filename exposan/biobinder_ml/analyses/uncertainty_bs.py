# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 19:57:01 2025

@author: aliah
"""

import os
import numpy as np
import pandas as pd
from datetime import datetime
from qsdsan import Model, PowerUtility
from chaospy import Uniform, Triangle, Normal
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


# Directory to save results
results_dir = 'results'
os.makedirs(results_dir, exist_ok=True)



def make_model(config_name, config_kwargs):
    sys = create_system(**config_kwargs)
    model = Model(system=sys)
    # --- Simulation and Metrics ---
    def simulate_and_return_metrics():
        sys.simulate()
        biobinder = sys.flowsheet.stream.biobinder
        MSP = sys.TEA.solve_price(biobinder)
        GWP = sys.LCA.get_allocated_impacts(
            streams=(biobinder,), operation_only=True, annual=True
        )['GWP'] / (biobinder.F_mass * sys.LCA.system.operating_hours)
        print(f" Simulated → MSP: {MSP}, GWP: {GWP}, flowrate: {biobinder.F_mass}")
        return MSP, GWP

    model.specification = lambda: sys.simulate()

    msp_metric = model.metric(
    name='MSP ($/kg)',
    units='USD/kg',
    )(lambda: simulate_and_return_metrics()[0])

    gwp_metric = model.metric(
    name='GWP (kg CO2e/kg)',
    units='kg CO2e/kg',
     )(lambda: simulate_and_return_metrics()[1])

    print(f" Registered metrics: {[m.name for m in model.metrics]}")


    # --- General Parameters ---
    model.parameter(
        name='Biofuel_price',
        element='stream.biofuel',
        kind='isolated', units='USD/kg',
        baseline=1.132979577,
        distribution=Uniform(1.019681619, 1.246277535),
    )(lambda i: setattr(sys.flowsheet.stream.biofuel, 'price', i))

    model.parameter(
        name='N_Fertilizer_Price',
        element='stream.recovered_N',
        kind='isolated', units='USD/kg N',
        baseline=0.9,
        distribution=Triangle(0.85, 0.9, 1.19),
    )(lambda i: setattr(sys.flowsheet.stream.recovered_N, 'price', i))

    model.parameter(
        name='Natural_Gas_Price',
        element='stream.natural_gas',
        kind='isolated', units='USD/kg',
        baseline=0.14,
        distribution=Normal(0.14, 0.0211),
    )(lambda i: setattr(sys.flowsheet.stream.natural_gas, 'price', i))
    
    model.parameter(
        name='Uptime Ratio',
        element='TEA',
        kind='isolated', units='ratio',
        baseline=0.9,
        distribution=Uniform(0.8, 1.0),
    )(lambda i: setattr(sys.flowsheet.unit.HTL, 'ratio', i))

    model.parameter(
        name='Electricity_Price',
        element='PowerUtility',
        kind='isolated', units='USD/kWh',
        baseline=0.074,
        distribution=Triangle(0.03, 0.074, 0.1059),
    )(lambda i: setattr(PowerUtility, 'price', i))

    # --- Scenario-Specific Parameters ---
    if 'CHCU' in config_name:
        model.parameter(
            name='Feedstock_Transportation_Cost',
            element='unit.FeedstockTrans',
            kind='isolated', units='USD/km/kg',
            baseline=0.053367347,
            distribution=Uniform(0.048030612, 0.058704082),
        )(lambda i: setattr(sys.flowsheet.unit.FeedstockTrans, 'cost', i))
    elif 'DHCU' in config_name:
        model.parameter(
            name='Biocrude_Transportation_Cost',
            element='unit.BiocrudeTrans',
            kind='isolated', units='USD/km/kg',
            baseline=0.098632868,
            distribution=Uniform(0.088769582, 0.108496155),
        )(lambda i: setattr(sys.flowsheet.unit.BiocrudeTrans, 'cost', i))

    model.parameter(
        name='IRR', element='TEA', kind='coupled', units='fraction',
        baseline=0.1, distribution=Uniform(0.05, 0.15),
    )(lambda i: setattr(sys.TEA, 'IRR', i))

    model.parameter(
        name='income_tax', element='TEA', kind='coupled', units='fraction',
        baseline=0.21, distribution=Uniform(0.15, 0.35),
    )(lambda i: setattr(sys.TEA, 'income_tax', i))

    # --- EC Unit Parameters ---
    if config['name'] in ["CHCU_EC", "DHCU_EC"]:
        for unit_name in ['HTL_EC', 'Upgrading_EC']:
            try:
                unit = getattr(sys.flowsheet.unit, unit_name)
                model.parameter(
                    name=f'{unit_name}_EO_voltage',
                    element=unit,
                    kind='coupled',
                    units='V',
                    baseline=5,
                    distribution=Uniform(2.5, 10),
                )(lambda i, u=unit: setattr(u, 'EO_voltage', i))
                model.parameter(
                    name=f'{unit_name}_ED_voltage',
                    element=unit,
                    kind='coupled',
                    units='V',
                    baseline=30,
                    distribution=Uniform(2.5, 50),
                )(lambda i, u=unit: setattr(u, 'ED_voltage', i))
                model.parameter(
                    name=f'{unit_name}_electrode_cost',
                    element=unit,
                    kind='coupled',
                    units='USD/m2',
                    baseline=40000,
                    distribution=Triangle(225, 40000, 80000),
                )(lambda i, u=unit: setattr(u, 'electrode_cost', i))
                model.parameter(
                    name=f'{unit_name}_AEM_cost',
                    element=unit,
                    kind='coupled',
                    units='USD/m2',
                    baseline=170,
                    distribution=Triangle(0, 170, 250),
                )(lambda i, u=unit: setattr(u, 'Anion_exchange_membrane_cost', i))
                model.parameter(
                    name=f'{unit_name}_CEM_cost',
                    element=unit,
                    kind='coupled',
                    units='USD/m2',
                    baseline=190,
                    distribution=Triangle(0, 190, 250),
                )(lambda i, u=unit: setattr(u, 'Cation_exchange_membrane_cost', i))
                model.parameter(
                    name=f'{unit_name}_COD_removal',
                    element=unit,
                    kind='coupled',
                    units='fraction',
                    baseline=0.95,
                    distribution=Uniform(0.9, 1),
                )(lambda i, u=unit: setattr(u, 'COD_removal', i))
                model.parameter(
                    name=f'{unit_name}_N_recovery',
                    element=unit,
                    kind='coupled',
                    units='fraction',
                    baseline=0.8,
                    distribution=Uniform(0.72, 0.88),
                )(lambda i, u=unit: setattr(u, 'N_recovery', i))
            except AttributeError:
                print(f"Unit '{unit_name}' not found in flowsheet. Skipping EC parameters.")
    print(f"Registered {len(model.parameters)} parameters for configuration: {config_name}")
    return model


# --- Run for all configs ---
configs = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "CHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": False, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central, "decentralized_HTL": True, "decentralized_upgrading": False, "skip_EC": False, "generate_H2": False, "EC_config": None}},
]

for config in configs:
    name = config['name']
    kwargs = config['config_kwargs']
    print(f"\n--- Running uncertainty analysis for {name} ---")
    model = make_model(name, kwargs)

    samples = model.sample(N=5, rule='L')
    model.load_samples(samples)

    failed_samples = []
    results = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")

    def exception_hook(e, sample):
        print(f"Simulation failed at sample {len(failed_samples) + 1} of {len(samples)}")
        print(f"   Error: {repr(e)}")
        import traceback
        traceback.print_exc()
        failed = dict(zip(model.parameters.keys(), sample))
        failed['Error'] = repr(e)
        failed_samples.append(failed)
        return [np.nan] * len(model.metrics)

    for i, sample in enumerate(model._samples, 1):
        print(f"▶️ Running sample {i} of {len(model._samples)}...")
        try:
            result = model._evaluate_sample(sample)
        except Exception as e:
            result = exception_hook(e, sample)
        MSP, GWP = result
        print(f" MSP: {MSP if not np.isnan(MSP) else 'nan'}, GWP: {GWP if not np.isnan(GWP) else 'nan'}")
        results.append(result)

    # Inject results back into the table 
    model.table.iloc[:, -len(model.metrics):] = results

    df = model.table
    df.columns = [' '.join(col).strip() if isinstance(col, tuple) else col for col in df.columns]
    result_file = os.path.join(results_dir, f"uncertainty_results_{name}_{timestamp}.xlsx")
    df.to_excel(result_file, index=False)
    print(f"Results saved: {result_file}")

    if failed_samples:
        failed_file = os.path.join(results_dir, f"failed_samples_{name}_{timestamp}.xlsx")
        pd.DataFrame(failed_samples).to_excel(failed_file, index=False)
        print(f" Failed sample info saved: {failed_file}")
    else:
        print(f" All simulations for {name} succeeded.")