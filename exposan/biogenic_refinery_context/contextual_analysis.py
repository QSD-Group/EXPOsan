#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 26 10:56:55 2025

@author: stetsonrowles
"""


import os
import pandas as pd
from qsdsan.utils import load_data
from exposan.biogenic_refinery_context import (
    create_system,
    data_path,
    results_path,
)

# Configuration
INPUT_EXCEL = 'Flow + Biosolids Combined Data (Major-Facilities)_subset.xlsx'
SHEET_NAME = 'Combined Sheet (Major-Only) (2)'
TONS_COL = 'Amount of Biosolids Generated'
TONS_TO_KG = 1000

# Dry solids fractions (from _units.py context)
ASH_FRAC = 0.194
VM_FRAC = 0.569
FC_FRAC = 0.237
MOISTURE_FRAC = 0.863  # 86.3% moisture

# Paths
this_dir = os.path.dirname(__file__)
input_path = os.path.join(this_dir, INPUT_EXCEL)
df = pd.read_excel(input_path, sheet_name=SHEET_NAME)
df = df[['Facility Name', TONS_COL]].dropna()
df = df[df[TONS_COL] > 0]

def run_contextual_analysis():
    # Create all systems to compare
    sysA = create_system('A')
    all_syses = [sysA]
    all_results = {}

    # Initialize result dictionary
    for metric in ('cost', 'gwp', 'fec', 'mec', 'energy', 'biochar_ton_per_day'):
        for ID in (sys.ID for sys in all_syses):
            all_results[f'{ID}_{metric}'] = {}

    # Run simulation for each facility
    for idx, row in df.iterrows():
        community = row['Facility Name']
        dry_tons_per_year = row[TONS_COL]
        dry_mass_kg_d = dry_tons_per_year * TONS_TO_KG / 365
        wet_mass_kg_d = dry_mass_kg_d / (1 - MOISTURE_FRAC)

        for sys in all_syses:
            # Set biosolids stream composition
            biosolids = sys.flowsheet.stream.biosolids
            biosolids.empty()
            biosolids.imass['AshContent'] = ASH_FRAC * dry_mass_kg_d
            biosolids.imass['VolatileMatter'] = VM_FRAC * dry_mass_kg_d
            biosolids.imass['FixedCarbon'] = FC_FRAC * dry_mass_kg_d
            biosolids.imass['H2O'] = MOISTURE_FRAC * wet_mass_kg_d

            # Simulate system
            sys.simulate()

            # Normalize per ton of biochar
            biochar_kg_d = sys.flowsheet.stream.biochar.F_mass
            biochar_ton = biochar_kg_d / 1000 

            sys_id = sys.ID
            all_results[f'{sys_id}_cost'][community] = sys.cost_per_ton_biochar
            all_results[f'{sys_id}_gwp'][community] = sys.gwp_per_ton_biochar
            all_results[f'{sys_id}_biochar'][community] = sys.biochar_generated
            all_results[f'{sys_id}_seq_C'][community] = sys.sequesterable_carbon
            all_results[f'{sys_id}_drying'][community] = sys.drying_requirement


    # Save to Excel
    results = pd.DataFrame.from_dict(all_results)
    output_path = os.path.join(results_path, 'contextual_analysis_all_systems.xlsx')
    results.to_excel(output_path)
    print(f"Contextual analysis saved to: {output_path}")

if __name__ == '__main__':
    run_contextual_analysis()
