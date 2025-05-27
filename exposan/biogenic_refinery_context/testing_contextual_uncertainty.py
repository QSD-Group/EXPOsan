
import os
import pandas as pd
import numpy as np
from exposan.biogenic_refinery_context.models import create_model
from exposan.biogenic_refinery_context import results_path

# Define constants here since they are not exposed in __init__.py
TONS_TO_KG = 907.18474  # 1 ton = 907.18474 kg
ASH_FRAC = 0.194
VM_FRAC = 0.569
FC_FRAC = 0.237
MOISTURE_FRAC = 0.863

# File configuration
INPUT_EXCEL = 'Flow + Biosolids Combined Data (Major-Facilities)_subset.xlsx'
SHEET_NAME = 'Combined Sheet (Major-Only) (2)'
TONS_COL = 'Amount of Biosolids Generated'

# Load input Excel
this_dir = os.path.dirname(__file__)
input_path = os.path.join(this_dir, INPUT_EXCEL)
df = pd.read_excel(input_path, sheet_name=SHEET_NAME)
df = df[['Facility Name', TONS_COL]].dropna()
df = df[df[TONS_COL] > 0]

def run_monte_carlo(N=100):
    writer = pd.ExcelWriter(os.path.join(results_path, 'monte_carlo_results.xlsx'))

    for _, row in df.iterrows():
        facility = row['Facility Name']
        dry_tons_per_year = row[TONS_COL]

        # Convert to dry mass per day
        dry_mass_kg_d = dry_tons_per_year * TONS_TO_KG / 365
        wet_mass_kg_d = dry_mass_kg_d / (1 - MOISTURE_FRAC)

        # Create and configure model
        model = create_model(model_ID='A')
        system = model.system
        biosolids = system.flowsheet.stream.biosolids
        biosolids.empty()
        biosolids.imass['AshContent'] = ASH_FRAC * dry_mass_kg_d
        biosolids.imass['VolatileMatter'] = VM_FRAC * dry_mass_kg_d
        biosolids.imass['FixedCarbon'] = FC_FRAC * dry_mass_kg_d
        biosolids.imass['H2O'] = MOISTURE_FRAC * wet_mass_kg_d

        # Run sampling and evaluation
        breakpoint()
        model.load_default_parameters()
        model.sample(N, rule='L')
        model.evaluate()

        # Save metrics to Excel
        df_metrics = model.table
        sheetname = facility[:31]  # Max Excel sheet name length = 31
        df_metrics.to_excel(writer, sheet_name=sheetname)

    writer.close()
    print("Monte Carlo results saved to monte_carlo_results.xlsx")

if __name__ == '__main__':
    run_monte_carlo(N=100)
