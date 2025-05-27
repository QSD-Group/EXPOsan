
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
df = df[['NPDES ID2', TONS_COL]].dropna()
df = df[df[TONS_COL] > 0]
results = []
mpath = os.path.join (results_path,'biosolids_contextual_scenario_results.xlsx')

writer = pd.ExcelWriter(os.path.join(results_path, 'monte_carlo_results.xlsx'))

# Create and configure model
model = create_model(model_ID='A')
system = model.system
u = system.flowsheet.unit
# biosolids = system.flowsheet.stream.biosolids
u.A1.AC = ASH_FRAC
u.A1.FC = FC_FRAC
u.A1.MC = MOISTURE_FRAC

N = 100  # Number of Latin Hypercube samples
samples = model.sample(N, rule='L')
model.load_samples(samples)

    
for _, row in df.iterrows():
    facility = row['NPDES ID2']
    # dry_tons_per_year = row[TONS_COL]
    u.A1.biosolids_generated = tpy = row[TONS_COL]

# =============================================================================
#     # Convert to dry mass per day
#     dry_mass_kg_d = dry_tons_per_year * TONS_TO_KG / 365
#     wet_mass_kg_d = dry_mass_kg_d / (1 - MOISTURE_FRAC)
# 
#     # Update model with mass flow from facility sheet, update VM and MC based on biosolids treatment methods.
#     biosolids.F_mass = wet_mass_kg_d
#     biosolids.imass['AshContent'] = ASH_FRAC * dry_mass_kg_d # TODO: add IF statement to reduce based on treatment methods
#     biosolids.imass['VolatileMatter'] = VM_FRAC * dry_mass_kg_d
#     biosolids.imass['FixedCarbon'] = FC_FRAC * dry_mass_kg_d
#     biosolids.imass['H2O'] = MOISTURE_FRAC * wet_mass_kg_d # TODO: add IF statement to reduce based on treatment methods
# =============================================================================
    
    
    #TODO: also update scale_factor of each unit to use the TON_COL data

    # Simulate system with updated input
    system.simulate()

    # ===== Monte Carlo sampling =====    
    model.evaluate()
    
    # ===== Collect detailed metrics for each facility =====   
    metric_values = {}
    for metric in model.metrics:
        for col in model.table.columns:
            if metric.name in col[-1]:  # Match just the name part
                values = model.table[col]
                break
        else:
            raise KeyError(f"Metric '{metric.name}' not found in model.table.")
        
        mean = np.mean(values)
        lower = np.percentile(values, 5)
        upper = np.percentile(values, 95)
        metric_values[metric.name] = (mean, lower, upper)

    # Store only summarized results for each facility
    result = {'NPDES ID2': facility, 'Biosolids_tpy': tpy}
    for name, (mean, lower, upper) in metric_values.items():
        result[f'{name} (mean)'] = mean
        result[f'{name} (5th pct)'] = lower
        result[f'{name} (95th pct)'] = upper

    results.append(result)


    # Save detailed metrics to Excel
    df_metrics = model.table
    sheetname = facility[:31]  # Max Excel sheet name length = 31
    df_metrics.to_excel(writer, sheet_name=sheetname)

writer.close()

# Save summarized results to excel
export_data_df = pd.DataFrame(results)
export_data_df.to_excel(mpath)
    

