
import os
import pandas as pd
import numpy as np
from exposan.biogenic_refinery_context.models import create_model
from exposan.biogenic_refinery_context import results_path

# Define constants here since they are not exposed in __init__.py
TONS_TO_KG = 1000  # 1 ton = 1000 kg
ASH_FRAC = 0.34         #TODO: change to uniform distr.     np.random.uniform(0.34, 0.42)  
FC_FRAC = 0.0625         #TODO: change to uniform distr.     np.random.uniform(0.025, 0.1)
MOISTURE_FRAC = 0.863   # np.random.uniform(.9, 0.97)
VM_FRAC = 1 - ASH_FRAC - FC_FRAC

# File configuration
INPUT_EXCEL = 'Flow + Biosolids Combined Data (Major-Facilities).xlsx'
SHEET_NAME = 'Combined Sheet (Major-Only)'
TONS_COL = 'Amount of Biosolids Generated'

# Define all the columns we need
columns_to_keep = [
    'NPDES ID2',
    TONS_COL,
    'All PSRP',
    'All Other',
    'All PFRP',
    'All PTO'
]

# Load input Excel
this_dir = os.path.dirname(__file__)
input_path = os.path.join(this_dir, INPUT_EXCEL)
df = pd.read_excel(input_path, sheet_name=SHEET_NAME)

# Keep only the necessary columns
df = df[columns_to_keep].dropna(subset=[TONS_COL])
df = df[df[TONS_COL] > 0]

results = []
mpath = os.path.join (results_path,'biosolids_contextual_scenario_results.xlsx')

writer = pd.ExcelWriter(os.path.join(results_path, 'monte_carlo_results.xlsx'))

# Create and configure model
model = create_model(model_ID='A')
system = model.system
u = system.flowsheet.unit

u.A1.AC = ASH_FRAC
u.A1.FC = FC_FRAC
u.A1.MC = MOISTURE_FRAC
count = 0

N = 100  # Number of Latin Hypercube samples
samples = model.sample(N, rule='L')
model.load_samples(samples)

    
for _, row in df.iterrows():
    facility = row['NPDES ID2']
    # dry_tons_per_year = row[TONS_COL]
    
    u.A1.biosolids_generated = tpy = row[TONS_COL]
    u.A2.biosolids_generated = u.A3.biosolids_generated = u.A4.biosolids_generated = \
        u.A5.biosolids_generated = u.A6.biosolids_generated = u.A7.biosolids_generated = u.A8.biosolids_generated = tpy
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
    management_cols = ['All PSRP', 'All Other', 'All PFRP', 'All PTO']
    management_processes = []
    for col in management_cols:
        val = row.get(col, '')
        if pd.notna(val) and val.strip() != '':
        # If the value is a string, split it
            val_str = str(val).lower()
            processes = [proc.strip() for proc in val_str.split(',') if proc.strip()]
            management_processes.extend(processes)
            # Remove duplicates
    management_processes = list({proc for proc in management_processes})
    
    MC = MOISTURE_FRAC
    
    if 'thickening' in management_processes:
        MC = np.random.uniform(0.94, 0.99)
    if 'airdry' in management_processes:
        MC = 0.6 * MC
    if 'heatdry' in management_processes:
        MC = np.random.uniform(0.09, 0.11) 
      
    VM = VM_FRAC
    AC = ASH_FRAC
    FC = FC_FRAC
    
    if 'anaerobic' in management_processes and 'aerobic' not in management_processes:
        VSR = np.random.uniform(0.38, 0.5)
        #Adjust ratios based on volatile matter having been removed from the raw biosolids
        VMr = tpy*VM/(1 - VM + VSR * VM) #tons/yr, initial mass of volatile matter saw before digestion 
        ACi = AC
        FCi = FC
        AC = (ACi * (tpy+VMr-VSR*VMr)) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
        FC = (FCi * (tpy+VMr-VSR*VMr)) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
        VM = (VSR*VMr) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
    elif 'aerobic' in management_processes and 'anaerobic' not in management_processes:
        VSR = np.random.uniform(0.45, 0.5)
        #Adjust ratios based on volatile matter having been removed from the raw biosolids
        VMr = tpy*VM/(1 - VM + VSR * VM) #tons/yr, initial mass of volatile matter saw before digestion 
        ACi = AC
        FCi = FC
        AC = (ACi * (tpy+VMr-VSR*VMr)) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
        FC = (FCi * (tpy+VMr-VSR*VMr)) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
        VM = (VSR*VMr) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
    elif 'anaerobic' in management_processes and 'aerobic' in management_processes:
        VSR = np.random.uniform(0.50, 0.60)
        #Adjust ratios based on volatile matter having been removed from the raw biosolids
        VMr = tpy*VM/(1 - VM + VSR * VM) #tons/yr, initial mass of volatile matter saw before digestion 
        ACi = AC
        FCi = FC
        AC = (ACi * (tpy+VMr-VSR*VMr)) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
        FC = (FCi * (tpy+VMr-VSR*VMr)) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
        VM = (VSR*VMr) / (ACi * (tpy+VMr-VSR*VMr) + FCi * (tpy+VMr-VSR*VMr) + VSR*VMr)
    
    u.A1.MC = MC
    u.A1.VM_db = VM
    u.A1.AC = AC
    u.A1.FC = FC
    
    #TODO: also update scale_factor of each unit to use the TON_COL data

    # Simulate system with updated input
    system.simulate()

    # ===== Monte Carlo sampling =====    
    model.evaluate()
    
    count += 1
    print(count)
    
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
        median = np.median(values)
        lower = np.percentile(values, 5)
        upper = np.percentile(values, 95)
        metric_values[metric.name] = (mean, median, lower, upper)

    # Store only summarized results for each facility
    result = {'NPDES ID2': facility, 'Biosolids_tpy': tpy}
    for name, (mean, median, lower, upper) in metric_values.items():
        result[f'{name} (mean)'] = mean
        result[f'{name} (median)'] = median
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
    

