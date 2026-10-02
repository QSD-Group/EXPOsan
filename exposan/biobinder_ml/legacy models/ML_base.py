# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 16:34:24 2025

@author: aliah
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime
from tensorflow.keras.models import load_model
import joblib
from exposan.biobinder import create_system, HTL_yields, feedstock_composition, central_dry_flowrate as default_central

# Set file paths
base_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder\models"
model_path = os.path.join(base_dir, "HTL DL Model.h5")
x_scaler_path = os.path.join(base_dir, "X_scaler.pkl")
y_scaler_path = os.path.join(base_dir, "Y_scaler.pkl")
excel_path = os.path.join(base_dir, "Experimental data.xlsx")
output_file = os.path.join(
    base_dir, f"Experimental_vs_Predicted_MSP_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx"
)

# Load model and scalers
model = load_model(model_path)
X_scaler = joblib.load(x_scaler_path)
Y_scaler = joblib.load(y_scaler_path)

# Load experimental data
data = pd.read_excel(excel_path)

# Ensure consistent feedstock type coding
# data['Feedstock Type'] = data['Feedstock Type'].replace({
#     'Woody Biomass': 0, 'Algae': 1, 'Model Compounds': 2, 'Food Waste': 3,
#     'Human Waste': 4, 'Household/Kitchen Waste': 5, 'Municipal Solid Waste': 6,
#     'Animal Manure': 7, 'Others': 8
# })

input_cols = [
    'Feedstock Type',
    'Pre-processing',
    'Carbohydrates wt%',
    'Protein wt%',
    'Lipids wt%',
    'Ash wt%',
    'Moisture',
    'C%',
    'H%',
    'O%',
    'N%',
    'S%',
    'HHV Biomass',
    'Catalyst',
    'Reactor Type',
    'Reactor Volume (mL)',
    'Solid content (w/w) %',
    'Residence Time/ Holding Time (min)',
    'Temperature (C)',
    'Solvent'
]

# Scale input and predict yield
X = data[input_cols]
X_scaled = X_scaler.transform(X)
predicted_yields_scaled = model.predict(X_scaled)
predicted_yields = Y_scaler.inverse_transform(predicted_yields_scaled).flatten()

# Initialize results list
results = []

for idx, row in data.iterrows():
    
    print(f"🔁 Running {idx+1}/{len(data)}")
    # Update feedstock composition
    moisture = row['Moisture']/100
    ash = row['Ash wt%']/100
    lipid = row['Lipids wt%']/100
    protein = row['Protein wt%']/100
    carbohydrate = row['Carbohydrates wt%']/100
    
    feedstock_composition.update({
        'Water': moisture,
        'Lipids': lipid,
        'Proteins': protein,
        'Carbohydrates': carbohydrate,
        'Ash': ash,
    })

    # Experimental biocrude yield
    Y_bio_exp = row['Biocrude wt%']/100
    Y_gas = 0.1756
    Y_char = 0.0100
    Y_aq_exp = 1 - (Y_bio_exp + Y_gas + Y_char)
    HTL_yields.update({'biocrude': Y_bio_exp, 'gas': Y_gas, 'char': Y_char, 'aqueous': Y_aq_exp})

    try:
        sys_exp = create_system(flowsheet=None, central_dry_flowrate=default_central,
                                decentralized_HTL=False, decentralized_upgrading=False,
                                skip_EC=True, generate_H2=False, EC_config=None)
        sys_exp.simulate()
        MSP_exp = sys_exp.TEA.solve_price(sys_exp.flowsheet.stream.biobinder)
    except:
        MSP_exp = np.nan

    # Predicted biocrude yield
    Y_bio_pred = predicted_yields[idx]/100
    Y_aq_pred = 1 - (Y_bio_pred + Y_gas + Y_char)
    HTL_yields.update({'biocrude': Y_bio_pred, 'gas': Y_gas, 'char': Y_char, 'aqueous': Y_aq_pred})

    try:
        sys_pred = create_system(flowsheet=None, central_dry_flowrate=default_central,
                                 decentralized_HTL=False, decentralized_upgrading=False,
                                 skip_EC=True, generate_H2=False, EC_config=None)
        sys_pred.simulate()
        MSP_pred = sys_pred.TEA.solve_price(sys_pred.flowsheet.stream.biobinder)
    except:
        MSP_pred = np.nan

    results.append({
        'Feedstock Type': row['Feedstock Type'],
        'Temperature (C)': row['Temperature (C)'],
        'Lipids wt%': lipid,
        'Protein wt%': protein,
        'Carbohydrates wt%': carbohydrate,
        'Moisture': moisture,
        'Ash wt%': ash,
        'Biocrude Yield (Exp)': Y_bio_exp,
        'MSP ($/kg) Experimental': MSP_exp,
        'Biocrude Yield (Predicted)': Y_bio_pred,
        'MSP ($/kg) Predicted': MSP_pred,
    })

# Export results
results_df = pd.DataFrame(results)
results_df.to_excel(output_file, index=False)
print(f"✅ Results saved to {output_file}")
