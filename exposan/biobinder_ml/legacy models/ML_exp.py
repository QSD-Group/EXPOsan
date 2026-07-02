# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 13:19:36 2025

@author: aliah
"""
# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')
import os
import pandas as pd
import numpy as np
from datetime import datetime
from tensorflow.keras.models import load_model
import joblib
from exposan.biobinder import create_system, HTL_yields, feedstock_composition, central_dry_flowrate as default_central

# Setup file paths
base_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder\models"
model_path = os.path.join(base_dir, "HTL DL Model.h5")
x_scaler_path = os.path.join(base_dir, "X_scaler.pkl")
y_scaler_path = os.path.join(base_dir, "Y_scaler.pkl")
excel_path = os.path.join(base_dir, "Experimental data.xlsx")
output_path = os.path.join(base_dir, f"Experimental_vs_Predicted_MSP_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx")

# Load DL model and scalers
model = load_model(model_path)
X_scaler = joblib.load(x_scaler_path)
Y_scaler = joblib.load(y_scaler_path)

# Load experimental data
data = pd.read_excel(excel_path)

# Input features (match your training)
input_cols = [
    'Feedstock Type', 'Pre-processing', 'Carbohydrates wt%', 'Protein wt%',
    'Lipids wt%', 'Ash wt%', 'Moisture', 'C%', 'H%', 'O%', 'N%', 'S%',
    'HHV Biomass', 'Catalyst', 'Reactor Type', 'Reactor Volume (mL)',
    'Solid content (w/w) %', 'Residence Time/ Holding Time (min)',
    'Temperature (C)', 'Solvent'
]

# Predict DL yield
X = data[input_cols]
X_scaled = X_scaler.transform(X)
Y_scaled = model.predict(X_scaled)
predicted_yields = Y_scaler.inverse_transform(Y_scaled).flatten()

# Define search ranges for distillation cuts
crude_light_range = np.linspace(0.01, 0.1, 5)
crude_heavy_range = np.linspace(0.8, 0.98, 5)
oil_frac_light_range = np.linspace(0.4, 0.6, 5)

# Initialize results
results = []

# Loop over experimental rows
for idx, row in data.iterrows():
    print(f"\n🔁 Running {idx+1}/{len(data)}")

    # Update feedstock composition
    feedstock_composition.update({
        'Water': row['Moisture']/100,
        'Lipids': row['Lipids wt%']/100,
        'Proteins': row['Protein wt%']/100,
        'Carbohydrates': row['Carbohydrates wt%']/100,
        'Ash': row['Ash wt%']/100,
    })

    Y_gas, Y_char = 0.1756, 0.0100

    def simulate_with_yield(Y_bio, label=""):
        Y_aq = 1 - (Y_bio + Y_gas + Y_char)
        HTL_yields.update({'biocrude': Y_bio, 'gas': Y_gas, 'char': Y_char, 'aqueous': Y_aq})
        
        for cl in crude_light_range:
          for ch in crude_heavy_range:
            for ofl in oil_frac_light_range:
                crude_fracs = [cl, ch]
                oil_fracs = [ofl, 1 - ofl]
                try:
                    sys = create_system(
                        flowsheet=None, central_dry_flowrate=default_central,
                        decentralized_HTL=False, decentralized_upgrading=False,
                        skip_EC=True, generate_H2=False, EC_config=None,
                        crude_fracs=crude_fracs, oil_fracs=oil_fracs
                    )
                    sys.simulate()
                    MSP = sys.TEA.solve_price(sys.flowsheet.stream.biobinder)

                    # Print and return on successful configuration
                    print(f"✅ Success with crude_fracs={crude_fracs}, oil_fracs={oil_fracs}")
                    return MSP, crude_fracs, oil_fracs

                except Exception as e:
                    print(f"⚠️ Failed: crude_fracs={crude_fracs}, oil_fracs={oil_fracs}, error: {e}")
                    continue
        print("❌ No valid distillation configuration found for this yield.")
        return np.nan, None, None

   

    # Run simulation for experimental yield
    Y_bio_exp = row['Biocrude wt%']
    MSP_exp, crude_exp, oil_exp = simulate_with_yield(Y_bio_exp, label="Experimental")

    # Run simulation for predicted yield
    Y_bio_pred = predicted_yields[idx]/100
    MSP_pred, crude_pred, oil_pred = simulate_with_yield(Y_bio_pred, label="Predicted")

    # Store result
    results.append({
        'Feedstock Type': row['Feedstock Type'],
        'Temperature (C)': row['Temperature (C)'],
        'Lipids wt%': row['Lipids wt%'],
        'Protein wt%': row['Protein wt%'],
        'Carbohydrates wt%': row['Carbohydrates wt%'],
        'Moisture': row['Moisture'],
        'Ash wt%': row['Ash wt%'],
        'Biocrude Yield (Exp)': Y_bio_exp,
        'MSP ($/kg) Experimental': MSP_exp,
        'crude_fracs_exp': crude_exp,
        'oil_fracs_exp': oil_exp,
        'Biocrude Yield (Predicted)': Y_bio_pred,
        'MSP ($/kg) Predicted': MSP_pred,
        'crude_fracs_pred': crude_pred,
        'oil_fracs_pred': oil_pred
    })

# Export results
results_df = pd.DataFrame(results)
results_df.to_excel(output_path, index=False)
print(f"\n✅ Results saved to {output_path}")
