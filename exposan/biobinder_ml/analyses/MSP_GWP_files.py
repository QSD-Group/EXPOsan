# -*- coding: utf-8 -*-
"""
Created on Sun Jun  1 20:47:47 2025

@author: aliah
"""

import os
import pandas as pd
import numpy as np

# Folder to search
input_dir = "results"  # Update this if your Excel files are in a different folder
output_file = "valid_msp_gwp_files.xlsx"

# List to store valid file names
valid_files = []

# Scan all .xlsx files in the folder
for file in os.listdir(input_dir):
    if file.endswith(".xlsx"):
        file_path = os.path.join(input_dir, file)
        try:
            df = pd.read_excel(file_path)

            # Check required columns exist
            if "MSP ($/kg)" in df.columns and "GWP (kg CO2e/kg)" in df.columns:
                # Drop rows where either column is NaN
                filtered = df[["MSP ($/kg)", "GWP (kg CO2e/kg)"]].dropna()
                if not filtered.empty:
                    valid_files.append(file)
        except Exception as e:
            print(f"❌ Failed to process {file}: {e}")

# Save to Excel
if valid_files:
    df_out = pd.DataFrame(valid_files, columns=["Valid Files"])
    df_out.to_excel(os.path.join(input_dir, output_file), index=False)
    print(f"✅ Saved list of valid files to '{output_file}' in folder '{input_dir}'")
else:
    print("⚠️ No valid files found.")
