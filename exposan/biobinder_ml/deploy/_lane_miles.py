# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aliahmad1331@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
import glob
import pandas as pd
import pyogrio

# =====================================================================
# 1. FILE PATH CONFIGURATION
# =====================================================================
HPMS_DIR = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\HPMS2024"   #https://data-usdot.opendata.arcgis.com/datasets/5e6a977c2d7c4ec1bdc82e684d3384f2/about
OUTPUT_CSV = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\county_lane_miles_summary.csv"

def main():
    print("--- Initializing High-Speed Tabular HPMS Pipeline ---")
    
    # Locate Geodatabase
    gdb_paths = glob.glob(os.path.join(HPMS_DIR, "*.gdb"))
    if not gdb_paths:
        raise FileNotFoundError(f"Could not locate a valid .gdb folder inside {HPMS_DIR}")
    gdb_path = gdb_paths[0]
    
    all_layers = [l[0] for l in pyogrio.list_layers(gdb_path)]
    print(f"Successfully connected to database. Processing {len(all_layers)} state layers...")
    
    aggregated_results = []
    
    for layer_name in all_layers:
        # Extract the state abbreviation directly from layer name (e.g., 'AK' from 'HPMS_FULL_AK_2024')
        state_parts = layer_name.split('_')
        state_abbr = state_parts[2] if len(state_parts) >= 3 else "UNKNOWN"
        
        print(f"Processing State: {state_abbr} ({layer_name})...")
        
        try:
            # Read ONLY the tabular data columns, completely bypassing the heavy spatial geometries
            df = pyogrio.read_dataframe(
                gdb_path,
                layer=layer_name,
                columns=['through_lanes', 'sectionlength', 'f_system', 'county_id'],
                read_geometry=False  # This makes processing 100x faster and eliminates RAM limits
            )
            
            if df.empty:
                continue
            
            # --- Data Cleaning ---
            df['through_lanes'] = pd.to_numeric(df['through_lanes'], errors='coerce').fillna(2).replace(0, 2)
            df['sectionlength'] = pd.to_numeric(df['sectionlength'], errors='coerce').fillna(0)
            
            # Ensure County IDs are padded strings (e.g., 2013 instead of 2013.0)
            df['county_id'] = pd.to_numeric(df['county_id'], errors='coerce').fillna(0).astype(int).astype(str)
            
            # Calculate segment lane-miles
            df['lane_miles'] = df['sectionlength'] * df['through_lanes']
            
            # Group metrics rapidly inside tabular memory
            summary = df.groupby(['county_id', 'f_system'])['lane_miles'].sum().reset_index()
            summary['state_abbr'] = state_abbr
            
            aggregated_results.append(summary)
            
        except Exception as e:
            print(f"  [Warning] Skipped layer {layer_name} due to error: {str(e)}")
            continue

    # =====================================================================
    # COMPILE AND SAVE OUTPUT MATRIX
    # =====================================================================
    if aggregated_results:
        print("\n--- Compiling Final National Aggregated Summary Table ---")
        final_df = pd.concat(aggregated_results, ignore_index=True)
        
        # Recalculate final groupings to catch cross-boundary edge cases
        final_df = final_df.groupby(['state_abbr', 'county_id', 'f_system'])['lane_miles'].sum().reset_index()
        
        # Clean column layout rules
        final_df = final_df.rename(columns={'county_id': 'county_fips', 'f_system': 'functional_system'})
        final_df = final_df.sort_values(by=['state_abbr', 'county_fips', 'functional_system'])
        
        # Export final matrix layout directly to disk
        final_df.to_csv(OUTPUT_CSV, index=False)
        print(f"\nSuccess! Custom County-level matrix exported to:\n{OUTPUT_CSV}")
        print(f"Total entries logged: {len(final_df)} matrix rows across the US.")
    else:
        print("\nPipeline failed. No data matrices could be extracted.")

if __name__ == "__main__":
    main()