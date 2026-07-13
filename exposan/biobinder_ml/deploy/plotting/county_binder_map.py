# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aliahmad1331@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd
import plotly.express as px
import numpy as np  # Added for capping the values safely

# --- 1. CONFIGURATION ---
input_csv = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\county_asphalt_binder_demand_projections.csv"
output_html_map = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\us_asphalt_demand_map_2023.html"

# Publicly hosted GeoJSON boundary file containing official US County FIPS vector outlines
GEOJSON_URL = "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"

print("Loading data matrix...")
df = pd.read_csv(input_csv)

# CHANGED HERE: Clean and pad FIPS properly to match the GeoJSON indexing keys perfectly
df['county_fips'] = df['county_fips'].astype(str).str.split('.').str[0].str.zfill(5)


# --- 2. OUTLIER HANDLING FOR VISUAL GRAPHICS ---
# Outlier counties with massive lane-mile footprints can compress the visual color scale.
color_cap = df['Binder_Demand_2023'].quantile(0.95)
print(f"95th Percentile Cap: {color_cap:,.1f} Tons (Prevents rendering wash-out)")

# CHANGED HERE: Instead of letting range_color hide values above the cap, 
# clamp them to the max color value so they stay colored on the map.
df['Display_Demand_2023'] = np.clip(df['Binder_Demand_2023'], 0, color_cap)


# --- 3. GENERATE CHOROPLETH MAP ---
print("Rendering interactive US County Choropleth Map...")
fig = px.choropleth(
    df, 
    geojson=GEOJSON_URL, 
    locations='county_fips',            
    color='Display_Demand_2023',         # CHANGED HERE: Use the clipped column for color
    color_continuous_scale="viridis",   
    scope="usa",                        
    # CHANGED HERE: Custom hover data allows you to still see the actual uncapped tonnage on hover
    hover_data={'county_fips': False, 'Display_Demand_2023': False, 'Binder_Demand_2023': ':, .1f'},
    labels={'Binder_Demand_2023': 'Actual Demand (Tons/Yr)'},
    title='<b> Annual Paving Asphalt Binder Demand by County</b>'
)

# FIXED: Changed colorcontinuousaxis_colorbar_title_text to coloraxis_colorbar_title_text
fig.update_layout(
    margin={"r":0, "t":50, "l":0, "b":0},
    title_font_size=16,
    coloraxis_colorbar_title_text="Short Tons / Year"
)

# --- 4. EXPORT MAP ---
# Save as a lightweight interactive HTML file (Drag-and-drop into a browser or PPT)
fig.write_html(output_html_map)
print(f"Map successfully generated and exported to:\n{output_html_map}")