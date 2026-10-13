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

input_csv = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\county_asphalt_binder_demand_projections_updated.csv"
output_html_map = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\us_asphalt_demand_map_2024.html"

GEOJSON_URL = "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"

df = pd.read_csv(input_csv, dtype={"county_fips": str})
df["county_fips"] = (
    df["county_fips"].str.replace(r"\.0$", "", regex=True).str.zfill(5)
)

demand_col = "Binder_Demand_2024"
if demand_col not in df.columns:
    raise KeyError(
        f"{demand_col} is missing. Rerun the updated _asphalt_req.py first."
    )

# Cap the displayed color at the 95th percentile.
# The actual, uncapped demand remains available on hover.
color_cap = df[demand_col].quantile(0.95)
df["Display_Demand_2024"] = df[demand_col].clip(lower=0, upper=color_cap)

print(f"95th percentile color cap: {color_cap:,.1f} short tons/year")

fig = px.choropleth(
    df,
    geojson=GEOJSON_URL,
    locations="county_fips",
    color="Display_Demand_2024",
    color_continuous_scale="viridis",
    scope="usa",
    hover_data={
        "county_fips": False,
        "Display_Demand_2024": False,
        demand_col: ":,.1f",
    },
    labels={
        "Display_Demand_2024": "Demand (short tons/year)",
        demand_col: "Actual demand (short tons/year)",
    },
    title="Annualized Asphalt Binder Demand by County (2024 Road Inventory)",
)

fig.update_layout(
    margin={"r": 0, "t": 50, "l": 0, "b": 0},
    title_font_size=16,
    coloraxis_colorbar_title_text="Short tons/year",
)

fig.write_html(output_html_map)
print(f"Map exported to: {output_html_map}")