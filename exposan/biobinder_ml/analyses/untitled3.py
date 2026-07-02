# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 14:34:04 2025

@author: aliah
"""

# ==============================================================
#  Food Waste Disposal Savings Map Export (300 dpi PNG + HTML)
#  Based on composition-adjusted MSW Type 10 dataset
#  Author: Ali Ahmad (Rutgers University)
#  ==============================================================
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import os as os

# --------------------------------------------------------------
# Load your dataset (update path if needed)
# --------------------------------------------------------------
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

df = pd.read_excel(os.path.join(output_dir, "MSW_Type10_By_State_CompositionBased_FoodFee1268.xlsx"))
# timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# --------------------------------------------------------------
# Compute savings metrics if not already in the file
# --------------------------------------------------------------
if "Food Disposal Cost (Million $, baseline)" not in df.columns:
    df["Food Disposal Cost (Million $, baseline)"] = (
        0.5 * df["Food (wet tons)"] * df["Fee Used ($/ton)"] / 1e6
    )

if "Food Disposal Cost (Million $, revised)" not in df.columns:
    df["Food Disposal Cost (Million $, revised)"] = (
        0.5 * df["Food (wet tons)"] * df["Revised Food Fee ($/wet ton)"] / 1e6
    )

df["Food Savings (Million $)"] = (
    df["Food Disposal Cost (Million $, baseline)"]
    - df["Food Disposal Cost (Million $, revised)"]
)
df["Food Savings (%)"] = (
    100
    * df["Food Savings (Million $)"]
    / df["Food Disposal Cost (Million $, baseline)"].replace(0, float("nan"))
)

# --------------------------------------------------------------
# Compute national summary values
# --------------------------------------------------------------
national_savings_billion = df["Food Savings (Million $)"].sum() / 1000
revised_fee = df["Revised Food Fee ($/wet ton)"].iloc[0]

# --------------------------------------------------------------
# Create choropleth map
# --------------------------------------------------------------
fig = px.choropleth(
    df,
    locations="State",
    locationmode="USA-states",
    color="Food Savings (Million $)",
    hover_name="State",
    hover_data={
        "Food Savings (Million $)": True,
        "Food Savings (%)": True,
        "Food Disposal Cost (Million $, baseline)": True,
        "Food Disposal Cost (Million $, revised)": True,
    },
    scope="usa",
    color_continuous_scale="Greens",
    title=(
        f"Food Waste Disposal Savings by State (50% Landfilled, Revised Fee = ${revised_fee:.2f}/wet ton)<br>"
        f"National Savings ≈ ${national_savings_billion:.2f} B per year"
    ),
    labels={"Food Savings (Million $)": "Savings (Million $)"},
)

# --------------------------------------------------------------
# Add state abbreviations overlay (for small states)
# --------------------------------------------------------------
state_coords = {
    'AL':[32.8,-86.8],'AK':[61.4,-152.4],'AZ':[33.7,-111.4],'AR':[34.97,-92.37],
    'CA':[36.1,-119.7],'CO':[39.06,-105.3],'CT':[41.6,-72.75],'DE':[39.3,-75.5],
    'FL':[27.76,-81.68],'GA':[33.04,-83.64],'HI':[21.09,-157.49],'ID':[44.24,-114.47],
    'IL':[40.34,-88.98],'IN':[39.84,-86.25],'IA':[42.01,-93.21],'KS':[38.52,-96.72],
    'KY':[37.66,-84.67],'LA':[31.16,-91.86],'ME':[44.69,-69.38],'MD':[39.06,-76.8],
    'MA':[42.23,-71.53],'MI':[43.32,-84.53],'MN':[45.69,-93.9],'MS':[32.74,-89.67],
    'MO':[38.45,-92.28],'MT':[46.92,-110.45],'NE':[41.12,-98.26],'NV':[38.31,-117.05],
    'NH':[43.45,-71.56],'NJ':[40.29,-74.52],'NM':[34.84,-106.24],'NY':[42.16,-74.94],
    'NC':[35.63,-79.8],'ND':[47.52,-99.78],'OH':[40.38,-82.76],'OK':[35.56,-96.92],
    'OR':[44.57,-122.07],'PA':[40.59,-77.2],'RI':[41.68,-71.51],'SC':[33.85,-80.94],
    'SD':[44.29,-99.43],'TN':[35.74,-86.69],'TX':[31.05,-97.56],'UT':[40.15,-111.86],
    'VT':[44.04,-72.71],'VA':[37.76,-78.16],'WA':[47.4,-121.49],'WV':[38.49,-80.95],
    'WI':[44.26,-89.61],'WY':[42.75,-107.3],'DC':[38.89,-77.02]
}
latitudes  = [v[0] for v in state_coords.values()]
longitudes = [v[1] for v in state_coords.values()]
state_abbr = list(state_coords.keys())

fig.add_trace(go.Scattergeo(
    locationmode='USA-states',
    lon=longitudes, lat=latitudes,
    text=state_abbr, mode='text',
    textfont=dict(color='black', size=9, family="Arial Black"),
    showlegend=False
))

# --------------------------------------------------------------
# Export high-DPI PNG and interactive HTML
# --------------------------------------------------------------
pio.write_image(
    fig,
    "FoodWaste_Savings_Map_withLabels_300dpi.png",
    format="png",
    scale=4,   # 4× base DPI ≈ 300–400 dpi
    width=2400, height=1400
)

fig.write_html("FoodWaste_Savings_Map_withLabels_Interactive.html")

print("✅ Export complete!")
print("Saved files:")
print(" - FoodWaste_Savings_Map_withLabels_300dpi.png (high-DPI PNG)")
print(" - FoodWaste_Savings_Map_withLabels_Interactive.html (interactive map)")
