#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 08:12:41 2023

@author: jiananfeng
"""

import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt
import geopy.distance, googlemaps
from exposan.htl.geospatial_HTL_systems import create_spatial_system
from exposan.htl import _m3perh_to_MGD
from qsdsan.utils import palettes
from warnings import filterwarnings

# read data
folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
#%%

sludge = pd.read_excel(folder + 'WWTP/sludge_catagorized_Seiple.xlsx', 'datasheet')

sludge['solid_fate'] = sludge['BASELINE:Process Code'].str[-1].astype(int)

# =============================================================================
# from Seiple et al. 2020
# 1: AD, dewatering
# 2: Aerobic digestion, dewatering
# 3: Dewatering, lime stabilization
# 4: AD, dewatering, thermal drying
# 5: Dewatering, multiple-hearth incineration (MHI)
# 6: Dewatering, fluidized bed incineration (FBI)
# 7: Dewatering only
# 8: Lagoon/constructed wetland
# =============================================================================

sludge = gpd.GeoDataFrame(sludge, crs='EPSG:4269',
                          geometry=gpd.points_from_xy(x=sludge.Longitude,
                                                      y=sludge.Latitude))

refinery = pd.read_excel(folder + 'refinery/petroleum_refineries_EIA.xlsx')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

for excluded in ('Alaska',
                 'Hawaii',
                 'Puerto Rico',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

sludge = sludge.to_crs(crs='EPSG:3857')
refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')


sludge = gpd.sjoin(sludge, US)
sludge = sludge.drop(['index_right'], axis=1)
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

# select WWTPs that are within 100 km of refineries, one WWTP is only assigned once to its nearest refinery

WWTP_within = sludge.sjoin_nearest(refinery, max_distance=150000, distance_col='distance')

# we have a problem here, that is CRS 3857 is not accurate at all, especially more far away from the equator, therefore I am thinking use a larger distance here, for example,
# 150000 (tested, when use 150000, all actual distance is smaller than 100000) as the max_distance and then use geopy.distance.geodesic to recalculate the distance
# and filter out those that are actually longer than 100000

WWTP_within['WWTP_location'] = list(zip(WWTP_within.Latitude_left, WWTP_within.Longitude_left))
WWTP_within['refinery_location'] = list(zip(WWTP_within.Latitude_right, WWTP_within.Longitude_right))

linear_distance = []
for i in range(len(WWTP_within)):
    linear_distance.append(geopy.distance.geodesic(WWTP_within['WWTP_location'].iloc[i], WWTP_within['refinery_location'].iloc[i]).km)

WWTP_within['linear_distance'] = linear_distance

WWTP_within = WWTP_within[WWTP_within['linear_distance'] <= 100] # select out WWTPs that are within 100 km of refineries

gmaps = googlemaps.Client(key='AIzaSyDBplFVcnCohVfZ9UuIhXPcd3N8elInQiM')
real_distance = []
for i in range(len(WWTP_within)):
    try:
        print(i)
        real_distance.append(gmaps.distance_matrix(WWTP_within['WWTP_location'].iloc[i], WWTP_within['refinery_location'].iloc[i], mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000)
    except KeyError:
        print('--------------------------------')
        real_distance.append('nan')
        
WWTP_within['real_distance'] = real_distance

# based on state of the WWTPs, determine electricity price and GHG

elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

WWTP_within = WWTP_within.merge(elec, left_on='ST', right_on='state')

# based on solid fate, determine sludge amount and biochemical composition

solid_fate = pd.read_excel(folder + 'WWTP/solid_fate.xlsx')

WWTP_within = WWTP_within.merge(solid_fate, left_on='solid_fate', right_on='category')

solid_amount = pd.read_excel(folder + 'WWTP/sludge_uncategorized_Seiple.xlsx')

WWTP_within = WWTP_within.merge(solid_amount, left_on=('FACILITY','ST','CITY'), right_on=('Facility Name','STATE','CITY'))

# WWTP_within.to_excel(folder + 'HTL_spatial_model_input.xlsx')

#%%
# this part is for WRRFs visualization

# read data

sludge = pd.read_excel(folder + 'WWTP/sludge_catagorized_Seiple.xlsx', 'datasheet')
sludge = gpd.GeoDataFrame(sludge, crs='EPSG:4269',
                          geometry=gpd.points_from_xy(x=sludge.Longitude,
                                                      y=sludge.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

for excluded in ('Alaska',
                 'Hawaii',
                 'Puerto Rico',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

sludge = sludge.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')

sludge = gpd.sjoin(sludge, US)
sludge = sludge.drop(['index_right'], axis=1)

# visualization of data
Guest = palettes['Guest']

fig, ax = plt.subplots(figsize=(30, 30))
ax.tick_params(top=False, bottom=False, left=False, right=False,
               labelleft=False, labelbottom=False)
ax.set_frame_on(False)

US.plot(ax=ax, color='w', edgecolor='k', linewidth=5)

# sludge[sludge['Influent Flow (MMGal/d)'] <= 10].plot(ax=ax, color=Guest.gray.HEX, markersize=5, alpha=0.5)
# sludge[(sludge['Influent Flow (MMGal/d)'] > 10) & (sludge['Influent Flow (MMGal/d)'] <= 100)].plot(ax=ax, color=Guest.red.HEX, markersize=100, edgecolor='k', linewidth=1.5, alpha=0.8)
# sludge[sludge['Influent Flow (MMGal/d)'] > 100].plot(ax=ax, color=Guest.green.HEX, markersize=500, edgecolor='k', linewidth=2)

# if you want to show all WRRFs together with the same symbols, comment the codes above and run the following codes
sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=10)

#%%
# this part is for oil refinery visualization

# read data

refinery = pd.read_excel(folder + 'refinery/petroleum_refineries_EIA.xlsx')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

for excluded in ('Alaska',
                 'Hawaii',
                 'Puerto Rico',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')

refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

# visualization of data
Guest = palettes['Guest']

fig, ax = plt.subplots(figsize=(30, 30))
ax.tick_params(top=False, bottom=False, left=False, right=False,
               labelleft=False, labelbottom=False)
ax.set_frame_on(False)

b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX

US.plot(ax=ax,
        color=[r, b, g, b, b, r, g, b, o,
               b, g, y, r, g, r, y, r, b,
               b, g, o, o, g, o, b, g, y,
               g, b, o, g, b, b, y, b, b,
               b, b, b, b, g, g, g, y, g,
               r, g, g, b],
        alpha=0.4,
        edgecolor='k',
        linewidth=0)

US.plot(ax=ax,
        color='none',
        edgecolor='k',
        linewidth=5)

refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=1000, edgecolor='k', linewidth=5)

#%%
# this part is for WRRFs+oil refineries visualization (just select states)

# read data

sludge = pd.read_excel(folder + 'WWTP/sludge_catagorized_Seiple.xlsx', 'datasheet')
sludge = gpd.GeoDataFrame(sludge, crs='EPSG:4269',
                          geometry=gpd.points_from_xy(x=sludge.Longitude,
                                                      y=sludge.Latitude))

refinery = pd.read_excel(folder + 'refinery/petroleum_refineries_EIA.xlsx')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

US = US.loc[US['NAME'].isin(('Pennsylvania',))]

sludge = sludge.to_crs(crs='EPSG:3857')
refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')

sludge = gpd.sjoin(sludge, US)
sludge = sludge.drop(['index_right'], axis=1)
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

# visualization of data
Guest = palettes['Guest']

fig, ax = plt.subplots(figsize=(30, 30))
ax.tick_params(top=False, bottom=False, left=False, right=False,
               labelleft=False, labelbottom=False)
ax.set_frame_on(False)

US.plot(ax=ax, color=Guest.blue.HEX, alpha=0.4, edgecolor='k', linewidth=0)
US.plot(ax=ax, color='none', edgecolor='k', linewidth=10)
sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=200)
refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=5000, edgecolor='k', linewidth=10)

#%%
# this part is for WRRFs+oil refineries visualization

# read data

sludge = pd.read_excel(folder + 'WWTP/sludge_catagorized_Seiple.xlsx', 'datasheet')
sludge = gpd.GeoDataFrame(sludge, crs='EPSG:4269',
                          geometry=gpd.points_from_xy(x=sludge.Longitude,
                                                      y=sludge.Latitude))

refinery = pd.read_excel(folder + 'refinery/petroleum_refineries_EIA.xlsx')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

for excluded in ('Alaska',
                 'Hawaii',
                 'Puerto Rico',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

sludge = sludge.to_crs(crs='EPSG:3857')
refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')

sludge = gpd.sjoin(sludge, US)
sludge = sludge.drop(['index_right'], axis=1)
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

# visualization of data
Guest = palettes['Guest']

fig, ax = plt.subplots(figsize=(30, 30))
ax.tick_params(top=False, bottom=False, left=False, right=False,
               labelleft=False, labelbottom=False)
ax.set_frame_on(False)

US.plot(ax=ax, color='w', edgecolor='k')
sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=20)
refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=300, edgecolor='k', linewidth=1.5)

#%%
# this part is for electricity price and carbon intensity visualization

# read data

electricity = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')


US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

for excluded in ('Alaska',
                 'Hawaii',
                 'Puerto Rico',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

electricity = pd.merge(electricity, US, left_on='name', right_on='NAME')

electricity = gpd.GeoDataFrame(electricity, crs='EPSG:4269')

# visualization of data
Guest = palettes['Guest']


fig, ax = plt.subplots(figsize=(30, 30))
ax.tick_params(top=False, bottom=False, left=False, right=False,
               labelleft=False, labelbottom=False)
ax.set_frame_on(False)

# electricity.plot('price (10-year median)', ax=ax, cmap='Oranges', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35})

electricity.plot('GHG (10-year median)', ax=ax, cmap='Blues', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35})


#%%
# this part is the analysis for WRRFs and oil refineries co-location

WWTP_within = pd.read_excel(folder + 'HTL_spatial_model_input.xlsx')

result = WWTP_within[['site_id']].drop_duplicates()

for distance in np.linspace(0, 100, 101):
    WWTP_within_distance = WWTP_within[WWTP_within['linear_distance']<=distance]
    WWTP_within_distance = WWTP_within_distance.groupby('site_id').sum('Influent Flow (MMGal/d)')
    WWTP_within_distance = WWTP_within_distance[['Influent Flow (MMGal/d)']]
    WWTP_within_distance = WWTP_within_distance.rename(columns={'Influent Flow (MMGal/d)': int(distance)})
    WWTP_within_distance.reset_index(inplace=True)
    if len(WWTP_within_distance) > 0:
        result = result.merge(WWTP_within_distance, how='left', on='site_id')

result = result.fillna(0)

refinery = pd.read_excel(folder + 'refinery/petroleum_refineries_EIA.xlsx')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME', 'geometry']]

for excluded in ('Alaska',
                 'Hawaii',
                 'Puerto Rico',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')

refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

if len(result) < len(refinery):
    result_index = pd.Index(result.site_id)
    refinery_index = pd.Index(refinery.site_id)
    refineries_left_id = refinery_index.difference(result_index).values

result = result.set_index('site_id')
for item in refineries_left_id:
    result.loc[item] = [0]*len(result.columns)

result = result.merge(refinery, how='left', on='site_id')

result.to_excel(folder + 'results/100km_results_MGD_vs_distance.xlsx')

#%%
# this part is the analysis for CO2 abatement cost
# run analysis for every row

filterwarnings('ignore')

WWTP_within = pd.read_excel(folder + 'HTL_spatial_model_input.xlsx')

CO2_reduction = []
USD_decarbonization = []
oil_BPD = []

# for i in range(len(WWTP_within)):
    
for i in range(308, 309): # open multiple consoles to speed up analysis
    waste_GHG_value = 800
    solid_reduction_value = WWTP_within.iloc[i]['reduction_factor_ave']
    sys, barrel = create_spatial_system(waste_price=400,
                                        waste_GHG=waste_GHG_value,
                                        size=WWTP_within.iloc[i]['Influent Flow (MMGal/d)'],
                                        distance=WWTP_within.iloc[i]['real_distance'],
                                        solid_fate=WWTP_within.iloc[i]['solid_fate'],
                                        ww_2_dry_sludge_ratio=WWTP_within.iloc[i]['TOTAL SOLIDS TONS_DAY']/WWTP_within.iloc[i]['Influent Flow (MMGal/d)'],
                                        solid_reduction=WWTP_within.iloc[i]['reduction_factor_ave'],
                                        state=WWTP_within.iloc[i]['state'],
                                        elec_GHG=WWTP_within.iloc[i]['GHG (10-year median)'])

    lca = sys.LCA
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    WWTP = unit.WWTP
    raw_wastewater = stream.raw_wastewater

    kg_CO2_per_ton_dry_sludge = lca.get_total_impacts(exclude=(raw_wastewater,))['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)/lca.lifetime

    CO2_reduction_result = WWTP_within.iloc[i]['TOTAL SOLIDS TONS_DAY']*(1-solid_reduction_value)*7920/24*30*(waste_GHG_value-kg_CO2_per_ton_dry_sludge)

    USD_per_ton_CO2_reduction = -sys.TEA.NPV/CO2_reduction_result*1000 # save money: the results are negative; spend money: the results are positive
    
    CO2_reduction.append(CO2_reduction_result)
    USD_decarbonization.append(USD_per_ton_CO2_reduction)
    oil_BPD.append(barrel)
    
    if USD_per_ton_CO2_reduction < 0:
        print('HERE WE GO!')
    
    if i%5 == 0:
        print(i) # check the progress
        
result = {'CO2_reduction': CO2_reduction,
          'USD_decarbonization': USD_decarbonization,
          'oil_BPD': oil_BPD}
        
result = pd.DataFrame(result)