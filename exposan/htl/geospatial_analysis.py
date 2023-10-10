#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 08:12:41 2023

@author: jiananfeng
"""

from exposan.htl.geospatial_HTL_systems import create_spatial_system
from exposan.htl import _m3perh_to_MGD
from qsdsan.utils import palettes
import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import geopy.distance, googlemaps
from datetime import date
from warnings import filterwarnings

# read data
folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'

sludge = pd.read_excel(folder + 'WWTP/sludge_catagorized_Seiple.xlsx', 'datasheet')

sludge['solid_fate'] = sludge['BASELINE:Process Code'].str[-1].astype(int)

# from Seiple et al. 2020
# 1: AD, dewatering
# 2: Aerobic digestion, dewatering
# 3: Dewatering, lime stabilization
# 4: AD, dewatering, thermal drying
# 5: Dewatering, multiple-hearth incineration (MHI)
# 6: Dewatering, fluidized bed incineration (FBI)
# 7: Dewatering only
# 8: Lagoon/constructed wetland

sludge = gpd.GeoDataFrame(sludge, crs='EPSG:4269',
                          geometry=gpd.points_from_xy(x=sludge.Longitude,
                                                      y=sludge.Latitude))

refinery = pd.read_excel(folder + 'refinery/petroleum_refineries_EIA.xlsx')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
# confirmed from the wesbite (next line): this is no change in US map since June 1, 1995. So, using 2018 US map data is OK.
# https://en.wikipedia.org/wiki/Territorial_evolution_of_the_United_States#1946%E2%80%93present_(Decolonization)
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

# if only select states, uncomment the next line of code
# US = US.loc[US['NAME'].isin(('Pennsylvania',))]

sludge = gpd.sjoin(sludge, US)
sludge = sludge.drop(['index_right'], axis=1)
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')
solid_fate = pd.read_excel(folder + 'WWTP/solid_fate.xlsx')
solid_amount = pd.read_excel(folder + 'WWTP/sludge_uncategorized_Seiple.xlsx')

electricity = pd.merge(elec, US, left_on='name', right_on='NAME')
electricity = gpd.GeoDataFrame(electricity)

# color palette
Guest = palettes['Guest']
b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX

def set_plot(figure_size=(30,30)):
    global fig, ax
    fig, ax = plt.subplots(figsize=(30, 30))
    ax.tick_params(top=False, bottom=False, left=False, right=False,
                   labelleft=False, labelbottom=False)
    ax.set_frame_on(False)

#%%
# WRRFs visualization

set_plot()

US.plot(ax=ax, color='w', edgecolor='k', linewidth=5) # for figures in SI, use default linewidth

less_than_10 = sludge['Influent Flow (MMGal/d)'] <= 10
sludge[less_than_10].plot(ax=ax, color=Guest.gray.HEX, markersize=5, alpha=0.5)

between_10_and_100 = (sludge['Influent Flow (MMGal/d)'] > 10) & (sludge['Influent Flow (MMGal/d)'] <= 100)
sludge[between_10_and_100].plot(ax=ax, color=Guest.red.HEX, markersize=100, edgecolor='k', linewidth=1.5, alpha=0.8)

more_than_100 = sludge['Influent Flow (MMGal/d)'] > 100
sludge[more_than_100].plot(ax=ax, color=Guest.green.HEX, markersize=500, edgecolor='k', linewidth=2)

# if you want to show all WRRFs together with the same symbols, comment the codes above and uncomment the next line of code
# sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=10)

#%%
# oil refinery visualization

set_plot()

US.plot(ax=ax,
        color=[r, b, g, b, b, r, g, b, o, b, g, y, r, g, r, y, r,
               b, b, g, o, o, g, o, b, g, y, g, b, o, g, b, b, y,
               b, b, b, b, b, b, g, g, g, y, g, r, g, g, b],
        alpha=0.4,
        edgecolor='k',
        linewidth=0)

US.plot(ax=ax, color='none', edgecolor='k', linewidth=5) # for figures in SI, use default linewidth

refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=1000, edgecolor='k', linewidth=5)

#%%
# WRRFs+oil refineries visualization (just select states)

set_plot()

US.plot(ax=ax, color=Guest.blue.HEX, alpha=0.4, edgecolor='k', linewidth=0)
US.plot(ax=ax, color='none', edgecolor='k', linewidth=10)

sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=200)

refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=5000, edgecolor='k', linewidth=10)

#%%
# WRRFs+oil refineries visualization

set_plot()

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3.5)
sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=15)
refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=500, edgecolor='k', linewidth=3.5)

#%%
# electricity price and carbon intensity visualization

set_plot()

# electricity.plot('price (10-year median)', ax=ax, cmap='Oranges', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35})

electricity.plot('GHG (10-year median)', ax=ax, cmap='Blues', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35})

#%%
# transporation distance calculation

# if want select WWTPs that are within certain ranges of refineries, set max_distance in sjoin_nearest.
# here we do not set a max distance

WWTP_within = sludge.sjoin_nearest(refinery, max_distance=None, distance_col='distance')

# if set max_distance = 100000, we have a problem here, that is CRS 3857 is not accurate at all,
# especially more far away from the equator, therefore can use a larger distance here,
# for example, 150000 (tested, when use 150000, all actual distances are smaller than 100000)
# as the max_distance and then use geopy.distance.geodesic to recalculate the distance and
# filter out those that are actually longer than 100000

WWTP_within['WWTP_location'] = list(zip(WWTP_within.Latitude_left, WWTP_within.Longitude_left))
WWTP_within['refinery_location'] = list(zip(WWTP_within.Latitude_right, WWTP_within.Longitude_right))

linear_distance = []
for i in range(len(WWTP_within)):
    linear_distance.append(geopy.distance.geodesic(WWTP_within['WWTP_location'].iloc[i], WWTP_within['refinery_location'].iloc[i]).km)

WWTP_within['linear_distance_km'] = linear_distance

# if we set the max_distance, then uncomment the next line (replace 100 km when necessary)
# WWTP_within = WWTP_within[WWTP_within['linear_distance'] <= 100]

gmaps = googlemaps.Client(key='XXX') # !!! get a google API key, and do not upload to GitHub
real_distance = []
for i in range(len(WWTP_within)):
    try:
        print(i)
        real_distance.append(gmaps.distance_matrix(WWTP_within['WWTP_location'].iloc[i],
                                                   WWTP_within['refinery_location'].iloc[i],
                                                   mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000)
    except KeyError:
        print('--------------------------------')
        real_distance.append('nan')
        
WWTP_within['real_distance_km'] = real_distance

# based on state of the WWTPs, determine electricity price and GHG

WWTP_within = WWTP_within.merge(elec, left_on='ST', right_on='state')

# based on solid fate, determine sludge amount and biochemical composition

WWTP_within = WWTP_within.merge(solid_fate, left_on='solid_fate', right_on='category')

WWTP_within = WWTP_within.merge(solid_amount, left_on=('FACILITY','ST','CITY'), right_on=('Facility Name','STATE','CITY'))

WWTP_within.to_excel(folder + f'HTL_geospatial_model_input_{date.today()}.xlsx') # this will be the input for the future analysis

#%%
# WRRFs GHG map (after filter out some WRRFs)

WWTP_within = pd.read_excel(folder + 'HTL_geospatial_model_input_final.xlsx')

WWTP_within = gpd.GeoDataFrame(WWTP_within, crs='EPSG:4269',
                               geometry=gpd.points_from_xy(x=WWTP_within.Longitude_left,
                                                           y=WWTP_within.Latitude_left))

WWTP_within = WWTP_within.to_crs(crs='EPSG:3857')

set_plot()

US.plot(ax=ax, color='w', edgecolor='k')

# use k-ton/yr

less_than_10 = WWTP_within['30_years_emission_ton_CO2']/30/1000 <= 10
WWTP_within[less_than_10].plot(ax=ax, color=Guest.gray.HEX, markersize=5, alpha=0.5)

between_10_and_100 = (WWTP_within['30_years_emission_ton_CO2']/30/1000 > 10) & (WWTP_within['30_years_emission_ton_CO2']/30/1000 <= 100)
WWTP_within[between_10_and_100].plot(ax=ax, color=Guest.red.HEX, markersize=100, edgecolor='k', linewidth=1.5, alpha=0.8)

more_than_100 = WWTP_within['30_years_emission_ton_CO2']/30/1000 > 100
WWTP_within[more_than_100].plot(ax=ax, color=Guest.green.HEX, markersize=500, edgecolor='k', linewidth=2)

#%%
# cumulative WRRFs capacity vs distances

# remember to use the correct file
WWTP_within = pd.read_excel(folder + 'HTL_geospatial_model_input_final.xlsx')

result = WWTP_within[['site_id']].drop_duplicates()

for distance in np.linspace(0, 100, 1001):
    WWTP_within_distance = WWTP_within[WWTP_within['linear_distance'] <= distance]
    WWTP_within_distance = WWTP_within_distance.groupby('site_id').sum('Influent Flow (MMGal/d)')
    WWTP_within_distance = WWTP_within_distance[['Influent Flow (MMGal/d)']]
    WWTP_within_distance = WWTP_within_distance.rename(columns={'Influent Flow (MMGal/d)': int(distance)})
    WWTP_within_distance.reset_index(inplace=True)
    if len(WWTP_within_distance) > 0:
        result = result.merge(WWTP_within_distance, how='left', on='site_id')

result = result.fillna(0)

if len(result) < len(refinery):
    result_index = pd.Index(result.site_id)
    refinery_index = pd.Index(refinery.site_id)
    refineries_left_id = refinery_index.difference(result_index).values

result = result.set_index('site_id')
for item in refineries_left_id:
    result.loc[item] = [0]*len(result.columns)

result = result.merge(refinery, how='left', on='site_id')

result.to_excel(folder + f'results/MGD_vs_distance_{date.today()}.xlsx')

#%%

# this part is the analysis for CO2 abatement cost
# run analysis for every row

filterwarnings('ignore')

# remember to use the correct file
final_WWTPs = pd.read_excel(folder + 'HTL_geospatial_model_input_final.xlsx')
print(len(final_WWTPs))

facility = []
city = []
state = []
CO2_reduction = []
sludge_CO2_reduction_ratio = []
WRRF_CO2_reduction_ratio = []
USD_decarbonization = []
oil_BPD = []
    
for i in range(14100, len(final_WWTPs)):
# for i in [0, 76, 718]:
        
    sys, barrel = create_spatial_system(waste_price=400,
                                        waste_GHG=final_WWTPs.iloc[i]['sludge_management_kg_per_ton'],
                                        size=final_WWTPs.iloc[i]['Influent Flow (MMGal/d)'],
                                        distance=final_WWTPs.iloc[i]['real_distance_km'],
                                        solid_fate=final_WWTPs.iloc[i]['solid_fate'],
                                        ww_2_dry_sludge_ratio=final_WWTPs.iloc[i]['BASELINE:SOLIDS (dry kg/y):Disposed (Biosolids)']/1000/365/final_WWTPs.iloc[i]['Influent Flow (MMGal/d)'],
                                        # ww_2_dry_sludge_ratio: how much metric ton/day sludge can be produced by 1 MGD of ww
                                        state=final_WWTPs.iloc[i]['state'],
                                        elec_GHG=final_WWTPs.iloc[i]['GHG (10-year median)'])
    
    lca = sys.LCA
    
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    WWTP = unit.WWTP
    raw_wastewater = stream.raw_wastewater

    kg_CO2_per_ton_dry_sludge = lca.get_total_impacts(exclude=(raw_wastewater,))['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)/lca.lifetime

    CO2_reduction_result = final_WWTPs.iloc[i]['BASELINE:SOLIDS (dry kg/y):Disposed (Biosolids)']/1000*30*(final_WWTPs.iloc[i]['sludge_management_kg_per_ton']-kg_CO2_per_ton_dry_sludge)
    
    sludge_CO2_reduction_ratio_result = CO2_reduction_result/final_WWTPs.iloc[i]['sludge_management_kg_CO2_per_day_AD_included']/365/30
    
    WRRF_CO2_reduction_ratio_result = CO2_reduction_result/final_WWTPs.iloc[i]['30_years_emission_ton_CO2']/1000

    # make sure CO2_reduction_result is positive if you want to calculate USD_per_ton_CO2_reduction
    
    if CO2_reduction_result > 0:
        USD_per_ton_CO2_reduction = -sys.TEA.NPV/CO2_reduction_result*1000 # save money: the results are negative; spend money: the results are positive
    else:
        USD_per_ton_CO2_reduction = np.nan
    
    facility.append(final_WWTPs.iloc[i]['FACILITY'])
    city.append(final_WWTPs.iloc[i]['CITY_x'])
    state.append(final_WWTPs.iloc[i]['ST'])
    CO2_reduction.append(CO2_reduction_result)
    sludge_CO2_reduction_ratio.append(sludge_CO2_reduction_ratio_result)
    WRRF_CO2_reduction_ratio.append(WRRF_CO2_reduction_ratio_result)
    USD_decarbonization.append(USD_per_ton_CO2_reduction)
    oil_BPD.append(barrel)
    
    if USD_per_ton_CO2_reduction < 0:
        print('HERE WE GO!')
    
    if i%5 == 0:
        print(i) # check progress
        
result = {'facility': facility,
          'city': city,
          'state': state,
          'CO2_reduction': CO2_reduction,
          'sludge_CO2_reduction_ratio': sludge_CO2_reduction_ratio,
          'WRRF_CO2_reduction_ratio': WRRF_CO2_reduction_ratio,
          'USD_decarbonization': USD_decarbonization,
          'oil_BPD': oil_BPD}
        
result = pd.DataFrame(result)

result.to_excel(folder + f'results/decarbonization_{date.today()}_{i}.xlsx')

#%%

# this part is to make plant-level decarbonization cost vs decarbonization ratio figure

decarbonization_result = pd.read_excel(folder + 'results/decarbonization_results_10_7_2023/carbonization_summary.xlsx')

WWTP_within = pd.read_excel(folder + 'HTL_geospatial_model_input_final.xlsx')

decarbonization_result = decarbonization_result.merge(WWTP_within[['FACILITY','CITY_x','Influent Flow (MMGal/d)', 'category']], how='left', left_on=['facility','city'], right_on=['FACILITY','CITY_x'])

decarbonization_result.loc[decarbonization_result['category'].isin((1, 2, 4)) ,'AD'] = b # blue

decarbonization_result.loc[~decarbonization_result['category'].isin((1, 2, 4)) ,'AD'] = r # red

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 500]

decarbonization_result_more = decarbonization_result[decarbonization_result['USD_decarbonization'] > -2000]

decarbonization_result_less = decarbonization_result[decarbonization_result['USD_decarbonization'] < -2000]

fig, ax_more = plt.subplots(figsize = (15, 10))

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.xticks(fontname = 'Arial')
plt.yticks(fontname = 'Arial')

ax_more = plt.gca()
ax_more.set_xlim([-3, 53])
ax_more.set_ylim([-1625, 625])
ax_more.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
ax_more.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False)

ax_more_right = ax_more.twinx()
ax_more_right.set_ylim(ax_more.get_ylim())
ax_more_right.tick_params(direction='in', length=15, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_more_top = ax_more.twiny()
ax_more_top.set_xlim(ax_more.get_xlim())
ax_more_top.tick_params(direction='in', length=15, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 55, 5))

ax_more.scatter(x = decarbonization_result_more['WRRF_CO2_reduction_ratio']*100,
                y = decarbonization_result_more['USD_decarbonization'],
                s = decarbonization_result_more['Influent Flow (MMGal/d)']*4,
                c = decarbonization_result_more['AD'],
                linewidths = 0,
                alpha = 0.7)

ax_more.scatter(x = decarbonization_result_more['WRRF_CO2_reduction_ratio']*100,
                y = decarbonization_result_more['USD_decarbonization'],
                s = decarbonization_result_more['Influent Flow (MMGal/d)']*4,
                color = 'none',
                linewidths = 2,
                edgecolors = 'k')

ax_less = fig.add_axes([0.5, 0.2, 0.36, 0.32])   

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.xticks(fontname = 'Arial')
plt.yticks(fontname = 'Arial')

ax_less.set_xlim([-0.5, 3.5])
ax_less.set_ylim([-40000, 0])
ax_less.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
ax_less.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False)

ax_less_right = ax_less.twinx()
ax_less_right.set_ylim(ax_less.get_ylim())
plt.yticks(np.arange(-35000, 5000, 10000))
ax_less_right.tick_params(direction='in', length=15, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_less_top = ax_less.twiny()
ax_less_top.set_xlim(ax_less.get_xlim())
ax_less_top.tick_params(direction='in', length=15, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-35000, 5000, 10000))

ax_less.scatter(x = decarbonization_result_less['WRRF_CO2_reduction_ratio']*100,
            y = decarbonization_result_less['USD_decarbonization'],
            s = decarbonization_result_less['Influent Flow (MMGal/d)']*4,
            c = decarbonization_result_less['AD'],
            linewidths = 0,
            alpha = 0.7)

ax_less.scatter(x = decarbonization_result_less['WRRF_CO2_reduction_ratio']*100,
            y = decarbonization_result_less['USD_decarbonization'],
            s = decarbonization_result_less['Influent Flow (MMGal/d)']*4,
            color = 'none',
            linewidths = 2,
            edgecolors = 'k')