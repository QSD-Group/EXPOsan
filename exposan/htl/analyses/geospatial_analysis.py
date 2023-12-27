#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Sun Jun 11 08:12:41 2023

@author: jiananfeng
'''

import geopy.distance, googlemaps
import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt, matplotlib.ticker as mtick
from exposan.htl.geospatial_HTL_systems import create_geospatial_system
from qsdsan.utils import palettes
from datetime import date
from warnings import filterwarnings

folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'

# color palette
Guest = palettes['Guest']
b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX
a = Guest.gray.HEX

def set_plot(figure_size=(30, 30)):
    global fig, ax
    fig, ax = plt.subplots(figsize=figure_size)
    ax.tick_params(top=False, bottom=False, left=False, right=False,
                   labelleft=False, labelbottom=False)
    ax.set_frame_on(False)

#%% read data and data pre-processing

WRRF = pd.read_excel(folder + 'AMO_results_12262023.xlsx')

WRRF['total_sludge_amount_kg_per_year'] = WRRF[['landfill','land_application',
                                                'incineration','other_management']].sum(axis=1)

# all emission data in kg CO2 eq/day
WRRF['total_emission'] = WRRF[['CO2_Emission','CH4_Emission','N2O_Emission',
                               'E_Emission','NG_Emission','NC_Emission',
                               'Biosolids_Emission']].sum(axis=1)

treatment_trains = np.array(['B1','B1E','B2','B3','B4','B5','B6',
                             'C1','C2','C3','C5','C6','D1','D2',
                             'D3','D5','D6','E2','E2P','F1','G1',
                             'G1E','G2','G3','G5','G6','H1','I1',
                             'I1E','I2','I3','I5','I6','N1','N2',
                             'O1','O1E','O2','O3','O5','O6',
                             'LAGOON_AER','LAGOON_ANAER','LAGOON_FAC'], dtype=object)

TT_indentifier = WRRF[treatment_trains].apply(lambda x: x > 0)
WRRF['treatment_train'] = TT_indentifier.apply(lambda x: list(treatment_trains[x.values]), axis=1)

WRRF = WRRF[['FACILITY','CITY','STATE','CWNS_NUM','FACILITY_CODE','LATITUDE',
             'LONGITUDE','FLOW_2022_MGD','treatment_train',
             'sludge_anaerobic_digestion','sludge_aerobic_digestion','landfill',
             'land_application','incineration','other_management',
             'total_sludge_amount_kg_per_year','CH4_Emission','N2O_Emission',
             'CO2_Emission','Biosolids_Emission','E_Emission','NG_Emission',
             'NC_Emission','total_emission']]

WRRF = WRRF.rename({'FACILITY':'facility',
                    'CITY':'city',
                    'STATE':'state',
                    'CWNS_NUM':'CWNS',
                    'FACILITY_CODE':'facility_code',
                    'LATITUDE':'latitude',
                    'LONGITUDE':'longitude',
                    'FLOW_2022_MGD':'flow_2022_MGD',
                    'other_management':'other_sludge_management',
                    'CH4_Emission':'CH4_emission',
                    'N2O_Emission':'N2O_emission',
                    'CO2_Emission':'CO2_emission',          
                    'Biosolids_Emission':'biosolids_emission',
                    'E_Emission':'electricity_emission',
                    'NG_Emission':'natural_gas_emission',
                    'NC_Emission':'non_combustion_CO2_emission'}, axis=1)

WRRF = WRRF.sort_values(by='flow_2022_MGD', ascending=False)

WRRF = gpd.GeoDataFrame(WRRF, crs='EPSG:4269',
                        geometry=gpd.points_from_xy(x=WRRF.longitude,
                                                    y=WRRF.latitude))

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
    
WRRF = WRRF.to_crs(crs='EPSG:3857')
refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')

# NOTE 1: if only select states, uncomment the next line of code
# US = US.loc[US['NAME'].isin(('Pennsylvania',))]

WRRF = gpd.sjoin(WRRF, US)
WRRF = WRRF.drop(['index_right'], axis=1)
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

electricity = pd.merge(elec, US, left_on='name', right_on='NAME')
electricity = gpd.GeoDataFrame(electricity)

#%% WRRFs visualization

set_plot()

# !!! for plot in SI, can use default linewidth for the US map
US.plot(ax=ax, color='w', edgecolor='k', linewidth=5)

WRRF = WRRF.sort_values(by='flow_2022_MGD', ascending=False)

WRRF_flow = WRRF['flow_2022_MGD']

# we use 'markersize' here (in .plot()) but 's' later (in .scatter())
# they are not the same
# markersize is proportional to length
# s is proportional to area
# see https://stackoverflow.com/questions/14827650/pyplot-scatter-plot-marker-size

less_than_10 = WRRF_flow <= 10
WRRF[less_than_10].plot(ax=ax, color=Guest.gray.HEX, markersize=WRRF.loc[less_than_10, WRRF_flow.name]**0.5*25, alpha=0.5)

between_10_and_100 = (WRRF_flow > 10) & (WRRF_flow <= 100)
WRRF[between_10_and_100].plot(ax=ax, color=Guest.red.HEX, markersize=WRRF.loc[between_10_and_100, WRRF_flow.name]**0.5*25, edgecolor='k', linewidth=1.5, alpha=0.7)

more_than_100 = WRRF_flow > 100
WRRF[more_than_100].plot(ax=ax, color=Guest.green.HEX, markersize=WRRF.loc[more_than_100, WRRF_flow.name]**0.5*25, edgecolor='k', linewidth=2, alpha=0.9)

# if you want to show all WRRFs together with the same symbols, comment the codes above and uncomment the next line of code
# sludge.plot(ax=ax, color=Guest.gray.HEX, markersize=10)

#%% oil refinery visualization

set_plot()

US.plot(ax=ax,
        color=[r, b, g, b, b, r, g, b, o, b, g, y, r, g, r, y, r,
               b, b, g, o, o, g, o, b, g, y, g, b, o, g, b, b, y,
               b, b, b, b, b, b, g, g, g, y, g, r, g, g, b],
        alpha=0.4,
        edgecolor='k',
        linewidth=0)

# !!! for plot in SI, can use default linewidth for the US map
US.plot(ax=ax, color='none', edgecolor='k', linewidth=5)

refinery['total_capacity'] = refinery[[i for i in refinery.columns if i[-4:] == 'Mbpd']].sum(axis=1)

refinery = refinery.sort_values(by='total_capacity', ascending=False)

refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=refinery['total_capacity']**0.5*25, edgecolor='k', linewidth=2, alpha=0.9)

#%% WRRFs+oil refineries visualization (just select states, search 'NOTE 1')

if len(set(US['NAME'])) == 1:

    set_plot()
    
    US.plot(ax=ax, color=Guest.blue.HEX, alpha=0.4, edgecolor='k', linewidth=0)
    US.plot(ax=ax, color='none', edgecolor='k', linewidth=10)
    
    WRRF.plot(ax=ax, color=Guest.gray.HEX, markersize=200)
    
    refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=5000, edgecolor='k', linewidth=10)

#%% WRRFs+oil refineries visualization (all)

set_plot()

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3.5)
WRRF.plot(ax=ax, color=Guest.gray.HEX, markersize=5)
refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=500, edgecolor='k', linewidth=3.5)

#%% electricity price and carbon intensity visualization

set_plot()

electricity.plot('price (10-year median)', ax=ax, cmap='Oranges', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35})

# electricity.plot('GHG (10-year median)', ax=ax, cmap='Blues', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35})

#%% transporation distance calculation

# if want select WRRFs that are within certain ranges of refineries, set max_distance in sjoin_nearest.
# here we do not set a max distance

WRRF_input = WRRF.sjoin_nearest(refinery, max_distance=None, distance_col='distance')

# if set max_distance = 100000, we have a problem here, that is CRS 3857 is not accurate at all,
# especially more far away from the equator, therefore we need use a larger distance here,
# for example, 150000 (tested, when use 150000, all actual distances are smaller than 100000)
# as the max_distance and then use geopy.distance.geodesic to recalculate the distance and
# filter out those that are actually longer than 100000

WRRF_input['WRRF_location'] = list(zip(WRRF_input.latitude, WRRF_input.longitude))
WRRF_input['refinery_location'] = list(zip(WRRF_input.Latitude, WRRF_input.Longitude))

linear_distance = []
for i in range(len(WRRF_input)):
    linear_distance.append(geopy.distance.geodesic(WRRF_input['WRRF_location'].iloc[i], WRRF_input['refinery_location'].iloc[i]).km)

WRRF_input['linear_distance_km'] = linear_distance

# if we set the max_distance, then uncomment the next line (replace 100 km when necessary)
# WRRF_input = WRRF_input[WRRF_input['linear_distance'] <= 100]

# we have previously calculated the distance using Google Maps API and the old dataset (based on Seiple et al.),
# now, we are using the AMO dataset, and to avoid do this again (which may result in a warning from Google...),
# we first merge these two datasets and only calculate distances that have not been calculated before.
# note there are 9 WRRFs in the old datasets whose locations were adjusted to enable distance calculation,
# they will be removed after getting the new HTL geospatial model input file.
# when calculate new distances, if there are WRRFs that cannot get a distance, remove them as well.

old_dataset_based_on_Seiple = pd.read_excel(folder + 'HTL_geospatial_model_input_final_old.xlsx')

old_dataset_based_on_Seiple = old_dataset_based_on_Seiple[['UID','FACILITY_CODE','site_id','real_distance_km']]

assert old_dataset_based_on_Seiple[['UID','FACILITY_CODE','site_id']].duplicated().sum() == 0

assert WRRF_input[['CWNS','facility_code','site_id']].duplicated().sum() == 0

WRRF_input = WRRF_input.merge(old_dataset_based_on_Seiple, how='left', left_on=['CWNS','facility_code','site_id'], right_on=['UID','FACILITY_CODE','site_id'])

WRRF_input = WRRF_input.drop(['UID','FACILITY_CODE'], axis=1)

WRRF_input_without_distance = WRRF_input[WRRF_input['real_distance_km'].isna()]

# !!! run <5000 datapoints every time !!! get a google API key !!! do not upload to GitHub
gmaps = googlemaps.Client(key='XXX')
real_distance = []
for i in WRRF_input_without_distance.index:
    try:
        print(i)
        WRRF_input['real_distance_km'].iloc[i] = gmaps.distance_matrix(WRRF_input['WRRF_location'].iloc[i],
                                                                       WRRF_input['refinery_location'].iloc[i],
                                                                       mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000
    except KeyError:
        print('--------------------------------')
        WRRF_input['real_distance_km'].iloc[i] = np.nan

# remove WRRFs whose distance to oil refineries cannot be calculated using Google Maps API
# there are 2 WRRFs that are removed here, they are:
# AVALON WWRF (0.915 MGD)
# MID-VALLEY WRP NO.4 (5.5 MGD)
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# there is a total of 11 WRRFs that are removed
# 9 of them were included in the old datasets and were manually removed, they are:
# RICKARDSVILLE WWTP (0.025 MGD)
# MACKINAC ISLAND STP (0.895 MGD)
# Ewell-Rhodes Point WWTP (0.0425 MGD)
# TYLERTON WWTP (0.013 MGD)
# TANGIER STP (0.075 MGD)
# WARM SPRINGS REHAB WWTF (0.05 MGD)
# BIG PINE WWTF (0.0545 MGD)
# NORTH HAVEN, WWTF (0.04 MGD)
# Eagle Nest, Village of (0.09 MGD)
WRRF_input = WRRF_input[~WRRF_input['facility'].isin(['RICKARDSVILLE WWTP','MACKINAC ISLAND STP',
                                                      'Ewell-Rhodes Point WWTP','TYLERTON WWTP',
                                                      'TANGIER STP','WARM SPRINGS REHAB WWTF',
                                                      'BIG PINE WWTF','NORTH HAVEN, WWTF',
                                                      'Eagle Nest, Village of'])]

WRRF_input.to_excel(folder + f'HTL_geospatial_model_input_{date.today()}.xlsx') # this will be the input for the future analysis

#%% travel distance box plot # !!! need to run 2 times to make the label font and size right

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

fig, ax = plt.subplots(figsize = (5, 8))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.yticks(fontname = 'Arial')

ax = plt.gca()
ax.set_ylim([-50, 850])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

bp = plt.boxplot(WRRF_input['real_distance_km'], showfliers=False, widths=0.5, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=b, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

#%% WRRFs GHG map

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

WRRF_input = WRRF_input.sort_values(by='total_emission', ascending=False)

WRRF_input = gpd.GeoDataFrame(WRRF_input, crs='EPSG:4269',
                               geometry=gpd.points_from_xy(x=WRRF_input.longitude,
                                                           y=WRRF_input.latitude))

WRRF_input = WRRF_input.to_crs(crs='EPSG:3857')

set_plot()

US.plot(ax=ax, color='w', edgecolor='k', linewidth=5)

# use tonne/day

WRRF_GHG_tonne_per_day = WRRF_input['total_emission']/1000

less_than_50 = WRRF_GHG_tonne_per_day <= 50
WRRF_input[less_than_50].plot(ax=ax, color=Guest.gray.HEX, markersize=WRRF_input.loc[less_than_50, WRRF_GHG_tonne_per_day.name]**0.5/3, alpha=0.5)

between_50_and_500 = (WRRF_GHG_tonne_per_day > 50) & (WRRF_GHG_tonne_per_day <= 500)
WRRF_input[between_50_and_500].plot(ax=ax, color=Guest.red.HEX, markersize=WRRF_input.loc[between_50_and_500, WRRF_GHG_tonne_per_day.name]**0.5/3, edgecolor='k', linewidth=1.5, alpha=0.7)

more_than_500 = WRRF_GHG_tonne_per_day > 500
WRRF_input[more_than_500].plot(ax=ax, color=Guest.green.HEX, markersize=WRRF_input.loc[more_than_500, WRRF_GHG_tonne_per_day.name]**0.5/3, edgecolor='k', linewidth=2, alpha=0.9)

#%% # WRRFs sludge management GHG map

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

WRRF_input = WRRF_input.sort_values(by='biosolids_emission', ascending=False)

WRRF_input = gpd.GeoDataFrame(WRRF_input, crs='EPSG:4269',
                               geometry=gpd.points_from_xy(x=WRRF_input.longitude,
                                                           y=WRRF_input.latitude))

WRRF_input = WRRF_input.to_crs(crs='EPSG:3857')

set_plot()

US.plot(ax=ax, color='w', edgecolor='k', linewidth=5)

# use tonne/day

sludge_GHG_tonne_per_day = WRRF_input['biosolids_emission']/1000

less_than_50 = sludge_GHG_tonne_per_day <= 50
WRRF_input[less_than_50].plot(ax=ax, color=Guest.gray.HEX, markersize=WRRF_input.loc[less_than_50, sludge_GHG_tonne_per_day.name]**0.5/3, alpha=0.5)

between_50_and_500 = (sludge_GHG_tonne_per_day > 50) & (sludge_GHG_tonne_per_day <= 500)
WRRF_input[between_50_and_500].plot(ax=ax, color=Guest.red.HEX, markersize=WRRF_input.loc[between_50_and_500, sludge_GHG_tonne_per_day.name]**0.5/3, edgecolor='k', linewidth=1.5, alpha=0.7)

more_than_500 = sludge_GHG_tonne_per_day > 500
WRRF_input[more_than_500].plot(ax=ax, color=Guest.green.HEX, markersize=WRRF_input.loc[more_than_500, sludge_GHG_tonne_per_day.name]**0.5/3, edgecolor='k', linewidth=2, alpha=0.9)

#%% cumulative WRRFs capacity vs distances

# remember to use the correct file
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

result = WRRF_input[['site_id']].drop_duplicates()

max_distance = 1000
for distance in np.linspace(0, max_distance, max_distance+1):
    WRRF_input_distance = WRRF_input[WRRF_input['linear_distance_km'] <= distance]
    WRRF_input_distance = WRRF_input_distance.groupby('site_id').sum('flow_2022_MGD')
    WRRF_input_distance = WRRF_input_distance[['flow_2022_MGD']]
    WRRF_input_distance = WRRF_input_distance.rename(columns={'flow_2022_MGD': int(distance)})
    WRRF_input_distance.reset_index(inplace=True)
    if len(WRRF_input_distance) > 0:
        result = result.merge(WRRF_input_distance, how='left', on='site_id')

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

#%% make the plot of cumulative WRRFs capacity vs distances (data preparation)

# import file for the cumulative figure (CF)
CF_input = pd.read_excel(folder + 'results/MGD_vs_distance_2023-12-26.xlsx')

CF_input[0] = 0

CF_input = CF_input[[0, *CF_input.columns[2:-1]]]

CF_input = CF_input.iloc[:, 0:max_distance+6]

CF_input.drop(['Company','Corp','Site'], axis=1, inplace=True)

CF_input.sort_values(by='PADD', inplace=True)

CF_input = CF_input[['PADD','State', *CF_input.columns[0:-2]]]

CF_input = CF_input.transpose()

#%% make the plot of cumulative WRRFs capacity vs distances (separated WRRF) # !!! need to run 2 times to make the label font and size right # TODO: not sure whether this problem happens to facility-level and regional level plots

fig = plt.figure(figsize=(20, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

PADD_1 = (CF_input.loc['PADD',:].isin([1,])).sum()

PADD_2 = (CF_input.loc['PADD',:].isin([1,2])).sum()

PADD_3 = (CF_input.loc['PADD',:].isin([1,2,3])).sum()

PADD_4 = (CF_input.loc['PADD',:].isin([1,2,3,4])).sum()

PADD_5 = (CF_input.loc['PADD',:].isin([1,2,3,4,5])).sum()

def set_cf_plot():
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    plt.xticks(fontname = 'Arial')
    plt.yticks(fontname = 'Arial')

def add_region(position, left_label, start_region, end_region, color):
    ax = fig.add_subplot(gs[0, position])
    
    set_cf_plot()
    
    for i in range(start_region, end_region):
        CF_input.iloc[2:, i].plot(ax=ax, color=color, linewidth=2)
        
    ax.set_xlim([0, max_distance])
    ax.set_ylim([0, 8000])
    
    if position == 4:
        ax.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=False, right=False, labelleft=left_label, labelbottom=True)
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        plt.yticks(np.arange(0, 8000, 1000))
        ax_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
    if position == 0:
        ax.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False, labelleft=left_label, labelbottom=True)
        plt.xticks(np.arange(0, max_distance*1.2, max_distance*0.2))
    else:
        ax.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=False, right=False, labelleft=left_label, labelbottom=True)
        plt.xticks(np.arange(max_distance*0.2, max_distance*1.2, max_distance*0.2))
    
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    plt.xticks(np.arange(max_distance*0.2, max_distance*1.2, max_distance*0.2))
    ax_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=False, labelcolor='none')

add_region(0, True, 0, PADD_1, b)
add_region(1, False, PADD_1, PADD_2, g)
add_region(2, False, PADD_2, PADD_3, r)
add_region(3, False, PADD_3, PADD_4, o)
add_region(4, False, PADD_4, PADD_5, y)
#%% make the plot of cumulative WRRFs capacity vs distances (regional total) # !!! need to run 2 times to make the label font and size right # TODO: not sure whether this problem happens to facility-level and regional level plots

fig, ax = plt.subplots(figsize=(5, 10))

set_cf_plot()

CF_input.iloc[2:, 0:PADD_1].sum(axis=1).plot(ax=ax, color=b, linewidth=2)

CF_input.iloc[2:, PADD_1:PADD_2].sum(axis=1).plot(ax=ax, color=g, linewidth=2)

CF_input.iloc[2:, PADD_2:PADD_3].sum(axis=1).plot(ax=ax, color=r, linewidth=2)

CF_input.iloc[2:, PADD_3:PADD_4].sum(axis=1).plot(ax=ax, color=o, linewidth=2)

CF_input.iloc[2:, PADD_4:PADD_5].sum(axis=1).plot(ax=ax, color=y, linewidth=2)

ax.set_xlim([0, max_distance])
ax.set_ylim([0, 20000])

ax.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
plt.yticks(np.arange(0, 20000, 1000))
ax_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(max_distance*0.2, max_distance*1.2, max_distance*0.2))
ax_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=False, labelcolor='none')

#%% CO2 abatement cost analysis

filterwarnings('ignore')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

print(len(WRRF_input))

# sludge disposal cost in $/kg sludge
# assume the cost for other_sludge_management is the weighted average of the other three, ratio based on Peccia and Westerhoff. 2015. https://pubs.acs.org/doi/full/10.1021/acs.est.5b01931
sludge_disposal_cost = {'landfill': 0.413,
                        'land_application': 0.606,
                        'incineration': 0.441,
                        'other_sludge_management': 0.413*0.3 + 0.606*0.55 + 0.441*0.15}

# sludge emission factor in kg CO2 eq/kg sludge
# assume the GHG for other_sludge_management is the weighted average of the other three, ratio based on Peccia and Westerhoff. 2015. https://pubs.acs.org/doi/full/10.1021/acs.est.5b01931
sludge_emission_factor = {'landfill': 1.36,
                          'land_application': 0.353,
                          'incineration': 0.519,
                          'other_sludge_management': 1.36*0.3 + 0.353*0.55 + 0.519*0.15}

WRRF_input['waste_cost'] = sum(WRRF_input[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

WRRF_input['waste_GHG'] =  sum(WRRF_input[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

CWNS = []
facility_code = []
CO2_reduction = []
sludge_CO2_reduction_ratio = []
WRRF_CO2_reduction_ratio = []
NPV = []
USD_decarbonization = []
oil_BPD = []
    
for i in range(10000, len(WRRF_input)): # !!! run in different consoles to speed up
# for i in range(0, 2):
    
    sys, barrel = create_geospatial_system(waste_cost=WRRF_input.iloc[i]['waste_cost'],
                                           waste_GHG=WRRF_input.iloc[i]['waste_GHG'],
                                           size=WRRF_input.iloc[i]['flow_2022_MGD'],
                                           distance=WRRF_input.iloc[i]['real_distance_km'],
                                           anaerobic_digestion=WRRF_input.iloc[i]['sludge_anaerobic_digestion'],
                                           aerobic_digestion=WRRF_input.iloc[i]['sludge_aerobic_digestion'],
                                           ww_2_dry_sludge_ratio=WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']/1000/365/WRRF_input.iloc[i]['flow_2022_MGD'],
                                           # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                           state=WRRF_input.iloc[i]['state'],
                                           elec_GHG=float(elec[elec['state']==WRRF_input.iloc[i]['state']]['GHG (10-year median)']))
    
    lca = sys.LCA
    
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    WRRF = unit.WWTP
    raw_wastewater = stream.feedstock_assumed_in_wastewater
    
    CO2_reduction_result = -lca.get_total_impacts()['GlobalWarming']
  
    sludge_CO2_reduction_ratio_result = CO2_reduction_result/WRRF_input.iloc[i]['biosolids_emission']/365/30
    
    WRRF_CO2_reduction_ratio_result = CO2_reduction_result/WRRF_input.iloc[i]['total_emission']/365/30

    # make sure CO2_reduction_result is positive if you want to calculate USD_per_tonne_CO2_reduction
    
    if CO2_reduction_result > 0:
        USD_per_tonne_CO2_reduction = -sys.TEA.NPV/CO2_reduction_result*1000 # 2020$/tonne, save money: the results are negative; spend money: the results are positive
    else:
        USD_per_tonne_CO2_reduction = np.nan
    
    CWNS.append(WRRF_input.iloc[i]['CWNS'])
    facility_code.append(WRRF_input.iloc[i]['facility_code'])
    CO2_reduction.append(CO2_reduction_result)
    sludge_CO2_reduction_ratio.append(sludge_CO2_reduction_ratio_result)
    WRRF_CO2_reduction_ratio.append(WRRF_CO2_reduction_ratio_result)
    NPV.append(sys.TEA.NPV)
    USD_decarbonization.append(USD_per_tonne_CO2_reduction)
    oil_BPD.append(barrel)
    
    if USD_per_tonne_CO2_reduction < 0:
        print('HERE WE GO!')
    
    if i%5 == 0:
        print(i) # check progress
    
result = {'CWNS': CWNS,
          'facility_code': facility_code,
          'CO2_reduction': CO2_reduction,
          'sludge_CO2_reduction_ratio': sludge_CO2_reduction_ratio,
          'WRRF_CO2_reduction_ratio': WRRF_CO2_reduction_ratio,
          'NPV': NPV,
          'USD_decarbonization': USD_decarbonization,
          'oil_BPD': oil_BPD}
        
result = pd.DataFrame(result)

result.to_excel(folder + f'results/decarbonization_{date.today()}_{i}.xlsx')

#%% merge the results and the input

input_data = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

output_result = pd.read_excel(folder + 'results/decarbonization_result_2023-12-26.xlsx')

assert input_data[['CWNS','facility_code']].duplicated().sum() == 0

assert output_result[['CWNS','facility_code']].duplicated().sum() == 0

integrated_result = input_data.merge(output_result, how='left', on=['CWNS','facility_code'])

integrated_result.to_excel(folder + f'results/integrated_decarbonization_result_{date.today()}.xlsx')

#%% facility-level decarbonization cost vs decarbonization ratio figure

# TODO: pick up from here after getting new results

decarbonization_result = pd.read_excel(folder + 'results/decarbonization_results_10_7_2023/carbonization_summary.xlsx')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

# TODO: may remove the merge since we may add more parameters eailier (line 388)
decarbonization_result = decarbonization_result.merge(WRRF_input[['FACILITY','CITY_x','flow_2022_MGD','category']], how='left', left_on=['facility','city'], right_on=['FACILITY','CITY_x']) # TODO: change to use 2022 MGD data

decarbonization_result.loc[decarbonization_result['category'].isin((1, 2, 4)) ,'AD'] = b # blue # TODO: update based on AMO data instead of Seiple data

decarbonization_result.loc[~decarbonization_result['category'].isin((1, 2, 4)) ,'AD'] = r # red # TODO: update based on AMO data instead of Seiple data

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 500] # TODO: decide the upper threshold, 0, 100, or other values?

decarbonization_result_more = decarbonization_result[decarbonization_result['USD_decarbonization'] > -2000]

decarbonization_result_less = decarbonization_result[decarbonization_result['USD_decarbonization'] < -2000] # TODO: discuss with Jeremy, how to deal with these data (very negative decarbonization cost due to low decarbonization ratio)

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

plt.xticks(np.arange(0, 55, 5))

ax_more_right = ax_more.twinx()
ax_more_right.set_ylim(ax_more.get_ylim())
ax_more_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_more_top = ax_more.twiny()
ax_more_top.set_xlim(ax_more.get_xlim())
plt.xticks(np.arange(0, 55, 5))
ax_more_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_more.scatter(x = decarbonization_result_more['WRRF_CO2_reduction_ratio']*100,
                y = decarbonization_result_more['USD_decarbonization'],
                s = decarbonization_result_more['flow_2022_MGD']*4, # TODO: change to use 2022 MGD data
                c = decarbonization_result_more['AD'],
                linewidths = 0,
                alpha = 0.7)

ax_more.scatter(x = decarbonization_result_more['WRRF_CO2_reduction_ratio']*100,
                y = decarbonization_result_more['USD_decarbonization'],
                s = decarbonization_result_more['flow_2022_MGD']*4, # TODO: change to use 2022 MGD data
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

plt.yticks(np.arange(-35000, 5000, 10000))

ax_less_right = ax_less.twinx()
ax_less_right.set_ylim(ax_less.get_ylim())
plt.yticks(np.arange(-35000, 5000, 10000))
ax_less_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_less_top = ax_less.twiny()
ax_less_top.set_xlim(ax_less.get_xlim())
ax_less_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_less.scatter(x = decarbonization_result_less['WRRF_CO2_reduction_ratio']*100,
                y = decarbonization_result_less['USD_decarbonization'],
                s = decarbonization_result_less['flow_2022_MGD']*4, # TODO: change to use 2022 MGD data
                c = decarbonization_result_less['AD'],
                linewidths = 0,
                alpha = 0.7)

ax_less.scatter(x = decarbonization_result_less['WRRF_CO2_reduction_ratio']*100,
                y = decarbonization_result_less['USD_decarbonization'],
                s = decarbonization_result_less['flow_2022_MGD']*4, # TODO: change to use 2022 MGD data 
                color = 'none',
                linewidths = 2,
                edgecolors = 'k')

#%%

# regional level analysis

decarbonization_result = pd.read_excel(folder + 'results/decarbonization_results_10_7_2023/carbonization_summary.xlsx')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

# TODO: may remove the merge since we may add more parameters eailier (line 388)
decarbonization_result = decarbonization_result.merge(WRRF_input[['FACILITY','CITY_x','PADD']], how='left', left_on=['facility','city'], right_on=['FACILITY','CITY_x'])

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 1000] # TODO: decide the upper threshold, 0, 100, or other values? # TODO: discuss with Jeremy, how to deal with very negative decarbonization cost (due to low decarbonization ratio)

PADD_1 = decarbonization_result[decarbonization_result['PADD'] == 1]
PADD_1.sort_values(by='USD_decarbonization', inplace=True)
PADD_1['cummulative_oil_BPD'] = PADD_1['oil_BPD'].cumsum()

PADD_2 = decarbonization_result[decarbonization_result['PADD'] == 2]
PADD_2.sort_values(by='USD_decarbonization', inplace=True)
PADD_2['cummulative_oil_BPD'] = PADD_2['oil_BPD'].cumsum()

PADD_3 = decarbonization_result[decarbonization_result['PADD'] == 3]
PADD_3.sort_values(by='USD_decarbonization', inplace=True)
PADD_3['cummulative_oil_BPD'] = PADD_3['oil_BPD'].cumsum()

PADD_4 = decarbonization_result[decarbonization_result['PADD'] == 4]
PADD_4.sort_values(by='USD_decarbonization', inplace=True)
PADD_4['cummulative_oil_BPD'] = PADD_4['oil_BPD'].cumsum()

PADD_5 = decarbonization_result[decarbonization_result['PADD'] == 5]
PADD_5.sort_values(by='USD_decarbonization', inplace=True)
PADD_5['cummulative_oil_BPD'] = PADD_5['oil_BPD'].cumsum()

fig, ax_more = plt.subplots(figsize = (15, 10))

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.xticks(fontname = 'Arial')
plt.yticks(fontname = 'Arial')

ax_more.set_xlim([0, 5000])
ax_more.set_ylim([-1500, 1000])

plt.axvspan(xmin=0, xmax=5000, ymin=0, ymax=0.6, facecolor=a, alpha=0.3)
plt.axvspan(xmin=0, xmax=5000, ymin=0, ymax=0.6, facecolor='none', linewidth=2, edgecolor='k')

def plot_line(ax, data, color):
    ax.plot((0, *data['cummulative_oil_BPD']), (0, *data['USD_decarbonization']),
             color=color, marker='x', markersize=10, markeredgewidth=2, linewidth=2)
    
plot_line(ax_more, PADD_1, b)
plot_line(ax_more, PADD_2, g)
plot_line(ax_more, PADD_3, r)
plot_line(ax_more, PADD_4, o)
plot_line(ax_more, PADD_5, y)

ax_more.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False)

ax_more_right = ax_more.twinx()
ax_more_right.set_ylim(ax_more.get_ylim())
ax_more_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_more_top = ax_more.twiny()
ax_more_top.set_xlim(ax_more.get_xlim())
ax_more_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_less = fig.add_axes([0.6, 0.2, 0.25, 0.3]) 

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.xticks(fontname = 'Arial')
plt.yticks(fontname = 'Arial')

ax_less.set_xlim([0, 1500])
ax_less.set_ylim([-33500, -1500])

plt.axvspan(xmin=0, xmax=1500, facecolor=a, alpha=0.3)

plot_line(ax_less, PADD_1, b)
plot_line(ax_less, PADD_2, g)
plot_line(ax_less, PADD_3, r)
plot_line(ax_less, PADD_4, o)
plot_line(ax_less, PADD_5, y)

plt.xticks(np.arange(0, 2000, 500))
plt.yticks(np.arange(-30000, 0, 5000))

ax_less.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False)

ax_less_right = ax_less.twinx()
ax_less_right.set_ylim(ax_less.get_ylim())
plt.yticks(np.arange(-30000, 0, 5000))
ax_less_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

ax_less_top = ax_less.twiny()
ax_less_top.set_xlim(ax_less.get_xlim())
plt.xticks(np.arange(0, 2000, 500))
ax_less_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')

#%%

# national level analysis

decarbonization_result = pd.read_excel(folder + 'results/decarbonization_results_10_7_2023/carbonization_summary.xlsx')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

decarbonization_result = decarbonization_result.merge(WRRF_input[['FACILITY','CITY_x','site_id',
                                                                   'AD_Mbpd','Vdist_Mbpd','CaDis_Mbpd',
                                                                   'HyCrk_Mbpd','VRedu_Mbpd','CaRef_Mbpd',
                                                                   'Isal_Mbpd','HDS_Mbpd','Cokin_Mbpd',
                                                                   'Asph_Mbpd','30_years_emission_tonne_CO2',
                                                                   'sludge_management_kg_CO2_per_day_AD_included']],
                                                      how='left', left_on=['facility','city'], right_on=['FACILITY','CITY_x'])

total_CO2_emission = decarbonization_result.sum(axis=0)['30_years_emission_tonne_CO2']

total_sludge_CO2_emission = decarbonization_result.sum(axis=0)['sludge_management_kg_CO2_per_day_AD_included']*365*30/1000

oil = decarbonization_result[['site_id','AD_Mbpd','Vdist_Mbpd','CaDis_Mbpd','HyCrk_Mbpd','VRedu_Mbpd','CaRef_Mbpd','Isal_Mbpd','HDS_Mbpd','Cokin_Mbpd','Asph_Mbpd']]

oil = oil.drop_duplicates(subset='site_id')

oil = oil.drop(labels='site_id', axis=1)

total_oil = oil.sum().sum()

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 500] # TODO: decide the upper threshold, 0, 100, or other values?

reduced_CO2_emission = decarbonization_result.sum(axis=0)['CO2_reduction']/1000

added_oil = decarbonization_result.sum(axis=0)['oil_BPD']/1000000

national_CO2_reduction_ratio = reduced_CO2_emission/total_CO2_emission

national_CO2_recuction_ratio_sludge_management = reduced_CO2_emission/total_sludge_CO2_emission

national_oil_production_ratio = added_oil/total_oil

print(f'National decarbonization ratio of wastewater treatment sector is {national_CO2_reduction_ratio*100:.2f}%')

print(f'National decarbonization ratio of sludge management is {national_CO2_recuction_ratio_sludge_management*100:.2f}%')

print(f'National increase ratio of crude oil production is {national_oil_production_ratio*100:.5f}%')