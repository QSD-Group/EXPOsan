#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Sun Jun 11 08:12:41 2023

@author: jiananfeng
'''

# TODO 1: consider adding transportation cost and carbon intensity as contextual parameters
# TODO 1: see if https://www.bls.gov/oes/tables.htm (from Abby) has some information
# TODO 2: update electricity CI to match up balancing areas instead of states

#%% initialization
import geopy.distance, googlemaps, random
import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt, matplotlib.colors as colors, matplotlib.ticker as mtick
from matplotlib.mathtext import _mathtext as mathtext
from matplotlib.patches import Rectangle
from colorpalette import Color
from exposan.htl import create_geospatial_system, create_geospatial_model
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
p = Guest.purple.HEX

# TODO: consider adding these dark colors to qsdsan
db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

# !!! for figures, run 2 times to make the ticks and font size right
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

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)
# US.plot(ax=ax, color='w', edgecolor='k', linewidth=6) # for all WRRFs together with the same symbols

WRRF = WRRF.sort_values(by='flow_2022_MGD', ascending=False)

WRRF_flow = WRRF['flow_2022_MGD']

# both 'markersize' (in gpd.plot()) and 's' (in plt.scatter()) are proportional to area

less_than_10 = WRRF_flow <= 10
WRRF[less_than_10].plot(ax=ax, color=Guest.gray.HEX, markersize=WRRF.loc[less_than_10, WRRF_flow.name]*10, edgecolor='none', alpha=0.5)

between_10_and_100 = (WRRF_flow > 10) & (WRRF_flow <= 100)
WRRF[between_10_and_100].plot(ax=ax, color=Guest.red.HEX, markersize=WRRF.loc[between_10_and_100, WRRF_flow.name]*10, edgecolor='k', linewidth=2, alpha=0.7)

more_than_100 = WRRF_flow > 100
WRRF[more_than_100].plot(ax=ax, color=Guest.green.HEX, markersize=WRRF.loc[more_than_100, WRRF_flow.name]*10, edgecolor='k', linewidth=2, alpha=0.9)

# if you want to show all WRRFs together with the same symbols, comment the codes above and uncomment the next line of code
# WRRF.plot(ax=ax, color=Guest.gray.HEX, markersize=10)

#%% oil refinery visualization

set_plot()

US.plot(ax=ax,
        color=[r, b, g, b, b, r, g, b, o, b, g, y, r, g, r, y, r,
               b, b, g, o, o, g, o, b, g, y, g, b, o, g, b, b, y,
               b, b, b, b, b, b, g, g, g, y, g, r, g, g, b],
        alpha=0.3,
        edgecolor='k',
        linewidth=0)

# US.plot(ax=ax, color='none', edgecolor='k', linewidth=0) # for TOC

US.plot(ax=ax, color='none', edgecolor='k', linewidth=3)

# US.plot(ax=ax, color='none', edgecolor='k', linewidth=6) # for all oil refineries together with the same symbols

refinery['total_capacity'] = refinery[[i for i in refinery.columns if i[-4:] == 'Mbpd']].sum(axis=1)

refinery = refinery.sort_values(by='total_capacity', ascending=False)

refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=refinery['total_capacity']*4.5, edgecolor='k', linewidth=3, alpha=0.9)

# if you want to show all oil refineries together with the same symbols, comment the codes above and uncomment the next line of code
# refinery.plot(ax=ax, color=Guest.orange.HEX, markersize=1000, edgecolor='k', linewidth=6)

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

norm = colors.TwoSlopeNorm(vmin=4, vcenter=10, vmax=16)
electricity.plot('price (10-year median)', ax=ax, linewidth=3, cmap='Oranges', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35}, norm=norm)

# norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)
# electricity.plot('GHG (10-year median)', ax=ax, linewidth=3, cmap='Blues', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35}, norm=norm)

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

#%% travel distance box plot

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

fig, ax = plt.subplots(figsize = (5, 8))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()
ax.set_ylim([-50, 850])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

ax.set_ylabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=35)

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
    
fig.savefig('/Users/jiananfeng/Desktop/distance.png', transparent=True, bbox_inches='tight')

#%% travel distance box plot (per region)

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

WRRF_input.loc[WRRF_input['state'].isin(['CT','DC','DE','FL','GA','MA','MD','ME','NC','NH','NJ','NY','PA','RI','SC','VA','VT','WV']), 'WRRF_PADD'] = 1
WRRF_input.loc[WRRF_input['state'].isin(['IA','IL','IN','KS','KY','MI','MN','MO','ND','NE','OH','OK','SD','TN','WI']), 'WRRF_PADD'] = 2
WRRF_input.loc[WRRF_input['state'].isin(['AL','AR','LA','MS','NM','TX']), 'WRRF_PADD'] = 3
WRRF_input.loc[WRRF_input['state'].isin(['CO','ID','MT','UT','WY']), 'WRRF_PADD'] = 4
WRRF_input.loc[WRRF_input['state'].isin(['AZ','CA','NV','OR','WA']), 'WRRF_PADD'] = 5

fig = plt.figure(figsize=(20, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

def set_cf_plot():
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

def add_region(position, region, color):
    ax = fig.add_subplot(gs[0, position])
    
    set_cf_plot()
    
    ax = plt.gca()
    ax.set_ylim([-100, 1500])
    ax.set_xlabel(region, fontname='Arial', fontsize=30, labelpad=15)
    
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)
        ax.set_ylabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=35)
    
    elif position == 4:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
    
    else:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
    
    bp = ax.boxplot(WRRF_input[WRRF_input['WRRF_PADD'] == position+1]['real_distance_km'], showfliers=True, widths=0.5, patch_artist=True)
    
    for box in bp['boxes']:
        box.set(color='k', facecolor=color, linewidth=3)

    for whisker in bp['whiskers']:
        whisker.set(color='k', linewidth=3)

    for median in bp['medians']:
        median.set(color='k', linewidth=3)
        
    for cap in bp['caps']:
        cap.set(color='k', linewidth=3)
        
    for flier in bp['fliers']:
        flier.set(marker='o', markersize=7, markerfacecolor=color, markeredgewidth=1.5)

# TODO: there are still ticks on the left axis of the leftmost figure, use PS for now
add_region(0, 'East Coast', b)
add_region(1, 'Midwest', g)
add_region(2, 'Gulf Coast', r)
add_region(3, 'Rocky Mountain', o)
add_region(4, 'West Coast', y)

#%% WRRFs GHG map

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

WRRF_input = WRRF_input.sort_values(by='total_emission', ascending=False)

WRRF_input = gpd.GeoDataFrame(WRRF_input, crs='EPSG:4269',
                               geometry=gpd.points_from_xy(x=WRRF_input.longitude,
                                                           y=WRRF_input.latitude))

WRRF_input = WRRF_input.to_crs(crs='EPSG:3857')

set_plot()

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)

# use tonne/day

WRRF_GHG_tonne_per_day = WRRF_input['total_emission']/1000

less_than_50 = WRRF_GHG_tonne_per_day <= 50
WRRF_input[less_than_50].plot(ax=ax, color=Guest.gray.HEX, markersize=WRRF_input.loc[less_than_50, WRRF_GHG_tonne_per_day.name]**0.5, alpha=0.5)

between_50_and_500 = (WRRF_GHG_tonne_per_day > 50) & (WRRF_GHG_tonne_per_day <= 500)
WRRF_input[between_50_and_500].plot(ax=ax, color=Guest.red.HEX, markersize=WRRF_input.loc[between_50_and_500, WRRF_GHG_tonne_per_day.name]**0.5, edgecolor='k', linewidth=1.5, alpha=0.7)

more_than_500 = WRRF_GHG_tonne_per_day > 500
WRRF_input[more_than_500].plot(ax=ax, color=Guest.green.HEX, markersize=WRRF_input.loc[more_than_500, WRRF_GHG_tonne_per_day.name]**0.5, edgecolor='k', linewidth=2, alpha=0.9)

#%% WRRFs sludge management GHG map

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
WRRF_input[less_than_50].plot(ax=ax, color=Guest.gray.HEX, markersize=WRRF_input.loc[less_than_50, sludge_GHG_tonne_per_day.name]**0.5, alpha=0.5)

between_50_and_500 = (sludge_GHG_tonne_per_day > 50) & (sludge_GHG_tonne_per_day <= 500)
WRRF_input[between_50_and_500].plot(ax=ax, color=Guest.red.HEX, markersize=WRRF_input.loc[between_50_and_500, sludge_GHG_tonne_per_day.name]**0.5, edgecolor='k', linewidth=1.5, alpha=0.7)

more_than_500 = sludge_GHG_tonne_per_day > 500
WRRF_input[more_than_500].plot(ax=ax, color=Guest.green.HEX, markersize=WRRF_input.loc[more_than_500, sludge_GHG_tonne_per_day.name]**0.5, edgecolor='k', linewidth=2, alpha=0.9)

#%% cumulative WRRFs capacity vs distances (data processing)

# remember to use the correct file
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

result = WRRF_input[['site_id']].drop_duplicates()

max_distance = 1500 # this will cover all WRRFs
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

result.to_excel(folder + f'results/MGD_vs_distance_{max_distance}_km.xlsx')

#%% make the plot of cumulative WRRFs capacity vs distances (data preparation)

# import file for the cumulative figure (CF)
CF_input = pd.read_excel(folder + 'results/MGD_vs_distance_1500_km.xlsx')

CF_input[0] = 0

CF_input = CF_input[[0, *CF_input.columns[2:-1]]]

max_distance = 1500

CF_input = CF_input.iloc[:, 0:max_distance+6]

CF_input.drop(['Company','Corp','Site'], axis=1, inplace=True)

CF_input.sort_values(by='PADD', inplace=True)

CF_input = CF_input[['PADD','State', *CF_input.columns[0:-2]]]

CF_input = CF_input.transpose()

max_distance_on_the_plot = 1500

CF_input = CF_input[0:max_distance_on_the_plot+3]

#%% make the plot of cumulative WRRFs capacity vs distances (separated oil refinery)

PADD_1 = (CF_input.loc['PADD',:].isin([1,])).sum()
PADD_2 = (CF_input.loc['PADD',:].isin([1,2])).sum()
PADD_3 = (CF_input.loc['PADD',:].isin([1,2,3])).sum()
PADD_4 = (CF_input.loc['PADD',:].isin([1,2,3,4])).sum()
PADD_5 = (CF_input.loc['PADD',:].isin([1,2,3,4,5])).sum()

fig = plt.figure(figsize=(20, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

def set_cf_plot():
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

def add_region(position, start_region, end_region, color):
    ax = fig.add_subplot(gs[0, position])
    
    set_cf_plot()
    
    for i in range(start_region, end_region):
        CF_input.iloc[2:, i].plot(ax=ax, color=color, linewidth=3)
        
    ax.set_xlim([0, max_distance])
    ax.set_ylim([0, 7000])
    
    if position == 4:
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        plt.yticks(np.arange(0, 7000, 1000))
        ax_right.tick_params(direction='in', length=7.5, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
    if position == 0:
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=True, right=False, labelleft=True, labelbottom=True)
        plt.xticks(np.arange(0, max_distance*1.2, max_distance*0.2))
        ax.set_ylabel(r'$\mathbf{Cumulative\ WRRFs\ capacity}$ [MGD]', fontname='Arial', fontsize=35)
    else:
        if position == 2:
            ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=35)
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        plt.xticks(np.arange(max_distance*0.2, max_distance*1.2, max_distance*0.2))
    
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        # label.set_ha('right')
    
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    plt.xticks(np.arange(max_distance*0.2, max_distance*1.2, max_distance*0.2))
    ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

add_region(0, 0, PADD_1, b)
add_region(1, PADD_1, PADD_2, g)
add_region(2, PADD_2, PADD_3, r)
add_region(3, PADD_3, PADD_4, o)
add_region(4, PADD_4, PADD_5, y)
#%% make the plot of cumulative WRRFs capacity vs distances (regional total)

PADD_1 = (CF_input.loc['PADD',:].isin([1,])).sum()
PADD_2 = (CF_input.loc['PADD',:].isin([1,2])).sum()
PADD_3 = (CF_input.loc['PADD',:].isin([1,2,3])).sum()
PADD_4 = (CF_input.loc['PADD',:].isin([1,2,3,4])).sum()
PADD_5 = (CF_input.loc['PADD',:].isin([1,2,3,4,5])).sum()

fig, ax = plt.subplots(figsize=(12, 10))

def set_cf_plot():
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

set_cf_plot()

CF_input.iloc[2:, 0:PADD_1].sum(axis=1).plot(ax=ax, color=b, linewidth=3)
CF_input.iloc[2:, PADD_1:PADD_2].sum(axis=1).plot(ax=ax, color=g, linewidth=3)
CF_input.iloc[2:, PADD_2:PADD_3].sum(axis=1).plot(ax=ax, color=r, linewidth=3)
CF_input.iloc[2:, PADD_3:PADD_4].sum(axis=1).plot(ax=ax, color=o, linewidth=3)
CF_input.iloc[2:, PADD_4:PADD_5].sum(axis=1).plot(ax=ax, color=y, linewidth=3)

ax.set_xlim([0, max_distance])
ax.set_ylim([0, 18000])

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=35)
ax.set_ylabel(r'$\mathbf{Cumulative\ WRRFs\ capacity}$ [MGD]', fontname='Arial', fontsize=35)

ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
plt.yticks(np.arange(0, 20000, 2000))
ax_right.tick_params(direction='in', length=7.5, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 1750, 250))
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

#%% CO2 abatement cost analysis

filterwarnings('ignore')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')
elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

# if just want to see the two plants in Urbana-Champaign:
# WRRF_input = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

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

for i in range(0, len(WRRF_input)): # !!! run in different consoles to speed up: 0, 5000, 10000
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
    raw_wastewater = stream.sludge_assumed_in_wastewater
    
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

output_result_1 = pd.read_excel(folder + 'results/decarbonization_biocrude_baseline_data/decarbonization_2023-12-31_4999.xlsx')
output_result_2 = pd.read_excel(folder + 'results/decarbonization_biocrude_baseline_data/decarbonization_2023-12-31_9999.xlsx')
output_result_3 = pd.read_excel(folder + 'results/decarbonization_biocrude_baseline_data/decarbonization_2023-12-31_14952.xlsx')

output_result = pd.concat([output_result_1, output_result_2, output_result_3])

assert input_data[['CWNS','facility_code']].duplicated().sum() == 0

assert output_result[['CWNS','facility_code']].duplicated().sum() == 0

integrated_result = input_data.merge(output_result, how='left', on=['CWNS','facility_code'])

integrated_result.to_excel(folder + f'results/integrated_decarbonization_result_{date.today()}.xlsx')

#%% read decarbonization data

decarbonization_result = pd.read_excel(folder + 'results/integrated_decarbonization_result_2023-12-31.xlsx')

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 0]

#%% decarbonization map

decarbonization_map = decarbonization_result.sort_values(by='total_emission', ascending=False)

decarbonization_map = gpd.GeoDataFrame(decarbonization_map, crs='EPSG:4269',
                                       geometry=gpd.points_from_xy(x=decarbonization_map.longitude,
                                                                   y=decarbonization_map.latitude))

decarbonization_map = decarbonization_map.to_crs(crs='EPSG:3857')

set_plot(figure_size=(15, 10))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=2)

# use tonne/day

decarbonization_map['CO2_reduction_tonne_per_day'] = decarbonization_map['CO2_reduction']/30/365/1000

WRRF_GHG_reduction_tonne_per_day = decarbonization_map['CO2_reduction_tonne_per_day']

less_than_50 = WRRF_GHG_reduction_tonne_per_day <= 5
decarbonization_map[less_than_50].plot(ax=ax, color=Guest.gray.HEX, markersize=decarbonization_map.loc[less_than_50, WRRF_GHG_reduction_tonne_per_day.name]*10, alpha=0.5)

between_50_and_500 = (WRRF_GHG_reduction_tonne_per_day > 5) & (WRRF_GHG_reduction_tonne_per_day <= 50)
decarbonization_map[between_50_and_500].plot(ax=ax, color=Guest.red.HEX, markersize=decarbonization_map.loc[between_50_and_500, WRRF_GHG_reduction_tonne_per_day.name]*10, edgecolor='k', linewidth=1.5, alpha=0.7)

more_than_500 = WRRF_GHG_reduction_tonne_per_day > 50
decarbonization_map[more_than_500].plot(ax=ax, color=Guest.green.HEX, markersize=decarbonization_map.loc[more_than_500, WRRF_GHG_reduction_tonne_per_day.name]*10, edgecolor='k', linewidth=1.5, alpha=0.9)

#%% facility level decarbonizaiton ratio and biocrude production

# import decarbonization_result from #%% read decarbonization data

fig = plt.figure(figsize=(15, 10))

gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)

def facility_plot(position, color):
    ax = fig.add_subplot(gs[0, position])
    
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    ax = plt.gca()
    ax.set_xlim([-6, 36])
    ax.set_ylim([0, 700])
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    
    if position == 0:
        data = decarbonization_result[(decarbonization_result['sludge_anaerobic_digestion'] != 1) & (decarbonization_result['sludge_aerobic_digestion'] != 1)]
        
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=True, right=False)
        
        plt.xticks(np.arange(0, 40, 10))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 40, 10))
        ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
        
        ax.set_ylabel(r'$\mathbf{Biocrude\ production}$ [BPD]', fontname='Arial', fontsize=35)
        
    if position == 1:
        data = decarbonization_result[decarbonization_result['sludge_aerobic_digestion'] == 1]

        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 40, 10))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 40, 10))
        ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')
        
        ax.set_xlabel(r'$\mathbf{Decarbonization\ potential}$', fontname='Arial', fontsize=35)
        
    if position == 2:
        data = decarbonization_result[decarbonization_result['sludge_anaerobic_digestion'] == 1]
        
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 40, 10))
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=7.5, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 40, 10))
        ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
    
    ax.scatter(x=data['WRRF_CO2_reduction_ratio']*100,
               y=data['oil_BPD'],
               s=data['flow_2022_MGD']*2,
               c=color,
               linewidths=0,
               alpha=0.9)
    
    ax.scatter(x=data['WRRF_CO2_reduction_ratio']*100,
               y=data['oil_BPD'],
               s=data['flow_2022_MGD']*2,
               c='none',
               linewidths=2,
               edgecolors='k')
    
facility_plot(0, p)
facility_plot(1, y)
facility_plot(2, b)

#%% facility level decarbonizaiton amount and biocrude production

# import decarbonization_result from #%% read decarbonization data

fig = plt.figure(figsize=(16, 8))

gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)

def facility_plot(position, color):
    
    ax = fig.add_subplot(gs[0, position])
    
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    ax = plt.gca()
    ax.set_xlim([-15, 165])
    ax.set_ylim([0, 600])
    
    if position == 0:
        data = decarbonization_result[(decarbonization_result['sludge_anaerobic_digestion'] != 1) & (decarbonization_result['sludge_aerobic_digestion'] != 1)]
        
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=True, right=False)
        
        plt.xticks(np.arange(0, 180, 30))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 180, 30))
        ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
        
        ax.set_ylabel(r'$\mathbf{Biocrude\ production}$ [BPD]', fontname='Arial', fontsize=35)
        
    if position == 1:
        data = decarbonization_result[decarbonization_result['sludge_aerobic_digestion'] == 1]

        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 180, 30))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 180, 30))
        ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')
        
        ax.set_xlabel(r'$\mathbf{Decarbonization\ amount}$ [tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=35)

        mathtext.FontConstantsBase.sup1 = 0.35
        
    if position == 2:
        data = decarbonization_result[decarbonization_result['sludge_anaerobic_digestion'] == 1]
    
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 180, 30))
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=7.5, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
        ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 180, 30))
        ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
    
    ax.scatter(x=data['CO2_reduction']/30/365/1000,
               y=data['oil_BPD'],
               s=data['flow_2022_MGD']*2,
               c=color,
               linewidths=0,
               alpha=0.9)
    
    ax.scatter(x=data['CO2_reduction']/30/365/1000,
               y=data['oil_BPD'],
               s=data['flow_2022_MGD']*2,
               c='none',
               linewidths=2,
               edgecolors='k')
    
facility_plot(0, p)
facility_plot(1, y)
facility_plot(2, b)
#%% find WRRFs with max decarbonization ratio, decarbonization amount, and biocrude production

no_digestion_WRRFs = decarbonization_result[(decarbonization_result['sludge_anaerobic_digestion'] != 1) & (decarbonization_result['sludge_aerobic_digestion'] != 1)]
print(no_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(no_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(no_digestion_WRRFs.sort_values('oil_BPD', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])

aerobic_digestion_WRRFs = decarbonization_result[decarbonization_result['sludge_aerobic_digestion'] == 1]
print(aerobic_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(aerobic_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(aerobic_digestion_WRRFs.sort_values('oil_BPD', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])

anaerobic_digestion_WRRFs = decarbonization_result[decarbonization_result['sludge_anaerobic_digestion'] == 1]
print(anaerobic_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(anaerobic_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(anaerobic_digestion_WRRFs.sort_values('oil_BPD', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])

#%% facility level decarbonizaiton amount vs sludge amount

decarbonization_vs_sludge = decarbonization_result[['total_sludge_amount_kg_per_year','CO2_reduction']]

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.set_xlim((0, 350))
ax.set_ylim((0, 125))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Sludge}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Decarbonization}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 400, 50))
plt.yticks(np.arange(0, 150, 25))

ax_right = ax.twinx()
ax_right.set_ylim((0, 125))
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 400, 50))
plt.yticks(np.arange(0, 150, 25))

ax_top = ax.twiny()
ax_top.set_xlim((0, 350))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 400, 50))
plt.yticks(np.arange(0, 150, 25))

ax.scatter(x=decarbonization_vs_sludge['total_sludge_amount_kg_per_year']/1000/365,
           y=decarbonization_vs_sludge['CO2_reduction']/30/365/1000,
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level decarbonizaiton amount vs distance

decarbonization_vs_sludge = decarbonization_result[['real_distance_km','CO2_reduction']]

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.set_xlim((0, 1200))
ax.set_ylim((0, 125))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Decarbonization}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 150, 25))

ax_right = ax.twinx()
ax_right.set_ylim((0, 125))
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 150, 25))

ax_top = ax.twiny()
ax_top.set_xlim((0, 1200))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 150, 25))

ax.scatter(x=decarbonization_vs_sludge['real_distance_km'],
           y=decarbonization_vs_sludge['CO2_reduction']/30/365/1000,
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level biocrude production vs sludge amount

decarbonization_vs_sludge = decarbonization_result[['total_sludge_amount_kg_per_year','oil_BPD']]

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.set_xlim((0, 350))
ax.set_ylim((0, 600))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Sludge}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Biocrude}$ [BPD]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 400, 50))
plt.yticks(np.arange(0, 700, 100))

ax_right = ax.twinx()
ax_right.set_ylim((0, 125))
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 400, 50))
plt.yticks(np.arange(0, 700, 100))

ax_top = ax.twiny()
ax_top.set_xlim((0, 350))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 400, 50))
plt.yticks(np.arange(0, 700, 100))

ax.scatter(x=decarbonization_vs_sludge['total_sludge_amount_kg_per_year']/1000/365,
           y=decarbonization_vs_sludge['oil_BPD'],
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level biocrude production vs distance

decarbonization_vs_sludge = decarbonization_result[['real_distance_km','oil_BPD']]

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.set_xlim((0, 1200))
ax.set_ylim((0, 600))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Biocrude}$ [BPD]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 700, 100))

ax_right = ax.twinx()
ax_right.set_ylim((0, 125))
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 700, 100))

ax_top = ax.twiny()
ax_top.set_xlim((0, 350))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 700, 100))

ax.scatter(x=decarbonization_vs_sludge['real_distance_km'],
           y=decarbonization_vs_sludge['oil_BPD'],
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level uncertainty analysis

filterwarnings('ignore')

decarbonization_result = pd.read_excel(folder + 'results/integrated_decarbonization_result_2023-12-31.xlsx')
elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]

decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 0]

print(len(decarbonization_result))

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

decarbonization_result['waste_cost'] = sum(decarbonization_result[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/decarbonization_result['total_sludge_amount_kg_per_year']*1000

decarbonization_result['waste_GHG'] =  sum(decarbonization_result[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/decarbonization_result['total_sludge_amount_kg_per_year']*1000

geo_uncertainty_decarbonization = pd.DataFrame()
geo_uncertainty_biocrude = pd.DataFrame()

for i in range(0, len(decarbonization_result)): # !!! run in different consoles to speed up: 0, 70, 140, 210, 280, 350, 420, 490
# for i in range(0, 2):
    
    sys, barrel = create_geospatial_system(waste_cost=decarbonization_result.iloc[i]['waste_cost'],
                                           waste_GHG=decarbonization_result.iloc[i]['waste_GHG'],
                                           size=decarbonization_result.iloc[i]['flow_2022_MGD'],
                                           distance=decarbonization_result.iloc[i]['real_distance_km'],
                                           anaerobic_digestion=decarbonization_result.iloc[i]['sludge_anaerobic_digestion'],
                                           aerobic_digestion=decarbonization_result.iloc[i]['sludge_aerobic_digestion'],
                                           ww_2_dry_sludge_ratio=decarbonization_result.iloc[i]['total_sludge_amount_kg_per_year']/1000/365/decarbonization_result.iloc[i]['flow_2022_MGD'],
                                           # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                           state=decarbonization_result.iloc[i]['state'],
                                           elec_GHG=float(elec[elec['state']==decarbonization_result.iloc[i]['state']]['GHG (10-year median)']))
    
    if decarbonization_result.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_values=[0.374, 0.468, 0.562, 'aerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612]
        
    elif decarbonization_result.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612]
        
    else:
        sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
        sludge_lipid_values=[0.080, 0.206, 0.308]
        sludge_protein_values=[0.380, 0.456, 0.485]
    
    model = create_geospatial_model(system=sys,
                                    sludge_ash=sludge_ash_values,
                                    sludge_lipid=sludge_lipid_values,
                                    sludge_protein=sludge_protein_values,
                                    raw_wastewater_price_baseline=sys.flowsheet.stream.sludge_assumed_in_wastewater.price,
                                    biocrude_and_transportation_price=[sys.flowsheet.stream.biocrude.price/6.80*4.21, # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                                                       sys.flowsheet.stream.biocrude.price,
                                                                       sys.flowsheet.stream.biocrude.price/6.80*11.9], # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                    electricity_cost=[elec[elec['state']==decarbonization_result.iloc[i]['state']]['price (10-year 5th)'].iloc[0]/100,
                                                      elec[elec['state']==decarbonization_result.iloc[i]['state']]['price (10-year median)'].iloc[0]/100,
                                                      elec[elec['state']==decarbonization_result.iloc[i]['state']]['price (10-year 95th)'].iloc[0]/100],
                                    electricity_GHG=[0.67848*elec[elec['state']==decarbonization_result.iloc[i]['state']]['GHG (10-year 5th)']/elec[elec['state']==decarbonization_result.iloc[i]['state']]['GHG (10-year median)'],
                                                     0.67848,
                                                     0.67848*elec[elec['state']==decarbonization_result.iloc[i]['state']]['GHG (10-year 95th)']/elec[elec['state']==decarbonization_result.iloc[i]['state']]['GHG (10-year median)']])
    
    kwargs = {'N':1000, 'rule':'L', 'seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    geo_uncertainty_decarbonization[decarbonization_result.iloc[i]['facility']] = results[('Geospatial','Decarbonization amount [tonne_per_day]')]
    geo_uncertainty_biocrude[decarbonization_result.iloc[i]['facility']] = results[('Geospatial','Biocrude production [BPD]')]
    
    if i%5 == 0:
        print(i) # check progress
    
geo_uncertainty_decarbonization.to_excel(folder + f'results/regional_decarbonization_uncertainty_{date.today()}_{i}.xlsx')
geo_uncertainty_biocrude.to_excel(folder + f'results/regional_biocrude_uncertainty_{date.today()}_{i}.xlsx')

#%% regional uncertainty (data preparation)

# decarbonization
regional_decarbonization_1 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_69.xlsx')
regional_decarbonization_2 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_139.xlsx')
regional_decarbonization_3 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_209.xlsx')
regional_decarbonization_4 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_279.xlsx')
regional_decarbonization_5 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_349.xlsx')
regional_decarbonization_6 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_419.xlsx')
regional_decarbonization_7 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_489.xlsx')
regional_decarbonization_8 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2023-12-31_559.xlsx')

regional_decarbonization_1.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_2.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_3.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_4.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_5.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_6.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_7.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_8.drop('Unnamed: 0', axis=1, inplace=True)

regional_decarbonization = pd.concat([regional_decarbonization_1, regional_decarbonization_2], axis=1)
regional_decarbonization = pd.concat([regional_decarbonization, regional_decarbonization_3], axis=1)
regional_decarbonization = pd.concat([regional_decarbonization, regional_decarbonization_4], axis=1)
regional_decarbonization = pd.concat([regional_decarbonization, regional_decarbonization_5], axis=1)
regional_decarbonization = pd.concat([regional_decarbonization, regional_decarbonization_6], axis=1)
regional_decarbonization = pd.concat([regional_decarbonization, regional_decarbonization_7], axis=1)
regional_decarbonization = pd.concat([regional_decarbonization, regional_decarbonization_8], axis=1)

regional_decarbonization.to_excel(folder + f'results/integrated_regional_decarbonization_uncertainty_{date.today()}.xlsx')

# biocrude
regional_biocrude_1 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_69.xlsx')
regional_biocrude_2 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_139.xlsx')
regional_biocrude_3 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_209.xlsx')
regional_biocrude_4 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_279.xlsx')
regional_biocrude_5 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_349.xlsx')
regional_biocrude_6 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_419.xlsx')
regional_biocrude_7 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_489.xlsx')
regional_biocrude_8 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2023-12-31_559.xlsx')

regional_biocrude_1.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_2.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_3.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_4.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_5.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_6.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_7.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_8.drop('Unnamed: 0', axis=1, inplace=True)

regional_biocrude = pd.concat([regional_biocrude_1, regional_biocrude_2], axis=1)
regional_biocrude = pd.concat([regional_biocrude, regional_biocrude_3], axis=1)
regional_biocrude = pd.concat([regional_biocrude, regional_biocrude_4], axis=1)
regional_biocrude = pd.concat([regional_biocrude, regional_biocrude_5], axis=1)
regional_biocrude = pd.concat([regional_biocrude, regional_biocrude_6], axis=1)
regional_biocrude = pd.concat([regional_biocrude, regional_biocrude_7], axis=1)
regional_biocrude = pd.concat([regional_biocrude, regional_biocrude_8], axis=1)

regional_biocrude.to_excel(folder + f'results/integrated_regional_biocrude_uncertainty_{date.today()}.xlsx')

#%% regional uncertainty (build and run model)

# import PADD information
WRRF_PADD = pd.read_excel(folder + 'results/integrated_decarbonization_result_2023-12-31.xlsx')

WRRF_PADD = WRRF_PADD[WRRF_PADD['USD_decarbonization'].notna()]

WRRF_PADD = WRRF_PADD[WRRF_PADD['USD_decarbonization'] <= 0]

WRRF_PADD.loc[WRRF_PADD['state'].isin(['CT','DC','DE','FL','GA','MA','MD','ME','NC','NH','NJ','NY','PA','RI','SC','VA','VT','WV']), 'WRRF_PADD'] = 1
WRRF_PADD.loc[WRRF_PADD['state'].isin(['IA','IL','IN','KS','KY','MI','MN','MO','ND','NE','OH','OK','SD','TN','WI']), 'WRRF_PADD'] = 2
WRRF_PADD.loc[WRRF_PADD['state'].isin(['AL','AR','LA','MS','NM','TX']), 'WRRF_PADD'] = 3
WRRF_PADD.loc[WRRF_PADD['state'].isin(['CO','ID','MT','UT','WY']), 'WRRF_PADD'] = 4
WRRF_PADD.loc[WRRF_PADD['state'].isin(['AZ','CA','NV','OR','WA']), 'WRRF_PADD'] = 5

WRRF_PADD = WRRF_PADD[['facility','WRRF_PADD']]

# decarbonization
regional_decarbonization_uncertainty = pd.read_excel(folder + 'results/integrated_regional_decarbonization_uncertainty_2024-01-01.xlsx')

regional_decarbonization_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)

regional_decarbonization_uncertainty = regional_decarbonization_uncertainty.transpose()

regional_decarbonization_uncertainty.reset_index(inplace=True, names='facility')

regional_decarbonization_uncertainty = regional_decarbonization_uncertainty.merge(WRRF_PADD, how='left', on='facility')

regional_decarbonization_uncertainty_PADD_1 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 1].iloc[:,1:-1]
regional_decarbonization_uncertainty_PADD_2 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 2].iloc[:,1:-1]
regional_decarbonization_uncertainty_PADD_3 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 3].iloc[:,1:-1]
regional_decarbonization_uncertainty_PADD_4 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 4].iloc[:,1:-1]
regional_decarbonization_uncertainty_PADD_5 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 5].iloc[:,1:-1]

PADD_1_decarbonization_uncertainty = []
PADD_2_decarbonization_uncertainty = []
PADD_3_decarbonization_uncertainty = []
PADD_4_decarbonization_uncertainty = []
PADD_5_decarbonization_uncertainty = []

for (region, result) in [(regional_decarbonization_uncertainty_PADD_1, PADD_1_decarbonization_uncertainty),
                         (regional_decarbonization_uncertainty_PADD_2, PADD_2_decarbonization_uncertainty),
                         (regional_decarbonization_uncertainty_PADD_3, PADD_3_decarbonization_uncertainty),
                         (regional_decarbonization_uncertainty_PADD_4, PADD_4_decarbonization_uncertainty),
                         (regional_decarbonization_uncertainty_PADD_5, PADD_5_decarbonization_uncertainty)]:
    for i in range(0,1000): # Monte Carlo simulation and Latin hypercube sampling
        if i%100 == 0:
            print(i) # check progress
        regional_total_decarbonization = 0
        for j in range(0, len(region)):
            regional_total_decarbonization += random.uniform(region.iloc[j,:].quantile(i/1000), region.iloc[j,:].quantile((i+1)/1000))
        result.append(regional_total_decarbonization)
            
integrated_regional_decarbonization_result = pd.DataFrame({'PADD_1':PADD_1_decarbonization_uncertainty,
                                                           'PADD_2':PADD_2_decarbonization_uncertainty,
                                                           'PADD_3':PADD_3_decarbonization_uncertainty,
                                                           'PADD_4':PADD_4_decarbonization_uncertainty,
                                                           'PADD_5':PADD_5_decarbonization_uncertainty})
        
integrated_regional_decarbonization_result.to_excel(folder + f'results/integrated_regional_total_decarbonization_uncertainty_{date.today()}.xlsx')

# biocrude
regional_biocrude_uncertainty = pd.read_excel(folder + 'results/integrated_regional_biocrude_uncertainty_2024-01-01.xlsx')

regional_biocrude_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)

regional_biocrude_uncertainty = regional_biocrude_uncertainty.transpose()

regional_biocrude_uncertainty.reset_index(inplace=True, names='facility')

regional_biocrude_uncertainty = regional_biocrude_uncertainty.merge(WRRF_PADD, how='left', on='facility')

regional_biocrude_uncertainty_PADD_1 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 1].iloc[:,1:-1]
regional_biocrude_uncertainty_PADD_2 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 2].iloc[:,1:-1]
regional_biocrude_uncertainty_PADD_3 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 3].iloc[:,1:-1]
regional_biocrude_uncertainty_PADD_4 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 4].iloc[:,1:-1]
regional_biocrude_uncertainty_PADD_5 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 5].iloc[:,1:-1]

PADD_1_biocrude_uncertainty = []
PADD_2_biocrude_uncertainty = []
PADD_3_biocrude_uncertainty = []
PADD_4_biocrude_uncertainty = []
PADD_5_biocrude_uncertainty = []

for (region, result) in [(regional_biocrude_uncertainty_PADD_1, PADD_1_biocrude_uncertainty),
                         (regional_biocrude_uncertainty_PADD_2, PADD_2_biocrude_uncertainty),
                         (regional_biocrude_uncertainty_PADD_3, PADD_3_biocrude_uncertainty),
                         (regional_biocrude_uncertainty_PADD_4, PADD_4_biocrude_uncertainty),
                         (regional_biocrude_uncertainty_PADD_5, PADD_5_biocrude_uncertainty)]:
    for i in range(0,1000): # Monte Carlo simulation and Latin hypercube sampling
        if i%100 == 0:
            print(i) # check progress
        regional_total_biocrude = 0
        for j in range(0, len(region)):
            regional_total_biocrude += random.uniform(region.iloc[j,:].quantile(i/1000), region.iloc[j,:].quantile((i+1)/1000))
        result.append(regional_total_biocrude)
            
integrated_regional_biocrude_result = pd.DataFrame({'PADD_1':PADD_1_biocrude_uncertainty,
                                                    'PADD_2':PADD_2_biocrude_uncertainty,
                                                    'PADD_3':PADD_3_biocrude_uncertainty,
                                                    'PADD_4':PADD_4_biocrude_uncertainty,
                                                    'PADD_5':PADD_5_biocrude_uncertainty})
        
integrated_regional_biocrude_result.to_excel(folder + f'results/integrated_regional_total_biocrude_uncertainty_{date.today()}.xlsx')

#%% regional uncertainty (visualization)

decarbonization = pd.read_excel(folder + 'results/integrated_regional_total_decarbonization_uncertainty_2024-01-01.xlsx')
biocrude = pd.read_excel(folder + 'results/integrated_regional_total_biocrude_uncertainty_2024-01-01.xlsx')

fig, ax = plt.subplots(figsize = (12, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax.set_xlim([0, 2400])
ax.set_ylim([0, 12000])

ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 2800, 400))
plt.yticks(np.arange(0, 14000, 2000))

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 2800, 400))
plt.yticks(np.arange(0, 14000, 2000))

ax.set_xlabel(r'$\mathbf{Decarbonization\ amount}$ [tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=35)
ax.set_ylabel(r'$\mathbf{Biocrude\ production}$ [BPD]', fontname='Arial', fontsize=35)

mathtext.FontConstantsBase.sup1 = 0.35

def add_rectangle(region, color, edgecolor):
    rectangle_fill = Rectangle((decarbonization[region].quantile(0.05), biocrude[region].quantile(0.05)),
                               decarbonization[region].quantile(0.95) - decarbonization[region].quantile(0.05),
                               biocrude[region].quantile(0.95) - biocrude[region].quantile(0.05),
                               fc=color, alpha=0.8)
    ax.add_patch(rectangle_fill)
    
    rectangle_edge = Rectangle((decarbonization[region].quantile(0.05), biocrude[region].quantile(0.05)),
                               decarbonization[region].quantile(0.95) - decarbonization[region].quantile(0.05),
                               biocrude[region].quantile(0.95) - biocrude[region].quantile(0.05),
                               color=edgecolor, lw=3, fc='none', alpha=1)
    ax.add_patch(rectangle_edge)
    
def add_line(region, color):
    plt.plot([decarbonization[region].quantile(0.25), decarbonization[region].quantile(0.75)],
             [biocrude[region].quantile(0.5), biocrude[region].quantile(0.5)],
             lw=3, color=color, solid_capstyle='round', zorder=1)
    plt.plot([decarbonization[region].quantile(0.5), decarbonization[region].quantile(0.5)],
             [biocrude[region].quantile(0.25), biocrude[region].quantile(0.75)],
             lw=3, color=color, solid_capstyle='round', zorder=1)

def add_point(region, edgecolor):
    ax_top.scatter(x=decarbonization[region].quantile(0.5),
                   y=biocrude[region].quantile(0.5),
                   marker='s',
                   s=150,
                   c='w',
                   linewidths=3,
                   alpha=1,
                   edgecolor=edgecolor)
    
add_rectangle('PADD_1', b, db)
add_line('PADD_1', db)

add_rectangle('PADD_2', g, dg)
add_line('PADD_2', dg)

add_rectangle('PADD_3', r, dr)
add_line('PADD_3', dr)

add_rectangle('PADD_4', o, do)
add_line('PADD_4', do)

add_rectangle('PADD_5', y, dy)
add_line('PADD_5', dy)

add_point('PADD_1', db)
add_point('PADD_2', dg)
add_point('PADD_3', dr)
# add_point('PADD_4', do)
add_point('PADD_5', dy)

#%% national uncertainty

# decarbonization
national_decarbonization_uncertainty = pd.read_excel(folder + 'results/integrated_regional_decarbonization_uncertainty_2024-01-01.xlsx')

national_decarbonization_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)

national_decarbonization_uncertainty = national_decarbonization_uncertainty.transpose()

national_decarbonization_uncertainty.reset_index(inplace=True, names='facility')

national_decarbonization_uncertainty = national_decarbonization_uncertainty.iloc[:,1:]

national_decarbonization_uncertainty_result = []

for i in range(0,1000): # Monte Carlo simulation and Latin hypercube sampling
    if i%100 == 0:
        print(i) # check progress
    national_total_decarbonization = 0
    for j in range(0, len(national_decarbonization_uncertainty)):
        national_total_decarbonization += random.uniform(national_decarbonization_uncertainty.iloc[j,:].quantile(i/1000), national_decarbonization_uncertainty.iloc[j,:].quantile((i+1)/1000))
    national_decarbonization_uncertainty_result.append(national_total_decarbonization)
    
# biocrude
national_biocrude_uncertainty = pd.read_excel(folder + 'results/integrated_regional_biocrude_uncertainty_2024-01-01.xlsx')

national_biocrude_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)

national_biocrude_uncertainty = national_biocrude_uncertainty.transpose()

national_biocrude_uncertainty.reset_index(inplace=True, names='facility')

national_biocrude_uncertainty = national_biocrude_uncertainty.iloc[:,1:]

national_biocrude_uncertainty_result = []

for i in range(0,1000): # Monte Carlo simulation and Latin hypercube sampling
    if i%100 == 0:
        print(i) # check progress
    national_total_biocrude = 0
    for j in range(0, len(national_biocrude_uncertainty)):
        national_total_biocrude += random.uniform(national_biocrude_uncertainty.iloc[j,:].quantile(i/1000), national_biocrude_uncertainty.iloc[j,:].quantile((i+1)/1000))
    national_biocrude_uncertainty_result.append(national_total_biocrude)
            
national_uncertainty_result = pd.DataFrame({'decarbonization':national_decarbonization_uncertainty_result,
                                            'biocrude':national_biocrude_uncertainty_result})

national_uncertainty_result.to_excel(folder + f'results/integrated_national_uncertainty_{date.today()}.xlsx')

#%% sludge transportation (Urbana-Champaign)

# two separated WRRFs
filterwarnings('ignore')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')
elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

# if just want to see the two plants in Urbana-Champaign:
# WRRF_input = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

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

Urbana_Champaign = WRRF_input[WRRF_input['CWNS'].isin([17000112001,17000112002])]

CU_uncertainty_NPV = pd.DataFrame()
CU_uncertainty_biocrude = pd.DataFrame()
CU_uncertainty_decarbonization = pd.DataFrame()
CU_uncertainty_decarbonization_cost = pd.DataFrame()

for i in range(0, 2):
    
    sys, barrel = create_geospatial_system(waste_cost=Urbana_Champaign.iloc[i]['waste_cost'],
                                           waste_GHG=Urbana_Champaign.iloc[i]['waste_GHG'],
                                           size=Urbana_Champaign.iloc[i]['flow_2022_MGD'],
                                           distance=Urbana_Champaign.iloc[i]['real_distance_km'],
                                           anaerobic_digestion=Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'],
                                           aerobic_digestion=Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'],
                                           ww_2_dry_sludge_ratio=Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365/Urbana_Champaign.iloc[i]['flow_2022_MGD'],
                                           # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                           state=Urbana_Champaign.iloc[i]['state'],
                                           elec_GHG=float(elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['GHG (10-year median)']))
    
    if Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_values=[0.374, 0.468, 0.562, 'aerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612]
        
    elif Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612]
        
    else:
        sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
        sludge_lipid_values=[0.080, 0.206, 0.308]
        sludge_protein_values=[0.380, 0.456, 0.485]
    
    model = create_geospatial_model(system=sys,
                                    sludge_ash=sludge_ash_values,
                                    sludge_lipid=sludge_lipid_values,
                                    sludge_protein=sludge_protein_values,
                                    raw_wastewater_price_baseline=sys.flowsheet.stream.sludge_assumed_in_wastewater.price,
                                    biocrude_and_transportation_price=[sys.flowsheet.stream.biocrude.price/6.80*4.21, # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                                                       sys.flowsheet.stream.biocrude.price,
                                                                       sys.flowsheet.stream.biocrude.price/6.80*11.9], # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                    electricity_cost=[elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['price (10-year 5th)'].iloc[0]/100,
                                                      elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['price (10-year median)'].iloc[0]/100,
                                                      elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['price (10-year 95th)'].iloc[0]/100],
                                    electricity_GHG=[0.67848*elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['GHG (10-year 5th)']/elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['GHG (10-year median)'],
                                                     0.67848,
                                                     0.67848*elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['GHG (10-year 95th)']/elec[elec['state']==Urbana_Champaign.iloc[i]['state']]['GHG (10-year median)']])
    
    kwargs = {'N':1000, 'rule':'L', 'seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    CU_uncertainty_NPV[Urbana_Champaign.iloc[i]['facility']] = results[('Geospatial','NPV [$]')]
    CU_uncertainty_biocrude[Urbana_Champaign.iloc[i]['facility']] = results[('Geospatial','Biocrude production [BPD]')]
    CU_uncertainty_decarbonization[Urbana_Champaign.iloc[i]['facility']] = results[('Geospatial','Decarbonization amount [tonne_per_day]')]
    CU_uncertainty_decarbonization_cost[Urbana_Champaign.iloc[i]['facility']] = -results[('Geospatial','NPV [$]')]/results[('Geospatial','Decarbonization amount [tonne_per_day]')]/365/30

# after sludge transportation
CU_combined_uncertainty_NPV = pd.DataFrame()
CU_combined_uncertainty_biocrude = pd.DataFrame()
CU_combined_uncertainty_decarbonization = pd.DataFrame()
CU_combined_uncertainty_decarbonization_cost = pd.DataFrame()

average_cost = (Urbana_Champaign.iloc[0]['waste_cost']*Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] +\
                Urbana_Champaign.iloc[1]['waste_cost']*Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']) /\
               (Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] + Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year'])

average_GHG = (Urbana_Champaign.iloc[0]['waste_GHG']*Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] +\
               Urbana_Champaign.iloc[1]['waste_GHG']*Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']) /\
              (Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] + Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year'])
               
total_size = Urbana_Champaign.iloc[0]['flow_2022_MGD']+Urbana_Champaign.iloc[1]['flow_2022_MGD']

distance_to_refinery = Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001]['real_distance_km'].iloc[0]

total_ash = 0
total_lipid = 0
total_protein = 0
total_sludge = Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year']+Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']

for i in range(len(Urbana_Champaign)):

    if Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_value=0.468
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
        
    elif Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_value=0.414
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
        
    else:
        sludge_ash_value=0.231
        sludge_lipid_value=0.206
        sludge_protein_value=0.456

    total_ash += sludge_ash_value*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']
    total_lipid += sludge_lipid_value*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)
    total_protein += sludge_protein_value*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)

average_ash_dw = total_ash/total_sludge
average_lipid_afdw = total_lipid/(total_sludge-total_ash)
average_protein_afdw = total_protein/(total_sludge-total_ash)

# !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion'] # set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

average_ww_2_dry_sludge_ratio = (Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] +\
                                 Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']) /\
                                 1000/365/(Urbana_Champaign.iloc[0]['flow_2022_MGD']+Urbana_Champaign.iloc[1]['flow_2022_MGD'])

state = Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001]['state'].iloc[0]

elec_GHG = float(elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['GHG (10-year median)'])

sys, barrel = create_geospatial_system(waste_cost=average_cost,
                                       waste_GHG=average_GHG,
                                       size=total_size,
                                       distance=distance_to_refinery,
                                       sludge_transportation=1,
                                       sludge_distance=16.2*Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']/(Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year']+Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']), # from Google Maps
                                       average_sludge_dw_ash=sludge_ash_values[1],
                                       average_sludge_afdw_lipid=sludge_lipid_values[1],
                                       average_sludge_afdw_protein=sludge_protein_values[1],
                                       anaerobic_digestion=None,
                                       aerobic_digestion=None,
                                       ww_2_dry_sludge_ratio=average_ww_2_dry_sludge_ratio,
                                       # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                       state=state,
                                       elec_GHG=elec_GHG)

model = create_geospatial_model(system=sys,
                                sludge_ash=sludge_ash_values,
                                sludge_lipid=sludge_lipid_values,
                                sludge_protein=sludge_protein_values,
                                raw_wastewater_price_baseline=sys.flowsheet.stream.sludge_assumed_in_wastewater.price,
                                biocrude_and_transportation_price=[sys.flowsheet.stream.biocrude.price/6.80*4.21, # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                                                   sys.flowsheet.stream.biocrude.price,
                                                                   sys.flowsheet.stream.biocrude.price/6.80*11.9], # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                electricity_cost=[elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['price (10-year 5th)'].iloc[0]/100,
                                                  elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['price (10-year median)'].iloc[0]/100,
                                                  elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['price (10-year 95th)'].iloc[0]/100],
                                electricity_GHG=[0.67848*elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['GHG (10-year 5th)']/elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['GHG (10-year median)'],
                                                  0.67848,
                                                  0.67848*elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['GHG (10-year 95th)']/elec[elec['state']==Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001].iloc[0]['state']]['GHG (10-year median)']])

kwargs = {'N':1000, 'rule':'L', 'seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()

idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
results = model.table.iloc[:, idx:]

CU_combined_uncertainty_NPV['combined'] = results[('Geospatial','NPV [$]')]
CU_combined_uncertainty_biocrude['combined'] = results[('Geospatial','Biocrude production [BPD]')]
CU_combined_uncertainty_decarbonization['combined'] = results[('Geospatial','Decarbonization amount [tonne_per_day]')]
CU_combined_uncertainty_decarbonization_cost['combined'] = -results[('Geospatial','NPV [$]')]/results[('Geospatial','Decarbonization amount [tonne_per_day]')]/365/30

CU_NPV_results = pd.concat([CU_uncertainty_NPV, CU_combined_uncertainty_NPV], axis=1)
CU_biocrude_results = pd.concat([CU_uncertainty_biocrude, CU_combined_uncertainty_biocrude], axis=1)
CU_decarbonization_results = pd.concat([CU_uncertainty_decarbonization, CU_combined_uncertainty_decarbonization], axis=1)
CU_decarbonization_cost_results = pd.concat([CU_uncertainty_decarbonization_cost, CU_combined_uncertainty_decarbonization_cost], axis=1)

CU_NPV_results.to_excel(folder + f'results/Urbana-Champaign/NPV_{date.today()}.xlsx')
CU_biocrude_results.to_excel(folder + f'results/Urbana-Champaign/biocrude_{date.today()}.xlsx')
CU_decarbonization_results.to_excel(folder + f'results/Urbana-Champaign/decarbonization_{date.today()}.xlsx')
CU_decarbonization_cost_results.to_excel(folder + f'results/Urbana-Champaign/decarbonization_cost_{date.today()}.xlsx')

fig = plt.figure(figsize=(10, 10))

gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)

def set_cf_plot():
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 38
    plt.rcParams['ytick.labelsize'] = 38
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

def add_region(position, xlabel, color):
    ax = fig.add_subplot(gs[0, position])
    
    set_cf_plot()
    
    ax = plt.gca()
    ax.set_ylim([-20, 10])
    ax.set_xlabel(xlabel, fontname='Arial', fontsize=38, labelpad=15)
    
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)
        ax.set_ylabel(r'$\mathbf{Net\ present\ value}$ [MM\$]', fontname='Arial', fontsize=45)
    
    elif position == 2:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
    
    else:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
    
    bp = ax.boxplot(CU_NPV_results.iloc[:,position]/1000000, showfliers=False, widths=0.7, patch_artist=True)
    
    for box in bp['boxes']:
        box.set(color='k', facecolor=color, linewidth=3)

    for whisker in bp['whiskers']:
        whisker.set(color='k', linewidth=3)

    for median in bp['medians']:
        median.set(color='k', linewidth=3)
        
    for cap in bp['caps']:
        cap.set(color='k', linewidth=3)
        
    for flier in bp['fliers']:
        flier.set(marker='o', markersize=7, markerfacecolor=color, markeredgewidth=1.5)

# TODO: remove the middle y-axes using PS
# TODO: add y=0 line in PowerPoint
add_region(0, 'Northeast', o)
add_region(1, 'Southwest', o)
add_region(2, 'Combined', b)

#%% sludge transportation (satellite, SA) (data preparation)

# import decarbonization_result from #%% read decarbonization data
Minnesota = decarbonization_result[decarbonization_result['state']=='MN']

Minnesota.sort_values(by='total_sludge_amount_kg_per_year', ascending=False, inplace=True)

center_WRRF = Minnesota.iloc[0,:]

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')

# get all WRRFs in MN, WI, IA
MN_WI_IA = WRRF_input[WRRF_input['state'].isin(['MN','WI','IA'])]

# select WRRFs that are within 200 km (linear distance) from the center WRRF

MN_WI_IA_distance = []
for i in range(len(MN_WI_IA)):
    MN_WI_IA_distance.append(geopy.distance.geodesic((MN_WI_IA['latitude'].iloc[i], MN_WI_IA['longitude'].iloc[i]), (center_WRRF.latitude, center_WRRF.longitude)).km)

MN_WI_IA['linear_distance_to_center_WRRF_km'] = MN_WI_IA_distance

MN_WI_IA = MN_WI_IA[MN_WI_IA['linear_distance_to_center_WRRF_km'] <= 200]

# remove WRRFs that are already decarbonized without sludge transportation
MN_WI_IA_decarbonized = decarbonization_result[decarbonization_result['state'].isin(['MN','WI','IA'])]

MN_WI_IA_transportation = MN_WI_IA[~MN_WI_IA['facility_code'].isin(list(MN_WI_IA_decarbonized['facility_code']))]

# calculate the travel distance using Google Maps API
# !!! get a google API key !!! do not upload to GitHub
gmaps = googlemaps.Client(key='XXX')
real_distance = []
for i in range(len(MN_WI_IA_transportation.index)):
    try:
        print(i)
        real_distance.append(gmaps.distance_matrix((MN_WI_IA_transportation['latitude'].iloc[i], MN_WI_IA_transportation['longitude'].iloc[i]),
                                                   (center_WRRF.latitude, center_WRRF.longitude),
                                                   mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000)
    except KeyError:
        print('--------------------------------')
        WRRF_input['real_distance_km'].iloc[i] = np.nan

MN_WI_IA_transportation['real_distance_to_center_WRRF_km'] = real_distance

MN_WI_IA_transportation.to_excel(folder + f'results/Minnesota/satellite_WRRFs_{date.today()}.xlsx')

#%% sludge transportation (satellite, SA) (center WRRF uncertainty)

filterwarnings('ignore')

Minnesota = decarbonization_result[decarbonization_result['state']=='MN']

Minnesota.sort_values(by='total_sludge_amount_kg_per_year', ascending=False, inplace=True)

elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

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

Minnesota['waste_cost'] = sum(Minnesota[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/Minnesota['total_sludge_amount_kg_per_year']*1000

Minnesota['waste_GHG'] =  sum(Minnesota[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/Minnesota['total_sludge_amount_kg_per_year']*1000

sys, barrel = create_geospatial_system(waste_cost=Minnesota.iloc[0]['waste_cost'],
                                       waste_GHG=Minnesota.iloc[0]['waste_GHG'],
                                       size=Minnesota.iloc[0]['flow_2022_MGD'],
                                       distance=Minnesota.iloc[0]['real_distance_km'],
                                       anaerobic_digestion=Minnesota.iloc[0]['sludge_anaerobic_digestion'],
                                       aerobic_digestion=Minnesota.iloc[0]['sludge_aerobic_digestion'],
                                       ww_2_dry_sludge_ratio=Minnesota.iloc[0]['total_sludge_amount_kg_per_year']/1000/365/Minnesota.iloc[0]['flow_2022_MGD'],
                                       # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                       state=Minnesota.iloc[0]['state'],
                                       elec_GHG=float(elec[elec['state']==Minnesota.iloc[0]['state']]['GHG (10-year median)']))

if Minnesota.iloc[0]['sludge_aerobic_digestion'] == 1:
    sludge_ash_values=[0.374, 0.468, 0.562, 'aerobic_digestion']
    sludge_lipid_values=[0.154, 0.193, 0.232]
    sludge_protein_values=[0.408, 0.510, 0.612]
    
elif Minnesota.iloc[0]['sludge_anaerobic_digestion'] == 1:
    sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
    sludge_lipid_values=[0.154, 0.193, 0.232]
    sludge_protein_values=[0.408, 0.510, 0.612]
    
else:
    sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
    sludge_lipid_values=[0.080, 0.206, 0.308]
    sludge_protein_values=[0.380, 0.456, 0.485]

model = create_geospatial_model(system=sys,
                                sludge_ash=sludge_ash_values,
                                sludge_lipid=sludge_lipid_values,
                                sludge_protein=sludge_protein_values,
                                raw_wastewater_price_baseline=sys.flowsheet.stream.sludge_assumed_in_wastewater.price,
                                biocrude_and_transportation_price=[sys.flowsheet.stream.biocrude.price/6.80*4.21, # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                                                   sys.flowsheet.stream.biocrude.price,
                                                                   sys.flowsheet.stream.biocrude.price/6.80*11.9], # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                electricity_cost=[elec[elec['state']==Minnesota.iloc[0]['state']]['price (10-year 5th)'].iloc[0]/100,
                                                  elec[elec['state']==Minnesota.iloc[0]['state']]['price (10-year median)'].iloc[0]/100,
                                                  elec[elec['state']==Minnesota.iloc[0]['state']]['price (10-year 95th)'].iloc[0]/100],
                                electricity_GHG=[0.67848*elec[elec['state']==Minnesota.iloc[0]['state']]['GHG (10-year 5th)']/elec[elec['state']==Minnesota.iloc[0]['state']]['GHG (10-year median)'],
                                                 0.67848,
                                                 0.67848*elec[elec['state']==Minnesota.iloc[0]['state']]['GHG (10-year 95th)']/elec[elec['state']==Minnesota.iloc[0]['state']]['GHG (10-year median)']])

kwargs = {'N':1000, 'rule':'L', 'seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()

idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
results = model.table.iloc[:, idx:]

biocrude_center = results[('Geospatial','Biocrude production [BPD]')]
decarbonization_center = results[('Geospatial','Decarbonization amount [tonne_per_day]')]

#%% sludge transportation (satellite, SA) (all WRRFs uncertainty)

MN_WI_IA_transportation = pd.read_excel(folder + 'results/Minnesota/satellite_WRRFs_2024-01-04.xlsx')

filterwarnings('ignore')

print(len(MN_WI_IA_transportation))

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

MN_WI_IA_transportation['waste_cost'] = sum(MN_WI_IA_transportation[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/MN_WI_IA_transportation['total_sludge_amount_kg_per_year']*1000

MN_WI_IA_transportation['waste_GHG'] =  sum(MN_WI_IA_transportation[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/MN_WI_IA_transportation['total_sludge_amount_kg_per_year']*1000

total_cost_except_center = sum(MN_WI_IA_transportation.iloc[i]['waste_cost']*MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year'] for i in range(len(MN_WI_IA_transportation)))
total_GHG_except_center = sum(MN_WI_IA_transportation.iloc[i]['waste_GHG']*MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year'] for i in range(len(MN_WI_IA_transportation)))
total_sludge_except_center = sum(MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year'] for i in range(len(MN_WI_IA_transportation)))

total_cost = total_cost_except_center + Minnesota.iloc[0]['waste_cost']*Minnesota.iloc[0]['total_sludge_amount_kg_per_year']
total_GHG = total_GHG_except_center + Minnesota.iloc[0]['waste_GHG']*Minnesota.iloc[0]['total_sludge_amount_kg_per_year']
total_sludge = total_sludge_except_center + Minnesota.iloc[0]['total_sludge_amount_kg_per_year']

average_cost = total_cost/total_sludge

average_GHG = total_GHG/total_sludge

total_size = sum(MN_WI_IA_transportation.iloc[i]['flow_2022_MGD'] for i in range(len(MN_WI_IA_transportation))) + Minnesota.iloc[0]['flow_2022_MGD']

distance_to_refinery = Minnesota.iloc[0]['real_distance_km']

sludge_transporataion_distance = sum(MN_WI_IA_transportation.iloc[i]['real_distance_to_center_WRRF_km']*MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year'] for i in range(len(MN_WI_IA_transportation)))/total_sludge

total_ash_except_center = 0
total_lipid_except_center = 0
total_protein_except_center = 0

for i in range(len(MN_WI_IA_transportation)):

    if MN_WI_IA_transportation.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_value=0.468
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
        
    elif MN_WI_IA_transportation.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_value=0.414
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
        
    else:
        sludge_ash_value=0.231
        sludge_lipid_value=0.206
        sludge_protein_value=0.456

    total_ash_except_center += sludge_ash_value*MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year']
    total_lipid_except_center += sludge_lipid_value*MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)
    total_protein_except_center += sludge_protein_value*MN_WI_IA_transportation.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)

# manually checked: the center WRRF does not have anaerobic or aerobic digestion

total_ash = total_ash_except_center + Minnesota.iloc[0]['total_sludge_amount_kg_per_year']*0.231
total_lipid = total_lipid_except_center + Minnesota.iloc[0]['total_sludge_amount_kg_per_year']*(1-0.231)*0.206
total_protein = total_protein_except_center + Minnesota.iloc[0]['total_sludge_amount_kg_per_year']*(1-0.231)*0.456

average_ash_dw = total_ash/total_sludge
average_lipid_afdw = total_lipid/(total_sludge-total_ash)
average_protein_afdw = total_protein/(total_sludge-total_ash)

# !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion'] # set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

average_ww_2_dry_sludge_ratio = total_sludge/1000/365/total_size

state = 'MN'

elec_GHG = float(elec[elec['state']==state]['GHG (10-year median)'])

sys, barrel = create_geospatial_system(waste_cost=average_cost,
                                       waste_GHG=average_GHG,
                                       size=total_size,
                                       distance=distance_to_refinery,
                                       sludge_transportation=1,
                                       sludge_distance=sludge_transporataion_distance,
                                       average_sludge_dw_ash=sludge_ash_values[1],
                                       average_sludge_afdw_lipid=sludge_lipid_values[1],
                                       average_sludge_afdw_protein=sludge_protein_values[1],
                                       anaerobic_digestion=None,
                                       aerobic_digestion=None,
                                       ww_2_dry_sludge_ratio=average_ww_2_dry_sludge_ratio,
                                       # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                       state=state,
                                       elec_GHG=elec_GHG)

model = create_geospatial_model(system=sys,
                                sludge_ash=sludge_ash_values,
                                sludge_lipid=sludge_lipid_values,
                                sludge_protein=sludge_protein_values,
                                raw_wastewater_price_baseline=sys.flowsheet.stream.sludge_assumed_in_wastewater.price,
                                biocrude_and_transportation_price=[sys.flowsheet.stream.biocrude.price/6.80*4.21, # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                                                   sys.flowsheet.stream.biocrude.price,
                                                                   sys.flowsheet.stream.biocrude.price/6.80*11.9], # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                electricity_cost=[elec[elec['state']==state]['price (10-year 5th)'].iloc[0]/100,
                                                  elec[elec['state']==state]['price (10-year median)'].iloc[0]/100,
                                                  elec[elec['state']==state]['price (10-year 95th)'].iloc[0]/100],
                                electricity_GHG=[0.67848*elec[elec['state']==state]['GHG (10-year 5th)']/elec[elec['state']==state]['GHG (10-year median)'],
                                                  0.67848,
                                                  0.67848*elec[elec['state']==state]['GHG (10-year 95th)']/elec[elec['state']==state]['GHG (10-year median)']])

kwargs = {'N':1000, 'rule':'L', 'seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()

idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
all_WRRFs_results = model.table.iloc[:, idx:]

biocrude_all = all_WRRFs_results[('Geospatial','Biocrude production [BPD]')]
decarbonization_all = all_WRRFs_results[('Geospatial','Decarbonization amount [tonne_per_day]')]

SA_results = pd.DataFrame({'center_decarbonization': decarbonization_center,
                           'all_decarbonization': decarbonization_all,
                           'center_biocrude': biocrude_center,
                           'all_biocrude': biocrude_all})

SA_results.to_excel(folder + f'results/Minnesota/center_vs_all_{date.today()}.xlsx')

#%% sludge transportation (satellite, SA) (visualization)

SA_results = pd.read_excel(folder + f'results/Minnesota/center_vs_all_2024-01-05.xlsx')

SA_results.drop('Unnamed: 0', axis=1, inplace=True)

fig, ax = plt.subplots(figsize = (10, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax.set_xlim([20, 120])
ax.set_ylim([300, 800])

ax.tick_params(direction='inout', length=15, width=3, bottom=True, top=False, left=True, right=False, pad=6)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(20, 140, 20))
plt.yticks(np.arange(300, 900, 100))

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(20, 140, 20))
plt.yticks(np.arange(300, 900, 100))

ax.set_xlabel(r'$\mathbf{Decarbonization\ amount}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Biocrude\ production}$ [BPD]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

def add_rectangle(item, color, edgecolor):
    rectangle_fill = Rectangle((SA_results[f'{item}_decarbonization'].quantile(0.05), SA_results[f'{item}_biocrude'].quantile(0.05)),
                                SA_results[f'{item}_decarbonization'].quantile(0.95) - SA_results[f'{item}_decarbonization'].quantile(0.05),
                                SA_results[f'{item}_biocrude'].quantile(0.95) - SA_results[f'{item}_biocrude'].quantile(0.05),
                                fc=color, alpha=0.8)
    ax.add_patch(rectangle_fill)
    
    rectangle_edge = Rectangle((SA_results[f'{item}_decarbonization'].quantile(0.05), SA_results[f'{item}_biocrude'].quantile(0.05)),
                                SA_results[f'{item}_decarbonization'].quantile(0.95) - SA_results[f'{item}_decarbonization'].quantile(0.05),
                                SA_results[f'{item}_biocrude'].quantile(0.95) - SA_results[f'{item}_biocrude'].quantile(0.05),
                                color=edgecolor, lw=3, fc='none', alpha=1)
    ax.add_patch(rectangle_edge)
    
def add_line(item, color):
    plt.plot([SA_results[f'{item}_decarbonization'].quantile(0.25), SA_results[f'{item}_decarbonization'].quantile(0.75)],
             [SA_results[f'{item}_biocrude'].quantile(0.5), SA_results[f'{item}_biocrude'].quantile(0.5)],
             lw=3, color=color, solid_capstyle='round', zorder=1)
    plt.plot([SA_results[f'{item}_decarbonization'].quantile(0.5), SA_results[f'{item}_decarbonization'].quantile(0.5)],
             [SA_results[f'{item}_biocrude'].quantile(0.25), SA_results[f'{item}_biocrude'].quantile(0.75)],
             lw=3, color=color, solid_capstyle='round', zorder=1)

def add_point(item, edgecolor):
    ax_top.scatter(x=SA_results[f'{item}_decarbonization'].quantile(0.5),
                   y=SA_results[f'{item}_biocrude'].quantile(0.5),
                   marker='s',
                   s=200,
                   c='w',
                   linewidths=3,
                   alpha=1,
                   edgecolor=edgecolor)
    
add_rectangle('center', o, do)
add_line('center', do)

add_rectangle('all', b, db)
add_line('all', db)

add_point('center', do)
add_point('all', db)

#%% sludge transportation (heat map, HM)

filterwarnings('ignore')

WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2023-12-26.xlsx')
elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')

# if just want to see the two plants in Urbana-Champaign:
# WRRF_input = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

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

waste_cost = sludge_disposal_cost['other_sludge_management']*1000

waste_GHG =  sludge_emission_factor['other_sludge_management']*1000

# use data for sludge + biosolids (as we did in the HTL model paper)
HM_sludge_ash_values = [0.174,0.257,0.414, 'digestion'] # set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
HM_sludge_lipid_values = [0.08,0.204,0.308]
HM_sludge_protein_values = [0.38,0.463,0.51]

HM_elec_GHG = elec['GHG (10-year median)'].mean()

ternary_results_dict = {'sludge_amount':[], 'sludge_transportation_distance':[],
                        'NPV_5th':[], 'NPV_50th':[], 'NPV_95th':[],
                        'biocrude_production_5th':[], 'biocrude_production_50th':[], 'biocrude_production_95th':[],
                        'decarbonization_amount_5th':[], 'decarbonization_amount_50th':[], 'decarbonization_amount_95th':[],
                        'decarbonization_cost_5th':[], 'decarbonization_cost_50th':[], 'decarbonization_cost_95th':[]}

ternary_results = pd.DataFrame(ternary_results_dict)

get_quantiles = lambda data, quantiles=(0.05, 0.5, 0.95): [data.quantile(q) for q in quantiles]

for size in np.linspace(2, 20, 10):
    for sludge_distance in np.linspace(20, 200, 10):
        
        print('\n\n', f'sludge amount: {size} metric tonne/day\n', f'sludge travel distance: {sludge_distance} km\n')

        sys, barrel = create_geospatial_system(waste_cost=waste_cost,
                                               waste_GHG=waste_GHG,
                                               size=size,
                                               distance=207, # 207 km is the median travel distance from WRRFs to oil refineries
                                               sludge_transportation=1,
                                               sludge_distance=sludge_distance,
                                               average_sludge_dw_ash=HM_sludge_ash_values[1],
                                               average_sludge_afdw_lipid=HM_sludge_lipid_values[1],
                                               average_sludge_afdw_protein=HM_sludge_protein_values[1],
                                               anaerobic_digestion=None,
                                               aerobic_digestion=None,
                                               ww_2_dry_sludge_ratio=1,
                                               # ww_2_dry_sludge_ratio: how much metric tonne/day sludge can be produced by 1 MGD of ww
                                               state='average',
                                               elec_GHG=HM_elec_GHG)
        
        model = create_geospatial_model(system=sys,
                                        sludge_ash=HM_sludge_ash_values,
                                        sludge_lipid=HM_sludge_lipid_values,
                                        sludge_protein=HM_sludge_protein_values,
                                        raw_wastewater_price_baseline=sys.flowsheet.stream.sludge_assumed_in_wastewater.price,
                                        biocrude_and_transportation_price=[sys.flowsheet.stream.biocrude.price/6.80*4.21, # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                                                           sys.flowsheet.stream.biocrude.price,
                                                                           sys.flowsheet.stream.biocrude.price/6.80*11.9], # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
                                        electricity_cost=[elec['price (10-year 5th)'].mean()/100,
                                                          elec['price (10-year median)'].mean()/100,
                                                          elec['price (10-year 95th)'].mean()/100],
                                        electricity_GHG=[0.67848*elec['GHG (10-year 5th)'].mean()/elec['GHG (10-year median)'].mean(),
                                                          0.67848,
                                                          0.67848*elec['GHG (10-year 95th)'].mean()/elec['GHG (10-year median)'].mean()])
        
        kwargs = {'N':1000, 'rule':'L', 'seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        HM_NPV = model.table['Geospatial']['NPV [$]'].dropna()
        HM_biocrude = model.table['Geospatial']['Biocrude production [BPD]'].dropna()
        HM_decarbonization = model.table['Geospatial']['Decarbonization amount [tonne_per_day]'].dropna()
        HM_decarbonization_cost = -model.table['Geospatial']['NPV [$]'].dropna()/model.table['Geospatial']['Decarbonization amount [tonne_per_day]'].dropna()/365/30
        
        ternary_results.loc[len(ternary_results.index)] = (
            [size, sludge_distance,] +
            get_quantiles(HM_NPV) +
            get_quantiles(HM_biocrude) +
            get_quantiles(HM_decarbonization) +
            get_quantiles(HM_decarbonization_cost)
            )
        
ternary_results.to_excel(folder + f'results/heat_map_{size}_tonne_per_day_{sludge_distance}_km_{date.today()}.xlsx')

#%% sludge transportation (heat map, HM) visualization (decarbonization) 

# HM = pd.read_excel(folder + 'results/heat_map_20.0_tonne_per_day_1000.0_km_2024-01-03.xlsx')
HM = pd.read_excel(folder + 'results/heat_map_20.0_tonne_per_day_200.0_km_2024-01-05.xlsx')

HM_decarbonization = HM[['sludge_amount','sludge_transportation_distance','decarbonization_amount_50th']]

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, pad=6)

ax.set_xlabel(r'$\mathbf{Average\ distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Total\ sludge}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

# plt.xticks(np.arange(100, 1150, 150))
plt.xticks(np.arange(20, 220, 30))
plt.yticks(np.arange(2, 22, 3))

ax_right = ax.twinx()
ax_right.set_ylim((2, 20))
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

# plt.xticks(np.arange(100, 1150, 150))
plt.xticks(np.arange(20, 220, 30))
plt.yticks(np.arange(2, 22, 3))

ax_top = ax.twiny()
# ax_top.set_xlim((100, 1000))
ax_top.set_xlim((20, 200))
# plt.xticks(np.arange(max_distance*0.2, max_distance*1.2, max_distance*0.2))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

# plt.xticks(np.arange(100, 1150, 150))
plt.xticks(np.arange(20, 220, 30))
plt.yticks(np.arange(2, 22, 3))

try:
    color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [r, o, y, g, b])
except ValueError:
    pass

x = np.array(HM_decarbonization['sludge_transportation_distance'])
y = np.array(HM_decarbonization['sludge_amount'])
z = np.array(HM_decarbonization['decarbonization_amount_50th'])

fills = ax.tricontourf(x, y, z, levels=10000, 
                      cmap=color_map_Guest)
fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(x, y, z, levels=7, linewidths=3, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=38)

#%% sludge transportation (heat map, HM) visualization (NPV) 

# HM = pd.read_excel(folder + 'results/heat_map_20.0_tonne_per_day_1000.0_km_2024-01-03.xlsx')
HM = pd.read_excel(folder + 'results/heat_map_20.0_tonne_per_day_200.0_km_2024-01-05.xlsx')

HM_NPV = HM[['sludge_amount','sludge_transportation_distance','NPV_50th']]

HM_NPV['NPV_50th'] = HM_NPV['NPV_50th']/1000000

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, pad=6)

ax.set_xlabel(r'$\mathbf{Average\ distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Total\ sludge}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

# plt.xticks(np.arange(100, 1150, 150))
plt.xticks(np.arange(20, 220, 30))
plt.yticks(np.arange(2, 22, 3))

ax_right = ax.twinx()
ax_right.set_ylim((2, 20))
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

# plt.xticks(np.arange(100, 1150, 150))
plt.xticks(np.arange(20, 220, 30))
plt.yticks(np.arange(2, 22, 3))

ax_top = ax.twiny()
# ax_top.set_xlim((100, 1000))
ax_top.set_xlim((20, 200))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=7.5, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

# plt.xticks(np.arange(100, 1150, 150))
plt.xticks(np.arange(20, 220, 30))
plt.yticks(np.arange(2, 22, 3))

try:
    color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [r, o, y, g, b])
except ValueError:
    pass

x = np.array(HM_NPV['sludge_transportation_distance'])
y = np.array(HM_NPV['sludge_amount'])
z = np.array(HM_NPV['NPV_50th'])

fills = ax.tricontourf(x, y, z, levels=10000, 
                      cmap=color_map_Guest)
fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(x, y, z, levels=7, linewidths=3, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=38)