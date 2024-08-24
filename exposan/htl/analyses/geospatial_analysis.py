#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Sun Jun 11 08:12:41 2023

@author: jiananfeng
'''

# TODO 1: consider adding transportation cost and carbon intensity as contextual parameters
# see if https://www.bls.gov/oes/tables.htm (from Abby) has some information
# TODO 2: update electricity CI to match up balancing areas instead of states
# TODO 3: refactor/reformat all codes

# confirmed: both 'markersize' (in gpd.plot()) and 's' (in plt.scatter()) are proportional to area

#%% initialization
import geopy.distance, googlemaps, random
import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt, matplotlib.colors as colors, matplotlib.ticker as mtick
from matplotlib.mathtext import _mathtext as mathtext
from matplotlib.patches import Rectangle
from colorpalette import Color
from exposan.htl import create_geospatial_system, create_geospatial_model
from qsdsan.utils import auom, palettes
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

db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

# kg/L
water_density = 1

_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_MMgal_to_L = auom('gal').conversion_factor('L')*1000000

# !!! when creating this input file, change year to 2022 in IEDO code for electricity GHG calculation
WRRF = pd.read_excel(folder + 'HTL_geospatial_input_08202024.xlsx')

assert WRRF.duplicated(subset='CWNS_NUM').sum() == 0

WRRF['total_sludge_amount_kg_per_year'] = WRRF[['landfill','land_application',
                                                'incineration']].sum(axis=1)

WRRF['biosolids_emission'] = WRRF[['LF_CH4','LA_N2O']].sum(axis=1)

# all emission data in kg CO2 eq/day
WRRF['total_emission'] = WRRF[['LF_CH4','LA_N2O','electricity_emission',
                               'onsite_NG_emission','upstream_NG_emission',
                               'CH4_emission','N2O_emission','CO2_emission']].sum(axis=1)

treatment_trains = np.array(['LAGOON_AER','LAGOON_ANAER','LAGOON_FAC',
                             'LAGOON_UNCATEGORIZED','C1','C2','C3','C5',
                             'C6','B1','B1E','B2','B3','B4','B5','B6',
                             'D1','D2','D3','D5','D6','E2','E2P','F1',
                             'I1','I1E','I2','I3','I5','I6','G1','G1E',
                             'G2','G3','G5','G6','H1','N1','N2','O1',
                             'O1E','O2','O3','O5','O6'], dtype=object)

TT_indentifier = WRRF[treatment_trains].apply(lambda x: x > 0)
WRRF['treatment_train'] = TT_indentifier.apply(lambda x: list(treatment_trains[x.values]), axis=1)

WRRF = WRRF[['FACILITY','CITY','STATE','CWNS_NUM','FACILITY_ID','LATITUDE',
             'LONGITUDE','FLOW_2022_MGD_FINAL','balancing_area','kg_CO2_kWh',
             'treatment_train','landfill','land_application','incineration',
             'total_sludge_amount_kg_per_year','CH4_emission','N2O_emission',
             'CO2_emission','electricity_emission','onsite_NG_emission',
             'upstream_NG_emission','biosolids_emission','total_emission']]

# assume LAGOON_UNCATEGORIZED has AeD (resulting in higher ash content than AD) to be conservative
TT_w_AeD = ['B2','C2','D2','E2','E2P','G2','I2','N2','O2','LAGOON_AER','LAGOON_UNCATEGORIZED']
TT_w_AD = ['B1','B1E','B4','C1','D1','F1','G1','G1E','H1','I1','I1E','N1','O1','O1E','LAGOON_ANAER','LAGOON_FAC']

# TODO: check this, and how this have been used?
# to be conservative (AeD has higher ash content than AD):
# if a WRRF has AeD (regardless of AD), assume all sludge composition follows the AeD sludge composition,
# if no AeD but has AD, then use AD composition
# if none, then use sludge composition
WRRF.loc[WRRF['treatment_train'].apply(lambda x: len([i for i in TT_w_AeD if i in x]) > 0), 'sludge_aerobic_digestion'] = 1
WRRF.loc[WRRF['treatment_train'].apply(lambda x: len([i for i in TT_w_AD if i in x]) > 0), 'sludge_anaerobic_digestion'] = 1
WRRF.fillna({'sludge_aerobic_digestion': 0}, inplace=True)
WRRF.fillna({'sludge_anaerobic_digestion': 0}, inplace=True)

WRRF = WRRF.rename({'FACILITY':'facility',
                    'CITY':'city',
                    'STATE':'state',
                    'CWNS_NUM':'CWNS',
                    'FACILITY_ID':'facility_ID',
                    'LATITUDE':'latitude',
                    'LONGITUDE':'longitude',
                    'FLOW_2022_MGD_FINAL':'flow_2022_MGD'}, axis=1)

# TODO: check all 4269 in this file and IEDO analysis (IEDO_GHG_inventory_visualization_JF.py and balancing_area.py)
WRRF = WRRF.sort_values(by='flow_2022_MGD', ascending=False)

WRRF = gpd.GeoDataFrame(WRRF, crs='EPSG:4269',
                        geometry=gpd.points_from_xy(x=WRRF.longitude,
                                                    y=WRRF.latitude))

refinery = pd.read_csv(folder + 'refinery/petroleum_refineries_EIA_06062024.csv')
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
# confirmed from the wesbite (next line): there is no change in US map since June 1, 1995. So, using 2018 US map data is OK.
# https://en.wikipedia.org/wiki/Territorial_evolution_of_the_United_States#1946%E2%80%93present_(Decolonization)
US = US[['NAME','geometry']]

for excluded in ('Alaska',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'Hawaii',
                 'Puerto Rico',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

balnc_area = gpd.read_file(folder + 'lpreg2/lpreg2.shp')

WRRF = WRRF.to_crs(crs='EPSG:3857')
refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')
balnc_area = balnc_area.to_crs(crs='EPSG:3857')

# a small number of WRRFs fall outside the boundary of the contiguous U.S. a little bit
# do not use sjoin for WRRFs, as this will remove the those WRRFs
# WRRF = gpd.sjoin(WRRF, US)
# WRRF = WRRF.drop(['index_right'], axis=1)

# we can use sjoin for refinery, this only removes 6 refineries in Alaska and Hawaii
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

# this one is the same as the one in exposan/htl/data
elec_price = pd.read_excel(folder + 'state_elec_price_2022.xlsx', 'elec_price_2022')
elec_price = elec_price.merge(US, how='right', left_on='name', right_on='NAME')
elec_price = gpd.GeoDataFrame(elec_price)

elec_GHG = WRRF[['balancing_area','kg_CO2_kWh']].copy()
elec_GHG.drop_duplicates(inplace=True)
elec_GHG = elec_GHG.merge(balnc_area, how='left', left_on='balancing_area', right_on='PCA_REG')
elec_GHG = gpd.GeoDataFrame(elec_GHG)

state_PADD = {'Alabama': 3,
              'Arizona': 5,
              'Arkansas': 3,
              'California': 5,
              'Colorado': 4,
              'Connecticut': 1,
              'Delaware': 1,
              'District of Columbia': 1,
              'Florida': 1,
              'Georgia': 1,
              'Idaho': 4,
              'Illinois': 2,
              'Indiana': 2,
              'Iowa': 2,
              'Kansas': 2,
              'Kentucky': 2,
              'Louisiana': 3,
              'Michigan': 2,
              'Maine': 1,
              'Maryland': 1,
              'Washington': 5,
              'Massachusetts': 1,
              'Minnesota': 2,
              'Mississippi': 3,
              'Missouri': 2,
              'Montana': 4,
              'Nebraska': 2,
              'Nevada': 5,
              'New Hampshire': 1,
              'New Jersey': 1,
              'New Mexico': 3,
              'New York': 1,
              'North Carolina': 1,
              'North Dakota': 2,
              'Ohio': 2,
              'Oklahoma': 2,
              'Oregon': 5,
              'Pennsylvania': 1,
              'Rhode Island': 1,
              'South Carolina': 1,
              'South Dakota': 2,
              'Tennessee': 2,
              'Texas': 3,
              'Utah': 4,
              'Vermont': 1,
              'Virginia': 1,
              'West Virginia': 1,
              'Wisconsin': 2,
              'Wyoming': 4}

PADD_color = {1: b,
              2: g,
              3: r,
              4: o,
              5: y}

# sludge disposal cost in $/kg sludge (Peccia and Westerhoff. 2015, https://pubs.acs.org/doi/full/10.1021/acs.est.5b01931)
sludge_disposal_cost = {'landfill': 0.413,
                        'land_application': 0.606,
                        'incineration': 0.441}

# sludge emission factor in kg CO2 eq/kg sludge (values from IEDO)
sludge_emission_factor = {'landfill': 5.65/1000*29.8,
                          'land_application': 0.05*0.01*44/28*273,
                          'incineration': 0}

#%% WRRFs visualization

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)
# US.plot(ax=ax, color='w', edgecolor='k', linewidth=6) # for all WRRFs together with the same symbols

WRRF = WRRF.sort_values(by='flow_2022_MGD', ascending=False)

WRRF_flow = WRRF['flow_2022_MGD']

more_than_100 = WRRF_flow > 100
WRRF[more_than_100].plot(ax=ax, color=dg, markersize=WRRF.loc[more_than_100, WRRF_flow.name]*10, edgecolor='k', linewidth=1.5, alpha=1)

between_10_and_100 = (WRRF_flow > 10) & (WRRF_flow <= 100)
WRRF[between_10_and_100].plot(ax=ax, color=g, markersize=WRRF.loc[between_10_and_100, WRRF_flow.name]*10, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_10 = WRRF_flow <= 10
WRRF[less_than_10].plot(ax=ax, color=g, markersize=WRRF.loc[less_than_10, WRRF_flow.name]*10, edgecolor='none', alpha=0.1)

ax.set_aspect(1)

ax.set_axis_off()

# comment out the code above and uncomment the following line to show all WRRFs together with the same symbols
# WRRF.plot(ax=ax, color=a, markersize=10)

#%% oil refinery visualization

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax,
        color=[PADD_color[state_PADD[i]] for i in list(US.NAME)],
        alpha=0.3,
        edgecolor='k',
        linewidth=0)

# for TOC
# US.plot(ax=ax, color='none', edgecolor='k', linewidth=0)

US.plot(ax=ax, color='none', edgecolor='k', linewidth=3)

 # for all oil refineries together with the same symbols
# US.plot(ax=ax, color='none', edgecolor='k', linewidth=6)

refinery['total_capacity'] = refinery[['Atmos. Crude Dist','Vacuum Dist','Catalytic Cracking',
                                       'Hydro Cracking','Thermal Cracking, Visbreaking',
                                       'Catalytic Recorming','Alkylates, Isomerization',
                                       'Desulfurization','Fluid and Delayed Coking',
                                       'Asphalt and Road Oil']].sum(axis=1)

refinery = refinery.sort_values(by='total_capacity', ascending=False)

refinery.plot(ax=ax, color=o, markersize=refinery['total_capacity']*4.5, edgecolor='k', linewidth=1.5, alpha=1)

ax.set_aspect(1)

ax.set_axis_off()

# comment out the code above and uncomment the following line to show all oil refineries together with the same symbols
# refinery.plot(ax=ax, color=o, markersize=1000, edgecolor='k', linewidth=6)

#%% WRRFs+oil refineries visualization (all)

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)
WRRF.plot(ax=ax, color=a, markersize=5)
refinery.plot(ax=ax, color=o, markersize=500, edgecolor='k', linewidth=1.5)

ax.set_aspect(1)

ax.set_axis_off()

#%% electricity price visualization

fig, ax = plt.subplots(figsize=(30, 30))

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['font.sans-serif'] = 'Arial'

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

mathtext.FontConstantsBase.sup1 = 0.35

assert elec_price.price.min() > 6, 'adjust the colormap range'
assert elec_price.price.max() < 18, 'adjust the colormap range'

norm = colors.TwoSlopeNorm(vmin=6, vcenter=12, vmax=18)
elec_price.plot('price', ax=ax, linewidth=1.5, cmap='Oranges', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35}, norm=norm)

fig.axes[1].set_ylabel('$\mathbf{Electricity\ price}$ [cent·${kWh^{-1}}$]',
                        fontname='Arial', fontsize=28)
fig.axes[1].tick_params(length=15, width=1.5)
# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% electricity carbon intensity visualization

fig, ax = plt.subplots(figsize=(30, 30))

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['font.sans-serif'] = 'Arial'

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

mathtext.FontConstantsBase.sup1 = 0.35

assert elec_GHG.kg_CO2_kWh.min() > 0, 'adjust the colormap range'
assert elec_GHG.kg_CO2_kWh.max() < 1.2, 'adjust the colormap range'

norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.6, vmax=1.2)
elec_GHG.plot('kg_CO2_kWh', ax=ax, linewidth=1.5, cmap='Blues', edgecolor='k', legend=True, legend_kwds={'shrink': 0.35}, norm=norm)

fig.axes[1].set_ylabel('$\mathbf{Electricity\ carbon\ intensity}$ [kg ${CO_2}$ eq·${kWh^{-1}}$]',
                        fontname='Arial', fontsize=28)
fig.axes[1].tick_params(length=15, width=1.5)
# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% transporation distance calculation

# no need to set a max distance since the transportation of biocrude was not a key/limiting driver
WRRF_input = WRRF.sjoin_nearest(refinery, max_distance=None, distance_col='distance')

WRRF_input['WRRF_location'] = list(zip(WRRF_input.latitude, WRRF_input.longitude))
WRRF_input['refinery_location'] = list(zip(WRRF_input.Latitude, WRRF_input.Longitude))

# =============================================================================
# previous code to calculate linear distance and set up real distance
# linear_distance = []
# for i in range(len(WRRF_input)):
#     linear_distance.append(geopy.distance.geodesic(WRRF_input['WRRF_location'].iloc[i], WRRF_input['refinery_location'].iloc[i]).km)
# 
# WRRF_input['linear_distance_km'] = linear_distance
# WRRF_input['real_distance_km'] = np.nan
# =============================================================================

distance_inventory = pd.read_excel(folder + 'distance_inventory.xlsx')

# match using WRRF ID ('CWNS') and oil refinery ID ('Site ID')
WRRF_input = WRRF_input.merge(distance_inventory, how='left', on=['CWNS','Site ID'])

missing_distance = []
for i in WRRF_input.index:
    if pd.isna(WRRF_input.loc[i, 'linear_distance_km']):
        missing_distance.append(i)

if len(missing_distance) == 0:
    WRRF_input.to_excel(folder + f'HTL_geospatial_model_input_{date.today()}_test.xlsx')
else:
    # !!! get a google API key !!! do not upload to GitHub
    gmaps = googlemaps.Client(key='XXX')
    
    for i in missing_distance:
        WRRF_input.loc[i, 'linear_distance_km'] = geopy.distance.geodesic(WRRF_input.loc[i, 'WRRF_location'],
                                                                          WRRF_input.loc[i, 'refinery_location']).km
        
        try:
            print(i)
            WRRF_input.loc[i, 'real_distance_km'] = gmaps.distance_matrix(WRRF_input.loc[i, 'WRRF_location'],
                                                                          WRRF_input.loc[i, 'refinery_location'],
                                                                          mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000
        except KeyError:
            print('--------------------------------')
            WRRF_input.loc[i, 'real_distance_km'] = np.nan
    
    # input for following analyses
    WRRF_input.to_excel(folder + f'HTL_geospatial_model_input_{date.today()}.xlsx')

#%% travel distance box plot

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# TODO: mention just the number of WRRFs included/removed in the main manuscript/SI
print(f"{WRRF_input['real_distance_km'].notna().sum()} WRRFs included")
print(f"{WRRF_input['real_distance_km'].isna().sum()} WRRFs excluded")
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

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
ax.set_ylim([0, 800])
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

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

WRRF_input.loc[WRRF_input['state'].isin(['CT','DC','DE','FL','GA','MA','MD','ME','NC','NH','NJ','NY','PA','RI','SC','VA','VT','WV']), 'WRRF_PADD'] = 1
WRRF_input.loc[WRRF_input['state'].isin(['IA','IL','IN','KS','KY','MI','MN','MO','ND','NE','OH','OK','SD','TN','WI']), 'WRRF_PADD'] = 2
WRRF_input.loc[WRRF_input['state'].isin(['AL','AR','LA','MS','NM','TX']), 'WRRF_PADD'] = 3
WRRF_input.loc[WRRF_input['state'].isin(['CO','ID','MT','UT','WY']), 'WRRF_PADD'] = 4
WRRF_input.loc[WRRF_input['state'].isin(['AZ','CA','NV','OR','WA']), 'WRRF_PADD'] = 5

fig = plt.figure(figsize=(20, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

def add_region(position, region, color):
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
    ax.set_ylim([0, 1400])
    ax.set_xlabel(region, fontname='Arial', fontsize=30, labelpad=15)
    
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)
        ax.set_ylabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=35)
    
    elif position == 4:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, length=0, labelcolor='none')
        
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

add_region(0, 'East Coast', b)
add_region(1, 'Midwest', g)
add_region(2, 'Gulf Coast', r)
add_region(3, 'Rocky Mountain', o)
add_region(4, 'West Coast', y)

#%% WRRFs GHG map

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

WRRF_input = WRRF_input.sort_values(by='total_emission', ascending=False)

WRRF_input = gpd.GeoDataFrame(WRRF_input, crs='EPSG:4269',
                               geometry=gpd.points_from_xy(x=WRRF_input.longitude,
                                                           y=WRRF_input.latitude))

WRRF_input = WRRF_input.to_crs(crs='EPSG:3857')

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)

# tonne/day
WRRF_GHG_tonne_per_day = WRRF_input['total_emission']/1000

more_than_500 = WRRF_GHG_tonne_per_day > 500
WRRF_input[more_than_500].plot(ax=ax, color=dp, markersize=WRRF_input.loc[more_than_500, WRRF_GHG_tonne_per_day.name]/350, edgecolor='k', linewidth=2, alpha=1)

between_50_and_500 = (WRRF_GHG_tonne_per_day > 50) & (WRRF_GHG_tonne_per_day <= 500)
WRRF_input[between_50_and_500].plot(ax=ax, color=p, markersize=WRRF_input.loc[between_50_and_500, WRRF_GHG_tonne_per_day.name]/350, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_50 = WRRF_GHG_tonne_per_day <= 50
WRRF_input[less_than_50].plot(ax=ax, color=p, markersize=WRRF_input.loc[less_than_50, WRRF_GHG_tonne_per_day.name]/350, alpha=0.1)

ax.set_aspect(1)

ax.set_axis_off()

#%% WRRFs sludge management GHG map

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

WRRF_input = WRRF_input.sort_values(by='biosolids_emission', ascending=False)

WRRF_input = gpd.GeoDataFrame(WRRF_input, crs='EPSG:4269',
                               geometry=gpd.points_from_xy(x=WRRF_input.longitude,
                                                           y=WRRF_input.latitude))

WRRF_input = WRRF_input.to_crs(crs='EPSG:3857')

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)

# tonne/day
sludge_GHG_tonne_per_day = WRRF_input['biosolids_emission']/1000

more_than_20 = sludge_GHG_tonne_per_day > 20
WRRF_input[more_than_20].plot(ax=ax, color=dr, markersize=WRRF_input.loc[more_than_20, sludge_GHG_tonne_per_day.name]/10, edgecolor='k', linewidth=2, alpha=1)

between_2_and_20 = (sludge_GHG_tonne_per_day > 2) & (sludge_GHG_tonne_per_day <= 20)
WRRF_input[between_2_and_20].plot(ax=ax, color=r, markersize=WRRF_input.loc[between_2_and_20, sludge_GHG_tonne_per_day.name]/10, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_2 = sludge_GHG_tonne_per_day <= 2
WRRF_input[less_than_2].plot(ax=ax, color=r, markersize=WRRF_input.loc[less_than_2, sludge_GHG_tonne_per_day.name]/10, alpha=0.1)

ax.set_aspect(1)

ax.set_axis_off()

#%% cumulative WRRFs capacity vs distances (data processing)

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

result = WRRF_input[['Site ID']].drop_duplicates()

max_distance = 1500
assert max_distance*1.2 > max_distance > WRRF_input.real_distance_km.max(), 'update max_distance'

for distance in np.linspace(0, max_distance, max_distance+1):
    WRRF_input_distance = WRRF_input[WRRF_input['real_distance_km'] <= distance]
    WRRF_input_distance = WRRF_input_distance.groupby('Site ID').sum('flow_2022_MGD')
    WRRF_input_distance = WRRF_input_distance[['flow_2022_MGD']]
    WRRF_input_distance = WRRF_input_distance.rename(columns={'flow_2022_MGD': int(distance)})
    WRRF_input_distance.reset_index(inplace=True)
    if len(WRRF_input_distance) > 0:
        result = result.merge(WRRF_input_distance, how='left', on='Site ID')

result = result.fillna(0)

if len(result) < len(refinery):
    result_index = pd.Index(result['Site ID'])
    refinery_index = pd.Index(refinery['Site ID'])
    refineries_left_id = refinery_index.difference(result_index).values

result = result.set_index('Site ID')
for item in refineries_left_id:
    result.loc[item] = [0]*len(result.columns)

result = result.merge(refinery, how='left', on='Site ID')

result.to_excel(folder + f'results/MGD_vs_real_distance_{max_distance}_km.xlsx')

#%% make the plot of cumulative WRRFs capacity vs distances (data preparation)

# !!! update the file if necessary
CF_input = pd.read_excel(folder + 'results/MGD_vs_real_distance_1500_km.xlsx')

CF_input[0] = 0

max_distance_plot = 1500
assert max_distance*1.2 > max_distance > WRRF_input.real_distance_km.max(), 'update max_distance'

CF_input = CF_input[['PADD', 'State', *list(range(0, max_distance_plot+1, 1))]]

CF_input.sort_values(by='PADD', inplace=True)

CF_input = CF_input.transpose()

#%% make the plot of cumulative WRRFs capacity vs distances (separated oil refinery)

PADD_1 = (CF_input.loc['PADD',:].isin([1,])).sum()
PADD_2 = (CF_input.loc['PADD',:].isin([1,2])).sum()
PADD_3 = (CF_input.loc['PADD',:].isin([1,2,3])).sum()
PADD_4 = (CF_input.loc['PADD',:].isin([1,2,3,4])).sum()
PADD_5 = (CF_input.loc['PADD',:].isin([1,2,3,4,5])).sum()

fig = plt.figure(figsize=(20, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

def add_region(position, start_region, end_region, color):
    ax = fig.add_subplot(gs[0, position])
    
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    for i in range(start_region, end_region):
        CF_input.iloc[2:, i].plot(ax=ax, color=color, linewidth=3)
        
    ax.set_xlim([0, max_distance_plot])
    ax.set_ylim([0, 7000])
    
    if position == 4:
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        plt.yticks(np.arange(0, 7000, 1000))
        ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, labelleft=True, labelbottom=True)
        plt.xticks(np.arange(0, max_distance_plot*1.2, max_distance_plot*0.2))
        ax.set_ylabel(r'$\mathbf{Cumulative\ WRRFs\ capacity}$ [MGD]', fontname='Arial', fontsize=35)
    else:
        if position == 2:
            ax.set_xlabel(r'$\mathbf{Travel\ distance}$ [km]', fontname='Arial', fontsize=35)
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        plt.xticks(np.arange(max_distance_plot*0.2, max_distance_plot*1.2, max_distance_plot*0.2))
    
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        # label.set_ha('right')
    
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    plt.xticks(np.arange(max_distance_plot*0.2, max_distance_plot*1.2, max_distance_plot*0.2))
    ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

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

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

CF_input.iloc[2:, 0:PADD_1].sum(axis=1).plot(ax=ax, color=b, linewidth=3)
CF_input.iloc[2:, PADD_1:PADD_2].sum(axis=1).plot(ax=ax, color=g, linewidth=3)
CF_input.iloc[2:, PADD_2:PADD_3].sum(axis=1).plot(ax=ax, color=r, linewidth=3)
CF_input.iloc[2:, PADD_3:PADD_4].sum(axis=1).plot(ax=ax, color=o, linewidth=3)
CF_input.iloc[2:, PADD_4:PADD_5].sum(axis=1).plot(ax=ax, color=y, linewidth=3)

ax.set_xlim([0, max_distance])
ax.set_ylim([0, 20000])

ax.set_xlabel(r'$\mathbf{Travel\ distance}$ [km]', fontname='Arial', fontsize=35)
ax.set_ylabel(r'$\mathbf{Cumulative\ WRRFs\ capacity}$ [MGD]', fontname='Arial', fontsize=35)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
plt.yticks(np.arange(0, 20000, 2000))
ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 1750, 250))
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

#%% CO2 abatement cost analysis

filterwarnings('ignore')

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# if just want to see the two plants in Urbana-Champaign:
# WRRF_input = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

print(len(WRRF_input))

WRRF_input['waste_cost'] = sum(WRRF_input[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000
WRRF_input['waste_GHG'] =  sum(WRRF_input[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

CWNS = []
facility_ID = []
CO2_reduction = []
sludge_CO2_reduction_ratio = []
WRRF_CO2_reduction_ratio = []
saves = []
USD_decarbonization = []
oil_BPD = []

# !!! run in different consoles to speed up: 0, 5000, 10000
for i in range(0, len(WRRF_input)):
# for i in range(0, 2):
    
    sys, barrel = create_geospatial_system(size=WRRF_input.iloc[i]['flow_2022_MGD'],
                                           sludge_transportation=0,
                                           sludge_distance=100,
                                           biocrude_distance=WRRF_input.iloc[i]['real_distance_km'],
                                           anaerobic_digestion=WRRF_input.iloc[i]['sludge_anaerobic_digestion'],
                                           aerobic_digestion=WRRF_input.iloc[i]['sludge_aerobic_digestion'],
                                           ww_2_dry_sludge_ratio=WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']/1000/365/WRRF_input.iloc[i]['flow_2022_MGD'],
                                           state=WRRF_input.iloc[i]['state'],
                                           elec_GHG=WRRF_input.iloc[i]['kg_CO2_kWh'])
    
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    WWTP = unit.WWTP
    raw_wastewater = stream.raw_wastewater
    tea = sys.TEA
    lca = sys.LCA
    
    # tonne
    sludge_tonne = raw_wastewater.F_vol*_m3perh_to_MGD*WWTP.ww_2_dry_sludge*(sys.operating_hours/24)*lca.lifetime
    
    try:
        # $/tonne
        sludge_cost = -tea.solve_price(raw_wastewater)*water_density*_MMgal_to_L/WWTP.ww_2_dry_sludge
    except RuntimeError:
        # TODO: double check this
        # due to high gross receipts tax (higher than net income, it might be impossible to break even)
        # TODO: there might be a problem: when solving price, the income tax in fixed, however, the raw_wastewater price could change net income, which can change income tax
        print('-------RUNTIME ERROR-------')
        
    try:
        saving_result = sludge_tonne*(WRRF_input.iloc[i]['waste_cost'] - sludge_cost)
    except NameError:
        saving_result = np.nan
        
    # kg CO2 eq/tonne
    sludge_CI = lca.get_total_impacts(exclude=(raw_wastewater,))['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)/lca.lifetime
    
    CO2_reduction_result = sludge_tonne*(WRRF_input.iloc[i]['waste_GHG'] - sludge_CI)
    
    try:
        sludge_CO2_reduction_ratio_result = CO2_reduction_result/WRRF_input.iloc[i]['biosolids_emission']/365/30
    # some WRRFs have incineration and the biosolids_emission is 0
    except FloatingPointError:
        sludge_CO2_reduction_ratio_result = np.nan
    
    WRRF_CO2_reduction_ratio_result = CO2_reduction_result/WRRF_input.iloc[i]['total_emission']/365/30
    
    # make sure CO2_reduction_result is positive if you want to calculate USD_per_tonne_CO2_reduction
    if CO2_reduction_result > 0:
        # save money: the results are negative; spend money: the results are positive
        USD_per_tonne_CO2_reduction = -saving_result/CO2_reduction_result*1000
    else:
        USD_per_tonne_CO2_reduction = np.nan
    
    CWNS.append(WRRF_input.iloc[i]['CWNS'])
    facility_ID.append(WRRF_input.iloc[i]['facility_ID'])
    CO2_reduction.append(CO2_reduction_result)
    sludge_CO2_reduction_ratio.append(sludge_CO2_reduction_ratio_result)
    WRRF_CO2_reduction_ratio.append(WRRF_CO2_reduction_ratio_result)
    saves.append(saving_result)
    USD_decarbonization.append(USD_per_tonne_CO2_reduction)
    oil_BPD.append(barrel)
    
    if USD_per_tonne_CO2_reduction < 0:
        print('HERE WE GO!')
    
    # check progress
    if i%50 == 0:
        print(i)
    
result = {'CWNS': CWNS,
          'facility_ID': facility_ID,
          'CO2_reduction': CO2_reduction,
          'sludge_CO2_reduction_ratio': sludge_CO2_reduction_ratio,
          'WRRF_CO2_reduction_ratio': WRRF_CO2_reduction_ratio,
          'saves': saves,
          'USD_decarbonization': USD_decarbonization,
          'oil_BPD': oil_BPD}
        
result = pd.DataFrame(result)

result.to_excel(folder + f'results/decarbonization_{date.today()}_{i}.xlsx')

#%% merge the results and the input

# !!! update the input file if necessary
input_data = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
input_data = input_data.dropna(subset='real_distance_km')

# !!! update these files if necessary
output_result_1 = pd.read_excel(folder + 'results/decarbonization_biocrude_baseline_data/decarbonization_2024-08-22_4999.xlsx')
output_result_2 = pd.read_excel(folder + 'results/decarbonization_biocrude_baseline_data/decarbonization_2024-08-22_9999.xlsx')
output_result_3 = pd.read_excel(folder + 'results/decarbonization_biocrude_baseline_data/decarbonization_2024-08-22_15858.xlsx')

output_result = pd.concat([output_result_1, output_result_2, output_result_3])

integrated_result = input_data.merge(output_result, how='left', on=['CWNS'])

integrated_result.to_excel(folder + f'results/integrated_decarbonization_result_{date.today()}.xlsx')

#%% read decarbonization data

# !!! update the file here if necessary
decarbonization_result = pd.read_excel(folder + 'results/integrated_decarbonization_result_2024-08-22.xlsx')
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 0]

#%% decarbonization map

# TODO: make 2 or 3 decarbonization maps; each for one type of WRRF (if 2: digestion & no digestion; if 3: AD & AeD & no digestion)

decarbonization_map = decarbonization_result.sort_values(by='total_emission', ascending=False).copy()
decarbonization_map = gpd.GeoDataFrame(decarbonization_map, crs='EPSG:4269',
                                       geometry=gpd.points_from_xy(x=decarbonization_map.longitude,
                                                                   y=decarbonization_map.latitude))
decarbonization_map = decarbonization_map.to_crs(crs='EPSG:3857')

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)

# tonne/day
decarbonization_map['CO2_reduction_tonne_per_day'] = decarbonization_map['CO2_reduction']/30/365/1000

WRRF_GHG_reduction_tonne_per_day = decarbonization_map['CO2_reduction_tonne_per_day']

more_than_20 = WRRF_GHG_reduction_tonne_per_day > 20
decarbonization_map[more_than_20].plot(ax=ax, color=dg, markersize=decarbonization_map.loc[more_than_20, WRRF_GHG_reduction_tonne_per_day.name]*100, edgecolor='k', linewidth=1.5, alpha=1)

between_2_and_20 = (WRRF_GHG_reduction_tonne_per_day > 2) & (WRRF_GHG_reduction_tonne_per_day <= 20)
decarbonization_map[between_2_and_20].plot(ax=ax, color=g, markersize=decarbonization_map.loc[between_2_and_20, WRRF_GHG_reduction_tonne_per_day.name]*100, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_2 = WRRF_GHG_reduction_tonne_per_day <= 2
decarbonization_map[less_than_2].plot(ax=ax, color=g, markersize=decarbonization_map.loc[less_than_2, WRRF_GHG_reduction_tonne_per_day.name]*100, alpha=0.1)

ax.set_aspect(1)

ax.set_axis_off()

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
    ax.set_xlim([-6, 26])
    ax.set_ylim([0, 500])
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    
    if position == 0:
        data = decarbonization_result[(decarbonization_result['sludge_anaerobic_digestion'] != 1) & (decarbonization_result['sludge_aerobic_digestion'] != 1)].copy()
        
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)
        
        plt.xticks(np.arange(0, 30, 10))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 30, 10))
        ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
        
        ax.set_ylabel(r'$\mathbf{Biocrude\ production}$ [BPD]', fontname='Arial', fontsize=35)
    
    # TODO: some WRRFs have both AD and AeD
    if position == 1:
        data = decarbonization_result[decarbonization_result['sludge_aerobic_digestion'] == 1].copy()

        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 30, 10))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 30, 10))
        ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')
        
        ax.set_xlabel(r'$\mathbf{Decarbonization\ potential}$', fontname='Arial', fontsize=35)
    
    if position == 2:
        data = decarbonization_result[decarbonization_result['sludge_anaerobic_digestion'] == 1].copy()
        
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 30, 10))
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 30, 10))
        ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
    
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
    ax.set_xlim([-15, 105])
    ax.set_ylim([0, 500])
    
    if position == 0:
        data = decarbonization_result[(decarbonization_result['sludge_anaerobic_digestion'] != 1) & (decarbonization_result['sludge_aerobic_digestion'] != 1)].copy()
        
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)
        
        plt.xticks(np.arange(0, 120, 30))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 120, 30))
        ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
        
        ax.set_ylabel(r'$\mathbf{Biocrude\ production}$ [BPD]', fontname='Arial', fontsize=35)
    
    # TODO: some WRRFs have both AD and AeD
    if position == 1:
        data = decarbonization_result[decarbonization_result['sludge_aerobic_digestion'] == 1].copy()

        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 120, 30))
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 120, 30))
        ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')
        
        ax.set_xlabel(r'$\mathbf{Decarbonization\ amount}$ [tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=35)

        mathtext.FontConstantsBase.sup1 = 0.35
    
    if position == 2:
        data = decarbonization_result[decarbonization_result['sludge_anaerobic_digestion'] == 1].copy()
    
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        plt.xticks(np.arange(0, 120, 30))
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        plt.xticks(np.arange(0, 120, 30))
        ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
    
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

no_digestion_WRRFs = decarbonization_result[(decarbonization_result['sludge_anaerobic_digestion'] != 1) & (decarbonization_result['sludge_aerobic_digestion'] != 1)].copy()
print(no_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(no_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(no_digestion_WRRFs.sort_values('oil_BPD', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])

# TODO: some WRRFs have both AD and AeD
aerobic_digestion_WRRFs = decarbonization_result[decarbonization_result['sludge_aerobic_digestion'] == 1].copy()
print(aerobic_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(aerobic_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(aerobic_digestion_WRRFs.sort_values('oil_BPD', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])

anaerobic_digestion_WRRFs = decarbonization_result[decarbonization_result['sludge_anaerobic_digestion'] == 1].copy()
print(anaerobic_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(anaerobic_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])
print(anaerobic_digestion_WRRFs.sort_values('oil_BPD', ascending=False).iloc[0,][['facility','city','flow_2022_MGD']])

#%% facility level decarbonizaiton amount vs sludge amount

decarbonization_vs_sludge = decarbonization_result[['total_sludge_amount_kg_per_year','CO2_reduction']].copy()

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.set_xlim((0, 300))
ax.set_ylim((0, 60))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Sludge}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Decarbonization}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 350, 50))
plt.yticks(np.arange(0, 70, 10))

ax_right = ax.twinx()
ax_right.set_ylim((0, 60))
ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 350, 50))
plt.yticks(np.arange(0, 70, 10))

ax_top = ax.twiny()
ax_top.set_xlim((0, 350))
# TODO: why length was 20/10 but 7.5 here
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 350, 50))
plt.yticks(np.arange(0, 70, 10))

ax.scatter(x=decarbonization_vs_sludge['total_sludge_amount_kg_per_year']/1000/365,
           y=decarbonization_vs_sludge['CO2_reduction']/30/365/1000,
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level decarbonizaiton amount vs distance

decarbonization_vs_sludge = decarbonization_result[['real_distance_km','CO2_reduction']].copy()

fig, ax = plt.subplots(figsize=(10, 10))

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
ax.set_ylim((0, 60))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Decarbonization}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 70, 10))

ax_right = ax.twinx()
ax_right.set_ylim((0, 60))
ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 70, 10))

ax_top = ax.twiny()
ax_top.set_xlim((0, 1200))
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 70, 10))

ax.scatter(x=decarbonization_vs_sludge['real_distance_km'],
           y=decarbonization_vs_sludge['CO2_reduction']/30/365/1000,
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level biocrude production vs sludge amount

decarbonization_vs_sludge = decarbonization_result[['total_sludge_amount_kg_per_year','oil_BPD']].copy()

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()

ax.set_xlim((0, 250))
ax.set_ylim((0, 500))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Sludge}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Biocrude}$ [BPD]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(0, 600, 100))

ax_right = ax.twinx()
ax_right.set_ylim((0, 500))
ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(0, 600, 100))

ax_top = ax.twiny()
ax_top.set_xlim((0, 250))
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(0, 600, 100))

ax.scatter(x=decarbonization_vs_sludge['total_sludge_amount_kg_per_year']/1000/365,
           y=decarbonization_vs_sludge['oil_BPD'],
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level biocrude production vs distance

decarbonization_vs_sludge = decarbonization_result[['real_distance_km','oil_BPD']].copy()

fig, ax = plt.subplots(figsize=(10, 10))

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
ax.set_ylim((0, 500))

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Biocrude}$ [BPD]', fontname='Arial', fontsize=45)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 600, 100))

ax_right = ax.twinx()
ax_right.set_ylim((0, 500))
ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 600, 100))

ax_top = ax.twiny()
ax_top.set_xlim((0, 1200))
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 1400, 200))
plt.yticks(np.arange(0, 600, 100))

ax.scatter(x=decarbonization_vs_sludge['real_distance_km'],
           y=decarbonization_vs_sludge['oil_BPD'],
           s=100,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% facility level uncertainty analysis

filterwarnings('ignore')

# !!! update the file here if necessary
decarbonization_result = pd.read_excel(folder + 'results/integrated_decarbonization_result_2024-08-22.xlsx')
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 0]

print(len(decarbonization_result))

decarbonization_result['waste_cost'] = sum(decarbonization_result[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/decarbonization_result['total_sludge_amount_kg_per_year']*1000
decarbonization_result['waste_GHG'] =  sum(decarbonization_result[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/decarbonization_result['total_sludge_amount_kg_per_year']*1000

geo_uncertainty_decarbonization = pd.DataFrame()
geo_uncertainty_biocrude = pd.DataFrame()

# !!! run in different consoles to speed up: 0, 80, 160, 240, 320, 400, 480
for i in range(0, len(decarbonization_result)):
# for i in range(0, 80):
    
    sys, barrel = create_geospatial_system(size=decarbonization_result.iloc[i]['flow_2022_MGD'],
                                           sludge_transportation=0,
                                           sludge_distance=100,
                                           biocrude_distance=decarbonization_result.iloc[i]['real_distance_km'],
                                           anaerobic_digestion=decarbonization_result.iloc[i]['sludge_anaerobic_digestion'],
                                           aerobic_digestion=decarbonization_result.iloc[i]['sludge_aerobic_digestion'],
                                           ww_2_dry_sludge_ratio=decarbonization_result.iloc[i]['total_sludge_amount_kg_per_year']/1000/365/decarbonization_result.iloc[i]['flow_2022_MGD'],
                                           state=decarbonization_result.iloc[i]['state'],
                                           elec_GHG=decarbonization_result.iloc[i]['kg_CO2_kWh'])
    
    # TODO: update if necessary; what if AD and AeD both exist or other co-existence scenarios?
    if decarbonization_result.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_values=[0.345, 0.436, 0.523, 'aerobic_digestion']
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
                                    sludge_protein=sludge_protein_values)
    
    kwargs = {'N':1000, 'rule':'L', 'seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    # kg CO2 eq/day
    geo_uncertainty_decarbonization[decarbonization_result.iloc[i]['CWNS']] = (decarbonization_result.iloc[i]['waste_GHG'] - results[('Geospatial','Sludge CI [kg CO2/tonne dry sludge]')])*decarbonization_result.iloc[i]['total_sludge_amount_kg_per_year']/1000/365
    # BPD
    geo_uncertainty_biocrude[decarbonization_result.iloc[i]['CWNS']] = results[('Geospatial','Biocrude production [BPD]')]
    
    # check progress
    if i%5 == 0:
        print(i)
    
geo_uncertainty_decarbonization.to_excel(folder + f'results/regional_decarbonization_uncertainty_{date.today()}_{i}.xlsx')
geo_uncertainty_biocrude.to_excel(folder + f'results/regional_biocrude_uncertainty_{date.today()}_{i}.xlsx')

#%% regional uncertainty (data preparation)

# TODO: continue from here

# decarbonization
# !!! update these files if necessary
regional_decarbonization_1 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_79.xlsx')
regional_decarbonization_2 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_159.xlsx')
regional_decarbonization_3 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_239.xlsx')
regional_decarbonization_4 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_319.xlsx')
regional_decarbonization_5 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_399.xlsx')
regional_decarbonization_6 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_479.xlsx')
regional_decarbonization_7 = pd.read_excel(folder + 'results/decarbonization_uncertainty_data/regional_decarbonization_uncertainty_2024-08-22_580.xlsx')

regional_decarbonization_1.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_2.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_3.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_4.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_5.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_6.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_7.drop('Unnamed: 0', axis=1, inplace=True)

regional_decarbonization = pd.concat([regional_decarbonization_1,
                                      regional_decarbonization_2,
                                      regional_decarbonization_3,
                                      regional_decarbonization_4,
                                      regional_decarbonization_5,
                                      regional_decarbonization_6,
                                      regional_decarbonization_7],
                                     axis=1)

regional_decarbonization.to_excel(folder + f'results/integrated_regional_decarbonization_uncertainty_{date.today()}.xlsx')

# biocrude
# !!! update these files if necessary
regional_biocrude_1 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_79.xlsx')
regional_biocrude_2 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_159.xlsx')
regional_biocrude_3 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_239.xlsx')
regional_biocrude_4 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_319.xlsx')
regional_biocrude_5 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_399.xlsx')
regional_biocrude_6 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_479.xlsx')
regional_biocrude_7 = pd.read_excel(folder + 'results/biocrude_uncertainty_data/regional_biocrude_uncertainty_2024-08-22_580.xlsx')

regional_biocrude_1.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_2.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_3.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_4.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_5.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_6.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_7.drop('Unnamed: 0', axis=1, inplace=True)

regional_biocrude = pd.concat([regional_biocrude_1,
                               regional_biocrude_2,
                               regional_biocrude_3,
                               regional_biocrude_4,
                               regional_biocrude_5,
                               regional_biocrude_6,
                               regional_biocrude_7],
                              axis=1)

regional_biocrude.to_excel(folder + f'results/integrated_regional_biocrude_uncertainty_{date.today()}.xlsx')

#%% regional uncertainty (build and run model)

# import PADD information
# !!! update the file here if necessary
WRRF_PADD = pd.read_excel(folder + 'results/integrated_decarbonization_result_2024-08-22.xlsx')
WRRF_PADD = WRRF_PADD[WRRF_PADD['USD_decarbonization'].notna()]
WRRF_PADD = WRRF_PADD[WRRF_PADD['USD_decarbonization'] <= 0]

WRRF_PADD.loc[WRRF_PADD['state'].isin(['CT','DC','DE','FL','GA','MA','MD','ME','NC','NH','NJ','NY','PA','RI','SC','VA','VT','WV']), 'WRRF_PADD'] = 1
WRRF_PADD.loc[WRRF_PADD['state'].isin(['IA','IL','IN','KS','KY','MI','MN','MO','ND','NE','OH','OK','SD','TN','WI']), 'WRRF_PADD'] = 2
WRRF_PADD.loc[WRRF_PADD['state'].isin(['AL','AR','LA','MS','NM','TX']), 'WRRF_PADD'] = 3
WRRF_PADD.loc[WRRF_PADD['state'].isin(['CO','ID','MT','UT','WY']), 'WRRF_PADD'] = 4
WRRF_PADD.loc[WRRF_PADD['state'].isin(['AZ','CA','NV','OR','WA']), 'WRRF_PADD'] = 5

assert WRRF_PADD.WRRF_PADD.isna().sum() == 0

WRRF_PADD = WRRF_PADD[['CWNS','WRRF_PADD']]

# decarbonization
# !!! update the file here if necessary
regional_decarbonization_uncertainty = pd.read_excel(folder + 'results/integrated_regional_decarbonization_uncertainty_2024-08-22.xlsx')

regional_decarbonization_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_uncertainty = regional_decarbonization_uncertainty.transpose()
regional_decarbonization_uncertainty.reset_index(inplace=True, names='CWNS')
regional_decarbonization_uncertainty = regional_decarbonization_uncertainty.merge(WRRF_PADD, how='left', on='CWNS')

regional_decarbonization_uncertainty_PADD_1 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 1].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_2 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 2].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_3 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 3].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_4 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 4].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_5 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 5].loc[:, 0:999]

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
    # Monte Carlo simulation and Latin hypercube sampling
    for i in range(0, 1000):
        # check progress
        if i%100 == 0:
            print(i)
        regional_total_decarbonization = 0
        for j in range(0, len(region)):
            regional_total_decarbonization += random.uniform(region.iloc[j,:].quantile(i/1000), region.iloc[j,:].quantile((i+1)/1000))
        result.append(regional_total_decarbonization)
            
integrated_regional_decarbonization_result = pd.DataFrame({'PADD_1': PADD_1_decarbonization_uncertainty,
                                                           'PADD_2': PADD_2_decarbonization_uncertainty,
                                                           'PADD_3': PADD_3_decarbonization_uncertainty,
                                                           'PADD_4': PADD_4_decarbonization_uncertainty,
                                                           'PADD_5': PADD_5_decarbonization_uncertainty})
        
integrated_regional_decarbonization_result.to_excel(folder + f'results/integrated_regional_total_decarbonization_uncertainty_{date.today()}.xlsx')

# biocrude
# !!! update the file here if necessary
regional_biocrude_uncertainty = pd.read_excel(folder + 'results/integrated_regional_biocrude_uncertainty_2024-08-22.xlsx')

regional_biocrude_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)
regional_biocrude_uncertainty = regional_biocrude_uncertainty.transpose()
regional_biocrude_uncertainty.reset_index(inplace=True, names='CWNS')
regional_biocrude_uncertainty = regional_biocrude_uncertainty.merge(WRRF_PADD, how='left', on='CWNS')

regional_biocrude_uncertainty_PADD_1 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 1].loc[:, 0:999]
regional_biocrude_uncertainty_PADD_2 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 2].loc[:, 0:999]
regional_biocrude_uncertainty_PADD_3 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 3].loc[:, 0:999]
regional_biocrude_uncertainty_PADD_4 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 4].loc[:, 0:999]
regional_biocrude_uncertainty_PADD_5 = regional_biocrude_uncertainty[regional_biocrude_uncertainty['WRRF_PADD'] == 5].loc[:, 0:999]

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
    # Monte Carlo simulation and Latin hypercube sampling
    for i in range(0,1000):
        # check progress
        if i%100 == 0:
            print(i)
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

# TODO: update the code/figure during execution, if necessary

# !!! update these files if necessary
decarbonization = pd.read_excel(folder + 'results/integrated_regional_total_decarbonization_uncertainty_2024-08-22.xlsx')
# tonne CO2 eq/day
decarbonization = decarbonization/1000
# BPD
biocrude = pd.read_excel(folder + 'results/integrated_regional_total_biocrude_uncertainty_2024-08-22.xlsx')


fig, ax = plt.subplots(figsize = (10, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax.set_xlim([0, 1000])
ax.set_ylim([0, 12000])

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1200, 200))
plt.yticks(np.arange(0, 14000, 2000))

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1200, 200))
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
# TODO: add point for PADD 4 or nor
add_point('PADD_4', do)
add_point('PADD_5', dy)

#%% national uncertainty

# decarbonization
# !!! update the file here if necessary
national_decarbonization_uncertainty = pd.read_excel(folder + 'results/integrated_regional_decarbonization_uncertainty_2024-08-22.xlsx')

national_decarbonization_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)
national_decarbonization_uncertainty = national_decarbonization_uncertainty.transpose()
national_decarbonization_uncertainty.reset_index(inplace=True, names='CWNS')
national_decarbonization_uncertainty = national_decarbonization_uncertainty.loc[:, 0:999]

national_decarbonization_uncertainty_result = []

# Monte Carlo simulation and Latin hypercube sampling
for i in range(0,1000):
    # check progress
    if i%100 == 0:
        print(i)
    national_total_decarbonization = 0
    for j in range(0, len(national_decarbonization_uncertainty)):
        national_total_decarbonization += random.uniform(national_decarbonization_uncertainty.iloc[j,:].quantile(i/1000), national_decarbonization_uncertainty.iloc[j,:].quantile((i+1)/1000))
    national_decarbonization_uncertainty_result.append(national_total_decarbonization)

# biocrude
# !!! update the file here if necessary
national_biocrude_uncertainty = pd.read_excel(folder + 'results/integrated_regional_biocrude_uncertainty_2024-08-22.xlsx')

national_biocrude_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)
national_biocrude_uncertainty = national_biocrude_uncertainty.transpose()
national_biocrude_uncertainty.reset_index(inplace=True, names='CWNS')
national_biocrude_uncertainty = national_biocrude_uncertainty.loc[:, 0:999]

national_biocrude_uncertainty_result = []

# Monte Carlo simulation and Latin hypercube sampling
for i in range(0,1000):
    # check progress
    if i%100 == 0:
        print(i)
    national_total_biocrude = 0
    for j in range(0, len(national_biocrude_uncertainty)):
        national_total_biocrude += random.uniform(national_biocrude_uncertainty.iloc[j,:].quantile(i/1000), national_biocrude_uncertainty.iloc[j,:].quantile((i+1)/1000))
    national_biocrude_uncertainty_result.append(national_total_biocrude)
            
national_uncertainty_result = pd.DataFrame({'decarbonization':national_decarbonization_uncertainty_result,
                                            'biocrude':national_biocrude_uncertainty_result})

national_uncertainty_result.to_excel(folder + f'results/integrated_national_uncertainty_{date.today()}.xlsx')

#%% national chord diagram (biocrude between states)

# import only if needed
from d3blocks import D3Blocks

biocrude_transportation = pd.read_excel(folder + 'results/integrated_decarbonization_result_2024-08-22.xlsx')
biocrude_transportation = biocrude_transportation[biocrude_transportation['USD_decarbonization'].notna()]
biocrude_transportation = biocrude_transportation[biocrude_transportation['USD_decarbonization'] <= 0]
biocrude_transportation = biocrude_transportation.merge(elec_price[['state','name']], how='left', on='state')
biocrude_transportation['source'] = biocrude_transportation['name']
biocrude_transportation['target'] = biocrude_transportation['State']
biocrude_transportation['weight'] = biocrude_transportation['oil_BPD']
biocrude_transportation = biocrude_transportation[['source','target','weight']]

d3 = D3Blocks()

ordering = sorted(state_PADD.items(), key=lambda x: x[1])   
ordering = [i[0] for i in ordering]

d3.chord(biocrude_transportation,
         ordering=ordering,
         color='source',
         fontsize=17.5,
         arrowhead=10,
         filepath=folder+'results/interstate_biocrude_transportation.html')

d3.node_properties['color'] = d3.node_properties['label'].apply(lambda x: PADD_color[state_PADD[x]])
d3.edge_properties['color'] = d3.edge_properties['source'].apply(lambda x: PADD_color[state_PADD[x]])

d3.show()

#%% sludge transportation (Urbana-Champaign, two separated WRRFs)

filterwarnings('ignore')

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

WRRF_input['waste_cost'] = sum(WRRF_input[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000
WRRF_input['waste_GHG'] =  sum(WRRF_input[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

Urbana_Champaign = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

CU_uncertainty_saving = pd.DataFrame()
CU_uncertainty_decarbonization = pd.DataFrame()
CU_uncertainty_biocrude = pd.DataFrame()

for i in range(0, 2):
    sys, barrel = create_geospatial_system(size=Urbana_Champaign.iloc[i]['flow_2022_MGD'],
                                           sludge_transportation=0,
                                           sludge_distance=100,
                                           biocrude_distance=Urbana_Champaign.iloc[i]['real_distance_km'],
                                           anaerobic_digestion=Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'],
                                           aerobic_digestion=Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'],
                                           ww_2_dry_sludge_ratio=Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365/Urbana_Champaign.iloc[i]['flow_2022_MGD'],
                                           state=Urbana_Champaign.iloc[i]['state'],
                                           elec_GHG=Urbana_Champaign.iloc[i]['kg_CO2_kWh'])
    
    # TODO: update if necessary; what if AD and AeD both exist or other co-existence scenarios?
    if Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_values=[0.345, 0.436, 0.523, 'aerobic_digestion']
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
                                    sludge_protein=sludge_protein_values)
    
    kwargs = {'N':1000, 'rule':'L', 'seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    # $/day
    CU_uncertainty_saving[Urbana_Champaign.iloc[i]['CWNS']] = (Urbana_Champaign.iloc[i]['waste_cost'] - results[('Geospatial','Sludge management price [$/tonne dry sludge]')])*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365
    # kg CO2 eq/day
    CU_uncertainty_decarbonization[Urbana_Champaign.iloc[i]['CWNS']] = (Urbana_Champaign.iloc[i]['waste_GHG'] - results[('Geospatial','Sludge CI [kg CO2/tonne dry sludge]')])*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365
    # BPD
    CU_uncertainty_biocrude[Urbana_Champaign.iloc[i]['CWNS']] = results[('Geospatial','Biocrude production [BPD]')]

#%% sludge transportation (Urbana-Champaign, after sludge transportation)

CU_combined_uncertainty_saving = pd.DataFrame()
CU_combined_uncertainty_decarbonization = pd.DataFrame()
CU_combined_uncertainty_biocrude = pd.DataFrame()

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
        sludge_ash_value=0.436
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
# set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion']
sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

average_ww_2_dry_sludge_ratio = (Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] +\
                                 Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']) /\
                                 1000/365/(Urbana_Champaign.iloc[0]['flow_2022_MGD']+Urbana_Champaign.iloc[1]['flow_2022_MGD'])

assert Urbana_Champaign.iloc[0]['state'] == Urbana_Champaign.iloc[1]['state']
assert Urbana_Champaign.iloc[0]['kg_CO2_kWh'] == Urbana_Champaign.iloc[1]['kg_CO2_kWh']

# sludge_distance from Google Maps
sys, barrel = create_geospatial_system(size=total_size,
                                       sludge_transportation=1,
                                       sludge_distance=16.2*Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']/(Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year']+Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']),
                                       biocrude_distance=distance_to_refinery,
                                       average_sludge_dw_ash=sludge_ash_values[1],
                                       average_sludge_afdw_lipid=sludge_lipid_values[1],
                                       average_sludge_afdw_protein=sludge_protein_values[1],
                                       anaerobic_digestion=None,
                                       aerobic_digestion=None,
                                       ww_2_dry_sludge_ratio=average_ww_2_dry_sludge_ratio,
                                       state=Urbana_Champaign.iloc[0]['state'],
                                       elec_GHG=Urbana_Champaign.iloc[0]['kg_CO2_kWh'])

model = create_geospatial_model(system=sys,
                                sludge_ash=sludge_ash_values,
                                sludge_lipid=sludge_lipid_values,
                                sludge_protein=sludge_protein_values)

kwargs = {'N':1000, 'rule':'L', 'seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()

idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
results = model.table.iloc[:, idx:]

# $/day
CU_combined_uncertainty_saving['combined'] = (average_cost - results[('Geospatial','Sludge management price [$/tonne dry sludge]')])*total_sludge/1000/365
# kg CO2 eq/day
CU_combined_uncertainty_decarbonization['combined'] = (average_GHG - results[('Geospatial','Sludge CI [kg CO2/tonne dry sludge]')])*total_sludge/1000/365
# BPD
CU_combined_uncertainty_biocrude['combined'] = results[('Geospatial','Biocrude production [BPD]')]

CU_saving_results = pd.concat([CU_uncertainty_saving, CU_combined_uncertainty_saving], axis=1)
CU_decarbonization_results = pd.concat([CU_uncertainty_decarbonization, CU_combined_uncertainty_decarbonization], axis=1)
CU_biocrude_results = pd.concat([CU_uncertainty_biocrude, CU_combined_uncertainty_biocrude], axis=1)

CU_saving_results.to_excel(folder + f'results/Urbana-Champaign/saving_{date.today()}.xlsx')
CU_decarbonization_results.to_excel(folder + f'results/Urbana-Champaign/decarbonization_{date.today()}.xlsx')
CU_biocrude_results.to_excel(folder + f'results/Urbana-Champaign/biocrude_{date.today()}.xlsx')

#%% Urbana-Champaign visualization

CU_saving_results = pd.read_excel(folder + 'results/Urbana-Champaign/saving_2024-08-23.xlsx')
CU_saving_results.drop('Unnamed: 0', axis=1, inplace=True)

fig = plt.figure(figsize=(10, 10))

gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)

def add_region(position, xlabel, color):
    ax = fig.add_subplot(gs[0, position])
    
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 38
    plt.rcParams['ytick.labelsize'] = 38
    
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    ax = plt.gca()
    ax.set_ylim([-3, 1])
    ax.set_xlabel(xlabel, fontname='Arial', fontsize=38, labelpad=15)
    
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)
        ax.set_ylabel(r'$\mathbf{Saving}$ [k-\$·day${^{-1}}$]', fontname='Arial', fontsize=45)
        
        mathtext.FontConstantsBase.sup1 = 0.35
    
    elif position == 2:
        ax.tick_params(direction='inout', length=0, labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
    
    else:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
    
    bp = ax.boxplot(CU_saving_results.iloc[:,position]/1000, showfliers=False, widths=0.7, patch_artist=True)
    
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
        
    ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0, 0], color='k', linestyle='--', linewidth=3)

add_region(0, 'Northeast', o)
add_region(1, 'Southwest', o)
add_region(2, 'Combined', b)

#%% sludge transportation (satellite, SA) (data preparation)

# !!! update the file here if necessary
decarbonization_result = pd.read_excel(folder + 'results/integrated_decarbonization_result_2024-08-22.xlsx')
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 0]

Minnesota = decarbonization_result[decarbonization_result['state']=='MN'].copy()
Minnesota.sort_values(by='total_sludge_amount_kg_per_year', ascending=False, inplace=True)

center_WRRF = Minnesota.iloc[0,:]

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# TODO: replace manual check with assertations, if necessary
# manually confirmed: only WRRFs in MN, WI, and IA can be within 200 km (linear distance) from the center WRRF
# get all WRRFs in MN, WI, IA
MN_WI_IA = WRRF_input[WRRF_input['state'].isin(['MN','WI','IA'])].copy()

# select WRRFs that are within 200 km (linear distance) from the center WRRF
MN_WI_IA_distance = []
for i in range(len(MN_WI_IA)):
    MN_WI_IA_distance.append(geopy.distance.geodesic((MN_WI_IA['latitude'].iloc[i], MN_WI_IA['longitude'].iloc[i]), (center_WRRF.latitude, center_WRRF.longitude)).km)

MN_WI_IA['linear_distance_to_center_WRRF_km'] = MN_WI_IA_distance

MN_WI_IA = MN_WI_IA[MN_WI_IA['linear_distance_to_center_WRRF_km'] <= 200]

# remove WRRFs that are already decarbonized without sludge transportation
MN_WI_IA_decarbonized = decarbonization_result[decarbonization_result['state'].isin(['MN','WI','IA'])]

MN_WI_IA_transportation = MN_WI_IA[~MN_WI_IA['CWNS'].isin(list(MN_WI_IA_decarbonized['CWNS']))].copy()

# !!! get a google API key !!! do not upload to GitHub
gmaps = googlemaps.Client(key='XXX')

for i in MN_WI_IA_transportation.index:
    try:
        print(i)
        MN_WI_IA_transportation.loc[i, 'real_distance_to_center_WRRF_km'] = gmaps.distance_matrix((MN_WI_IA_transportation.loc[i, 'latitude'], MN_WI_IA_transportation.loc[i, 'longitude']),
                                                                                                  (center_WRRF.latitude, center_WRRF.longitude),
                                                                                                  mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000
    except KeyError:
        print('--------------------------------')
        MN_WI_IA_transportation.loc[i, 'real_distance_to_center_WRRF_km'] = np.nan

MN_WI_IA_transportation.to_excel(folder + f'results/Minnesota/satellite_WRRFs_{date.today()}.xlsx')

#%% sludge transportation (satellite, SA) (center WRRF uncertainty)

filterwarnings('ignore')

# !!! update the file here if necessary
decarbonization_result = pd.read_excel(folder + 'results/integrated_decarbonization_result_2024-08-22.xlsx')
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'].notna()]
decarbonization_result = decarbonization_result[decarbonization_result['USD_decarbonization'] <= 0]

Minnesota = decarbonization_result[decarbonization_result['state']=='MN'].copy()
Minnesota.sort_values(by='total_sludge_amount_kg_per_year', ascending=False, inplace=True)

center_WRRF = Minnesota.iloc[0,:]

sys, barrel = create_geospatial_system(size=center_WRRF['flow_2022_MGD'],
                                       sludge_transportation=0,
                                       sludge_distance=100,
                                       biocrude_distance=center_WRRF['real_distance_km'],
                                       anaerobic_digestion=center_WRRF['sludge_anaerobic_digestion'],
                                       aerobic_digestion=center_WRRF['sludge_aerobic_digestion'],
                                       ww_2_dry_sludge_ratio=center_WRRF['total_sludge_amount_kg_per_year']/1000/365/center_WRRF['flow_2022_MGD'],
                                       state=center_WRRF['state'],
                                       elec_GHG=center_WRRF['kg_CO2_kWh'])

# TODO: update if necessary; what if AD and AeD both exist or other co-existence scenarios?
if center_WRRF['sludge_aerobic_digestion'] == 1:
    sludge_ash_values=[0.345, 0.436, 0.523, 'aerobic_digestion']
    sludge_lipid_values=[0.154, 0.193, 0.232]
    sludge_protein_values=[0.408, 0.510, 0.612]
elif center_WRRF['sludge_anaerobic_digestion'] == 1:
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
                                sludge_protein=sludge_protein_values)

kwargs = {'N':1000, 'rule':'L', 'seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()

idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
results = model.table.iloc[:, idx:]

# $/day
saving_center = (center_WRRF['waste_cost'] - results[('Geospatial','Sludge management price [$/tonne dry sludge]')])*center_WRRF['total_sludge_amount_kg_per_year']/1000/365
# kg CO2 eq/day
decarbonization_center = (center_WRRF['waste_GHG'] - results[('Geospatial','Sludge CI [kg CO2/tonne dry sludge]')])*center_WRRF['total_sludge_amount_kg_per_year']/1000/365
# BPD
biocrude_center = results[('Geospatial','Biocrude production [BPD]')]

#%% sludge transportation (satellite, SA) (all WRRFs uncertainty)

filterwarnings('ignore')

# !!! update the file here if necessary
MN_WI_IA_transportation = pd.read_excel(folder + 'results/Minnesota/satellite_WRRFs_2024-08-24.xlsx')

print(len(MN_WI_IA_transportation))

MN_WI_IA_transportation['waste_cost'] = sum(MN_WI_IA_transportation[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/MN_WI_IA_transportation['total_sludge_amount_kg_per_year']*1000
MN_WI_IA_transportation['waste_GHG'] =  sum(MN_WI_IA_transportation[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/MN_WI_IA_transportation['total_sludge_amount_kg_per_year']*1000

# the units here do not matter
total_cost_except_center = (MN_WI_IA_transportation['waste_cost']*MN_WI_IA_transportation['total_sludge_amount_kg_per_year']).sum(axis=0)
total_GHG_except_center = (MN_WI_IA_transportation['waste_GHG']*MN_WI_IA_transportation['total_sludge_amount_kg_per_year']).sum(axis=0)
total_sludge_except_center = MN_WI_IA_transportation['total_sludge_amount_kg_per_year'].sum(axis=0)

total_cost = total_cost_except_center + center_WRRF['waste_cost']*center_WRRF['total_sludge_amount_kg_per_year']
total_GHG = total_GHG_except_center + center_WRRF['waste_GHG']*center_WRRF['total_sludge_amount_kg_per_year']
total_sludge = total_sludge_except_center + center_WRRF['total_sludge_amount_kg_per_year']

average_cost = total_cost/total_sludge
average_GHG = total_GHG/total_sludge

total_size = MN_WI_IA_transportation['flow_2022_MGD'].sum(axis=0) + center_WRRF['flow_2022_MGD']

distance_to_refinery = center_WRRF['real_distance_km']

sludge_transporataion_distance = (MN_WI_IA_transportation['real_distance_to_center_WRRF_km']*MN_WI_IA_transportation['total_sludge_amount_kg_per_year']).sum(axis=0)/total_sludge

total_ash_except_center = 0
total_lipid_except_center = 0
total_protein_except_center = 0

for i in range(len(MN_WI_IA_transportation)):
    # TODO: update if necessary; what if AD and AeD both exist or other co-existence scenarios?
    if MN_WI_IA_transportation.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_value=0.436
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

# the center WRRF does not have anaerobic or aerobic digestion
assert center_WRRF.sludge_aerobic_digestion == center_WRRF.sludge_anaerobic_digestion == 0

total_ash = total_ash_except_center + center_WRRF['total_sludge_amount_kg_per_year']*0.231
total_lipid = total_lipid_except_center + center_WRRF['total_sludge_amount_kg_per_year']*(1-0.231)*0.206
total_protein = total_protein_except_center + center_WRRF['total_sludge_amount_kg_per_year']*(1-0.231)*0.456

average_ash_dw = total_ash/total_sludge
average_lipid_afdw = total_lipid/(total_sludge-total_ash)
average_protein_afdw = total_protein/(total_sludge-total_ash)

# !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
# set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion']
sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

average_ww_2_dry_sludge_ratio = total_sludge/1000/365/total_size

sys, barrel = create_geospatial_system(size=total_size,
                                       sludge_transportation=1,
                                       sludge_distance=sludge_transporataion_distance,
                                       biocrude_distance=distance_to_refinery,
                                       average_sludge_dw_ash=sludge_ash_values[1],
                                       average_sludge_afdw_lipid=sludge_lipid_values[1],
                                       average_sludge_afdw_protein=sludge_protein_values[1],
                                       anaerobic_digestion=None,
                                       aerobic_digestion=None,
                                       ww_2_dry_sludge_ratio=average_ww_2_dry_sludge_ratio,
                                       state=center_WRRF['state'],
                                       elec_GHG=center_WRRF['kg_CO2_kWh'])

model = create_geospatial_model(system=sys,
                                sludge_ash=sludge_ash_values,
                                sludge_lipid=sludge_lipid_values,
                                sludge_protein=sludge_protein_values)

kwargs = {'N':1000, 'rule':'L', 'seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()

idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
all_WRRFs_results = model.table.iloc[:, idx:]

# $/day
saving_all = (average_cost - all_WRRFs_results[('Geospatial','Sludge management price [$/tonne dry sludge]')])*total_sludge/1000/365
# kg CO2 eq/day
decarbonization_all = (average_GHG - all_WRRFs_results[('Geospatial','Sludge CI [kg CO2/tonne dry sludge]')])*total_sludge/1000/365
# BPD
biocrude_all = all_WRRFs_results[('Geospatial','Biocrude production [BPD]')]

SA_results = pd.DataFrame({'center_saving': saving_center,
                           'all_saving': saving_all,
                           'center_decarbonization': decarbonization_center,
                           'all_decarbonization': decarbonization_all,
                           'center_biocrude': biocrude_center,
                           'all_biocrude': biocrude_all})

SA_results.to_excel(folder + f'results/Minnesota/center_vs_all_{date.today()}.xlsx')

#%% sludge transportation (satellite, SA) (visualization)

# !!! update the file here if necessary
SA_results = pd.read_excel(folder + 'results/Minnesota/center_vs_all_2024-08-24.xlsx')
SA_results.drop('Unnamed: 0', axis=1, inplace=True)
# convert decarbonization from kg CO2 eq/day to tonne CO2 eq/day
SA_results['center_decarbonization'] /= 1000
SA_results['all_decarbonization'] /= 1000

fig, ax = plt.subplots(figsize = (10, 10))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax.set_xlim([-20, 50])
ax.set_ylim([200, 700])

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, pad=6)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(-20, 60, 10))
plt.yticks(np.arange(200, 800, 100))

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(-20, 60, 10))
plt.yticks(np.arange(200, 800, 100))

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

ax.plot([0, 0], [200, 700], color='k', linestyle='--', linewidth=3)

#%% sludge transportation (heat map, HM)

filterwarnings('ignore')

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2024-08-20.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

print(len(WRRF_input))

# TODO: continue from here

# TODO: check the following notes
# TODO: is this the average sludge/biosolids compositions? is there a way to also considered lime stabilization?
# use data for sludge + biosolids (as we did in the HTL model paper)
# set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
HM_sludge_ash_values = [0.174,0.257,0.414, 'digestion']
HM_sludge_lipid_values = [0.08,0.204,0.308]
HM_sludge_protein_values = [0.38,0.463,0.51]

# TODO: use nation average?
HM_elec_GHG = elec['GHG (10-year median)'].mean()

# TODO: plot saving and decarbonization amount?
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
        
        # TODO: check parameters
        # TODO: how do we get 207? can we replace this with a formula, if possible?
        # 207 km is the median travel distance from WRRFs to oil refineries
        sys, barrel = create_geospatial_system(size=size,
                                               sludge_transportation=1,
                                               sludge_distance=sludge_distance,
                                               biocrude_distance=207,
                                               average_sludge_dw_ash=HM_sludge_ash_values[1],
                                               average_sludge_afdw_lipid=HM_sludge_lipid_values[1],
                                               average_sludge_afdw_protein=HM_sludge_protein_values[1],
                                               anaerobic_digestion=None,
                                               aerobic_digestion=None,
                                               ww_2_dry_sludge_ratio=1,
                                               state='US',
                                               elec_GHG=HM_elec_GHG)
        
        # TODO: check parameters
        # the distribution of biocrude transportation cost was assumed to be proportional to the distribution of biocrude price (triangular: 0.2384/0.3847/0.6703, or 4.21/6.80/11.9)
        model = create_geospatial_model(system=sys,
                                        sludge_ash=HM_sludge_ash_values,
                                        sludge_lipid=HM_sludge_lipid_values,
                                        sludge_protein=HM_sludge_protein_values)
        
        kwargs = {'N':1000, 'rule':'L', 'seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        # TODO: confirm 'Geospatial' is a tab
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

# TODO: update the code/figure during execution, if necessary

# TODO: update the file
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
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

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

# TODO: update the code/figure during execution, if necessary

# TODO: update the file
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
ax_top.tick_params(direction='in', length=15, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

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