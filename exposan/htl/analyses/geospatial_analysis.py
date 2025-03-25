#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Sun Jun 11 08:12:41 2023

@author: jiananfeng

Note the word 'sludge' in this file refers to either sludge or biosolids.
'''

# confirmed: both 'markersize' (in gpd.plot()) and 's' (in plt.scatter()) are proportional to area

#%% initialization

import geopy.distance, googlemaps
import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt, matplotlib.colors as colors, scipy.stats as stats, qsdsan as qs, networkx as nx
from math import floor
from matplotlib.mathtext import _mathtext as mathtext
from matplotlib.patches import Rectangle
from colorpalette import Color
from exposan.htl import create_geospatial_system, create_geospatial_model
from qsdsan.utils import auom, palettes
from datetime import date
from warnings import filterwarnings
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.cluster import AgglomerativeClustering
from scipy.linalg import cholesky
from scipy.spatial import KDTree
from numba import njit

mathtext.FontConstantsBase.sup1 = 0.35

# TODO: update file paths later
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
_oil_barrel_to_L = auom('oil_barrel').conversion_factor('L')

WRRF = pd.read_excel(folder + 'HTL_geospatial_input_2025-02-10.xlsx')

assert WRRF.duplicated(subset='CWNS_NUM').sum() == 0

WRRF['total_sludge_amount_kg_per_year'] = WRRF[['landfill','land_application',
                                                'incineration']].sum(axis=1)

WRRF['biosolids_emission'] = WRRF[['LF_CH4','LA_N2O']].sum(axis=1)

# all emission data in kg CO2 eq/day
WRRF['total_emission'] = WRRF[['LF_CH4','LA_N2O','electricity_emission',
                               'onsite_NG_emission','upstream_NG_emission',
                               'CH4_emission','N2O_emission','CO2_emission']].sum(axis=1)

treatment_trains = np.array(['LAGOON_AER','LAGOON_ANAER','LAGOON_FAC',
                             'LAGOON_UNCATEGORIZED','C1','C1E','C2','C3',
                             'C5','C6','B1','B1E','B2','B3','B4','B5',
                             'B6','D1','D1E','D2','D3','D5','D6','E2',
                             'E2P','F1','F1E','I1','I1E','I2','I3','I5',
                             'I6','G1','G1E','G2','G3','G5','G6','H1',
                             'H1E','N1','N1E','N2','O1','O1E','O2','O3',
                             'O5','O6'], dtype=object)

TT_indentifier = WRRF[treatment_trains].apply(lambda x: x > 0)
WRRF['treatment_train'] = TT_indentifier.apply(lambda x: list(treatment_trains[x.values]), axis=1)

# !!! 'total_emission' here is based on the IEDO grid CI (2020 standard scenario), but for HTL-based systems, use 2021 standard scenario, which may have minimal impact on WRRF life cycle GHG reduction calculation
WRRF = WRRF[['FACILITY','CITY','STATE','CWNS_NUM','LATITUDE','LONGITUDE',
             'FLOW_2022_MGD_FINAL','balancing_area','treatment_train',
             'landfill','land_application','incineration',
             'total_sludge_amount_kg_per_year','CH4_emission','N2O_emission',
             'CO2_emission','electricity_emission','onsite_NG_emission',
             'upstream_NG_emission','biosolids_emission','total_emission']]

# assume LAGOON_FAC and LAGOON_UNCATEGORIZED has AeD (resulting in higher ash content than AD) to be conservative
TT_w_AeD = ['B2','C2','D2','E2','E2P','G2','I2','N2','O2','LAGOON_AER','LAGOON_FAC','LAGOON_UNCATEGORIZED']
TT_w_AD = ['B1','B1E','B4','C1','C1E','D1','D1E','F1','F1E','G1','G1E','H1','H1E','I1','I1E','N1','N1E','O1','O1E','LAGOON_ANAER']

# to be conservative (AeD has higher ash content than AD):
# if a WRRF has AeD (regardless of AD), assume all sludge composition follows the AeD sludge composition,
# if no AeD but has AD, then use AD composition
# if none, then use sludge composition
WRRF.loc[WRRF['treatment_train'].apply(lambda x: len([i for i in TT_w_AeD if i in x]) > 0), 'sludge_aerobic_digestion'] = 1
WRRF.loc[WRRF['treatment_train'].apply(lambda x: len([i for i in TT_w_AD if i in x]) > 0), 'sludge_anaerobic_digestion'] = 1
WRRF.fillna({'sludge_aerobic_digestion': 0,
             'sludge_anaerobic_digestion': 0},
            inplace=True)

elec_CI = pd.read_csv(folder + 'StdScen21_MidCase_annual_balancingArea.csv')
elec_CI = elec_CI[['r','kg_CO2e_kWh']]
elec_CI['r'] = elec_CI['r'].apply(lambda x: x[1:])
elec_CI['r'] = pd.to_numeric(elec_CI['r'])

WRRF = WRRF.merge(elec_CI, how='left', left_on='balancing_area', right_on='r')

assert WRRF.kg_CO2e_kWh.isna().sum() == 0

WRRF = WRRF.sort_values(by='FLOW_2022_MGD_FINAL', ascending=False)

WRRF = gpd.GeoDataFrame(WRRF, crs='EPSG:4269',
                        geometry=gpd.points_from_xy(x=WRRF.LONGITUDE,
                                                    y=WRRF.LATITUDE))

refinery = pd.read_csv(folder + 'petroleum_refineries_EIA_2024-06-06.csv')
refinery['capacity'] = refinery[['Atmos. Crude Dist','Vacuum Dist','Catalytic Cracking',
                                 'Hydro Cracking','Thermal Cracking, Visbreaking',
                                 'Catalytic Recorming','Alkylates, Isomerization',
                                 'Desulfurization','Fluid and Delayed Coking',
                                 'Asphalt and Road Oil']].sum(axis=1)
# remove oil refineries with a 0 MBPD production
refinery = refinery[refinery['capacity'] > 0]
refinery = gpd.GeoDataFrame(refinery, crs='EPSG:4269',
                            geometry=gpd.points_from_xy(x=refinery.Longitude,
                                                        y=refinery.Latitude))

# confirmed from the wesbite (next line): there is no change in US map since June 1, 1995. So, using 2018 US map data is OK.
# https://en.wikipedia.org/wiki/Territorial_evolution_of_the_United_States#1946%E2%80%93present_(Decolonization)
US = gpd.read_file(folder + 'US/cb_2018_us_state_500k.shp')
US = US[['NAME','geometry']]

for excluded in ('Alaska',
                 'American Samoa',
                 'Commonwealth of the Northern Mariana Islands',
                 'Guam',
                 'Hawaii',
                 'Puerto Rico',
                 'United States Virgin Islands'):
    US = US.loc[US['NAME'] != excluded]

US_county = gpd.read_file('/Users/jiananfeng/Desktop/PhD_CEE/ARPA-E_proposal/preliminary_analysis/cb_2018_us_county_500k/cb_2018_us_county_500k.shp')

US_county = US_county[US_county['STATEFP'].isin(['01','04','05','06','08','09','10','11',
                                                 '12','13','16','17','18','19','20','21',
                                                 '22','23','24','25','26','27','28','29',
                                                 '30','31','32','33','34','35','36','37',
                                                 '38','39','40','41','42','44','45','46',
                                                 '47','48','49','50','51','53','54','55','56'])]

state_ID = {'01':'AL',
            '04':'AZ',
            '05':'AR',
            '06':'CA',
            '08':'CO',
            '09':'CT',
            '10':'DE',
            '11':'DC',
            '12':'FL',
            '13':'GA',
            '16':'ID',
            '17':'IL',
            '18':'IN',
            '19':'IA',
            '20':'KS',
            '21':'KY',
            '22':'LA',
            '23':'ME',
            '24':'MD',
            '25':'MA',
            '26':'MI',
            '27':'MN',
            '28':'MS',
            '29':'MO',
            '30':'MT',
            '31':'NE',
            '32':'NV',
            '33':'NH',
            '34':'NJ',
            '35':'NM',
            '36':'NY',
            '37':'NC',
            '38':'ND',
            '39':'OH',
            '40':'OK',
            '41':'OR',
            '42':'PA',
            '44':'RI',
            '45':'SC',
            '46':'SD',
            '47':'TN',
            '48':'TX',
            '49':'UT',
            '50':'VT',
            '51':'VA',
            '53':'WA',
            '54':'WV',
            '55':'WI',
            '56':'WY'}

US_county['STATE'] = US_county['STATEFP'].apply(lambda x: state_ID[x])

balnc_area = gpd.read_file(folder + 'lpreg2/lpreg2.shp')

WRRF = WRRF.to_crs(crs='EPSG:3857')
refinery = refinery.to_crs(crs='EPSG:3857')
US = US.to_crs(crs='EPSG:3857')
US_county = US_county.to_crs(crs='EPSG:3857')
balnc_area = balnc_area.to_crs(crs='EPSG:3857')

US_county_labor_cost = gpd.read_file('/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/county_labor_cost_2022_processed.geojson')

WRRF = WRRF.sjoin_nearest(US_county_labor_cost)
WRRF = WRRF.drop(['index_right','STATE_right'], axis=1)

assert WRRF.quotient.isna().sum() == 0

WRRF = WRRF.rename({'FACILITY':'facility',
                    'CITY':'city',
                    'STATE_left':'state',
                    'CWNS_NUM':'CWNS',
                    'LATITUDE':'latitude',
                    'LONGITUDE':'longitude',
                    'FLOW_2022_MGD_FINAL':'flow_2022_MGD_final',
                    'Area\nCode':'county_code',
                    'NAME':'county',
                    'quotient':'wage_quotient'}, axis=1)

# a small number of WRRFs fall outside the boundary of the contiguous U.S. a little bit
# do not use sjoin for WRRFs, as this will remove those WRRFs
# [do not uncomment] WRRF = gpd.sjoin(WRRF, US)
# [do not uncomment] WRRF = WRRF.drop(['index_right'], axis=1)

# sjoin refinery and US
# confirmed: this only removes refineries in Alaska and Hawaii
refinery = gpd.sjoin(refinery, US)
refinery = refinery.drop(['index_right'], axis=1)

farm_fertilizer = pd.read_excel('/Users/jiananfeng/Desktop/PhD_CEE/ARPA-E_proposal/preliminary_analysis/N-P_from_fertilizer_1950-2017-july23-2020.xlsx','farm')
nonfarm_fertilizer = pd.read_excel('/Users/jiananfeng/Desktop/PhD_CEE/ARPA-E_proposal/preliminary_analysis/N-P_from_fertilizer_1950-2017-july23-2020.xlsx','nonfarm')

N_farm = farm_fertilizer[['CountyName','State','farmfertN-kg-2017']]
N_nonfarm = nonfarm_fertilizer[['CountyName','State','nonffertN-kg-2017']]

N = N_farm.merge(N_nonfarm, on=['CountyName','State'], how='inner')
N['total'] = N['farmfertN-kg-2017'] + N['nonffertN-kg-2017']

P_farm = farm_fertilizer[['CountyName','State','farmfertP-kg-2017']]
P_nonfarm = nonfarm_fertilizer[['CountyName','State','nonffertP-kg-2017']]

P = P_farm.merge(P_nonfarm, on=['CountyName','State'], how='inner')
P['total'] = P['farmfertP-kg-2017'] + P['nonffertP-kg-2017']

elec_price = pd.read_excel(folder + 'state_elec_price_2022.xlsx', 'elec_price_2022')
elec_price = elec_price.merge(US, how='right', left_on='name', right_on='NAME')
elec_price = gpd.GeoDataFrame(elec_price)

elec_GHG = WRRF[['balancing_area','kg_CO2e_kWh']].copy()
elec_GHG.drop_duplicates(inplace=True)
elec_GHG = elec_GHG.merge(balnc_area, how='right', left_on='balancing_area', right_on='PCA_REG')
elec_GHG = gpd.GeoDataFrame(elec_GHG)

US_county_labor_cost = gpd.read_file(folder + 'county_labor_cost_2022_processed.geojson')

income_tax = pd.read_excel(folder + 'state_corporate_income_tax_2022.xlsx', 'tax_2022')
income_tax = income_tax.merge(US, how='right', left_on='name', right_on='NAME')
income_tax = gpd.GeoDataFrame(income_tax)

state_ID = {'Alabama':'AL',
            'Arizona':'AZ',
            'Arkansas':'AR',
            'California':'CA',
            'Colorado':'CO',
            'Connecticut':'CT',
            'Delaware':'DE',
            'District of Columbia':'DC',
            'Florida':'FL',
            'Georgia':'GA',
            'Idaho':'ID',
            'Illinois':'IL',
            'Indiana':'IN',
            'Iowa':'IA',
            'Kansas':'KS',
            'Kentucky':'KY',
            'Louisiana':'LA',
            'Maine':'ME',
            'Maryland':'MD',
            'Massachusetts':'MA',
            'Michigan':'MI',
            'Minnesota':'MN',
            'Mississippi':'MS',
            'Missouri':'MO',
            'Montana':'MT',
            'Nebraska':'NE',
            'Nevada':'NV',
            'New Hampshire':'NH',
            'New Jersey':'NJ',
            'New Mexico':'NM',
            'New York':'NY',
            'North Carolina':'NC',
            'North Dakota':'ND',
            'Ohio':'OH',
            'Oklahoma':'OK',
            'Oregon':'OR',
            'Pennsylvania':'PA',
            'Rhode Island':'RI',
            'South Carolina':'SC',
            'South Dakota':'SD',
            'Tennessee':'TN',
            'Texas':'TX',
            'Utah':'UT',
            'Vermont':'VT',
            'Virginia':'VA',
            'Washington':'WA',
            'West Virginia':'WV',
            'Wisconsin':'WI',
            'Wyoming':'WY'}

ID_state = dict((value, key) for key, value in state_ID.items())

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
              'Maine': 1,
              'Maryland': 1,
              'Massachusetts': 1,
              'Michigan': 2,
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
              'Washington': 5,
              'West Virginia': 1,
              'Wisconsin': 2,
              'Wyoming': 4}

PADD_color = {1: b,
              2: g,
              3: r,
              4: o,
              5: y}

# the purchase of crude oil is assumed at WRRFs, therefore, use the first purchase price of the location of WRRFs
# at the same time, assume WRRFs will be responsible for the transportation of the crude oil
crude_oil_price_data = pd.read_excel(folder + 'crude_oil_price_2022.xlsx')

# state dominant nitrogen fertilizer form
state_nitrogen_fertilizer = {'AL':'UAN',
                             'AZ':'UAN',
                             'AR':'UAN',
                             'CA':'UAN',
                             'CO':'UAN',
                             'CT':'UAN',
                             'DE':'UAN',
                             'DC':'UAN',
                             'FL':'UAN',
                             'GA':'UAN',
                             'ID':'urea',
                             'IL':'NH3',
                             'IN':'NH3',
                             'IA':'NH3',
                             'KS':'UAN',
                             'KY':'UAN',
                             'LA':'UAN',
                             'MI':'NH3',
                             'ME':'UAN',
                             'MD':'UAN',
                             'MA':'UAN',
                             'MN':'NH3',
                             'MS':'UAN',
                             'MO':'NH3',
                             'MT':'urea',
                             'NE':'urea',
                             'NV':'UAN',
                             'NH':'UAN',
                             'NJ':'UAN',
                             'NM':'UAN',
                             'NY':'UAN',
                             'NC':'UAN',
                             'ND':'urea',
                             'OH':'NH3',
                             'OK':'UAN',
                             'OR':'urea',
                             'PA':'UAN',
                             'RI':'UAN',
                             'SC':'UAN',
                             'SD':'urea',
                             'TN':'UAN',
                             'TX':'UAN',
                             'UT':'UAN',
                             'VT':'UAN',
                             'VA':'UAN',
                             'WA':'urea',
                             'WV':'UAN',
                             'WI':'NH3',
                             'WY':'urea'}

WRRF['nitrogen_fertilizer'] = WRRF['state'].apply(lambda x: state_nitrogen_fertilizer[x])

fertilizer_price_uncertainty = {'AL': {'DAP':[588, 972.5, 1357],
                                       'urea':[575, 850, 1125],
                                       'UAN':[560, 705, 850]},
                                'IA': {'DAP':[830, 1000, 1170],
                                       'anhydrous_ammonia':[1065, 1386.5, 1708],
                                       'urea':[650, 867, 1084],
                                       'UAN':[532, 658.5, 785]},
                                'IL': {'DAP':[800, 962.5, 1125],
                                       'anhydrous_ammonia':[1085, 1417.5, 1750],
                                       'urea':[700, 875, 1050],
                                       'UAN':[545, 702.5, 860]},
                                'NC': {'DAP':[609, 995, 1381],
                                       'urea':[447, 855, 1263],
                                       'UAN':[350, 605, 860]},
                                'OK': {'DAP':[800, 1022, 1244],
                                       'anhydrous_ammonia':[1100, 1312.5, 1525],
                                       'urea':[700, 882.5, 1065],
                                       'UAN':[415, 607.5, 800]},
                                'SC': {'DAP':[820, 1002.5, 1185],
                                       'urea':[650, 895, 1140],
                                       'UAN':[450, 587.5, 725]}}

# sludge disposal cost in $/kg sludge (Peccia and Westerhoff. 2015, https://pubs.acs.org/doi/full/10.1021/acs.est.5b01931)
sludge_disposal_cost = {'landfill': 0.413,
                        'land_application': 0.606,
                        'incineration': 0.441}

# TODO: are these just fugutive emissions? is it a fair comparison between these values and the emission from the HTL-based system, which includes chemical, transportation, etc.
# sludge emission factor in kg CO2 eq/kg sludge (values from IEDO)
sludge_emission_factor = {'landfill': 5.65/1000*29.8,
                          'land_application': 0.05*0.01*44/28*273,
                          'incineration': 0}

#%% hauling contribution to GHG emissions from landfill and land application

# kg CO2 eq/tonne/km
hauling_CI = 0.13004958

sludge_dw = 0.2

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_xlim(0, 80)
ax.set_ylim(0, 240)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Hauling\ distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Management\ CI}$' + '\n[kg CO${_2}$ eq·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 90, 10))
plt.yticks(np.arange(0, 280, 40))

ax_top = ax.twiny()
ax_top.set_xlim(0, 80)
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 90, 10))
plt.yticks(np.arange(0, 280, 40))

ax_right = ax.twinx()
ax_right.set_ylim(0, 240)
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 90, 10))
plt.yticks(np.arange(0, 280, 40))

ax.plot([0, 80],
        [sludge_emission_factor['landfill']*1000, sludge_emission_factor['landfill']*1000],
        c='k',
        linewidth=3)

ax.plot([0, 80],
        [sludge_emission_factor['land_application']*1000, sludge_emission_factor['land_application']*1000],
        c='k',
        linewidth=3)

ax.plot([0, 80],
        [0*hauling_CI/sludge_dw, 80*hauling_CI/sludge_dw],
        c='k',
        linewidth=3)

ax.plot([8, 8],
        [0, 240],
        c=db,
        linestyle='--',
        linewidth=3)

ax.plot([13, 13],
        [0, 240],
        c=dr,
        linestyle='--',
        linewidth=3)

ax.scatter(8,
           8*hauling_CI/sludge_dw,
           c=b,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

ax.scatter(8,
           sludge_emission_factor['landfill']*1000,
           c=b,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

ax.scatter(8,
           sludge_emission_factor['land_application']*1000,
           c=b,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

ax.scatter(13,
           13*hauling_CI/sludge_dw,
           c=r,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

ax.scatter(13,
           sludge_emission_factor['landfill']*1000,
           c=r,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

ax.scatter(13,
           sludge_emission_factor['land_application']*1000,
           c=r,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

#%% WRRFs visualization

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)
# for all WRRFs together with the same symbols
# US.plot(ax=ax, color='w', edgecolor='k', linewidth=6)

WRRF = WRRF.sort_values(by='flow_2022_MGD_final', ascending=False)

WRRF_flow = WRRF['flow_2022_MGD_final']

more_than_100 = WRRF_flow >= 100
WRRF[more_than_100].plot(ax=ax, color=dg, markersize=WRRF.loc[more_than_100,'flow_2022_MGD_final']*10, edgecolor='k', linewidth=1.5, alpha=1)

between_10_and_100 = (WRRF_flow >= 10) & (WRRF_flow < 100)
WRRF[between_10_and_100].plot(ax=ax, color=g, markersize=WRRF.loc[between_10_and_100,'flow_2022_MGD_final']*10, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_10 = WRRF_flow < 10
WRRF[less_than_10].plot(ax=ax, color=g, markersize=WRRF.loc[less_than_10,'flow_2022_MGD_final']*10, edgecolor='none', alpha=0.2)

# comment out the code above and uncomment the following line to show all WRRFs together with the same symbol
# WRRF.plot(ax=ax, color=a, markersize=10)

rectangle_edge = Rectangle((-13320000, 2880000), 1280000, 680000,
                           color='k', lw=3, fc='none', alpha=1)
ax.add_patch(rectangle_edge)

# note the size of legends does not match exactly with the size of WRRFs
ax.scatter(x=-13140000, y=3420000, marker='o', s=200*10, c=dg, linewidths=3,
           alpha=1, edgecolor='k')
ax.scatter(x=-13140000, y=3220000, marker='o', s=50*10, c=g, linewidths=3,
           alpha=0.5, edgecolor='k')
ax.scatter(x=-13140000, y=3020000, marker='o', s=10*10, c=g, linewidths=3,
           alpha=0.2, edgecolor='none')

plt.figtext(0.261, 0.37, '≥ 100 MGD', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})
plt.figtext(0.261, 0.3485, '10 to 100 MGD', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})
plt.figtext(0.261, 0.327, '< 10 MGD', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})

ax.set_aspect(1)

ax.set_axis_off()

#%% oil refinery visualization

# note oil refinery production numbers reflect the production at some certain time points, not necessarily the capacity of refineries

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax,
        color=[PADD_color[state_PADD[i]] for i in list(US.NAME)],
        alpha=0.3,
        edgecolor='k',
        linewidth=0)

US.plot(ax=ax, color='none', edgecolor='k', linewidth=3)
# for TOC
# US.plot(ax=ax, color='none', edgecolor='k', linewidth=0)
# for all oil refineries together with the same symbols
# US.plot(ax=ax, color='none', edgecolor='k', linewidth=6)

refinery.plot(ax=ax, color=o, markersize=1000, edgecolor='k', linewidth=1.5, alpha=1)

# comment out the code above and uncomment the following line to show all oil refineries together with the same symbols
# refinery.plot(ax=ax, color=o, markersize=1000, edgecolor='k', linewidth=6)

ax.set_aspect(1)

ax.set_axis_off()

#%% WRRFs+oil refineries visualization (all)

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)
WRRF.plot(ax=ax, color=a, markersize=5)
refinery.plot(ax=ax, color=o, markersize=500, edgecolor='k', linewidth=1.5)

ax.set_aspect(1)

ax.set_axis_off()

#%% N needs

US_county_N = US_county.merge(N, left_on=['NAME','STATE'], right_on=['CountyName','State'], how='left')
US_county_N = US_county_N[['NAME','STATE','total','geometry']]

US_county_N['total'] = US_county_N['total'].fillna(0)

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', b, db])

US_county_N.plot(column='total', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=0.5)

fig.axes[1].set_ylabel('$\mathbf{Nitrogen}$ [kg]', fontname='Arial', fontsize=41)
fig.axes[1].tick_params(length=7.5, width=1.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% P needs

US_county_P = US_county.merge(P, left_on=['NAME','STATE'], right_on=['CountyName','State'], how='left')
US_county_P = US_county_P[['NAME','STATE','total','geometry']]

US_county_P['total'] = US_county_P['total'].fillna(0)

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', r, dr])

US_county_P.plot(column='total', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=0.5)

fig.axes[1].set_ylabel('$\mathbf{Phosphorus}$ [kg]', fontname='Arial', fontsize=41)
fig.axes[1].tick_params(length=7.5, width=1.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% electricity price visualization

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', o, do])

elec_price.plot(column='price', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=3)

fig.axes[1].set_ylabel('$\mathbf{Electricity\ price}$ [cent·${kWh^{-1}}$]', fontname='Arial', fontsize=35)
fig.axes[1].tick_params(length=10, width=3)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% electricity carbon intensity visualization

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', b, db])

elec_GHG.plot(column='kg_CO2e_kWh', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=3)

fig.axes[1].set_ylabel('$\mathbf{Electricity\ carbon\ intensity}$' + '\n[kg ${CO_2}$ eq·${kWh^{-1}}$]', fontname='Arial', fontsize=35, linespacing=0.8)
fig.axes[1].tick_params(length=10, width=3)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% labor wage visualization

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', p, dp])

US_county_labor_cost.plot(column='quotient', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=0.5)

fig.axes[1].set_ylabel('$\mathbf{Relative\ labor\ wage}$ [%]', fontname='Arial', fontsize=35)
fig.axes[1].tick_params(length=10, width=3)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% tax rate visualization (just for demonstration, not for model)

income_tax['tax_percentage'] = income_tax['tax']*100

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', g, dg])

income_tax.plot(column='tax_percentage', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=3)

fig.axes[1].set_ylabel('$\mathbf{Corporate\ income\ tax}$ [%]', fontname='Arial', fontsize=35)
fig.axes[1].tick_params(length=10, width=3)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

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
# # code to generate the inventory
# 
# # !!! get a google API key
# # !!! do not upload to GitHub
# gmaps = googlemaps.Client(key='XXX')
# 
# linear_distance = []
# real_distance = []
# 
# for i in range(len(WRRF_input)):
#     linear_distance.append(geopy.distance.geodesic(WRRF_input['WRRF_location'].iloc[i], WRRF_input['refinery_location'].iloc[i]).km)
#     
#     try:
#         print(i)
#         real_distance.append(gmaps.distance_matrix(WRRF_input['WRRF_location'].iloc[i], WRRF_input['refinery_location'].iloc[i],
#                                                    mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000)
#     except KeyError:
#         print('--------------------------------')
#         real_distance.append(np.nan)
# 
# WRRF_input['linear_distance_km'] = linear_distance
# WRRF_input['real_distance_km'] = real_distance
# 
# distance_inventory = WRRF_input[['CWNS','Site ID','linear_distance_km','real_distance_km']]
# 
# distance_inventory.to_excel(folder + f'distance_inventory_{date.today()}.xlsx')
# =============================================================================

distance_inventory = pd.read_excel(folder + 'distance_inventory_2025-02-10.xlsx')

# match using WRRF ID ('CWNS') and oil refinery ID ('Site ID')
WRRF_input = WRRF_input.merge(distance_inventory, how='left', on=['CWNS','Site ID'])

missing_distance = []
for i in WRRF_input.index:
    if pd.isna(WRRF_input.loc[i,'linear_distance_km']):
        missing_distance.append(i)

if len(missing_distance) == 0:
    WRRF_input.to_excel(folder + f'HTL_geospatial_model_input_{date.today()}.xlsx')
else:
    # !!! get a google API key
    # !!! do not upload to GitHub
    gmaps = googlemaps.Client(key='XXX')
    
    for i in missing_distance:
        WRRF_input.loc[i,'linear_distance_km'] = geopy.distance.geodesic(WRRF_input.loc[i,'WRRF_location'],
                                                                         WRRF_input.loc[i,'refinery_location']).km
        
        try:
            print(i)
            WRRF_input.loc[i,'real_distance_km'] = gmaps.distance_matrix(WRRF_input.loc[i,'WRRF_location'],
                                                                         WRRF_input.loc[i,'refinery_location'],
                                                                         mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000
        except KeyError:
            print('--------------------------------')
            WRRF_input.loc[i,'real_distance_km'] = np.nan
    
    distance_inventory = WRRF_input[['CWNS','Site ID','linear_distance_km','real_distance_km']]
    
    distance_inventory.to_excel(folder + f'distance_inventory_{date.today()}.xlsx')
    
    # input for following analyses
    WRRF_input.to_excel(folder + f'HTL_geospatial_model_input_{date.today()}.xlsx')

#%% travel distance box plot

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
print(f"{WRRF_input['real_distance_km'].notna().sum()} WRRFs included")
print(f"{WRRF_input['real_distance_km'].isna().sum()} WRRFs excluded")
print(WRRF_input[WRRF_input['real_distance_km'].isna()])
print(WRRF_input[WRRF_input['real_distance_km'].isna()]['flow_2022_MGD_final'])
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')
print(f"{((WRRF_input.sludge_anaerobic_digestion == 1) & (WRRF_input.sludge_aerobic_digestion == 0)).sum()} WRRFs just have AD")
print(f"{((WRRF_input.sludge_anaerobic_digestion == 0) & (WRRF_input.sludge_aerobic_digestion == 1)).sum()} WRRFs just have AeD")
print(f"{((WRRF_input.sludge_anaerobic_digestion == 1) & (WRRF_input.sludge_aerobic_digestion == 1)).sum()} WRRFs have both AD and AeD")
print(f"{((WRRF_input.sludge_anaerobic_digestion == 0) & (WRRF_input.sludge_aerobic_digestion == 0)).sum()} WRRFs have neither AD nor AeD")

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize = (5, 8))

ax = plt.gca()
ax.set_ylim(0, 800)
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
    
ax_right.scatter(x=1,
                 y=WRRF_input['real_distance_km'].mean(),
                 marker='*',
                 s=600,
                 c='w',
                 linewidths=3,
                 alpha=1,
                 edgecolor='k',
                 zorder=2)

# uncomment for outliers
# for flier in bp['fliers']:
#     flier.set(marker='o', markersize=7, markerfacecolor='k', markeredgewidth=1.5)
    
# fig.savefig('/Users/jiananfeng/Desktop/distance.png', transparent=True, bbox_inches='tight')

#%% travel distance box plot (per region)

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# results grouped by different regions (not by PADD regions, just they are the same as PADD regions)
WRRF_input.loc[WRRF_input['state'].isin(['CT','DC','DE','FL','GA','MA','MD','ME','NC','NH','NJ','NY','PA','RI','SC','VA','VT','WV']),'WRRF_PADD'] = 1
WRRF_input.loc[WRRF_input['state'].isin(['IA','IL','IN','KS','KY','MI','MN','MO','ND','NE','OH','OK','SD','TN','WI']),'WRRF_PADD'] = 2
WRRF_input.loc[WRRF_input['state'].isin(['AL','AR','LA','MS','NM','TX']),'WRRF_PADD'] = 3
WRRF_input.loc[WRRF_input['state'].isin(['CO','ID','MT','UT','WY']),'WRRF_PADD'] = 4
WRRF_input.loc[WRRF_input['state'].isin(['AZ','CA','NV','OR','WA']),'WRRF_PADD'] = 5

fig = plt.figure(figsize=(20, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

def add_region(position, region, color):
    ax = fig.add_subplot(gs[0, position])
    
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    ax = plt.gca()
    ax.set_ylim(0, 1200)
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
    
    bp = ax.boxplot(WRRF_input[WRRF_input['WRRF_PADD'] == position+1]['real_distance_km'], showfliers=False, widths=0.5, patch_artist=True)
    
    for box in bp['boxes']:
        box.set(color='k', facecolor=color, linewidth=3)

    for whisker in bp['whiskers']:
        whisker.set(color='k', linewidth=3)

    for median in bp['medians']:
        median.set(color='k', linewidth=3)
        
    for cap in bp['caps']:
        cap.set(color='k', linewidth=3)
        
    ax.scatter(x=1,
               y=WRRF_input[WRRF_input['WRRF_PADD'] == position+1]['real_distance_km'].mean(),
               marker='*',
               s=600,
               c='w',
               linewidths=3,
               alpha=1,
               edgecolor='k',
               zorder=3)
        
    # for flier in bp['fliers']:
    #     flier.set(marker='o', markersize=7, markerfacecolor=color, markeredgewidth=1.5)

add_region(0, 'East Coast', b)
add_region(1, 'Midwest', g)
add_region(2, 'Gulf Coast', r)
add_region(3, 'Rocky Mountain', o)
add_region(4, 'West Coast', y)

#%% WRRFs GHG map

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
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
WRRF_input[more_than_500].plot(ax=ax, color=dp, markersize=WRRF_input.loc[more_than_500,'total_emission']/350, edgecolor='k', linewidth=2, alpha=1)

between_50_and_500 = (WRRF_GHG_tonne_per_day > 50) & (WRRF_GHG_tonne_per_day <= 500)
WRRF_input[between_50_and_500].plot(ax=ax, color=p, markersize=WRRF_input.loc[between_50_and_500,'total_emission']/350, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_50 = WRRF_GHG_tonne_per_day <= 50
WRRF_input[less_than_50].plot(ax=ax, color=p, markersize=WRRF_input.loc[less_than_50,'total_emission']/350, alpha=0.1)

ax.set_aspect(1)

ax.set_axis_off()

#%% WRRFs sludge management GHG map

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
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
WRRF_input[more_than_20].plot(ax=ax, color=dr, markersize=WRRF_input.loc[more_than_20,'biosolids_emission']/10, edgecolor='k', linewidth=2, alpha=1)

between_2_and_20 = (sludge_GHG_tonne_per_day > 2) & (sludge_GHG_tonne_per_day <= 20)
WRRF_input[between_2_and_20].plot(ax=ax, color=r, markersize=WRRF_input.loc[between_2_and_20,'biosolids_emission']/10, edgecolor='k', linewidth=1.5, alpha=0.5)

less_than_2 = sludge_GHG_tonne_per_day <= 2
WRRF_input[less_than_2].plot(ax=ax, color=r, markersize=WRRF_input.loc[less_than_2,'biosolids_emission']/10, alpha=0.1)

ax.set_aspect(1)

ax.set_axis_off()

#%% cumulative WRRFs capacity vs distances (data processing)

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

result = WRRF_input[['Site ID']].drop_duplicates()

max_distance = 1500
assert max_distance > WRRF_input.real_distance_km.max(), 'update max_distance'

for distance in np.linspace(0, max_distance, max_distance+1):
    WRRF_input_distance = WRRF_input[WRRF_input['real_distance_km'] <= distance]
    WRRF_input_distance = WRRF_input_distance.groupby('Site ID').sum('flow_2022_MGD_final')
    WRRF_input_distance = WRRF_input_distance[['flow_2022_MGD_final']]
    WRRF_input_distance = WRRF_input_distance.rename(columns={'flow_2022_MGD_final': int(distance)})
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

result.to_excel(folder + f'results/distance/MGD_vs_real_distance_{max_distance}_km.xlsx')

#%% make the plot of cumulative WRRFs capacity vs distances (data preparation)

# !!! update the file if necessary
CF_input = pd.read_excel(folder + 'results/distance/MGD_vs_real_distance_1500_km.xlsx')

CF_input[0] = 0

max_distance_plot = 1500
assert max_distance_plot > WRRF_input.real_distance_km.max(), 'update max_distance'

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
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    ax.set_xlim(0, max_distance_plot)
    ax.set_ylim(0, 7000)
    
    plt.xticks(np.arange(0, max_distance_plot*1.2, max_distance_plot*0.2))
    plt.yticks(np.arange(0, 8000, 1000))
    
    for i in range(start_region, end_region):
        CF_input.iloc[2:, i].plot(ax=ax, color=color, linewidth=3)
    
    if position == 4:
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True, pad=0)
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        
        ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, labelleft=True, labelbottom=True, pad=0)
        ax.set_ylabel(r'$\mathbf{Cumulative\ WRRFs\ capacity}$ [MGD]', fontname='Arial', fontsize=35)
    else:
        if position == 2:
            ax.set_xlabel(r'$\mathbf{Travel\ distance}$ [km]', fontname='Arial', fontsize=35)
        ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True, pad=0)
        plt.xticks(np.arange(max_distance_plot*0.2, max_distance_plot*1.2, max_distance_plot*0.2))
    
    for label in ax.get_xticklabels():
        label.set_rotation(45)
    
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    plt.xticks(np.arange(max_distance_plot*0.2, max_distance_plot*1.2, max_distance_plot*0.2))
    ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

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

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(11, 10))

ax.set_xlim(0, max_distance_plot)
ax.set_ylim(0, 20000)

plt.xticks(np.arange(0, max_distance_plot*1.2, max_distance_plot*0.2))
plt.yticks(np.arange(0, 24000, 4000))

CF_input.iloc[2:, 0:PADD_1].sum(axis=1).plot(ax=ax, color=b, linewidth=3)
CF_input.iloc[2:, PADD_1:PADD_2].sum(axis=1).plot(ax=ax, color=g, linewidth=3)
CF_input.iloc[2:, PADD_2:PADD_3].sum(axis=1).plot(ax=ax, color=r, linewidth=3)
CF_input.iloc[2:, PADD_3:PADD_4].sum(axis=1).plot(ax=ax, color=o, linewidth=3)
CF_input.iloc[2:, PADD_4:PADD_5].sum(axis=1).plot(ax=ax, color=y, linewidth=3)

ax.set_xlabel(r'$\mathbf{Travel\ distance}$ [km]', fontname='Arial', fontsize=35)
ax.set_ylabel(r'$\mathbf{Cumulative\ WRRFs\ capacity}$ [MGD]', fontname='Arial', fontsize=35)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
plt.yticks(np.arange(0, 24000, 4000))
ax_right.tick_params(direction='in', length=10, width=3, bottom=True, top=False, left=False, right=True, labelcolor='none')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, max_distance_plot*1.2, max_distance_plot*0.2))
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

#%% CO2 abatement cost analysis

filterwarnings('ignore')

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# if just want to see the two plants in Urbana-Champaign:
# WRRF_input = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

print(len(WRRF_input))

# $/tonne
WRRF_input['waste_cost'] = sum(WRRF_input[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000
# kg CO2 eq/tonne
WRRF_input['waste_GHG'] =  sum(WRRF_input[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

CWNS = []
cost = []
CI = []
saving = []
CO2_reduction = []
sludge_CO2_reduction_ratio = []
WRRF_CO2_reduction_ratio = []
USD_decarbonization = []
income_tax = []

# confirmed: the baseline results from the model are the same as directly calculated from the system
# (since the baseline values in the model are the same as parameters in the system)
# !!! run in different consoles to speed up: 0, 5000, 10000, len(WRRF_input)
for i in range(0, len(WRRF_input)):
    sys = create_geospatial_system(size=WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']/1000/365,
                                   sludge_transportation=0,
                                   sludge_distance=100,
                                   biocrude_distance=WRRF_input.iloc[i]['real_distance_km'],
                                   anaerobic_digestion=WRRF_input.iloc[i]['sludge_anaerobic_digestion'],
                                   aerobic_digestion=WRRF_input.iloc[i]['sludge_aerobic_digestion'],
                                   ww_2_dry_sludge_ratio=1,
                                   state=WRRF_input.iloc[i]['state'],
                                   nitrogen_fertilizer=WRRF_input.iloc[i]['nitrogen_fertilizer'],
                                   elec_GHG=WRRF_input.iloc[i]['kg_CO2e_kWh'],
                                   wage_adjustment=WRRF_input.iloc[i]['wage_quotient']/100)
    
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    WWTP = unit.WWTP
    BiocrudeTank = unit.BiocrudeTank
    raw_wastewater = stream.raw_wastewater
    biocrude = stream.biocrude
    DAP = stream.DAP
    tea = sys.TEA
    lca = sys.LCA
    
    # tonne
    sludge_tonne = raw_wastewater.F_vol*_m3perh_to_MGD*WWTP.ww_2_dry_sludge*(sys.operating_hours/24)*lca.lifetime
    
    try:
        # $/tonne
        sludge_cost = -tea.solve_price(raw_wastewater)*water_density*_MMgal_to_L/WWTP.ww_2_dry_sludge
    except RuntimeError:
        # due to high gross receipts tax (higher than net income, it might be impossible to break even)
        print('-------RUNTIME ERROR-------')
        
    try:
        saving_result = sludge_tonne*(WRRF_input.iloc[i]['waste_cost'] - sludge_cost)
    except NameError:
        saving_result = np.nan
        
    # kg CO2 eq/tonne
    sludge_CI = lca.get_total_impacts(operation_only=True,
                                      exclude=(raw_wastewater,),
                                      annual=True)['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)
    
    CO2_reduction_result = sludge_tonne*(WRRF_input.iloc[i]['waste_GHG'] - sludge_CI)
    
    try:
        sludge_CO2_reduction_ratio_result = CO2_reduction_result/WRRF_input.iloc[i]['biosolids_emission']/(sys.operating_hours/24)/lca.lifetime
    # some WRRFs have incineration and the biosolids_emission is 0
    except FloatingPointError:
        sludge_CO2_reduction_ratio_result = np.nan
    
    WRRF_CO2_reduction_ratio_result = CO2_reduction_result/WRRF_input.iloc[i]['total_emission']/(sys.operating_hours/24)/lca.lifetime
    
    # make sure CO2_reduction_result is positive to calculate USD_per_tonne_CO2_reduction
    if CO2_reduction_result > 0:
        # save money: the results are negative; spend money: the results are positive
        USD_per_tonne_CO2_reduction = -saving_result/CO2_reduction_result*1000
    else:
        USD_per_tonne_CO2_reduction = np.nan
    
    CWNS.append(WRRF_input.iloc[i]['CWNS'])
    try:
        cost.append(sludge_cost)
    except NameError:
        cost.append(np.nan)
    CI.append(sludge_CI)
    saving.append(saving_result)
    CO2_reduction.append(CO2_reduction_result)
    sludge_CO2_reduction_ratio.append(sludge_CO2_reduction_ratio_result)
    WRRF_CO2_reduction_ratio.append(WRRF_CO2_reduction_ratio_result)
    USD_decarbonization.append(USD_per_tonne_CO2_reduction)
    income_tax.append(tea.income_tax)
    
    if USD_per_tonne_CO2_reduction < 0:
        print('HERE WE GO!')
    
    # check progress
    if i%50 == 0:
        print(i)

result = {'CWNS': CWNS,
          'cost': cost,
          'CI': CI,
          'saving': saving,
          'CO2_reduction': CO2_reduction,
          'sludge_CO2_reduction_ratio': sludge_CO2_reduction_ratio,
          'WRRF_CO2_reduction_ratio': WRRF_CO2_reduction_ratio,
          'USD_decarbonization': USD_decarbonization,
          'income_tax': income_tax}
        
result = pd.DataFrame(result)

result.to_excel(folder + f'results/baseline/baseline_{date.today()}_{i}.xlsx')

#%% merge the results and the input

# !!! update the input file if necessary
input_data = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
input_data = input_data.dropna(subset='real_distance_km')

# !!! update these files if necessary
output_result_1 = pd.read_excel(folder + 'results/baseline/baseline_2025-02-16_4999.xlsx')
output_result_2 = pd.read_excel(folder + 'results/baseline/baseline_2025-02-16_9999.xlsx')
output_result_3 = pd.read_excel(folder + 'results/baseline/baseline_2025-02-16_15858.xlsx')

output_result = pd.concat([output_result_1, output_result_2, output_result_3])

assert len(input_data) == len(output_result)

integrated_result = input_data.merge(output_result, how='left', on='CWNS')

integrated_result.to_excel(folder + f'results/baseline/integrated_baseline_{date.today()}.xlsx')

#%% decarbonization map (preparation)

# !!! update the file here if necessary
decarbonization_map = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
decarbonization_map = decarbonization_map[decarbonization_map['USD_decarbonization'].notna()]
decarbonization_map = decarbonization_map[decarbonization_map['USD_decarbonization'] <= 0]
decarbonization_map = decarbonization_map.sort_values(by='total_emission', ascending=False).copy()
decarbonization_map = gpd.GeoDataFrame(decarbonization_map, crs='EPSG:4269',
                                       geometry=gpd.points_from_xy(x=decarbonization_map.longitude,
                                                                   y=decarbonization_map.latitude))
decarbonization_map = decarbonization_map.to_crs(crs='EPSG:3857')

decarbonization_map['CO2_reduction_tonne_per_day'] = decarbonization_map['CO2_reduction']/30/365/1000

def plot_map(dataset, color):    
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30

    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    fig, ax = plt.subplots(figsize=(30, 30))

    US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)

    dataset.plot(ax=ax, color=color, markersize=dataset['CO2_reduction_tonne_per_day']*200, edgecolor='k', linewidth=1.5, alpha=1)
    
    color_1 = color_2 = color_3 = color_4 = 'none'
    
    max_size = (dataset['CO2_reduction_tonne_per_day']*200).max()
    min_size = (dataset['CO2_reduction_tonne_per_day']*200).min()
    
    if max_size > 50*200:
        raise ValueError('add another layer of legend')
    elif max_size > 25*200:
        color_1 = color
    elif max_size > 10*200:
        color_2 = color
    elif max_size > 1*200:
        color_3 = color
    else:
        color_4 = color
        
    if min_size > 50*200:
        color_1 = 'w'
        raise ValueError('add another layer of legend')
    elif min_size > 25*200:
        color_2 = 'w'
    elif min_size > 10*200:
        color_3 = 'w'
    elif min_size > 1*200:
        color_4 = 'w'
    
    rectangle_edge = Rectangle((-13900000, 2848000), 2150000, 750000,
                               color='k', lw=3, fc='none', alpha=1)
    ax.add_patch(rectangle_edge)
    
    ax.scatter(x=-13580000, y=3100000, marker='o', s=50*200, c=color_1, linewidths=3,
               alpha=1, edgecolor='k')
    ax.scatter(x=-13580000, y=3100000, marker='o', s=25*200, c=color_2, linewidths=3,
               alpha=1, edgecolor='k')
    ax.scatter(x=-13580000, y=3100000, marker='o', s=10*200, c=color_3, linewidths=3,
               alpha=1, edgecolor='k')
    ax.scatter(x=-13580000, y=3100000, marker='o', s=1*200, c=color_4, linewidths=3,
               alpha=1, edgecolor='k')
    
    plt.figtext(0.163, 0.372, '[tonne ${CO_2}$ eq·${day^{-1}}$]', fontname='Arial', fontdict={'fontsize': 50,'color':'k','fontweight':'bold'})
    plt.figtext(0.23, 0.345, '1st: 50  2nd: 25', fontname='Arial', fontdict={'fontsize': 50,'color':'k','style':'italic'})
    plt.figtext(0.23, 0.319, '3rd:10  4th: 1', fontname='Arial', fontdict={'fontsize': 50,'color':'k','style':'italic'})
    
    ax.set_aspect(1)
    
    ax.set_axis_off()

#%% decarbonization map (all)

all_map = decarbonization_map.copy()
plot_map(all_map, g)

#%% decarbonization map (AD only)

AD_map = decarbonization_map[(decarbonization_map['sludge_anaerobic_digestion'] == 1) & (decarbonization_map['sludge_aerobic_digestion'] == 0)].copy()
plot_map(AD_map, b)

#%% decarbonization map (AeD only)

AeD_map = decarbonization_map[(decarbonization_map['sludge_anaerobic_digestion'] == 0) & (decarbonization_map['sludge_aerobic_digestion'] == 1)].copy()
plot_map(AeD_map, y)

#%% decarbonization map (both AD and AeD)

AD_and_AeD_map = decarbonization_map[(decarbonization_map['sludge_anaerobic_digestion'] == 1) & (decarbonization_map['sludge_aerobic_digestion'] == 1)].copy()
plot_map(AD_and_AeD_map, g)

#%% decarbonization map (AD or AeD)

AD_or_AeD_map = decarbonization_map[(decarbonization_map['sludge_anaerobic_digestion'] == 1) | (decarbonization_map['sludge_aerobic_digestion'] == 1)].copy()
plot_map(AD_or_AeD_map, b)

#%% decarbonization map (no AD and no AeD)

none_map = decarbonization_map[(decarbonization_map['sludge_anaerobic_digestion'] == 0) & (decarbonization_map['sludge_aerobic_digestion'] == 0)].copy()
plot_map(none_map, r)

#%% cumulative GHG reduction

# !!! update the file here if necessary
decarbonization_map = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
decarbonization_map = decarbonization_map[decarbonization_map['USD_decarbonization'].notna()]
decarbonization_map = decarbonization_map[decarbonization_map['USD_decarbonization'] <= 0]
decarbonization_map = decarbonization_map.sort_values(by='total_emission', ascending=False).copy()
decarbonization_map = gpd.GeoDataFrame(decarbonization_map, crs='EPSG:4269',
                                       geometry=gpd.points_from_xy(x=decarbonization_map.longitude,
                                                                   y=decarbonization_map.latitude))
decarbonization_map = decarbonization_map.to_crs(crs='EPSG:3857')

decarbonization_map['CO2_reduction_tonne_per_day'] = decarbonization_map['CO2_reduction']/30/365/1000

sorted_data = decarbonization_map.sort_values(by='CO2_reduction_tonne_per_day', ascending=False).reset_index(drop=True)
sorted_data['cumulative_emissions'] = sorted_data['CO2_reduction_tonne_per_day'].cumsum()

sorted_data['facility_rank'] = sorted_data.index + 1

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(11, 10))

ax = plt.gca()
ax.set_xlim(0, 600)
ax.set_ylim(0, 1400)

ax.tick_params(direction='inout', length=20, width=3,
                bottom=True, top=False, left=True, right=False, pad=0)

plt.xticks(np.arange(0, 700, 100))
plt.yticks(np.arange(0, 1600, 200))

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 700, 100))

ax_top.tick_params(direction='in', length=10, width=3,
                    bottom=False, top=True, left=False, right=False,
                    labelcolor='none')

ax_bottom = ax.twinx()
ax_bottom.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0, 700, 100))

ax_bottom.tick_params(direction='in', length=10, width=3,
                      bottom=False, top=False, left=False, right=True,
                      labelcolor='none')

plt.plot(sorted_data['facility_rank'],
         sorted_data['cumulative_emissions'],
         linewidth=0,
         marker='o',
         color='k',
         markersize=10)

percentile_10 = floor(len(sorted_data)/10) - 1
percentile_50 = floor(len(sorted_data)/2) - 1

# 10th percentile of the facility number
plt.plot(sorted_data['facility_rank'].iloc[percentile_10],
         sorted_data['cumulative_emissions'].iloc[percentile_10],
         linewidth=0,
         marker='o',
         color=r,
         markeredgecolor='k',
         markeredgewidth=3,
         markersize=20)

# 50th percentile of the facility number
plt.plot(sorted_data['facility_rank'].iloc[percentile_50],
         sorted_data['cumulative_emissions'].iloc[percentile_50],
         linewidth=0,
         marker='o',
         color=r,
         markeredgecolor='k',
         markeredgewidth=3,
         markersize=20)

plt.plot([0, sorted_data['facility_rank'].iloc[percentile_10]],
         [sorted_data['cumulative_emissions'].iloc[percentile_10], sorted_data['cumulative_emissions'].iloc[percentile_10]],
         lw=3, color='k', linestyle='--', solid_capstyle='round', zorder=0)

plt.plot([0, sorted_data['facility_rank'].iloc[percentile_50]],
         [sorted_data['cumulative_emissions'].iloc[percentile_50], sorted_data['cumulative_emissions'].iloc[percentile_50]],
         lw=3, color='k', linestyle='--', solid_capstyle='round', zorder=0)

ax.set_xlabel('$\mathbf{Facitity\ number}$',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{Cumulative\ GHG\ reduction}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

print(sorted_data['cumulative_emissions'].iloc[percentile_10]/sorted_data['cumulative_emissions'].iloc[-1])
print(sorted_data['cumulative_emissions'].iloc[percentile_50]/sorted_data['cumulative_emissions'].iloc[-1])

#%% GHG reduction boxplots

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize = (7, 10))

ax = plt.gca()
ax.set_ylim(0, 100)
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

ax.set_xticklabels(['solids','WRRF'])
plt.yticks(np.arange(0, 120, 20))

ax.set_ylabel(r'$\mathbf{GHG\ reduction}$ [%]', fontname='Arial', fontsize=45)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

bp = plt.boxplot([decarbonization_map['sludge_CO2_reduction_ratio'].dropna()*100, decarbonization_map['WRRF_CO2_reduction_ratio'].dropna()*100], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=b, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=3)
    
ax_right.scatter(x=1,
                 y=decarbonization_map['sludge_CO2_reduction_ratio'].mean()*100,
                 marker='*',
                 s=600,
                 c='w',
                 linewidths=3,
                 alpha=1,
                 edgecolor='k',
                 zorder=3)

ax_right.scatter(x=2,
                 y=decarbonization_map['WRRF_CO2_reduction_ratio'].mean()*100,
                 marker='*',
                 s=600,
                 c='w',
                 linewidths=3,
                 alpha=1,
                 edgecolor='k',
                 zorder=3)

# uncomment for outliers
# for flier in bp['fliers']:
#     flier.set(marker='o', markersize=7, markerfacecolor='k', markeredgewidth=1.5)
    
# fig.savefig('/Users/jiananfeng/Desktop/distance.png', transparent=True, bbox_inches='tight')

#%% GHG reduction zoom-in boxplot

plt.rcParams['axes.linewidth'] = 6
plt.rcParams['xtick.labelsize'] = 75
plt.rcParams['ytick.labelsize'] = 75

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize = (5, 10))

ax = plt.gca()
ax.set_ylim(0, 6)
ax.tick_params(direction='inout', length=40, width=6, labelbottom=False, bottom=False, top=False, left=True, right=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=20, width=6, bottom=False, top=True, left=False, right=True, labelcolor='none')

bp = plt.boxplot(decarbonization_map['WRRF_CO2_reduction_ratio'].dropna()*100, showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=b, linewidth=6)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=6)

for median in bp['medians']:
    median.set(color='k', linewidth=6)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=6)
    
ax_right.scatter(x=1,
                 y=decarbonization_map['WRRF_CO2_reduction_ratio'].mean()*100,
                 marker='*',
                 s=2200,
                 c='w',
                 linewidths=6,
                 alpha=1,
                 edgecolor='k',
                 zorder=3)

# uncomment for outliers
# for flier in bp['fliers']:
#     flier.set(marker='o', markersize=7, markerfacecolor='k', markeredgewidth=1.5)
    
# fig.savefig('/Users/jiananfeng/Desktop/distance.png', transparent=True, bbox_inches='tight')

#%% cost vs sludge amount

# !!! update the file here if necessary
saving_vs_sludge = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
saving_vs_sludge = saving_vs_sludge[saving_vs_sludge['USD_decarbonization'].notna()]
saving_vs_sludge = saving_vs_sludge[saving_vs_sludge['USD_decarbonization'] <= 0]

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_xlim(0, 250)
ax.set_ylim(-200, 700)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Solids}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Management\ cost}$' + '\n[\$·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(-200, 800, 100))

ax_top = ax.twiny()
ax_top.set_xlim(0, 250)
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(-200, 800, 100))

ax_right = ax.twinx()
ax_right.set_ylim(-200, 700)
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(-200, 800, 100))

ax.scatter(x=saving_vs_sludge['total_sludge_amount_kg_per_year']/1000/365,
           y=saving_vs_sludge['cost'],
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% CI vs sludge amount

# !!! update the file here if necessary
decarbonization_vs_sludge = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
decarbonization_vs_sludge = decarbonization_vs_sludge[decarbonization_vs_sludge['USD_decarbonization'].notna()]
decarbonization_vs_sludge = decarbonization_vs_sludge[decarbonization_vs_sludge['USD_decarbonization'] <= 0]

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

assert (decarbonization_vs_sludge['total_sludge_amount_kg_per_year']/1000/365).max() <= 250
assert (decarbonization_vs_sludge['CO2_reduction']/30/365/1000).max() <= 50

ax.set_xlim(0, 250)
ax.set_ylim(-100, 200)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Solids}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Management\ CI}$' + '\n[kg CO${_2}$ eq·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(-100, 250, 50))

ax_top = ax.twiny()
ax_top.set_xlim(0, 250)
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(-100, 250, 50))

ax_right = ax.twinx()
ax_right.set_ylim(-100, 200)
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 300, 50))
plt.yticks(np.arange(-100, 250, 50))

ax.scatter(x=decarbonization_vs_sludge['total_sludge_amount_kg_per_year']/1000/365,
           y=decarbonization_vs_sludge['CI'],
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% cost vs distance

# !!! update the file here if necessary
saving_vs_distance = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
saving_vs_distance = saving_vs_distance[saving_vs_distance['USD_decarbonization'].notna()]
saving_vs_distance = saving_vs_distance[saving_vs_distance['USD_decarbonization'] <= 0]

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

assert (saving_vs_distance['real_distance_km']).max() <= 1400
assert (saving_vs_sludge['saving']/30/365/1000).max() <= 180

ax.set_xlim(0, 1400)
ax.set_ylim(-200, 700)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Management\ cost}$' + '\n[\$·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 1600, 200))
plt.yticks(np.arange(-200, 800, 100))

ax_top = ax.twiny()
ax_top.set_xlim(0, 1400)
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 1600, 200))
plt.yticks(np.arange(-200, 800, 100))

ax_right = ax.twinx()
ax_right.set_ylim(-200, 700)
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1600, 200))
plt.yticks(np.arange(-200, 800, 100))

ax.scatter(x=saving_vs_distance['real_distance_km'],
           y=saving_vs_distance['cost'],
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% CI vs distance

# !!! update the file here if necessary
decarbonization_vs_distance = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
decarbonization_vs_distance = decarbonization_vs_distance[decarbonization_vs_distance['USD_decarbonization'].notna()]
decarbonization_vs_distance = decarbonization_vs_distance[decarbonization_vs_distance['USD_decarbonization'] <= 0]

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

assert (decarbonization_vs_distance['real_distance_km']).max() <= 1400
assert (decarbonization_vs_sludge['CO2_reduction']/30/365/1000).max() <= 50

ax.set_xlim(0, 1400)
ax.set_ylim(-100, 200)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Management\ CI}$' + '\n[kg CO${_2}$ eq·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 1600, 200))
plt.yticks(np.arange(-100, 250, 50))

ax_top = ax.twiny()
ax_top.set_xlim(0, 1400)
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 1600, 200))
plt.yticks(np.arange(-100, 250, 50))

ax_right = ax.twinx()
ax_right.set_ylim(-100, 200)
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 1600, 200))
plt.yticks(np.arange(-100, 250, 50))

ax.scatter(x=decarbonization_vs_distance['real_distance_km'],
           y=decarbonization_vs_distance['CI'],
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% cost vs degestion or not

# !!! update the file here if necessary
digestion_or_not = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
digestion_or_not = digestion_or_not[digestion_or_not['USD_decarbonization'].notna()]
digestion_or_not = digestion_or_not[digestion_or_not['USD_decarbonization'] <= 0]

digestion = digestion_or_not.loc[(digestion_or_not['sludge_anaerobic_digestion'] == 1) | (digestion_or_not['sludge_aerobic_digestion'] == 1), 'cost']
no_digestion = digestion_or_not.loc[(digestion_or_not['sludge_anaerobic_digestion'] == 0) & (digestion_or_not['sludge_aerobic_digestion'] == 0), 'cost']

digestion_or_not_plot = pd.DataFrame({'digestion': digestion,
                                      'no_digestion': no_digestion})

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_ylim(-200, 700)

ax.tick_params(direction='inout', length=20, width=3,
               bottom=False, top=False, left=True, right=False, pad=0)

ax.set_ylabel(r'$\mathbf{Management\ cost}$' + '\n[\$·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

ax.set_xticklabels(['digestion','no digestion'])
plt.yticks(np.arange(-200, 800, 100))

ax_right = ax.twinx()
ax_right.tick_params(direction='in', length=10, width=3,
                     bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-200, 800, 100))

bp = ax.boxplot([digestion_or_not_plot['digestion'].dropna(),
                 digestion_or_not_plot['no_digestion'].dropna()],
                showfliers=False, widths=0.7, patch_artist=True)

for box in bp['boxes']:
    box.set(color='none', facecolor='none', linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='none', linewidth=3)

for median in bp['medians']:
    median.set(color='none', linewidth=3)
    
for cap in bp['caps']:
    cap.set(color='none', linewidth=3)

scatter_y = digestion_or_not_plot['digestion'].dropna()
scatter_x = np.random.normal(1, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

scatter_y = digestion_or_not_plot['no_digestion'].dropna()
scatter_x = np.random.normal(2, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% CI vs degestion or not

# !!! update the file here if necessary
digestion_or_not = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
digestion_or_not = digestion_or_not[digestion_or_not['USD_decarbonization'].notna()]
digestion_or_not = digestion_or_not[digestion_or_not['USD_decarbonization'] <= 0]

digestion = digestion_or_not.loc[(digestion_or_not['sludge_anaerobic_digestion'] == 1) | (digestion_or_not['sludge_aerobic_digestion'] == 1), 'CI']
no_digestion = digestion_or_not.loc[(digestion_or_not['sludge_anaerobic_digestion'] == 0) & (digestion_or_not['sludge_aerobic_digestion'] == 0), 'CI']

digestion_or_not_plot = pd.DataFrame({'digestion': digestion,
                                      'no_digestion': no_digestion})

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_ylim(-100, 200)

ax.tick_params(direction='inout', length=20, width=3,
               bottom=False, top=False, left=True, right=False, pad=0)

ax.set_ylabel(r'$\mathbf{Management\ CI}$' + '\n[kg CO${_2}$ eq·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

ax.set_xticklabels(['digestion','no digestion'])
plt.yticks(np.arange(-100, 250, 50))

ax_right = ax.twinx()
ax_right.tick_params(direction='in', length=10, width=3,
                     bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-100, 250, 50))

bp = ax.boxplot([digestion_or_not_plot['digestion'].dropna(),
                 digestion_or_not_plot['no_digestion'].dropna()],
                showfliers=False, widths=0.7, patch_artist=True)

for box in bp['boxes']:
    box.set(color='none', facecolor='none', linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='none', linewidth=3)

for median in bp['medians']:
    median.set(color='none', linewidth=3)
    
for cap in bp['caps']:
    cap.set(color='none', linewidth=3)

scatter_y = digestion_or_not_plot['digestion'].dropna()
scatter_x = np.random.normal(1, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

scatter_y = digestion_or_not_plot['no_digestion'].dropna()
scatter_x = np.random.normal(2, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% cost vs nitrogen fertilizer type

# !!! update the file here if necessary
nitrogen_fertilizer_type = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
nitrogen_fertilizer_type = nitrogen_fertilizer_type[nitrogen_fertilizer_type['USD_decarbonization'].notna()]
nitrogen_fertilizer_type = nitrogen_fertilizer_type[nitrogen_fertilizer_type['USD_decarbonization'] <= 0]

anhydrous_ammonia = nitrogen_fertilizer_type.loc[nitrogen_fertilizer_type['nitrogen_fertilizer'] == 'NH3', 'cost']
urea = nitrogen_fertilizer_type.loc[nitrogen_fertilizer_type['nitrogen_fertilizer'] == 'urea', 'cost']
UAN = nitrogen_fertilizer_type.loc[nitrogen_fertilizer_type['nitrogen_fertilizer'] == 'UAN', 'cost']

nitrogen_fertilizer_type_plot = pd.DataFrame({'anhydrous_ammonia': anhydrous_ammonia,
                                              'urea': urea,
                                              'UAN': UAN})

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_ylim(-200, 700)

ax.tick_params(direction='inout', length=20, width=3,
               bottom=False, top=False, left=True, right=False, pad=0)

ax.set_ylabel(r'$\mathbf{Management\ cost}$' + '\n[\$·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

ax.set_xticklabels(['anhydrous\nammonia','urea','UAN'], linespacing=0.8)
plt.yticks(np.arange(-200, 800, 100))

ax_right = ax.twinx()
ax_right.tick_params(direction='in', length=10, width=3,
                     bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-200, 800, 100))

bp = ax.boxplot([nitrogen_fertilizer_type_plot['anhydrous_ammonia'].dropna(),
                 nitrogen_fertilizer_type_plot['urea'].dropna(),
                 nitrogen_fertilizer_type_plot['UAN'].dropna()],
                showfliers=False, widths=0.7, patch_artist=True)

for box in bp['boxes']:
    box.set(color='none', facecolor='none', linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='none', linewidth=3)

for median in bp['medians']:
    median.set(color='none', linewidth=3)
    
for cap in bp['caps']:
    cap.set(color='none', linewidth=3)

scatter_y = nitrogen_fertilizer_type_plot['anhydrous_ammonia'].dropna()
scatter_x = np.random.normal(1, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

scatter_y = nitrogen_fertilizer_type_plot['urea'].dropna()
scatter_x = np.random.normal(2, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

scatter_y = nitrogen_fertilizer_type_plot['UAN'].dropna()
scatter_x = np.random.normal(3, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% CI vs nitrogen fertilizer type

# !!! update the file here if necessary
nitrogen_fertilizer_type = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
nitrogen_fertilizer_type = nitrogen_fertilizer_type[nitrogen_fertilizer_type['USD_decarbonization'].notna()]
nitrogen_fertilizer_type = nitrogen_fertilizer_type[nitrogen_fertilizer_type['USD_decarbonization'] <= 0]

anhydrous_ammonia = nitrogen_fertilizer_type.loc[nitrogen_fertilizer_type['nitrogen_fertilizer'] == 'NH3', 'CI']
urea = nitrogen_fertilizer_type.loc[nitrogen_fertilizer_type['nitrogen_fertilizer'] == 'urea', 'CI']
UAN = nitrogen_fertilizer_type.loc[nitrogen_fertilizer_type['nitrogen_fertilizer'] == 'UAN', 'CI']

nitrogen_fertilizer_type_plot = pd.DataFrame({'anhydrous_ammonia': anhydrous_ammonia,
                                              'urea': urea,
                                              'UAN': UAN})

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_ylim(-100, 200)

ax.tick_params(direction='inout', length=20, width=3,
               bottom=False, top=False, left=True, right=False, pad=0)

ax.set_ylabel(r'$\mathbf{Management\ CI}$' + '\n[kg CO${_2}$ eq·tonne${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

ax.set_xticklabels(['anhydrous\nammonia','urea','UAN'], linespacing=0.8)
plt.yticks(np.arange(-100, 250, 50))

ax_right = ax.twinx()
ax_right.tick_params(direction='in', length=10, width=3,
                     bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-100, 250, 50))

bp = ax.boxplot([nitrogen_fertilizer_type_plot['anhydrous_ammonia'].dropna(),
                 nitrogen_fertilizer_type_plot['urea'].dropna(),
                 nitrogen_fertilizer_type_plot['UAN'].dropna()],
                showfliers=False, widths=0.7, patch_artist=True)

for box in bp['boxes']:
    box.set(color='none', facecolor='none', linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='none', linewidth=3)

for median in bp['medians']:
    median.set(color='none', linewidth=3)
    
for cap in bp['caps']:
    cap.set(color='none', linewidth=3)

scatter_y = nitrogen_fertilizer_type_plot['anhydrous_ammonia'].dropna()
scatter_x = np.random.normal(1, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

scatter_y = nitrogen_fertilizer_type_plot['urea'].dropna()
scatter_x = np.random.normal(2, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

scatter_y = nitrogen_fertilizer_type_plot['UAN'].dropna()
scatter_x = np.random.normal(3, 0.12, size=len(scatter_y))
ax.scatter(scatter_x, scatter_y,
           s=300,
           c=a,
           linewidths=2,
           edgecolors='k')

#%% relative importance of inter-WRRF contextual parameters to cost

# !!! update the file here if necessary
relative_importance = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
relative_importance = relative_importance[relative_importance['USD_decarbonization'].notna()]
relative_importance = relative_importance[relative_importance['USD_decarbonization'] <= 0]
relative_importance['digestion'] = relative_importance[['sludge_aerobic_digestion','sludge_anaerobic_digestion']].values.max(axis=1)
relative_importance.loc[relative_importance['digestion'] == 1, 'digestion'] = 'Y'
relative_importance.loc[relative_importance['digestion'] == 0, 'digestion'] = 'N'

X_ML = relative_importance[['total_sludge_amount_kg_per_year','real_distance_km','digestion','nitrogen_fertilizer']]
y_ML = relative_importance['cost']

num_features = ['total_sludge_amount_kg_per_year','real_distance_km']
cat_features = ['digestion','nitrogen_fertilizer']

preprocessor = ColumnTransformer([('num', StandardScaler(), num_features),
                                  ('cat', OneHotEncoder(), cat_features)])

# create a RandomForest pipeline
pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                           ('regressor', RandomForestRegressor(n_estimators=100))])

# split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_ML, y_ML, test_size=0.2, random_state=42)

# fit the model
pipeline.fit(X_train, y_train)

# get feature importances
importances = pipeline.named_steps['regressor'].feature_importances_

numerical_columns = num_features
one_hot_columns = pipeline.named_steps['preprocessor'].transformers_[1][1].get_feature_names_out(cat_features)
feature_names =  numerical_columns + list(one_hot_columns)

feature_importance = pd.DataFrame({'Feature': feature_names,
                                   'Importance': importances})

print(feature_importance)

#%% relative importance of inter-WRRF contextual parameters to CI

# !!! update the file here if necessary
relative_importance = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
relative_importance = relative_importance[relative_importance['USD_decarbonization'].notna()]
relative_importance = relative_importance[relative_importance['USD_decarbonization'] <= 0]
relative_importance['digestion'] = relative_importance[['sludge_aerobic_digestion','sludge_anaerobic_digestion']].values.max(axis=1)
relative_importance.loc[relative_importance['digestion'] == 1, 'digestion'] = 'Y'
relative_importance.loc[relative_importance['digestion'] == 0, 'digestion'] = 'N'

X_ML = relative_importance[['total_sludge_amount_kg_per_year','real_distance_km','digestion','nitrogen_fertilizer']]
y_ML = relative_importance['CI']

num_features = ['total_sludge_amount_kg_per_year','real_distance_km']
cat_features = ['digestion','nitrogen_fertilizer']

preprocessor = ColumnTransformer([('num', StandardScaler(), num_features),
                                  ('cat', OneHotEncoder(), cat_features)])

# create a RandomForest pipeline
pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                           ('regressor', RandomForestRegressor(n_estimators=100))])

# split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_ML, y_ML, test_size=0.2, random_state=42)

# fit the model
pipeline.fit(X_train, y_train)

# get feature importances
importances = pipeline.named_steps['regressor'].feature_importances_

numerical_columns = num_features
one_hot_columns = pipeline.named_steps['preprocessor'].transformers_[1][1].get_feature_names_out(cat_features)
feature_names =  numerical_columns + list(one_hot_columns)

feature_importance = pd.DataFrame({'Feature': feature_names,
                                   'Importance': importances})

print(feature_importance)

#%% test required Spearman sample number

def spearman_sample_size(rho, alpha, power):
    fisher_z = 0.5*np.log((1+rho)/(1-rho))
    # two-tailed test
    z_alpha = stats.norm.ppf(1-alpha/2)
    z_beta = stats.norm.ppf(power)
    n = ((z_alpha+z_beta)/fisher_z)**2+3
    return int(np.ceil(n))

# demonstrate 1000 samples are enough when power = 0.8 and rho = 0.2 (weak correlation)
# adjusted alpha for all tests using Bonferroni correction
# may have 50-100 parameters (based on geospatial_models.py)
# may have 1-15 metrics (based on geospatial_models.py)
for parameter_number in range(50, 101):
    for metric_number in range(1, 16):
        alpha_adjusted = 0.05/parameter_number/metric_number
        assert spearman_sample_size(0.2, alpha_adjusted, 0.8) < 1000

#%% qualified facility level uncertainty and sensitivity analyses

filterwarnings('ignore')

# !!! update the file here if necessary
qualified_facility = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
qualified_facility = qualified_facility[qualified_facility['USD_decarbonization'].notna()]
qualified_facility = qualified_facility[qualified_facility['USD_decarbonization'] <= 0]

print(len(qualified_facility))

# $/tonne
qualified_facility['waste_cost'] = sum(qualified_facility[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/qualified_facility['total_sludge_amount_kg_per_year']*1000
# kg CO2 eq/tonne
qualified_facility['waste_GHG'] =  sum(qualified_facility[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/qualified_facility['total_sludge_amount_kg_per_year']*1000

geo_uncertainty_cost = pd.DataFrame()
geo_uncertainty_CI = pd.DataFrame()
geo_uncertainty_saving = pd.DataFrame()
geo_uncertainty_decarbonization = pd.DataFrame()
geo_uncertainty_sludge_N = pd.DataFrame()
geo_uncertainty_sludge_P = pd.DataFrame()
geo_uncertainty_HNO3_N = pd.DataFrame()
geo_uncertainty_biocrude = pd.DataFrame()
geo_uncertainty_DAP = pd.DataFrame()
geo_uncertainty_anhydrous_ammonia = pd.DataFrame()
geo_uncertainty_urea = pd.DataFrame()
geo_uncertainty_UAN = pd.DataFrame()
cost_spearman_r = pd.DataFrame()
cost_spearman_p = pd.DataFrame()
CI_spearman_r = pd.DataFrame()
CI_spearman_p = pd.DataFrame()

# !!! run in different consoles to speed up: 0, 80, 160, 240, 320, 400, 480, len(qualified_facility)
for i in range(0, len(qualified_facility)):
    sys = create_geospatial_system(size=qualified_facility.iloc[i]['total_sludge_amount_kg_per_year']/1000/365,
                                   sludge_transportation=0,
                                   sludge_distance=100,
                                   biocrude_distance=qualified_facility.iloc[i]['real_distance_km'],
                                   anaerobic_digestion=qualified_facility.iloc[i]['sludge_anaerobic_digestion'],
                                   aerobic_digestion=qualified_facility.iloc[i]['sludge_aerobic_digestion'],
                                   ww_2_dry_sludge_ratio=1,
                                   state=qualified_facility.iloc[i]['state'],
                                   nitrogen_fertilizer=qualified_facility.iloc[i]['nitrogen_fertilizer'],
                                   elec_GHG=qualified_facility.iloc[i]['kg_CO2e_kWh'],
                                   wage_adjustment=qualified_facility.iloc[i]['wage_quotient']/100)
    
    # if AD and AeD both exist, assume AeD (due to higher ash content) to be conservative
    if qualified_facility.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_values=[0.349, 0.436, 0.523, 'aerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612]
    elif qualified_facility.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612] 
    else:
        sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
        sludge_lipid_values=[0.080, 0.206, 0.308]
        sludge_protein_values=[0.380, 0.456, 0.485]
    
    crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']==qualified_facility.iloc[i]['state']]['2022_min'].iloc[0],
                              crude_oil_price_data[crude_oil_price_data['state']==qualified_facility.iloc[i]['state']]['2022_average'].iloc[0],
                              crude_oil_price_data[crude_oil_price_data['state']==qualified_facility.iloc[i]['state']]['2022_max'].iloc[0]]
    
    try:
        DAP_price_values = fertilizer_price_uncertainty[qualified_facility.iloc[i]['state']]['DAP']
    except KeyError:
        DAP_price_values = [588, 984.5, 1381]
    
    try:
        anhydrous_ammonia_price_values = fertilizer_price_uncertainty[qualified_facility.iloc[i]['state']]['anhydrous_ammonia']
    except KeyError:
        anhydrous_ammonia_price_values = [1065, 1407.5, 1750]
    
    try:
        urea_price_values = fertilizer_price_uncertainty[qualified_facility.iloc[i]['state']]['urea']
    except KeyError:
        urea_price_values = [447, 855, 1263]
    
    try:
        UAN_price_values = fertilizer_price_uncertainty[qualified_facility.iloc[i]['state']]['UAN']
    except KeyError:
        UAN_price_values = [350, 605, 860]
    
    model = create_geospatial_model(system=sys,
                                    sludge_ash=sludge_ash_values,
                                    sludge_lipid=sludge_lipid_values,
                                    sludge_protein=sludge_protein_values,
                                    crude_oil_price=crude_oil_price_values,
                                    DAP_price=DAP_price_values,
                                    anhydrous_ammonia_price=anhydrous_ammonia_price_values,
                                    urea_price=urea_price_values,
                                    UAN_price=UAN_price_values)
    
    kwargs = {'N':1000,'rule':'L','seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    # note three N fertilizers have different numbers of parameters, need to separate
    r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', nan_policy='omit')
    cost_r_list = [i for i in r_df['Geospatial','Sludge management cost [$/tonne dry sludge]']]
    cost_p_list = [i for i in p_df['Geospatial','Sludge management cost [$/tonne dry sludge]']]
    CI_r_list = [i for i in r_df['Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]']]
    CI_p_list = [i for i in p_df['Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]']]
    
    for data in [cost_r_list, cost_p_list, CI_r_list, CI_p_list]:
        assert len(data) <= 100
        nan_to_add = 100 - len(data)
        data.extend([np.nan]*nan_to_add)
    
    # $/tonne
    geo_uncertainty_cost[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Sludge management cost [$/tonne dry sludge]')]
    # kg CO2 eq/tonne
    geo_uncertainty_CI[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')]
    # $/day
    geo_uncertainty_saving[qualified_facility.iloc[i]['CWNS']] = (qualified_facility.iloc[i]['waste_cost'] -\
                                                                  results[('Geospatial','Sludge management cost [$/tonne dry sludge]')])*\
                                                                 qualified_facility.iloc[i]['total_sludge_amount_kg_per_year']/1000/365
    # kg CO2 eq/day
    geo_uncertainty_decarbonization[qualified_facility.iloc[i]['CWNS']] = (qualified_facility.iloc[i]['waste_GHG'] -\
                                                                             results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')])*\
                                                                            qualified_facility.iloc[i]['total_sludge_amount_kg_per_year']/1000/365
    # tonne/year
    geo_uncertainty_sludge_N[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Sludge N [tonne/year]')]
    # tonne/year
    geo_uncertainty_HNO3_N[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','HNO3 N [tonne/year]')]
    # tonne/year
    geo_uncertainty_sludge_P[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Sludge P [tonne/year]')]
    # BPD
    geo_uncertainty_biocrude[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Biocrude production [BPD]')]
    # tonne/year
    geo_uncertainty_DAP[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','DAP production [tonne/year]')]
    # tonne/year
    geo_uncertainty_anhydrous_ammonia[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Anhydrous ammonia production [tonne/year]')]
    # tonne/year
    geo_uncertainty_urea[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','Urea production [tonne/year]')]
    # tonne/year
    geo_uncertainty_UAN[qualified_facility.iloc[i]['CWNS']] = results[('Geospatial','UAN production [tonne/year]')]
    
    cost_spearman_r[qualified_facility.iloc[i]['CWNS']] = cost_r_list
    cost_spearman_p[qualified_facility.iloc[i]['CWNS']] = cost_p_list
    CI_spearman_r[qualified_facility.iloc[i]['CWNS']] = CI_r_list
    CI_spearman_p[qualified_facility.iloc[i]['CWNS']] = CI_p_list
    
    # check progress
    if i%5 == 0:
        print(i)

geo_uncertainty_cost.to_excel(folder + f'results/qualified_facility/before_integration/cost_dollar_per_tonne_{date.today()}_{i}.xlsx')
geo_uncertainty_CI.to_excel(folder + f'results/qualified_facility/before_integration/CI_kg_CO2_eq_per_day_{date.today()}_{i}.xlsx')
geo_uncertainty_saving.to_excel(folder + f'results/qualified_facility/before_integration/saving_dollar_per_day_{date.today()}_{i}.xlsx')
geo_uncertainty_decarbonization.to_excel(folder + f'results/qualified_facility/before_integration/decarbonization_kg_CO2_eq_per_day_{date.today()}_{i}.xlsx')
geo_uncertainty_sludge_N.to_excel(folder + f'results/qualified_facility/before_integration/sludge_N_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_HNO3_N.to_excel(folder + f'results/qualified_facility/before_integration/HNO3_N_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_sludge_P.to_excel(folder + f'results/qualified_facility/before_integration/sludge_P_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_biocrude.to_excel(folder + f'results/qualified_facility/before_integration/biocrude_BPD_{date.today()}_{i}.xlsx')
geo_uncertainty_DAP.to_excel(folder + f'results/qualified_facility/before_integration/DAP_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_anhydrous_ammonia.to_excel(folder + f'results/qualified_facility/before_integration/anhydrous_ammonia_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_urea.to_excel(folder + f'results/qualified_facility/before_integration/urea_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_UAN.to_excel(folder + f'results/qualified_facility/before_integration/UAN_tonne_per_year_{date.today()}_{i}.xlsx')

cost_spearman_r.to_excel(folder + f'results/qualified_facility/before_integration/cost_spearman_r_{date.today()}_{i}.xlsx')
cost_spearman_p.to_excel(folder + f'results/qualified_facility/before_integration/cost_spearman_p_{date.today()}_{i}.xlsx')
CI_spearman_r.to_excel(folder + f'results/qualified_facility/before_integration/CI_spearman_r_{date.today()}_{i}.xlsx')
CI_spearman_p.to_excel(folder + f'results/qualified_facility/before_integration/CI_spearman_p_{date.today()}_{i}.xlsx')

#%% qualified facility uncertainty - integration

# !!! update these files if necessary
def combine_file(name):
    file_1 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_79.xlsx')
    file_2 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_159.xlsx')
    file_3 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_239.xlsx')
    file_4 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_319.xlsx')
    file_5 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_399.xlsx')
    file_6 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_479.xlsx')
    file_7 = pd.read_excel(folder + f'results/qualified_facility/before_integration/{name}_2025-02-16_554.xlsx')

    for file in [file_1, file_2, file_3, file_4, file_5, file_6, file_7]:
        file.drop('Unnamed: 0', axis=1, inplace=True)

    integrated_file = pd.concat([file_1, file_2, file_3, file_4, file_5, file_6, file_7], axis=1)
    
    integrated_file.to_excel(folder + f'results/qualified_facility/integrated_{name}_{date.today()}.xlsx')

for name in ['cost_dollar_per_tonne','CI_kg_CO2_eq_per_day','saving_dollar_per_day',
             'decarbonization_kg_CO2_eq_per_day','sludge_N_tonne_per_year',
             'HNO3_N_tonne_per_year','sludge_P_tonne_per_year','biocrude_BPD','DAP_tonne_per_year',
             'anhydrous_ammonia_tonne_per_year','urea_tonne_per_year','UAN_tonne_per_year']:
    combine_file(name)

#%% biocrude transportation Chord diagram

# import only if needed
from d3blocks import D3Blocks

biocrude_transportation = pd.read_excel(folder + 'results/qualified_facility/integrated_biocrude_BPD_2025-02-18.xlsx')
biocrude_transportation.drop('Unnamed: 0', axis=1, inplace=True)
biocrude_transportation = biocrude_transportation.transpose()
biocrude_transportation['5th'] = biocrude_transportation[0:999].quantile(q=0.05, axis=1)
biocrude_transportation['50th'] = biocrude_transportation[0:999].quantile(q=0.5, axis=1)
biocrude_transportation['95th'] = biocrude_transportation[0:999].quantile(q=0.95, axis=1)
biocrude_transportation.reset_index(names='CWNS', inplace=True)
biocrude_transportation = biocrude_transportation[['CWNS','5th','50th','95th']]

WRRF_all = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')

biocrude_transportation = biocrude_transportation.merge(WRRF_all, how='left', on='CWNS')

biocrude_transportation['source'] = biocrude_transportation['state']
biocrude_transportation['target'] = biocrude_transportation['State'].apply(lambda x: state_ID[x])
biocrude_transportation['weight'] = biocrude_transportation['50th']
biocrude_transportation = biocrude_transportation[['source','target','weight']]

d3 = D3Blocks()

# ordering = sorted(state_PADD.items(), key=lambda x: x[1])   
# ordering = [state_ID[i[0]] for i in ordering]

ordering = ['CT','DE','DC','FL','ME','GA','RI','NJ','NH','NY','VT','PA','MD',
            'MA','NC','SC','VA','WV','IL','IN','IA','KS','KY','NE','MI','MN',
            'MO','OH','ND','OK','SD','TN','WI','AR','AL','MS','TX','NM','LA',
            'ID','CO','MT','UT','WY','AZ','NV','CA','OR','WA']

d3.chord(biocrude_transportation,
         ordering=ordering,
         color='source',
         figsize=[1200,1200],
         fontsize=32,
         arrowhead=10,
         filepath=folder+f'results/interstate_biocrude_transportation_{date.today()}.html',
         save_button=False)

d3.node_properties['color'] = d3.node_properties['label'].apply(lambda x: PADD_color[state_PADD[ID_state[x]]])
d3.edge_properties['color'] = d3.edge_properties['source'].apply(lambda x: PADD_color[state_PADD[ID_state[x]]])

d3.show()

#%% N and P offsets - data preparation

def get_copula_sum(item):
    if type(item) == str:
        production = pd.read_excel(folder + f'results/qualified_facility/{item}.xlsx')
        production.drop('Unnamed: 0', axis=1, inplace=True)
    else:
        production = item
    
    # generate samples
    np.random.seed(3221)
    
    # combine datasets
    data_samples = np.column_stack([np.array([production[i]]).reshape(len(production),) for i in production.columns])
    
    # apply Gaussian Copula
    # generate independent normal samples
    normal_samples = np.random.randn(len(data_samples), len(production.columns))
    
    # define correlation matrix for copula
    rho = 0.8
    cov_matrix = np.full((len(production.columns), len(production.columns)), rho)
    np.fill_diagonal(cov_matrix, 1)
    
    # Cholesky decomposition
    L = cholesky(cov_matrix, lower=True)
    dependent_samples = normal_samples @ L.T
    
    # map back to dependent uniform [0, 1] space
    dependent_uniform = stats.norm.cdf(dependent_samples)
    
    # transform back to original distributions using quantiles
    sum_results = sum(np.quantile(np.array([production[production.columns[i]]]).reshape(len(production),), dependent_uniform[:, i]) for i in range(len(production.columns)))
    
    return sum_results

# !!! update the file here if necessary
P_offset = get_copula_sum('integrated_DAP_tonne_per_year_2025-02-18')/132.06*30.973762/P['total'].sum()*1000*100
# !!! update the file here if necessary
N_offset = (get_copula_sum('integrated_anhydrous_ammonia_tonne_per_year_2025-02-18')/17.031*14.0067 +\
            get_copula_sum('integrated_urea_tonne_per_year_2025-02-18')/60.06*2*14.0067 +\
            get_copula_sum('integrated_UAN_tonne_per_year_2025-02-18')*0.3)/N['total'].sum()*1000*100

all_facility = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
selected_facility = all_facility[all_facility['USD_decarbonization'].notna()]
selected_facility = selected_facility[selected_facility['USD_decarbonization'] <= 0]

# estimate 83% coverage
P_ratio = selected_facility['total_sludge_amount_kg_per_year'].sum()/all_facility['total_sludge_amount_kg_per_year'].sum()/0.83

# estimate 83% coverage
def N_ratio(fertilizer):
    selected_N = selected_facility[selected_facility['nitrogen_fertilizer'] == fertilizer]['total_sludge_amount_kg_per_year'].sum()
    all_N = all_facility[all_facility['nitrogen_fertilizer'] == fertilizer]['total_sludge_amount_kg_per_year'].sum()
    return selected_N/all_N/0.83

P_offset_future = P_offset/P_ratio
# !!! update the file here if necessary
N_offset_future = (get_copula_sum('integrated_anhydrous_ammonia_tonne_per_year_2025-02-18')/N_ratio('NH3')/17.031*14.0067 +\
                   get_copula_sum('integrated_urea_tonne_per_year_2025-02-18')/N_ratio('urea')/60.06*2*14.0067 +\
                   get_copula_sum('integrated_UAN_tonne_per_year_2025-02-18')/N_ratio('UAN')*0.3)/N['total'].sum()*1000*100

#%% N and P offsets - visualization

offset = pd.DataFrame({'N_near_term': pd.Series(N_offset),
                       'N_future': pd.Series(N_offset_future),
                       'P_near_term': pd.Series(P_offset),
                       'P_future': pd.Series(P_offset_future)})

fig = plt.figure(figsize=(12.5, 10))

gs = fig.add_gridspec(1, 5, hspace=0, wspace=0)

def add_region(position, color):  
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['hatch.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    ax = fig.add_subplot(gs[0, position])
    
    ax = plt.gca()
    ax.set_ylim(0, 6)
    
    if position == 0:
        ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)
        ax.set_ylabel(r'$\mathbf{Offset}$ [%]', fontname='Arial', fontsize=35)
        ax.spines[['right']].set_visible(False)
    
    elif position == 3:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, length=0, labelcolor='none')
        
        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
        ax.spines[['left']].set_visible(False)
        ax_right.spines[['left']].set_visible(False)
    
    else:
        ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
        ax.spines[['left','right']].set_visible(False)
    
    bp = ax.boxplot(offset.iloc[:,position].dropna(), showfliers=False, widths=0.75, patch_artist=True)
    
    if position == 0 or position == 2:
        for box in bp['boxes']:
            box.set(color='k', facecolor=color, linewidth=3)
        
    if position == 1 or position == 3:
        for box in bp['boxes']:
            box.set(color='k', facecolor=color, linewidth=3, hatch = '//',  edgecolor='w')

    for whisker in bp['whiskers']:
        whisker.set(color='k', linewidth=3)

    for median in bp['medians']:
        median.set(color='k', linewidth=3)
        
    for cap in bp['caps']:
        cap.set(color='k', linewidth=3)
        
    if position == 1 or position == 3:  
        bp_2 = ax.boxplot(offset.iloc[:,position].dropna(), showfliers=False, widths=0.75, patch_artist=True)    
        for box in bp_2['boxes']:
            box.set(color='k', facecolor='none', linewidth=3)

        for median in bp_2['medians']:
            median.set(color='k', linewidth=3)
        
    ax.scatter(x=1,
               y=offset.iloc[:,position].dropna().mean(),
               marker='*',
               s=600,
               c='w',
               linewidths=3,
               alpha=1,
               edgecolor='k',
               zorder=3)
    
    # for flier in bp['fliers']:
    #     flier.set(marker='o', markersize=7, markerfacecolor=color, markeredgewidth=1.5)

add_region(0, b)
add_region(1, b)
add_region(2, r)
add_region(3, r)

#%% sampled facility level uncertainty and sensitivity analyses (data preparation)

# !!! update the file here if necessary
sampled_facility = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')

print(len(sampled_facility))

sampled_facility['dollar_per_kWh'] = sampled_facility['state'].apply(lambda x: elec_price[elec_price['state'] == x]['price'].iloc[0]/100)
sampled_facility['crude_oil_dollar_per_barrel'] = sampled_facility['state'].apply(lambda x: crude_oil_price_data[crude_oil_price_data['state'] == x]['2022_average'].iloc[0])

sampled_facility['DAP_price'] = sampled_facility['state'].apply(lambda x: fertilizer_price_uncertainty[x]['DAP'][1] if x in ['AL','IA','IL','NC','OK','SC'] else 984.5)
sampled_facility['anhydrous_ammonia_price'] = sampled_facility['state'].apply(lambda x: fertilizer_price_uncertainty[x]['anhydrous_ammonia'][1] if x in ['IA','IL','OK'] else 1407.5)
sampled_facility['urea_price'] = sampled_facility['state'].apply(lambda x: fertilizer_price_uncertainty[x]['urea'][1] if x in ['AL','IA','IL','NC','OK','SC'] else 855)
sampled_facility['UAN_price'] = sampled_facility['state'].apply(lambda x: fertilizer_price_uncertainty[x]['UAN'][1] if x in ['AL','IA','IL','NC','OK','SC'] else 605)

sampled_facility.loc[sampled_facility['nitrogen_fertilizer']=='NH3', 'nitrogen_fertilizer'] = 1
sampled_facility.loc[sampled_facility['nitrogen_fertilizer']=='urea', 'nitrogen_fertilizer'] = 2
sampled_facility.loc[sampled_facility['nitrogen_fertilizer']=='UAN', 'nitrogen_fertilizer'] = 3

# exclude income_tax since it is decided by the included parameters (even the location can be represented by dollar_per_kWh)
sampled_facility = sampled_facility[['total_sludge_amount_kg_per_year',
                                     'sludge_aerobic_digestion',
                                     'sludge_anaerobic_digestion',
                                     'dollar_per_kWh',
                                     'kg_CO2e_kWh',
                                     'crude_oil_dollar_per_barrel',
                                     'DAP_price',
                                     'anhydrous_ammonia_price',
                                     'urea_price',
                                     'UAN_price',
                                     'wage_quotient',
                                     'nitrogen_fertilizer',
                                     'real_distance_km']]

# scale the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(sampled_facility)

# 1000 clusters
k = 1000

# initialize the AgglomerativeClustering model and fit it to the data
agg_cluster = AgglomerativeClustering(n_clusters=k)
agg_cluster.fit(X_scaled)

# get the cluster labels
labels = agg_cluster.labels_

sampled_facility['label'] = labels

sampled_facility = sampled_facility.groupby('label').median()

sampled_facility.reset_index(inplace=True)

sampled_facility['nitrogen_fertilizer'] = sampled_facility['nitrogen_fertilizer'].apply(lambda x: round(x))

sampled_facility.loc[sampled_facility['nitrogen_fertilizer']==1, 'nitrogen_fertilizer'] = 'NH3'
sampled_facility.loc[sampled_facility['nitrogen_fertilizer']==2, 'nitrogen_fertilizer'] = 'urea'
sampled_facility.loc[sampled_facility['nitrogen_fertilizer']==3, 'nitrogen_fertilizer'] = 'UAN'

sampled_facility.to_excel(folder + f'results/sampled_facility/sampled_facility_{date.today()}.xlsx')

#%% sampled facility level uncertainty and sensitivity analyses - part 1

filterwarnings('ignore')

# !!! update the file here if necessary
sampled_facility = pd.read_excel(folder + 'results/sampled_facility/sampled_facility_2025-02-16.xlsx')

print(len(sampled_facility))

geo_uncertainty_cost = pd.DataFrame()
geo_uncertainty_CI = pd.DataFrame()
geo_uncertainty_saving = pd.DataFrame()
geo_uncertainty_decarbonization = pd.DataFrame()
geo_uncertainty_sludge_N = pd.DataFrame()
geo_uncertainty_sludge_P = pd.DataFrame()
geo_uncertainty_HNO3_N = pd.DataFrame()
geo_uncertainty_biocrude = pd.DataFrame()
geo_uncertainty_DAP = pd.DataFrame()
geo_uncertainty_anhydrous_ammonia = pd.DataFrame()
geo_uncertainty_urea = pd.DataFrame()
geo_uncertainty_UAN = pd.DataFrame()
cost_spearman_r = pd.DataFrame()
cost_spearman_p = pd.DataFrame()
CI_spearman_r = pd.DataFrame()
CI_spearman_p = pd.DataFrame()

# !!! run in different consoles to speed up: 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, len(sampled_facility)
for i in range(0, len(sampled_facility)):
    sys = create_geospatial_system(size=sampled_facility.iloc[i]['total_sludge_amount_kg_per_year']/1000/365,
                                   sludge_transportation=0,
                                   sludge_distance=100,
                                   biocrude_distance=sampled_facility.iloc[i]['real_distance_km'],
                                   anaerobic_digestion=sampled_facility.iloc[i]['sludge_anaerobic_digestion'],
                                   aerobic_digestion=sampled_facility.iloc[i]['sludge_aerobic_digestion'],
                                   ww_2_dry_sludge_ratio=1,
                                   state='US',
                                   nitrogen_fertilizer=sampled_facility.iloc[i]['nitrogen_fertilizer'],
                                   elec_price=sampled_facility.iloc[i]['dollar_per_kWh'],
                                   elec_GHG=sampled_facility.iloc[i]['kg_CO2e_kWh'],
                                   wage_adjustment=sampled_facility.iloc[i]['wage_quotient']/100)
    
    # if AD and AeD both exist, assume AeD (due to higher ash content) to be conservative
    if sampled_facility.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_values=[0.349, 0.436, 0.523, 'aerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612]
    elif sampled_facility.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
        sludge_lipid_values=[0.154, 0.193, 0.232]
        sludge_protein_values=[0.408, 0.510, 0.612] 
    else:
        sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
        sludge_lipid_values=[0.080, 0.206, 0.308]
        sludge_protein_values=[0.380, 0.456, 0.485]
    
    crude_oil_price_values = [sampled_facility.iloc[i]['crude_oil_dollar_per_barrel']*0.8,
                              sampled_facility.iloc[i]['crude_oil_dollar_per_barrel'],
                              sampled_facility.iloc[i]['crude_oil_dollar_per_barrel']*1.2]
    
    DAP_price_values = [sampled_facility.iloc[i]['DAP_price']*0.8,
                        sampled_facility.iloc[i]['DAP_price'],
                        sampled_facility.iloc[i]['DAP_price']*1.2]
    
    anhydrous_ammonia_price_values = [sampled_facility.iloc[i]['anhydrous_ammonia_price']*0.8,
                                      sampled_facility.iloc[i]['anhydrous_ammonia_price'],
                                      sampled_facility.iloc[i]['anhydrous_ammonia_price']*1.2]

    urea_price_values = [sampled_facility.iloc[i]['urea_price']*0.8,
                         sampled_facility.iloc[i]['urea_price'],
                         sampled_facility.iloc[i]['urea_price']*1.2]

    UAN_price_values = [sampled_facility.iloc[i]['UAN_price']*0.8,
                        sampled_facility.iloc[i]['UAN_price'],
                        sampled_facility.iloc[i]['UAN_price']*1.2]
    
    model = create_geospatial_model(system=sys,
                                    sludge_ash=sludge_ash_values,
                                    sludge_lipid=sludge_lipid_values,
                                    sludge_protein=sludge_protein_values,
                                    crude_oil_price=crude_oil_price_values,
                                    DAP_price=DAP_price_values,
                                    anhydrous_ammonia_price=anhydrous_ammonia_price_values,
                                    urea_price=urea_price_values,
                                    UAN_price=UAN_price_values)
    
    kwargs = {'N':1000,'rule':'L','seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    # note three N fertilizers have different numbers of parameters, need to separate
    r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', nan_policy='omit')
    cost_r_list = [i for i in r_df['Geospatial','Sludge management cost [$/tonne dry sludge]']]
    cost_p_list = [i for i in p_df['Geospatial','Sludge management cost [$/tonne dry sludge]']]
    CI_r_list = [i for i in r_df['Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]']]
    CI_p_list = [i for i in p_df['Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]']]
    
    for data in [cost_r_list, cost_p_list, CI_r_list, CI_p_list]:
        assert len(data) <= 100
        nan_to_add = 100 - len(data)
        data.extend([np.nan]*nan_to_add)
    
    # $/tonne
    geo_uncertainty_cost[sampled_facility.iloc[i]['label']] = results[('Geospatial','Sludge management cost [$/tonne dry sludge]')]
    # kg CO2 eq/tonne
    geo_uncertainty_CI[sampled_facility.iloc[i]['label']] = results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')]
    # tonne/year
    geo_uncertainty_sludge_N[sampled_facility.iloc[i]['label']] = results[('Geospatial','Sludge N [tonne/year]')]
    # tonne/year
    geo_uncertainty_HNO3_N[sampled_facility.iloc[i]['label']] = results[('Geospatial','HNO3 N [tonne/year]')]
    # tonne/year
    geo_uncertainty_sludge_P[sampled_facility.iloc[i]['label']] = results[('Geospatial','Sludge P [tonne/year]')]
    # BPD
    geo_uncertainty_biocrude[sampled_facility.iloc[i]['label']] = results[('Geospatial','Biocrude production [BPD]')]
    # tonne/year
    geo_uncertainty_DAP[sampled_facility.iloc[i]['label']] = results[('Geospatial','DAP production [tonne/year]')]
    # tonne/year
    geo_uncertainty_anhydrous_ammonia[sampled_facility.iloc[i]['label']] = results[('Geospatial','Anhydrous ammonia production [tonne/year]')]
    # tonne/year
    geo_uncertainty_urea[sampled_facility.iloc[i]['label']] = results[('Geospatial','Urea production [tonne/year]')]
    # tonne/year
    geo_uncertainty_UAN[sampled_facility.iloc[i]['label']] = results[('Geospatial','UAN production [tonne/year]')]
    
    cost_spearman_r[sampled_facility.iloc[i]['label']] = cost_r_list
    cost_spearman_p[sampled_facility.iloc[i]['label']] = cost_p_list
    CI_spearman_r[sampled_facility.iloc[i]['label']] = CI_r_list
    CI_spearman_p[sampled_facility.iloc[i]['label']] = CI_p_list
    
    # check progress
    if i%5 == 0:
        print(i)

geo_uncertainty_cost.to_excel(folder + f'results/sampled_facility/before_integration/cost_dollar_per_tonne_{date.today()}_{i}.xlsx')
geo_uncertainty_CI.to_excel(folder + f'results/sampled_facility/before_integration/CI_kg_CO2_eq_per_day_{date.today()}_{i}.xlsx')
geo_uncertainty_saving.to_excel(folder + f'results/sampled_facility/before_integration/saving_dollar_per_day_{date.today()}_{i}.xlsx')
geo_uncertainty_decarbonization.to_excel(folder + f'results/sampled_facility/before_integration/decarbonization_kg_CO2_eq_per_day_{date.today()}_{i}.xlsx')
geo_uncertainty_sludge_N.to_excel(folder + f'results/sampled_facility/before_integration/sludge_N_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_HNO3_N.to_excel(folder + f'results/sampled_facility/before_integration/HNO3_N_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_sludge_P.to_excel(folder + f'results/sampled_facility/before_integration/sludge_P_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_biocrude.to_excel(folder + f'results/sampled_facility/before_integration/biocrude_BPD_{date.today()}_{i}.xlsx')
geo_uncertainty_DAP.to_excel(folder + f'results/sampled_facility/before_integration/DAP_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_anhydrous_ammonia.to_excel(folder + f'results/sampled_facility/before_integration/anhydrous_ammonia_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_urea.to_excel(folder + f'results/sampled_facility/before_integration/urea_tonne_per_year_{date.today()}_{i}.xlsx')
geo_uncertainty_UAN.to_excel(folder + f'results/sampled_facility/before_integration/UAN_tonne_per_year_{date.today()}_{i}.xlsx')

cost_spearman_r.to_excel(folder + f'results/sampled_facility/before_integration/cost_spearman_r_{date.today()}_{i}.xlsx')
cost_spearman_p.to_excel(folder + f'results/sampled_facility/before_integration/cost_spearman_p_{date.today()}_{i}.xlsx')
CI_spearman_r.to_excel(folder + f'results/sampled_facility/before_integration/CI_spearman_r_{date.today()}_{i}.xlsx')
CI_spearman_p.to_excel(folder + f'results/sampled_facility/before_integration/CI_spearman_p_{date.today()}_{i}.xlsx')

#%% sampled facility level uncertainty and sensitivity analyses - part 2

# !!! update these files if necessary
def combine_file(name):
    file_1 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_99.xlsx')
    file_2 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_199.xlsx')
    file_3 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_299.xlsx')
    file_4 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_399.xlsx')
    file_5 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_499.xlsx')
    file_6 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_599.xlsx')
    file_7 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_699.xlsx')
    file_8 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_799.xlsx')
    file_9 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_899.xlsx')
    file_10 = pd.read_excel(folder + f'results/sampled_facility/before_integration/{name}_2025-02-16_999.xlsx')

    for file in [file_1, file_2, file_3, file_4, file_5, file_6, file_7, file_8, file_9, file_10]:
        file.drop('Unnamed: 0', axis=1, inplace=True)

    integrated_file = pd.concat([file_1, file_2, file_3, file_4, file_5, file_6, file_7, file_8, file_9, file_10], axis=1)
    
    integrated_file.to_excel(folder + f'results/sampled_facility/integrated_{name}_{date.today()}.xlsx')

for name in ['cost_spearman_r','cost_spearman_p','CI_spearman_r','CI_spearman_p']:
    combine_file(name)

#%% sampled facility level uncertainty and sensitivity analyses - part 3

def plot_sensitivity(data_type):
    data = pd.read_excel(folder + f'results/sampled_facility/integrated_{data_type}_2025-02-17.xlsx')
    data.drop('Unnamed: 0', axis=1, inplace=True)
    
    len_1, len_2, len_3 = sorted(set([len(data[i].dropna()) for i in range(0, data.shape[1])]))
    
    anhydrous_ammonia_data = data[[i for i in range(0, data.shape[1]) if len(data[i].dropna()) == len_1]].dropna().transpose()
    urea_data = data[[i for i in range(0, data.shape[1]) if len(data[i].dropna()) == len_2]].dropna().transpose()
    UAN_data = data[[i for i in range(0, data.shape[1]) if len(data[i].dropna()) == len_3]].dropna().transpose()
    
    def get_parameters(nitrogen_fertilizer):
        sys = create_geospatial_system(size=1, ww_2_dry_sludge_ratio=1, nitrogen_fertilizer=nitrogen_fertilizer)
        model = create_geospatial_model(system=sys,
                                        sludge_ash=[0.1, 0.2, 0.3, 'aerobic_digestion'],
                                        sludge_lipid=[0.1, 0.2, 0.3],
                                        sludge_protein=[0.1, 0.2, 0.3],
                                        crude_oil_price=[1, 2, 3],
                                        DAP_price=[1, 2, 3],
                                        anhydrous_ammonia_price=[1, 2, 3],
                                        urea_price=[1, 2, 3],
                                        UAN_price=[1, 2, 3])
        return model.parameters
    
    anhydrous_ammonia_columns = get_parameters('NH3')
    urea_columns = get_parameters('urea')
    UAN_columns = get_parameters('UAN')
    
    anhydrous_ammonia_data.columns=[anhydrous_ammonia_columns[i].name for i in range(0, len(anhydrous_ammonia_columns))]
    urea_data.columns=[urea_columns[i].name for i in range(0, len(urea_columns))]
    UAN_data.columns=[UAN_columns[i].name for i in range(0, len(UAN_columns))]
    
    anhydrous_ammonia_data.reset_index(inplace=True, names='label')
    urea_data.reset_index(inplace=True, names='label')
    UAN_data.reset_index(inplace=True, names='label')
    
    sampled_facility = pd.read_excel(folder + 'results/sampled_facility/sampled_facility_2025-02-16.xlsx')
    sampled_facility = sampled_facility[['label','total_sludge_amount_kg_per_year']]
    
    anhydrous_ammonia_data = anhydrous_ammonia_data.merge(sampled_facility, how='left', on='label')
    urea_data = urea_data.merge(sampled_facility, how='left', on='label')
    UAN_data = UAN_data.merge(sampled_facility, how='left', on='label')
    
    integrated_data = pd.concat([anhydrous_ammonia_data, urea_data, UAN_data])
    
    for i in integrated_data.columns.drop(['label','total_sludge_amount_kg_per_year']):    
        plt.rcParams['axes.linewidth'] = 2
        plt.rcParams['xtick.labelsize'] = 24
        plt.rcParams['ytick.labelsize'] = 24
    
        plt.rcParams.update({'mathtext.fontset': 'custom'})
        plt.rcParams.update({'mathtext.default': 'regular'})
        plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
        
        fig, ax = plt.subplots(figsize = (5, 2.5))
        
        ax = plt.gca()
        
        low_limit = (integrated_data['total_sludge_amount_kg_per_year']/1000/365).min()
        high_limit = (integrated_data['total_sludge_amount_kg_per_year']/1000/365).max()
        
        ax.set_xlim(low_limit, high_limit)
        if data_type[-1] == 'r':
            ax.set_ylim(-1, 1)
        else:
            ax.set_ylim(0, 1)
        
        plt.xscale('log')
        plt.xticks([0.01, 0.1, 1, 10, 100])
        if data_type[-1] == 'r':
            plt.yticks(np.arange(-1, 1.5, 0.5))
        else:
            plt.yticks(np.arange(0, 1.2, 0.2))
        
        ax.tick_params(direction='inout', length=15, width=2, bottom=True, top=False, left=True, right=False)
        
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        ax_top.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')
        
        plt.xscale('log')
        plt.xticks([0.01, 0.1, 1, 10, 100])
        if data_type[-1] == 'r':
            plt.yticks(np.arange(-1, 1.5, 0.5))
        else:
            plt.yticks(np.arange(0, 1.2, 0.2))

        ax_right = ax.twinx()
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(direction='in', length=7.5, width=2, bottom=False, top=True, left=False, right=True, labelcolor='none')
        
        plt.xscale('log')
        plt.xticks([0.01, 0.1, 1, 10, 100])
        if data_type[-1] == 'r':
            plt.yticks(np.arange(-1, 1.5, 0.5))
        else:
            plt.yticks(np.arange(0, 1.2, 0.2))
        
        ax.set_xlabel(r'$\mathbf{Wastewater\ solids\ quantity}$' + '\n[tonne·day${^{-1}}$]', fontname='Arial', fontsize=24, linespacing=0.8)
        if (data_type[0:4] == 'cost') & (data_type[-1] == 'r'):
            ax.set_ylabel(r"$\mathbf{Spearman's\ \rho}$", fontname='Arial', fontsize=24)
        elif (data_type[0:4] == 'cost') & (data_type[-1] == 'p'):
                ax.set_ylabel(r"$\mathbf{Spearman's\ p}$", fontname='Arial', fontsize=24)
        elif (data_type[0:2] == 'CI') & (data_type[-1] == 'r'):
                ax.set_ylabel(r"$\mathbf{Spearman's\ \rho}$", fontname='Arial', fontsize=24)
        elif (data_type[0:2] == 'CI') & (data_type[-1] == 'p'):
                ax.set_ylabel(r"$\mathbf{Spearman's\ p}$", fontname='Arial', fontsize=24)
        else:
            raise ValueError(f'Check name: {data_type}.')

        mathtext.FontConstantsBase.sup1 = 0.35
    
        if data_type[-1] == 'r':
            color = b
            edgecolor=db
        else:
            color=r
            edgecolor=dr
        
        plt.scatter(integrated_data['total_sludge_amount_kg_per_year']/1000/365, integrated_data[i], color=color, s=50, edgecolors=edgecolor, linewidths=1.5)
        
        if data_type[-1] == 'r':
            plt.plot([0, 250],[-0.2, -0.2], lw=2, color='k', linestyle='--')
            plt.plot([0, 250],[0.2, 0.2], lw=2, color='k', linestyle='--')
        else:
            plt.plot([0, 250],[0.05, 0.05], lw=2, color='k', linestyle='--')
        
        # add _ only
        add_underscore = ['N 2 P','HTLaqueous C slope','TOC TC','P acid recovery ratio','P syn recovery ratio',
                          'HTL TIC factor','CHG TIC factor','CHP unit TIC','IRR','CHG catalyst price',
                          'DAP price','CHG catalyst global warming','H2SO4 global warming','DAP global warming',
                          'CO2 recovery','MEA to CO2','UANSyn product loss ratio','HNO3 price','UAN water price',
                          'UAN price','HNO3 global warming','UAN water global warming','UAN global warming']
        
        # lowercase and add _
        lowercase_add_underscore = ['Sludge dw ash','Sludge afdw lipid','Sludge afdw protein','Lipid 2 C',
                                    'Protein 2 C','Carbo 2 C','Lipid 2 H','Protein 2 H','Carbo 2 H',
                                    'Protein 2 N','Sludge wet density','Sludge distance',
                                    'Enforced heating transfer coefficient','Lipid 2 biocrude',
                                    'Protein 2 biocrude','Carbo 2 biocrude','Protein 2 gas','Carbo 2 gas',
                                    'Biocrude C slope','Biocrude C intercept','Biocrude N slope',
                                    'Biocrude H slope','Biocrude H intercept','Hydrochar C slope',
                                    'Biocrude moisture content','Hydrochar P recovery ratio',
                                    'Crude oil density','Crude oil HHV','Biocrude wet density',
                                    'Biocrude distance','Catalyst lifetime','Gas C 2 total C','Acid vol',
                                    'Electricity price','Wage adjustment','Biocrude price',
                                    'Water steam price','Anhydrous ammonia price','Natural gas price',
                                    'Ash disposal price','Cooling tower makeup water price',
                                    'Cooling tower chemicals price','Sludge transportation price',
                                    'Biocrude transportation price','Electricity global warming',
                                    'Deionized water global warming','Sludge trucking global warming',
                                    'Biocrude trucking global warming','Natural gas global warming',
                                    'Biocrude global warming','Ash disposal global warming',
                                    'Anhydrous ammonia global warming','Makeup MEA price',
                                    'Makeup water price','Urea price','Makeup MEA global warming',
                                    'Makeup water global warming','Urea global warming']
        
        # manually update
        manually_update = {'WHSV':'WHSV',
                           'Influent p H':'influent_pH',
                           '5% h2so4 price':'5%_H2SO4_price',
                           'Na OH price':'NaOH_price',
                           'Na OH global warming':'NaOH_global_warming',
                           'Urea syn nh3 to co2 ratio':'UreaSyn_NH3_to_CO2_ratio',
                           'Urea syn co2 to urea efficiency':'UreaSyn_CO2_to_urea_efficiency',
                           'Urea syn product loss ratio':'UreaSyn_product_loss_ratio',
                           'UANSyn nh3 to co2 ratio':'UANSyn_NH3_to_CO2_ratio',
                           'UANSyn co2 to urea efficiency':'UANSyn_CO2_to_urea_efficiency'}
        
        if i in add_underscore:
            title = i.replace(' ', '_')
        elif i in lowercase_add_underscore:
            title = i.replace(' ', '_')
            title = title[0].lower() + title[1:]
        elif i in manually_update.keys():
            title = manually_update[i]
        else:
            raise ValueError(f'Check parameters name of {i}.')
        
        plt.title(title, fontdict={'font':'Arial','size':24,'fontweight':'bold'}, pad=15, y=1)
        
    return integrated_data

#%% cost_spearman_p visualization

cost_spearman_p = plot_sensitivity('cost_spearman_p')

#%% cost_spearman_r visualization

cost_spearman_r = plot_sensitivity('cost_spearman_r')

#%% select key cost drivers

cost_spearman_p_selected = (cost_spearman_p < 0.05).sum()

cost_spearman_r_selected = (abs(cost_spearman_r) > 0.2).sum()

key_cost_drivers = (cost_spearman_p_selected > 100) & (cost_spearman_r_selected > 100)

cost_mask = key_cost_drivers == True

key_cost_drivers_only = key_cost_drivers[cost_mask]

print(key_cost_drivers_only)

#%% CI_spearman_p visualization

CI_spearman_p = plot_sensitivity('CI_spearman_p')

#%% CI_spearman_r visualization

CI_spearman_r = plot_sensitivity('CI_spearman_r')

#%% select key CI drivers

CI_spearman_p_selected = (CI_spearman_p < 0.05).sum()

CI_spearman_r_selected = (abs(CI_spearman_r) > 0.2).sum()

key_CI_drivers = (CI_spearman_p_selected > 100) & (CI_spearman_r_selected > 100)

CI_mask = key_CI_drivers == True

key_CI_drivers_only = key_CI_drivers[CI_mask]

print(key_CI_drivers_only)

#%% cost vs. IRR

filterwarnings('ignore')

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

print(len(WRRF_input))

total_sludge = WRRF_input['total_sludge_amount_kg_per_year'].sum(axis=0)

total_ash = 0
total_lipid = 0
total_protein = 0
total_sludge = WRRF_input['total_sludge_amount_kg_per_year'].sum()

for i in range(len(WRRF_input)):
    if WRRF_input.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_value=0.436
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
    elif WRRF_input.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_value=0.414
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
    else:
        sludge_ash_value=0.231
        sludge_lipid_value=0.206
        sludge_protein_value=0.456

    total_ash += sludge_ash_value*WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']
    total_lipid += sludge_lipid_value*WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)
    total_protein += sludge_protein_value*WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)

average_ash_dw = total_ash/total_sludge
average_lipid_afdw = total_lipid/(total_sludge-total_ash)
average_protein_afdw = total_protein/(total_sludge-total_ash)

# !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
# set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion']
sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

cost_vs_IRR = pd.DataFrame()

for size in [WRRF_input['total_sludge_amount_kg_per_year'].min()/1000/365, 0.01, 0.1, 1, 10, 100, WRRF_input['total_sludge_amount_kg_per_year'].max()/1000/365]:
    for IRR in [0, 0.01, 0.02, 0.03, 0.04 ,0.05]:
        print('\n\n', f'size: {size} metric tonne/day\n', f'IRR: {IRR}\n')
        
        sys = create_geospatial_system(size=size,
                                       biocrude_distance=(WRRF_input['real_distance_km']*WRRF_input['total_sludge_amount_kg_per_year']).sum()/WRRF_input['total_sludge_amount_kg_per_year'].sum(),
                                       average_sludge_dw_ash=sludge_ash_values[1],
                                       average_sludge_afdw_lipid=sludge_lipid_values[1],
                                       average_sludge_afdw_protein=sludge_protein_values[1],
                                       anaerobic_digestion=None,
                                       aerobic_digestion=None,
                                       ww_2_dry_sludge_ratio=1,
                                       state='US',
                                       nitrogen_fertilizer='UAN',
                                       elec_GHG=(WRRF_input['kg_CO2e_kWh']*WRRF_input['total_sludge_amount_kg_per_year']).sum()/WRRF_input['total_sludge_amount_kg_per_year'].sum(),
                                       wage_adjustment=(WRRF_input.iloc[i]['wage_quotient']/100*WRRF_input['total_sludge_amount_kg_per_year']).sum()/WRRF_input['total_sludge_amount_kg_per_year'].sum())
        
        sys.TEA.IRR = IRR
        
        crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']=='US']['2022_min'].iloc[0],
                                  crude_oil_price_data[crude_oil_price_data['state']=='US']['2022_average'].iloc[0],
                                  crude_oil_price_data[crude_oil_price_data['state']=='US']['2022_max'].iloc[0]]
        
        DAP_price_values = [588, 984.5, 1381]
        anhydrous_ammonia_price_values = [1065, 1407.5, 1750]
        urea_price_values = [447, 855, 1263]
        UAN_price_values = [350, 605, 860]
                
        model = create_geospatial_model(system=sys,
                                        sludge_ash=sludge_ash_values,
                                        sludge_lipid=sludge_lipid_values,
                                        sludge_protein=sludge_protein_values,
                                        crude_oil_price=crude_oil_price_values,
                                        DAP_price=DAP_price_values,
                                        anhydrous_ammonia_price=anhydrous_ammonia_price_values,
                                        urea_price=urea_price_values,
                                        UAN_price=UAN_price_values,
                                        exclude_IRR=True)
        
        kwargs = {'N':1000,'rule':'L','seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        idx = len(model.parameters)
        parameters = model.table.iloc[:, :idx]
        results = model.table.iloc[:, idx:]
        
        if size < 0.01:
            size_name = 'min'
        elif size > 100:
            size_name = 'max'
        else:
            size_name = size
        
        # $/day
        cost_vs_IRR[f'{size_name}_{IRR}'] = results[('Geospatial','Sludge management cost [$/tonne dry sludge]')]
        
cost_vs_IRR.to_excel(folder + f'results/cost_vs_IRR/cost_vs_IRR_{date.today()}.xlsx')

#%% regional uncertainty (build and run model)

# import PADD information
# !!! update the file here if necessary
WRRF_PADD = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
WRRF_PADD = WRRF_PADD[WRRF_PADD['USD_decarbonization'].notna()]
WRRF_PADD = WRRF_PADD[WRRF_PADD['USD_decarbonization'] <= 0]

WRRF_PADD.loc[WRRF_PADD['state'].isin(['CT','DC','DE','FL','GA','MA','MD','ME','NC','NH','NJ','NY','PA','RI','SC','VA','VT','WV']),'WRRF_PADD'] = 1
WRRF_PADD.loc[WRRF_PADD['state'].isin(['IA','IL','IN','KS','KY','MI','MN','MO','ND','NE','OH','OK','SD','TN','WI']),'WRRF_PADD'] = 2
WRRF_PADD.loc[WRRF_PADD['state'].isin(['AL','AR','LA','MS','NM','TX']),'WRRF_PADD'] = 3
WRRF_PADD.loc[WRRF_PADD['state'].isin(['CO','ID','MT','UT','WY']),'WRRF_PADD'] = 4
WRRF_PADD.loc[WRRF_PADD['state'].isin(['AZ','CA','NV','OR','WA']),'WRRF_PADD'] = 5

assert WRRF_PADD.WRRF_PADD.isna().sum() == 0

WRRF_PADD = WRRF_PADD[['CWNS','WRRF_PADD']]

# saving
# !!! update the file here if necessary
regional_saving_uncertainty = pd.read_excel(folder + 'results/qualified_facility/integrated_saving_dollar_per_day_2025-02-18.xlsx')

regional_saving_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)
regional_saving_uncertainty = regional_saving_uncertainty.transpose()
regional_saving_uncertainty.reset_index(inplace=True, names='CWNS')
regional_saving_uncertainty = regional_saving_uncertainty.merge(WRRF_PADD, how='left', on='CWNS')

regional_saving_uncertainty_PADD_1 = regional_saving_uncertainty[regional_saving_uncertainty['WRRF_PADD'] == 1].loc[:, 0:999]
regional_saving_uncertainty_PADD_2 = regional_saving_uncertainty[regional_saving_uncertainty['WRRF_PADD'] == 2].loc[:, 0:999]
regional_saving_uncertainty_PADD_3 = regional_saving_uncertainty[regional_saving_uncertainty['WRRF_PADD'] == 3].loc[:, 0:999]
regional_saving_uncertainty_PADD_4 = regional_saving_uncertainty[regional_saving_uncertainty['WRRF_PADD'] == 4].loc[:, 0:999]
regional_saving_uncertainty_PADD_5 = regional_saving_uncertainty[regional_saving_uncertainty['WRRF_PADD'] == 5].loc[:, 0:999]

regional_saving_result = pd.DataFrame({'PADD_1': get_copula_sum(regional_saving_uncertainty_PADD_1.transpose()),
                                       'PADD_2': get_copula_sum(regional_saving_uncertainty_PADD_2.transpose()),
                                       'PADD_3': get_copula_sum(regional_saving_uncertainty_PADD_3.transpose()),
                                       'PADD_4': get_copula_sum(regional_saving_uncertainty_PADD_4.transpose()),
                                       'PADD_5': get_copula_sum(regional_saving_uncertainty_PADD_5.transpose())})

regional_saving_result.to_excel(folder + f'results/regional/saving_{date.today()}.xlsx')

# decarbonization
# !!! update the file here if necessary
regional_decarbonization_uncertainty = pd.read_excel(folder + 'results/qualified_facility/integrated_decarbonization_kg_CO2_eq_per_day_2025-02-18.xlsx')

regional_decarbonization_uncertainty.drop('Unnamed: 0', axis=1, inplace=True)
regional_decarbonization_uncertainty = regional_decarbonization_uncertainty.transpose()
regional_decarbonization_uncertainty.reset_index(inplace=True, names='CWNS')
regional_decarbonization_uncertainty = regional_decarbonization_uncertainty.merge(WRRF_PADD, how='left', on='CWNS')

regional_decarbonization_uncertainty_PADD_1 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 1].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_2 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 2].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_3 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 3].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_4 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 4].loc[:, 0:999]
regional_decarbonization_uncertainty_PADD_5 = regional_decarbonization_uncertainty[regional_decarbonization_uncertainty['WRRF_PADD'] == 5].loc[:, 0:999]

regional_decarbonization_result = pd.DataFrame({'PADD_1': get_copula_sum(regional_decarbonization_uncertainty_PADD_1.transpose()),
                                                'PADD_2': get_copula_sum(regional_decarbonization_uncertainty_PADD_2.transpose()),
                                                'PADD_3': get_copula_sum(regional_decarbonization_uncertainty_PADD_3.transpose()),
                                                'PADD_4': get_copula_sum(regional_decarbonization_uncertainty_PADD_4.transpose()),
                                                'PADD_5': get_copula_sum(regional_decarbonization_uncertainty_PADD_5.transpose())})

regional_decarbonization_result.to_excel(folder + f'results/regional/decarbonization_{date.today()}.xlsx')

#%% regional uncertainty (visualization)

# saving
# !!! update these files if necessary
saving = pd.read_excel(folder + 'results/regional/saving_2025-02-20.xlsx')
# MM$/day
saving = saving/1000000

# decarbonization
# !!! update these files if necessary
decarbonization = pd.read_excel(folder + 'results/regional/decarbonization_2025-02-20.xlsx')
# tonne CO2 eq/day
decarbonization = decarbonization/1000

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize = (11, 10))

ax.set_xlim(0, 700)
ax.set_ylim(0, 1.8)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 800, 100))
plt.yticks(np.arange(0, 2, 0.2))

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 800, 100))
plt.yticks(np.arange(0, 2, 0.2))

ax_top.set_zorder(ax.get_zorder()+1)

ax.set_xlabel(r'$\mathbf{GHG\ reduction}$ [tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=35)
ax.set_ylabel(r'$\mathbf{Saving}$ [MM\$·day${^{-1}}$]', fontname='Arial', fontsize=35)

mathtext.FontConstantsBase.sup1 = 0.35

def add_rectangle(region, color, edgecolor):
    rectangle_fill = Rectangle((decarbonization[region].quantile(0.05), saving[region].quantile(0.05)),
                               decarbonization[region].quantile(0.95) - decarbonization[region].quantile(0.05),
                               saving[region].quantile(0.95) - saving[region].quantile(0.05),
                               fc=color, alpha=0.8, zorder=0)
    ax.add_patch(rectangle_fill)
    
    rectangle_edge = Rectangle((decarbonization[region].quantile(0.05), saving[region].quantile(0.05)),
                               decarbonization[region].quantile(0.95) - decarbonization[region].quantile(0.05),
                               saving[region].quantile(0.95) - saving[region].quantile(0.05),
                               color=edgecolor, lw=3, fc='none', alpha=1, zorder=0)
    ax.add_patch(rectangle_edge)
    
def add_line(region, color):
    plt.plot([decarbonization[region].quantile(0.25), decarbonization[region].quantile(0.75)],
             [saving[region].quantile(0.5), saving[region].quantile(0.5)],
             lw=3, color=color, solid_capstyle='round', zorder=1)
    plt.plot([decarbonization[region].quantile(0.5), decarbonization[region].quantile(0.5)],
             [saving[region].quantile(0.25), saving[region].quantile(0.75)],
             lw=3, color=color, solid_capstyle='round', zorder=1)

def add_point(region, edgecolor):
    ax_top.scatter(x=decarbonization[region].mean(),
               y=saving[region].mean(),
               marker='*',
               s=500,
               c='w',
               linewidths=3,
               alpha=1,
               edgecolor=edgecolor,
               zorder=2)

add_rectangle('PADD_2', g, dg)
add_line('PADD_2', dg)
add_point('PADD_2', dg)

add_rectangle('PADD_1', b, db)
add_line('PADD_1', db)
add_point('PADD_1', db)

add_rectangle('PADD_5', y, dy)
add_line('PADD_5', dy)
add_point('PADD_5', dy)

add_rectangle('PADD_3', r, dr)
add_line('PADD_3', dr)
add_point('PADD_3', dr)

add_rectangle('PADD_4', o, do)
add_line('PADD_4', do)
add_point('PADD_4', do)

#%% national uncertainty

# saving
# !!! update the file here if necessary
get_copula_sum('integrated_saving_dollar_per_day_2025-02-18')

# decarbonization
# !!! update the file here if necessary
get_copula_sum('integrated_decarbonization_kg_CO2_eq_per_day_2025-02-18')

#%% sludge transportation (heat map, HM) data preparation - part 1

filterwarnings('ignore')

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')

print(len(WRRF_input))

# $/tonne
WRRF_input['waste_cost'] = sum(WRRF_input[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000
# kg CO2 eq/tonne
WRRF_input['waste_GHG'] =  sum(WRRF_input[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

# the units here do not matter
total_cost = (WRRF_input['waste_cost']*WRRF_input['total_sludge_amount_kg_per_year']).sum(axis=0)
total_GHG = (WRRF_input['waste_GHG']*WRRF_input['total_sludge_amount_kg_per_year']).sum(axis=0)
total_sludge = WRRF_input['total_sludge_amount_kg_per_year'].sum(axis=0)

average_cost = total_cost/total_sludge
average_GHG = total_GHG/total_sludge

total_ash = 0
total_lipid = 0
total_protein = 0
total_sludge = WRRF_input['total_sludge_amount_kg_per_year'].sum()

for i in range(len(WRRF_input)):
    if WRRF_input.iloc[i]['sludge_aerobic_digestion'] == 1:
        sludge_ash_value=0.436
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
    elif WRRF_input.iloc[i]['sludge_anaerobic_digestion'] == 1:
        sludge_ash_value=0.414
        sludge_lipid_value=0.193
        sludge_protein_value=0.510
    else:
        sludge_ash_value=0.231
        sludge_lipid_value=0.206
        sludge_protein_value=0.456

    total_ash += sludge_ash_value*WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']
    total_lipid += sludge_lipid_value*WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)
    total_protein += sludge_protein_value*WRRF_input.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)

average_ash_dw = total_ash/total_sludge
average_lipid_afdw = total_lipid/(total_sludge-total_ash)
average_protein_afdw = total_protein/(total_sludge-total_ash)

# !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
# set HM_sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
HM_sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion']
HM_sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
HM_sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

HM_results_dict = {'sludge_amount':[], 'sludge_transportation_distance':[],
                   'saving_5th':[], 'saving_50th':[], 'saving_95th':[],
                   'decarbonization_5th':[], 'decarbonization_50th':[], 'decarbonization_95th':[]}

HM_results = pd.DataFrame(HM_results_dict)

get_quantiles = lambda data, quantiles=(0.05, 0.5, 0.95): [data.quantile(q) for q in quantiles]

for size in np.linspace(20, 200, 10):
    for sludge_distance in np.linspace(20, 200, 10):
        print('\n\n', f'sludge amount: {size} metric tonne/day\n', f'sludge travel distance: {sludge_distance} km\n')
        
        # !!! replace nitrogen_fertilizer the other 2 nitrogen fertilizers
        sys = create_geospatial_system(size=size,
                                       sludge_transportation=1,
                                       sludge_distance=sludge_distance,
                                       biocrude_distance=(WRRF_input['real_distance_km']*WRRF_input['total_sludge_amount_kg_per_year']).sum()/WRRF_input['total_sludge_amount_kg_per_year'].sum(),
                                       average_sludge_dw_ash=HM_sludge_ash_values[1],
                                       average_sludge_afdw_lipid=HM_sludge_lipid_values[1],
                                       average_sludge_afdw_protein=HM_sludge_protein_values[1],
                                       anaerobic_digestion=None,
                                       aerobic_digestion=None,
                                       ww_2_dry_sludge_ratio=1,
                                       state='US',
                                       nitrogen_fertilizer='NH3',
                                       elec_GHG=(WRRF_input['kg_CO2e_kWh']*WRRF_input['total_sludge_amount_kg_per_year']).sum()/WRRF_input['total_sludge_amount_kg_per_year'].sum(),
                                       wage_adjustment=(WRRF_input.iloc[i]['wage_quotient']/100*WRRF_input['total_sludge_amount_kg_per_year']).sum()/WRRF_input['total_sludge_amount_kg_per_year'].sum())
        
        crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']=='US']['2022_min'].iloc[0],
                                  crude_oil_price_data[crude_oil_price_data['state']=='US']['2022_average'].iloc[0],
                                  crude_oil_price_data[crude_oil_price_data['state']=='US']['2022_max'].iloc[0]]
        
        DAP_price_values = [588, 984.5, 1381]
        anhydrous_ammonia_price_values = [1065, 1407.5, 1750]
        urea_price_values = [447, 855, 1263]
        UAN_price_values = [350, 605, 860]
        
        model = create_geospatial_model(system=sys,
                                        sludge_ash=HM_sludge_ash_values,
                                        sludge_lipid=HM_sludge_lipid_values,
                                        sludge_protein=HM_sludge_protein_values,
                                        crude_oil_price=crude_oil_price_values,
                                        DAP_price=DAP_price_values,
                                        anhydrous_ammonia_price=anhydrous_ammonia_price_values,
                                        urea_price=urea_price_values,
                                        UAN_price=UAN_price_values)
        
        kwargs = {'N':1000,'rule':'L','seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        idx = len(model.parameters)
        parameters = model.table.iloc[:, :idx]
        results = model.table.iloc[:, idx:]
        
        # $/day
        HM_saving = (average_cost - results[('Geospatial','Sludge management cost [$/tonne dry sludge]')])*size
        # kg CO2 eq/day
        HM_decarbonization = (average_GHG - results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')])*size
        
        HM_results.loc[len(HM_results.index)] = ([size, sludge_distance,] +
                                                 get_quantiles(HM_saving) +
                                                 get_quantiles(HM_decarbonization))

HM_results.to_excel(folder + f'results/heat_map/anhydrous_ammonia/heat_map_{size}_tonne_per_day_{sludge_distance}_km_{date.today()}.xlsx')

#%% sludge transportation (heat map, HM) data preparation - part 2

def plot_heat_map(nitrogen_fertilizer, solids_quantity, item):
    # !!! update the input file if necessary
    HM = pd.read_excel(folder + f'results/heat_map/{nitrogen_fertilizer}/heat_map_{solids_quantity}_tonne_per_day_200.0_km_2025-02-19.xlsx')
    
    HM_saving = HM[['sludge_amount','sludge_transportation_distance',item]]
    
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.labelsize'] = 38
    plt.rcParams['ytick.labelsize'] = 38
    
    plt.rcParams.update({'mathtext.fontset': 'custom'})
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
    fig, ax = plt.subplots(figsize=(12.5, 10))
    
    ax = plt.gca()
    
    ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, pad=6)
    
    ax.set_xlabel(r'$\mathbf{Average\ distance}$ [km]', fontname='Arial', fontsize=45)
    ax.set_ylabel(r'$\mathbf{Total\ solids}$ [tonne·day${^{-1}}$]', fontname='Arial', fontsize=45)
    
    mathtext.FontConstantsBase.sup1 = 0.35
    
    plt.xticks(np.arange(20, 220, 30))
    plt.yticks(np.arange((solids_quantity)*0.1, (solids_quantity)*1.1, (solids_quantity)*0.9/6))
    
    ax_top = ax.twiny()
    ax_top.set_xlim(20, 200)
    ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')
    
    plt.xticks(np.arange(20, 220, 30))
    plt.yticks(np.arange((solids_quantity)*0.1, (solids_quantity)*1.1, (solids_quantity)*0.9/6))
    
    ax_right = ax.twinx()
    ax_right.set_ylim((solids_quantity)*0.1, (solids_quantity)*1.1)
    ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')
    
    plt.xticks(np.arange(20, 220, 30))
    plt.yticks(np.arange((solids_quantity)*0.1, (solids_quantity)*1.1, (solids_quantity)*0.9/6))
    
    try:
        color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [r, o, y, g, b][::-1])
    except ValueError:
        pass
    
    X = np.array(HM_saving['sludge_transportation_distance'])
    Y = np.array(HM_saving['sludge_amount'])
    Z = np.array(HM_saving[item]/1000)
    
    fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)
    
    fig.colorbar(fills, ax=ax)
    
    fig.delaxes(fig.axes[3])
    
    lines = ax.tricontour(X, Y, Z, levels=7, linewidths=3, linestyles='solid', colors='k')
    
    ax.clabel(lines, lines.levels, inline=True, fontsize=38)


#%% sludge transportation (heat map, HM) visualization - anhydrous_ammonia, 200.0 tonne per day, saving_50th

plot_heat_map('anhydrous_ammonia', 200.0, 'saving_50th')

#%% sludge transportation (heat map, HM) visualization - anhydrous_ammonia, 200.0 tonne per day, decarbonization_50th

plot_heat_map('anhydrous_ammonia', 200.0, 'decarbonization_50th')

#%% sludge transportation (heat map, HM) visualization - urea, 200.0 tonne per day, saving_50th

plot_heat_map('urea', 200.0, 'saving_50th')

#%% sludge transportation (heat map, HM) visualization - urea, 200.0 tonne per day, saving_50th

plot_heat_map('urea', 200.0, 'decarbonization_50th')

#%% sludge transportation (heat map, HM) visualization - UAN, 200.0 tonne per day, saving_50th

plot_heat_map('UAN', 200.0, 'saving_50th')

#%% sludge transportation (heat map, HM) visualization - UAN, 200.0 tonne per day, saving_50th

plot_heat_map('UAN', 200.0, 'decarbonization_50th')

#%% sludge transportation (heat map, HM) visualization - anhydrous_ammonia, 20.0 tonne per day, saving_50th

plot_heat_map('anhydrous_ammonia', 20.0, 'saving_50th')

#%% sludge transportation (heat map, HM) visualization - anhydrous_ammonia, 20.0 tonne per day, decarbonization_50th

plot_heat_map('anhydrous_ammonia', 20.0, 'decarbonization_50th')

#%% sludge transportation (heat map, HM) visualization - urea, 20.0 tonne per day, saving_50th

plot_heat_map('urea', 20.0, 'saving_50th')

#%% sludge transportation (heat map, HM) visualization - urea, 20.0 tonne per day, saving_50th

plot_heat_map('urea', 20.0, 'decarbonization_50th')

#%% sludge transportation (heat map, HM) visualization - UAN, 20.0 tonne per day, saving_50th

plot_heat_map('UAN', 20.0, 'saving_50th')

#%% sludge transportation (heat map, HM) visualization - UAN, 20.0 tonne per day, saving_50th

plot_heat_map('UAN', 20.0, 'decarbonization_50th')

#%% future coverage (data preparation for coverage visualization)

distance_x = []
coverage_y = []

for distance_threshold in np.linspace(0, 30, 151):
    print(distance_threshold)
    
    # max distance to consider WWTPs as neighbors
    distance_threshold_km = distance_threshold
    solids_threshold = 8
    # /1.3 to convert the real distance to linear distance, 1.3 is the average ratio from WRRF-oil refinery transportation
    distance_constant = 80/1.3
    
    WRRF_coverage = []
    for i in range(len(WRRF)):
        WRRF_x = WRRF.iloc[i].latitude
        WRRF_y = WRRF.iloc[i].longitude
        solids = WRRF.iloc[i].total_sludge_amount_kg_per_year/1000/365
        WRRF_coverage.append({'id': i, 'pos': (WRRF_x, WRRF_y), 'solids': solids})
    
    G = nx.Graph()
    for facility in WRRF_coverage:
        G.add_node(facility['id'], pos=facility['pos'], solids=facility['solids'])
    
    # convert the distance threshold from km to degrees.
    # roughly, 1 degree latitude ≈ 111 km.
    degree_threshold = distance_threshold_km / 111.0
    
    points = np.array([facility['pos'] for facility in WRRF_coverage])
    tree = KDTree(points)
    # query candidate pairs using the approximate threshold in degrees
    candidate_pairs = tree.query_pairs(r=degree_threshold)
    print(f"Candidate pairs from KDTree: {len(candidate_pairs)}")
    
    @njit
    def haversine(lat1, lon1, lat2, lon2):
        """
        calculate the great-circle distance between two points on Earth (in km)
        using the haversine formula
        """
        # Earth radius in km
        R = 6371.0
        # convert decimal degrees to radians
        lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat/2.0)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
        c = 2*np.arcsin(np.sqrt(a))
        return R*c
    
    times = 0   
    
    for i, j in candidate_pairs:
        if times%500000 == 0:
            print(times)
        times += 1
        
        pos_i = points[i]
        pos_j = points[j]
        geo_distance = haversine(pos_i[0], pos_i[1], pos_j[0], pos_j[1])
        if geo_distance < distance_threshold_km:
            G.add_edge(i, j, weight=geo_distance)
    
    # print(f"Total edges added: {G.number_of_edges()}")
    
    # cluster identification
    clusters = list(nx.connected_components(G))
    # print(f"Found {len(clusters)} clusters.")
    
    # hub identification
    hubs = []
    for cluster in clusters:
        cluster_nodes = list(cluster)
        
        total_solids = sum(G.nodes[node]['solids'] for node in cluster_nodes)
        
        # compute the centroid as the candidate hub position
        xs = [G.nodes[node]['pos'][0] for node in cluster_nodes]
        ys = [G.nodes[node]['pos'][1] for node in cluster_nodes]
        centroid = (np.mean(xs), np.mean(ys))
        
        total_distance = sum(haversine(G.nodes[node]['pos'][0],
                                       G.nodes[node]['pos'][1],
                                       centroid[0],
                                       centroid[1])
                             for node in cluster_nodes)
        
        if total_solids > solids_threshold and total_distance < total_solids * distance_constant:
            hubs.append({'hub_pos': centroid, 
                         'cluster': cluster_nodes, 
                         'total_solids': total_solids, 
                         'total_distance': total_distance})
        #     print(f"Cluster {cluster_nodes}: total_solids = {total_solids:.2f}, "
        #           f"total_distance = {total_distance:.2f} -> HUB accepted at {centroid}")
        # else:
        #     print(f"Cluster {cluster_nodes}: Does not meet hub criteria (solids = {total_solids:.2f}, "
        #           f"distance = {total_distance:.2f}).")
    
    # calculate coverage
    total_solids_all = sum(facility['solids'] for facility in WRRF_coverage)
    covered_solids = sum(hub['total_solids'] for hub in hubs)
    percentage_covered = covered_solids / total_solids_all * 100
    
    distance_x.append(distance_threshold)
    coverage_y.append(percentage_covered)
    
#%% future coverage (coverage visualization)

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38
plt.rcParams['font.sans-serif'] = 'Arial'

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

fig, ax = plt.subplots(figsize=(10, 10))

ax = plt.gca()

ax.set_xlim(0, 30)
ax.set_ylim(0, 100)

ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.set_xlabel(r'$\mathbf{Max.\ clustering\ distance}$ [km]', fontname='Arial', fontsize=45)
ax.set_ylabel(r'$\mathbf{Coverage}$ [%]', fontname='Arial', fontsize=45, linespacing=0.8)

mathtext.FontConstantsBase.sup1 = 0.35

plt.xticks(np.arange(0, 35, 5))
plt.yticks(np.arange(0, 120, 20))

ax_top = ax.twiny()
ax_top.set_xlim(0, 30)
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 35, 5))
plt.yticks(np.arange(0, 120, 20))

ax_right = ax.twinx()
ax_right.set_ylim(0, 100)
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

plt.xticks(np.arange(0, 35, 5))
plt.yticks(np.arange(0, 120, 20))

ax.plot(distance_x,
        coverage_y,
        c='k',
        linewidth=3)

ax.plot([0, distance_x[coverage_y.index(max(coverage_y))], distance_x[coverage_y.index(max(coverage_y))]],
        [max(coverage_y), max(coverage_y), 0],
        c='k',
        linestyle='--',
        linewidth=3)

ax.scatter(distance_x[coverage_y.index(max(coverage_y))],
           max(coverage_y),
           c=r,
           s=500,
           linewidths=3,
           edgecolors='k',
           zorder=2)

#%% future coverage (data preparation for map visualization)

# max distance to consider WWTPs as neighbors
distance_threshold_km = 16.2
solids_threshold = 8
# /1.3 to convert the real distance to linear distance, 1.3 is the average ratio from WRRF-oil refinery transportation
distance_constant = 80/1.3

WRRF_coverage = []
for i in range(len(WRRF)):
    WRRF_coverage.append({'id': i,
                          'pos': (WRRF.iloc[i].latitude,
                                           WRRF.iloc[i].longitude),
                          'solids': WRRF.iloc[i].total_sludge_amount_kg_per_year/1000/365})

G = nx.Graph()
for facility in WRRF_coverage:
    G.add_node(facility['id'], pos=facility['pos'], solids=facility['solids'])

# convert the distance threshold from km to degrees.
# roughly, 1 degree latitude ≈ 111 km.
degree_threshold = distance_threshold_km / 111.0

points = np.array([facility['pos'] for facility in WRRF_coverage])
tree = KDTree(points)
# query candidate pairs using the approximate threshold in degrees
candidate_pairs = tree.query_pairs(r=degree_threshold)
print(f"Candidate pairs from KDTree: {len(candidate_pairs)}")

@njit
def haversine(lat1, lon1, lat2, lon2):
    """
    calculate the great-circle distance between two points on Earth (in km)
    using the haversine formula
    """
    # Earth radius in km
    R = 6371.0
    # convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
    c = 2*np.arcsin(np.sqrt(a))
    return R*c

times = 0   

for i, j in candidate_pairs:
    if times%500000 == 0:
        print(times)
    times += 1
    
    pos_i = points[i]
    pos_j = points[j]
    geo_distance = haversine(pos_i[0], pos_i[1], pos_j[0], pos_j[1])
    if geo_distance < distance_threshold_km:
        G.add_edge(i, j, weight=geo_distance)

# print(f"Total edges added: {G.number_of_edges()}")

# cluster identification
clusters = list(nx.connected_components(G))
# print(f"Found {len(clusters)} clusters.")

# hub identification
hubs = []
for cluster in clusters:
    cluster_nodes = list(cluster)
    
    total_solids = sum(G.nodes[node]['solids'] for node in cluster_nodes)
    
    # compute the centroid as the candidate hub position
    xs = [G.nodes[node]['pos'][0] for node in cluster_nodes]
    ys = [G.nodes[node]['pos'][1] for node in cluster_nodes]
    centroid = (np.mean(xs), np.mean(ys))
    
    total_distance = sum(haversine(G.nodes[node]['pos'][0],
                                   G.nodes[node]['pos'][1],
                                   centroid[0],
                                   centroid[1])
                         for node in cluster_nodes)
    
    if total_solids > solids_threshold and total_distance < total_solids * distance_constant:
        hubs.append({'hub_pos': centroid, 
                     'cluster': cluster_nodes, 
                     'total_solids': total_solids, 
                     'total_distance': total_distance})
    #     print(f"Cluster {cluster_nodes}: total_solids = {total_solids:.2f}, "
    #           f"total_distance = {total_distance:.2f} -> HUB accepted at {centroid}")
    # else:
    #     print(f"Cluster {cluster_nodes}: Does not meet hub criteria (solids = {total_solids:.2f}, "
    #           f"distance = {total_distance:.2f}).")

#%% future coverage (map visualization)

pos_dict = nx.get_node_attributes(G, 'pos')

WRRF_cluster = pd.DataFrame()

WRRF_cluster['index'] = pos_dict.keys()
WRRF_cluster['latitide'] = [i[0] for i in pos_dict.values()]
WRRF_cluster['longitude'] = [i[1] for i in pos_dict.values()]

covered_WRRF_index = [k for j in [i['cluster'] for i in hubs] for k in j]

WRRF_cluster.loc[WRRF_cluster['index'].isin(covered_WRRF_index), 'color'] = b
WRRF_cluster.loc[~WRRF_cluster['index'].isin(covered_WRRF_index), 'color'] = a

WRRF_cluster = gpd.GeoDataFrame(WRRF_cluster, crs='EPSG:4269',
                                geometry=gpd.points_from_xy(x=WRRF_cluster.longitude,
                                                            y=WRRF_cluster.latitide))
WRRF_cluster = WRRF_cluster.to_crs(crs='EPSG:3857')

new_hub_position = pd.DataFrame()

new_hub_position['latitide'] = [j[0] for j in [i['hub_pos'] for i in hubs if len(i['cluster']) != 1]]
new_hub_position['longitude'] = [j[1] for j in [i['hub_pos'] for i in hubs if len(i['cluster']) != 1]]

new_hub_position = gpd.GeoDataFrame(new_hub_position, crs='EPSG:4269',
                                    geometry=gpd.points_from_xy(x=new_hub_position.longitude,
                                                                y=new_hub_position.latitide))
new_hub_position = new_hub_position.to_crs(crs='EPSG:3857')


old_hub_position = pd.DataFrame()

old_hub_position['latitide'] = [j[0] for j in [i['hub_pos'] for i in hubs if len(i['cluster']) == 1]]
old_hub_position['longitude'] = [j[1] for j in [i['hub_pos'] for i in hubs if len(i['cluster']) == 1]]

old_hub_position = gpd.GeoDataFrame(old_hub_position, crs='EPSG:4269',
                                    geometry=gpd.points_from_xy(x=old_hub_position.longitude,
                                                                y=old_hub_position.latitide))
old_hub_position = old_hub_position.to_crs(crs='EPSG:3857')

fig, ax = plt.subplots(figsize=(30, 30))

US.plot(ax=ax, color='w', edgecolor='k', linewidth=3)

WRRF_cluster_N = WRRF_cluster[WRRF_cluster['color'] == a]

WRRF_cluster_Y = WRRF_cluster[WRRF_cluster['color'] == b]

WRRF_cluster_N.plot(ax=ax, color=WRRF_cluster_N['color'], markersize=30, edgecolor='k', linewidth=0.5, alpha=0.3)

WRRF_cluster_Y.plot(ax=ax, color=WRRF_cluster_Y['color'], markersize=30, edgecolor='k', linewidth=0.5, alpha=1)

new_hub_position.plot(ax=ax, color=y, marker='*', markersize=600, edgecolor='k', linewidth=1.5, alpha=1)

old_hub_position.plot(ax=ax, color=r, marker='*', markersize=600, edgecolor='k', linewidth=1.5, alpha=1)

rectangle_edge = Rectangle((-13730000, 2930000), 1980000, 630000,
                           color='k', lw=3, fc='none', alpha=1)
ax.add_patch(rectangle_edge)

# note the size of legends does not match exactly with the size of WRRFs
ax.scatter(x=-13640000, y=3450000, marker='o', s=50, c=b, linewidths=0.5,
           alpha=1, edgecolor='k')
ax.scatter(x=-13640000, y=3310000, marker='o', s=50, c=a, linewidths=0.5,
           alpha=0.5, edgecolor='k')
ax.scatter(x=-13640000, y=3170000, marker='*', s=600, c=y, linewidths=1.5,
           alpha=1, edgecolor='k')
ax.scatter(x=-13640000, y=3030000, marker='*', s=600, c=r, linewidths=1.5,
           alpha=1, edgecolor='k')

plt.figtext(0.2, 0.3751, 'covered WRRFs', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})
plt.figtext(0.2, 0.3596, 'uncovered WRRFs', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})
plt.figtext(0.2, 0.3441, 'hubs covering multiple WRRFs', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})
plt.figtext(0.2, 0.3286, 'hubs covering one WRRF', fontname='Arial', fontdict={'fontsize': 30,'color':'k','style':'italic'})

ax.set_aspect(1)

ax.set_axis_off()

#%% writing results

qualified_facility_biocrude = pd.read_excel(folder + 'results/qualified_facility/integrated_biocrude_BPD_2025-02-18.xlsx')
qualified_facility_biocrude.drop('Unnamed: 0', axis=1, inplace=True)
qualified_facility_biocrude_max = pd.DataFrame(qualified_facility_biocrude.max())
qualified_facility_biocrude_max.reset_index(names='CWNS', inplace=True)
qualified_facility_biocrude_max.rename(columns={0:'BPD'}, inplace=True)

# !!! update the input file if necessary
WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# removal WRRFs with no real_distance_km
WRRF_input = WRRF_input.dropna(subset='real_distance_km')
WRRF_input = WRRF_input[['CWNS','Site ID','capacity']]

BPD_capacity = qualified_facility_biocrude_max.merge(WRRF_input, how='left', on='CWNS')
BPD_capacity = BPD_capacity.groupby(by=['Site ID','capacity']).sum()
BPD_capacity.reset_index(inplace=True)
BPD_capacity['ratio'] = BPD_capacity['BPD']/BPD_capacity['capacity']/1000000
print(f'The highest blending ratio in the near-term scenario is {BPD_capacity["ratio"].max()*100: .3f}%.')

#%% line plot of all WRRFs (superseded)

# # !!! update the file here if necessary
# line_plot = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
# # line_plot = line_plot[line_plot['USD_decarbonization'].notna()]
# # line_plot = line_plot[line_plot['USD_decarbonization'] <= 0]

# line_plot = line_plot.copy()

# line_plot['dollar_per_kWh'] = line_plot['state'].apply(lambda x: elec_price[elec_price['state'] == x]['price'].iloc[0]/100)
# line_plot['crude_oil_dollar_per_barrel'] = line_plot['state'].apply(lambda x: crude_oil_price_data[crude_oil_price_data['state'] == x]['2022_average'].iloc[0])

# line_plot['DAP_price'] = line_plot['state'].apply(lambda x: fertilizer_price_uncertainty[x]['DAP'][1] if x in ['AL','IA','IL','NC','OK','SC'] else 984.5)
# line_plot['anhydrous_ammonia_price'] = line_plot['state'].apply(lambda x: fertilizer_price_uncertainty[x]['anhydrous_ammonia'][1] if x in ['IA','IL','OK'] else 1407.5)
# line_plot['urea_price'] = line_plot['state'].apply(lambda x: fertilizer_price_uncertainty[x]['urea'][1] if x in ['AL','IA','IL','NC','OK','SC'] else 855)
# line_plot['UAN_price'] = line_plot['state'].apply(lambda x: fertilizer_price_uncertainty[x]['UAN'][1] if x in ['AL','IA','IL','NC','OK','SC'] else 605)

# line_plot = line_plot[['state','flow_2022_MGD_final','total_sludge_amount_kg_per_year',
#                        'sludge_aerobic_digestion','sludge_anaerobic_digestion','landfill',
#                        'land_application','incineration','DAP_price',
#                        'anhydrous_ammonia_price','urea_price','UAN_price',
#                        'dollar_per_kWh','crude_oil_dollar_per_barrel','kg_CO2e_kWh',
#                        'wage_quotient','nitrogen_fertilizer','real_distance_km',
#                        'income_tax','saving','CO2_reduction',
#                        'USD_decarbonization']]

# line_plot.loc[(line_plot['sludge_anaerobic_digestion'] == 1) | (line_plot['sludge_aerobic_digestion'] == 1), 'digestion'] = 1
# line_plot['digestion'] = line_plot['digestion'].fillna(0)

# line_plot.loc[(line_plot['incineration'] != 0) & (line_plot['landfill'] == 0) & (line_plot['land_application'] == 0), 'incineration_only'] = 1
# line_plot['incineration_only'] = line_plot['incineration_only'].fillna(0)

# line_plot['unit_saving'] = line_plot['saving']/line_plot['total_sludge_amount_kg_per_year']

# plt.rcParams['axes.linewidth'] = 3
# plt.rcParams['xtick.labelsize'] = 38
# plt.rcParams['ytick.labelsize'] = 38

# plt.rcParams.update({'mathtext.fontset': 'custom'})
# plt.rcParams.update({'mathtext.default': 'regular'})
# plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

# fig, ax = plt.subplots(figsize=(20, 10))

# ax = plt.gca()

# for i in range(len(line_plot)):
#     if line_plot['USD_decarbonization'].iloc[i] > 0:
#         plt.plot([0, 1, 2, 3],
#                  [line_plot['digestion'].iloc[i],
#                   line_plot['nitrogen_fertilizer'].iloc[i]],
#                   (line_plot['total_sludge_amount_kg_per_year'].iloc[i] - line_plot['total_sludge_amount_kg_per_year'].min())/\
#                       (line_plot['total_sludge_amount_kg_per_year'].max()-line_plot['total_sludge_amount_kg_per_year'].min()),
#                   (line_plot['real_distance_km'].iloc[i] - line_plot['real_distance_km'].min())/\
#                       (line_plot['real_distance_km'].max()-line_plot['real_distance_km'].min()),
#                   (line_plot['dollar_per_kWh'].iloc[i] - line_plot['dollar_per_kWh'].min())/\
#                       (line_plot['dollar_per_kWh'].max()-line_plot['dollar_per_kWh'].min()),
#                   (line_plot['kg_CO2e_kWh'].iloc[i] - line_plot['kg_CO2e_kWh'].min())/\
#                       (line_plot['kg_CO2e_kWh'].max()-line_plot['kg_CO2e_kWh'].min()),
#                   (line_plot['wage_quotient'].iloc[i] - line_plot['wage_quotient'].min())/\
#                       (line_plot['wage_quotient'].max()-line_plot['wage_quotient'].min()),
#                  lw=5*((line_plot['saving'].iloc[i] - line_plot['saving'].min())/\
#                        (line_plot['saving'].max()-line_plot['saving'].min())),
#                  color='k', solid_capstyle='round', zorder=0,
#                  alpha=((line_plot['saving'].iloc[i] - line_plot['saving'].min())/\
#                         (line_plot['saving'].max()-line_plot['saving'].min()))**0.5)

#%% statistic tests (superseded)

# import seaborn as sns

# test = line_plot[line_plot['USD_decarbonization'].notna()]

# test_1 = test[test['USD_decarbonization'] <= 0]
# test_2 = test[test['USD_decarbonization'] > 0]

# plt.rcParams['axes.linewidth'] = 3
# plt.rcParams['xtick.labelsize'] = 38
# plt.rcParams['ytick.labelsize'] = 38

# plt.rcParams.update({'mathtext.fontset': 'custom'})
# plt.rcParams.update({'mathtext.default': 'regular'})
# plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

# fig, ax = plt.subplots(figsize=(20, 10))

# sns.distplot(np.log(test['total_sludge_amount_kg_per_year']), hist = False, kde = True,
#              kde_kws = {'shade': True, 'linewidth': 3})

# sns.distplot(test_1['wage_quotient'], hist = False, kde = True,
#                   kde_kws = {'shade': True, 'linewidth': 3})

# sns.distplot(test_2['wage_quotient'], hist = False, kde = True,
#                   kde_kws = {'shade': True, 'linewidth': 3})

# from scipy import stats

# stats.kstest(test_1['wage_quotient'].dropna(), test_2['wage_quotient'].dropna())

# stats.ttest_ind(test_1['dollar_per_kWh'].dropna(), test_2['dollar_per_kWh'].dropna())

# from scipy.stats import ranksums

# ranksums(test_1['real_distance_km'].dropna(), test_2['real_distance_km'].dropna())

#%% PCA-part 1 (superseded)

# from sklearn.decomposition import PCA
# from sklearn.preprocessing import PowerTransformer, StandardScaler, OneHotEncoder
# from sklearn.compose import ColumnTransformer

# line_plot.loc[(line_plot['USD_decarbonization'].notna()) & (line_plot['USD_decarbonization'] <=0), 'qualified'] = 1

# line_plot.loc[line_plot['incineration'] > 0, 'incineration_yn'] = 1
# line_plot.loc[line_plot['landfill'] > 0, 'landfill_yn'] = 1
# line_plot.loc[line_plot['land_application'] > 0, 'land_application_yn'] = 1

# line_plot.fillna({'qualified': 0,
#                   'incineration_yn': 0,
#                   'landfill_yn': 0,
#                   'land_application_yn': 0},
#                  inplace=True)

# transformer = PowerTransformer(method='yeo-johnson')

# transformed_feature = transformer.fit_transform(line_plot[['flow_2022_MGD_final','total_sludge_amount_kg_per_year',
#                                                            'incineration','landfill','land_application']])

# line_plot[['flow_2022_MGD_final','total_sludge_amount_kg_per_year',
#            'incineration','landfill','land_application']] = transformed_feature

# import seaborn as sns

# corr_test = line_plot[['total_sludge_amount_kg_per_year','digestion','landfill',
#                        'land_application','incineration','dollar_per_kWh',
#                        'crude_oil_dollar_per_barrel','kg_CO2e_kWh','wage_quotient',
#                        'nitrogen_fertilizer','real_distance_km','income_tax']]

# plt.figure(figsize=(10, 5))
# # example: correlation heatmap
# correlation_matrix = corr_test.corr()
# sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm')
# plt.show()

# # identify numerical and categorical columns
# num_features = ['flow_2022_MGD_final',
#                 'total_sludge_amount_kg_per_year',
#                 'landfill',
#                 'land_application',
#                 'incineration',
#                 'dollar_per_kWh',
#                 'kg_CO2e_kWh',
#                 'wage_quotient',
#                 'crude_oil_dollar_per_barrel',
#                 'real_distance_km',
#                 'income_tax']

# cat_features = ['digestion',
#                 'nitrogen_fertilizer',
#                 'incineration_yn',
#                 'landfill_yn',
#                 'land_application_yn']

# # preprocessing: standardizing numerical data and one-hot encoding categorical data
# preprocessor = ColumnTransformer([('num', StandardScaler(), num_features),
#                                   ('cat', OneHotEncoder(), cat_features)])

# # apply transformation
# processed_data = preprocessor.fit_transform(line_plot)

# # apply PCA
# pca = PCA(n_components=5)
# pca_result = pca.fit_transform(processed_data)

# # convert PCA results into a DataFrame
# pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])

# line_plot.loc[line_plot['unit_saving'] <= line_plot['unit_saving'].median(), 'qualified_qualified'] = 1
# line_plot.loc[line_plot['unit_saving'] > line_plot['unit_saving'].median(), 'qualified_qualified'] = 0

# # add continuous target for color mapping
# pca_df['Target'] = line_plot['qualified']

# # explained variance ratio
# explained_variance = pca.explained_variance_ratio_

# pca_df.sort_values('Target', ascending=True, inplace=True)

# # plot PCA result with a color gradient for the continuous target
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
# sc = ax.scatter(pca_df['PC1'], pca_df['PC2'], pca_df['PC3'], c=pca_df['Target'], edgecolors='k', alpha=0.5)

# # customize plot
# plt.colorbar(sc, label='Target Value')  # Add color scale legend
# plt.xlabel(f'PC1 ({explained_variance[0]:.2%} Variance)')
# plt.ylabel(f'PC2 ({explained_variance[1]:.2%} Variance)')
# ax.set_zlabel(f'PC2 ({explained_variance[1]:.2%} Variance)')
# plt.title('PCA with Continuous Target Coloring')
# plt.grid(True)
# plt.show()

#%% PCA-part 2 (superseded)

# # get PCA loadings (i.e., contributions of each original feature to each principal component)
# loadings = pd.DataFrame(pca.components_, columns=(num_features+\
#                                                   list(preprocessor.named_transformers_['cat'].get_feature_names_out(cat_features))),
#                         index=[f'PC{i+1}' for i in range(pca.n_components_)])

# print("PCA Loadings (Feature Contributions to Each PC):")
# print(loadings)

#%% PCA-part 3 (superseded)

# import seaborn as sns

# # plot feature contributions to PC1
# plt.figure(figsize=(10, 5))
# sns.barplot(x=loadings.columns, y=loadings.loc['PC1'], palette='viridis')
# plt.xticks(rotation=45, ha='right')
# plt.ylabel('Feature Contribution')
# plt.title('Feature Contributions to PC1')
# plt.grid(True)
# plt.show()

#%% PCA-part 4 (superseded)

# # correlation of PCs with original features
# pca_correlation = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(pca.n_components_)])
# pca_correlation[num_features + cat_features] = line_plot[num_features + cat_features]

# # display correlation matrix
# print(pca_correlation.corr().round(2))

#%% PCA-part 5 (superseded)

# from sklearn.preprocessing import MinMaxScaler

# # scale loadings for better visualization
# scaler = MinMaxScaler(feature_range=(-1, 1))
# # transpose to align features
# scaled_loadings = scaler.fit_transform(loadings.T)

# # create biplot
# plt.figure(figsize=(20,20))
# sc = plt.scatter(pca_df['PC1'], pca_df['PC2'], c=pca_df['Target'], cmap='viridis', edgecolors='k', alpha=0.8)
# plt.colorbar(sc, label='Target Value')

# # plot arrows for feature contributions
# for i, feature in enumerate(loadings.columns):
#     plt.arrow(0, 0, scaled_loadings[i, 0], scaled_loadings[i, 1], color='red', alpha=0.7, head_width=0.05)
#     plt.text(s=feature, x=scaled_loadings[i, 0], y=scaled_loadings[i, 1], color='w', fontsize=10)

# plt.xlabel(f'PC1 ({explained_variance[0]:.2%} Variance)')
# plt.ylabel(f'PC2 ({explained_variance[1]:.2%} Variance)')
# plt.title('PCA Biplot (Projection + Feature Contributions)')
# plt.grid(True)
# plt.show()

#%% PCA-part 6 (superseded)

# import seaborn as sns
# import matplotlib.pyplot as plt

# plt.figure(figsize=(15,8))
# sns.heatmap(loadings, annot=True, cmap='coolwarm', center=0)
# plt.title("PCA Feature Loadings")
# plt.show()

#%% sludge transportation (Urbana-Champaign, two separated WRRFs) (superseded)

# filterwarnings('ignore')

# # !!! update the input file if necessary
# WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# # removal WRRFs with no real_distance_km
# WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# # $/tonne
# WRRF_input['waste_cost'] = sum(WRRF_input[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000
# # kg CO2 eq/tonne
# WRRF_input['waste_GHG'] =  sum(WRRF_input[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/WRRF_input['total_sludge_amount_kg_per_year']*1000

# Urbana_Champaign = WRRF_input[WRRF_input['CWNS'].isin([17000112001, 17000112002])]

# CU_uncertainty_saving = pd.DataFrame()
# CU_uncertainty_decarbonization = pd.DataFrame()

# for i in range(0, 2):
#     sys = create_geospatial_system(size=Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365,
#                                    sludge_transportation=0,
#                                    sludge_distance=100,
#                                    biocrude_distance=Urbana_Champaign.iloc[i]['real_distance_km'],
#                                    anaerobic_digestion=Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'],
#                                    aerobic_digestion=Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'],
#                                    ww_2_dry_sludge_ratio=1,
#                                    state=Urbana_Champaign.iloc[i]['state'],
#                                    nitrogen_fertilizer=Urbana_Champaign.iloc[i]['nitrogen_fertilizer'],
#                                    elec_GHG=Urbana_Champaign.iloc[i]['kg_CO2e_kWh'],
#                                    wage_adjustment=Urbana_Champaign.iloc[i]['wage_quotient']/100)
    
#     # if AD and AeD both exist, assume AeD (due to higher ash content) to be conservative
#     if Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'] == 1:
#         sludge_ash_values=[0.349, 0.436, 0.523, 'aerobic_digestion']
#         sludge_lipid_values=[0.154, 0.193, 0.232]
#         sludge_protein_values=[0.408, 0.510, 0.612]
#     elif Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'] == 1:
#         sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
#         sludge_lipid_values=[0.154, 0.193, 0.232]
#         sludge_protein_values=[0.408, 0.510, 0.612]
#     else:
#         sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
#         sludge_lipid_values=[0.080, 0.206, 0.308]
#         sludge_protein_values=[0.380, 0.456, 0.485]

#     crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']==Urbana_Champaign.iloc[i]['state']]['2022_min'].iloc[0],
#                               crude_oil_price_data[crude_oil_price_data['state']==Urbana_Champaign.iloc[i]['state']]['2022_average'].iloc[0],
#                               crude_oil_price_data[crude_oil_price_data['state']==Urbana_Champaign.iloc[i]['state']]['2022_max'].iloc[0]]

#     try:
#         DAP_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[i]['state']]['DAP']
#     except KeyError:
#         DAP_price_values = [588, 984.5, 1381]
    
#     try:
#         anhydrous_ammonia_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[i]['state']]['anhydrous_ammonia']
#     except KeyError:
#         anhydrous_ammonia_price_values = [1065, 1407.5, 1750]
    
#     try:
#         urea_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[i]['state']]['urea']
#     except KeyError:
#         urea_price_values = [447, 855, 1263]
    
#     try:
#         UAN_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[i]['state']]['UAN']
#     except KeyError:
#         UAN_price_values = [350, 605, 860]

#     model = create_geospatial_model(system=sys,
#                                     sludge_ash=sludge_ash_values,
#                                     sludge_lipid=sludge_lipid_values,
#                                     sludge_protein=sludge_protein_values,
#                                     crude_oil_price=crude_oil_price_values,
#                                     DAP_price=DAP_price_values,
#                                     anhydrous_ammonia_price=anhydrous_ammonia_price_values,
#                                     urea_price=urea_price_values,
#                                     UAN_price=UAN_price_values)
    
#     kwargs = {'N':1000,'rule':'L','seed':3221}
#     samples = model.sample(**kwargs)
#     model.load_samples(samples)
#     model.evaluate()
    
#     idx = len(model.parameters)
#     parameters = model.table.iloc[:, :idx]
#     results = model.table.iloc[:, idx:]
    
#     # $/day
#     CU_uncertainty_saving[Urbana_Champaign.iloc[i]['CWNS']] = (Urbana_Champaign.iloc[i]['waste_cost'] - results[('Geospatial','Sludge management cost [$/tonne dry sludge]')])*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365
#     # kg CO2 eq/day
#     CU_uncertainty_decarbonization[Urbana_Champaign.iloc[i]['CWNS']] = (Urbana_Champaign.iloc[i]['waste_GHG'] - results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')])*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']/1000/365

#%% sludge transportation (Urbana-Champaign, after sludge transportation) (superseded)

# # note the sludge is transported from 17000112002 to 17000112001

# CU_combined_uncertainty_saving = pd.DataFrame()
# CU_combined_uncertainty_decarbonization = pd.DataFrame()

# # $/tonne
# average_cost = (Urbana_Champaign.iloc[0]['waste_cost']*Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] +\
#                 Urbana_Champaign.iloc[1]['waste_cost']*Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']) /\
#                (Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] + Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year'])
# # kg CO2 eq/tonne
# average_GHG = (Urbana_Champaign.iloc[0]['waste_GHG']*Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] +\
#                Urbana_Champaign.iloc[1]['waste_GHG']*Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']) /\
#               (Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year'] + Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year'])

# distance_to_refinery = Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112001]['real_distance_km'].iloc[0]

# total_ash = 0
# total_lipid = 0
# total_protein = 0
# total_sludge = Urbana_Champaign.iloc[0]['total_sludge_amount_kg_per_year']+Urbana_Champaign.iloc[1]['total_sludge_amount_kg_per_year']
# total_size = total_sludge/1000/365

# for i in range(len(Urbana_Champaign)):
#     if Urbana_Champaign.iloc[i]['sludge_aerobic_digestion'] == 1:
#         sludge_ash_value=0.436
#         sludge_lipid_value=0.193
#         sludge_protein_value=0.510
#     elif Urbana_Champaign.iloc[i]['sludge_anaerobic_digestion'] == 1:
#         sludge_ash_value=0.414
#         sludge_lipid_value=0.193
#         sludge_protein_value=0.510
#     else:
#         sludge_ash_value=0.231
#         sludge_lipid_value=0.206
#         sludge_protein_value=0.456

#     total_ash += sludge_ash_value*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']
#     total_lipid += sludge_lipid_value*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)
#     total_protein += sludge_protein_value*Urbana_Champaign.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)

# average_ash_dw = total_ash/total_sludge
# average_lipid_afdw = total_lipid/(total_sludge-total_ash)
# average_protein_afdw = total_protein/(total_sludge-total_ash)

# # !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
# # set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
# sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion']
# sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
# sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

# assert Urbana_Champaign.iloc[0]['state'] == Urbana_Champaign.iloc[1]['state']
# assert Urbana_Champaign.iloc[0]['kg_CO2e_kWh'] == Urbana_Champaign.iloc[1]['kg_CO2e_kWh']
# assert Urbana_Champaign.iloc[0]['nitrogen_fertilizer'] == Urbana_Champaign.iloc[1]['nitrogen_fertilizer']

# # 16.2 km is the distance between these two WRRFs; the distance was manually obtained from Google Maps
# sys = create_geospatial_system(size=total_size,
#                                sludge_transportation=1,
#                                sludge_distance=16.2*Urbana_Champaign[Urbana_Champaign['CWNS'] == 17000112002]['total_sludge_amount_kg_per_year'].iloc[0]/total_sludge,
#                                biocrude_distance=distance_to_refinery,
#                                average_sludge_dw_ash=sludge_ash_values[1],
#                                average_sludge_afdw_lipid=sludge_lipid_values[1],
#                                average_sludge_afdw_protein=sludge_protein_values[1],
#                                anaerobic_digestion=None,
#                                aerobic_digestion=None,
#                                ww_2_dry_sludge_ratio=1,
#                                state=Urbana_Champaign.iloc[0]['state'],
#                                nitrogen_fertilizer=Urbana_Champaign.iloc[0]['nitrogen_fertilizer'],
#                                elec_GHG=Urbana_Champaign.iloc[0]['kg_CO2e_kWh'],
#                                wage_adjustment=Urbana_Champaign.iloc[i]['wage_quotient']/100)

# crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']==Urbana_Champaign.iloc[0]['state']]['2022_min'].iloc[0],
#                           crude_oil_price_data[crude_oil_price_data['state']==Urbana_Champaign.iloc[0]['state']]['2022_average'].iloc[0],
#                           crude_oil_price_data[crude_oil_price_data['state']==Urbana_Champaign.iloc[0]['state']]['2022_max'].iloc[0]]

# try:
#     DAP_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[0]['state']]['DAP']
# except KeyError:
#     DAP_price_values = [588, 984.5, 1381]

# try:
#     anhydrous_ammonia_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[0]['state']]['anhydrous_ammonia']
# except KeyError:
#     anhydrous_ammonia_price_values = [1065, 1407.5, 1750]

# try:
#     urea_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[0]['state']]['urea']
# except KeyError:
#     urea_price_values = [447, 855, 1263]

# try:
#     UAN_price_values = fertilizer_price_uncertainty[Urbana_Champaign.iloc[0]['state']]['UAN']
# except KeyError:
#     UAN_price_values = [350, 605, 860]

# model = create_geospatial_model(system=sys,
#                                 sludge_ash=sludge_ash_values,
#                                 sludge_lipid=sludge_lipid_values,
#                                 sludge_protein=sludge_protein_values,
#                                 crude_oil_price=crude_oil_price_values,
#                                 DAP_price=DAP_price_values,
#                                 anhydrous_ammonia_price=anhydrous_ammonia_price_values,
#                                 urea_price=urea_price_values,
#                                 UAN_price=UAN_price_values)

# kwargs = {'N':1000,'rule':'L','seed':3221}
# samples = model.sample(**kwargs)
# model.load_samples(samples)
# model.evaluate()

# idx = len(model.parameters)
# parameters = model.table.iloc[:, :idx]
# results = model.table.iloc[:, idx:]

# # $/day
# CU_combined_uncertainty_saving['combined'] = (average_cost - results[('Geospatial','Sludge management cost [$/tonne dry sludge]')])*total_sludge/1000/365
# # kg CO2 eq/day
# CU_combined_uncertainty_decarbonization['combined'] = (average_GHG - results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')])*total_sludge/1000/365

# CU_saving_results = pd.concat([CU_uncertainty_saving, CU_combined_uncertainty_saving], axis=1)
# CU_decarbonization_result = pd.concat([CU_uncertainty_decarbonization, CU_combined_uncertainty_decarbonization], axis=1)

# CU_saving_results.to_excel(folder + f'results/Urbana-Champaign/saving_{date.today()}.xlsx')
# CU_decarbonization_result.to_excel(folder + f'results/Urbana-Champaign/decarbonization_{date.today()}.xlsx')

#%% Urbana-Champaign visualization (superseded)

# CU_saving_results = pd.read_excel(folder + 'results/Urbana-Champaign/saving_2025-02-17.xlsx')
# CU_saving_results.drop('Unnamed: 0', axis=1, inplace=True)

# fig = plt.figure(figsize=(10, 10))

# gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)

# def add_region(position, xlabel, color):
#     ax = fig.add_subplot(gs[0, position])
    
#     plt.rcParams['axes.linewidth'] = 3
#     plt.rcParams['xtick.labelsize'] = 38
#     plt.rcParams['ytick.labelsize'] = 38
    
#     plt.xticks(fontname='Arial')
#     plt.yticks(fontname='Arial')
    
#     plt.rcParams.update({'mathtext.fontset': 'custom'})
#     plt.rcParams.update({'mathtext.default': 'regular'})
#     plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
    
#     ax = plt.gca()
#     ax.set_ylim(-3000, 1000)
#     ax.set_xlabel(xlabel, fontname='Arial', fontsize=38, labelpad=15)
    
#     if position == 0:
#         ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)
#         ax.set_ylabel(r'$\mathbf{Saving}$ [\$·day${^{-1}}$]', fontname='Arial', fontsize=45)
        
#         mathtext.FontConstantsBase.sup1 = 0.35
    
#     elif position == 2:
#         ax.tick_params(direction='inout', length=0, labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
        
#         ax_right = ax.twinx()
#         ax_right.set_ylim(ax.get_ylim())
#         ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')
    
#     else:
#         ax.tick_params(direction='inout', labelbottom=False, bottom=False, top=False, left=False, right=False, labelcolor='none')
    
#     bp = ax.boxplot(CU_saving_results.iloc[:,position], showfliers=False, widths=0.7, patch_artist=True)
    
#     for box in bp['boxes']:
#         box.set(color='k', facecolor=color, linewidth=3)
    
#     for whisker in bp['whiskers']:
#         whisker.set(color='k', linewidth=3)
    
#     for median in bp['medians']:
#         median.set(color='k', linewidth=3)
    
#     for cap in bp['caps']:
#         cap.set(color='k', linewidth=3)
    
#     for flier in bp['fliers']:
#         flier.set(marker='o', markersize=7, markerfacecolor=color, markeredgewidth=1.5)

# add_region(0, 'Northeast', o)
# add_region(1, 'Southwest', o)
# add_region(2, 'Combined', b)

#%% sludge transportation (satellite, SA) (data preparation) (superseded)

# # !!! update the file here if necessary
# satellite_data_preparation = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
# satellite_data_preparation = satellite_data_preparation[satellite_data_preparation['USD_decarbonization'].notna()]
# satellite_data_preparation = satellite_data_preparation[satellite_data_preparation['USD_decarbonization'] <= 0]

# center_WRRF_data_preparation = satellite_data_preparation[satellite_data_preparation['facility']=='MWRDGC, Stickney Treatment Plant']

# # !!! update the input file if necessary
# WRRF_input = pd.read_excel(folder + 'HTL_geospatial_model_input_2025-02-10.xlsx')
# # removal WRRFs with no real_distance_km
# WRRF_input = WRRF_input.dropna(subset='real_distance_km')

# # select WRRFs that are within 200 km (linear distance) from the center WRRF
# linear_distance = []
# for i in range(len(WRRF_input)):
#     linear_distance.append(geopy.distance.geodesic((WRRF_input['latitude'].iloc[i], WRRF_input['longitude'].iloc[i]), (center_WRRF_data_preparation.latitude[0], center_WRRF_data_preparation.longitude[0])).km)

# WRRF_input['linear_distance_to_center_WRRF_km'] = linear_distance

# WRRF_input = WRRF_input[WRRF_input['linear_distance_to_center_WRRF_km'] <= 200]

# # remove WRRFs that are already decarbonized without sludge transportation
# need_transportation = WRRF_input[~WRRF_input['CWNS'].isin(list(satellite_data_preparation['CWNS']))].copy()

# need_transportation.sort_values(by='total_sludge_amount_kg_per_year', ascending=False, inplace=True)

# need_transportation['WRRF_location'] = list(zip(need_transportation.latitude, need_transportation.longitude))
# center_WRRF_location = (center_WRRF_data_preparation.latitude[0], center_WRRF_data_preparation.longitude[0])

# # =============================================================================
# # # code to generate the inventory
# # 
# # # !!! get a google API key
# # # !!! do not upload to GitHub
# # gmaps = googlemaps.Client(key='XXX')
# # 
# # real_distance = []
# # 
# # for i in range(len(need_transportation)):
# #     try:
# #         print(i)
# #         real_distance.append(gmaps.distance_matrix(need_transportation['WRRF_location'].iloc[i], center_WRRF_location,
# #                                                    mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000)
# #     except KeyError:
# #         print('--------------------------------')
# #         real_distance.append(np.nan)
# # 
# # need_transportation['real_distance_to_center_WRRF_km'] = real_distance
# # 
# # distance_inventory_satellite = need_transportation[['CWNS','linear_distance_to_center_WRRF_km','real_distance_to_center_WRRF_km']]
# # 
# # distance_inventory_satellite.to_excel(folder + 'satellite_distance_inventory.xlsx')
# # =============================================================================

# distance_inventory_satellite = pd.read_excel(folder + 'satellite_distance_inventory.xlsx')

# # match using WRRF ID ('CWNS')
# need_transportation = need_transportation.merge(distance_inventory_satellite, how='left', on=['CWNS'])

# missing_satellite_distance = []
# for i in need_transportation.index:
#     if pd.isna(need_transportation.loc[i,'linear_distance_to_center_WRRF_km_y']):
#         missing_satellite_distance.append(i)

# if len(missing_satellite_distance) == 0:
#     need_transportation.to_excel(folder + f'results/satellite/satellite_WRRFs_{date.today()}.xlsx')
# else:
#     # !!! get a google API key
#     # !!! do not upload to GitHub
#     gmaps = googlemaps.Client(key='XXX')
    
#     for i in missing_satellite_distance:
#         need_transportation.loc[i,'linear_distance_to_center_WRRF_km_test'] = geopy.distance.geodesic(need_transportation.loc[i,'WRRF_location'],
#                                                                                                       center_WRRF_data_preparation.loc[i,'WRRF_location']).km
        
#         try:
#             print(i)
#             need_transportation.loc[i,'real_distance_to_center_WRRF_km'] = gmaps.distance_matrix(need_transportation.loc[i,'WRRF_location'],
#                                                                                                  center_WRRF_data_preparation.loc[i,'WRRF_location'],
#                                                                                                  mode='driving')['rows'][0]['elements'][0]['distance']['value']/1000
#         except KeyError:
#             print('--------------------------------')
#             need_transportation.loc[i,'real_distance_to_center_WRRF_km'] = np.nan
    
#     distance_inventory_satellite = need_transportation[['CWNS','linear_distance_to_center_WRRF_km_test','real_distance_to_center_WRRF_km']]
    
#     distance_inventory_satellite.to_excel(folder + 'satellite_distance_inventory.xlsx')
    
#     need_transportation.to_excel(folder + f'results/satellite/satellite_WRRFs_{date.today()}.xlsx')

#%% sludge transportation (satellite, SA) (center WRRF uncertainty) (superseded)

# filterwarnings('ignore')

# # !!! update the file here if necessary
# satellite_result = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
# satellite_result = satellite_result[satellite_result['USD_decarbonization'].notna()]
# satellite_result = satellite_result[satellite_result['USD_decarbonization'] <= 0]

# center_WRRF = satellite_result[satellite_result['facility']=='MWRDGC, Stickney Treatment Plant']

# # $/tonne
# center_WRRF['waste_cost'] = sum(center_WRRF[i][0]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/center_WRRF['total_sludge_amount_kg_per_year'][0]*1000
# # kg CO2 eq/tonne
# center_WRRF['waste_GHG'] =  sum(center_WRRF[i][0]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/center_WRRF['total_sludge_amount_kg_per_year'][0]*1000

# sys = create_geospatial_system(size=center_WRRF['total_sludge_amount_kg_per_year'][0]/1000/365,
#                                sludge_transportation=0,
#                                sludge_distance=100,
#                                biocrude_distance=center_WRRF['real_distance_km'][0],
#                                anaerobic_digestion=center_WRRF['sludge_anaerobic_digestion'][0],
#                                aerobic_digestion=center_WRRF['sludge_aerobic_digestion'][0],
#                                ww_2_dry_sludge_ratio=1,
#                                state=center_WRRF['state'][0],
#                                nitrogen_fertilizer=center_WRRF['nitrogen_fertilizer'][0],
#                                elec_GHG=center_WRRF['kg_CO2e_kWh'][0],
#                                wage_adjustment=center_WRRF['wage_quotient'][0]/100)

# # if AD and AeD both exist, assume AeD (due to higher ash content) to be conservative
# if center_WRRF['sludge_aerobic_digestion'][0] == 1:
#     sludge_ash_values=[0.349, 0.436, 0.523, 'aerobic_digestion']
#     sludge_lipid_values=[0.154, 0.193, 0.232]
#     sludge_protein_values=[0.408, 0.510, 0.612]
# elif center_WRRF['sludge_anaerobic_digestion'][0] == 1:
#     sludge_ash_values=[0.331, 0.414, 0.497, 'anaerobic_digestion']
#     sludge_lipid_values=[0.154, 0.193, 0.232]
#     sludge_protein_values=[0.408, 0.510, 0.612]
# else:
#     sludge_ash_values=[0.174, 0.231, 0.308, 'no_digestion']
#     sludge_lipid_values=[0.080, 0.206, 0.308]
#     sludge_protein_values=[0.380, 0.456, 0.485]

# crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']==center_WRRF['state'][0]]['2022_min'].iloc[0],
#                           crude_oil_price_data[crude_oil_price_data['state']==center_WRRF['state'][0]]['2022_average'].iloc[0],
#                           crude_oil_price_data[crude_oil_price_data['state']==center_WRRF['state'][0]]['2022_max'].iloc[0]]

# try:
#     DAP_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['DAP']
# except KeyError:
#     DAP_price_values = [588, 984.5, 1381]

# try:
#     anhydrous_ammonia_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['anhydrous_ammonia']
# except KeyError:
#     anhydrous_ammonia_price_values = [1065, 1407.5, 1750]

# try:
#     urea_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['urea']
# except KeyError:
#     urea_price_values = [447, 855, 1263]

# try:
#     UAN_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['UAN']
# except KeyError:
#     UAN_price_values = [350, 605, 860]

# model = create_geospatial_model(system=sys,
#                                 sludge_ash=sludge_ash_values,
#                                 sludge_lipid=sludge_lipid_values,
#                                 sludge_protein=sludge_protein_values,
#                                 crude_oil_price=crude_oil_price_values,
#                                 DAP_price=DAP_price_values,
#                                 anhydrous_ammonia_price=anhydrous_ammonia_price_values,
#                                 urea_price=urea_price_values,
#                                 UAN_price=UAN_price_values)

# kwargs = {'N':1000,'rule':'L','seed':3221}
# samples = model.sample(**kwargs)
# model.load_samples(samples)
# model.evaluate()

# idx = len(model.parameters)
# parameters = model.table.iloc[:, :idx]
# results = model.table.iloc[:, idx:]

# # $/day
# saving_center = (center_WRRF['waste_cost'][0] - results[('Geospatial','Sludge management cost [$/tonne dry sludge]')])*center_WRRF['total_sludge_amount_kg_per_year'][0]/1000/365
# # kg CO2 eq/day
# decarbonization_center = (center_WRRF['waste_GHG'][0] - results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')])*center_WRRF['total_sludge_amount_kg_per_year'][0]/1000/365

#%% sludge transportation (satellite, SA) (all WRRFs uncertainty) (superseded)

# filterwarnings('ignore')

# # !!! update the file here if necessary
# satellite_WRRF = pd.read_excel(folder + 'results/satellite/satellite_WRRFs_2025-02-19.xlsx')

# print(len(satellite_WRRF))

# # $/tonne
# satellite_WRRF['waste_cost'] = sum(satellite_WRRF[i]*sludge_disposal_cost[i] for i in sludge_disposal_cost.keys())/satellite_WRRF['total_sludge_amount_kg_per_year']*1000
# # kg CO2 eq/tonne
# satellite_WRRF['waste_GHG'] =  sum(satellite_WRRF[i]*sludge_emission_factor[i] for i in sludge_emission_factor.keys())/satellite_WRRF['total_sludge_amount_kg_per_year']*1000

# # the units here do not matter
# total_cost_except_center = (satellite_WRRF['waste_cost']*satellite_WRRF['total_sludge_amount_kg_per_year']).sum(axis=0)
# total_GHG_except_center = (satellite_WRRF['waste_GHG']*satellite_WRRF['total_sludge_amount_kg_per_year']).sum(axis=0)
# total_sludge_except_center = satellite_WRRF['total_sludge_amount_kg_per_year'].sum(axis=0)

# total_cost = total_cost_except_center + center_WRRF['waste_cost'][0]*center_WRRF['total_sludge_amount_kg_per_year'][0]
# total_GHG = total_GHG_except_center + center_WRRF['waste_GHG'][0]*center_WRRF['total_sludge_amount_kg_per_year'][0]
# total_sludge = total_sludge_except_center + center_WRRF['total_sludge_amount_kg_per_year'][0]
# total_size = total_sludge/1000/365

# average_cost = total_cost/total_sludge
# average_GHG = total_GHG/total_sludge

# distance_to_refinery = center_WRRF['real_distance_km'][0]

# sludge_transporataion_distance = (satellite_WRRF['real_distance_to_center_WRRF_km']*satellite_WRRF['total_sludge_amount_kg_per_year']).sum(axis=0)/total_sludge

# total_ash_except_center = 0
# total_lipid_except_center = 0
# total_protein_except_center = 0

# for i in range(len(satellite_WRRF)):
#     # if AD and AeD both exist, assume AeD (due to higher ash content) to be conservative
#     if satellite_WRRF.iloc[i]['sludge_aerobic_digestion'] == 1:
#         sludge_ash_value=0.436
#         sludge_lipid_value=0.193
#         sludge_protein_value=0.510
#     elif satellite_WRRF.iloc[i]['sludge_anaerobic_digestion'] == 1:
#         sludge_ash_value=0.414
#         sludge_lipid_value=0.193
#         sludge_protein_value=0.510
#     else:
#         sludge_ash_value=0.231
#         sludge_lipid_value=0.206
#         sludge_protein_value=0.456

#     total_ash_except_center += sludge_ash_value*satellite_WRRF.iloc[i]['total_sludge_amount_kg_per_year']
#     total_lipid_except_center += sludge_lipid_value*satellite_WRRF.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)
#     total_protein_except_center += sludge_protein_value*satellite_WRRF.iloc[i]['total_sludge_amount_kg_per_year']*(1-sludge_ash_value)

# # the center WRRF has anaerobic digestion
# assert center_WRRF.sludge_aerobic_digestion[0] == 0
# assert center_WRRF.sludge_anaerobic_digestion[0] == 1

# total_ash = total_ash_except_center + center_WRRF['total_sludge_amount_kg_per_year']*0.414
# total_lipid = total_lipid_except_center + center_WRRF['total_sludge_amount_kg_per_year']*(1-0.414)*0.193
# total_protein = total_protein_except_center + center_WRRF['total_sludge_amount_kg_per_year']*(1-0.414)*0.510

# average_ash_dw = total_ash/total_sludge
# average_lipid_afdw = total_lipid/(total_sludge-total_ash)
# average_protein_afdw = total_protein/(total_sludge-total_ash)

# # !!! use uniform distribution (0.8x, 1x, 1.2x) for combined WRRFs
# # set sludge_ash_values[-1] = 'digestion', this does not necessarily mean there is digestion, but this will set the uncertainty distribution of ash, lipid, and protein to uniform in the model
# sludge_ash_values = [average_ash_dw*0.8, average_ash_dw, average_ash_dw*1.2, 'digestion']
# sludge_lipid_values = [average_lipid_afdw*0.8, average_lipid_afdw, average_lipid_afdw*1.2]
# sludge_protein_values = [average_protein_afdw*0.8, average_protein_afdw, average_protein_afdw*1.2]

# sys = create_geospatial_system(size=total_size,
#                                sludge_transportation=1,
#                                sludge_distance=sludge_transporataion_distance,
#                                biocrude_distance=distance_to_refinery,
#                                average_sludge_dw_ash=sludge_ash_values[1][0],
#                                average_sludge_afdw_lipid=sludge_lipid_values[1][0],
#                                average_sludge_afdw_protein=sludge_protein_values[1][0],
#                                anaerobic_digestion=None,
#                                aerobic_digestion=None,
#                                ww_2_dry_sludge_ratio=1,
#                                state=center_WRRF['state'][0],
#                                nitrogen_fertilizer=center_WRRF['nitrogen_fertilizer'][0],
#                                elec_GHG=center_WRRF['kg_CO2e_kWh'][0],
#                                wage_adjustment=center_WRRF['wage_quotient'][0]/100)

# crude_oil_price_values = [crude_oil_price_data[crude_oil_price_data['state']==center_WRRF['state'][0]]['2022_min'].iloc[0],
#                           crude_oil_price_data[crude_oil_price_data['state']==center_WRRF['state'][0]]['2022_average'].iloc[0],
#                           crude_oil_price_data[crude_oil_price_data['state']==center_WRRF['state'][0]]['2022_max'].iloc[0]]

# try:
#     DAP_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['DAP']
# except KeyError:
#     DAP_price_values = [588, 984.5, 1381]

# try:
#     anhydrous_ammonia_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['anhydrous_ammonia']
# except KeyError:
#     anhydrous_ammonia_price_values = [1065, 1407.5, 1750]

# try:
#     urea_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['urea']
# except KeyError:
#     urea_price_values = [447, 855, 1263]

# try:
#     UAN_price_values = fertilizer_price_uncertainty[center_WRRF['state'][0]]['UAN']
# except KeyError:
#     UAN_price_values = [350, 605, 860]

# model = create_geospatial_model(system=sys,
#                                 sludge_ash=sludge_ash_values,
#                                 sludge_lipid=sludge_lipid_values,
#                                 sludge_protein=sludge_protein_values,
#                                 crude_oil_price=crude_oil_price_values,
#                                 DAP_price=DAP_price_values,
#                                 anhydrous_ammonia_price=anhydrous_ammonia_price_values,
#                                 urea_price=urea_price_values,
#                                 UAN_price=UAN_price_values)

# kwargs = {'N':1000,'rule':'L','seed':3221}
# samples = model.sample(**kwargs)
# model.load_samples(samples)
# model.evaluate()

# idx = len(model.parameters)
# parameters = model.table.iloc[:, :idx]
# all_WRRFs_results = model.table.iloc[:, idx:]

# # $/day
# saving_all = (average_cost - all_WRRFs_results[('Geospatial','Sludge management cost [$/tonne dry sludge]')])*total_sludge/1000/365
# # kg CO2 eq/day
# decarbonization_all = (average_GHG - all_WRRFs_results[('Geospatial','Sludge CI [kg CO2 eq/tonne dry sludge]')])*total_sludge/1000/365

# SA_results = pd.DataFrame({'center_saving': saving_center,
#                            'all_saving': saving_all,
#                            'center_decarbonization': decarbonization_center,
#                            'all_decarbonization': decarbonization_all})

# SA_results.to_excel(folder + f'results/satellite/center_vs_all_{date.today()}.xlsx')

#%% sludge transportation (satellite, SA) (visualization) (superseded)

# # !!! update the file here if necessary
# SA_results = pd.read_excel(folder + 'results/satellite/center_vs_all_2025-02-19.xlsx')
# SA_results.drop('Unnamed: 0', axis=1, inplace=True)
# # convert decarbonization from kg CO2 eq/day to tonne CO2 eq/day
# SA_results['center_decarbonization'] /= 1000
# SA_results['all_decarbonization'] /= 1000

# SA_results['center_saving'] /= 1000
# SA_results['all_saving'] /= 1000

# plt.rcParams['axes.linewidth'] = 3
# plt.rcParams['xtick.labelsize'] = 38
# plt.rcParams['ytick.labelsize'] = 38

# plt.rcParams.update({'mathtext.fontset': 'custom'})
# plt.rcParams.update({'mathtext.default': 'regular'})
# plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

# fig, ax = plt.subplots(figsize = (10, 10))

# ax.set_xlim(-50, 70)
# ax.set_ylim(0, 600)

# ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False, pad=6)

# ax_top = ax.twiny()
# ax_top.set_xlim(ax.get_xlim())
# ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

# plt.xticks(np.arange(-50, 90, 20))
# plt.yticks(np.arange(0, 700, 100))

# ax_right = ax.twinx()
# ax_right.set_ylim(ax.get_ylim())
# ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=True, labelcolor='none')

# plt.xticks(np.arange(-50, 90, 20))
# plt.yticks(np.arange(0, 700, 100))

# ax.set_xlabel(r'$\mathbf{Decarbonization\ amount}$' + '\n[tonne CO${_2}$ eq·day${^{-1}}$]', fontname='Arial', fontsize=45, linespacing=0.8)
# ax.set_ylabel(r'$\mathbf{Saving}$ [k-\$·day${^{-1}}$]', fontname='Arial', fontsize=45)

# mathtext.FontConstantsBase.sup1 = 0.35

# def add_rectangle(item, color, edgecolor):
#     rectangle_fill = Rectangle((SA_results[f'{item}_decarbonization'].quantile(0.05), SA_results[f'{item}_saving'].quantile(0.05)),
#                                 SA_results[f'{item}_decarbonization'].quantile(0.95) - SA_results[f'{item}_decarbonization'].quantile(0.05),
#                                 SA_results[f'{item}_saving'].quantile(0.95) - SA_results[f'{item}_saving'].quantile(0.05),
#                                 fc=color, alpha=0.8)
#     ax.add_patch(rectangle_fill)
    
#     rectangle_edge = Rectangle((SA_results[f'{item}_decarbonization'].quantile(0.05), SA_results[f'{item}_saving'].quantile(0.05)),
#                                 SA_results[f'{item}_decarbonization'].quantile(0.95) - SA_results[f'{item}_decarbonization'].quantile(0.05),
#                                 SA_results[f'{item}_saving'].quantile(0.95) - SA_results[f'{item}_saving'].quantile(0.05),
#                                 color=edgecolor, lw=3, fc='none', alpha=1)
#     ax.add_patch(rectangle_edge)
    
# def add_line(item, color):
#     ax.plot([SA_results[f'{item}_decarbonization'].quantile(0.25), SA_results[f'{item}_decarbonization'].quantile(0.75)],
#              [SA_results[f'{item}_saving'].quantile(0.5), SA_results[f'{item}_saving'].quantile(0.5)],
#              lw=3, color=color, solid_capstyle='round')
#     ax.plot([SA_results[f'{item}_decarbonization'].quantile(0.5), SA_results[f'{item}_decarbonization'].quantile(0.5)],
#              [SA_results[f'{item}_saving'].quantile(0.25), SA_results[f'{item}_saving'].quantile(0.75)],
#              lw=3, color=color, solid_capstyle='round')

# def add_point(item, edgecolor):
#     ax_right.scatter(x=SA_results[f'{item}_decarbonization'].quantile(0.5),
#                      y=SA_results[f'{item}_saving'].quantile(0.5),
#                      marker='s',
#                      s=200,
#                      c='w',
#                      linewidths=3,
#                      alpha=1,
#                      edgecolor=edgecolor)
    
# add_rectangle('center', o, do)
# add_line('center', do)
# add_point('center', do)

# add_rectangle('all', b, db)
# add_line('all', db)
# add_point('all', db)

#%% writing results - part 1 (superseded)

# # !!! update the file here if necessary
# WRRF_finder = pd.read_excel(folder + 'results/baseline/integrated_baseline_2025-02-16.xlsx')
# WRRF_finder = WRRF_finder[WRRF_finder['USD_decarbonization'].notna()]
# WRRF_finder = WRRF_finder[WRRF_finder['USD_decarbonization'] <= 0]

# digestion_WRRFs = WRRF_finder[(WRRF_finder['sludge_anaerobic_digestion'] == 1) | (WRRF_finder['sludge_aerobic_digestion'] == 1)].copy()
# print(digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD_final']])
# print(digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD_final']])

# no_digestion_WRRFs = WRRF_finder[(WRRF_finder['sludge_anaerobic_digestion'] == 0) & (WRRF_finder['sludge_aerobic_digestion'] == 0)].copy()
# print(no_digestion_WRRFs.sort_values('CO2_reduction', ascending=False).iloc[0,][['facility','city','flow_2022_MGD_final']])
# print(no_digestion_WRRFs.sort_values('WRRF_CO2_reduction_ratio', ascending=False).iloc[0,][['facility','city','flow_2022_MGD_final']])