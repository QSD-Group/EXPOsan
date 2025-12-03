#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 08:54:11 2024

@author: jiananfeng
"""

#%% initialization

import pandas as pd, geopandas as gpd, numpy as np, matplotlib.pyplot as plt, matplotlib.colors as colors
from colorpalette import Color
from qsdsan.utils import palettes

gallon_to_liter = 3.78541

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

US_county = gpd.read_file('/Users/jiananfeng/Desktop/PhD_CEE/NSF/HTL_geospatial/cb_2018_us_county_500k/cb_2018_us_county_500k.shp')

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

US_county['STCOFIPS'] = US_county['STATEFP'] + US_county['COUNTYFP']
US_county['STCOFIPS'] = US_county['STCOFIPS'].astype('int64')

US_county = US_county.to_crs(crs='EPSG:3857')

farm_fertilizer = pd.read_excel('/Users/jiananfeng/Desktop/PhD_CEE/NSF/HTL_geospatial/N-P_from_fertilizer_1950-2017-july23-2020.xlsx','farm')
nonfarm_fertilizer = pd.read_excel('/Users/jiananfeng/Desktop/PhD_CEE/NSF/HTL_geospatial/N-P_from_fertilizer_1950-2017-july23-2020.xlsx','nonfarm')

N_farm = farm_fertilizer[['STCOFIPS','farmfertN-kg-2017']]
N_nonfarm = nonfarm_fertilizer[['STCOFIPS','nonffertN-kg-2017']]

N = N_farm.merge(N_nonfarm, on='STCOFIPS', how='inner')
N['total'] = N['farmfertN-kg-2017'] + N['nonffertN-kg-2017']

P_farm = farm_fertilizer[['STCOFIPS','farmfertP-kg-2017']]
P_nonfarm = nonfarm_fertilizer[['STCOFIPS','nonffertP-kg-2017']]

P = P_farm.merge(P_nonfarm, on='STCOFIPS', how='inner')
P['total'] = P['farmfertP-kg-2017'] + P['nonffertP-kg-2017']

#%% N needs

US_county_N = US_county.merge(N, on='STCOFIPS', how='left')
US_county_N = US_county_N[['STCOFIPS','NAME','STATE','total','geometry']]

US_county_N['total'] = US_county_N['total'].fillna(0)

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', b, db])

fig, ax = plt.subplots(figsize=(30, 30))

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 35
plt.rcParams['ytick.labelsize'] = 35
plt.rcParams['font.sans-serif'] = 'Arial'

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

ax = plt.gca()

US_county_N.plot(ax=ax, column='total', legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=0.5)

fig.axes[1].set_yticks(np.arange(0, 70000000, 10000000))
fig.axes[1].set_ylabel('$\mathbf{Nitrogen}$ [kg]', fontname='Arial', fontsize=41)
fig.axes[1].tick_params(length=7.5, width=1.5)

ax.set_aspect(1)

ax.set_axis_off()

#%% P needs

US_county_P = US_county.merge(P, on='STCOFIPS', how='left')
US_county_P = US_county_P[['STCOFIPS','NAME','STATE','total','geometry']]

US_county_P['total'] = US_county_P['total'].fillna(0)

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', r, dr])

fig, ax = plt.subplots(figsize=(30, 30))

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 35
plt.rcParams['ytick.labelsize'] = 35
plt.rcParams['font.sans-serif'] = 'Arial'

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})
plt.rcParams.update({'figure.max_open_warning': 100})

ax = plt.gca()

US_county_P.plot(ax=ax, column='total', legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, edgecolor='k', linewidth=0.5)

fig.axes[1].set_yticks(np.arange(0, 10000000, 1000000))
fig.axes[1].set_ylabel('$\mathbf{Phosphorus}$ [kg]', fontname='Arial', fontsize=41)
fig.axes[1].tick_params(length=7.5, width=1.5)

ax.set_aspect(1)

ax.set_axis_off()