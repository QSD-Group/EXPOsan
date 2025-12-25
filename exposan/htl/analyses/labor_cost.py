#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

#%% initialization

import os, pandas as pd, geopandas as gpd
from colorpalette import Color
from qsdsan.utils import palettes

folder = os.path.dirname(os.path.dirname(__file__))

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

US_county = gpd.read_file(folder + '/data/cb_2018_us_county_500k/cb_2018_us_county_500k.shp')

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

US_county['Area\nCode'] = US_county['STATEFP'] + US_county['COUNTYFP']

US_county = US_county.to_crs(crs='EPSG:3857')

labor_cost = pd.read_excel(folder + '/data/county_labor_cost_2022.xlsx','US_St_Cn_MSA')

US_average_labor_cost = labor_cost[(labor_cost['St'] == 'US') & (labor_cost['Ownership'] == 'Local Government')]['Annual Average Pay']

state_labor_cost = labor_cost[(labor_cost['Cnty'] == 0) & (labor_cost['Ownership'] == 'Local Government')][['St','Annual Average Pay']]
state_labor_cost.rename(columns={'Annual Average Pay': 'State Annual Average Pay'}, inplace=True)

labor_cost = labor_cost[labor_cost['St'].isin(['01','04','05','06','08','09','10','11',
                                               '12','13','16','17','18','19','20','21',
                                               '22','23','24','25','26','27','28','29',
                                               '30','31','32','33','34','35','36','37',
                                               '38','39','40','41','42','44','45','46',
                                               '47','48','49','50','51','53','54','55','56'])]
labor_cost = labor_cost[labor_cost['Ownership'] == 'Local Government']

#%% labor wage

US_county_labor_cost = US_county.merge(labor_cost, on=['Area\nCode'], how='inner')

US_county_labor_cost = US_county_labor_cost.merge(state_labor_cost, on='St', how='left')

US_county_labor_cost.loc[US_county_labor_cost['Annual Average Pay'] == 0, 'Annual Average Pay'] = US_county_labor_cost['State Annual Average Pay']

US_county_labor_cost['quotient'] = US_county_labor_cost['Annual Average Pay']/float(US_average_labor_cost.iloc[0])*100

US_county_labor_cost = US_county_labor_cost[['Area\nCode','NAME','STATE','quotient','geometry']]

US_county_labor_cost.to_file(folder + '/data/county_labor_cost_2022_processed.geojson')