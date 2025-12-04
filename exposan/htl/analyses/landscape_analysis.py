#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# 'country' refers to both countries and territories

#%% initialization

import numpy as np, pandas as pd, geopandas as gpd, matplotlib.pyplot as plt, matplotlib.colors as colors, matplotlib.ticker as mtick, chaospy as cp, json
from colorpalette import Color
from scipy.stats import spearmanr, linregress
from matplotlib.colors import to_hex
from matplotlib.mathtext import _mathtext as mathtext
from matplotlib.gridspec import GridSpec
from chaospy import distributions as shape
from warnings import filterwarnings
from datetime import date
from qsdsan.utils import auom
from exposan.htl import (create_C1_system, create_C2_system, create_C3_system,
                         create_C4_system, create_C5_system, create_C6_system,
                         create_C7_system, create_C8_system, create_C9_system,
                         create_C10_system, create_C11_system, create_C12_system,
                         create_C13_system, create_C14_system, create_C15_system,
                         create_C16_system, create_C17_system, create_C18_system,
                         create_C19_system, create_C20_system, create_C21_system,
                         create_C22_system, create_C23_system, create_C24_system,
                         create_C25_system, create_T1_system, create_T2_system,
                         create_T3_system, create_T4_system, create_T5_system,
                         create_T6_system, create_T7_system, create_T8_system,
                         create_T9_system, create_T10_system, create_T11_system,
                         create_T12_system, create_T13_system, create_T14_system,
                         create_T15_system)

folder = '/Users/jiananfeng/Desktop/UIUC_PhD/PhD_CEE/NSF/HTL_landscape/'

_ton_to_tonne = auom('ton').conversion_factor('tonne')
_MMgal_to_m3 = auom('gal').conversion_factor('m3')*1000000
# 1 ton dry WWRS/MG wastewater, MOP 8, converted to tonne dry WWRS/MG
ww_2_dry_solids = 1*_ton_to_tonne

mathtext.FontConstantsBase.sup1 = 0.35

# color palette
b = Color('blue', (96, 193, 207)).HEX
g = Color('green', (121, 191, 130)).HEX
r = Color('red', (237, 88, 111)).HEX
o = Color('orange', (249, 143, 96)).HEX
y = Color('yellow', (243, 195, 84)).HEX
a = Color('gray', (144, 145, 142)).HEX
p = Color('purple', (162, 128, 185)).HEX

lb = to_hex((96/256*0.5 + 1*0.5, 193/256*0.5 + 1*0.5, 207/256*0.5 + 1*0.5))
lg = to_hex((121/256*0.5 + 1*0.5, 191/256*0.5 + 1*0.5, 130/256*0.5 + 1*0.5))
lr = to_hex((237/256*0.5 + 1*0.5, 88/256*0.5 + 1*0.5, 111/256*0.5 + 1*0.5))
lo = to_hex((249/256*0.5 + 1*0.5, 143/256*0.5 + 1*0.5, 96/256*0.5 + 1*0.5))
ly = to_hex((243/256*0.5 + 1*0.5, 195/256*0.5 + 1*0.5, 84/256*0.5 + 1*0.5))
la = to_hex((144/256*0.5 + 1*0.5, 145/256*0.5 + 1*0.5, 142/256*0.5 + 1*0.5))
lp = to_hex((162/256*0.5 + 1*0.5, 128/256*0.5 + 1*0.5, 185/256*0.5 + 1*0.5))

db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

# =============================================================================
# word
# =============================================================================
world = gpd.read_file(folder + 'analyses/World Bank Official Boundaries - Admin 0_all_layers/WB_GAD_ADM0_complete.shp')

# =============================================================================
# global north / south countries and regions
# =============================================================================
# this classification may not be accurate for all countries but it is accurate for countries included in the figures
global_north_countries_regions = ['United States','Russia','Japan','Turkey','Germany',
                                  'United Kingdom','France','Italy','South Korea',
                                  'Spain','Canada','Ukraine','Poland','Uzbekistan',
                                  'Australia','Kazakhstan','Romania','Netherlands',
                                  'Belgium','Sweden','Portugal','Azerbaijan','Greece',
                                  'Hungary','Israel','Austria','Belarus','Switzerland',
                                  'Kyrgyzstan','Bulgaria','Serbia','Denmark','Finland',
                                  'Norway','Slovakia','Ireland','New Zealand','Croatia',
                                  'Georgia','Moldova','Armenia','Lithuania','Albania',
                                  'Slovenia','Latvia','North Macedonia','Cyprus',
                                  'Estonia','Luxembourg','Montenegro','Malta','Iceland',
                                  'Andorra','Liechtenstein','Monaco','San Marino',
                                  'Palau','Tuvalu','Vatican City']

# this classification may not be accurate for all countries but it is accurate for countries included in the figures
global_south_countries_regions = ['India','China','Mexico','Indonesia','Pakistan',
                                  'Nigeria','Brazil','Bangladesh','Ethiopia','Egypt',
                                  'Philippines','DR Congo','Vietnam','Iran','Thailand',
                                  'Tanzania','South Africa','Kenya','Myanmar','Colombia',
                                  'Sudan','Uganda','Algeria','Iraq','Argentina',
                                  'Afghanistan','Yemen','Angola','Morocco','Malaysia',
                                  'Mozambique','Ghana','Peru','Saudi Arabia','Madagascar',
                                  'Ivory Coast','Cameroon','Nepal','Venezuela',
                                  'Niger','North Korea','Syria','Mali','Burkina Faso',
                                  'Sri Lanka','Malawi','Zambia','Chad','Chile',
                                  'Somalia','Senegal','Guatemala','Ecuador','Cambodia',
                                  'Zimbabwe','Guinea','Benin','Rwanda','Burundi',
                                  'Bolivia','Tunisia','South Sudan','Haiti','Jordan',
                                  'Dominican Republic','United Arab Emirates','Honduras',
                                  'Cuba','Tajikistan','Papua New Guinea','Togo',
                                  'Sierra Leone','Laos','Turkmenistan','Libya',
                                  'Paraguay','Nicaragua','Republic of the Congo',
                                  'El Salvador','Singapore','Lebanon','Liberia',
                                  'Palestine','Central African Republic','Oman',
                                  'Mauritania','Costa Rica','Kuwait','Panama','Eritrea',
                                  'Mongolia','Uruguay','Bosnia and Herzegovina',
                                  'Qatar','Namibia','Jamaica','Gambia','Gabon',
                                  'Botswana','Lesotho','Guinea-Bissau','Equatorial Guinea',
                                  'Bahrain','Trinidad and Tobago','Timor-Leste',
                                  'Mauritius','Eswatini','Djibouti','Fiji','Comoros',
                                  'Solomon Islands','Guyana','Bhutan','Suriname',
                                  'Maldives','Cape Verde','Brunei','Belize','Bahamas',
                                  'Vanuatu','Barbados','Sao Tome and Principe',
                                  'Samoa','Saint Lucia','Kiribati','Seychelles',
                                  'Grenada','Micronesia','Tonga','Saint Vincent and the Grenadines',
                                  'Antigua and Barbuda','Dominica','Saint Kitts and Nevis',
                                  'Marshall Islands','Nauru']

# =============================================================================
# WRRF
# =============================================================================
WRRF = pd.read_csv(folder + 'analyses/HydroWASTE_v10/HydroWASTE_v10.csv', encoding='latin-1')
# WASTE_DIS in m3/day
WRRF = WRRF[WRRF['WASTE_DIS'] != 0]
WRRF = WRRF[~WRRF['STATUS'].isin(['Closed','Decommissioned','Non-Operational'])]
WRRF['dry_solids_tonne_per_day'] = WRRF['WASTE_DIS']/_MMgal_to_m3*ww_2_dry_solids
# only keep WRRFs with ≥ 5 dry tonne WWRS/day for more feasible thermochemical technology deployment (5 dry tonne/day is from the HTL geospatial paper)
WRRF_filtered = WRRF[WRRF['dry_solids_tonne_per_day'] >= 5]
WRRF_filtered.reset_index(inplace=True)
print(f"{WRRF_filtered['dry_solids_tonne_per_day'].sum()/WRRF['dry_solids_tonne_per_day'].sum()*100:.1f}% global WWRS are included." )
print(f"{len(WRRF_filtered)} WRRFs are included." )
print(f"{len(set(WRRF_filtered['COUNTRY']))} countires or regions ('Taiwan' listed separately) are included." )

WRRF_filtered = gpd.GeoDataFrame(WRRF_filtered, crs='EPSG:4269',
                                 geometry=gpd.points_from_xy(x=WRRF_filtered.LON_WWTP,
                                                             y=WRRF_filtered.LAT_WWTP))

HDI = pd.read_excel(folder + 'analyses/HDR25_Statistical_Annex_HDI_Table.xlsx')

# =============================================================================
# electricity price
# =============================================================================
electricity_price = pd.read_excel(folder + 'analyses/P_Electric Prices by Country.xlsx')
electricity_price.rename(columns={'Country Name':'country',
                                  'Country Code':'country_code',
                                  'Time':'year',
                                  'Getting electricity: Price of electricity (US cents per kWh) (DB16-20 methodology) [IC.ELC.PRI.KH.DB1619]':'US_cents_per_kWh'},
                         inplace=True)
electricity_price = electricity_price[electricity_price['US_cents_per_kWh'] != '..']
electricity_price = electricity_price[electricity_price['US_cents_per_kWh'].notna()]
electricity_price = electricity_price.loc[electricity_price.groupby('country')['year'].idxmax()].reset_index(drop=True)
electricity_price = electricity_price[['country','country_code','year','US_cents_per_kWh']]

# =============================================================================
# electricity CI
# =============================================================================
electriicty_CI = pd.read_excel(folder + 'analyses/low_voltage_electricity_CI.xlsx')

# =============================================================================
# labor cost
# =============================================================================
labor_cost = pd.read_excel(folder + 'analyses/income.xlsx')

# =============================================================================
# price level index (PLI)
# =============================================================================
PLI = pd.read_excel(folder + 'analyses/P_Data_Extract_From_World_Development_Indicators.xlsx')
year_cols = ['1990 [YR1990]','2000 [YR2000]','2015 [YR2015]','2016 [YR2016]','2017 [YR2017]',
             '2018 [YR2018]','2019 [YR2019]','2020 [YR2020]','2021 [YR2021]','2022 [YR2022]',
             '2023 [YR2023]','2024 [YR2024]']
PLI[year_cols] = PLI[year_cols].replace('..', np.nan)
PLI['PLI'] = PLI[year_cols].ffill(axis=1).iloc[:, -1]
PLI = PLI[PLI['PLI'].notna()]
PLI.reset_index(inplace=True)
PLI = PLI[['Country Name','Country Code','PLI']]
PLI.rename(columns={'Country Name':'country',
                    'Country Code':'country_code'},
           inplace=True)

#%% data processing

# # =============================================================================
# # electricity price
# # =============================================================================
# for i in set(WRRF_filtered['CNTRY_ISO']):
#     if i not in list(electricity_price['country_code']):
#         print(set(WRRF_filtered[WRRF_filtered['CNTRY_ISO'] == i]['COUNTRY']))
#         print(i)

# # https://www.voronoiapp.com/energy/Whats-the-Average-Cost-of-1-kWh-Electricity-around-the-World--3398
# electricity_price.loc[len(electricity_price)] = ['Monaco', 'MCO', 9999, 18.0]
# electricity_price.loc[len(electricity_price)] = ['British Virgin Islands', 'VGB', 9999, 22.5]
# electricity_price.loc[len(electricity_price)] = ['Aruba', 'ABW', 9999, 17.8]
# electricity_price.loc[len(electricity_price)] = ['Cuba', 'CUB', 9999, 38.6]
# electricity_price.loc[len(electricity_price)] = ['French Polynesia', 'PYF', 9999, 25.9]
# electricity_price.loc[len(electricity_price)] = ['Curaçao', 'CUW', 9999, 41.9]
# # https://www.globalpetrolprices.com/Macao/electricity_prices/
# electricity_price.loc[len(electricity_price)] = ['Macao', 'MAC', 2025, 15.9]
# # https://www.nkeconwatch.com/category/communications/page/6/ (6.5 2012 euro cent, converted to 2012 US cent)
# electricity_price.loc[len(electricity_price)] = ['North Korea', 'PRK', 2012, 8.4]
# # https://www.ceicdata.com/en/turkmenistan/environmental-environmental-policy-taxes-and-transfers-non-oecd-member-annual/tm-industry-electricity-price-usd-per-kwh
# electricity_price.loc[len(electricity_price)] = ['Turkmenistan', 'TKM', 2020, 1.5]

# electricity_price.drop_duplicates(inplace=True)

# assert len([i for i in set(WRRF_filtered['CNTRY_ISO']) if i not in list(electricity_price['country_code'])]) == 0

# electricity_price.to_excel(folder + f'analyses/electricity_price_{date.today()}.xlsx')

# # =============================================================================
# # electricity CI
# # =============================================================================
# for i in set(WRRF_filtered['CNTRY_ISO']):
#     if i not in list(electriicty_CI['country_ISO_A3']):
#         print(set(WRRF_filtered[WRRF_filtered['CNTRY_ISO'] == i]['COUNTRY']))
#         print(i)

# # Europe without Switzerland
# electriicty_CI.loc[len(electriicty_CI)] = ['Monaco', 'MC', 'MCO', 0.335772839657998]
# # Africa
# electriicty_CI.loc[len(electriicty_CI)] = ['Equatorial Guinea', 'GQ', 'GNQ', 0.84866548695154]
# electriicty_CI.loc[len(electriicty_CI)] = ['Mali', 'ML', 'MLI', 0.84866548695154]
# electriicty_CI.loc[len(electriicty_CI)] = ['Guinea', 'GN', 'GIN', 0.84866548695154]
# electriicty_CI.loc[len(electriicty_CI)] = ['Seychelles', 'SC', 'SYC', 0.84866548695154]
# # Latin America and the Caribbean
# electriicty_CI.loc[len(electriicty_CI)] = ['British Virgin Islands', 'VG', 'VGB', 0.35047605957495]
# electriicty_CI.loc[len(electriicty_CI)] = ['Aruba', 'AW', 'ABW', 0.35047605957495]
# electriicty_CI.loc[len(electriicty_CI)] = ['Bahamas', 'BS', 'BHS', 0.35047605957495]
# electriicty_CI.loc[len(electriicty_CI)] = ['Barbados', 'BB', 'BRB', 0.35047605957495]
# electriicty_CI.loc[len(electriicty_CI)] = ['Antigua and Barbuda', 'AG', 'ATG', 0.35047605957495]
# electriicty_CI.loc[len(electriicty_CI)] = ['Belize', 'BZ', 'BLZ', 0.35047605957495]

# # Middle East
# electriicty_CI.loc[len(electriicty_CI)] = ['Palestina', 'PS', 'PSE', 0.907394533030651]

# # Asia
# electriicty_CI.loc[len(electriicty_CI)] = ['Macao', 'MO', 'MAC', 0.876728132896613]
# electriicty_CI.loc[len(electriicty_CI)] = ['Laos', 'LA', 'LAO', 0.876728132896613]
# electriicty_CI.loc[len(electriicty_CI)] = ['Afghanistan', 'AF', 'AFG', 0.876728132896613]

# # Global
# electriicty_CI.loc[len(electriicty_CI)] = ['Papua New Guinea', 'PG', 'PNG', 0.691007559959689]
# electriicty_CI.loc[len(electriicty_CI)] = ['Fiji', 'FJ', 'FJI', 0.691007559959689]
# electriicty_CI.loc[len(electriicty_CI)] = ['French Polynesia', 'PF', 'PYF', 0.691007559959689]

# electriicty_CI.drop_duplicates(inplace=True)

# assert len([i for i in set(WRRF_filtered['CNTRY_ISO']) if i not in list(electriicty_CI['country_ISO_A3'])]) == 0

# electriicty_CI.to_excel(folder + f'analyses/electricity_CI_{date.today()}.xlsx')

# # =============================================================================
# # labor cost
# # =============================================================================
# for i in set(WRRF_filtered['CNTRY_ISO']):
#     if i not in list(labor_cost['country_code']):
#         print(set(WRRF_filtered[WRRF_filtered['CNTRY_ISO'] == i]['COUNTRY']))
#         print(i)

# # use the minimum value in the available dataset
# labor_cost.loc[len(labor_cost)] = ['North Korea', 'PRK', 190]

# labor_cost.drop_duplicates(inplace=True)

# assert len([i for i in set(WRRF_filtered['CNTRY_ISO']) if i not in list(labor_cost['country_code'])]) == 0

# labor_cost.to_excel(folder + f'analyses/labor_cost_{date.today()}.xlsx')

# # =============================================================================
# # PLI
# # =============================================================================
# for i in set(WRRF_filtered['CNTRY_ISO']):
#     if i not in list(PLI['country_code']):
#         print(set(WRRF_filtered[WRRF_filtered['CNTRY_ISO'] == i]['COUNTRY']))
#         print(i)

# # use France data
# PLI.loc[len(PLI)] = ['Monaco', 'MCO', 0.753]
# # use the average of China, Hong Kong SAR, China, Macao SAR, China, South Korea, and Japan
# PLI.loc[len(PLI)] = ['Taiwan', 'TWN', 0.604]
# # use the average of Caribbean countries (Dominica, Haiti, and Trinidad and Tobago)
# PLI.loc[len(PLI)] = ['British Virgin Islands', 'VGB', 0.566]
# # use the average of Latin American countries (Belize, Costa Rica, El Salvador, Guatemala, Honduras,
# #                                              Mexico, Nicaragua, Panama, Argentina, Bolivia, Brazil,
# #                                              Chile, Colombia, Ecuador, Guyana, Paraguay, Peru, Suriname,
# #                                              Uruguay, Dominica, Dominican Republic, Haiti, Trinidad and Tobago)
# PLI.loc[len(PLI)] = ['Cuba', 'CUB', 0.464]
# # use the average of Pacific Islands countries (Fiji, Kiribati, Marshall Islands, Nauru, Palau,
# #                                               Papua New Guinea, Samoa, Solomon Islands, Tonga,
# #                                               Tuvalu, Vanuatu)
# PLI.loc[len(PLI)] = ['French Polynesia', 'PYF', 0.764]
# # use the minimum value in the available dataset
# PLI.loc[len(PLI)] = ['North Korea', 'PRK', 0.125]

# PLI.drop_duplicates(inplace=True)

# assert len([i for i in set(WRRF_filtered['CNTRY_ISO']) if i not in list(PLI['country_code'])]) == 0

# PLI.to_excel(folder + f'analyses/PLI_{date.today()}.xlsx')

#%% labor costs comparison

labor_cost_comparison = pd.read_excel(folder + 'analyses/labor_cost_comparison.xlsx')
labor_cost_comparison = labor_cost_comparison[~labor_cost_comparison['ILO'].isna()]

pearson_r = np.corrcoef(labor_cost_comparison['Worlddata'], labor_cost_comparison['ILO'])[0, 1]
print("Pearson's r:", pearson_r.round(2))

spearman_rho, spearman_p = spearmanr(labor_cost_comparison['Worlddata'], labor_cost_comparison['ILO'])
print("Spearman's rho:", spearman_rho.round(2))

#%% MPs concentration

MPs = pd.read_excel(folder + 'analyses/EC_data.xlsx','MPs_summary')

print('\n' + str([len(MPs[i].dropna()) for i in MPs.columns]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim([10**-7, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{C}$ [particles·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([MPs[i].dropna() for i in MPs.columns],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=dr, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lr, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/MPs_concentration.pdf', transparent=True, bbox_inches='tight')

#%% MPs concentration (unique values)

MPs = pd.read_excel(folder + 'analyses/EC_data.xlsx','MPs_summary')

print('\n' + str([len(MPs[i].dropna().drop_duplicates()) for i in MPs.columns]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim([10**-7, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{MPs}$\n[particles·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([MPs[i].dropna().drop_duplicates() for i in MPs.columns],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=dr, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lr, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/MPs_concentration_unique.pdf', transparent=True, bbox_inches='tight')

#%% MPs capture

textile_mass_flow = shape.Uniform(398000, 597000)
vehicle_tires_mass_flow = shape.Uniform(1130000, 1690000)
city_dust_mass_flow = shape.Uniform(520000, 780000)
road_markings_mass_flow = shape.Uniform(472000, 708000)
personal_care_product_mass_flow = shape.Uniform(20000, 30000)
marine_coatings_mass_flow = shape.Uniform(40000, 60000)

textile_to_WRRF = shape.Uniform(0.414, 0.620)
vehicle_tires_to_runoff = shape.Uniform(0.72, 1)
vehicle_tires_runoff_to_sewer = shape.Uniform(0.548, 0.822)
city_dust_to_runoff = shape.Uniform(0.456, 0.684)
road_markings_to_runoff = shape.Uniform(0.72, 1)
road_markings_runoff_to_sewer = shape.Uniform(0.188, 0.282)
personal_care_product_to_WRRF = shape.Uniform(0.398, 0.596)

combined_sewer_fraction = shape.Uniform(0.254, 0.380)

WWRS_removal = shape.Triangle(0.806, 0.952, 0.998)

joint = cp.distributions.J(textile_mass_flow, vehicle_tires_mass_flow,
                           city_dust_mass_flow, road_markings_mass_flow, personal_care_product_mass_flow,
                           marine_coatings_mass_flow, textile_to_WRRF, vehicle_tires_to_runoff,
                           vehicle_tires_runoff_to_sewer, city_dust_to_runoff, road_markings_to_runoff,
                           road_markings_runoff_to_sewer, personal_care_product_to_WRRF,
                           combined_sewer_fraction, WWRS_removal)

sample = joint.sample(100000)

MPs_WWRS_MC = pd.DataFrame()

MPs_textile = sample[0]*sample[6]
MPs_vehicle_tires = sample[1]*sample[7]*sample[8]*sample[13]
MPs_city_dust = sample[2]*sample[9]*sample[13]
MPs_road_markings = sample[3]*sample[10]*sample[11]*sample[13]
MPs_personal_care_products = sample[4]*sample[12]
MPs_marine_coatings = sample[5]*0
MPs_total = sample[0] + sample[1] + sample[2] + sample[3] + sample[4] + sample[5]


MPs_WWRS_MC = (MPs_textile + MPs_vehicle_tires + MPs_city_dust + MPs_road_markings +\
               MPs_personal_care_products + MPs_marine_coatings)*sample[14]/MPs_total

fig, ax = plt.subplots(figsize=(2.5, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 1])
ax.set_ylim([0, 40])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 50, 10), fontname='Arial')

ax.set_ylabel('$\mathbf{Capture}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 50, 10), fontname='Arial')

ax.bar(0.5,
       np.quantile(MPs_WWRS_MC, 0.5)*100,
       yerr=np.array([[(np.quantile(MPs_WWRS_MC, 0.5) - np.quantile(MPs_WWRS_MC, 0.05))*100], [(np.quantile(MPs_WWRS_MC, 0.95) - np.quantile(MPs_WWRS_MC, 0.5))*100]]),
       error_kw=dict(capsize=15, lw=3, capthick=3),
       width=0.7,
       color=dr,
       edgecolor='k',
       linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/MPs_capture.pdf', transparent=True, bbox_inches='tight')

#%% MPs removal

MPs_removal = pd.read_excel(folder + 'analyses/EC_removal_data.xlsx','MPs_summary')

conventional_MPs_removal = MPs_removal[MPs_removal['process'] == 'C']
wet_thermochemical_MPs_removal = MPs_removal[MPs_removal['process'] == 'wet_T']
dry_thermochemical_MPs_removal = MPs_removal[MPs_removal['process'] == 'dry_T']

fig = plt.figure(figsize=(20.16, 5))
gs = GridSpec(1, 2, width_ratios=[4, 1.25], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 700])
ax.set_ylim([0, 100])
ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.yaxis.set_major_formatter(mtick.PercentFormatter())

plt.xticks(np.arange(0, 800, 100), fontname='Arial')
plt.yticks(np.arange(0, 120, 20), fontname='Arial')

# ax.set_xlabel('$\mathbf{Temperature}$ [°C]',
#               fontname='Arial',
#               fontsize=45,
#               labelpad=13)

ax.set_ylabel('$\mathbf{Removal}$',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 800, 100), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.scatter(conventional_MPs_removal['temperature'], conventional_MPs_removal['reduction']*100, marker='o', color=lr, s=1000, edgecolor='k', linewidth=3)
ax.scatter(wet_thermochemical_MPs_removal['temperature'], wet_thermochemical_MPs_removal['reduction']*100, marker='^', color=r, s=1000, edgecolor='k', linewidth=3)
ax.scatter(dry_thermochemical_MPs_removal['temperature'], dry_thermochemical_MPs_removal['reduction']*100, marker='s', color=dr, s=1000, edgecolor='k', linewidth=3)

ax_box = fig.add_subplot(gs[0, 1], sharey=ax)

data = [conventional_MPs_removal['reduction']*100, wet_thermochemical_MPs_removal['reduction']*100, dry_thermochemical_MPs_removal['reduction']*100]
bp = ax_box.boxplot(data, vert=True, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=lr, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=r, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=dr, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

ax_box.axis('off')

plt.savefig('/Users/jiananfeng/Desktop/MPs_removal.pdf', transparent=True, bbox_inches='tight')

#%% PhACs concentration

PhACs = pd.read_excel(folder + 'analyses/pharms-uba_v3_2021_0.xlsx')
set(PhACs['Matrix'])
assert (PhACs['MEC standardized'].isna()).sum() == 0
set(PhACs['Statistics'])
PhACs = PhACs[PhACs['Statistics'].isin(['Average','Max','Mean','Median','Min',
                                        'Single Value','Single value',
                                        'Single value of a composite sample',
                                        'average','max','median','min','single Value',
                                        'time-weighted-average concentration'])]

PhACs = PhACs.fillna('NaN')

# remove average or mean if min or max is present
group_cols = [c for c in PhACs.columns if c not in ['ID','Statistics','Detection','MEC original','MEC standardized']]

def filter_stats(group):
    stats = set(group['Statistics'])
    if ('Max' in stats) or ('Min' in stats) or ('max' in stats) or ('min' in stats):
        return group[~group['Statistics'].isin(['Average','Mean','average','time-weighted-average concentration'])]
    return group

PhACs = PhACs.groupby(group_cols, group_keys=False).apply(filter_stats)

# WWRS
WWRS_PhACs = PhACs[PhACs['Matrix'].isin(['Dissolved activated sludge',
                                         'WWTP biosolid',
                                         'WWTP dehydrated sludge',
                                         'WWTP digested sludge',
                                         'WWTP primary sludge',
                                         'WWTP secondary sludge',
                                         'WWTP sludge',
                                         'sewage sludge'])]
set(WWRS_PhACs['Unit standard'])
WWRS_PhACs = WWRS_PhACs[WWRS_PhACs['Unit standard'].isin(['mg/kg dry-weight','µg/kg dry-weight'])]
WWRS_PhACs.loc[WWRS_PhACs['Unit standard'] == 'mg/kg dry-weight', 'MEC standardized'] *= 1000

# soil
soil_PhACs = PhACs[PhACs['Matrix'] == 'Soil']
set(soil_PhACs['Unit standard'])
soil_PhACs = soil_PhACs[soil_PhACs['Unit standard'] == 'mg/kg dry-weight']
soil_PhACs['MEC standardized'] *= 1000

# sediment
sediment_PhACs = PhACs[PhACs['Matrix'].isin(['Sediment - Aquaculture',
                                             'Sediment - Estuary',
                                             'Sediment - Lagoon',
                                             'Sediment - Lake',
                                             'Sediment - River/Stream',
                                             'Sediment - Sea or Ocean',
                                             'Sediment - unspecific'])]
set(sediment_PhACs['Unit standard'])
sediment_PhACs = sediment_PhACs[sediment_PhACs['Unit standard'] == 'mg/kg dry-weight']
sediment_PhACs['MEC standardized'] *= 1000

# air
air_PhACs = pd.read_excel(folder + 'analyses/EC_data.xlsx','PhACs_summary')

# water
water_PhACs = PhACs[PhACs['Matrix'].isin(['Drinking Water',
                                          'Groundwater',
                                          'Infiltration water',
                                          'Rain',
                                          'Raw Water - Drinking Water Treatment Plant',
                                          'Reclaimed Water',
                                          'Reservoir drainage',
                                          'Riverbank filtration',
                                          'Soil Water',
                                          'Spring water',
                                          'Surface Water - Aquaculture',
                                          'Surface Water - Estuary',
                                          'Surface Water - Lake',
                                          'Surface Water - Pond',
                                          'Surface Water - River/Stream',
                                          'Surface Water - Sea or Ocean',
                                          'Surface Water - unspecific',
                                          'Tap Water',
                                          'Well Water (untreated)'])]
set(water_PhACs['Unit standard'])
# just remove 'mg/kg dry-weight' and assume others are roughly to be the same as per volume
water_PhACs = water_PhACs[water_PhACs['Unit standard'].isin(['mg/kg','mg/kg moist-weight','µg/L'])]
water_PhACs.loc[water_PhACs['Unit standard'].isin(['mg/kg','mg/kg moist-weight']), 'MEC standardized'] *= 1000

#%% PhACs concentration visualization

print('\n' + str([len(i) for i in [WWRS_PhACs, soil_PhACs, sediment_PhACs, air_PhACs, water_PhACs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim([10**-7, 10**3])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{C}$ [ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil*','sediment*','air','water*'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([WWRS_PhACs['MEC standardized'], soil_PhACs['MEC standardized'], sediment_PhACs['MEC standardized'], air_PhACs['air'], water_PhACs['MEC standardized']],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=do, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lo, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PhACs_concentration.pdf', transparent=True, bbox_inches='tight')

#%% PhACs concentration visualization (unique values)

print('\n' + str([len(i['MEC standardized'].drop_duplicates()) for i in [WWRS_PhACs, soil_PhACs, sediment_PhACs]]))
print(str(len(air_PhACs['air'].drop_duplicates())))
print(str(len(water_PhACs['MEC standardized'].drop_duplicates())) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim([10**-7, 10**3])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PhACs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil*','sediment*','air','water*'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([WWRS_PhACs['MEC standardized'].drop_duplicates(), soil_PhACs['MEC standardized'].drop_duplicates(), sediment_PhACs['MEC standardized'].drop_duplicates(), air_PhACs['air'].drop_duplicates(), water_PhACs['MEC standardized'].drop_duplicates()],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=do, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lo, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PhACs_concentration_unique.pdf', transparent=True, bbox_inches='tight')

#%% PhACs capture

unused_PhACs = shape.Uniform(0.15, 0.98)
take_back = shape.Triangle(0.136, 0.169, 0.203)
toilet_trash_ratio = shape.Triangle(0.319, 0.399, 0.479)
toilet_to_WRRF = shape.Uniform(0.398, 0.596)

human_use = shape.Uniform(0.784, 1)
human_excretion = shape.Triangle(0.02, 0.39, 0.862)
animal_excretion = shape.Triangle(0.02, 0.39, 0.862)
removal_rate = shape.Triangle(0.15, 0.784, 1)
sorption_biotransformation_ratio = shape.Triangle(0.168, 2.91, 7)

joint = cp.distributions.J(unused_PhACs, take_back, toilet_trash_ratio,
                           toilet_to_WRRF, human_use, human_excretion,
                           animal_excretion, removal_rate, sorption_biotransformation_ratio)

sample = joint.sample(100000)

PhACs_WWRS_MC = pd.DataFrame()

unused_PhACs_to_WRRF = sample[0]*(1 - sample[1])*sample[2]/(sample[2] + 1)*sample[3]
human_use_to_WRRF = (1 - sample[0])*sample[4]*sample[5]*sample[3]
sorption_ratio = sample[7]*sample[8]/(sample[8] + 1)

unused_PhACs_to_environment = sample[0]*(1 - sample[1])*sample[2]/(sample[2] + 1)*(1 - sample[3])
human_use_to_environment = (1 - sample[0])*sample[4]*sample[5]*(1 - sample[3])
animal_use_to_environment = (1 - sample[0])*(1 - sample[4])*sample[6]
WRRF_to_environment_ratio = (1 - sample[7]) + sample[7]*sample[8]/(sample[8] + 1)

PhACs_WWRS = (unused_PhACs_to_WRRF + human_use_to_WRRF)*sorption_ratio

PhACs_total = unused_PhACs_to_environment + human_use_to_environment + animal_use_to_environment +\
              unused_PhACs_to_WRRF*WRRF_to_environment_ratio + human_use_to_WRRF*WRRF_to_environment_ratio

PhACs_WWRS_MC = PhACs_WWRS/PhACs_total

fig, ax = plt.subplots(figsize=(2.5, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 1])
ax.set_ylim([0, 40])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 50, 10), fontname='Arial')

ax.set_ylabel('$\mathbf{Capture}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 50, 10), fontname='Arial')

ax.bar(0.5,
       np.quantile(PhACs_WWRS_MC, 0.5)*100,
       yerr=np.array([[(np.quantile(PhACs_WWRS_MC, 0.5) - np.quantile(PhACs_WWRS_MC, 0.05))*100], [(np.quantile(PhACs_WWRS_MC, 0.95) - np.quantile(PhACs_WWRS_MC, 0.5))*100]]),
       error_kw=dict(capsize=15, lw=3, capthick=3),
       width=0.7,
       color=do,
       edgecolor='k',
       linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PhACs_capture.pdf', transparent=True, bbox_inches='tight')

#%% PhACs removal

PhACs_removal = pd.read_excel(folder + 'analyses/EC_removal_data.xlsx','PhACs_summary')

conventional_PhACs_removal = PhACs_removal[PhACs_removal['process'] == 'C']
wet_thermochemical_PhACs_removal = PhACs_removal[PhACs_removal['process'] == 'wet_T']
dry_thermochemical_PhACs_removal = PhACs_removal[PhACs_removal['process'] == 'dry_T']

fig = plt.figure(figsize=(20.16, 5))
gs = GridSpec(1, 2, width_ratios=[4, 1.25], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 700])
ax.set_ylim([-200, 100])
ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

ax.yaxis.set_major_formatter(mtick.PercentFormatter())

plt.xticks(np.arange(0, 800, 100), fontname='Arial')
plt.yticks(np.arange(-200, 150, 50), fontname='Arial')

# ax.set_xlabel('$\mathbf{Temperature}$ [°C]',
#               fontname='Arial',
#               fontsize=45,
#               labelpad=13)

ax.set_ylabel('$\mathbf{Removal}$',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 800, 100), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-200, 150, 50), fontname='Arial')

ax.scatter(conventional_PhACs_removal['temperature'], conventional_PhACs_removal['reduction']*100, marker='o', color=lo, s=1000, edgecolor='k', linewidth=3)
ax.scatter(wet_thermochemical_PhACs_removal['temperature'], wet_thermochemical_PhACs_removal['reduction']*100, marker='^', color=o, s=1000, edgecolor='k', linewidth=3)
ax.scatter(dry_thermochemical_PhACs_removal['temperature'], dry_thermochemical_PhACs_removal['reduction']*100, marker='s', color=do, s=1000, edgecolor='k', linewidth=3)

ax_box = fig.add_subplot(gs[0, 1], sharey=ax)

data = [conventional_PhACs_removal['reduction']*100, wet_thermochemical_PhACs_removal['reduction']*100, dry_thermochemical_PhACs_removal['reduction']*100]
bp = ax_box.boxplot(data, vert=True, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=lo, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=o, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=do, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

ax_box.axis('off')

plt.savefig('/Users/jiananfeng/Desktop/PhACs_removal.pdf', transparent=True, bbox_inches='tight')

#%% ARGs concentration

# manually collected data
ARGs_0 = pd.read_excel(folder + 'analyses/EC_data.xlsx','ARGs_summary')

ARGs_1 = pd.read_excel(folder + 'analyses/ARGs_data/ARG Abundance.xlsx')
# confirmed all types are ARG types but not mobile gene element types
set(ARGs_1['Type'])
ARGs_1['Absolute Abundance (copies per g/L sample)'] = ARGs_1['Absolute Abundance (copies per g/L sample)'].apply(pd.to_numeric, errors='coerce')
ARGs_1 = ARGs_1[~ARGs_1['Absolute Abundance (copies per g/L sample)'].isna()]
set(ARGs_1['Habitat'])

ARGs_2 = pd.read_excel(folder + 'analyses/ARGs_data/ARGs_NORMAN.xlsx')
set(ARGs_2['Gene name'])
# just remove all intI genes, assume all others are ARGs
ARGs_2 = ARGs_2[~ARGs_2['Gene name'].isin(['intI1','intI1_1','intI1_2','intI1_3','intI1_4',
                                           'intI2','intI2_2','intI3_1','intI3_2','intl1'])]
ARGs_2['gene copy number/mL of sample]'] = ARGs_2['gene copy number/mL of sample]'].apply(pd.to_numeric, errors='coerce')
ARGs_2 = ARGs_2[~ARGs_2['gene copy number/mL of sample]'].isna()]
set(ARGs_2['Sample matrix'])
set(ARGs_2[ARGs_2['Sample matrix'] == 'Other']['Other'])

# WWRS
WWRS_ARGs_1 = ARGs_0['WWRS']
WWRS_ARGs_2 = pd.read_excel(folder + 'analyses/ARGs_data/ARGs_WWRS_Harrison.xlsx')
set(WWRS_ARGs_2['Gene'])

WWRS_ARGs_2 = WWRS_ARGs_2[~WWRS_ARGs_2['Gene'].isin(['16S rRNA','IncQ','intI1'])]
set(WWRS_ARGs_2['Units'])
WWRS_ARGs_2 = WWRS_ARGs_2[~WWRS_ARGs_2['Units'].isna()]
WWRS_ARGs_2 = WWRS_ARGs_2[~WWRS_ARGs_2['Units'].isin(['Log 10 (gene copies/g sample)',
                                                      'Log Concentration (copies/g)',
                                                      'Total ARGs (log10(copies/mL)',
                                                      'copies/g soil',
                                                      'gcn/mL sludge',
                                                      'ln (ARG copies)',
                                                      'log (copies/g)',
                                                      'log copies/g',
                                                      'log10 (copies)/mL',
                                                      'log10 (copies/g sample)',
                                                      'log10 (gene copies/g wet sludge)',
                                                      'ng/g',
                                                      'qnrA copies/g volatile suspended solids'])]
set(WWRS_ARGs_2['Units'])
WWRS_ARGs_2.loc[WWRS_ARGs_2['Units'].isin(['Log10 gene copies/g DS','log (ARG copies/g dry basis)','log 10 (gene copies/g DS)','log10 (gene copies/g DS)']),'Conc Before post-treatment'] =\
    WWRS_ARGs_2.loc[WWRS_ARGs_2['Units'].isin(['Log10 gene copies/g DS','log (ARG copies/g dry basis)','log 10 (gene copies/g DS)','log10 (gene copies/g DS)']),'Conc Before post-treatment'].apply(lambda x: 10**x)
WWRS_ARGs_2.loc[WWRS_ARGs_2['Units'].isin(['Log10 gene copies/g DS','log (ARG copies/g dry basis)','log 10 (gene copies/g DS)','log10 (gene copies/g DS)']),'Conc After post-treatment'] =\
    WWRS_ARGs_2.loc[WWRS_ARGs_2['Units'].isin(['Log10 gene copies/g DS','log (ARG copies/g dry basis)','log 10 (gene copies/g DS)','log10 (gene copies/g DS)']),'Conc After post-treatment'].apply(lambda x: 10**x)

WWRS_ARGs = list(WWRS_ARGs_1.dropna()) + list(WWRS_ARGs_2['Conc Before post-treatment']) + list(WWRS_ARGs_2['Conc After post-treatment'])

# soil
soil_ARGs = list(ARGs_0['soil'].dropna())

# sediment
sediment_ARGs = list(ARGs_0['sediment'].dropna())

# air
air_ARGs = list(ARGs_0['air'].dropna())

# water
water_ARGs_0 = list(ARGs_0['ocean'].dropna())

water_ARGs_1 = ARGs_1[ARGs_1['Habitat'] == 'Water']
water_ARGs_1['Absolute Abundance (copies per g/L sample)'] *= (1/1000)

water_ARGs_2 = ARGs_2[ARGs_2['Sample matrix'].isin(['Groundwater','Other','Surface water - Lake water',
                                                    'Surface water - Other','Surface water - Reservoirs',
                                                    'Surface water - River water'])]

water_ARGs = water_ARGs_0 + list(water_ARGs_1['Absolute Abundance (copies per g/L sample)']) + list(water_ARGs_2['gene copy number/mL of sample]'])

#%% ARGs concentration visualization

print('\n' + str([len(i) for i in [WWRS_ARGs, soil_ARGs, sediment_ARGs, air_ARGs, water_ARGs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim([10**-3, 10**17])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-3, 10**1, 10**5, 10**9, 10**13, 10**17], fontname='Arial')

ax.set_ylabel('$\mathbf{C}$ [copies·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-3, 10**1, 10**5, 10**9, 10**13, 10**17], fontname='Arial')

bp = ax.boxplot([WWRS_ARGs, soil_ARGs, sediment_ARGs, air_ARGs, water_ARGs],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=dp, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lp, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/ARGs_concentration.pdf', transparent=True, bbox_inches='tight')

#%% ARGs concentration visualization (unique values)

print('\n' + str([len(set(i)) for i in [WWRS_ARGs, soil_ARGs, sediment_ARGs, air_ARGs, water_ARGs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_ylim([10**-3, 10**17])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-3, 10**1, 10**5, 10**9, 10**13, 10**17], fontname='Arial')

ax.set_ylabel('$\mathbf{ARGs}$\n[copies·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-3, 10**1, 10**5, 10**9, 10**13, 10**17], fontname='Arial')

bp = ax.boxplot([list(set(WWRS_ARGs)), list(set(soil_ARGs)), list(set(sediment_ARGs)), list(set(air_ARGs)), list(set(water_ARGs))],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=dp, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lp, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/ARGs_concentration_unique.pdf', transparent=True, bbox_inches='tight')

#%% ARGs capture

human_use = shape.Uniform(0.224, 0.336)
human_excretion = shape.Triangle(0.02, 0.39, 0.862)
excretion_to_WRRF = shape.Uniform(0.398, 0.596)
animal_excretion = shape.Triangle(0.02, 0.39, 0.862)
WWRS_capture_rate = shape.Uniform(0.90, 0.95)

joint = cp.distributions.J(human_use, human_excretion, excretion_to_WRRF,
                           animal_excretion, WWRS_capture_rate)

sample = joint.sample(100000)

ARGs_WWRS_MC = pd.DataFrame()

ARGs_WRRS = sample[0]*sample[1]*sample[2]*sample[4]

# assume ARGs can not be degraded during wastewater treatment processes
ARGs_human_to_environment = sample[0]*sample[1]
ARGs_animal_to_environment = (1 - sample[0])*sample[3]

ARGs_WWRS_MC = ARGs_WRRS/(ARGs_human_to_environment + ARGs_animal_to_environment)

fig, ax = plt.subplots(figsize=(2.5, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 1])
ax.set_ylim([0, 40])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 50, 10), fontname='Arial')

ax.set_ylabel('$\mathbf{Capture}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 50, 10), fontname='Arial')

ax.bar(0.5,
       np.quantile(ARGs_WWRS_MC, 0.5)*100,
       yerr=np.array([[(np.quantile(ARGs_WWRS_MC, 0.5) - np.quantile(ARGs_WWRS_MC, 0.05))*100], [(np.quantile(ARGs_WWRS_MC, 0.95) - np.quantile(ARGs_WWRS_MC, 0.5))*100]]),
       error_kw=dict(capsize=15, lw=3, capthick=3),
       width=0.7,
       color=dp,
       edgecolor='k',
       linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/ARGs_capture.pdf', transparent=True, bbox_inches='tight')

#%% ARGs removal

ARGs_removal = pd.read_excel(folder + 'analyses/EC_removal_data.xlsx','ARGs_summary')

conventional_ARGs_removal = ARGs_removal[ARGs_removal['process'] == 'C']
wet_thermochemical_ARGs_removal = ARGs_removal[ARGs_removal['process'] == 'wet_T']
dry_thermochemical_ARGs_removal = ARGs_removal[ARGs_removal['process'] == 'dry_T']

fig = plt.figure(figsize=(20.16, 5))
gs = GridSpec(1, 2, width_ratios=[4, 1.25], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 700])
ax.set_ylim([-2, 5])
ax.tick_params(direction='inout', length=20, width=3, bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 800, 100), fontname='Arial')
plt.yticks(np.arange(-2, 6, 1), fontname='Arial')

ax.set_xlabel('$\mathbf{Temperature}$ [°C]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax.set_ylabel('$\mathbf{Log\ removal}$',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.tick_params(direction='in', length=10, width=3, bottom=False, top=True, left=False, right=False, labelcolor='none')

plt.xticks(np.arange(0, 800, 100), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(-2, 6, 1), fontname='Arial')

ax.scatter(conventional_ARGs_removal['temperature'], conventional_ARGs_removal['reduction'], marker='o', color=lp, s=1000, edgecolor='k', linewidth=3)
ax.scatter(wet_thermochemical_ARGs_removal['temperature'], wet_thermochemical_ARGs_removal['reduction'], marker='^', color=p, s=1000, edgecolor='k', linewidth=3)
ax.scatter(dry_thermochemical_ARGs_removal['temperature'], dry_thermochemical_ARGs_removal['reduction'], marker='s', color=dp, s=1000, edgecolor='k', linewidth=3)

ax_box = fig.add_subplot(gs[0, 1], sharey=ax)

data = [conventional_ARGs_removal['reduction'], wet_thermochemical_ARGs_removal['reduction'], dry_thermochemical_ARGs_removal['reduction']]
bp = ax_box.boxplot(data, vert=True, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=lp, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=p, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=dp, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

ax_box.axis('off')

plt.savefig('/Users/jiananfeng/Desktop/ARGs_removal.pdf', transparent=True, bbox_inches='tight')

#%% PFAS concentration data processing 1 - EU - data cleaning to reduce the file size

# EU_PFAS_raw = pd.read_csv(folder + 'analyses/PFAS_data/EU_PFAS_raw.csv')
# EU_PFAS = EU_PFAS_raw[['pfas_values','unit','matrix']]
# EU_PFAS = EU_PFAS.dropna(subset='matrix')
# EU_PFAS = EU_PFAS[EU_PFAS['matrix'].isin(['Drinking water','Groundwater','Rainwater','Sea water','Sediment','Soil','Surface water'])]
# EU_PFAS['pfas_values'] = EU_PFAS['pfas_values'].apply(json.loads)

# def extract_totals(pfas_list):
#     PFOS_conc, PFOS_unit, PFOA_conc, PFOA_unit = None, None, None, None
    
#     for item in pfas_list:
#         if 'isomer' in item:
#             continue
        
#         if item['substance'] == 'PFOA':
#             if 'value' in item.keys():
#                 val = float(item.get('value'))
#                 PFOA_conc, PFOA_unit = val, item['unit']
#         elif item['substance'] == 'PFOS':
#             if 'value' in item.keys():
#                 val = float(item.get('value'))
#                 PFOS_conc, PFOS_unit = val, item['unit']

#     return pd.Series({'PFOA_conc': PFOA_conc,
#                       'PFOA_unit': PFOA_unit,
#                       'PFOS_conc': PFOS_conc,
#                       'PFOS_unit': PFOS_unit})

# EU_PFAS_results = EU_PFAS['pfas_values'].apply(extract_totals)

# EU_PFAS = pd.concat([EU_PFAS, EU_PFAS_results], axis=1)

# EU_PFOA = EU_PFAS[['PFOA_conc','PFOA_unit','matrix']]
# EU_PFOA = EU_PFOA.dropna(subset='PFOA_conc')

# EU_PFOS = EU_PFAS[['PFOS_conc','PFOS_unit','matrix']]
# EU_PFOS = EU_PFOS.dropna(subset='PFOS_conc')

# EU_PFOA.to_excel(folder + 'analyses/PFAS_data/EU_PFOA.xlsx')
# EU_PFOS.to_excel(folder + 'analyses/PFAS_data/EU_PFOS.xlsx')

#%% PFOA concentration

# PFOA results in ng/g
WWRS_PFOA = []
soil_PFOA = []
sediment_PFOA = []
air_PFOA = []
water_PFOA = []

# =============================================================================
# CN_1_PFOA
# =============================================================================

CN_1_PFOA = pd.read_excel(folder + 'analyses/PFAS_data/CN_1_PFOA.xlsx')
set(CN_1_PFOA['Pathway'])
CN_1_PFOA = CN_1_PFOA[CN_1_PFOA['Pathway'].isin(['air','others','sediment','soil','water'])]
CN_1_PFOA = CN_1_PFOA[((CN_1_PFOA['Pathway'] == 'others') & (CN_1_PFOA['Detailed_pathway'].str.contains('sludge|Sludge'))) | (CN_1_PFOA['Pathway'] != 'others')]
# remove entries containing 'sludge' but are not WWRS
set(CN_1_PFOA[(CN_1_PFOA['Pathway'] == 'others') & (CN_1_PFOA['Detailed_pathway'].str.contains('sludge'))]['Detailed_pathway'])
CN_1_PFOA = CN_1_PFOA[~CN_1_PFOA['Detailed_pathway'].isin(['drinking water sludge','surface sediment&sludge'])]
CN_1_PFOA = CN_1_PFOA[['Pathway','Unit','Min','Max','Mean','Median']]
CN_1_PFOA.replace('nd', 0, inplace=True)
for column in ['Min','Max','Mean','Median']:
    CN_1_PFOA[column] = CN_1_PFOA[column].apply(pd.to_numeric, errors='coerce')

CN_1_PFOA_WWRS = CN_1_PFOA[CN_1_PFOA['Pathway'] == 'others']
CN_1_PFOA_soil = CN_1_PFOA[CN_1_PFOA['Pathway'] == 'soil']
CN_1_PFOA_sediment = CN_1_PFOA[CN_1_PFOA['Pathway'] == 'sediment']
CN_1_PFOA_air = CN_1_PFOA[CN_1_PFOA['Pathway'] == 'air']
CN_1_PFOA_water = CN_1_PFOA[CN_1_PFOA['Pathway'] == 'water']

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['Min']) or not pd.isna(row['Max']):
        if not pd.isna(row['Min']):
            values.append(row['Min'])
        if not pd.isna(row['Max']):
            values.append(row['Max'])
    else:
        if not pd.isna(row['Mean']):
            values.append(row['Mean'])
    
    if not pd.isna(row['Median']):
        values.append(row['Median'])
    return values

# WWRS
set(CN_1_PFOA_WWRS['Unit'])
CN_1_PFOA_WWRS = CN_1_PFOA_WWRS[CN_1_PFOA_WWRS['Unit'].str.contains('dw')]
CN_1_PFOA_WWRS['valid_data'] = CN_1_PFOA_WWRS.apply(add_data, axis=1)

WWRS_PFOA += [x for sublist in CN_1_PFOA_WWRS['valid_data'] for x in sublist]

# soil
set(CN_1_PFOA_soil['Unit'])
CN_1_PFOA_soil = CN_1_PFOA_soil[CN_1_PFOA_soil['Unit'].str.contains('dw')]
CN_1_PFOA_soil.loc[CN_1_PFOA_soil['Unit'] == 'pg/g dw', ['Min','Max','Mean','Median']] *= (1/1000)
CN_1_PFOA_soil['valid_data'] = CN_1_PFOA_soil.apply(add_data, axis=1)

soil_PFOA += [x for sublist in CN_1_PFOA_soil['valid_data'] for x in sublist]

# sediment
set(CN_1_PFOA_sediment['Unit'])
CN_1_PFOA_sediment = CN_1_PFOA_sediment[CN_1_PFOA_sediment['Unit'].str.contains('dw')]
CN_1_PFOA_sediment.loc[CN_1_PFOA_sediment['Unit'] == 'pg/g dw', ['Min','Max','Mean','Median']] *= (1/1000)
CN_1_PFOA_sediment['valid_data'] = CN_1_PFOA_sediment.apply(add_data, axis=1)

sediment_PFOA += [x for sublist in CN_1_PFOA_sediment['valid_data'] for x in sublist]

# air
set(CN_1_PFOA_air['Unit'])
# assume the density of air is 1.225 kg/m3
CN_1_PFOA_air.loc[CN_1_PFOA_air['Unit'] == 'ng/L', ['Min','Max','Mean','Median']] *= (1000/1.225/1000)
CN_1_PFOA_air.loc[CN_1_PFOA_air['Unit'] == 'ng/m3', ['Min','Max','Mean','Median']] *= (1/1.225/1000)
CN_1_PFOA_air.loc[CN_1_PFOA_air['Unit'] == 'pg/L', ['Min','Max','Mean','Median']] *= (1000/1.225/1000/1000)
CN_1_PFOA_air.loc[CN_1_PFOA_air['Unit'] == 'pg/m3', ['Min','Max','Mean','Median']] *= (1/1.225/1000/1000)
CN_1_PFOA_air['valid_data'] = CN_1_PFOA_air.apply(add_data, axis=1)

air_PFOA += [x for sublist in CN_1_PFOA_air['valid_data'] for x in sublist]

# water
set(CN_1_PFOA_water['Unit'])
# assume the density of water is 1000 kg/m3
CN_1_PFOA_water.loc[CN_1_PFOA_water['Unit'].isin(['ng/L','ng/l','pg/mL']), ['Min','Max','Mean','Median']] *= (1/1000)
CN_1_PFOA_water.loc[CN_1_PFOA_water['Unit'] == 'pg/L', ['Min','Max','Mean','Median']] *= (1/1000/1000)
# PFOA MW: 414.07 g/mpl
CN_1_PFOA_water.loc[CN_1_PFOA_water['Unit'] == 'pmol/L', ['Min','Max','Mean','Median']] *= (414.07/1000/1000)
CN_1_PFOA_water['valid_data'] = CN_1_PFOA_water.apply(add_data, axis=1)

water_PFOA += [x for sublist in CN_1_PFOA_water['valid_data'] for x in sublist]

# =============================================================================
# CN_2_PFAS
# =============================================================================

CN_2_PFAS = pd.read_excel(folder + 'analyses/PFAS_data/CN_2_PFAS.xlsx')
CN_2_PFOA = CN_2_PFAS[CN_2_PFAS['ECs_type'] == 'PFOA']
set(CN_2_PFOA['Ecs_Pathway'])
CN_2_PFOA = CN_2_PFOA[CN_2_PFOA['Ecs_Pathway'].isin(['Sediment','Soil','air','sediment','sludge','soil','water'])]
CN_2_PFOA = CN_2_PFOA[['Ecs_Pathway','Ecs_unit','ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']]
for column in ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']:
    CN_2_PFOA[column] = CN_2_PFOA[column].apply(pd.to_numeric, errors='coerce')

CN_2_PFOA_WWRS = CN_2_PFOA[CN_2_PFOA['Ecs_Pathway'] == 'sludge']
CN_2_PFOA_soil = CN_2_PFOA[CN_2_PFOA['Ecs_Pathway'].isin(['Soil','soil'])]
CN_2_PFOA_sediment = CN_2_PFOA[CN_2_PFOA['Ecs_Pathway'].isin(['Sediment','sediment'])]
CN_2_PFOA_air = CN_2_PFOA[CN_2_PFOA['Ecs_Pathway'] == 'air']
CN_2_PFOA_water = CN_2_PFOA[CN_2_PFOA['Ecs_Pathway'] == 'water']

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['ECs_conc_min']) or not pd.isna(row['ECs_conc_max']):
        if not pd.isna(row['ECs_conc_min']):
            values.append(row['ECs_conc_min'])
        if not pd.isna(row['ECs_conc_max']):
            values.append(row['ECs_conc_max'])
    else:
        if not pd.isna(row['ECs_conc_mean']):
            values.append(row['ECs_conc_mean'])
    
    if not pd.isna(row['Ecs_conc_median']):
        values.append(row['Ecs_conc_median'])
    return values

# WWRS
set(CN_2_PFOA_WWRS['Ecs_unit'])

# soil
set(CN_2_PFOA_soil['Ecs_unit'])
CN_2_PFOA_soil = CN_2_PFOA_soil[CN_2_PFOA_soil['Ecs_unit'].str.contains('dw')]
CN_2_PFOA_soil.loc[CN_2_PFOA_soil['Ecs_unit'] == 'pg/g dw', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000)
CN_2_PFOA_soil['valid_data'] = CN_2_PFOA_soil.apply(add_data, axis=1)

soil_PFOA += [x for sublist in CN_2_PFOA_soil['valid_data'] for x in sublist]

# sediment
set(CN_2_PFOA_sediment['Ecs_unit'])
CN_2_PFOA_sediment = CN_2_PFOA_sediment[CN_2_PFOA_sediment['Ecs_unit'].str.contains('dw')]
CN_2_PFOA_sediment.loc[CN_2_PFOA_sediment['Ecs_unit'].isin(['ng/kg  dw','pg/g dw']), ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000)
CN_2_PFOA_sediment['valid_data'] = CN_2_PFOA_sediment.apply(add_data, axis=1)

sediment_PFOA += [x for sublist in CN_2_PFOA_sediment['valid_data'] for x in sublist]

# air
set(CN_2_PFOA_air['Ecs_unit'])
# assume the density of air is 1.225 kg/m3
CN_2_PFOA_air.loc[CN_2_PFOA_air['Ecs_unit'] == 'ng/L', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1000/1.225/1000)
CN_2_PFOA_air.loc[CN_2_PFOA_air['Ecs_unit'] == 'ng/m3', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1.225/1000)
CN_2_PFOA_air.loc[CN_2_PFOA_air['Ecs_unit'] == 'pg/m3', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1.225/1000/1000)
CN_2_PFOA_air['valid_data'] = CN_2_PFOA_air.apply(add_data, axis=1)

air_PFOA += [x for sublist in CN_2_PFOA_air['valid_data'] for x in sublist]

# water
set(CN_2_PFOA_water['Ecs_unit'])
# assume the density of water is 1000 kg/m3
CN_2_PFOA_water.loc[CN_2_PFOA_water['Ecs_unit'] == 'ng/L', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000)
CN_2_PFOA_water.loc[CN_2_PFOA_water['Ecs_unit'] == 'pg/L', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000/1000)
CN_2_PFOA_water['valid_data'] = CN_2_PFOA_water.apply(add_data, axis=1)

water_PFOA += [x for sublist in CN_2_PFOA_water['valid_data'] for x in sublist]

# =============================================================================
# US_PFAS
# =============================================================================

US_PFAS = pd.read_excel(folder + 'analyses/PFAS_data/US_PFAS.xlsx')
# PFOA CAS number: 335-67-1
US_PFOA = US_PFAS[US_PFAS['CAS Number'] == '335-67-1']
set(US_PFOA['Environmental Media Name'])
US_PFOA = US_PFOA[US_PFOA['Environmental Media Name'].isin(['Air','Other','Sediment','Soil','Water'])]
set(US_PFOA['Activity Media Subdivision Name'])
US_PFOA = US_PFOA[US_PFOA['Environmental Media Name'].isin(['Air','Sediment','Soil','Water'])]
US_PFOA = US_PFOA[['Environmental Media Name','Result Measure Value','Result Unit of Measure','Activity Comment']]
US_PFOA['Activity Comment'] = US_PFOA['Activity Comment'].fillna(' ')
US_PFOA.replace('-', 0, inplace=True)
US_PFOA['Result Measure Value'] = US_PFOA['Result Measure Value'].apply(pd.to_numeric, errors='coerce')

US_PFOA_soil = US_PFOA[US_PFOA['Environmental Media Name'] == 'Soil']
US_PFOA_sediment = US_PFOA[US_PFOA['Environmental Media Name'] == 'Sediment']
US_PFOA_air = US_PFOA[US_PFOA['Environmental Media Name'] == 'Air']
US_PFOA_water = US_PFOA[US_PFOA['Environmental Media Name'] == 'Water']

# soil
set(US_PFOA_soil['Activity Comment'])
set(US_PFOA_soil['Result Unit of Measure'])

# sediment
set(US_PFOA_sediment['Activity Comment'])
set(US_PFOA_sediment['Result Unit of Measure'])

# air
set(US_PFOA_air['Result Unit of Measure'])
US_PFOA_air = US_PFOA_air[US_PFOA_air['Result Unit of Measure'] == 'pg/L']
# assume the density of air is 1.225 kg/m3
US_PFOA_air.loc[US_PFOA_air['Result Unit of Measure'] == 'pg/L', 'Result Measure Value'] *= (1000/1.225/1000/1000)

air_PFOA += list(US_PFOA_air['Result Measure Value'])

# water
set(US_PFOA_water['Result Unit of Measure'])
US_PFOA_water = US_PFOA_water[US_PFOA_water['Result Unit of Measure'].isin(['ng/L','ug/L','ug/kg'])]
# assume the density of water is 1000 kg/m3
US_PFOA_water.loc[US_PFOA_water['Result Unit of Measure'] == 'ng/L', 'Result Measure Value'] *= (1/1000)

water_PFOA += list(US_PFOA_water['Result Measure Value'])

# =============================================================================
# EU_PFOA
# =============================================================================

EU_PFOA = pd.read_excel(folder + 'analyses/PFAS_data/EU_PFOA.xlsx')
set(EU_PFOA['matrix'])
EU_PFOA['PFOA_conc'] = EU_PFOA['PFOA_conc'].apply(pd.to_numeric, errors='coerce')

EU_PFOA_soil = EU_PFOA[EU_PFOA['matrix'] == 'Soil']
EU_PFOA_sediment = EU_PFOA[EU_PFOA['matrix'] == 'Sediment']
EU_PFOA_water = EU_PFOA[EU_PFOA['matrix'].isin(['Drinking water','Groundwater','Rainwater','Sea water','Surface water'])]

# soil
set(EU_PFOA_soil['PFOA_unit'])
EU_PFOA_soil = EU_PFOA_soil[EU_PFOA_soil['PFOA_unit'].str.contains('dw')]
EU_PFOA_soil['PFOA_conc'] *= (1/1000)

soil_PFOA += list(EU_PFOA_soil['PFOA_conc'])

# sediment
set(EU_PFOA_sediment['PFOA_unit'])
EU_PFOA_sediment = EU_PFOA_sediment[EU_PFOA_sediment['PFOA_unit'].str.contains('dw')]
EU_PFOA_sediment['PFOA_conc'] *= (1/1000)

sediment_PFOA += list(EU_PFOA_sediment['PFOA_conc'])

# water
set(EU_PFOA_water['PFOA_unit'])
# assume the density of water is 1000 kg/m3
EU_PFOA_water['PFOA_conc'] *= (1/1000)

water_PFOA += list(EU_PFOA_water['PFOA_conc'])

# =============================================================================
# additional WWRS PFOA
# =============================================================================

additional_WWRS_PFOA = pd.read_excel(folder + 'analyses/EC_data.xlsx','PFAS_summary')

WWRS_PFOA += list(additional_WWRS_PFOA['PFOA'].dropna())

#%% PFOA concentrations visualization

print('\n' + str([len(i) for i in [WWRS_PFOA, soil_PFOA, sediment_PFOA, air_PFOA, water_PFOA]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-7, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PFOA}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([WWRS_PFOA, soil_PFOA, sediment_PFOA, air_PFOA, water_PFOA],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PFOA.pdf', transparent=True, bbox_inches='tight')

#%% PFOA concentrations visualization (unique values)

print('\n' + str([len(set(i)) for i in [WWRS_PFOA, soil_PFOA, sediment_PFOA, air_PFOA, water_PFOA]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-7, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PFOA}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([list(set(WWRS_PFOA)), list(set(soil_PFOA)), list(set(sediment_PFOA)), list(set(air_PFOA)), list(set(water_PFOA))],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PFOA_unique.pdf', transparent=True, bbox_inches='tight')

#%% PFOS concentration

# PFOS results in ng/g
WWRS_PFOS = []
soil_PFOS = []
sediment_PFOS = []
air_PFOS = []
water_PFOS = []

# =============================================================================
# CN_1_PFOS
# =============================================================================

CN_1_PFOS = pd.read_excel(folder + 'analyses/PFAS_data/CN_1_PFOS.xlsx')
set(CN_1_PFOS['Pathway'])
CN_1_PFOS = CN_1_PFOS[CN_1_PFOS['Pathway'].isin(['air','others','sediment','soil','water'])]
CN_1_PFOS = CN_1_PFOS[((CN_1_PFOS['Pathway'] == 'others') & (CN_1_PFOS['Detailed_pathway'].str.contains('sludge|Sludge'))) | (CN_1_PFOS['Pathway'] != 'others')]
# remove entries containing 'sludge' but are not WWRS
set(CN_1_PFOS[(CN_1_PFOS['Pathway'] == 'others') & (CN_1_PFOS['Detailed_pathway'].str.contains('sludge'))]['Detailed_pathway'])
CN_1_PFOS = CN_1_PFOS[~CN_1_PFOS['Detailed_pathway'].isin(['drinking water sludge','surface sediment&sludge'])]
CN_1_PFOS = CN_1_PFOS[['Pathway','Unit','Min','Max','Mean','Median']]
CN_1_PFOS.replace('nd', 0, inplace=True)
for column in ['Min','Max','Mean','Median']:
    CN_1_PFOS[column] = CN_1_PFOS[column].apply(pd.to_numeric, errors='coerce')

CN_1_PFOS_WWRS = CN_1_PFOS[CN_1_PFOS['Pathway'] == 'others']
CN_1_PFOS_soil = CN_1_PFOS[CN_1_PFOS['Pathway'] == 'soil']
CN_1_PFOS_sediment = CN_1_PFOS[CN_1_PFOS['Pathway'] == 'sediment']
CN_1_PFOS_air = CN_1_PFOS[CN_1_PFOS['Pathway'] == 'air']
CN_1_PFOS_water = CN_1_PFOS[CN_1_PFOS['Pathway'] == 'water']

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['Min']) or not pd.isna(row['Max']):
        if not pd.isna(row['Min']):
            values.append(row['Min'])
        if not pd.isna(row['Max']):
            values.append(row['Max'])
    else:
        if not pd.isna(row['Mean']):
            values.append(row['Mean'])
    
    if not pd.isna(row['Median']):
        values.append(row['Median'])
    return values

# WWRS
set(CN_1_PFOS_WWRS['Unit'])
CN_1_PFOS_WWRS = CN_1_PFOS_WWRS[CN_1_PFOS_WWRS['Unit'].str.contains('dw')]
CN_1_PFOS_WWRS['valid_data'] = CN_1_PFOS_WWRS.apply(add_data, axis=1)

WWRS_PFOS += [x for sublist in CN_1_PFOS_WWRS['valid_data'] for x in sublist]

# soil
set(CN_1_PFOS_soil['Unit'])
CN_1_PFOS_soil = CN_1_PFOS_soil[CN_1_PFOS_soil['Unit'].str.contains('dw')]
CN_1_PFOS_soil.loc[CN_1_PFOS_soil['Unit'] == 'pg/g dw', ['Min','Max','Mean','Median']] *= (1/1000)
CN_1_PFOS_soil['valid_data'] = CN_1_PFOS_soil.apply(add_data, axis=1)

soil_PFOS += [x for sublist in CN_1_PFOS_soil['valid_data'] for x in sublist]

# sediment
set(CN_1_PFOS_sediment['Unit'])
CN_1_PFOS_sediment = CN_1_PFOS_sediment[CN_1_PFOS_sediment['Unit'].str.contains('dw')]
CN_1_PFOS_sediment.loc[CN_1_PFOS_sediment['Unit'] == 'pg/g dw', ['Min','Max','Mean','Median']] *= (1/1000)
CN_1_PFOS_sediment['valid_data'] = CN_1_PFOS_sediment.apply(add_data, axis=1)

sediment_PFOS += [x for sublist in CN_1_PFOS_sediment['valid_data'] for x in sublist]

# air
set(CN_1_PFOS_air['Unit'])
# assume the density of air is 1.225 kg/m3
CN_1_PFOS_air.loc[CN_1_PFOS_air['Unit'] == 'ng/L', ['Min','Max','Mean','Median']] *= (1000/1.225/1000)
CN_1_PFOS_air.loc[CN_1_PFOS_air['Unit'] == 'ng/m3', ['Min','Max','Mean','Median']] *= (1/1.225/1000)
CN_1_PFOS_air.loc[CN_1_PFOS_air['Unit'] == 'pg/L', ['Min','Max','Mean','Median']] *= (1000/1.225/1000/1000)
CN_1_PFOS_air.loc[CN_1_PFOS_air['Unit'] == 'pg/m3', ['Min','Max','Mean','Median']] *= (1/1.225/1000/1000)
CN_1_PFOS_air['valid_data'] = CN_1_PFOS_air.apply(add_data, axis=1)

air_PFOS += [x for sublist in CN_1_PFOS_air['valid_data'] for x in sublist]

# water
set(CN_1_PFOS_water['Unit'])
# assume the density of water is 1000 kg/m3
CN_1_PFOS_water.loc[CN_1_PFOS_water['Unit'].isin(['ng/L','ng/l','pg/mL']), ['Min','Max','Mean','Median']] *= (1/1000)
CN_1_PFOS_water.loc[CN_1_PFOS_water['Unit'] == 'pg/L', ['Min','Max','Mean','Median']] *= (1/1000/1000)
# PFOS MW: 500.13 g/mpl
CN_1_PFOS_water.loc[CN_1_PFOS_water['Unit'] == 'pmol/L', ['Min','Max','Mean','Median']] *= (500.13/1000/1000)
CN_1_PFOS_water['valid_data'] = CN_1_PFOS_water.apply(add_data, axis=1)

water_PFOS += [x for sublist in CN_1_PFOS_water['valid_data'] for x in sublist]

# =============================================================================
# CN_2_PFAS
# =============================================================================

CN_2_PFAS = pd.read_excel(folder + 'analyses/PFAS_data/CN_2_PFAS.xlsx')
CN_2_PFOS = CN_2_PFAS[CN_2_PFAS['ECs_type'] == 'PFOS']
set(CN_2_PFOS['Ecs_Pathway'])
CN_2_PFOS = CN_2_PFOS[CN_2_PFOS['Ecs_Pathway'].isin(['air','sediment','sediment ','sludge','soil','water'])]
CN_2_PFOS = CN_2_PFOS[['Ecs_Pathway','Ecs_unit','ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']]
for column in ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']:
    CN_2_PFOS[column] = CN_2_PFOS[column].apply(pd.to_numeric, errors='coerce')

CN_2_PFOS_WWRS = CN_2_PFOS[CN_2_PFOS['Ecs_Pathway'] == 'sludge']
CN_2_PFOS_soil = CN_2_PFOS[CN_2_PFOS['Ecs_Pathway'] == 'soil']
CN_2_PFOS_sediment = CN_2_PFOS[CN_2_PFOS['Ecs_Pathway'].isin(['sediment','sediment '])]
CN_2_PFOS_air = CN_2_PFOS[CN_2_PFOS['Ecs_Pathway'] == 'air']
CN_2_PFOS_water = CN_2_PFOS[CN_2_PFOS['Ecs_Pathway'] == 'water']

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['ECs_conc_min']) or not pd.isna(row['ECs_conc_max']):
        if not pd.isna(row['ECs_conc_min']):
            values.append(row['ECs_conc_min'])
        if not pd.isna(row['ECs_conc_max']):
            values.append(row['ECs_conc_max'])
    else:
        if not pd.isna(row['ECs_conc_mean']):
            values.append(row['ECs_conc_mean'])
    
    if not pd.isna(row['Ecs_conc_median']):
        values.append(row['Ecs_conc_median'])
    return values

# WWRS
set(CN_2_PFOS_WWRS['Ecs_unit'])

# soil
set(CN_2_PFOS_soil['Ecs_unit'])
CN_2_PFOS_soil = CN_2_PFOS_soil[CN_2_PFOS_soil['Ecs_unit'].str.contains('dw')]
CN_2_PFOS_soil.loc[CN_2_PFOS_soil['Ecs_unit'] == 'pg/g dw', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000)
CN_2_PFOS_soil['valid_data'] = CN_2_PFOS_soil.apply(add_data, axis=1)

soil_PFOS += [x for sublist in CN_2_PFOS_soil['valid_data'] for x in sublist]

# sediment
set(CN_2_PFOS_sediment['Ecs_unit'])
CN_2_PFOS_sediment = CN_2_PFOS_sediment[CN_2_PFOS_sediment['Ecs_unit'].str.contains('dw')]
CN_2_PFOS_sediment.loc[CN_2_PFOS_sediment['Ecs_unit'] == 'ng/kg  dw', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000)
CN_2_PFOS_sediment['valid_data'] = CN_2_PFOS_sediment.apply(add_data, axis=1)

sediment_PFOS += [x for sublist in CN_2_PFOS_sediment['valid_data'] for x in sublist]

# air
set(CN_2_PFOS_air['Ecs_unit'])
# assume the density of air is 1.225 kg/m3
CN_2_PFOS_air.loc[CN_2_PFOS_air['Ecs_unit'] == 'pg/m3', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1.225/1000/1000)
CN_2_PFOS_air['valid_data'] = CN_2_PFOS_air.apply(add_data, axis=1)

air_PFOS += [x for sublist in CN_2_PFOS_air['valid_data'] for x in sublist]

# water
set(CN_2_PFOS_water['Ecs_unit'])
# assume the density of water is 1000 kg/m3
CN_2_PFOS_water.loc[CN_2_PFOS_water['Ecs_unit'] == 'ng/L', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000)
CN_2_PFOS_water.loc[CN_2_PFOS_water['Ecs_unit'] == 'pg/L', ['ECs_conc_min','ECs_conc_max','ECs_conc_mean','Ecs_conc_median']] *= (1/1000/1000)
CN_2_PFOS_water['valid_data'] = CN_2_PFOS_water.apply(add_data, axis=1)

water_PFOS += [x for sublist in CN_2_PFOS_water['valid_data'] for x in sublist]

# =============================================================================
# US_PFAS
# =============================================================================

US_PFAS = pd.read_excel(folder + 'analyses/PFAS_data/US_PFAS.xlsx')
# PFOS CAS number: 1763-23-1
US_PFOS = US_PFAS[US_PFAS['CAS Number'] == '1763-23-1']
set(US_PFOS['Environmental Media Name'])
US_PFOS = US_PFOS[US_PFOS['Environmental Media Name'].isin(['Air','Sediment','Water'])]
US_PFOS = US_PFOS[['Environmental Media Name','Result Measure Value','Result Unit of Measure','Activity Comment']]
US_PFOS['Activity Comment'] = US_PFOS['Activity Comment'].fillna(' ')
US_PFOS.replace('-', 0, inplace=True)
US_PFOS['Result Measure Value'] = US_PFOS['Result Measure Value'].apply(pd.to_numeric, errors='coerce')

US_PFOS_sediment = US_PFOS[US_PFOS['Environmental Media Name'] == 'Sediment']
US_PFOS_air = US_PFOS[US_PFOS['Environmental Media Name'] == 'Air']
US_PFOS_water = US_PFOS[US_PFOS['Environmental Media Name'] == 'Water']

# sediment
set(US_PFOS_sediment['Activity Comment'])
set(US_PFOS_sediment['Result Unit of Measure'])

# air
set(US_PFOS_air['Result Unit of Measure'])
US_PFOS_air = US_PFOS_air[US_PFOS_air['Result Unit of Measure'] == 'pg/L']
# assume the density of air is 1.225 kg/m3
US_PFOS_air.loc[US_PFOS_air['Result Unit of Measure'] == 'pg/L', 'Result Measure Value'] *= (1000/1.225/1000/1000)

air_PFOS += list(US_PFOS_air['Result Measure Value'])

# water
set(US_PFOS_water['Result Unit of Measure'])
US_PFOS_water = US_PFOS_water[US_PFOS_water['Result Unit of Measure'].isin(['ng/L','ug/L','ug/kg'])]
# assume the density of water is 1000 kg/m3
US_PFOS_water.loc[US_PFOS_water['Result Unit of Measure'] == 'ng/L', 'Result Measure Value'] *= (1/1000)

water_PFOS += list(US_PFOS_water['Result Measure Value'])

# =============================================================================
# EU_PFOS
# =============================================================================

EU_PFOS = pd.read_excel(folder + 'analyses/PFAS_data/EU_PFOS.xlsx')
set(EU_PFOS['matrix'])
EU_PFOS['PFOS_conc'] = EU_PFOS['PFOS_conc'].apply(pd.to_numeric, errors='coerce')

EU_PFOS_soil = EU_PFOS[EU_PFOS['matrix'] == 'Soil']
EU_PFOS_sediment = EU_PFOS[EU_PFOS['matrix'] == 'Sediment']
EU_PFOS_water = EU_PFOS[EU_PFOS['matrix'].isin(['Drinking water','Groundwater','Rainwater','Sea water','Surface water'])]

# soil
set(EU_PFOS_soil['PFOS_unit'])
EU_PFOS_soil = EU_PFOS_soil[EU_PFOS_soil['PFOS_unit'].str.contains('dw')]
EU_PFOS_soil['PFOS_conc'] *= (1/1000)

soil_PFOS += list(EU_PFOS_soil['PFOS_conc'])

# sediment
set(EU_PFOS_sediment['PFOS_unit'])
EU_PFOS_sediment = EU_PFOS_sediment[EU_PFOS_sediment['PFOS_unit'].str.contains('dw')]
EU_PFOS_sediment['PFOS_conc'] *= (1/1000)

sediment_PFOS += list(EU_PFOS_sediment['PFOS_conc'])

# water
set(EU_PFOS_water['PFOS_unit'])
# assume the density of water is 1000 kg/m3
EU_PFOS_water['PFOS_conc'] *= (1/1000)

water_PFOS += list(EU_PFOS_water['PFOS_conc'])

# =============================================================================
# additional WWRS PFOS
# =============================================================================

additional_WWRS_PFOS = pd.read_excel(folder + 'analyses/EC_data.xlsx','PFAS_summary')

WWRS_PFOS += list(additional_WWRS_PFOS['PFOS'].dropna())

#%% PFOS concentrations visualization

print('\n' + str([len(i) for i in [WWRS_PFOS, soil_PFOS, sediment_PFOS, air_PFOS, water_PFOS]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-7, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PFOS}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([WWRS_PFOS, soil_PFOS, sediment_PFOS, air_PFOS, water_PFOS],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PFOS.pdf', transparent=True, bbox_inches='tight')

#%% PFOS concentrations visualization (unique values)

print('\n' + str([len(set(i)) for i in [WWRS_PFOS, soil_PFOS, sediment_PFOS, air_PFOS, water_PFOS]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-7, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PFOS}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([list(set(WWRS_PFOS)), list(set(soil_PFOS)), list(set(sediment_PFOS)), list(set(air_PFOS)), list(set(water_PFOS))],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PFOS_unique.pdf', transparent=True, bbox_inches='tight')

#%% PBDEs concentration data

CN_PBDEs = pd.read_excel(folder + 'analyses/POPs_data/exposure data/Standardized data used for analysis and human health risk assessment/PBDEs_analysis.xlsx')
set(CN_PBDEs['Pathway'])
CN_PBDEs = CN_PBDEs[CN_PBDEs['Pathway'].isin(['Water','air','others','sediment','soil','water'])]
CN_PBDEs = CN_PBDEs[((CN_PBDEs['Pathway'] == 'others') & (CN_PBDEs['Detailed_pathway'].str.contains('sludge|Sludge'))) | (CN_PBDEs['Pathway'] != 'others')]
# remove entries containing 'sludge' but are not WWRS
set(CN_PBDEs[(CN_PBDEs['Pathway'] == 'others') & (CN_PBDEs['Detailed_pathway'].str.contains('sludge'))]['Detailed_pathway'])
CN_PBDEs = CN_PBDEs[['Pathway','Unit','Min','Max','Mean','Median']]
CN_PBDEs.replace('nd', 0, inplace=True)
for column in ['Min','Max','Mean','Median']:
    CN_PBDEs[column] = CN_PBDEs[column].apply(pd.to_numeric, errors='coerce')

CN_PBDEs_WWRS = CN_PBDEs[CN_PBDEs['Pathway'] == 'others']
CN_PBDEs_soil = CN_PBDEs[CN_PBDEs['Pathway'] == 'soil']
CN_PBDEs_sediment = CN_PBDEs[CN_PBDEs['Pathway'] == 'sediment']
CN_PBDEs_air = CN_PBDEs[CN_PBDEs['Pathway'] == 'air']
CN_PBDEs_water = CN_PBDEs[CN_PBDEs['Pathway'].isin(['Water','water'])]

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['Min']) or not pd.isna(row['Max']):
        if not pd.isna(row['Min']):
            values.append(row['Min'])
        if not pd.isna(row['Max']):
            values.append(row['Max'])
    else:
        if not pd.isna(row['Mean']):
            values.append(row['Mean'])
    
    if not pd.isna(row['Median']):
        values.append(row['Median'])
    return values

# WWRS
set(CN_PBDEs_WWRS['Unit'])
CN_PBDEs_WWRS = CN_PBDEs_WWRS[CN_PBDEs_WWRS['Unit'].str.contains('dw')]
CN_PBDEs_WWRS[['Min','Max','Mean','Median']] *= (1/1000)
CN_PBDEs_WWRS['valid_data'] = CN_PBDEs_WWRS.apply(add_data, axis=1)

WWRS_PBDEs = [x for sublist in CN_PBDEs_WWRS['valid_data'] for x in sublist]

# soil
set(CN_PBDEs_soil['Unit'])
CN_PBDEs_soil = CN_PBDEs_soil[CN_PBDEs_soil['Unit'].str.contains('dw')]
CN_PBDEs_soil[['Min','Max','Mean','Median']] *= (1/1000)
CN_PBDEs_soil['valid_data'] = CN_PBDEs_soil.apply(add_data, axis=1)

soil_PBDEs = [x for sublist in CN_PBDEs_soil['valid_data'] for x in sublist]

# sediment
set(CN_PBDEs_sediment['Unit'])
CN_PBDEs_sediment = CN_PBDEs_sediment[CN_PBDEs_sediment['Unit'].str.contains('dw')]
CN_PBDEs_sediment[['Min','Max','Mean','Median']] *= (1/1000)
CN_PBDEs_sediment['valid_data'] = CN_PBDEs_sediment.apply(add_data, axis=1)

sediment_PBDEs = [x for sublist in CN_PBDEs_sediment['valid_data'] for x in sublist]

# air
set(CN_PBDEs_air['Unit'])
CN_PBDEs_air = CN_PBDEs_air[CN_PBDEs_air['Unit'].isin(['ng/L','ng/kg','ng/m3'])]
# assume the density of air is 1.225 kg/m3
CN_PBDEs_air.loc[CN_PBDEs_air['Unit'] == 'ng/L', ['Min','Max','Mean','Median']] *= (1000/1.225/1000)
CN_PBDEs_air.loc[CN_PBDEs_air['Unit'] == 'ng/m3', ['Min','Max','Mean','Median']] *= (1/1.225/1000)
CN_PBDEs_air.loc[CN_PBDEs_air['Unit'] == 'ng/kg', ['Min','Max','Mean','Median']] *= (1/1000)
CN_PBDEs_air['valid_data'] = CN_PBDEs_air.apply(add_data, axis=1)

air_PBDEs = [x for sublist in CN_PBDEs_air['valid_data'] for x in sublist]

# water
set(CN_PBDEs_water['Unit'])
CN_PBDEs_water = CN_PBDEs_water[CN_PBDEs_water['Unit'] != 'ng/kg dw']
# assume the density of water is 1000 kg/m3
CN_PBDEs_water[['Min','Max','Mean','Median']] *= (1/1000)
CN_PBDEs_water['valid_data'] = CN_PBDEs_water.apply(add_data, axis=1)

water_PBDEs = [x for sublist in CN_PBDEs_water['valid_data'] for x in sublist]

#%% PBDEs concentrations visualization

print('\n' + str([len(i) for i in [WWRS_PBDEs, soil_PBDEs, sediment_PBDEs, air_PBDEs, water_PBDEs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-8, 10**6])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4, 10**6], fontname='Arial')

ax.set_ylabel('$\mathbf{PBDEs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4, 10**6], fontname='Arial')

bp = ax.boxplot([WWRS_PBDEs, soil_PBDEs, sediment_PBDEs, air_PBDEs, water_PBDEs],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PBDEs.pdf', transparent=True, bbox_inches='tight')

#%% PBDEs concentrations visualization (unique values)

print('\n' + str([len(set(i)) for i in [WWRS_PBDEs, soil_PBDEs, sediment_PBDEs, air_PBDEs, water_PBDEs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-8, 10**6])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4, 10**6], fontname='Arial')

ax.set_ylabel('$\mathbf{PBDEs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4, 10**6], fontname='Arial')

bp = ax.boxplot([list(set(WWRS_PBDEs)), list(set(soil_PBDEs)), list(set(sediment_PBDEs)), list(set(air_PBDEs)), list(set(water_PBDEs))],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PBDEs_unique.pdf', transparent=True, bbox_inches='tight')

#%% PCBs concentration data

CN_PCBs = pd.read_excel(folder + 'analyses/POPs_data/exposure data/Standardized data used for analysis and human health risk assessment/PCBs_analysis.xlsx')
set(CN_PCBs['Pathway'])
CN_PCBs = CN_PCBs[CN_PCBs['Pathway'].isin(['Water','air','others','sediment','soil','water'])]
CN_PCBs = CN_PCBs[((CN_PCBs['Pathway'] == 'others') & (CN_PCBs['Detailed_pathway'].str.contains('sludge|Sludge'))) | (CN_PCBs['Pathway'] != 'others')]
# remove entries containing 'sludge' but are not WWRS
set(CN_PCBs[(CN_PCBs['Pathway'] == 'others') & (CN_PCBs['Detailed_pathway'].str.contains('sludge'))]['Detailed_pathway'])
CN_PCBs = CN_PCBs[CN_PCBs['Detailed_pathway'] != 'sludge, pond sludge']
CN_PCBs = CN_PCBs[['Pathway','Unit','Min','Max','Mean','Median']]
CN_PCBs.replace('nd', 0, inplace=True)
for column in ['Min','Max','Mean','Median']:
    CN_PCBs[column] = CN_PCBs[column].apply(pd.to_numeric, errors='coerce')

CN_PCBs_WWRS = CN_PCBs[CN_PCBs['Pathway'] == 'others']
CN_PCBs_soil = CN_PCBs[CN_PCBs['Pathway'] == 'soil']
CN_PCBs_sediment = CN_PCBs[CN_PCBs['Pathway'] == 'sediment']
CN_PCBs_air = CN_PCBs[CN_PCBs['Pathway'] == 'air']
CN_PCBs_water = CN_PCBs[CN_PCBs['Pathway'].isin(['Water','water'])]

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['Min']) or not pd.isna(row['Max']):
        if not pd.isna(row['Min']):
            values.append(row['Min'])
        if not pd.isna(row['Max']):
            values.append(row['Max'])
    else:
        if not pd.isna(row['Mean']):
            values.append(row['Mean'])
    
    if not pd.isna(row['Median']):
        values.append(row['Median'])
    return values

# WWRS
set(CN_PCBs_WWRS['Unit'])
CN_PCBs_WWRS = CN_PCBs_WWRS[CN_PCBs_WWRS['Unit'].str.contains('dw')]
CN_PCBs_WWRS[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCBs_WWRS['valid_data'] = CN_PCBs_WWRS.apply(add_data, axis=1)

WWRS_PCBs = [x for sublist in CN_PCBs_WWRS['valid_data'] for x in sublist]

# soil
set(CN_PCBs_soil['Unit'])
CN_PCBs_soil = CN_PCBs_soil[CN_PCBs_soil['Unit'].str.contains('dw')]
CN_PCBs_soil[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCBs_soil['valid_data'] = CN_PCBs_soil.apply(add_data, axis=1)

soil_PCBs = [x for sublist in CN_PCBs_soil['valid_data'] for x in sublist]

# sediment
set(CN_PCBs_sediment['Unit'])
CN_PCBs_sediment = CN_PCBs_sediment[CN_PCBs_sediment['Unit'].str.contains('dw')]
CN_PCBs_sediment[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCBs_sediment['valid_data'] = CN_PCBs_sediment.apply(add_data, axis=1)

sediment_PCBs = [x for sublist in CN_PCBs_sediment['valid_data'] for x in sublist]

# air
set(CN_PCBs_air['Unit'])
# assume the density of air is 1.225 kg/m3
CN_PCBs_air.loc[CN_PCBs_air['Unit'] == 'ng/L', ['Min','Max','Mean','Median']] *= (1000/1.225/1000)
CN_PCBs_air.loc[CN_PCBs_air['Unit'].isin(['ng/Nm3','ng/Nm3 TEQ','ng/m3','ng/m3 TEQ']), ['Min','Max','Mean','Median']] *= (1/1.225/1000)
CN_PCBs_air.loc[CN_PCBs_air['Unit'].isin(['ng/kg','ng/kg TEQ']), ['Min','Max','Mean','Median']] *= (1/1000)
CN_PCBs_air['valid_data'] = CN_PCBs_air.apply(add_data, axis=1)

air_PCBs = [x for sublist in CN_PCBs_air['valid_data'] for x in sublist]

# water
set(CN_PCBs_water['Unit'])
CN_PCBs_water = CN_PCBs_water[CN_PCBs_water['Unit'] != 'ng/kg dw']
# assume the density of water is 1000 kg/m3
CN_PCBs_water[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCBs_water['valid_data'] = CN_PCBs_water.apply(add_data, axis=1)

water_PCBs = [x for sublist in CN_PCBs_water['valid_data'] for x in sublist]

#%% PCBs concentrations visualization

print('\n' + str([len(i) for i in [WWRS_PCBs, soil_PCBs, sediment_PCBs, air_PCBs, water_PCBs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-9, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-9, 10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PCBs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-9, 10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([WWRS_PCBs, soil_PCBs, sediment_PCBs, air_PCBs, water_PCBs],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PCBs.pdf', transparent=True, bbox_inches='tight')

#%% PCBs concentrations visualization (unique values)

print('\n' + str([len(set(i)) for i in [WWRS_PCBs, soil_PCBs, sediment_PCBs, air_PCBs, water_PCBs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-9, 10**5])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-9, 10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

ax.set_ylabel('$\mathbf{PCBs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-9, 10**-7, 10**-5, 10**-3, 10**-1, 10**1, 10**3, 10**5], fontname='Arial')

bp = ax.boxplot([list(set(WWRS_PCBs)), list(set(soil_PCBs)), list(set(sediment_PCBs)), list(set(air_PCBs)), list(set(water_PCBs))],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PCBs_unique.pdf', transparent=True, bbox_inches='tight')

#%% PCDD&Fs concentration data

CN_PCDDFs = pd.read_excel(folder + 'analyses/POPs_data/exposure data/Standardized data used for analysis and human health risk assessment/PCDD&Fs_analysis.xlsx')
set(CN_PCDDFs['Pathway'])
CN_PCDDFs = CN_PCDDFs[CN_PCDDFs['Pathway'].isin(['air','others','sediment','soil','water'])]
CN_PCDDFs = CN_PCDDFs[((CN_PCDDFs['Pathway'] == 'others') & (CN_PCDDFs['Detailed_pathway'].str.contains('sludge|Sludge'))) | (CN_PCDDFs['Pathway'] != 'others')]
# remove entries containing 'sludge' but are not WWRS
set(CN_PCDDFs[(CN_PCDDFs['Pathway'] == 'others') & (CN_PCDDFs['Detailed_pathway'].str.contains('sludge'))]['Detailed_pathway'])
CN_PCDDFs = CN_PCDDFs[CN_PCDDFs['Detailed_pathway'] != 'sludge, pond sludge']
CN_PCDDFs = CN_PCDDFs[['Pathway','Detailed_pathway','Unit','Min','Max','Mean','Median']]
CN_PCDDFs.replace('nd', 0, inplace=True)
for column in ['Min','Max','Mean','Median']:
    CN_PCDDFs[column] = CN_PCDDFs[column].apply(pd.to_numeric, errors='coerce')

CN_PCDDFs_WWRS = CN_PCDDFs[(CN_PCDDFs['Pathway'] == 'others') & (~CN_PCDDFs['Detailed_pathway'].isin(['sludge amended soils-cauliflower soil',
                                                                                                      'sludge amended soils-melon soil',
                                                                                                      'sludge amended soils-wheat field soil']))]
CN_PCDDFs_soil = CN_PCDDFs[(CN_PCDDFs['Pathway'] == 'soil') | ((CN_PCDDFs['Pathway'] == 'others') & (CN_PCDDFs['Detailed_pathway'].isin(['sludge amended soils-cauliflower soil',
                                                                                                                                         'sludge amended soils-melon soil',
                                                                                                                                         'sludge amended soils-wheat field soil'])))]
CN_PCDDFs_sediment = CN_PCDDFs[CN_PCDDFs['Pathway'] == 'sediment']
CN_PCDDFs_air = CN_PCDDFs[CN_PCDDFs['Pathway'] == 'air']
CN_PCDDFs_water = CN_PCDDFs[CN_PCDDFs['Pathway'] == 'water']

# always collect min, max, and median; collect average only if none of min and max are present
def add_data(row):
    values = []
    if not pd.isna(row['Min']) or not pd.isna(row['Max']):
        if not pd.isna(row['Min']):
            values.append(row['Min'])
        if not pd.isna(row['Max']):
            values.append(row['Max'])
    else:
        if not pd.isna(row['Mean']):
            values.append(row['Mean'])
    
    if not pd.isna(row['Median']):
        values.append(row['Median'])
    return values

# WWRS
set(CN_PCDDFs_WWRS['Unit'])
CN_PCDDFs_WWRS = CN_PCDDFs_WWRS[CN_PCDDFs_WWRS['Unit'].str.contains('dw')]
CN_PCDDFs_WWRS[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCDDFs_WWRS['valid_data'] = CN_PCDDFs_WWRS.apply(add_data, axis=1)

WWRS_PCDDFs = [x for sublist in CN_PCDDFs_WWRS['valid_data'] for x in sublist]

# soil
set(CN_PCDDFs_soil['Unit'])
CN_PCDDFs_soil = CN_PCDDFs_soil[CN_PCDDFs_soil['Unit'].str.contains('dw')]
CN_PCDDFs_soil[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCDDFs_soil['valid_data'] = CN_PCDDFs_soil.apply(add_data, axis=1)

soil_PCDDFs = [x for sublist in CN_PCDDFs_soil['valid_data'] for x in sublist]

# sediment
set(CN_PCDDFs_sediment['Unit'])
CN_PCDDFs_sediment = CN_PCDDFs_sediment[CN_PCDDFs_sediment['Unit'].str.contains('dw')]
CN_PCDDFs_sediment[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCDDFs_sediment['valid_data'] = CN_PCDDFs_sediment.apply(add_data, axis=1)

sediment_PCDDFs = [x for sublist in CN_PCDDFs_sediment['valid_data'] for x in sublist]

# air
set(CN_PCDDFs_air['Unit'])
# assume the density of air is 1.225 kg/m3
CN_PCDDFs_air.loc[CN_PCDDFs_air['Unit'] == 'ng/dm3', ['Min','Max','Mean','Median']] *= (1000/1.225/1000)
CN_PCDDFs_air.loc[CN_PCDDFs_air['Unit'].isin(['ng/Nm3','ng/Nm3 TEQ','ng/m3','ng/m3 TEQ']), ['Min','Max','Mean','Median']] *= (1/1.225/1000)
CN_PCDDFs_air.loc[CN_PCDDFs_air['Unit'] == 'pg WHO2005-TEQ/m3', ['Min','Max','Mean','Median']] *= (1/1.225/1000/1000)
CN_PCDDFs_air['valid_data'] = CN_PCDDFs_air.apply(add_data, axis=1)

air_PCDDFs = [x for sublist in CN_PCDDFs_air['valid_data'] for x in sublist]

# water
# assume the density of water is 1000 kg/m3
set(CN_PCDDFs_water['Unit'])
CN_PCDDFs_water[['Min','Max','Mean','Median']] *= (1/1000)
CN_PCDDFs_water['valid_data'] = CN_PCDDFs_water.apply(add_data, axis=1)

water_PCDDFs = [x for sublist in CN_PCDDFs_water['valid_data'] for x in sublist]

#%% PCDD&Fs concentrations visualization

print('\n' + str([len(i) for i in [WWRS_PCDDFs, soil_PCDDFs, sediment_PCDDFs, air_PCDDFs, water_PCDDFs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-10, 10**2])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-10, 10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4], fontname='Arial')

ax.set_ylabel('$\mathbf{PCDD&Fs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-10, 10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4], fontname='Arial')

bp = ax.boxplot([WWRS_PCDDFs, soil_PCDDFs, sediment_PCDDFs, air_PCDDFs, water_PCDDFs],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PCDD&Fs.pdf', transparent=True, bbox_inches='tight')

#%% PCDD&Fs concentrations visualization (unique values)

print('\n' + str([len(set(i)) for i in [WWRS_PCDDFs, soil_PCDDFs, sediment_PCDDFs, air_PCDDFs, water_PCDDFs]]) + '\n')

fig, ax = plt.subplots(figsize=(15, 5))

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([0.5, 5.5])
ax.set_ylim([10**-10, 10**2])
ax.tick_params(direction='inout', length=20, width=3, bottom=False, top=False, left=True, right=False)

plt.yticks([10**-10, 10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4], fontname='Arial')

ax.set_ylabel('$\mathbf{PCDD&Fs}$\n[ng·${g^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_xlim(ax.get_xlim())
ax_right.set_yscale('log')
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

ax_right.set_xticklabels(['WWRS','soil','sediment','air','water'])
plt.yticks([10**-10, 10**-8, 10**-6, 10**-4, 10**-2, 10**0, 10**2, 10**4], fontname='Arial')

bp = ax.boxplot([list(set(WWRS_PCDDFs)), list(set(soil_PCDDFs)), list(set(sediment_PCDDFs)), list(set(air_PCDDFs)), list(set(water_PCDDFs))],
                whis=[5, 95], showfliers=False, widths=0.7, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=da, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=la, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

plt.savefig('/Users/jiananfeng/Desktop/PCDD&Fs_unique.pdf', transparent=True, bbox_inches='tight')

#%% country average

filterwarnings('ignore')

C_cost_mean = []
C_CI_mean = []
T_cost_FOAK_mean = []
T_CI_FOAK_mean = []
T_cost_NOAK_mean = []
T_CI_NOAK_mean = []

C_cost_median = []
C_CI_median = []
T_cost_FOAK_median = []
T_CI_FOAK_median = []
T_cost_NOAK_median = []
T_CI_NOAK_median = []

C_cost_min = []
C_CI_min = []
T_cost_FOAK_min = []
T_CI_FOAK_min = []
T_cost_NOAK_min = []
T_CI_NOAK_min = []

C_cost_max = []
C_CI_max = []
T_cost_FOAK_max = []
T_CI_FOAK_max = []
T_cost_NOAK_max = []
T_CI_NOAK_max = []

# run in different consoles to speed up
# do not have more than two consoles running together since the application memory is not enough
for i in range(2500, len(WRRF_filtered)):
    if i%100 == 0:
        print(i)
    
    C_TEA = []
    C_LCA = []
    for function in (create_C1_system, create_C2_system, create_C3_system,
                     create_C4_system, create_C5_system, create_C6_system,
                     create_C7_system, create_C8_system, create_C9_system,
                     create_C10_system, create_C11_system, create_C12_system,
                     create_C13_system, create_C14_system, create_C15_system,
                     create_C16_system, create_C17_system, create_C18_system,
                     create_C19_system, create_C20_system, create_C21_system,
                     create_C22_system, create_C23_system, create_C24_system,
                     create_C25_system):
        sys = function(country_code=WRRF_filtered.iloc[i]['CNTRY_ISO'], size=WRRF_filtered.iloc[i]['dry_solids_tonne_per_day'])
        
        C_TEA.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        C_LCA.append(sys.LCA.get_total_impacts(operation_only=True,
                                               exclude=(sys.flowsheet.raw_wastewater,),
                                               annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    C_cost_mean.append(np.mean(C_TEA))
    C_CI_mean.append(np.mean(C_LCA))
    C_cost_median.append(np.median(C_TEA))
    C_CI_median.append(np.median(C_LCA))
    C_cost_min.append(np.min(C_TEA))
    C_CI_min.append(np.min(C_LCA))
    C_cost_max.append(np.max(C_TEA))
    C_CI_max.append(np.max(C_LCA))
    
    T_FOAK_TEA = []
    T_FOAK_LCA = []
    for function in (create_T1_system, create_T2_system, create_T3_system,
                     create_T4_system, create_T5_system, create_T6_system,
                     create_T7_system, create_T8_system, create_T9_system,
                     create_T10_system, create_T11_system, create_T12_system,
                     create_T13_system, create_T14_system, create_T15_system):
        sys = function(country_code=WRRF_filtered.iloc[i]['CNTRY_ISO'], size=WRRF_filtered.iloc[i]['dry_solids_tonne_per_day'], FOAK=True)
        
        T_FOAK_TEA.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        T_FOAK_LCA.append(sys.LCA.get_total_impacts(operation_only=True,
                                                    exclude=(sys.flowsheet.raw_wastewater,),
                                                    annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    T_cost_FOAK_mean.append(np.mean(T_FOAK_TEA))
    T_CI_FOAK_mean.append(np.mean(T_FOAK_LCA))
    T_cost_FOAK_median.append(np.median(T_FOAK_TEA))
    T_CI_FOAK_median.append(np.median(T_FOAK_LCA))
    T_cost_FOAK_min.append(np.min(T_FOAK_TEA))
    T_CI_FOAK_min.append(np.min(T_FOAK_LCA))
    T_cost_FOAK_max.append(np.max(T_FOAK_TEA))
    T_CI_FOAK_max.append(np.max(T_FOAK_LCA))
    
    T_NOAK_TEA = []
    T_NOAK_LCA = []
    for function in (create_T1_system, create_T2_system, create_T3_system,
                     create_T4_system, create_T5_system, create_T6_system,
                     create_T7_system, create_T8_system, create_T9_system,
                     create_T10_system, create_T11_system, create_T12_system,
                     create_T13_system, create_T14_system, create_T15_system):
        sys = function(country_code=WRRF_filtered.iloc[i]['CNTRY_ISO'], size=WRRF_filtered.iloc[i]['dry_solids_tonne_per_day'], FOAK=False)
        
        T_NOAK_TEA.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        T_NOAK_LCA.append(sys.LCA.get_total_impacts(operation_only=True,
                                                    exclude=(sys.flowsheet.raw_wastewater,),
                                                    annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    T_cost_NOAK_mean.append(np.mean(T_NOAK_TEA))
    T_CI_NOAK_mean.append(np.mean(T_NOAK_LCA))
    T_cost_NOAK_median.append(np.median(T_NOAK_TEA))
    T_CI_NOAK_median.append(np.median(T_NOAK_LCA))
    T_cost_NOAK_min.append(np.min(T_NOAK_TEA))
    T_CI_NOAK_min.append(np.min(T_NOAK_LCA))
    T_cost_NOAK_max.append(np.max(T_NOAK_TEA))
    T_CI_NOAK_max.append(np.max(T_NOAK_LCA))

country_mean = pd.DataFrame({'C_cost_mean': C_cost_mean,
                             'C_CI_mean': C_CI_mean,
                             'T_cost_FOAK_mean': T_cost_FOAK_mean,
                             'T_CI_FOAK_mean': T_CI_FOAK_mean,
                             'T_cost_NOAK_mean': T_cost_NOAK_mean,
                             'T_CI_NOAK_mean': T_CI_NOAK_mean})

country_median = pd.DataFrame({'C_cost_median': C_cost_median,
                               'C_CI_median': C_CI_median,
                               'T_cost_FOAK_median': T_cost_FOAK_median,
                               'T_CI_FOAK_median': T_CI_FOAK_median,
                               'T_cost_NOAK_median': T_cost_NOAK_median,
                               'T_CI_NOAK_median': T_CI_NOAK_median})

country_min = pd.DataFrame({'C_cost_min': C_cost_min,
                            'C_CI_min': C_CI_min,
                            'T_cost_FOAK_min': T_cost_FOAK_min,
                            'T_CI_FOAK_min': T_CI_FOAK_min,
                            'T_cost_NOAK_min': T_cost_NOAK_min,
                            'T_CI_NOAK_min': T_CI_NOAK_min})

country_max = pd.DataFrame({'C_cost_max': C_cost_max,
                            'C_CI_max': C_CI_max,
                            'T_cost_FOAK_max': T_cost_FOAK_max,
                            'T_CI_FOAK_max': T_CI_FOAK_max,
                            'T_cost_NOAK_max': T_cost_NOAK_max,
                            'T_CI_NOAK_max': T_CI_NOAK_max})

country_mean.to_excel(folder + f'results/country_mean_{date.today()}_{i}.xlsx')

country_median.to_excel(folder + f'results/country_median_{date.today()}_{i}.xlsx')

country_min.to_excel(folder + f'results/country_min_{date.today()}_{i}.xlsx')

country_max.to_excel(folder + f'results/country_max_{date.today()}_{i}.xlsx')

#%% merge country results - using 5 dry tonne/day as the cutoff size

# =============================================================================
# mean
# =============================================================================
country_mean_1 = pd.read_excel(folder + 'results/country_mean_2025-11-14_2499.xlsx')
country_mean_2 = pd.read_excel(folder + 'results/country_mean_2025-11-14_4726.xlsx')

integrated_country_mean = pd.concat([country_mean_1, country_mean_2])
integrated_country_mean.reset_index(inplace=True)
integrated_country_mean = integrated_country_mean[['C_cost_mean','C_CI_mean','T_cost_FOAK_mean','T_CI_FOAK_mean','T_cost_NOAK_mean','T_CI_NOAK_mean']]

WRRF_mean = pd.concat([WRRF_filtered, integrated_country_mean], axis=1)

WRRF_mean['C_cost_times_mass_flow'] = WRRF_mean['C_cost_mean']*WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['C_CI_times_mass_flow'] = WRRF_mean['C_CI_mean']*WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_cost_FOAK_times_mass_flow'] = WRRF_mean['T_cost_FOAK_mean']*WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_CI_FOAK_times_mass_flow'] = WRRF_mean['T_CI_FOAK_mean']*WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_cost_NOAK_times_mass_flow'] = WRRF_mean['T_cost_NOAK_mean']*WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_CI_NOAK_times_mass_flow'] = WRRF_mean['T_CI_NOAK_mean']*WRRF_mean['dry_solids_tonne_per_day']

WRRF_mean = WRRF_mean[['COUNTRY','dry_solids_tonne_per_day','C_cost_times_mass_flow','C_CI_times_mass_flow',
                       'T_cost_FOAK_times_mass_flow','T_CI_FOAK_times_mass_flow','T_cost_NOAK_times_mass_flow','T_CI_NOAK_times_mass_flow']]

# the world map data include 'Taiwan' as part of 'China'
WRRF_mean.loc[WRRF_mean['COUNTRY'] == 'Taiwan', 'COUNTRY'] = 'China'

WRRF_mean = WRRF_mean.groupby('COUNTRY').sum()

WRRF_mean['C_cost_weighted_average'] = WRRF_mean['C_cost_times_mass_flow']/WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['C_CI_weighted_average'] = WRRF_mean['C_CI_times_mass_flow']/WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_cost_FOAK_weighted_average'] = WRRF_mean['T_cost_FOAK_times_mass_flow']/WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_CI_FOAK_weighted_average'] = WRRF_mean['T_CI_FOAK_times_mass_flow']/WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_cost_NOAK_weighted_average'] = WRRF_mean['T_cost_NOAK_times_mass_flow']/WRRF_mean['dry_solids_tonne_per_day']
WRRF_mean['T_CI_NOAK_weighted_average'] = WRRF_mean['T_CI_NOAK_times_mass_flow']/WRRF_mean['dry_solids_tonne_per_day']

WRRF_mean = WRRF_mean[['C_cost_weighted_average','C_CI_weighted_average','T_cost_FOAK_weighted_average',
                       'T_CI_FOAK_weighted_average', 'T_cost_NOAK_weighted_average', 'T_CI_NOAK_weighted_average']]

WRRF_mean['carbon_credit_needed'] = (WRRF_mean['T_cost_NOAK_weighted_average'] - WRRF_mean['C_cost_weighted_average'])/(WRRF_mean['C_CI_weighted_average'] - WRRF_mean['T_CI_NOAK_weighted_average'])*1000

WRRF_mean.reset_index(inplace=True)

WRRF_mean = WRRF_mean.merge(WRRF_filtered[['COUNTRY','CNTRY_ISO']].drop_duplicates(), how='left', on='COUNTRY')

world_mean = world.merge(WRRF_mean, how='left', left_on='ISO_A3', right_on='CNTRY_ISO')

assert len(WRRF_mean) == len(world_mean[world_mean['CNTRY_ISO'].notna()]['ISO_A3'].drop_duplicates())

# world_mean.to_excel(folder + f'results/country_mean_{date.today()}_summary.xlsx')

# =============================================================================
# median
# =============================================================================
country_median_1 = pd.read_excel(folder + 'results/country_median_2025-11-14_2499.xlsx')
country_median_2 = pd.read_excel(folder + 'results/country_median_2025-11-14_4726.xlsx')

integrated_country_median = pd.concat([country_median_1, country_median_2])
integrated_country_median.reset_index(inplace=True)
integrated_country_median = integrated_country_median[['C_cost_median','C_CI_median','T_cost_FOAK_median','T_CI_FOAK_median','T_cost_NOAK_median','T_CI_NOAK_median']]

WRRF_median = pd.concat([WRRF_filtered, integrated_country_median], axis=1)

WRRF_median['C_cost_times_mass_flow'] = WRRF_median['C_cost_median']*WRRF_median['dry_solids_tonne_per_day']
WRRF_median['C_CI_times_mass_flow'] = WRRF_median['C_CI_median']*WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_cost_FOAK_times_mass_flow'] = WRRF_median['T_cost_FOAK_median']*WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_CI_FOAK_times_mass_flow'] = WRRF_median['T_CI_FOAK_median']*WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_cost_NOAK_times_mass_flow'] = WRRF_median['T_cost_NOAK_median']*WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_CI_NOAK_times_mass_flow'] = WRRF_median['T_CI_NOAK_median']*WRRF_median['dry_solids_tonne_per_day']

WRRF_median = WRRF_median[['COUNTRY','dry_solids_tonne_per_day','C_cost_times_mass_flow','C_CI_times_mass_flow',
                           'T_cost_FOAK_times_mass_flow','T_CI_FOAK_times_mass_flow','T_cost_NOAK_times_mass_flow','T_CI_NOAK_times_mass_flow']]

# the world map data include 'Taiwan' as part of 'China'
WRRF_median.loc[WRRF_median['COUNTRY'] == 'Taiwan', 'COUNTRY'] = 'China'

WRRF_median = WRRF_median.groupby('COUNTRY').sum()

WRRF_median['C_cost_weighted_average'] = WRRF_median['C_cost_times_mass_flow']/WRRF_median['dry_solids_tonne_per_day']
WRRF_median['C_CI_weighted_average'] = WRRF_median['C_CI_times_mass_flow']/WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_cost_FOAK_weighted_average'] = WRRF_median['T_cost_FOAK_times_mass_flow']/WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_CI_FOAK_weighted_average'] = WRRF_median['T_CI_FOAK_times_mass_flow']/WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_cost_NOAK_weighted_average'] = WRRF_median['T_cost_NOAK_times_mass_flow']/WRRF_median['dry_solids_tonne_per_day']
WRRF_median['T_CI_NOAK_weighted_average'] = WRRF_median['T_CI_NOAK_times_mass_flow']/WRRF_median['dry_solids_tonne_per_day']

WRRF_median = WRRF_median[['C_cost_weighted_average','C_CI_weighted_average','T_cost_FOAK_weighted_average',
                           'T_CI_FOAK_weighted_average', 'T_cost_NOAK_weighted_average', 'T_CI_NOAK_weighted_average']]

WRRF_median['carbon_credit_needed'] = (WRRF_median['T_cost_NOAK_weighted_average'] - WRRF_median['C_cost_weighted_average'])/(WRRF_median['C_CI_weighted_average'] - WRRF_median['T_CI_NOAK_weighted_average'])*1000

WRRF_median.reset_index(inplace=True)

WRRF_median = WRRF_median.merge(WRRF_filtered[['COUNTRY','CNTRY_ISO']].drop_duplicates(), how='left', on='COUNTRY')

world_median = world.merge(WRRF_median, how='left', left_on='ISO_A3', right_on='CNTRY_ISO')

assert len(WRRF_median) == len(world_median[world_median['CNTRY_ISO'].notna()]['ISO_A3'].drop_duplicates())

# world_median.to_excel(folder + f'results/country_median_{date.today()}_summary.xlsx')

# =============================================================================
# min
# =============================================================================
# note cost and CI minimums may not come from the same system (this is fine for identifying the possible range)
country_min_1 = pd.read_excel(folder + 'results/country_min_2025-11-14_2499.xlsx')
country_min_2 = pd.read_excel(folder + 'results/country_min_2025-11-14_4726.xlsx')

integrated_country_min = pd.concat([country_min_1, country_min_2])
integrated_country_min.reset_index(inplace=True)
integrated_country_min = integrated_country_min[['C_cost_min','C_CI_min','T_cost_FOAK_min','T_CI_FOAK_min','T_cost_NOAK_min','T_CI_NOAK_min']]

WRRF_min = pd.concat([WRRF_filtered, integrated_country_min], axis=1)

WRRF_min['C_cost_times_mass_flow'] = WRRF_min['C_cost_min']*WRRF_min['dry_solids_tonne_per_day']
WRRF_min['C_CI_times_mass_flow'] = WRRF_min['C_CI_min']*WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_cost_FOAK_times_mass_flow'] = WRRF_min['T_cost_FOAK_min']*WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_CI_FOAK_times_mass_flow'] = WRRF_min['T_CI_FOAK_min']*WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_cost_NOAK_times_mass_flow'] = WRRF_min['T_cost_NOAK_min']*WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_CI_NOAK_times_mass_flow'] = WRRF_min['T_CI_NOAK_min']*WRRF_min['dry_solids_tonne_per_day']

WRRF_min = WRRF_min[['COUNTRY','dry_solids_tonne_per_day','C_cost_times_mass_flow','C_CI_times_mass_flow',
                     'T_cost_FOAK_times_mass_flow','T_CI_FOAK_times_mass_flow','T_cost_NOAK_times_mass_flow','T_CI_NOAK_times_mass_flow']]

# the world map data include 'Taiwan' as part of 'China'
WRRF_min.loc[WRRF_min['COUNTRY'] == 'Taiwan', 'COUNTRY'] = 'China'

WRRF_min = WRRF_min.groupby('COUNTRY').sum()

WRRF_min['C_cost_weighted_average'] = WRRF_min['C_cost_times_mass_flow']/WRRF_min['dry_solids_tonne_per_day']
WRRF_min['C_CI_weighted_average'] = WRRF_min['C_CI_times_mass_flow']/WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_cost_FOAK_weighted_average'] = WRRF_min['T_cost_FOAK_times_mass_flow']/WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_CI_FOAK_weighted_average'] = WRRF_min['T_CI_FOAK_times_mass_flow']/WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_cost_NOAK_weighted_average'] = WRRF_min['T_cost_NOAK_times_mass_flow']/WRRF_min['dry_solids_tonne_per_day']
WRRF_min['T_CI_NOAK_weighted_average'] = WRRF_min['T_CI_NOAK_times_mass_flow']/WRRF_min['dry_solids_tonne_per_day']

WRRF_min = WRRF_min[['C_cost_weighted_average','C_CI_weighted_average','T_cost_FOAK_weighted_average',
                     'T_CI_FOAK_weighted_average', 'T_cost_NOAK_weighted_average', 'T_CI_NOAK_weighted_average']]

WRRF_min['carbon_credit_needed'] = (WRRF_min['T_cost_NOAK_weighted_average'] - WRRF_min['C_cost_weighted_average'])/(WRRF_min['C_CI_weighted_average'] - WRRF_min['T_CI_NOAK_weighted_average'])*1000

WRRF_min.reset_index(inplace=True)

WRRF_min = WRRF_min.merge(WRRF_filtered[['COUNTRY','CNTRY_ISO']].drop_duplicates(), how='left', on='COUNTRY')

world_min = world.merge(WRRF_min, how='left', left_on='ISO_A3', right_on='CNTRY_ISO')

assert len(WRRF_min) == len(world_min[world_min['CNTRY_ISO'].notna()]['ISO_A3'].drop_duplicates())

# world_min.to_excel(folder + f'results/country_min_{date.today()}_summary.xlsx')

# =============================================================================
# max
# =============================================================================
# note cost and CI maximums may not come from the same system (this is fine for identifying the possible range)
country_max_1 = pd.read_excel(folder + 'results/country_max_2025-11-14_2499.xlsx')
country_max_2 = pd.read_excel(folder + 'results/country_max_2025-11-14_4726.xlsx')

integrated_country_max = pd.concat([country_max_1, country_max_2])
integrated_country_max.reset_index(inplace=True)
integrated_country_max = integrated_country_max[['C_cost_max','C_CI_max','T_cost_FOAK_max','T_CI_FOAK_max','T_cost_NOAK_max','T_CI_NOAK_max']]

WRRF_max = pd.concat([WRRF_filtered, integrated_country_max], axis=1)

WRRF_max['C_cost_times_mass_flow'] = WRRF_max['C_cost_max']*WRRF_max['dry_solids_tonne_per_day']
WRRF_max['C_CI_times_mass_flow'] = WRRF_max['C_CI_max']*WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_cost_FOAK_times_mass_flow'] = WRRF_max['T_cost_FOAK_max']*WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_CI_FOAK_times_mass_flow'] = WRRF_max['T_CI_FOAK_max']*WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_cost_NOAK_times_mass_flow'] = WRRF_max['T_cost_NOAK_max']*WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_CI_NOAK_times_mass_flow'] = WRRF_max['T_CI_NOAK_max']*WRRF_max['dry_solids_tonne_per_day']

WRRF_max = WRRF_max[['COUNTRY','dry_solids_tonne_per_day','C_cost_times_mass_flow','C_CI_times_mass_flow',
                     'T_cost_FOAK_times_mass_flow','T_CI_FOAK_times_mass_flow','T_cost_NOAK_times_mass_flow','T_CI_NOAK_times_mass_flow']]

# the world map data include 'Taiwan' as part of 'China'
WRRF_max.loc[WRRF_max['COUNTRY'] == 'Taiwan', 'COUNTRY'] = 'China'

WRRF_max = WRRF_max.groupby('COUNTRY').sum()

WRRF_max['C_cost_weighted_average'] = WRRF_max['C_cost_times_mass_flow']/WRRF_max['dry_solids_tonne_per_day']
WRRF_max['C_CI_weighted_average'] = WRRF_max['C_CI_times_mass_flow']/WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_cost_FOAK_weighted_average'] = WRRF_max['T_cost_FOAK_times_mass_flow']/WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_CI_FOAK_weighted_average'] = WRRF_max['T_CI_FOAK_times_mass_flow']/WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_cost_NOAK_weighted_average'] = WRRF_max['T_cost_NOAK_times_mass_flow']/WRRF_max['dry_solids_tonne_per_day']
WRRF_max['T_CI_NOAK_weighted_average'] = WRRF_max['T_CI_NOAK_times_mass_flow']/WRRF_max['dry_solids_tonne_per_day']

WRRF_max = WRRF_max[['C_cost_weighted_average','C_CI_weighted_average','T_cost_FOAK_weighted_average',
                     'T_CI_FOAK_weighted_average', 'T_cost_NOAK_weighted_average', 'T_CI_NOAK_weighted_average']]

WRRF_max['carbon_credit_needed'] = (WRRF_max['T_cost_NOAK_weighted_average'] - WRRF_max['C_cost_weighted_average'])/(WRRF_max['C_CI_weighted_average'] - WRRF_max['T_CI_NOAK_weighted_average'])*1000

WRRF_max.reset_index(inplace=True)

WRRF_max = WRRF_max.merge(WRRF_filtered[['COUNTRY','CNTRY_ISO']].drop_duplicates(), how='left', on='COUNTRY')

world_max = world.merge(WRRF_max, how='left', left_on='ISO_A3', right_on='CNTRY_ISO')

assert len(WRRF_max) == len(world_max[world_max['CNTRY_ISO'].notna()]['ISO_A3'].drop_duplicates())

# world_max.to_excel(folder + f'results/country_max_{date.today()}_summary.xlsx')

#%% world map visualization - C cost mean

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 54.5
plt.rcParams['ytick.labelsize'] = 54.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', b, db])

world.plot(ax=ax, color='none', edgecolor='k', hatch='//', linewidth=1)

world_mean.plot(column='C_cost_weighted_average', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=-200, vmax=1000)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{Cost}$ [\$·${tonne^{−1}}$]', fontname='Arial', fontsize=54.5)
fig.axes[1].tick_params(length=15, width=3)
fig.axes[1].set_yticks([-200, 0, 200, 400, 600, 800, 1000])
fig.axes[1].set_yticklabels(['-200','0','200','400','600','800','1000'], fontname='Arial', fontsize=54.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.0575, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

# plt.savefig('/Users/jiananfeng/Desktop/C_cost.png', transparent=True, bbox_inches='tight')
plt.savefig('/Users/jiananfeng/Desktop/C_cost.pdf', transparent=True, bbox_inches='tight')

#%% world map visualization - T NOAK cost mean

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 54.5
plt.rcParams['ytick.labelsize'] = 54.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', b, db])

world.plot(ax=ax, color='none', edgecolor='k', hatch='//', linewidth=1)

world_mean.plot(column='T_cost_NOAK_weighted_average', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=-200, vmax=1000)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{Cost}$ [\$·${tonne^{−1}}$]', fontname='Arial', fontsize=54.5)
fig.axes[1].tick_params(length=15, width=3)
fig.axes[1].set_yticks([-200, 0, 200, 400, 600, 800, 1000])
fig.axes[1].set_yticklabels(['-200','0','200','400','600','800','1000'], fontname='Arial', fontsize=54.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.0575, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

# plt.savefig('/Users/jiananfeng/Desktop/T_NOAK_cost.png', transparent=True, bbox_inches='tight')
plt.savefig('/Users/jiananfeng/Desktop/T_NOAK_cost.pdf', transparent=True, bbox_inches='tight')

#%% world map visualization - C CI mean

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 54.5
plt.rcParams['ytick.labelsize'] = 54.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', g, dg])

world.plot(ax=ax, color='none', edgecolor='k', hatch='//', linewidth=1)

world_mean.plot(column='C_CI_weighted_average', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=-300, vmax=900)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{CI}$ [kg CO${_{2}}$e·${tonne^{−1}}$]', fontname='Arial', fontsize=54.5, labelpad=16)
fig.axes[1].tick_params(length=15, width=3)
fig.axes[1].set_yticks([-300, -150, 0, 150, 300, 450, 600, 750, 900])
fig.axes[1].set_yticklabels(['-300','-150','0','150','300','450','600','750','900'], fontname='Arial', fontsize=54.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.0575, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

# plt.savefig('/Users/jiananfeng/Desktop/C_CI.png', transparent=True, bbox_inches='tight')
plt.savefig('/Users/jiananfeng/Desktop/C_CI.pdf', transparent=True, bbox_inches='tight')

#%% world map visualization - T NOAK CI mean

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 54.5
plt.rcParams['ytick.labelsize'] = 54.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', g, dg])

world.plot(ax=ax, color='none', edgecolor='k', hatch='//', linewidth=1)

world_mean.plot(column='T_CI_NOAK_weighted_average', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=-300, vmax=900)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{CI}$ [kg CO${_{2}}$e·${tonne^{−1}}$]', fontname='Arial', fontsize=54.5, labelpad=16)
fig.axes[1].tick_params(length=15, width=3)
fig.axes[1].set_yticks([-300, -150, 0, 150, 300, 450, 600, 750, 900])
fig.axes[1].set_yticklabels(['-300','-150','0','150','300','450','600','750','900'], fontname='Arial', fontsize=54.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.0575, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

# plt.savefig('/Users/jiananfeng/Desktop/T_NOAK_CI.png', transparent=True, bbox_inches='tight')
plt.savefig('/Users/jiananfeng/Desktop/T_NOAK_CI.pdf', transparent=True, bbox_inches='tight')

#%% world map visualization - carbon credit needed

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 54.5
plt.rcParams['ytick.labelsize'] = 54.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', ['w', r, dr])

world.plot(ax=ax, color='none', edgecolor='k', hatch='//', linewidth=1)

world_mean.plot(column='carbon_credit_needed', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=0, vmax=600)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{Carbon\ credit}$' + '\n[\$·${tonne\ CO_{2}e^{−1}}$]', fontname='Arial', fontsize=54.5)
fig.axes[1].tick_params(length=15, width=3)
fig.axes[1].set_yticks([0, 100, 200, 300, 400, 500, 600])
fig.axes[1].set_yticklabels(['<0','100','200','300','400','500','600'], fontname='Arial', fontsize=54.5)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.0575, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

ax.set_aspect(1)

ax.set_axis_off()

# plt.savefig('/Users/jiananfeng/Desktop/NOAK_carbon_credit.png', transparent=True, bbox_inches='tight')
plt.savefig('/Users/jiananfeng/Desktop/NOAK_carbon_credit.pdf', transparent=True, bbox_inches='tight')

#%% country order

country_order = WRRF_filtered['COUNTRY'].value_counts()

# 'China' and 'Taiwan' data are merged
country_order['China'] += country_order['Taiwan']
country_order = country_order[country_order.index != 'Taiwan']

country_order = country_order[country_order > 10].index.tolist()[::-1]

#%% cost ranges for C and T systems

cost_min = WRRF_min.copy()
cost_min = cost_min.set_index('COUNTRY').loc[country_order]

assert (cost_min.index != country_order).sum() == 0

print([i for i in cost_min.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)])

global_north_countries_regions.append('Hong Kong')
global_north_countries_regions.append('Czech Republic')

assert(len([i for i in cost_min.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)]) == 0)

cost_min.loc[cost_min.index.isin(global_north_countries_regions), 'color'] = la
cost_min.loc[cost_min.index.isin(global_south_countries_regions), 'color'] = y

cost_max = WRRF_max.copy()
cost_max = cost_max.set_index('COUNTRY').loc[country_order]

assert (cost_max.index != country_order).sum() == 0

print([i for i in cost_max.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)])

global_north_countries_regions.append('Hong Kong')
global_north_countries_regions.append('Czech Republic')

assert(len([i for i in cost_max.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)]) == 0)

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 48
plt.rcParams['ytick.labelsize'] = 48
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(11.5, 50))

ax = plt.gca()

ax.set_xlim([-300, 1200])
ax.set_ylim([0.25, len(country_order) + 0.75])

ax.tick_params(axis='x', direction='inout', length=30, width=3, pad=5)
ax.tick_params(axis='y', direction='inout', length=30, width=3, pad=5)

index = np.arange(1, len(country_order)+1, 1)

ax.barh(y=index+0.1875,
        width=cost_max['C_cost_weighted_average'] - cost_min['C_cost_weighted_average'],
        height=0.375,
        color=lb,
        edgecolor='k',
        linewidth=3,
        left=cost_min['C_cost_weighted_average'])

ax.barh(y=index-0.1875,
        width=cost_max['T_cost_NOAK_weighted_average'] - cost_min['T_cost_NOAK_weighted_average'],
        height=0.375,
        color=db,
        edgecolor='k',
        linewidth=3,
        left=cost_min['T_cost_NOAK_weighted_average'])

plt.xticks(np.arange(-300, 1500, 300), fontname='Arial')
plt.yticks(index, cost_max['CNTRY_ISO'], fontname='Arial')

ax.set_xlabel('$\mathbf{Cost}$ [\$·${tonne^{-1}}$]',
              fontname='Arial', fontsize=48)

ax_top = ax.twiny()

ax_top.set_xlim([-300, 1200])

ax_top.tick_params(axis='x', direction='inout', length=30, width=3, pad=-5)

plt.xticks(np.arange(-300, 1500, 300), fontname='Arial')
plt.yticks(index, cost_max['CNTRY_ISO'], fontname='Arial')

plt.savefig('/Users/jiananfeng/Desktop/cost_ranges.pdf', transparent=True, bbox_inches='tight')

#%% cost changes - data preparation

cost_change = WRRF_mean.copy()
cost_change['cost_difference'] = cost_change['T_cost_NOAK_weighted_average'] - cost_change['C_cost_weighted_average']
cost_change = cost_change.set_index('COUNTRY').loc[country_order]

assert (cost_change.index != country_order).sum() == 0

print([i for i in cost_change.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)])

global_north_countries_regions.append('Hong Kong')
global_north_countries_regions.append('Czech Republic')

assert(len([i for i in cost_change.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)]) == 0)

cost_change.loc[cost_change.index.isin(global_north_countries_regions), 'color'] = la
cost_change.loc[cost_change.index.isin(global_south_countries_regions), 'color'] = y

#%% cost changes - visualization

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 48
plt.rcParams['ytick.labelsize'] = 48
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(11.5, 50))

ax = plt.gca()

ax.set_xlim([-100, 300])
ax.set_ylim([0.25, len(country_order) + 0.75])

ax.plot([0, 0],
        [0, 51],
        c='k',
        linewidth=3)

ax.tick_params(axis='x', direction='inout', length=30, width=3, pad=5)
ax.tick_params(axis='y', direction='inout', length=30, width=3, pad=5)

index = np.arange(1, len(country_order)+1, 1)

ax.barh(y=index,
        width=cost_change['cost_difference'],
        height=0.75,
        color=cost_change['color'],
        edgecolor='k',
        linewidth=3)

plt.xticks(np.arange(-100, 400, 100), fontname='Arial')
plt.yticks(index, cost_change['CNTRY_ISO'], fontname='Arial')

ax.set_xlabel('$\mathbf{Cost\ change}$ [\$·${tonne^{-1}}$]',
              fontname='Arial', fontsize=48)

ax_top = ax.twiny()

ax_top.set_xlim([-100, 300])

ax_top.tick_params(axis='x', direction='inout', length=30, width=3, pad=-5)

plt.xticks(np.arange(-100, 400, 100), fontname='Arial')
plt.yticks(index, cost_change['CNTRY_ISO'], fontname='Arial')

plt.savefig('/Users/jiananfeng/Desktop/cost_changess.pdf', transparent=True, bbox_inches='tight')

#%% CI ranges for C and T systems

CI_min = WRRF_min.copy()
CI_min = CI_min.set_index('COUNTRY').loc[country_order]

assert (CI_min.index != country_order).sum() == 0

print([i for i in CI_min.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)])

global_north_countries_regions.append('Hong Kong')
global_north_countries_regions.append('Czech Republic')

assert(len([i for i in CI_min.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)]) == 0)

CI_min.loc[CI_min.index.isin(global_north_countries_regions), 'color'] = la
CI_min.loc[CI_min.index.isin(global_south_countries_regions), 'color'] = y

CI_max = WRRF_max.copy()
CI_max = CI_max.set_index('COUNTRY').loc[country_order]

assert (CI_max.index != country_order).sum() == 0

print([i for i in CI_max.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)])

global_north_countries_regions.append('Hong Kong')
global_north_countries_regions.append('Czech Republic')

assert(len([i for i in CI_max.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)]) == 0)

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 48
plt.rcParams['ytick.labelsize'] = 48
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(11.5, 50))

ax = plt.gca()

ax.set_xlim([-3, 2.5])
ax.set_ylim([0.25, len(country_order) + 0.75])

ax.tick_params(axis='x', direction='inout', length=30, width=3, pad=5)
ax.tick_params(axis='y', direction='inout', length=30, width=3, pad=5)

index = np.arange(1, len(country_order)+1, 1)

ax.barh(y=index+0.1875,
        width=(CI_max['C_CI_weighted_average'] - CI_min['C_CI_weighted_average'])/1000,
        height=0.375,
        color=dg,
        edgecolor='k',
        linewidth=3,
        left=CI_min['C_CI_weighted_average']/1000)

ax.barh(y=index-0.1875,
        width=(CI_max['T_CI_NOAK_weighted_average'] - CI_min['T_CI_NOAK_weighted_average'])/1000,
        height=0.375,
        color=lg,
        edgecolor='k',
        linewidth=3,
        left=CI_min['T_CI_NOAK_weighted_average']/1000)

plt.xticks(np.arange(-3, 4, 1), fontname='Arial')
plt.yticks(index, CI_max['CNTRY_ISO'], fontname='Arial')

ax.set_xlabel('$\mathbf{CI}$ [tonne CO${_{2}}$e·${tonne^{-1}}$]',
              fontname='Arial', fontsize=48)

ax_top = ax.twiny()

ax_top.set_xlim([-3, 2.5])

ax_top.tick_params(axis='x', direction='inout', length=30, width=3, pad=-5)

plt.xticks(np.arange(-3, 4, 1), fontname='Arial')
plt.yticks(index, CI_max['CNTRY_ISO'], fontname='Arial')

plt.savefig('/Users/jiananfeng/Desktop/CI_ranges.pdf', transparent=True, bbox_inches='tight')

#%% CI changes - data preparation

CI_change = WRRF_mean.copy()
CI_change['CI_difference'] = CI_change['T_CI_NOAK_weighted_average'] - CI_change['C_CI_weighted_average']
CI_change = CI_change.set_index('COUNTRY').loc[country_order]

assert (CI_change.index != country_order).sum() == 0

print([i for i in CI_change.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)])

global_north_countries_regions.append('Hong Kong')
global_north_countries_regions.append('Czech Republic')

assert(len([i for i in CI_change.index if (i not in global_north_countries_regions) and (i not in global_south_countries_regions)]) == 0)

CI_change.loc[CI_change.index.isin(global_north_countries_regions), 'color'] = la
CI_change.loc[CI_change.index.isin(global_south_countries_regions), 'color'] = y

#%% CI changes - visualization

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 48
plt.rcParams['ytick.labelsize'] = 48
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(11.5, 50))

ax = plt.gca()

ax.set_xlim([-1200, 0])
ax.set_ylim([0.25, len(country_order) + 0.75])

ax.plot([0, 0],
        [0, 51],
        c='k',
        linewidth=3)

ax.tick_params(axis='x', direction='inout', length=30, width=3, pad=5)
ax.tick_params(axis='y', direction='inout', length=30, width=3, pad=5)

index = np.arange(1, len(country_order)+1, 1)

ax.barh(y=index,
        width=CI_change['CI_difference'],
        height=0.75,
        color=CI_change['color'],
        edgecolor='k',
        linewidth=3)

plt.xticks(np.arange(-1200, 300, 300), fontname='Arial')
plt.yticks(index, CI_change['CNTRY_ISO'], fontname='Arial')

ax.set_xlabel('$\mathbf{CI\ change}$ [kg CO${_{2}}$e·${tonne^{-1}}$]',
              fontname='Arial', fontsize=48)

ax_top = ax.twiny()

ax_top.set_xlim([-1200, 0])

ax_top.tick_params(axis='x', direction='inout', length=30, width=3, pad=-5)

plt.xticks(np.arange(-1200, 300, 300), fontname='Arial')
plt.yticks(index, CI_change['CNTRY_ISO'], fontname='Arial')

plt.savefig('/Users/jiananfeng/Desktop/CI_changes.pdf', transparent=True, bbox_inches='tight')

#%% unit capital cost equation

sys_size = [1, 5, 25, 125, 625]
scaling_stream = []
purchase_cost = []
# replace systems and units as needed
for i in sys_size:
    sys = create_C1_system(size=i)
    scaling_stream.append(sys.flowsheet.Thickening.ins[0].F_mass/1000)
    purchase_cost.append(sys.flowsheet.Thickening.purchase_cost)

log_scaling_stream = np.log(scaling_stream)
log_purchase_cost = np.log(purchase_cost)

slope, intercept, r_value, p_value, std_err = linregress(log_scaling_stream, log_purchase_cost)

a = np.exp(intercept)
b = slope

print(f'fitted model: cost = {a:.3g} * size^{b:.3g}')
print(f'R² = {r_value**2:.3g}')

#%% unit capital cost BM

# replace systems and units as needed
sys = create_C1_system(size=10)
purchase_cost = sys.flowsheet.Thickening.purchase_cost
installed_cost = sys.flowsheet.Thickening.installed_cost

BM = installed_cost/purchase_cost

print(f'BM: {BM:.3g}')

#%% global equity

# country_average_solids = WRRF.groupby('COUNTRY').mean('dry_solids_tonne_per_day')

# country_average_solids_HDI = country_average_solids.merge(HDI, how='inner', left_on='COUNTRY', right_on='Country')

# plt.rcParams['axes.linewidth'] = 3
# plt.rcParams['hatch.linewidth'] = 3
# plt.rcParams['xtick.labelsize'] = 36
# plt.rcParams['ytick.labelsize'] = 36
# plt.rcParams['font.sans-serif'] = 'Arial'
# plt.rcParams['pdf.fonttype'] = 42
# plt.rcParams['ps.fonttype'] = 42

# plt.rcParams.update({'mathtext.fontset':'custom'})
# plt.rcParams.update({'mathtext.default':'regular'})
# plt.rcParams.update({'mathtext.bf':'Arial: bold'})

# fig, ax = plt.subplots(figsize=(12, 10))

# ax = plt.gca()

# ax.set_xlim([0.35, 1])
# ax.set_ylim([0, 90])

# ax.tick_params(direction='inout', length=20, width=3,
#                bottom=True, top=False, left=True, right=False, pad=0)

# plt.xticks(np.arange(0.4, 1.1, 0.1), fontname='Arial')
# plt.yticks(np.arange(0, 100, 20), fontname='Arial')

# ax_right = ax.twinx()
# ax_right.set_ylim(ax.get_ylim())
# plt.xticks(np.arange(0.4, 1.1, 0.1), fontname='Arial')
# plt.yticks(np.arange(0, 100, 20), fontname='Arial')

# ax_right.tick_params(direction='in', length=10, width=3,
#                      bottom=False, top=False, left=False, right=True,
#                      labelcolor='none')

# ax_top = ax.twiny()
# ax_top.set_xlim(ax.get_xlim())
# plt.xticks(np.arange(0.4, 1.1, 0.1), fontname='Arial')
# plt.yticks(np.arange(0, 100, 20), fontname='Arial')

# ax_top.tick_params(direction='in', length=10, width=3,
#                    bottom=False, top=True, left=False, right=False,
#                    labelcolor='none')

# ax.set_xlabel('$\mathbf{HDI}$',
#               fontname='Arial',
#               fontsize=45)

# ax.set_ylabel('$\mathbf{Average\ mass\ flow\ rate}$\n[dry tonne·${day^{-1}}$]',
#               fontname='Arial',
#               fontsize=45)

# # very high (≥0.8), high (0.7-0.8), medium (0.55-0.7), low (<0.55)
# country_average_solids_HDI['HDI_color'] = country_average_solids_HDI['HDI'].apply(lambda x: y if x < 0.55 else (o if x < 0.7 else (r if x < 0.8 else dr)))

# plt.scatter(country_average_solids_HDI['HDI'], country_average_solids_HDI['dry_solids_tonne_per_day'], color=country_average_solids_HDI['HDI_color'], edgecolor='k', linewidths=3, s=200)

# ax.plot([0.55, 0.55],
#         [0, 100],
#         c='k',
#         linestyle='--',
#         linewidth=3)

# ax.plot([0.7, 0.7],
#         [0, 100],
#         c='k',
#         linestyle='--',
#         linewidth=3)

# ax.plot([0.8, 0.8],
#         [0, 100],
#         c='k',
#         linestyle='--',
#         linewidth=3)

#%% writing

WRRF_mean['C_cost_weighted_average'].min()
WRRF_mean['C_cost_weighted_average'].max()

WRRF_mean['T_cost_NOAK_weighted_average'].min()
WRRF_mean['T_cost_NOAK_weighted_average'].max()

plt.boxplot([electricity_price[electricity_price['country'].isin(global_north_countries_regions)]['US_cents_per_kWh'], electricity_price[electricity_price['country'].isin(global_south_countries_regions)]['US_cents_per_kWh']])
plt.boxplot([labor_cost[labor_cost['country'].isin(global_north_countries_regions)]['average_annual_income_USD'], labor_cost[labor_cost['country'].isin(global_south_countries_regions)]['average_annual_income_USD']])
plt.boxplot([PLI[PLI['country'].isin(global_north_countries_regions)]['PLI'], PLI[PLI['country'].isin(global_south_countries_regions)]['PLI']])

cost_change = WRRF_mean.copy()
cost_change['cost_difference'] = cost_change['T_cost_NOAK_weighted_average'] - cost_change['C_cost_weighted_average']
cost_change[cost_change['COUNTRY'].isin(global_north_countries_regions)]['cost_difference'].mean()
cost_change[cost_change['COUNTRY'].isin(['India','Egypt','Peru','Thailand','Indonesia'])]

CI_change = WRRF_mean.copy()
CI_change['CI_difference'] = CI_change['T_CI_NOAK_weighted_average'] - CI_change['C_CI_weighted_average']
CI_change['CI_difference'].min()
CI_change['CI_difference'].max()

WRRF_mean['carbon_credit_needed'].min()
WRRF_mean['carbon_credit_needed'].max()

WRRF_mean[WRRF_mean['carbon_credit_needed']<0]