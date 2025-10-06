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

#%% housekeeping

# TODO: use Matplotlib.rcParams['pdf.fonttype'] = 42

#%% initialization

import numpy as np, pandas as pd, geopandas as gpd, matplotlib.pyplot as plt, matplotlib.colors as colors, chaospy as cp, json
from matplotlib.animation import FuncAnimation
from matplotlib.mathtext import _mathtext as mathtext
from colorpalette import Color
from matplotlib.colors import to_hex
from chaospy import distributions as shape
from scipy.stats import qmc
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

mathtext.FontConstantsBase.sup1 = 0.35

folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_landscape/'

_MMgal_to_m3 = auom('gal').conversion_factor('m3')*1000000
# TODO: this can be highly uncertain
# tonne/MG
ww_2_dry_solids = 1

# color palette
b = Color('blue', (96, 193, 207)).HEX
g = Color('green', (121, 191, 130)).HEX
r = Color('red', (237, 88, 111)).HEX
o = Color('orange', (249, 143, 96)).HEX
y = Color('yellow', (243, 195, 84)).HEX
a = Color('gray', (144, 145, 142)).HEX
p = Color('purple', (162, 128, 185)).HEX

lb = to_hex((96/256, 193/256, 207/256, 0.5), keep_alpha=True)
lg = to_hex((121/256, 191/256, 130/256, 0.5), keep_alpha=True)
lr = to_hex((237/256, 88/256, 111/256, 0.5), keep_alpha=True)
lo = to_hex((249/256, 143/256, 96/256, 0.5), keep_alpha=True)
ly = to_hex((243/256, 195/256, 84/256, 0.5), keep_alpha=True)
la = to_hex((144/256, 145/256, 142/256, 0.5), keep_alpha=True)
lp = to_hex((162/256, 128/256, 185/256, 0.5), keep_alpha=True)

db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

world = gpd.read_file(folder + 'analyses/World Bank Official Boundaries - Admin 0_all_layers/WB_GAD_ADM0_complete.shp')

WRRF = pd.read_csv(folder + 'analyses/HydroWASTE_v10/HydroWASTE_v10.csv', encoding='latin-1')
# WASTE_DIS in m3/day
WRRF = WRRF[WRRF['WASTE_DIS'] != 0]
WRRF = WRRF[~WRRF['STATUS'].isin(['Closed','Decommissioned','Non-Operational'])]
WRRF['dry_solids_tonne_per_day'] = WRRF['WASTE_DIS']/_MMgal_to_m3*ww_2_dry_solids
# TODO: add this to the manuscript or SI
# TODO: to demonstrate this, may look at the ITO dataset, see how much percentage of WRRF with the WWRF mass flow < 1 dry tonne/day is lagoon
# models will not be accurate when the WWRS mass flow rate is < 1 dry tonne/day (in reality, those WRRFs may use other WRRS treatment methods, e.g., lagoon, wetland, etc.)
# only keep WRRFs with ≥ 1 tonne dry solids per day
WRRF_filtered = WRRF[WRRF['dry_solids_tonne_per_day'] >= 1]
WRRF_filtered.reset_index(inplace=True)
# TODO: add this to the manuscript or SI
print(f"{WRRF_filtered['dry_solids_tonne_per_day'].sum()/WRRF['dry_solids_tonne_per_day'].sum()*100:.1f}% global WWRS are included." )
print(f"{len(set(WRRF_filtered['COUNTRY']))} countires are included." )

WRRF_filtered = gpd.GeoDataFrame(WRRF_filtered, crs='EPSG:4269',
                                 geometry=gpd.points_from_xy(x=WRRF_filtered.LON_WWTP,
                                                             y=WRRF_filtered.LAT_WWTP))

HDI = pd.read_excel(folder + 'analyses/HDR25_Statistical_Annex_HDI_Table.xlsx')

# TODO: add this to the manuscript or SI
# electricity price
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

# TODO: add this to the manuscript or SI
# =============================================================================
# code for processing ecoinvent data
# data = pd.read_excel('/Users/jiananfeng/Desktop/electricity_CI_more_digits.xlsx')
# data['country_code'] = data['country'].str.extract(r'\((.*?)\)')
# data['country'] = data['country'].str.replace(r"\s*\(.*?\)", "", regex=True)
# data.to_excel('/Users/jiananfeng/Desktop/electricity_CI_more_digits.xlsx')
# =============================================================================
# TODO: for countries w/o country data, use regional or global data
# electricity CI
# electriicty_CI = 

# TODO: add this to the manuscript or SI
# labor cost
# TODO: for countries with no data, consider (1) manually searching for the data, or (2) use the average value of the labor costs of countries with similar development stages
labor_cost = pd.read_excel(folder + 'analyses/EAR_4HRL_SEX_CUR_NB_A-20250917T1926.xlsx')
labor_cost = labor_cost[labor_cost['sex.label'] == 'Total']
labor_cost = labor_cost[labor_cost['classif1.label'] == 'Currency: U.S. dollars']
labor_cost = labor_cost.loc[labor_cost.groupby('ref_area.label')['time'].idxmax()].reset_index(drop=True)
# note obs_value can be 0 if the wage after being converted to the U.S. dollar is too low
labor_cost = labor_cost[labor_cost['obs_value'].notna()]
labor_cost = labor_cost[['ref_area.label','obs_value']]
labor_cost.rename(columns={'ref_area.label':'country',
                           'obs_value':'labor_cost'},
                  inplace=True)

# TODO: add this to the manuscript or SI
# price level index (PLI) - for chemical price adjustment
PLI = pd.read_excel(folder + 'analyses/P_Data_Extract_From_World_Development_Indicators.xlsx')
year_cols = ['1990 [YR1990]','2000 [YR2000]','2015 [YR2015]','2016 [YR2016]','2017 [YR2017]',
             '2018 [YR2018]','2019 [YR2019]','2020 [YR2020]','2021 [YR2021]','2022 [YR2022]',
             '2023 [YR2023]','2024 [YR2024]']
PLI[year_cols] = PLI[year_cols].replace('..', np.nan)
PLI['PLI'] = PLI[year_cols].ffill(axis=1).iloc[:, -1]
PLI = PLI[PLI['PLI'].notna()]
PLI = PLI[['Country Name','Country Code','PLI']]
PLI.rename(columns={'Country Name':'country',
                    'Country Code':'country_code'},
           inplace=True)

# TODO: index for capital costs (any others, like transportation and utilities, but these two are less important and may use other indexes or just the same globally - if other indexes are not very different, so these remains not important, like if capital cost in country A is 0.01% of U.S., in this case, if still keep transportation cost the U.S. value, the transporation may become the most impactful driver)

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

ax.set_ylabel('$\mathbf{C}$ [particle·${g^{-1}}$]',
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

# plt.savefig('/Users/jiananfeng/Desktop/MPs.pdf', transparent=True, bbox_inches='tight')

#%% MPs distribution

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
ax.set_ylim([0, 100])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.set_ylabel('$\mathbf{Passage}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.bar(0.5,
       100,
       width=0.7,
       color=lr,
       edgecolor='k',
       linewidth=3)

ax.bar(0.5,
       np.quantile(MPs_WWRS_MC, 0.5)*100,
       yerr=np.array([[(np.quantile(MPs_WWRS_MC, 0.5) - np.quantile(MPs_WWRS_MC, 0.05))*100], [(np.quantile(MPs_WWRS_MC, 0.95) - np.quantile(MPs_WWRS_MC, 0.5))*100]]),
       error_kw=dict(capsize=15, lw=3, capthick=3),
       width=0.7,
       color=dr,
       edgecolor='k',
       linewidth=3)

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
# TODO: include all kinds of water for now (except sewage-related water)
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
# TODO: just remove 'mg/kg dry-weight' and assume others are roughly to be the same as per volume
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

# plt.savefig('/Users/jiananfeng/Desktop/PhACs.pdf', transparent=True, bbox_inches='tight')

#%% PhACs distribution

unused_PhACs = shape.Uniform(0.15, 0.98)
take_back = shape.Triangle(0.136, 0.169, 0.203)
toilet_trash_ratio = shape.Triangle(0.319, 0.399, 0.479)
toilet_to_WRRF = shape.Uniform(0.398, 0.596)
# trash_to_combustion = shape.Uniform(0.152, 0.228)

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
ax.set_ylim([0, 100])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.set_ylabel('$\mathbf{Passage}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.bar(0.5,
       100,
       width=0.7,
       color=lo,
       edgecolor='k',
       linewidth=3)

ax.bar(0.5,
       np.quantile(PhACs_WWRS_MC, 0.5)*100,
       yerr=np.array([[(np.quantile(PhACs_WWRS_MC, 0.5) - np.quantile(PhACs_WWRS_MC, 0.05))*100], [(np.quantile(PhACs_WWRS_MC, 0.95) - np.quantile(PhACs_WWRS_MC, 0.5))*100]]),
       error_kw=dict(capsize=15, lw=3, capthick=3),
       width=0.7,
       color=do,
       edgecolor='k',
       linewidth=3)

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
    
bp['boxes'][0].set(color='k', facecolor=dg, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=lg, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=lg, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=lg, linewidth=3)
bp['boxes'][4].set(color='k', facecolor=lg, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

# plt.savefig('/Users/jiananfeng/Desktop/PhACs.pdf', transparent=True, bbox_inches='tight')

#%% ARGs distribution

human_use = shape.Uniform(0.224, 0.336)
human_excretion = shape.Triangle(0.02, 0.39, 0.862)
excretion_to_WRRF = shape.Uniform(0.398, 0.596)
animal_excretion = shape.Triangle(0.02, 0.39, 0.862)
WWRS_capture_rate = shape.Uniform(0.90, 0.95)

joint = cp.distributions.J(human_use, human_excretion, excretion_to_WRRF,
                           animal_excretion, WWRS_capture_rate)

sample = joint.sample(100000)

ARCs_WWRS_MC = pd.DataFrame()

ARCs_WRRS = sample[0]*sample[1]*sample[2]*sample[4]

# TODO: assume ARGs can not be degraded during wastewater treatment processes
ARCs_human_to_environment = sample[0]*sample[1]
ARCs_animal_to_environment = (1 - sample[0])*sample[3]

ARCs_WWRS_MC = ARCs_WRRS/(ARCs_human_to_environment + ARCs_animal_to_environment)

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
ax.set_ylim([0, 100])
ax.tick_params(direction='inout', length=20, width=3, labelbottom=False, bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.set_ylabel('$\mathbf{Passage}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=13)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(direction='in', length=10, width=3, bottom=False, top=False, left=False, right=True, labelcolor='none')

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax.bar(0.5,
       100,
       width=0.7,
       color=lg,
       edgecolor='k',
       linewidth=3)

ax.bar(0.5,
       np.quantile(ARCs_WWRS_MC, 0.5)*100,
       yerr=np.array([[(np.quantile(ARCs_WWRS_MC, 0.5) - np.quantile(ARCs_WWRS_MC, 0.05))*100], [(np.quantile(ARCs_WWRS_MC, 0.95) - np.quantile(ARCs_WWRS_MC, 0.5))*100]]),
       error_kw=dict(capsize=15, lw=3, capthick=3),
       width=0.7,
       color=dg,
       edgecolor='k',
       linewidth=3)

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

# plt.savefig('/Users/jiananfeng/Desktop/PFOA.pdf', transparent=True, bbox_inches='tight')

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

# plt.savefig('/Users/jiananfeng/Desktop/PFOS.pdf', transparent=True, bbox_inches='tight')

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
# TODO: check if it is reasonable to remove 'ng/kg dw''
CN_PBDEs_air = CN_PBDEs_air[CN_PBDEs_air['Unit'] != 'ng/kg dw']
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

# plt.savefig('/Users/jiananfeng/Desktop/PBDEs.pdf', transparent=True, bbox_inches='tight')

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

# plt.savefig('/Users/jiananfeng/Desktop/PCBs.pdf', transparent=True, bbox_inches='tight')

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

# plt.savefig('/Users/jiananfeng/Desktop/PCDD&Fs.pdf', transparent=True, bbox_inches='tight')

#%% country average

# TODO: as a kind of uncertainty analysis, show the minimum values on maps in the SI

# TODO: assign electricity price and CI (any other parameters?) for each country

filterwarnings('ignore')

C_cost_mean = []
C_CI_mean = []
T_cost_FOAK_mean = []
T_CI_FOAK_mean = []
T_cost_NOAK_mean = []
T_CI_NOAK_mean = []

# run in different consoles to speed up
# do not round all together since the application memory is not enough
# first round: 3000 and 6000
# second round: 9000 and 12000
# last round: 14000 and len(WRRF_filtered)
for i in range(0, len(WRRF_filtered)):
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
        sys = function(size=WRRF_filtered.iloc[i]['dry_solids_tonne_per_day'])
        
        C_TEA.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        C_LCA.append(sys.LCA.get_total_impacts(operation_only=True,
                                               exclude=(sys.flowsheet.raw_wastewater,),
                                               annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    C_cost_mean.append(np.median(C_TEA))
    C_CI_mean.append(np.median(C_LCA))
    
    T_FOAK_TEA = []
    T_FOAK_LCA = []
    for function in (create_T1_system, create_T2_system, create_T3_system,
                     create_T4_system, create_T5_system, create_T6_system,
                     create_T7_system, create_T8_system, create_T9_system,
                     create_T10_system, create_T11_system, create_T12_system,
                     create_T13_system, create_T14_system, create_T15_system):
        sys = function(size=WRRF_filtered.iloc[i]['dry_solids_tonne_per_day'], FOAK=True)
        
        T_FOAK_TEA.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        T_FOAK_LCA.append(sys.LCA.get_total_impacts(operation_only=True,
                                                    exclude=(sys.flowsheet.raw_wastewater,),
                                                    annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    T_cost_FOAK_mean.append(np.median(T_FOAK_TEA))
    T_CI_FOAK_mean.append(np.median(T_FOAK_LCA))
    
    T_NOAK_TEA = []
    T_NOAK_LCA = []
    for function in (create_T1_system, create_T2_system, create_T3_system,
                     create_T4_system, create_T5_system, create_T6_system,
                     create_T7_system, create_T8_system, create_T9_system,
                     create_T10_system, create_T11_system, create_T12_system,
                     create_T13_system, create_T14_system, create_T15_system):
        sys = function(size=WRRF_filtered.iloc[i]['dry_solids_tonne_per_day'], FOAK=False)
        
        T_NOAK_TEA.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        T_NOAK_LCA.append(sys.LCA.get_total_impacts(operation_only=True,
                                                    exclude=(sys.flowsheet.raw_wastewater,),
                                                    annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    T_cost_NOAK_mean.append(np.median(T_NOAK_TEA))
    T_CI_NOAK_mean.append(np.median(T_NOAK_LCA))

country_result = pd.DataFrame({'C_cost_mean': C_cost_mean,
                               'C_CI_mean': C_CI_mean,
                               'T_cost_FOAK_mean': T_cost_FOAK_mean,
                               'T_CI_FOAK_mean': T_CI_FOAK_mean,
                               'T_cost_NOAK_mean': T_cost_NOAK_mean,
                               'T_CI_NOAK_mean': T_CI_NOAK_mean,})

country_result.to_excel(folder + f'results/country_results_{date.today()}_{i}.xlsx')

#%% merge country results

country_result_1 = pd.read_excel(folder + 'results/country_results_2025-09-15_2999.xlsx')
country_result_2 = pd.read_excel(folder + 'results/country_results_2025-09-15_5999.xlsx')
country_result_3 = pd.read_excel(folder + 'results/country_results_2025-09-15_8999.xlsx')
country_result_4 = pd.read_excel(folder + 'results/country_results_2025-09-15_11999.xlsx')
country_result_5 = pd.read_excel(folder + 'results/country_results_2025-09-15_13999.xlsx')
country_result_6 = pd.read_excel(folder + 'results/country_results_2025-09-15_15964.xlsx')

integrated_country_result = pd.concat([country_result_1, country_result_2, country_result_3, country_result_4, country_result_5, country_result_6])
integrated_country_result.reset_index(inplace=True)
integrated_country_result = integrated_country_result[['C_cost_mean','C_CI_mean','T_cost_FOAK_mean','T_CI_FOAK_mean','T_cost_NOAK_mean','T_CI_NOAK_mean']]

WRRF_result = pd.concat([WRRF_filtered, integrated_country_result], axis=1)

WRRF_result['C_cost_times_mass_flow'] = WRRF_result['C_cost_mean']*WRRF_result['dry_solids_tonne_per_day']
WRRF_result['C_CI_times_mass_flow'] = WRRF_result['C_CI_mean']*WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_cost_FOAK_times_mass_flow'] = WRRF_result['T_cost_FOAK_mean']*WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_CI_FOAK_times_mass_flow'] = WRRF_result['T_CI_FOAK_mean']*WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_cost_NOAK_times_mass_flow'] = WRRF_result['T_cost_NOAK_mean']*WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_CI_NOAK_times_mass_flow'] = WRRF_result['T_CI_NOAK_mean']*WRRF_result['dry_solids_tonne_per_day']

WRRF_result = WRRF_result[['COUNTRY','dry_solids_tonne_per_day','C_cost_times_mass_flow','C_CI_times_mass_flow','T_cost_FOAK_times_mass_flow','T_CI_FOAK_times_mass_flow','T_cost_NOAK_times_mass_flow','T_CI_NOAK_times_mass_flow']]

WRRF_result = WRRF_result.groupby('COUNTRY').sum()

WRRF_result['C_cost_weighted_average'] = WRRF_result['C_cost_times_mass_flow']/WRRF_result['dry_solids_tonne_per_day']
WRRF_result['C_CI_weighted_average'] = WRRF_result['C_CI_times_mass_flow']/WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_cost_FOAK_weighted_average'] = WRRF_result['T_cost_FOAK_times_mass_flow']/WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_CI_FOAK_weighted_average'] = WRRF_result['T_CI_FOAK_times_mass_flow']/WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_cost_NOAK_weighted_average'] = WRRF_result['T_cost_NOAK_times_mass_flow']/WRRF_result['dry_solids_tonne_per_day']
WRRF_result['T_CI_NOAK_weighted_average'] = WRRF_result['T_CI_NOAK_times_mass_flow']/WRRF_result['dry_solids_tonne_per_day']

WRRF_result = WRRF_result[['C_cost_weighted_average','C_CI_weighted_average','T_cost_FOAK_weighted_average',
                           'T_CI_FOAK_weighted_average', 'T_cost_NOAK_weighted_average', 'T_CI_NOAK_weighted_average']]

WRRF_result.reset_index(inplace=True)

WRRF_result = WRRF_result.merge(WRRF_filtered[['COUNTRY','CNTRY_ISO']].drop_duplicates(), how='left', on='COUNTRY')

world_result = world.merge(WRRF_result, how='left', left_on='ISO_A3', right_on='CNTRY_ISO')

#%% world map visualization - C cost

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

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [b, g, y, o, r])

world_result.plot(column='C_cost_weighted_average', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=0, vmax=3500)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{Cost}$\n[\$·${tonne^{−1}}$]', fontname='Arial', fontsize=35)
fig.axes[1].tick_params(length=10, width=3)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% world map visualization - T NOAK cost

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

fig, ax = plt.subplots(figsize=(30, 30))

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [b, g, y, o, r])

world_result.plot(column='T_cost_NOAK_weighted_average', ax=ax, legend=True, legend_kwds={'shrink': 0.35}, cmap=color_map_Guest, vmin=0, vmax=3500)

world.plot(ax=ax, color='none', edgecolor='k', linewidth=1)

fig.axes[1].set_ylabel('$\mathbf{Cost}$\n[\$·${tonne^{−1}}$]', fontname='Arial', fontsize=35)
fig.axes[1].tick_params(length=10, width=3)

pos1 = fig.axes[1].get_position()
pos2 = [pos1.x0-0.035, pos1.y0, pos1.width, pos1.height] 
fig.axes[1].set_position(pos2)

# comment out the following line if the colorbar is needed
# fig.delaxes(fig.axes[1])

ax.set_aspect(1)

ax.set_axis_off()

#%% global equity

country_average_solids = WRRF.groupby('COUNTRY').mean('dry_solids_tonne_per_day')

country_average_solids_HDI = country_average_solids.merge(HDI, how='inner', left_on='COUNTRY', right_on='Country')

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

fig, ax = plt.subplots(figsize=(12, 10))

ax = plt.gca()

ax.set_xlim([0.35, 1])
ax.set_ylim([0, 90])

ax.tick_params(direction='inout', length=20, width=3,
               bottom=True, top=False, left=True, right=False, pad=0)

plt.xticks(np.arange(0.4, 1.1, 0.1), fontname='Arial')
plt.yticks(np.arange(0, 100, 20), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0.4, 1.1, 0.1), fontname='Arial')
plt.yticks(np.arange(0, 100, 20), fontname='Arial')

ax_right.tick_params(direction='in', length=10, width=3,
                     bottom=False, top=False, left=False, right=True,
                     labelcolor='none')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0.4, 1.1, 0.1), fontname='Arial')
plt.yticks(np.arange(0, 100, 20), fontname='Arial')

ax_top.tick_params(direction='in', length=10, width=3,
                   bottom=False, top=True, left=False, right=False,
                   labelcolor='none')

ax.set_xlabel('$\mathbf{HDI}$',
              fontname='Arial',
              fontsize=45)

ax.set_ylabel('$\mathbf{Average\ mass\ flow\ rate}$\n[dry tonne·${day^{-1}}$]',
              fontname='Arial',
              fontsize=45)

# very high (≥0.8), high (0.7-0.8), medium (0.55-0.7), low (<0.55)
country_average_solids_HDI['HDI_color'] = country_average_solids_HDI['HDI'].apply(lambda x: y if x < 0.55 else (o if x < 0.7 else (r if x < 0.8 else dr)))

plt.scatter(country_average_solids_HDI['HDI'], country_average_solids_HDI['dry_solids_tonne_per_day'], color=country_average_solids_HDI['HDI_color'], edgecolor='k', linewidths=3, s=200)

ax.plot([0.55, 0.55],
        [0, 100],
        c='k',
        linestyle='--',
        linewidth=3)

ax.plot([0.7, 0.7],
        [0, 100],
        c='k',
        linestyle='--',
        linewidth=3)

ax.plot([0.8, 0.8],
        [0, 100],
        c='k',
        linestyle='--',
        linewidth=3)

# TODO: update
# plt.savefig('/Users/jiananfeng/Desktop/Figure_XXX.pdf', transparent=True, bbox_inches='tight')

#%% WWTP random sampling

# TODO: adjust typology if including AD and AeD (any any units before them)
# all use uniform distribution (not Monte Carlo)
# for now, assume WRRFs are only responsible for the transportation of e.g., solids (and ash) to landfills, biocrude to oil refineries
# for other products, assume they can be sold at the gate of WRRFs
typology_bounds = {
    # tonne/day
    'solids_mass_flow': (0, 250),
    # dry weight ratio
    'solids_dw_ash': (0.15, 0.45),
    # ash free dry weight ratio
    'solids_afdw_lipid': (0.05, 0.35),
    # ash free dry weight ratio
    'solids_afdw_protein': (0.3, 0.55),
    # $/kWh
    'electricity_cost': (0.05, 0.2),
    # kg CO2 eq/kWh
    'electricity_CI': (0, 1.5),
    # km
    'landfill_distance': (0, 1500),
    # km
    'biofuel_distance': (0, 1500)
    }

typology_keys = list(typology_bounds.keys())

# Latin Hypercube Sampling
sampler = qmc.LatinHypercube(d=len(typology_keys), seed=3221)
# TODO: decide the sample size to balance the sample representativeness and the computation burden
LHS_unit = sampler.random(n=10000)

typology_mins = np.array([typology_bounds[key][0] for key in typology_keys])
typology_maxs = np.array([typology_bounds[key][1] for key in typology_keys])
LHS_samples = qmc.scale(LHS_unit, typology_mins, typology_maxs)

typology = {key: LHS_samples[:, i] for i, key in enumerate(typology_keys)}

# the biochemical composition parameters here not only affect HTL/HALT, but other through affecting the volatile solids content
lipid = typology['solids_afdw_lipid']
protein = typology['solids_afdw_protein']
carbohydrate = 1 - lipid - protein
# no negative carb
mask = carbohydrate >= 0
typology = {k: v[mask] for k, v in typology.items()}
typology['solids_afdw_carbohydrate'] = carbohydrate[mask]

for i in range(len(LHS_samples)):
    solids_mass_flow = typology['solids_mass_flow'][i]
    solids_dw_ash = typology['solids_dw_ash'][i]
    solids_afdw_lipid = typology['solids_afdw_lipid'][i]
    solids_afdw_protein = typology['solids_afdw_protein'][i]
    solids_afdw_carbohydrate = typology['solids_afdw_carbohydrate'][i]
    electricity_cost = typology['electricity_cost'][i]
    electricity_CI = typology['electricity_CI'][i]
    landfill_distance = typology['landfill_distance'][i]
    biofuel_distance = typology['biofuel_distance'][i]
    
    # TODO: compare the baseline cost and GHG results from different systems
    # TODO: for each typology, identify a winner
    # TODO: plot 5th and 95th percentiles to represent the opportunity space

#%% preliminary utility check

results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
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
                 create_T15_system):
    sys = function(size=10)
    print('\n' + sys.ID)
    print(round(sys.get_cooling_duty()))
    print(round(sys.get_heating_duty()))
    print(sys.heat_utilities)

#%% preliminary TEA check

TEA_results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
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
                 create_T15_system):
    sys = function(size=10)
    
    print('\n' + sys.ID)
    
    print(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    
    TEA_results.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)

for function in (create_T1_system, create_T2_system, create_T3_system,
                 create_T4_system, create_T5_system, create_T6_system,
                 create_T7_system, create_T8_system, create_T9_system,
                 create_T10_system, create_T11_system, create_T12_system,
                 create_T13_system, create_T14_system, create_T15_system):
    sys = function(size=10, FOAK=False)
    
    print('\n' + sys.ID + ' nth plant')
    
    print(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    
    TEA_results.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)

print(np.quantile(TEA_results[0:25], 0.05))
print(np.quantile(TEA_results[0:25], 0.5))
print(np.quantile(TEA_results[0:25], 0.95))

print(np.quantile(TEA_results[25:40], 0.05))
print(np.quantile(TEA_results[25:40], 0.5))
print(np.quantile(TEA_results[25:40], 0.95))

print(np.quantile(TEA_results[40:55], 0.05))
print(np.quantile(TEA_results[40:55], 0.5))
print(np.quantile(TEA_results[40:55], 0.95))

plt.boxplot((TEA_results[0:25], TEA_results[25:40], TEA_results[40:55]), whis=[5, 95], showfliers=False)

#%% preliminary LCA check

LCA_results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
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
                 create_T15_system):
    sys = function(size=10)
    
    print('\n' + sys.ID)
    
    print(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    LCA_results.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

for function in (create_T1_system, create_T2_system, create_T3_system,
                 create_T4_system, create_T5_system, create_T6_system,
                 create_T7_system, create_T8_system, create_T9_system,
                 create_T10_system, create_T11_system, create_T12_system,
                 create_T13_system, create_T14_system, create_T15_system):
    sys = function(size=10, FOAK=False)
    
    print('\n' + sys.ID + ' nth plant')
    
    print(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    LCA_results.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

print(np.quantile(LCA_results[0:25], 0.05))
print(np.quantile(LCA_results[0:25], 0.5))
print(np.quantile(LCA_results[0:25], 0.95))

print(np.quantile(LCA_results[25:40], 0.05))
print(np.quantile(LCA_results[25:40], 0.5))
print(np.quantile(LCA_results[25:40], 0.95))

print(np.quantile(LCA_results[40:55], 0.05))
print(np.quantile(LCA_results[40:55], 0.5))
print(np.quantile(LCA_results[40:55], 0.95))

# TODO: LCA_results[25:40] and LCA_results[40:55] may have differences after the parameter plant_performance_factor is implemented
plt.boxplot((LCA_results[0:25], LCA_results[25:40], LCA_results[40:55]), whis=[5, 95], showfliers=False)

#%% preliminary TEA and LCA animation

filterwarnings('ignore')

system_ID_list = []
system_type_list = []
size_list = []
TEA_result_list = []
LCA_result_list = []

for function in (create_C1_system, create_C2_system, create_C3_system,
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
                 create_T15_system):
    for size in np.linspace(1, 50, 50):
        sys = function(size=size)
        
        if size == 1:
            print(sys.ID)
        
        system_ID_list.append(sys.ID)
        system_type_list.append(sys.ID[7])
        size_list.append(size)
        TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                     exclude=(sys.flowsheet.raw_wastewater,),
                                                     annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

animation_data = pd.DataFrame({'system_ID': system_ID_list,
                               'system_type': system_type_list,
                               'size': size_list,
                               'cost': TEA_result_list,
                               'CI': LCA_result_list})

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

fig, ax = plt.subplots(figsize=(14, 12))

ax = plt.gca()
ax.set_xlim(0, 4000)
ax.set_ylim(-1000, 2500)

ax.tick_params(direction='inout', length=20, width=3,
                bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 4500, 500), fontname='Arial')
plt.yticks(np.arange(-1000, 3000, 500), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 4500, 500), fontname='Arial')

ax_top.tick_params(direction='in', length=10, width=3,
                   bottom=False, top=True, left=False, right=False,
                   labelcolor='none')

ax_bottom = ax.twinx()
ax_bottom.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0, 4500, 500), fontname='Arial')

ax_bottom.tick_params(direction='in', length=10, width=3,
                      bottom=False, top=False, left=False, right=True,
                      labelcolor='none')

ax.set_xlabel(r'$\mathbf{Cost}$ ' + '[\$·day${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{CI}$ ' + '[tonne CO${_2}$ eq·day${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

scatter_C = ax.scatter(x=[], y=[], s=500, color=a, linewidths=2, edgecolors='k')
scatter_T = ax.scatter(x=[], y=[], s=500, color=r, linewidths=2, edgecolors='k')

flow_scenarios = animation_data['size'].unique()

def update(frame):
    fr = flow_scenarios[frame]
    animation_data_fr = animation_data[animation_data['size'] == fr]
    danimation_data_C = animation_data_fr[animation_data_fr['system_ID'].str.contains('C')]
    animation_data_T = animation_data_fr[animation_data_fr['system_ID'].str.contains('T')]

    scatter_C.set_offsets(list(zip(danimation_data_C['cost'], danimation_data_C['CI'])))
    scatter_T.set_offsets(list(zip(animation_data_T['cost'], animation_data_T['CI'])))
    
    ax.set_title(f'{int(fr)} ' + 'dry tonne·day${^{-1}}$', fontname='Arial', fontsize=45)
    return scatter_C, scatter_T

ani = FuncAnimation(fig, update, frames=len(flow_scenarios), interval=800, blit=True)

fig.tight_layout()

ani.save('preliminary_TEA_LCA_results.mp4', writer="ffmpeg", fps=2, savefig_kwargs={'bbox_inches':'tight'})

plt.show()