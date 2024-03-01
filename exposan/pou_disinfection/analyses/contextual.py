#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, pandas as pd, geopandas as gpd
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from qsdsan import PowerUtility, ImpactItem
from qsdsan.utils import load_data
from exposan.pou_disinfection import (
    create_system,
    data_path,
    figures_path,
    ppl,
    results_path,
    update_number_of_householdsize,
    )

# Data for contextual analysis
water_path = os.path.join(data_path, '_raw_water.xlsx')
contextual_data = load_data(water_path, sheet='contextual')

# Function to update data and record results
def run_contextual_analysis():
    # Create systems
    sysA = create_system('A')
    sysB = create_system('B')
    sysC = create_system('C')
    sysD = create_system('D')
    
    all_syses = (sysA, sysB, sysC, sysD)
    for sys in all_syses:
        sys.units[0].Mg = 0 # assume all hardness is Ca
        
    e_item = ImpactItem.get_item('Electricity')
    
    all_results = {}
    # Sort by metric type
    for metric in ('cost', 'gwp'):
        for ID in (sys.ID for sys in all_syses):
            all_results[f'{ID}_{metric}'] = {}
    
    # # Sort by country
    # for ID in (sys.ID for sys in all_syses):
    #     all_results[f'{ID}_cost'] = {}
    #     all_results[f'{ID}_gwp'] = {}
    
    for community, data in contextual_data.iterrows():
        PowerUtility.price = data.e_cost
        e_item.CFs['GWP'] = data.e_gwp
        for sys in all_syses:
            update_number_of_householdsize(sys, household_size=data.household_size, ppl=ppl)
            RawWater = sys.units[0]
            RawWater.Ecoli = data.e_coli
            RawWater.turbidity = data.turbidity
            RawWater.Ca = data.hardness # assume all hardness is Ca
            sys.simulate()
            all_results[f'{sys.ID}_cost'][community] = sys.TEA.EAC/ppl
            all_results[f'{sys.ID}_gwp'][community] = sys.LCA.total_impacts['GWP']/(sys.LCA.lifetime*ppl)
    
    results = pd.DataFrame.from_dict(all_results)
    results.to_csv(os.path.join(results_path, 'contextual_analysis.csv'))


# %%

def plot_world_map():
    colors = ['w', colormaps['viridis'].colors[-1]]
    cmap = LinearSegmentedColormap.from_list('contextual', colors)
    
    # Natural Earth shape file, downloaded 2023-08-02
    # https://www.naturalearthdata.com/downloads/110m-cultural-vectors/
    world_folder = os.path.join(data_path, 'ne_110m_admin_0_countries')
    world_shape = gpd.read_file(os.path.join(world_folder, 'ne_110m_admin_0_countries.shp'))

    # Extract only useful columns
    world_trimmed = gpd.GeoDataFrame(geometry=world_shape['geometry'].copy())
    for column in ['NAME', 'ISO_A2', 'CONTINENT']:
        world_trimmed[column] = world_shape[column].copy()
    
    world_trimmed['latitude'] = world_shape['LABEL_Y'].copy()
    world_trimmed['longitude'] = world_shape['LABEL_X'].copy()
    
    # Color the country of interest
    world_trimmed = world_trimmed.set_index('ISO_A2')
    world_trimmed['used'] = 0
    world_trimmed.loc[world_trimmed.NAME.isin(contextual_data.country), 'used'] = 1
    
    ax = world_trimmed.plot(
        column='used',
        figsize=(25, 20),
        # legend=True,
        cmap=cmap,
        vmax=1,
        vmin=0,
        edgecolor='k',
        )
    
    fig = ax.figure
    ax.axis('off')
    fig.tight_layout(pad=0)
    fig.savefig(os.path.join(figures_path, 'world_map.png'), dpi=300)
    
    
# %%

if __name__ == '__main__':
    run_contextual_analysis()
    # plot_world_map()