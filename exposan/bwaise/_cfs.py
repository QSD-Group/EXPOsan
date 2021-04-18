#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# %%

# =============================================================================
# Database download, only run this once to download the database
# =============================================================================

# from bw2qsd import DataDownloader

# downloader = DataDownloader()
# downloader.download_ecoinvent() # apos 3.7.1



# %%

import os
c_path = os.path.realpath(__file__)
raw_path = os.path.join(os.path.split(c_path)[0], 'data/raw')

from bw2qsd import CFgetter

apos371 = CFgetter('apos371')
# ecoinvent version 3.7.1, at the point of substitution
apos371.load_database('ecoinvent_apos371')

# No need for Midpoint as Endpoint gives categorized results in addition to totals
apos371.load_indicators(add=True, method=('recipe'), method_exclude=('LT', 'obsolete'))
apos371.load_indicators(add=True, method=('recipe endpoint'), method_exclude='LT')
# For comparison with Trimmer et al.
apos371.load_indicators(add=True, method='TRACI')


# %%

import pandas as pd
all_acts = {}

def new_act(name):
    act = apos371.copy(name)
    all_acts[name] = act
    return act

def get_stats(df):
    df2 = df.copy()
    df2 = df2.append(df[1:].min(), ignore_index=True)
    df2 = df2.append(df[1:].mean(), ignore_index=True)
    df2 = df2.append(df[1:].median(), ignore_index=True)
    df2 = df2.append(df[1:].max(), ignore_index=True)
    df2.index = pd.Index((*df.index, 'min', 'mean', 'median', 'max'))
    df2.loc['mean'][0] = df2.loc['median'][0] = df2.loc['max'][0]
    return df2

brick = new_act('brick')
brick.load_activities('market brick', add=True, filter={'location': 'GLO'}, mask={'product': 'facility'}, limit=None)
brick.load_activities('market brick', add=True, filter={'location': 'RoW', 'product': 'brick'}, limit=None)

cement = new_act('cement')
cement.load_activities('market cement, unspecified', add=True, filter={'location': 'GLO'}, limit=None)
cement.load_activities('market cement, unspecified', add=True, filter={'location': 'RoW'}, limit=None)

concrete = new_act('concrete')
concrete.load_activities('market concrete, normal', add=True, filter={'location': 'GLO'}, limit=None)
concrete.load_activities('market concrete, normal', add=True, filter={'location': 'RoW'}, limit=None)

excavation = new_act('excavation')
excavation.load_activities('market excavation', add=True, limit=None)

gravel_sand = new_act('gravel_sand')
gravel_sand.load_activities('market gravel', add=True, filter={'product': 'gravel'}, mask={'name': 'infrastructure'}, limit=None)
gravel_sand.load_activities('market sand', add=True, filter={'product': 'sand', 'location': 'GLO'}, mask={'name': 'infrastructure'}, limit=None)
gravel_sand.load_activities('market sand', add=True, filter={'product': 'sand', 'location': 'RoW'}, mask={'name': 'infrastructure'}, limit=None)

hdpe_liner = new_act('hdpe_liner')
hdpe_liner.load_activities('hdpe', add=True, filter={'product': 'polyethylene', 'location': 'RoW'}, limit=None)
# hdpe_liner.load_activities('hdpe', add=True, filter={'product': 'polyethylene', 'location': 'GLO'}, limit=None) # None

reinforcing_steel = new_act('reinforcing_steel')
reinforcing_steel.load_activities('market reinforcing steel', add=True, filter={'product': 'steel'})
reinforcing_steel.load_activities('production reinforcing steel', add=True, filter={'product': 'reinforcing'})

steel = new_act('steel')
steel.load_activities('market steel, unalloyed', add=True, limit=None)
steel.load_activities('market steel, low-alloyed', add=True, limit=None)

stainless_steel = new_act('stainless_steel')
stainless_steel.load_activities('market steel, chromium steel 18/8', add=True, limit=None)

steel_rolling = new_act('steel_rolling')
steel_rolling.load_activities('market sheet rolling, chromium steel', add=True, limit=None)

wood = new_act('wood')
# wood.load_activities('market sawnwood', add=True, filter={'location': 'GLO'}, mask={'product': 'raw'}, limit=None) # None
wood.load_activities('market sawnwood', add=True, filter={'location': 'RoW'}, mask={'product': 'raw'}, limit=None)
to_remove = []
for act in wood.activities.keys():
    if ('azobe' in act) or ('paraná pine' in act) or ('lath' in act) or ('beam' in act) or ('board' in act):
        to_remove.append(act)
wood.remove('activity', to_remove)


# %%

for k, v in all_acts.items():
    get_stats(v.CFs).to_excel(os.path.join(raw_path, f'{k}_CFs.xlsx'))


# df = apos371.get_CF(path=os.path.join(raw_path, 'temp.xlsx'))

# Export impact indicators to be imported by QSDsan
# apos371.export_indicators(show=False, path=os.path.join(raw_path, 'temp.tsv'))
# apos371.export_indicators(show=False, path=os.path.join(raw_path, 'impact_indicators.tsv'))



# %%

# # Only run this at the very end to remove the outdated setup.pickle file
from bw2qsd import remove_setups_pickle
remove_setups_pickle()