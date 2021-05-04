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

import os, pickle
import pandas as pd
import qsdsan as qs
from bw2qsd import CFgetter, remove_setups_pickle
from bw2qsd.utils import format_name

c_path = os.path.dirname(__file__)
data_path = os.path.join(c_path, 'data')

__all__ = ('load_lca_data',)


# %%

# =============================================================================
# Database download, only run this once to download the database
# =============================================================================

def download_data():
    from bw2qsd import DataDownloader
    downloader = DataDownloader()
    downloader.download_ecoinvent() # apos 3.7.1


# %%

# =============================================================================
# Impact indicators
# =============================================================================

def create_indicators():
    apos371 = CFgetter('apos371')
    # ecoinvent version 3.7.1, at the point of substitution
    apos371.load_database('ecoinvent_apos371')

    # No need for Midpoint as Endpoint gives categorized results in addition to totals
    # apos371.load_indicators(add=True, method=('recipe'), method_exclude=('LT', 'obsolete'))
    apos371.load_indicators(add=True, method=('recipe endpoint'), method_exclude='LT')
    # For comparison with Trimmer et al.
    apos371.load_indicators(add=True, method='TRACI')

    ind_df_raw = apos371.export_indicators(show=False, path='')
    ind_df_processed = ind_df_raw.copy()

    for num, ind in ind_df_raw.iterrows():
        old_name = ind['indicator']
        if 'TRACI' in ind['method']:
            new_name = old_name
        else:
            kind = ind['method'].split(',A)')[0][-1] # egalitarian, hierarchist, individualist
            if old_name == 'Total':
                if ind['category'] == 'total':
                    new_name = f'{kind}_All' # total of all categories in E/H/I
                else:
                    new_name = f'{kind}_{format_name(ind["category"])}_Total' # category total
                    ind_df_processed.iloc[num]['category'] += ' total'
            else:
                new_name = f'{kind}_{old_name}'

        ind_df_processed.iloc[num]['indicator'] = new_name

    ind_df_processed.sort_values(by=['method', 'category', 'indicator'], inplace=True)

    qs.ImpactIndicator.load_indicators_from_file(ind_df_processed)
    indicators = qs.ImpactIndicator.get_all_indicators()

    return apos371, ind_df_processed, indicators


# %%

# =============================================================================
# Impact items
# =============================================================================

def select_items(database):
    all_acts = {}

    def new_act(database, all_acts, name):
        act = database.copy(name)
        all_acts[name] = act
        return act

    brick = new_act(database, all_acts, 'brick')
    brick.load_activities('market brick', add=True, filter={'location': 'GLO'}, mask={'product': 'facility'}, limit=None)
    brick.load_activities('market brick', add=True, filter={'location': 'RoW', 'product': 'brick'}, limit=None)

    cement = new_act(database, all_acts, 'cement')
    cement.load_activities('market cement, unspecified', add=True, filter={'location': 'GLO'}, limit=None)
    cement.load_activities('market cement, unspecified', add=True, filter={'location': 'RoW'}, limit=None)

    concrete = new_act(database, all_acts, 'concrete')
    concrete.load_activities('market concrete, normal', add=True, filter={'location': 'GLO'}, limit=None)
    concrete.load_activities('market concrete, normal', add=True, filter={'location': 'RoW'}, limit=None)

    excavation = new_act(database, all_acts, 'excavation')
    excavation.load_activities('market excavation', add=True, limit=None)

    gravel = new_act(database, all_acts, 'gravel')
    gravel.load_activities('market gravel', add=True, filter={'product': 'gravel'}, mask={'name': 'infrastructure'}, limit=None)

    hdpe_liner = new_act(database, all_acts, 'hdpe_liner')
    hdpe_liner.load_activities('hdpe', add=True, filter={'product': 'polyethylene', 'location': 'RoW'}, limit=None)
    # hdpe_liner.load_activities('hdpe', add=True, filter={'product': 'polyethylene', 'location': 'GLO'}, limit=None) # None

    # Not used
    # reinforcing_steel = new_act(database, all_acts, 'reinforcing_steel')
    # reinforcing_steel.load_activities('market reinforcing steel', add=True, filter={'product': 'steel'})
    # reinforcing_steel.load_activities('production reinforcing steel', add=True, filter={'product': 'reinforcing'})

    sand = new_act(database, all_acts, 'sand')
    sand.load_activities('market sand', add=True, filter={'product': 'sand', 'location': 'GLO'}, mask={'name': 'infrastructure'}, limit=None)
    sand.load_activities('market sand', add=True, filter={'product': 'sand', 'location': 'RoW'}, mask={'name': 'infrastructure'}, limit=None)

    stainless_steel = new_act(database, all_acts, 'stainless_steel')
    stainless_steel.load_activities('market steel, chromium steel 18/8', add=True, limit=None)

    steel_rolling = new_act(database, all_acts, 'steel_rolling')
    steel_rolling.load_activities('market sheet rolling, chromium steel', add=True, limit=None)

    steel = new_act(database, all_acts, 'steel')
    steel.load_activities('market steel, unalloyed', add=True, limit=None)
    steel.load_activities('market steel, low-alloyed', add=True, limit=None)

    trucking = new_act(database, all_acts, 'trucking')
    trucking.load_activities('market for transport, freight, lorry, all sizes', add=True, filter={'location': 'Row'}, limit=None)
    # trucking.load_activities('market for transport, freight, lorry, all sizes', add=False, filter={'location': 'Row'}, limit=None) # None

    wood = new_act(database, all_acts, 'wood')
    wood.load_activities('market sawnwood', add=True, filter={'location': 'RoW'}, mask={'product': 'raw'}, limit=None)
    # wood.load_activities('market sawnwood', add=True, filter={'location': 'GLO'}, mask={'product': 'raw'}, limit=None) # None
    to_remove = []
    for act in wood.activities.keys():
        if ('azobe' in act) or ('paraná pine' in act) or ('lath' in act) or ('beam' in act) or ('board' in act):
            to_remove.append(act)
    wood.remove('activity', to_remove)

    return all_acts


def get_stats(df, keep_raw_data=False):
    df2 = pd.DataFrame(columns=df.columns, dtype='float64') if not keep_raw_data else df.copy()
    df2 = df2.append(df[1:].min(), ignore_index=True)
    df2 = df2.append(df[1:].mean(), ignore_index=True)
    # df2 = df2.append(df[1:].median(), ignore_index=True)
    df2 = df2.append(df[1:].max(), ignore_index=True)
    df2.loc[-3:, ('-', '-', 'activity name')] = ('min', 'mean', 'max')
    functional_unit = df.loc[1, ('-', '-', 'functional unit')]
    df2.loc[1:, ('-', '-', 'functional unit')] = functional_unit
    return df2

def organize_cfs(all_acts):
    cf_dct = {}
    cf_dct['Brick'] = get_stats(all_acts['brick'].CFs)
    cf_dct['Cement'] = get_stats(all_acts['cement'].CFs)
    cf_dct['Concrete'] = get_stats(all_acts['concrete'].CFs)
    cf_dct['Excavation'] = get_stats(all_acts['excavation'].CFs)
    cf_dct['Gravel'] = get_stats(all_acts['gravel'].CFs)
    cf_dct['Plastic'] = get_stats(all_acts['hdpe_liner'].CFs)
    cf_dct['Sand'] = get_stats(all_acts['sand'].CFs)
    cf_dct['StainlessSteel'] = get_stats(all_acts['stainless_steel'].CFs)
    # Sum of stainless steel and steel rolling
    cf_dct['StainlessSteelSheet'] = get_stats(all_acts['stainless_steel'].CFs)
    cols = all_acts['stainless_steel'].CFs.columns[2:]
    cf_dct['StainlessSteelSheet'][cols] += get_stats(all_acts['steel_rolling'].CFs)[cols]
    cf_dct['Steel'] = get_stats(all_acts['steel'].CFs)
    cf_dct['Trucking'] = get_stats(all_acts['trucking'].CFs)
    cf_dct['Wood'] = get_stats(all_acts['wood'].CFs)
    return cf_dct

def create_items(ind_df_processed, cf_dct):
    items = []
    for item_ID, df in cf_dct.items():
        item = qs.ImpactItem(ID=item_ID,
                             functional_unit=df.loc[1, ('-', '-', 'functional unit')])

        for num in ind_df_processed.index:
            ind_ID = ind_df_processed.iloc[num]['indicator']
            ind_str = ind_df_processed.iloc[num]['full_name']
            ind_col = tuple(i for i in ind_str.split("'") if len(i)>2)
            item.add_indicator(ind_ID, CF_value=df[df.values=='mean'][ind_col].item())

        items.append(item)
    return items


# %%

# =============================================================================
# Run and save data
# =============================================================================

def get_cf_data():
    apos371, ind_df_processed, indicators = create_indicators()
    all_acts = select_items(apos371)
    cf_dct = organize_cfs(all_acts)

    # Only run this at the very end to remove the outdated setup.pickle file
    remove_setups_pickle()

    return ind_df_processed, all_acts, cf_dct

def save_cf_data():
    ind_df_processed, all_acts, cf_dct = get_cf_data()

    # path = data_path if data_path else os.path.join(os.path.split(c_path)[0], 'data')
    ind_file = 'indicators_new.tsv'
    raw_path = os.path.join(data_path, 'new_CFs_raw')

    # Impact indicators
    ind_df_processed.to_csv(os.path.join(data_path, ind_file), sep='\t')

    # Raw data
    if not os.path.isdir(raw_path):
        os.mkdir(raw_path)
    for k, v in all_acts.items():
        v.CFs.to_csv(os.path.join(raw_path, f'{k}_CFs.tsv'), sep='\t')

    # Organized data
    f = open(os.path.join(data_path, 'cf_dct.pckl'), 'wb')
    pickle.dump(cf_dct, f)
    f.close()


def load_lca_data(kind='original'):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    kind : str
        "original" loads the data from Trimmer et al.
        (TRACI, ecoinvent v3.2),
        "new" loads the data for ReCiPe and TRACI
        (ecoinvent 3.7.1, at the point of substitution).
    '''
    indicator_path = os.path.join(data_path, f'indicators_{kind}.tsv')
    indel_col = None if kind=='original' else 0
    ind_df_processed = pd.read_csv(indicator_path, sep='\t', index_col=indel_col)
    qs.ImpactIndicator.load_indicators_from_file(indicator_path)
    indicators = qs.ImpactIndicator.get_all_indicators()

    if kind == 'original':
        item_path = os.path.join(data_path, 'items_original.xlsx')
        qs.ImpactItem.load_items_from_excel(item_path)
        items = qs.ImpactItem.get_all_items()
    else:
        item_path = os.path.join(data_path, 'cf_dct.pckl')
        f = open(item_path, 'rb')
        cf_dct = pickle.load(f)
        f.close()
        items = create_items(ind_df_processed, cf_dct)

    return indicators, items