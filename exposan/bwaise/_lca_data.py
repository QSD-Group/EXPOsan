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

import os, sys, pickle
import pandas as pd
import qsdsan as qs

_ImpactItem_LOADED = False
lca_data_kind = 'original'

c_path = os.path.dirname(__file__)
data_path = os.path.join(c_path, 'data')

__all__ = ('get_cf_data', 'save_cf_data', 'load_lca_data',)


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

def create_indicators(replace=True):
    from bw2qsd import CFgetter
    from bw2qsd.utils import format_name
    apos371 = CFgetter('apos371')
    # ecoinvent version 3.7.1, at the point of substitution
    apos371.load_database('ecoinvent_apos371', print_msg=False)

    # Only include categorical total (not including total of categories)
    apos371.load_indicators(add=True,
                            method='recipe endpoint', method_exclude='LT',
                            category_exclude='total', indicator='total',
                            print_msg=False)
    # For comparison with Trimmer et al.
    apos371.load_indicators(add=True, method='TRACI', indicator='global warming',
                            print_msg=False)

    ind_df_raw = apos371.export_indicators(show=False, path='')
    ind_df_processed = ind_df_raw.copy()

    for num, ind in ind_df_raw.iterrows():
        old_name = ind['indicator']
        if 'TRACI' in ind['method']:
            new_name = old_name
        else:
            kind = ind['method'].split(',A)')[0][-1] # egalitarian, hierarchist, individualist
            new_name = f'{kind}_{format_name(ind["category"])}_Total' # categorical total

            # # Legacy code to include breakdown of all categories and total of categories
            # kind = ind['method'].split(',A)')[0][-1] # egalitarian, hierarchist, individualist
            # if old_name == 'Total':
            #     if ind['category'] == 'total':
            #         new_name = f'{kind}_All' # total of all categories in E/H/I
            #     else:
            #         new_name = f'{kind}_{format_name(ind["category"])}_Total' # categorical total
            #         ind_df_processed.iloc[num]['category'] += ' total'
            # else:
            #     new_name = f'{kind}_{old_name}'

        ind_df_processed.iloc[num]['indicator'] = new_name

    ind_df_processed.sort_values(by=['method', 'category', 'indicator'], inplace=True)

    if replace:
        for ind_ID in ind_df_processed.indicator:
            ind = qs.ImpactIndicator.get_indicator(ind_ID)
            if ind:
                stdout = sys.stdout
                sys.stdout = open(os.devnull, 'w')
                ind.deregister()
                sys.stdout = stdout

    qs.ImpactIndicator.load_from_file(ind_df_processed)
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

    # Catch all printouts (i.e., don't show them in the console)
    stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    # Construction
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

    # Only one result
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

    wood = new_act(database, all_acts, 'wood')
    wood.load_activities('market sawnwood', add=True, filter={'location': 'RoW'}, mask={'product': 'raw'}, limit=None)
    # wood.load_activities('market sawnwood', add=True, filter={'location': 'GLO'}, mask={'product': 'raw'}, limit=None) # None
    to_remove = []
    for act in wood.activities.keys():
        if ('azobe' in act) or ('paraná pine' in act) or ('lath' in act) or ('beam' in act) or ('board' in act):
            to_remove.append(act)
    wood.remove('activity', to_remove)

    # Transportation
    trucking = new_act(database, all_acts, 'trucking')
    trucking.load_activities('market for transport, freight, lorry, all sizes', add=True, filter={'location': 'Row'}, limit=None)
    # trucking.load_activities('market for transport, freight, lorry, all sizes', add=False, filter={'location': 'Row'}, limit=None) # None

    # Stream
    biogas = new_act(database, all_acts, 'biogas')
    biogas.load_activities('market group for liquefied petroleum gas', add=True, filter={'location': 'GLO'}, limit=None)
    biogas.load_activities('market for liquefied petroleum gas', add=True, filter={'location': 'RoW'}, limit=None)

    nitrogen = new_act(database, all_acts, 'nitrogen')
    nitrogen.load_activities('nitrogen fertiliser', add=True, limit=None)
    to_remove = [k for k in nitrogen.activities.keys() if (not 'as N' in k or 'organo' in k)]
    nitrogen.remove('activity', to_remove)

    phosphorus = new_act(database, all_acts, 'phosphorus')
    phosphorus.load_activities('phosphorus fertiliser', add=True, limit=None)
    to_remove = [k for k in phosphorus.activities.keys() if (not 'as P2O5' in k or 'organo' in k)]
    phosphorus.remove('activity', to_remove)

    potassium = new_act(database, all_acts, 'potassium')
    potassium.load_activities('potassium fertiliser', add=True, limit=None)
    to_remove = [k for k in potassium.activities.keys() if (not 'as K2O' in k or 'organo' in k)]
    potassium.remove('activity', to_remove)

    # Others
    electricity = new_act(database, all_acts, 'electricity')
    electricity.load_activities('market electricity', add=True, filter={'location': 'RAF'}, limit=None) # RAF is Africa
    electricity.load_activities('electricity, hydro', add=True, filter={'name': 'hydro'}, mask={'name': 'pumped'}, limit=None)

    # Restore printouts
    sys.stdout = stdout

    return all_acts


def get_stats(df, keep_raw_data=False):
    df2 = pd.DataFrame(columns=df.columns, dtype='float64') if not keep_raw_data else df.copy()
    df2 = df2.append(df[1:].min(), ignore_index=True)
    df2 = df2.append(df[1:].mean(), ignore_index=True)
    # df2 = df2.append(df[1:].median(), ignore_index=True) # use mean instead
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
    cf_dct['Wood'] = get_stats(all_acts['wood'].CFs)
    cf_dct['Trucking'] = get_stats(all_acts['trucking'].CFs)

    cf_dct['Biogas_item'] = get_stats(all_acts['biogas'].CFs)
    cf_dct['Biogas_item'][cols] *= -1 # credit

    cf_dct['N_item'] = get_stats(all_acts['nitrogen'].CFs)
    cf_dct['N_item'][cols] *= -1 # credit

    cf_dct['P_item'] = get_stats(all_acts['phosphorus'].CFs)
    cf_dct['P_item'][cols] /= -(30.973762*2/283.886) # convert from P2O5 to P

    cf_dct['K_item'] = get_stats(all_acts['potassium'].CFs)
    cf_dct['K_item'][cols] /= -(39.098*2/94.2) # convert from K2O to K

    cf_dct['E_item'] = get_stats(all_acts['electricity'].CFs)

    return cf_dct


def create_items(ind_df_processed, cf_dct, replace=True):
    items = []
    for item_ID, df in cf_dct.items():
        item = qs.ImpactItem.get_item(item_ID)
        if not (replace and item):
            if not ('item' in item_ID and item_ID!='E_item'):
                item = qs.ImpactItem(ID=item_ID,
                                     functional_unit=df.loc[1, ('-', '-', 'functional unit')])
            else:
                item = qs.StreamImpactItem(ID=item_ID)

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
    from bw2qsd import remove_setups_pickle
    apos371, ind_df_processed, indicators = create_indicators()
    all_acts = select_items(apos371)
    cf_dct = organize_cfs(all_acts)

    # Only run this at the very end to remove the outdated setup.pickle file
    remove_setups_pickle(False)

    return ind_df_processed, all_acts, cf_dct


def save_cf_data():
    ind_df_processed, all_acts, cf_dct = get_cf_data()

    ind_file = 'indicators_new.tsv'
    raw_path = os.path.join(data_path, 'CFs_new')

    # Impact indicators
    ind_df_processed.to_csv(os.path.join(data_path, ind_file), sep='\t')

    # Raw data
    if not os.path.isdir(raw_path):
        os.mkdir(raw_path)
    for k, v in all_acts.items():
        v.CFs.to_csv(os.path.join(raw_path, f'CFs_{k}.tsv'), sep='\t')

    # Organized data
    f = open(os.path.join(data_path, 'cf_dct.pckl'), 'wb')
    pickle.dump(cf_dct, f)
    f.close()


# %%

# =============================================================================
# Load data
# =============================================================================

def load_lca_data(kind, return_loaded=False):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    kind : str
        "original" loads the data from Trimmer et al.
        (TRACI, ecoinvent v3.2),
        "new" loads the data for ReCiPe and TRACI
        (ecoinvent 3.7.1, at the point of substitution).
    return_loaded : bool
        If True, will return the loaded indicators and items.
    '''
    indicator_path = os.path.join(data_path, f'indicators_{kind}.tsv')
    indel_col = None if kind=='original' else 0
    ind_df_processed = pd.read_csv(indicator_path, sep='\t', index_col=indel_col)
    qs.ImpactIndicator.load_from_file(indicator_path)
    indicators = qs.ImpactIndicator.get_all_indicators()

    global lca_data_kind
    lca_data_kind = kind

    if kind == 'original':
        item_path = os.path.join(data_path, 'items_original.xlsx')
        qs.ImpactItem.load_from_file(item_path)
        items = qs.ImpactItem.get_all_items()
    else:
        item_path = os.path.join(data_path, 'cf_dct.pckl')
        f = open(item_path, 'rb')
        cf_dct = pickle.load(f)
        f.close()
        items = create_items(ind_df_processed, cf_dct)

    global _ImpactItem_LOADED
    _ImpactItem_LOADED = True

    if return_loaded:
        return indicators, items