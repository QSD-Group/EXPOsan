# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd
from math import ceil
from qsdsan.utils import auom

__all__ = ('categorize_construction_impacts',
           'categorize_OM_impacts',
           'categorize_all_impacts')

#%%

def get_construction_impacts(lca, constructions=None, time=None, time_unit='hr'):
    '''
    Return all construction-related impacts for the given unit,
    normalized to a certain time frame.
    '''
    annualize = lca.annualize_construction
    if time is None:
        time = lca.lifetime_hr
    else:
        time = auom(time_unit).convert(float(time), 'hr')
    impacts = dict.fromkeys((i.ID for i in lca.indicators), 0.)
    if constructions is None: constructions = lca.construction_inventory
    for j in constructions:
        impact = j.impacts
        if j.lifetime is not None: # this equipment has a lifetime
            constr_lifetime = auom('yr').convert(j.lifetime, 'hr')
            ratio = ceil(time/constr_lifetime) if not annualize else time/constr_lifetime
        else: # equipment doesn't have a lifetime
            i = j.linked_unit
            if i is not None and i.lifetime and not isinstance(i.lifetime, dict): # unit has a uniform lifetime
                    constr_lifetime = auom('yr').convert(i.lifetime, 'hr')
                    ratio = ceil(time/constr_lifetime) if not annualize else time/constr_lifetime
            else: # no lifetime, assume just need one
                ratio = 1.
        for m, n in impact.items():
            if m not in impacts.keys():
                continue
            impacts[m] += n*ratio
    return impacts

def categorize_construction_impacts(lca, time=None, time_unit='hr'):
    inv = set(lca.construction_inventory)
    cats = {}
    cats['beads'] = [i for i in inv if '_beads_' in i.ID]
    inv -= set(cats['beads'])
    cats['vessel'] = [i for i in inv if '_R1_' in i.ID or '_R2_' in i.ID]
    inv -= set(cats['vessel'])
    cats['dm'] = [i for i in inv 
                  if ('_DMe_' in i.ID or '_DMs_' in i.ID)
                  and not 'VacPump' in i.ID]
    inv -= set(cats['dm'])
    cats['others'] = list(inv)
    cat_impacts = {
        k: get_construction_impacts(lca, v, time, time_unit) 
        for k,v in cats.items()
        }
    return cat_impacts

def categorize_OM_impacts(lca, time=None, time_unit='hr'):
    if not time:
        time = lca.lifetime_hr
    else:
        time = auom(time_unit).convert(float(time), 'hr')
    factor = time / lca.lifetime_hr
    cat_impacts = {}
    lca.refresh_other_items()
    other_dct = lca.other_items
    for i in ('electricity', 'heat_onsite'):
        itm, quantity = other_dct[i]['item'], other_dct[i]['quantity']
        cat_impacts[i] = {k: factor*quantity*v for k,v in itm.CFs.items()}
    inv = lca.stream_inventory
    cat_impacts['chemicals'] = lca.get_stream_impacts(
        stream_items=[i for i in inv if i.ID.startswith(('DMe_', 'DMs_'))],
        time=time, time_unit=time_unit
        )
    cat_impacts['biogas_offset'] = lca.get_stream_impacts(
        stream_items=[i for i in inv if i.ID.endswith('_NG_offset')],
        time=time, time_unit=time_unit
        )
    cat_impacts['fug_ch4'] = lca.get_stream_impacts(
        stream_items=[i for i in inv if i.ID.endswith('_fugitive_ch4')],
        time=time, time_unit=time_unit
        )
    return cat_impacts

def categorize_all_impacts(lca, time=None, time_unit='hr'):
    const = categorize_construction_impacts(lca, time, time_unit)
    om = categorize_OM_impacts(lca, time, time_unit)
    dct = {**const, **om}
    dct['total'] = lca.get_total_impacts(time=time, time_unit=time_unit)
    out = pd.DataFrame.from_dict(dct)
    return out