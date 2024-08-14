#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    sanunits as su,
    ImpactItem,
    System, TEA, LCA,
    )
import biosteam as bst
from qsdsan.utils import clear_lca_registries
from exposan.utils import (
    add_fugitive_items,
    get_decay_k,
    )
# from exposan.VR_toilet import (
#     _load_components,
#     _load_lca_data,
#     discount_rate,
#     get_toilet_user,
#     max_CH4_emission,
#     operator_daily_wage,
#     ppl,
#     price_dct,
#     )

# __all__= ('create_system',)

#%%
# =========================================================================
#Universal units and functions
# =========

def batch_create_streams(prefix= 'A', phases=('liq', 'sol')):
    #link impact item CH4 & N2O to their streams
    #Q: where do we get impact item from?
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item',set_as_source=True)
    WasteStream('CH4', phase= 'g', stream_impact_item= item)
    item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item',set_as_source=True)
    WasteStream('N2O', phase= 'g', stream_impact_item= item)
    
    for nutrient in ('N','P','K'):
        for phase in phases:
            item = ImpactItem.get_item(f'{nutrient}_item').copy(f'{phase}_{nutrient}_item', 
                                                                set_as_source=True)
            WasteStream(f'{phase}_{nutrient}', phase='l',
                        #price=price_dct[nutrient], 
                        stream_impact_item=item)

    def create_stream_with_impact_item(stream_ID='', item_ID='', dct_key=''):
        item_ID = item_ID or stream_ID+'_item'
        dct_key = dct_key or stream_ID
        item = ImpactItem.get_item(item_ID).copy(f'{prefix}_{item_ID}', set_as_source=True)
        WasteStream(f'{stream_ID}', phase='s',
                    #price=price_dct.get(dct_key) or 0., 
                    stream_impact_item=item)

    # create_stream_with_impact_item(stream_ID='MgOH2')
    # create_stream_with_impact_item(stream_ID='struvite')
    # create_stream_with_impact_item(stream_ID='salt')
    # create_stream_with_impact_item(stream_ID='HCl_acid')

# %%

# =========================================================================
# Unit parameters
# =========================================================================
default_ppl = 5
discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

emptying_fee = 0.15

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05


# %%
# =========================================================================
#Volume reduction toilet based on
#https://patentimages.storage.googleapis.com/73/c1/74/aa1f6a78b89957/WO2023288326A1.pdf
# =========================================================================

ww1 = WasteStream(ID='ww1',H2O= 1500, units = 'L/d')
# Q: how about components?

#### Human Inputs ####
A1 = su.Excretion('A1', outs=('urine','feces'))

A2 = su.MURT('A2',
             ins= (A1-0, A1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
             outs = ('mixed_waste','A2_CH4','A2_N2O'),
             N_user = 5, N_tot_user=default_ppl, lifetime = 10, 
             if_include_front_end=True, if_toilet_paper=True, 
             if_flushing=True, if_cleansing=False,
             if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
             CAPEX=500*max(1, default_ppl/100), OPEX_over_CAPEX=0.06,
             decay_k_COD=get_decay_k(),
             decay_k_N=get_decay_k(),
             max_CH4_emission=max_CH4_emission)

A3 = bst.RotaryVacuumFilter('A3',
                            ins = A2-0,
                            outs = ('A3_solids','A3_liquid'),
                            split = 0.1,
                            moisture_content = 0.99 #[0045]
                            )
A4 = su.BeltThickener('A4', 
                      ins = (A3-0,A3-1,A5-1),#Q: A5 can't be here?
                      outs = ('A4_sludge','A4_liquid'),
                      max_capacity = 0.1,#m3/h, capacity is too small and falls outside the given range
                      power_demand = 0.2, #kW, assumption, need data support
                      sludge_moisture = 0.96
                      )
A5 = su.ReclaimerUltrafiltration('A5',
                                 ins = A4-1,
                                 outs = ('A5_treated','A5_retentate'),
                                 ppl = 5
                                 )
A6 = bst.units.ReverseOsmosis('A6', #no attribute ReverseOsmosis?
                        ins = A5-0,
                        outs = ('A6_treated','A6_brine')
                        )



    
    
    

