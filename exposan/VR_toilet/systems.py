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
    Flowsheet, 
    main_flowsheet,
    WasteStream,
    sanunits as qsu,
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
# =============================================================================
# Universal units and functions
# =============================================================================

# !!! cmd 4 & cmd 1

def batch_create_streams(prefix= 'A', phases=('liq', 'sol')):
    #link impact item CH4 & N2O to their streams
    #!!! where do we get impact item from?
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
# TODO: how about components?

#### Human Inputs ####
A1 = qsu.Excretion('A1', outs=('urine','feces'))

recycle_fw = WasteStream('reverse_osmosis_treated')
A2 = qsu.MURT('A2',
             ins= (A1-0, A1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 
                   recycle_fw),
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
# make a recycle loop
recycle_uf = WasteStream('ultrafiltration_reject')
A4 = qsu.BeltThickener('A4', 
                      ins = (A3-0,A3-1,recycle_uf),
                      outs = ('A4_sludge','A4_liquid'),
                      max_capacity = 0.1,#m3/h, based on 5 users. 
                      #capacity is too small and falls outside the given range
                      power_demand = 0.2, #kW, assumption, need data support
                      sludge_moisture = 0.96
                      )

A5 = qsu.ReclaimerUltrafiltration('A5',
                                 ins = A4-1,
                                 outs = ('A5_treated',recycle_uf),
                                 ppl = 5
                                 )
# recycle_ro = WasteStream('reverse_osmosis_brine')
A6 = bst.wastewater.conventional.ReverseOsmosis('A6', 
                        ins = A5-0,
                        outs = (recycle_fw,'A6_brine')
                        )

A7 = qsu.BinaryDistillation('A7', 
                           ins = A6-1,
                           outs = ('A7_distillate','A7_bottoms_product'),
                           LHK = ('Water','NaCl'),#!!!not sure what the heavy key should be
                           #NH3 should be the lighter than water?
                           y_top = 0.99, #mostly water vapor in distillate
                           x_bot = 0.90, #90% water in brine
                           k=0 #!!!not sure what this should be.reflux ratio=0?
                           )
A8 = qsu.BiogenicRefineryGrinder('A8',
                                ins = A4-0,
                                outs = 'A8_grinded'
                                )
A9 = qsu.SludgePasteurization('A9', 
                             ins = A8-0,
                             outs = 'A9_treated',
                             )
A10 = 


    
    
    

