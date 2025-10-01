#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 16:01:43 2025

@author: blues
"""

# !!! check the units in the first tab of the LCI table, also check with Zack
# the values he used.


_impact_item_loaded = False
def _load_lca_data(reload=False):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    reload : bool
        Whether to force reload LCA data.
    '''
    global _impact_item_loaded
    if not _impact_item_loaded or reload:
        indicator_path = os.path.join(g2rt_data_path, 'impact_indicators.csv')
        qs.ImpactIndicator.load_from_file(indicator_path)

        item_path = os.path.join(g2rt_data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)
        
        price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct = update_resource_recovery_settings()

        # Impacts associated with streams and electricity
        def create_stream_impact_item(item_ID, dct_key=''):
            dct_key = dct_key or item_ID.rsplit('_item')[0] # `rstrip` will change "struvite_item" to "struv"
            StreamImpactItem(ID=item_ID, GWP=GWP_dct[dct_key],
                             H_Ecosystems=H_Ecosystems_dct[dct_key],
                             H_Health=H_Health_dct[dct_key],
                             H_Resources=H_Resources_dct[dct_key])

        create_stream_impact_item(item_ID='CH4_item')
        create_stream_impact_item(item_ID='N2O_item')
        create_stream_impact_item(item_ID='N_item')
        create_stream_impact_item(item_ID='P_item')
        create_stream_impact_item(item_ID='K_item')
        create_stream_impact_item(item_ID='H2O_item')
        create_stream_impact_item(item_ID='NO_item')
        create_stream_impact_item(item_ID='SO2_item')
        create_stream_impact_item(item_ID='NH3_item')
        ImpactItem(ID='e_item', functional_unit='kWh',
                   GWP=GWP_dct['Electricity'],
                   H_Ecosystems=H_Ecosystems_dct['Electricity'],
                   H_Health=H_Health_dct['Electricity'],
                   H_Resources=H_Resources_dct['Electricity'])

        _impact_item_loaded = True

    return _impact_item_loaded
