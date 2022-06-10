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


# %%

import qsdsan as qs
from qsdsan import ImpactItem, StreamImpactItem

__all__ = ('load_process_settings', )

exchange_rate = 3700 # UGX per USD

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25

# Energetic content of the biogas
# biogas_energy in kJ/mol (as CH4, 16 is the MW of CH4),
# LPG_energy in MJ/kg
get_biogas_factor = lambda biogas_energy=803, LPG_energy=50: biogas_energy/16/LPG_energy

price_dct = {
    'Electricity': 0.17,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*price_factor,
    'P': 3.983*price_factor,
    'K': 1.333*price_factor,
    'Biogas': 6500/exchange_rate*get_biogas_factor()
    }

GWP_dct = {
    'Electricity': 0.1135,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Biogas': -3*get_biogas_factor()
    }

# =============================================================================
# Prices and GWP CFs
# =============================================================================

def batch_create_stream_items(kind):
    if kind == 'original':
        for k, v in GWP_dct.items():
            if k == 'Electricity':
                ImpactItem(ID='E_item', functional_unit='kWh', GWP=v)
            else:
                StreamImpactItem(ID=f'{k}_item', GWP=v)
    elif kind == 'new':
        EcosystemQuality_factor = 29320 # pt/species/yr
        HumanHealth_factor = 436000 # pt/DALY

        E_factor = {
            # Global warming to (terrestrial+freshwater) ecosystem
            'GW2ECO': (2.5e-08+6.82e-13)*EcosystemQuality_factor,
            'GW2HH': 1.25e-05*HumanHealth_factor, # global warming to human health
            'OD2HH': 0.00134*HumanHealth_factor, # stratospheric ozone depletion to human health
                    }
        H_factor = {
            'GW2ECO': (2.8e-09+7.65e-14)*EcosystemQuality_factor,
            'GW2HH': 9.28e-07*HumanHealth_factor,
            'OD2HH': 0.000531*HumanHealth_factor,
            }
        I_factor = {
            'GW2ECO': (5.32e-10+1.45e-14)*EcosystemQuality_factor,
            'GW2HH': 8.12e-08*HumanHealth_factor,
            'OD2HH': 0.000237*HumanHealth_factor,
            }

        StreamImpactItem(ID='CH4_item',
                         E_EcosystemQuality_Total=E_factor['GW2ECO']*4.8,
                         E_HumanHealth_Total=E_factor['GW2HH']*4.8,
                         H_EcosystemQuality_Total=H_factor['GW2ECO']*34,
                         H_HumanHealth_Total=H_factor['GW2HH']*34,
                         I_EcosystemQuality_Total=I_factor['GW2ECO']*84,
                         I_HumanHealth_Total=I_factor['GW2HH']*84
                         )
        StreamImpactItem(ID='N2O_item',
                         E_EcosystemQuality_Total=E_factor['GW2ECO']*78.8,
                         # From climate change + ozone depletion
                         E_HumanHealth_Total=\
                             E_factor['GW2HH']*78.8+E_factor['OD2HH']*0.017,
                         H_EcosystemQuality_Total=H_factor['GW2ECO']*298,
                         H_HumanHealth_Total=\
                             H_factor['GW2HH']*298+H_factor['OD2HH']*0.011,
                         I_EcosystemQuality_Total=I_factor['GW2ECO']*264,
                         I_HumanHealth_Total=\
                             I_factor['GW2HH']*264+I_factor['OD2HH']*0.007
                         )
    else:
        raise ValueError(f'`kind` can only be "original" or "new", not "{kind}".')

    global _ImpactItem_LOADED
    _ImpactItem_LOADED = True


def load_process_settings():
    qs.currency = 'USD'
    qs.CEPCI = qs.CEPCI_by_year[2018]

    qs.PowerUtility.price = price_dct['Electricity']

    from exposan.bwaise._lca_data import lca_data_kind, load_lca_data, _ImpactItem_LOADED

    load_lca_data('original')
    if not _ImpactItem_LOADED: batch_create_stream_items(kind=lca_data_kind)

    ImpactItem.get_item('Concrete').price = price_dct['Concrete']
    ImpactItem.get_item('Steel').price = price_dct['Steel']


    '''
    try: # prevent from reloading
        Excavation = ImpactItem.get_item('Excavation')
        Excavation.indicators[0]
    except:
        load_lca_data('original')
    '''