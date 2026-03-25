#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import SanUnit,  WasteStream
from qsdsan.utils import auom
from qsdsan.sanunits import HXutility, SludgePump
from qsdsan.equipments import VerticalMixer
from biosteam import Stream
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch, pressure_vessel_material_factors as factors
from thermosteam.reaction import ParallelReaction
from qsdsan.sanunits import SludgeCentrifuge

# =============================================================================
# Elemental accounting utilities (for Sankey / flow tracing)
# =============================================================================

MW_P = 30.97
MW_PO4 = 94.97
MW_Fe = 55.845
MW_FePO4 = 150.82
MW_FePO4_2H2O = 186.85

def P_mass(stream):
    '''kg P per hour'''
    return (
        stream.imass['PO4'] * MW_P / MW_PO4
        + stream.imass['FePO4'] * MW_P / MW_FePO4
        + stream.imass['FePO4_2H2O'] * MW_P / MW_FePO4_2H2O
    )

def Fe_mass(stream):
    '''kg Fe per hour'''
    return (
        stream.imass['Fe2']
        + stream.imass['Fe3']
        + stream.imass['FePO4'] * MW_Fe / MW_FePO4
        + stream.imass['FePO4_2H2O'] * MW_Fe / MW_FePO4_2H2O
    )

class ElementFlowMixin:
    '''
    Provide elemental flow properties (Fe / P / C) for SanUnit.
    '''

    @property
    def element_in(self):
        '''Total elemental inflow to the unit [kg/hr].'''
        return {
            'Fe': sum(Fe_mass(s) for s in self.ins),
            'P':  sum(P_mass(s)  for s in self.ins)
        }

    @property
    def element_out(self):
        '''Total elemental outflow from the unit (kg/hr).'''
        return {
            'Fe': sum(Fe_mass(s) for s in self.outs),
            'P':  sum(P_mass(s)  for s in self.outs)
        }

    @property
    def element_distribution(self):
        '''
        Elemental distribution to each outlet.
        Returns: list of dicts, one per outlet.
        '''
        total_in = self.element_in
        dist = []

        for s in self.outs:
            dist.append({
                'stream': s.ID,
                'Fe': Fe_mass(s) / total_in['Fe'] if total_in['Fe'] > 0 else 0,
                'P':  P_mass(s)  / total_in['P']  if total_in['P']  > 0 else 0
            })
        return dist

_C_to_K = 273.15
_316_over_304 = factors['Stainless steel 316'] / factors['Stainless steel 304']
_ton_to_tonne = auom('ton').conversion_factor('tonne')

__all__ = (
    'AcidogenicFermenter',
    'SludgeCentrifugeWithElementFlow',
    'SelectivePrecipitation',
    'HeatDrying',
    'Sintering',
    'FePO4_recovery'
)

# =============================================================================
# AcidogenicFermenter
# =============================================================================

@cost('Recirculation flow rate', 'Recirculation pumps', kW=30, S=77.22216,
      cost=47200, n=0.8, BM=2.3, CE=522, N='Number of reactors')
@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3*_316_over_304, S=3785, n=0.6, BM=1.8)
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500*_316_over_304,
      S=3785, n=0.5, kW=22.371, BM=1.5, N='Number of reactors')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000*_316_over_304,
      S=3785, n=0.5, BM=1.5, N='Number of reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900*_316_over_304,
      S=20920000.0, n=0.7, BM=2.2, N='Number of reactors',
      magnitude=True)
class AcidogenicFermenter(ElementFlowMixin, SanUnit):
    '''
    Fermentation of sludge and food waste.
    Assume only pumps already use 316 SS; others need to be converted
    from 304 SS to 316 SS due to acidic conditions.
    See the reference (https://docs.nrel.gov/docs/fy11osti/47764.pdf) cited in
    :class:`biosteam.units.NRELBatchBioreactor` for more information.
    
    Parameters
    ----------
    ins : iterable
        sludge, food_waste.
    outs : iterable
        fermentate, gas.
    food_sludge_ratio : float
        Mass ratio of organics bewteen food waste and sludge, [-].
    food_waste_moisture : float
        Mositure content of food waste, [-].
    org_to_gas : float
        Mass ratio of total organics converted to gas, [-].
    org_to_vfa : float
        Mass ratio of total organics converted to volatile fatty acids (VFAs), [-].
    org_to_ethanol : float
        Mass ratio of total organics converted to ethanol, [-].
    org_to_residue : float
        Mass ratio of total organics in the residue after acidogenic fermentation, [-].
    VFA_ratio : float
        Mass fractions of individual VFAs, [-].
    Fe_reduction : float
        Fe3+ reduction ratio during acidogenic fermentation, [-].
    metal_release_ratio_X : float
        metal release ratio at different food_sludge_ratio (X), [-].
    P_release_ratio_X : float
        phosphorus release ratio at different food_sludge_ratio (X), [-].
    T : float
        Required temperature for fermentation, [K].
    
    See Also
    --------
    :class:`biosteam.units.NRELBatchBioreactor`
    '''
    _N_ins = 2
    _N_outs = 2
    
    auxiliary_unit_names=(
        'fermentate_pump',
    )
    
    _units = {
        'Reactor volume': 'm3',
        'Cycle time': 'hr',
        'Batch time': 'hr',
        'Loading time': 'hr',
        'Total dead time': 'hr',
        'Reactor duty': 'kJ/hr',
        'Recirculation flow rate': 'm3/hr'
    }
    
    # hydraulic retention time, [hr]
    # HRT = 5*24
    
    # cleaning and unloading, [hr]
    tau_cleaning = 3
    
    # number of fermenter, [-]
    N = 2
    
    # fraction of filled tank to total tank volume, [-]
    V_wf = 0.9
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 food_sludge_ratio=1,
                 HRT=132,
                 food_waste_moisture=0.74,
                 org_to_gas=0.05,
                 org_to_vfa=0.65,
                 org_to_ethanol=0.02,
                 org_to_residue=0.25,
                 VFA_ratio={'Ac': 0.5, 'Pr': 0.24, 'Bu': 0.23, 'Va': 0.02, 'Lac': 0.01},
                 Fe_reduction=0.98,
                 metal_release_ratio_0_0 = 0.00739,
                 metal_release_ratio_0_12 = 0.007535,
                 metal_release_ratio_0_24 = 0.006234,
                 metal_release_ratio_0_36 = 0.004403,
                 metal_release_ratio_0_48 = 0.0058,
                 metal_release_ratio_0_60 = 0.004596,
                 metal_release_ratio_0_72 = 0.006041,
                 metal_release_ratio_0_84 = 0.0077053,
                 metal_release_ratio_0_96 = 0.006668,
                 metal_release_ratio_0_108 = 0.004837,
                 metal_release_ratio_0_120 = 0.00397,              
                 metal_release_ratio_0_132 = 0.02124263,
                 metal_release_ratio_1_3_0 = 0.00739,
                 metal_release_ratio_1_3_12 = 0.105381,
                 metal_release_ratio_1_3_24 = 0.140121,
                 metal_release_ratio_1_3_36 = 0.16433,
                 metal_release_ratio_1_3_48 = 0.175555,
                 metal_release_ratio_1_3_60 = 0.180373,
                 metal_release_ratio_1_3_72 = 0.19049,
                 metal_release_ratio_1_3_84 = 0.186636,
                 metal_release_ratio_1_3_96 = 0.196753,
                 metal_release_ratio_1_3_108 = 0.173146,
                 metal_release_ratio_1_3_120 = 0.148576,
                 metal_release_ratio_1_3_132 = 0.2113664,
                 metal_release_ratio_2_3_0 = 0.00739,
                 metal_release_ratio_2_3_12 = 0.149463,
                 metal_release_ratio_2_3_24 = 0.290793,
                 metal_release_ratio_2_3_36 = 0.361517,
                 metal_release_ratio_2_3_48 = 0.40054,
                 metal_release_ratio_2_3_60 = 0.427037,
                 metal_release_ratio_2_3_72 = 0.443898,
                 metal_release_ratio_2_3_84 = 0.475695,
                 metal_release_ratio_2_3_96 = 0.488221,
                 metal_release_ratio_2_3_108 = 0.5099,
                 metal_release_ratio_2_3_120 = 0.514718,
                 metal_release_ratio_2_3_132 = 0.6067351,
                 metal_release_ratio_1_0 = 0.00739,
                 metal_release_ratio_1_12 = 0.175141,
                 metal_release_ratio_1_24 = 0.357759,
                 metal_release_ratio_1_36 = 0.454979,
                 metal_release_ratio_1_48 = 0.48533,
                 metal_release_ratio_1_60 = 0.516645,
                 metal_release_ratio_1_72 = 0.577829,
                 metal_release_ratio_1_84 = 0.627451,
                 metal_release_ratio_1_96 = 0.671292,
                 metal_release_ratio_1_108 = 0.728622,
                 metal_release_ratio_1_120 = 0.772462,
                 metal_release_ratio_1_132 = 0.8307559,
                 metal_release_ratio_4_3_132 = 0.85,
                 P_release_ratio_0_0 = 0.08391,
                 P_release_ratio_0_12 = 0.138184,
                 P_release_ratio_0_24 = 0.147972,
                 P_release_ratio_0_36 = 0.160052,
                 P_release_ratio_0_48 = 0.162135,
                 P_release_ratio_0_60 = 0.093407,
                 P_release_ratio_0_72 = 0.099447,
                 P_release_ratio_0_84 = 0.10809,
                 P_release_ratio_0_96 = 0.112151,
                 P_release_ratio_0_108 = 0.113921,
                 P_release_ratio_0_120 = 0.110693,
                 P_release_ratio_0_132 = 0.110693,
                 P_release_ratio_1_3_0 = 0.083931,
                 P_release_ratio_1_3_12 = 0.248252,
                 P_release_ratio_1_3_24 = 0.412365,
                 P_release_ratio_1_3_36 = 0.460266,
                 P_release_ratio_1_3_48 = 0.502336,
                 P_release_ratio_1_3_60 = 0.494005,
                 P_release_ratio_1_3_72 = 0.516498,
                 P_release_ratio_1_3_84 = 0.514832,
                 P_release_ratio_1_3_96 = 0.527744,
                 P_release_ratio_1_3_108 = 0.522746,
                 P_release_ratio_1_3_120 = 0.512749,
                 P_release_ratio_1_3_132 = 0.4431886,
                 P_release_ratio_2_3_0 = 0.083931,
                 P_release_ratio_2_3_12 = 0.284178,
                 P_release_ratio_2_3_24 = 0.544405,
                 P_release_ratio_2_3_36 = 0.634793,
                 P_release_ratio_2_3_48 = 0.612717,
                 P_release_ratio_2_3_60 = 0.62313,
                 P_release_ratio_2_3_72 = 0.651454,
                 P_release_ratio_2_3_84 = 0.650204,
                 P_release_ratio_2_3_96 = 0.667699,
                 P_release_ratio_2_3_108 = 0.655619,
                 P_release_ratio_2_3_120 = 0.69394,
                 P_release_ratio_2_3_132 = 0.7122673,
                 P_release_ratio_1_0 = 0.083931,
                 P_release_ratio_1_12 = 0.299277,
                 P_release_ratio_1_24 = 0.576062,
                 P_release_ratio_1_36 = 0.628128,
                 P_release_ratio_1_48 = 0.655619,
                 P_release_ratio_1_60 = 0.67228,
                 P_release_ratio_1_72 = 0.706019,
                 P_release_ratio_1_84 = 0.700188,
                 P_release_ratio_1_96 = 0.736426,
                 P_release_ratio_1_108 = 0.741424,
                 P_release_ratio_1_120 = 0.79724,
                 P_release_ratio_1_132 = 0.8230645,
                 P_release_ratio_4_3_132 = 0.83,
                 T=37+_C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.food_sludge_ratio = food_sludge_ratio
        self.food_waste_moisture = food_waste_moisture
        self.org_to_gas = org_to_gas
        self.org_to_vfa = org_to_vfa
        self.org_to_ethanol = org_to_ethanol
        self.org_to_residue = org_to_residue
        self.VFA_ratio = VFA_ratio
        self.HRT = HRT
        self.Fe_reduction = Fe_reduction
        self.metal_release_ratio_0_0 = metal_release_ratio_0_0
        self.metal_release_ratio_0_12 = metal_release_ratio_0_12
        self.metal_release_ratio_0_24 = metal_release_ratio_0_24
        self.metal_release_ratio_0_36 = metal_release_ratio_0_36
        self.metal_release_ratio_0_48 = metal_release_ratio_0_48
        self.metal_release_ratio_0_60 = metal_release_ratio_0_60
        self.metal_release_ratio_0_72 = metal_release_ratio_0_72
        self.metal_release_ratio_0_84 = metal_release_ratio_0_84
        self.metal_release_ratio_0_96 = metal_release_ratio_0_96
        self.metal_release_ratio_0_108 = metal_release_ratio_0_108
        self.metal_release_ratio_0_120 = metal_release_ratio_0_120
        self.metal_release_ratio_0_132 = metal_release_ratio_0_132
        self.metal_release_ratio_1_3_0 = metal_release_ratio_1_3_0
        self.metal_release_ratio_1_3_12 = metal_release_ratio_1_3_12
        self.metal_release_ratio_1_3_24 = metal_release_ratio_1_3_24
        self.metal_release_ratio_1_3_36 = metal_release_ratio_1_3_36
        self.metal_release_ratio_1_3_48 = metal_release_ratio_1_3_48
        self.metal_release_ratio_1_3_60 = metal_release_ratio_1_3_60
        self.metal_release_ratio_1_3_72 = metal_release_ratio_1_3_72
        self.metal_release_ratio_1_3_84 = metal_release_ratio_1_3_84
        self.metal_release_ratio_1_3_96 = metal_release_ratio_1_3_96
        self.metal_release_ratio_1_3_108 = metal_release_ratio_1_3_108
        self.metal_release_ratio_1_3_120 = metal_release_ratio_1_3_120
        self.metal_release_ratio_1_3_132 = metal_release_ratio_1_3_132
        self.metal_release_ratio_2_3_0 = metal_release_ratio_2_3_0
        self.metal_release_ratio_2_3_12 = metal_release_ratio_2_3_12
        self.metal_release_ratio_2_3_24 = metal_release_ratio_2_3_24
        self.metal_release_ratio_2_3_36 = metal_release_ratio_2_3_36
        self.metal_release_ratio_2_3_48 = metal_release_ratio_2_3_48
        self.metal_release_ratio_2_3_60 = metal_release_ratio_2_3_60
        self.metal_release_ratio_2_3_72 = metal_release_ratio_2_3_72
        self.metal_release_ratio_2_3_84 = metal_release_ratio_2_3_84
        self.metal_release_ratio_2_3_96 = metal_release_ratio_2_3_96
        self.metal_release_ratio_2_3_108 = metal_release_ratio_2_3_108
        self.metal_release_ratio_2_3_120 = metal_release_ratio_2_3_120
        self.metal_release_ratio_2_3_132 = metal_release_ratio_2_3_132
        self.metal_release_ratio_1_0 = metal_release_ratio_1_0
        self.metal_release_ratio_1_12 = metal_release_ratio_1_12
        self.metal_release_ratio_1_24 = metal_release_ratio_1_24
        self.metal_release_ratio_1_36 = metal_release_ratio_1_36
        self.metal_release_ratio_1_48 = metal_release_ratio_1_48
        self.metal_release_ratio_1_60 = metal_release_ratio_1_60
        self.metal_release_ratio_1_72 = metal_release_ratio_1_72
        self.metal_release_ratio_1_84 = metal_release_ratio_1_84
        self.metal_release_ratio_1_96 = metal_release_ratio_1_96
        self.metal_release_ratio_1_108 = metal_release_ratio_1_108
        self.metal_release_ratio_1_120 = metal_release_ratio_1_120
        self.metal_release_ratio_1_132 = metal_release_ratio_1_132
        self.metal_release_ratio_4_3_132 = metal_release_ratio_4_3_132
        self.P_release_ratio_0_0 = P_release_ratio_0_0
        self.P_release_ratio_0_12 = P_release_ratio_0_12
        self.P_release_ratio_0_24 = P_release_ratio_0_24
        self.P_release_ratio_0_36 = P_release_ratio_0_36
        self.P_release_ratio_0_48 = P_release_ratio_0_48
        self.P_release_ratio_0_60 = P_release_ratio_0_60
        self.P_release_ratio_0_72 = P_release_ratio_0_72
        self.P_release_ratio_0_84 = P_release_ratio_0_84
        self.P_release_ratio_0_96 = P_release_ratio_0_96
        self.P_release_ratio_0_108 = P_release_ratio_0_108
        self.P_release_ratio_0_120 = P_release_ratio_0_120
        self.P_release_ratio_0_132 = P_release_ratio_0_132
        self.P_release_ratio_1_3_0 = P_release_ratio_1_3_0
        self.P_release_ratio_1_3_12 = P_release_ratio_1_3_12
        self.P_release_ratio_1_3_24 = P_release_ratio_1_3_24
        self.P_release_ratio_1_3_36 = P_release_ratio_1_3_36
        self.P_release_ratio_1_3_48 = P_release_ratio_1_3_48
        self.P_release_ratio_1_3_60 = P_release_ratio_1_3_60
        self.P_release_ratio_1_3_72 = P_release_ratio_1_3_72
        self.P_release_ratio_1_3_84 = P_release_ratio_1_3_84
        self.P_release_ratio_1_3_96 = P_release_ratio_1_3_96
        self.P_release_ratio_1_3_108 = P_release_ratio_1_3_108
        self.P_release_ratio_1_3_120 = P_release_ratio_1_3_120
        self.P_release_ratio_1_3_132 = P_release_ratio_1_3_132
        self.P_release_ratio_2_3_0 = P_release_ratio_2_3_0
        self.P_release_ratio_2_3_12 = P_release_ratio_2_3_12
        self.P_release_ratio_2_3_24 = P_release_ratio_2_3_24
        self.P_release_ratio_2_3_36 = P_release_ratio_2_3_36
        self.P_release_ratio_2_3_48 = P_release_ratio_2_3_48
        self.P_release_ratio_2_3_60 = P_release_ratio_2_3_60
        self.P_release_ratio_2_3_72 = P_release_ratio_2_3_72
        self.P_release_ratio_2_3_84 = P_release_ratio_2_3_84
        self.P_release_ratio_2_3_96 = P_release_ratio_2_3_96
        self.P_release_ratio_2_3_108 = P_release_ratio_2_3_108
        self.P_release_ratio_2_3_120 = P_release_ratio_2_3_120
        self.P_release_ratio_2_3_132 = P_release_ratio_2_3_132
        self.P_release_ratio_1_0 = P_release_ratio_1_0
        self.P_release_ratio_1_12 = P_release_ratio_1_12
        self.P_release_ratio_1_24 = P_release_ratio_1_24
        self.P_release_ratio_1_36 = P_release_ratio_1_36
        self.P_release_ratio_1_48 = P_release_ratio_1_48
        self.P_release_ratio_1_60 = P_release_ratio_1_60
        self.P_release_ratio_1_72 = P_release_ratio_1_72
        self.P_release_ratio_1_84 = P_release_ratio_1_84
        self.P_release_ratio_1_96 = P_release_ratio_1_96
        self.P_release_ratio_1_108 = P_release_ratio_1_108
        self.P_release_ratio_1_120 = P_release_ratio_1_120
        self.P_release_ratio_1_132 = P_release_ratio_1_132
        self.P_release_ratio_4_3_132 = P_release_ratio_4_3_132
        self.T = T
        fermentate = self.outs[0].proxy(f'{ID}_fermentate')
        self.fermentate_pump = SludgePump(f'.{ID}_fermentate_pump', ins=fermentate, init_with=init_with)
    
    def _run(self):
        sludge, food_waste = self.ins
        fermentate, gas = self.outs
        
        sludge.phase = 'l'
        food_waste.phase = 'l'
        fermentate.phase = 'l'
        gas.phase = 'g'
        
        if self.food_sludge_ratio not in [0, 1/3, 2/3, 1, 4/3]:
            raise RuntimeError('food_sludge_ratio must be one of the follow: 0, 1/3, 2/3, 1, 4/3.')
        
        food_waste.imass['Org'] = sludge.imass['Org']*self.food_sludge_ratio
        food_waste.imass['H2O'] = food_waste.imass['Org']/(1-self.food_waste_moisture) * self.food_waste_moisture
        
        fermentate.mix_from(self.ins)
        
        org_total = fermentate.imass['Org']
        
        # gas production
        gas.imass['CO2'] = org_total*self.org_to_gas
        
        # VFAs and inorganics in fermentate
        vfa_mass = org_total * self.org_to_vfa
        
        fermentate.imass['Acetic_acid'] = vfa_mass * self.VFA_ratio['Ac']
        fermentate.imass['Propionic_acid'] = vfa_mass * self.VFA_ratio['Pr']
        fermentate.imass['Butyric_acid'] = vfa_mass * self.VFA_ratio['Bu']
        fermentate.imass['Valeric_acid'] = vfa_mass * self.VFA_ratio['Va']
        fermentate.imass['Lactic_acid'] = vfa_mass * self.VFA_ratio['Lac']
        
        fermentate.imass['Ethanol'] = org_total*self.org_to_ethanol
        
        fermentate.imass['Residue'] = org_total*self.org_to_residue
        
        if self.org_to_gas + self.org_to_vfa + self.org_to_ethanol + self.org_to_residue > 1:
            raise RuntimeError('org cannot be balanced')
        
        fermentate.imass['Org'] *= (1 - self.org_to_gas - self.org_to_vfa - self.org_to_ethanol - self.org_to_residue)

        fermentate.imass['Fe2'] = fermentate.imass['Fe3'] * self.Fe_reduction
        fermentate.imass['Fe3'] -= fermentate.imass['Fe2']
        
        # metal_P_release = {
        #     0: {'metal': self.metal_release_ratio_0, 'P': self.P_release_ratio_0},
        #     1/3: {'metal': self.metal_release_ratio_1_3, 'P': self.P_release_ratio_1_3},
        #     2/3: {'metal': self.metal_release_ratio_2_3, 'P': self.P_release_ratio_2_3},
        #     1: {'metal': self.metal_release_ratio_1, 'P': self.P_release_ratio_1},
        #     4/3: {'metal': self.metal_release_ratio_4_3, 'P': self.P_release_ratio_4_3}
        # }
        
        metal_P_release = {
            0: {
                0:{'metal': self.metal_release_ratio_0_0, 'P': self.P_release_ratio_0_0},
                12:{'metal': self.metal_release_ratio_0_12, 'P': self.P_release_ratio_0_12},
                24:{'metal': self.metal_release_ratio_0_24, 'P': self.P_release_ratio_0_24},
                36:{'metal': self.metal_release_ratio_0_36, 'P': self.P_release_ratio_0_36},
                48:{'metal': self.metal_release_ratio_0_48, 'P': self.P_release_ratio_0_48},
                60:{'metal': self.metal_release_ratio_0_60, 'P': self.P_release_ratio_0_60},
                72:{'metal': self.metal_release_ratio_0_72, 'P': self.P_release_ratio_0_72},
                84:{'metal': self.metal_release_ratio_0_84, 'P': self.P_release_ratio_0_84},
                96:{'metal': self.metal_release_ratio_0_96, 'P': self.P_release_ratio_0_96},
                108:{'metal': self.metal_release_ratio_0_108, 'P': self.P_release_ratio_0_108},           
                120: {'metal': self.metal_release_ratio_0_120, 'P': self.P_release_ratio_0_120},
                132:{'metal': self.metal_release_ratio_0_132, 'P': self.P_release_ratio_0_132},
                },
            1/3: {
                0:{'metal': self.metal_release_ratio_1_3_0, 'P': self.P_release_ratio_1_3_0},
                12:{'metal': self.metal_release_ratio_1_3_12, 'P': self.P_release_ratio_1_3_12},
                24:{'metal': self.metal_release_ratio_1_3_24, 'P': self.P_release_ratio_1_3_24},
                36:{'metal': self.metal_release_ratio_1_3_36, 'P': self.P_release_ratio_1_3_36},
                48:{'metal': self.metal_release_ratio_1_3_48, 'P': self.P_release_ratio_1_3_48},
                60:{'metal': self.metal_release_ratio_1_3_60, 'P': self.P_release_ratio_1_3_60},
                72:{'metal': self.metal_release_ratio_1_3_72, 'P': self.P_release_ratio_1_3_72},
                84:{'metal': self.metal_release_ratio_1_3_84, 'P': self.P_release_ratio_1_3_84},
                96:{'metal': self.metal_release_ratio_1_3_96, 'P': self.P_release_ratio_1_3_96},
                108:{'metal': self.metal_release_ratio_1_3_108, 'P': self.P_release_ratio_1_3_108},
                120:{'metal': self.metal_release_ratio_1_3_120, 'P': self.P_release_ratio_1_3_120},
                132:{'metal': self.metal_release_ratio_1_3_132, 'P': self.P_release_ratio_1_3_132},
                },
            2/3: {
                0:{'metal': self.metal_release_ratio_2_3_0, 'P': self.P_release_ratio_2_3_0},
                12:{'metal': self.metal_release_ratio_2_3_12, 'P': self.P_release_ratio_2_3_12},
                24:{'metal': self.metal_release_ratio_2_3_24, 'P': self.P_release_ratio_2_3_24},
                36:{'metal': self.metal_release_ratio_2_3_36, 'P': self.P_release_ratio_2_3_36},
                48:{'metal': self.metal_release_ratio_2_3_48, 'P': self.P_release_ratio_2_3_48},
                60:{'metal': self.metal_release_ratio_2_3_60, 'P': self.P_release_ratio_2_3_60},
                72:{'metal': self.metal_release_ratio_2_3_72, 'P': self.P_release_ratio_2_3_72},
                84:{'metal': self.metal_release_ratio_2_3_84, 'P': self.P_release_ratio_2_3_84},
                96:{'metal': self.metal_release_ratio_2_3_96, 'P': self.P_release_ratio_2_3_96},
                108:{'metal': self.metal_release_ratio_2_3_108, 'P': self.P_release_ratio_2_3_108},
                120:{'metal': self.metal_release_ratio_2_3_120, 'P': self.P_release_ratio_2_3_120},
                132:{'metal': self.metal_release_ratio_2_3_132, 'P': self.P_release_ratio_2_3_132},
                },
            1: {
                0:{'metal': self.metal_release_ratio_1_0, 'P': self.P_release_ratio_1_0},
                12:{'metal': self.metal_release_ratio_1_12, 'P': self.P_release_ratio_1_12},
                24:{'metal': self.metal_release_ratio_1_24, 'P': self.P_release_ratio_1_24},
                36:{'metal': self.metal_release_ratio_1_36, 'P': self.P_release_ratio_1_36},
                48:{'metal': self.metal_release_ratio_1_48, 'P': self.P_release_ratio_1_48},
                60:{'metal': self.metal_release_ratio_1_60, 'P': self.P_release_ratio_1_60},
                72:{'metal': self.metal_release_ratio_1_72, 'P': self.P_release_ratio_1_72},
                84:{'metal': self.metal_release_ratio_1_84, 'P': self.P_release_ratio_1_84},
                96:{'metal': self.metal_release_ratio_1_96, 'P': self.P_release_ratio_1_96},
                108:{'metal': self.metal_release_ratio_1_108, 'P': self.P_release_ratio_1_108},           
                120: {'metal': self.metal_release_ratio_1_120, 'P': self.P_release_ratio_1_120},
                132:{'metal': self.metal_release_ratio_1_132, 'P': self.P_release_ratio_1_132},
                },
            4/3: {'metal': self.metal_release_ratio_4_3_132, 'P': self.P_release_ratio_4_3_132}
        }
        
        metal_P_to_residue = 0
        
        for metal in ['Fe2','Fe3','Ca2','Mg2']:
            metal_P_to_residue += fermentate.imass[metal]*(1 - metal_P_release[self.food_sludge_ratio][self.HRT]['metal'])
            fermentate.imass[metal] *= metal_P_release[self.food_sludge_ratio][self.HRT]['metal']
        
        metal_P_to_residue += fermentate.imass['PO4']*(1 - metal_P_release[self.food_sludge_ratio][self.HRT]['P'])
        fermentate.imass['PO4'] *= metal_P_release[self.food_sludge_ratio][self.HRT]['P']
        
        fermentate.imass['Residue'] += metal_P_to_residue
        
        fermentate.T = gas.T = self.T
    
    def _design(self):
        sludge, food_waste = self.ins
        
        F_vol = sludge.F_vol + food_waste.F_vol
        
        D = self.design_results
        
        D.update(size_batch(F_vol, self.HRT, self.tau_cleaning, self.N, self.V_wf))
        D['Number of reactors'] = self.N
        D['Recirculation flow rate'] =  F_vol/self.N
        
        duty = self.Hnet
        
        D['Reactor duty'] = duty
        
        self.add_heat_utility(duty, self.T)
        
        self.fermentate_pump.simulate()

# =============================================================================
# SludgeCentrifugeWithElementFlow
# =============================================================================

class SludgeCentrifugeWithElementFlow(ElementFlowMixin, SludgeCentrifuge):
    '''SludgeCentrifuge with elemental flow @property.'''
    pass

# =============================================================================
# SelectivePrecipitation
# =============================================================================

@cost('Recirculation flow rate', 'Recirculation pumps', kW=30, S=77.22216,
      cost=47200, n=0.8, BM=2.3, CE=522, N='Number of reactors')
@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3*_316_over_304, S=3785, n=0.6, BM=1.8)
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500*_316_over_304,
      S=3785, n=0.5, kW=22.371, BM=1.5, N='Number of reactors')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000*_316_over_304,
      S=3785, n=0.5, BM=1.5, N='Number of reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900*_316_over_304,
      S=20920000.0, n=0.7, BM=2.2, N='Number of reactors',
      magnitude=True)
class SelectivePrecipitation(ElementFlowMixin, SanUnit):
    '''
    Selective precipitation of FePO4 by adding acid and oxidant.
    Assume only pumps already use 316 SS; others need to be converted
    from 304 SS to 316 SS due to acidic and oxidative conditions.
    See the reference (https://docs.nrel.gov/docs/fy11osti/47764.pdf) cited in
    :class:`biosteam.units.NRELBatchBioreactor` for more information.
    
    Parameters
    ----------
    ins : iterable
        supernatant, acid, oxidant.
    outs : iterable
        slurry.
    acid_dose : float
        Mass ratio of acid (measured as pure H2SO4) to the supernatant after
        acidogenic fermentation to adjust pH to 2, [-].
    oxidant_excess : float
        Excess oxidant ratio, [-].
    T : float
        Required temperature for FePO4 precipitation, [K].
    
    See Also
    --------
    :class:`biosteam.units.NRELBatchBioreactor
    '''
    _N_ins = 3
    _N_outs = 1
    
    auxiliary_unit_names=(
        'acid_pump',
        'oxidant_pump',
        'slurry_pump'
    )
    
    _units = {
        'Reactor volume': 'm3',
        'Cycle time': 'hr',
        'Batch time': 'hr',
        'Loading time': 'hr',
        'Total dead time': 'hr',
        'Reactor duty': 'kJ/hr',
        'Recirculation flow rate': 'm3/hr'
    }
    
    # hydraulic retention time, [hr]
    HRT = 6
    
    # cleaning and unloading, [hr]
    tau_cleaning = 1
    
    # number of precipitation tank, [-]
    N = 2
    
    # fraction of filled tank to total tank volume, [-]
    V_wf = 0.9
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                acid_dose=0.003,
                oxidant_excess=2,
                T=40+_C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_dose = acid_dose
        self.oxidant_excess = oxidant_excess
        self.T = T
        hx_in = Stream(ID=f'{ID}_hx_in')
        hx_out = Stream(ID=f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        self.verticle_mixer = VerticalMixer(ID='verticle_mixer', linked_unit=self)
        self.equipments = (self.verticle_mixer,)
        acid = self.ins[1].proxy(f'{ID}_acid')
        self.acid_pump = SludgePump(f'.{ID}_acid_pump', ins=acid, init_with=init_with)
        oxidant = self.ins[2].proxy(f'{ID}_oxidant')
        self.oxidant_pump = SludgePump(f'.{ID}_oxidant_pump', ins=oxidant, init_with=init_with)
        slurry = self.outs[0].proxy(f'{ID}_slurry')
        self.slurry_pump = SludgePump(f'.{ID}_slurry_pump', ins=slurry, init_with=init_with)
    
    def _run(self):
        supernatant, acid, oxidant = self.ins
        slurry = self.outs[0]
        
        supernatant.phase = 'l'
        acid.phase = 'l'
        oxidant.phase = 'l'
        slurry.phase = 'l'
        
        # acid dosing while adjusting the pH
        acid.imass['H2SO4'] = supernatant.F_mass * self.acid_dose
        # assume the volumetic ratio between H2SO4 and H2O is 1:1
        acid.ivol['H2O'] = acid.ivol['H2SO4']
        
        # H2O2 comsumption (H2O2 + 2Fe2+ -> 2Fe3+ + 2OH-)
        # double H2O2 to ensure complete oxidation
        oxidant.imass['H2O2'] = supernatant.imass['Fe2'] / 56 / 2 * 34 * self.oxidant_excess
        # assume the volumetic ratio between H2O2 and H2O is 3:7
        oxidant.ivol['H2O'] = oxidant.ivol['H2O2'] / 3 * 7
        
        slurry.mix_from(self.ins)
        
        # Fe2+ oxidation to Fe3+
        slurry.imass['Fe3'] += slurry.imass['Fe2']
        slurry.imass['Fe2'] = 0
        
        # FePO4 precipitation
        P_mol = slurry.imass['PO4'] / 95
        
        Fe_mol = slurry.imass['Fe3'] / 56
        Fe_reacted_mol = min(Fe_mol, P_mol)
        
        slurry.imass['PO4'] -= Fe_reacted_mol * 95
        slurry.imass['Fe3'] -= Fe_reacted_mol * 56
        slurry.imass['FePO4_2H2O'] = Fe_reacted_mol * 187
        
        slurry.imass['H2O2'] = 0
        slurry.imass['H2O'] = 0
        slurry.imass['H2O'] = supernatant.F_mass + acid.F_mass + oxidant.F_mass - slurry.F_mass
        
        slurry.T = self.T
    
    def _design(self):
        supernatant, acid, oxidant = self.ins
        
        F_vol = supernatant.F_vol + acid.F_vol + oxidant.F_vol
        
        D = self.design_results
        
        D.update(size_batch(F_vol, self.HRT, self.tau_cleaning, self.N, self.V_wf))
        D['Number of reactors'] = self.N
        D['Recirculation flow rate'] =  F_vol/self.N
        
        duty = self.Hnet
        
        D['Reactor duty'] = duty
        
        self.add_heat_utility(duty, self.T)
        
        self.acid_pump.simulate()
        self.oxidant_pump.simulate()
        self.slurry_pump.simulate()

# =============================================================================
# HeatDrying
# =============================================================================

# n: 0.7, assumed
# BM: 3.17, from biosteam/units/heat_exchange.py
@cost(ID='Dryer 1', basis='Half dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
@cost(ID='Dryer 2', basis='Half dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
class HeatDrying(ElementFlowMixin, SanUnit):
    '''
    Heat drying of sludge or biosolids.
    
    from [1]: 220000 2018$ CAPEX/dry tonne solids/day for a 80 dry tonne/day.
    # assume CAPEX / installed cost = 2 and installed cost / purchase cost (BM) = 3.17
    
    Parameters
    ----------
    ins : iterable
        feed, natural_gas.
    outs : iterable
        dried_solids, vapor.
    T : float
        Heat drying temperature, [K].
    unit_heat : float
        Energy for removing unit water from solids, [GJ/tonne water].
    natural_gas_HHV : float
        Higher heating value of natural gas, [MJ/m3].
    unit_electricity : float
        Electricity for heat drying, [kWh/dry tonne solids].
    
    See Also
    --------
    :class:`exposan.htl.landscape_sanunits.HeatDrying` under `pfas` branch
    
    References
    ----------
    [1] Hao, X.; Chen, Q.; van Loosdrecht, M. C. M.; Li, J.; Jiang, H.
        Sustainable Disposal of Excess Sludge: Incineration without
        Anaerobic Digestion. Water Research 2020, 170, 115298.
        https://doi.org/10.1016/j.watres.2019.115298.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _units = {'Half dry solids flow': 'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', T=105 + _C_to_K, unit_heat=4.5,
                 natural_gas_HHV=39, unit_electricity=214):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.T = T
        self.unit_heat = unit_heat
        self.unit_electricity = unit_electricity
        self.natural_gas_HHV = natural_gas_HHV
    
    def _run(self):
        feed, natural_gas = self.ins
        dried_solids, vapor = self.outs
        
        feed.phase = 'l'
        natural_gas.phase = 'g'
        dried_solids.phase = 's'
        vapor.phase = 'g'
        
        dried_solids.copy_like(feed)    
            
        for cmp in ['H2O', 'Acetic_acid', 'Propionic_acid', 'Butyric_acid',
                    'Valeric_acid', 'Lactic_acid', 'Ethanol']:
            vapor.imass[cmp] = feed.imass[cmp]
            dried_solids.imass[cmp] = 0
        
        # use natural gas for heat drying base on the BEAM*2024 model
        natural_gas.ivol['CH4'] = vapor.F_mass/1000*self.unit_heat*1000/self.natural_gas_HHV
        
        dried_solids.T = vapor.T = self.T
    
    def _design(self):
        self.design_results['Half dry solids flow'] = self.outs[0].F_mass/1000*24/2
        
        # kW
        self.add_power_utility(self.outs[0].F_mass/1000*self.unit_electricity)

# =============================================================================
# Sintering
# =============================================================================

# assume 2023 dollar
@cost('Product dry mass flow', 'Rotary kiln', S=400/_ton_to_tonne, cost=1250000, n=0.6,
      BM=1, CE=798)
class Sintering(ElementFlowMixin, SanUnit):
    '''
    Sintering of FePO4·2H2O-rich feed. Capital costs are based on 
    https://www.cementequipments.com/info/vertical-kiln-vs-rotary-kiln-a-cost-analysis-103128381.html (accessed 2025-12-26).
    Sintering contributes minimally to the capital cost of the entire system.
    
    Parameters
    ----------
    ins : iterable
        feed, natural_gas, air.
    outs : iterable
        product, gas.
    T : float
        Sintering precipitation, [K].
    combustion_eff : float
        Combustion efficiency to convert the energy embedded in feed to heat.
    natural_gas_HHV : float
        Higher heating value of natural gas, [MJ/m3].
    unit_electricity : float
        Electricity for sintering, [kWh/tonne feed].
    '''
    _N_ins = 3
    _N_outs = 2
    
    _units = {'Product dry mass flow': 'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 T=700 + _C_to_K, combustion_eff=0.8, natural_gas_HHV=39,
                 # 15-40 kWh/tonne feed, https://zhuanlan.zhihu.com/p/30646376322 (accessed 2025-12-23)
                 unit_electricity=30):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.T = T
        self.combustion_eff = combustion_eff
        self.natural_gas_HHV = natural_gas_HHV
        self.unit_electricity = unit_electricity
        # references for recovery calculation (assigned in systems.py)
        self.feedstock = None
    
    def _run(self):
        feed, natural_gas, air = self.ins
        product, vapor = self.outs
        
        feed.phase = 's'
        natural_gas.phase = 'g'
        air.phase = 'g'
        product.phase = 's'
        vapor.phase = 'g'
        
        product.T = vapor.T = self.T
        
        feed_copy = qs.WasteStream()
        
        feed_copy.copy_like(feed)
        
        feed_copy.imass['FePO4'] = feed_copy.imass['FePO4_2H2O']/187*151
        low_T_water_vapor = qs.WasteStream('low_T_water_vapor', H2O=feed_copy.imass['FePO4_2H2O']/187*36, T=feed_copy.T, units='kg/hr', phase='g')
        feed_copy.imass['FePO4_2H2O'] = 0
        
        rxns = []
        
        for cmp in [cmp for cmp in self.components if feed_copy.imass[cmp.ID] > 0]:
            if cmp.locked_state in ('l', 's') and (not cmp.organic or cmp.degradability=='Undegradable'):
                continue
            
            rxn = cmp.get_combustion_reaction()
            
            if rxn:
                rxns.append(rxn)
        
        combustion_rxns = self.combustion_reactions = ParallelReaction(rxns)
        
        vapor.copy_flow(feed_copy)
        combustion_rxns.force_reaction(vapor.mol)
        
        vapor.imass['H2O'] += feed.imass['FePO4_2H2O']/187*36 
        
        product.copy_flow(vapor, IDs=('Fe3', 'PO4', 'Ca2', 'Mg2', 'FePO4'), remove=True)
        
        if vapor.imass['O2'] < 0:
            air.imass['O2'] = -vapor.imass['O2']
            air.imol['N2'] = vapor.imol['N2'] = air.imol['O2']/0.21*0.79
            vapor.imass['O2'] = 0
        
        # from https://www.sciencedirect.com/science/article/pii/S0040603115001641?via%3Dihub:
        # the apparent activation energies were around 93 kJ/mol (likely the removal of the first H2O) and
        # 73 kJ/mol (likely the removal of the second H2O)
        # for endothermic reactions of removing H2O, the activation energy is higher than the total heat required (enthalpy change)
        # so the required heat to remove H2O will be less than (93 + 73)/2 = 83 kJ/mol
        # to be conservative, use 90 kJ/mol H2O
        # kJ/hr
        required_heat = product.H + product.HHV + vapor.H + 90/18*1000*feed.imass['FePO4_2H2O']/187*36
        provided_heat = feed_copy.H + feed_copy.HHV + air.H + low_T_water_vapor.H
        natural_gas_heat = (required_heat - provided_heat)/self.combustion_eff
        
        # 50-120 m3/tonne feed, https://zhuanlan.zhihu.com/p/30646376322 (accessed 2025-12-23)
        # the current calculation is around 56 m3/tonne
        natural_gas.ivol['CH4'] = natural_gas_heat/1000/self.natural_gas_HHV
    
    def _design(self):
        self.design_results['Product dry mass flow'] = self.outs[0].F_mass/1000*24
        
        # kW
        self.add_power_utility(self.ins[0].F_mass/1000*self.unit_electricity)
        
    @property
    def Fe_recovery(self):
        '''
        Fe recovery to final product (dimensionless).
        Defined as Fe in product / Fe in feedstock.
        '''
        if self.feedstock is None:
            raise AttributeError('Sintering.feedstock is not set. Assign it in systems.py.')
        Fe_in = Fe_mass(self.feedstock)
        Fe_out = Fe_mass(self.outs[0])  # product
        return Fe_out / Fe_in if Fe_in > 0 else 0

    @property
    def P_recovery(self):
        '''
        P recovery to final product (dimensionless).
        Defined as P in product / P in feedstock.
        '''
        if self.feedstock is None:
            raise AttributeError('Sintering.feedstock is not set. Assign it in systems.py.')
        P_in = P_mass(self.feedstock)
        P_out = P_mass(self.outs[0])  # product
        return P_out / P_in if P_in > 0 else 0

# TODO: update this unit based on the changes in the units above later
# =============================================================================
# FePO4_recovery
# =============================================================================

class FePO4_recovery(SanUnit):
    '''
    A unit representing the FePO4 recovery system. This unit can be connected to
    water resource recovery facility benchmarking models in 'exposan.werf'.
    A `MetalDosage` unit is needed before the `PrimaryClarifier` unit.
    The 'IdealClarifier' unit after the `PrimaryClarifier` unit needs to be replaced
    with this unit. This unit only needs operational costs as the benchmarking
    models do not have capital costs and LCA.
    
    The unit has two inputs primary sludge, and food waste. FePO4 is recovered, along
    with cake, effluent stream, and a gas stream.
    
    Operational costs can be added using the `add_OPEX` function. Operational costs
    include the following items. Relevant data and calculation methods may be found
    in `exposan.werf.models`.
        (1) food waste (may be 0),
        (2) acid,
        (3) oxidant,
        (4) natural gas,
        (5) electricity,
        (6) heat, and
        (7) solids disposal.
    
    Parameters
    ----------
    ins : iterable
        sludge, food_waste.
    outs : iterable
        product, effluent.
    food_sludge_ratio : float
        Mass ratio of organics bewteen food waste and sludge, [-].
    food_waste_moisture : float
        Moisture content of food waste, [-].
    product_purity : float
        Mass purity of produced FePO4, [-].
    f_XS_SA : float,
        Fraction of food waste organics fermented to acetate
    f_XS_vfa : float
        Fraction of food waste organics fermented to VFAs and Lactate
    f_XS_eth : float
        Fraction of food waste organics fermented to ethanol
    f_XS_cake : float
        Fraction of food waste organics lost to residue
    f_XS_gas : float
        Fraction of food waste organics lost to CO2
    cake_moisture : float
        Moisture content of cake, [-]. 
    '''
    _N_ins = 2
    _N_outs = 4

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 isdynamic=False, food_sludge_ratio=1, food_waste_moisture=0.2,
                 product_purity=1, f_XS_SA = 0.325, f_XS_vfa = 0.325, f_XS_eth = 0.02, 
                 f_XS_cake = 0.25, f_XS_gas = 0.05, cake_moisture=0.8):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic)
        self.food_sludge_ratio = food_sludge_ratio
        self.food_waste_moisture = food_waste_moisture
        self.product_purity = product_purity
        self.f_XS_SA = f_XS_SA 
        self.f_XS_vfa = f_XS_vfa 
        self.f_XS_eth = f_XS_eth 
        self.f_XS_cake = f_XS_cake 
        self.f_XS_gas = f_XS_gas 
        self.cake_moisture = cake_moisture
        
    def _run(self):
        sludge, food_waste = self.ins
        product, effluent, cake, gas = self.outs
        self._solve_streams(sludge, food_waste, product, effluent, cake, gas)
        
    def _solve_streams(self, sludge, food_waste, product, effluent, cake, gas):
        sludge.phase = 'l'
        food_waste.phase = 'l'
        product.phase = 's'
        effluent.phase = 'l'
        cake.phase = 'l'
        gas.phase = 'g'

        cmps = self.thermo.chemicals
        self.cmps = cmps
        MW_P = 30.97  # g/mol
        IDs = list(cmps.IDs)

        if self.food_sludge_ratio not in [0, 1/3, 2/3, 1, 4/3]:
            raise RuntimeError(
                'food_sludge_ratio must be one of the follow: 0, 1/3, 2/3, 1, 4/3.'
            )

        # Reset outlet streams before solving
        for ws in (product, effluent, cake, gas):
            ws.empty()

        # Food waste defined from sludge X_S and chosen moisture
        food_waste.imass['X_S'] = sludge.imass['X_S'] * self.food_sludge_ratio
        food_waste.imass['H2O'] = (
            food_waste.imass['X_S'] / (1 - self.food_waste_moisture)
        ) * self.food_waste_moisture

        metal_P_release = {
            0:   {'metal': 2.124263, 'P': 11.0693},
            1/3: {'metal': 21.13664, 'P': 44.31886},
            2/3: {'metal': 60.67351, 'P': 71.22673},
            1:   {'metal': 83.07559, 'P': 82.30645},
            4/3: {'metal': 85, 'P': 83},
        }

        P_pct = metal_P_release[self.food_sludge_ratio]['P'] / 100
        Fe_pct = metal_P_release[self.food_sludge_ratio]['metal'] / 100

        # -----------------------------
        # Store original balances
        # -----------------------------
        iP = np.array(cmps.i_P, dtype=float)
        iCOD = np.array(cmps.i_COD, dtype=float)

        sludge_imass_initial = np.array([sludge.imass[i] for i in IDs], dtype=float)  # kg/hr
        P_initial = np.sum(iP * sludge_imass_initial)      # kg-P/hr
        COD_initial = np.sum(iCOD * sludge_imass_initial)  # kg-COD/hr

        P_selected_particulates = (
            np.array(cmps.x, dtype=float)
            * np.array(cmps.i_P, dtype=float)
            * np.array([sludge.imass[i] for i in cmps.IDs], dtype=float)
            / MW_P
        )

        P_selected_solubles = (
            np.array(cmps.s, dtype=float)
            * np.array(cmps.i_P, dtype=float)
            * np.array([sludge.imass[i] for i in cmps.IDs], dtype=float)
            / MW_P
        )

        P_selected_particulates[[IDs.index(i) for i in ('X_I', 'X_FePO4') if i in IDs]] = 0.0
        P_selected_solubles[[IDs.index(i) for i in ('S_I',) if i in IDs]] = 0.0

        dFePO4_mob = sludge.imol['X_FePO4'] * min(P_pct, Fe_pct)
        dFe_from_FeOH = sludge.imol['X_FeOH'] * Fe_pct

        P_released_total = (
            np.sum(P_selected_solubles)
            + np.sum(P_selected_particulates) * P_pct
            + dFePO4_mob
        )
        Fe_released_total = dFe_from_FeOH + dFePO4_mob

        product.imol['X_FePO4'] = min(P_released_total, Fe_released_total)

        # -----------------------------
        # POST-PROCESSING MASS BALANCE
        # -----------------------------
        sludge.imol['X_FePO4'] = max(0.0, sludge.imol['X_FePO4'] - dFePO4_mob)

        P_used_total_mol = product.imol['X_FePO4']
        f_used = 0.0 if P_released_total <= 0 else min(1.0, P_used_total_mol / P_released_total)

        imass0 = np.array([sludge.imass[i] for i in IDs], dtype=float)
        P_used_part_by_cmp = (P_selected_particulates * P_pct) * f_used

        mass_consumed_part = np.zeros_like(imass0)
        mask = iP > 0
        mass_consumed_part[mask] = P_used_part_by_cmp[mask] * MW_P / iP[mask]
        mass_consumed_part = np.minimum(mass_consumed_part, imass0)

        P_used_sol_by_cmp = P_selected_solubles * f_used

        mass_consumed_sol = np.zeros_like(imass0)
        mass_consumed_sol[mask] = P_used_sol_by_cmp[mask] * MW_P / iP[mask]
        mass_consumed_sol = np.minimum(mass_consumed_sol, imass0 - mass_consumed_part)

        new_imass = imass0.copy()
        new_imass -= mass_consumed_part
        new_imass -= mass_consumed_sol

        # COD transfer to acetate
        COD_lost_total = (iCOD * (mass_consumed_part + mass_consumed_sol)).sum()
        SA_idx = IDs.index('S_A')
        if iCOD[SA_idx] == 0:
            raise RuntimeError("i_COD of S_A is zero; cannot transfer COD to S_A.")
        new_imass[SA_idx] += COD_lost_total / iCOD[SA_idx]

        for i, cmp_id in enumerate(IDs):
            sludge.imass[cmp_id] = new_imass[i]

        # -----------------------------
        # BALANCE CHECKS
        # -----------------------------
        sludge_imass_final = np.array([sludge.imass[i] for i in IDs], dtype=float)

        P_final_sludge = np.sum(iP * sludge_imass_final)
        P_in_product = product.imol['X_FePO4'] * MW_P
        P_total_after = P_final_sludge + P_in_product

        COD_final = np.sum(iCOD * sludge_imass_final)

        P_err = 0.0 if P_initial == 0 else abs((P_total_after - P_initial) / P_initial)
        COD_err = 0.0 if COD_initial == 0 else abs((COD_final - COD_initial) / COD_initial)

        tol = 0.025

        if P_err > tol:
            raise RuntimeError(
                f"P balance error exceeds tolerance: {P_err*100:.2f}%\n"
                f"P_initial={P_initial:.6f}, P_final_total={P_total_after:.6f}"
            )

        if COD_err > tol:
            raise RuntimeError(
                f"COD balance error exceeds tolerance: {COD_err*100:.2f}%\n"
                f"COD_initial={COD_initial:.6f}, COD_final={COD_final:.6f}"
            )

        # Product impurity and inert distribution
        product.imass['X_I'] = (
            product.imass['X_FePO4'] / self.product_purity * (1 - self.product_purity)
        )

        if (sludge.imass['X_I'] - product.imass['X_I']) < -1e-9:
            raise ValueError(
                f"Negative X_I mass in cake: {sludge.imass['X_I'] - product.imass['X_I']:.6e} kg/hr."
                "Check mass balance or upstream removal."
            )

        cake.imass['X_I'] = max(sludge.imass['X_I'] - product.imass['X_I'], 0.0)

        sludge.imass['X_I'] -= (product.imass['X_I'] + cake.imass['X_I'])
        sludge.imass['X_I'] = max(sludge.imass['X_I'], 0.0)

        cake.imass['H2O'] = cake.imass['X_I'] / (1 - self.cake_moisture) * self.cake_moisture

        effluent.mix_from((sludge, food_waste))

        effluent.imass['H2O'] -= cake.imass['H2O']
        if effluent.imass['H2O'] < -1e-9:
            raise ValueError('Negative flow in effluent. Reduce moisture in cake')

        # Food waste / X_S conversion and partitioning
        effluent.imass['S_A'] += effluent.imass['X_S'] * self.f_XS_SA
        effluent.imass['S_F'] += effluent.imass['X_S'] * (self.f_XS_eth + self.f_XS_vfa)
        cake.imass['X_S'] += effluent.imass['X_S'] * self.f_XS_cake

        gas.imass['S_IC'] += (
            effluent.imass['X_S']
            * self.f_XS_gas
            * effluent.components.X_S.i_C
            * gas.components.S_IC.i_mass
        )

        effluent.imass['X_S'] -= effluent.imass['X_S'] * (
            self.f_XS_SA + self.f_XS_vfa + self.f_XS_eth + self.f_XS_cake + self.f_XS_gas
        )
        
        # =====================================================================
        
        # ('S_N2', 0.21009170299490135): minimal; assume no change
        # ('S_NH4', 0.302415334699883): minimal; assume no change
        # ('S_PO4', 0.03285492291999533): updated
        # ('S_F', 2.4797507673691084): increase from X_S
        # ('S_A', 0.10028377289623293): increase from fermentation
        # ('S_I', 0.1536277619843353): minimal; assume no change
        # ('S_IC', 0.9804279473095397): minimal; assume no change
        # ('S_K', 0.32680931576984656): minimal; assume no change
        # ('S_Mg', 0.5835880638747258): minimal; assume no change
        # ('X_I', 72.5977222466668): assume removed
        # ('X_S', 146.14937577251877): convert to S_F and S_A
        # ('S_Ca', 1.6340465788492327): minimal; assume no change
        # ('X_FeOH', 10.899540897515893): updated
        # ('X_FePO4', 10.118187033007603): updated
        # ('S_Na', 1.0154432311420232): minimal; assume no change
        # ('S_Cl', 4.96049854293517): minimal; assume no change
        # ('H2O', 11632.405004568363): assume no change
    
    # TODO: uncomment
    # def _cost(self):
    #     # TODO: operating expense per hour in addition to utility cost (if no utility is actually added to this unit, then also add the utility costs here)
    #     self.add_OPEX = {
    #     '': 
    #     }
    
    # ---------------------------------------------------------------------
    # Dynamic simulation functions
    # ---------------------------------------------------------------------

    def _stream_to_state(self, ws):
        """
        Convert a WasteStream to the standard QSDsan dynamic convention:
        [C1, C2, ..., Cn, Q]
        where C is mg/L and Q is m3/d.
        """
        cmps = self.thermo.chemicals
        IDs = cmps.IDs
        Q_d = max(ws.F_vol * 24.0, 0.0)  # m3/d

        state = np.zeros(len(IDs) + 1, dtype=float)
        if Q_d > 0:
            # kg/hr -> mg/L using C = mass / flow
            # C [mg/L] = imass[kg/hr] * 24 * 1000 / Q[m3/d]
            state[:-1] = np.array([ws.imass[i] for i in IDs], dtype=float) * 24.0 * 1000.0 / Q_d
        state[-1] = Q_d
        return state

    def _state_to_stream(self, y, ws, phase='l'):
        """
        Reconstruct a WasteStream from [C1...Cn,Q] using set_flow_by_concentration.
        Assumes concentration entries are mg/L and Q is m3/d.
        """
        cmps = self.thermo.chemicals
        IDs = cmps.IDs

        C = np.asarray(y[:-1], dtype=float).copy()
        Q_d = float(y[-1])

        C[C < 0] = 0.0
        if Q_d < 0:
            Q_d = 0.0

        ws.empty()
        ws.phase = phase

        conc = {ID: c for ID, c in zip(IDs, C) if c != 0.0}

        # QSDsan flow basis here is m3/hr
        ws.set_flow_by_concentration(
            flow_tot=Q_d / 24.0,
            concentrations=conc,
            units=('m3/hr', 'mg/L'),
        )
        return ws

    def _evaluate_outlet_states(self, y_ins):
        """
        Algebraically map inlet states -> outlet states by reconstructing temporary
        inlet streams and reusing the static solver.
        Returns concatenated outlet states:
        [product_state, effluent_state, cake_state, gas_state]
        """
        thermo = self.thermo

        sludge = WasteStream(f'{self.ID}_dyn_sludge', thermo=thermo)
        food_waste = WasteStream(f'{self.ID}_dyn_food', thermo=thermo)

        product = WasteStream(f'{self.ID}_dyn_product', thermo=thermo)
        effluent = WasteStream(f'{self.ID}_dyn_effluent', thermo=thermo)
        cake = WasteStream(f'{self.ID}_dyn_cake', thermo=thermo)
        gas = WasteStream(f'{self.ID}_dyn_gas', thermo=thermo)

        # Both inlets are treated as liquid-state vectors
        self._state_to_stream(y_ins[0], sludge, phase='l')
        self._state_to_stream(y_ins[1], food_waste, phase='l')

        self._solve_streams(sludge, food_waste, product, effluent, cake, gas)

        y_product = self._stream_to_state(product)
        y_effluent = self._stream_to_state(effluent)
        y_cake = self._stream_to_state(cake)
        y_gas = self._stream_to_state(gas)

        return np.concatenate((y_product, y_effluent, y_cake, y_gas))

    def _init_state(self):
        """
        Since this unit is purely algebraic, store all four outlet states in the
        unit state vector:
        [product, effluent, cake, gas]
        """
        self._run()

        y_product = self._stream_to_state(self.outs[0])
        y_effluent = self._stream_to_state(self.outs[1])
        y_cake = self._stream_to_state(self.outs[2])
        y_gas = self._stream_to_state(self.outs[3])

        self._state = np.concatenate((y_product, y_effluent, y_cake, y_gas))
        self._dstate = np.zeros_like(self._state)
        
    def _update_state(self):
        """
        Push slices of self._state to each outlet stream.state.
        """
        n = len(self.thermo.chemicals.IDs) + 1
        product, effluent, cake, gas = self.outs

        y_product = self._state[0:n]
        y_effluent = self._state[n:2*n]
        y_cake = self._state[2*n:3*n]
        y_gas = self._state[3*n:4*n]

        if product.state is None:
            product.state = y_product.copy()
        else:
            product.state[:] = y_product

        if effluent.state is None:
            effluent.state = y_effluent.copy()
        else:
            effluent.state[:] = y_effluent

        if cake.state is None:
            cake.state = y_cake.copy()
        else:
            cake.state[:] = y_cake

        if gas.state is None:
            gas.state = y_gas.copy()
        else:
            gas.state[:] = y_gas

    def _update_dstate(self):
        """
        Push slices of self._dstate to each outlet stream.dstate.
        """
        n = len(self.thermo.chemicals.IDs) + 1
        product, effluent, cake, gas = self.outs

        dy_product = self._dstate[0:n]
        dy_effluent = self._dstate[n:2*n]
        dy_cake = self._dstate[2*n:3*n]
        dy_gas = self._dstate[3*n:4*n]

        if product.dstate is None:
            product.dstate = dy_product.copy()
        else:
            product.dstate[:] = dy_product

        if effluent.dstate is None:
            effluent.dstate = dy_effluent.copy()
        else:
            effluent.dstate[:] = dy_effluent

        if cake.dstate is None:
            cake.dstate = dy_cake.copy()
        else:
            cake.dstate[:] = dy_cake

        if gas.dstate is None:
            gas.dstate = dy_gas.copy()
        else:
            gas.dstate[:] = dy_gas

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        """
        Pure algebraic unit. Outlet states are an instantaneous function of inlet states.
        Because the mapping is nonlinear and piecewise, computing dstate numerically from inlet dy_ins.
        """
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        _evaluate_outlet_states = self._evaluate_outlet_states

        def y_t(t, y_ins, dy_ins):
            # Main algebraic evaluation
            y0 = _evaluate_outlet_states(y_ins)
            _state[:] = y0

            # Numerical directional derivative:
            # y(t+dt) ≈ F(y_ins + dy_ins*dt)
            # yins​(t+dt) ≈ yins​(t) + y˙​ins​(t)*dt
            dt_num = 1e-6  # days # TODO: dt_num can be updated.
            y_ins_pert = np.asarray(y_ins, dtype=float) + np.asarray(dy_ins, dtype=float) * dt_num

            # avoid negative concentrations/flows in perturbed inlet states
            y_ins_pert[:, :-1] = np.maximum(y_ins_pert[:, :-1], 0.0)
            y_ins_pert[:, -1] = np.maximum(y_ins_pert[:, -1], 0.0)

            y1 = _evaluate_outlet_states(y_ins_pert)
            _dstate[:] = (y1 - y0) / dt_num

            _update_state()
            _update_dstate()

        self._AE = y_t

# components in the system:
# CO2: just TEA + LCA
# CH4: just TEA + LCA
# H2O: H2O
# Fe2: TBD
# Fe3: TBD
# PO4: S_PO4
# Ca2: S_Ca 
# Mg2: S_Mg
# Org: S_F
# Ac: S_A
# Pr: S_A
# Bu: S_A
# Va: S_A
# Lac: S_F or S_A
# Etoh: S_F or S_A
# Inert: X_I?
# Residue: X_I?
# O2: not needed
# N2: not needed
# SO2: just TEA + LCA
# H2SO4: OPEX + CI
# H2O2: OPEX + CI
# FePO4: X_FePO4

# ins of the system
# fe_sludge
# heat_drying_natural_gas
# food_waste
# sintering_natural_gas
# air
# acid
# oxidant

# outs of the system
# residue
# precipitation_supernatant
# heat_drying_vapor
# fermentation_gas
# product
# sintering_vapor