#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
    Saumitra Rai <raisaumitra9@gmail.com>
    
    Rishabh Puri <rp34@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import SanUnit,  WasteStream
from qsdsan.utils import auom
from qsdsan.unit_operations import SludgePump, SludgeCentrifuge, HXutility
from qsdsan.equipments import VerticalMixer
from biosteam import Stream
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch, pressure_vessel_material_factors as factors
from thermosteam.reaction import ParallelReaction

# =============================================================================
# Elemental accounting utilities (for Sankey / flow tracing)
# =============================================================================

MW_P = 30.97
MW_PO4 = 94.97
MW_Fe = 55.845
MW_FePO4 = 150.82
MW_FePO4_2H2O = 186.85

_C_to_K = 273.15
_316_over_304 = factors['Stainless steel 316']/factors['Stainless steel 304']
_ton_to_tonne = auom('ton').conversion_factor('tonne')

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

__all__ = (
    'ElementFlowMixin',
    'AcidogenicFermenter',
    'SludgeCentrifugeWithElementFlow',
    'SelectivePrecipitation',
    'HeatDrying',
    'Sintering',
    'FePO4_recovery'
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
            'P':  sum(P_mass(s) for s in self.ins)
        }

    @property
    def element_out(self):
        '''Total elemental outflow from the unit (kg/hr).'''
        return {
            'Fe': sum(Fe_mass(s) for s in self.outs),
            'P':  sum(P_mass(s) for s in self.outs)
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
                'Fe': Fe_mass(s)/total_in['Fe'] if total_in['Fe']> 0 else 0,
                'P':  P_mass(s)/total_in['P']  if total_in['P'] > 0 else 0
            })
        return dist

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
    fermentation_time : float
        reaction time, [hr].
    food_sludge_ratio : float
        Mass ratio of organics bewteen food waste and sludge, [-].
    food_waste_moisture : float
        Mositure content of food waste, [-].
    org_unconverted : float
        Mass ratio of total organics remained as organics (excluding VFAs and ethanol) after fermentation, [-].
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
        'Reactor volume':'m3',
        'Cycle time':'hr',
        'Batch time':'hr',
        'Loading time':'hr',
        'Total dead time':'hr',
        'Reactor duty':'kJ/hr',
        'Recirculation flow rate':'m3/hr'
    }
    
    # cleaning and unloading, [hr]
    tau_cleaning = 3
    
    # number of fermenter, [-]
    N = 2
    
    # fraction of filled tank to total tank volume, [-]
    V_wf = 0.9
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 fermentation_time=132,
                 food_sludge_ratio=1,
                 food_waste_moisture=0.74,
                 org_unconverted=0.03,
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
        self.fermentation_time = fermentation_time
        self.food_sludge_ratio = food_sludge_ratio
        self.food_waste_moisture = food_waste_moisture
        self.org_unconverted = org_unconverted
        self.org_to_gas = org_to_gas
        self.org_to_vfa = org_to_vfa
        self.org_to_ethanol = org_to_ethanol
        self.org_to_residue = org_to_residue
        self.VFA_ratio = VFA_ratio
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
        
        if self.fermentation_time not in np.arange(0, 144, 12):
            raise RuntimeError('fermentation_time must be one of the follow: 0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132.')
        
        if self.food_sludge_ratio not in [0, 1/3, 2/3, 1, 4/3]:
            raise RuntimeError('food_sludge_ratio must be one of the follow: 0, 1/3, 2/3, 1, 4/3.')
        
        food_waste.imass['Org'] = sludge.imass['Org']*self.food_sludge_ratio
        food_waste.imass['H2O'] = food_waste.imass['Org']/(1-self.food_waste_moisture) * self.food_waste_moisture
        
        fermentate.mix_from(self.ins)
        
        org_total = fermentate.imass['Org']
        
        org_input_sum = self.org_unconverted + self.org_to_gas + self.org_to_vfa + self.org_to_ethanol + self.org_to_residue
        
        fermentate.imass['Org'] = org_total*self.org_unconverted/org_input_sum
        
        # gas production
        gas.imass['CO2'] = org_total*self.org_to_gas/org_input_sum
        
        # VFAs and inorganics in fermentate
        vfa_mass = org_total*self.org_to_vfa/org_input_sum
        
        fermentate.imass['Acetic_acid'] = vfa_mass * self.VFA_ratio['Ac']
        fermentate.imass['Propionic_acid'] = vfa_mass * self.VFA_ratio['Pr']
        fermentate.imass['Butyric_acid'] = vfa_mass * self.VFA_ratio['Bu']
        fermentate.imass['Valeric_acid'] = vfa_mass * self.VFA_ratio['Va']
        fermentate.imass['Lactic_acid'] = vfa_mass * self.VFA_ratio['Lac']
        
        fermentate.imass['Ethanol'] = org_total*self.org_to_ethanol/org_input_sum
        
        fermentate.imass['Residue'] = org_total*self.org_to_residue/org_input_sum

        fermentate.imass['Fe2'] = fermentate.imass['Fe3'] * self.Fe_reduction
        fermentate.imass['Fe3'] -= fermentate.imass['Fe2']
        
        metal_P_release = {
            0: {
                0: {'metal': self.metal_release_ratio_0_0, 'P': self.P_release_ratio_0_0},
                12: {'metal': self.metal_release_ratio_0_12, 'P': self.P_release_ratio_0_12},
                24: {'metal': self.metal_release_ratio_0_24, 'P': self.P_release_ratio_0_24},
                36: {'metal': self.metal_release_ratio_0_36, 'P': self.P_release_ratio_0_36},
                48: {'metal': self.metal_release_ratio_0_48, 'P': self.P_release_ratio_0_48},
                60: {'metal': self.metal_release_ratio_0_60, 'P': self.P_release_ratio_0_60},
                72: {'metal': self.metal_release_ratio_0_72, 'P': self.P_release_ratio_0_72},
                84: {'metal': self.metal_release_ratio_0_84, 'P': self.P_release_ratio_0_84},
                96: {'metal': self.metal_release_ratio_0_96, 'P': self.P_release_ratio_0_96},
                108: {'metal': self.metal_release_ratio_0_108, 'P': self.P_release_ratio_0_108},           
                120: {'metal': self.metal_release_ratio_0_120, 'P': self.P_release_ratio_0_120},
                132: {'metal': self.metal_release_ratio_0_132, 'P': self.P_release_ratio_0_132},
                },
            1/3: {
                0: {'metal': self.metal_release_ratio_1_3_0, 'P': self.P_release_ratio_1_3_0},
                12: {'metal': self.metal_release_ratio_1_3_12, 'P': self.P_release_ratio_1_3_12},
                24: {'metal': self.metal_release_ratio_1_3_24, 'P': self.P_release_ratio_1_3_24},
                36: {'metal': self.metal_release_ratio_1_3_36, 'P': self.P_release_ratio_1_3_36},
                48: {'metal': self.metal_release_ratio_1_3_48, 'P': self.P_release_ratio_1_3_48},
                60: {'metal': self.metal_release_ratio_1_3_60, 'P': self.P_release_ratio_1_3_60},
                72: {'metal': self.metal_release_ratio_1_3_72, 'P': self.P_release_ratio_1_3_72},
                84: {'metal': self.metal_release_ratio_1_3_84, 'P': self.P_release_ratio_1_3_84},
                96: {'metal': self.metal_release_ratio_1_3_96, 'P': self.P_release_ratio_1_3_96},
                108: {'metal': self.metal_release_ratio_1_3_108, 'P': self.P_release_ratio_1_3_108},
                120: {'metal': self.metal_release_ratio_1_3_120, 'P': self.P_release_ratio_1_3_120},
                132: {'metal': self.metal_release_ratio_1_3_132, 'P': self.P_release_ratio_1_3_132},
                },
            2/3: {
                0: {'metal': self.metal_release_ratio_2_3_0, 'P': self.P_release_ratio_2_3_0},
                12: {'metal': self.metal_release_ratio_2_3_12, 'P': self.P_release_ratio_2_3_12},
                24: {'metal': self.metal_release_ratio_2_3_24, 'P': self.P_release_ratio_2_3_24},
                36: {'metal': self.metal_release_ratio_2_3_36, 'P': self.P_release_ratio_2_3_36},
                48: {'metal': self.metal_release_ratio_2_3_48, 'P': self.P_release_ratio_2_3_48},
                60: {'metal': self.metal_release_ratio_2_3_60, 'P': self.P_release_ratio_2_3_60},
                72: {'metal': self.metal_release_ratio_2_3_72, 'P': self.P_release_ratio_2_3_72},
                84: {'metal': self.metal_release_ratio_2_3_84, 'P': self.P_release_ratio_2_3_84},
                96: {'metal': self.metal_release_ratio_2_3_96, 'P': self.P_release_ratio_2_3_96},
                108: {'metal': self.metal_release_ratio_2_3_108, 'P': self.P_release_ratio_2_3_108},
                120: {'metal': self.metal_release_ratio_2_3_120, 'P': self.P_release_ratio_2_3_120},
                132: {'metal': self.metal_release_ratio_2_3_132, 'P': self.P_release_ratio_2_3_132},
                },
            1: {
                0: {'metal': self.metal_release_ratio_1_0, 'P': self.P_release_ratio_1_0},
                12: {'metal': self.metal_release_ratio_1_12, 'P': self.P_release_ratio_1_12},
                24: {'metal': self.metal_release_ratio_1_24, 'P': self.P_release_ratio_1_24},
                36: {'metal': self.metal_release_ratio_1_36, 'P': self.P_release_ratio_1_36},
                48: {'metal': self.metal_release_ratio_1_48, 'P': self.P_release_ratio_1_48},
                60: {'metal': self.metal_release_ratio_1_60, 'P': self.P_release_ratio_1_60},
                72: {'metal': self.metal_release_ratio_1_72, 'P': self.P_release_ratio_1_72},
                84: {'metal': self.metal_release_ratio_1_84, 'P': self.P_release_ratio_1_84},
                96: {'metal': self.metal_release_ratio_1_96, 'P': self.P_release_ratio_1_96},
                108: {'metal': self.metal_release_ratio_1_108, 'P': self.P_release_ratio_1_108},           
                120: {'metal': self.metal_release_ratio_1_120, 'P': self.P_release_ratio_1_120},
                132: {'metal': self.metal_release_ratio_1_132, 'P': self.P_release_ratio_1_132},
                },
            4/3: {
                132: {'metal': self.metal_release_ratio_4_3_132, 'P': self.P_release_ratio_4_3_132}
                }
        }
        
        metal_P_to_residue = 0
        
        for metal in ['Fe2','Fe3','Ca2','Mg2']:
            metal_P_to_residue += fermentate.imass[metal]*(1 - metal_P_release[self.food_sludge_ratio][self.fermentation_time]['metal'])
            fermentate.imass[metal] *= metal_P_release[self.food_sludge_ratio][self.fermentation_time]['metal']
        
        metal_P_to_residue += fermentate.imass['PO4']*(1 - metal_P_release[self.food_sludge_ratio][self.fermentation_time]['P'])
        fermentate.imass['PO4'] *= metal_P_release[self.food_sludge_ratio][self.fermentation_time]['P']
        
        fermentate.imass['Residue'] += metal_P_to_residue
        
        fermentate.T = gas.T = self.T
    
    def _design(self):
        sludge, food_waste = self.ins
        
        F_vol = sludge.F_vol + food_waste.F_vol
        
        D = self.design_results
        
        D.update(size_batch(F_vol, self.fermentation_time, self.tau_cleaning, self.N, self.V_wf))
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
    
    # reaction time, [hr]
    precipitation_time = 6
    
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
        
        # H2O2 comsumption (H2O2 + 2Fe2+ + 2H+ -> 2Fe3+ + 2H2O)
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
        
        D.update(size_batch(F_vol, self.precipitation_time, self.tau_cleaning, self.N, self.V_wf))
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
@cost(ID='Dryer 1', basis='Dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
@cost(ID='Dryer 2', basis='Dry solids flow', units='tonne/day',
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
    combustion_eff_HD : float
        Combustion efficiency to convert the energy embedded in feed to heat.
    natural_gas_HHV : float
        Higher heating value of natural gas, [MJ/m3].
    unit_electricity_HD : float
        Electricity for heat drying, [kWh/tonne dried solids].
    
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
    
    _units = {'Dry solids flow': 'tonne/day'}
    
    N = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', T=105 + _C_to_K, unit_heat=4.5,
                 combustion_eff_HD=0.8, natural_gas_HHV=39, unit_electricity_HD=214):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.T = T
        self.unit_heat = unit_heat
        self.combustion_eff_HD = combustion_eff_HD
        self.unit_electricity_HD = unit_electricity_HD
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
        natural_gas.ivol['CH4'] = vapor.F_mass/1000*self.unit_heat*1000/self.combustion_eff_HD/self.natural_gas_HHV
        
        dried_solids.T = vapor.T = self.T
    
    def _design(self):
        self.design_results['Dry solids flow'] = self.outs[0].F_mass/1000*24/self.N
        
        # kW
        self.add_power_utility(self.outs[0].F_mass/1000*self.unit_electricity_HD)

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
    combustion_eff_SI : float
        Combustion efficiency to convert the energy embedded in feed to heat.
    natural_gas_HHV : float
        Higher heating value of natural gas, [MJ/m3].
    unit_electricity_SI : float
        Electricity for sintering, [kWh/tonne feed].
    '''
    _N_ins = 3
    _N_outs = 2
    
    _units = {'Product dry mass flow': 'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 T=700 + _C_to_K, combustion_eff_SI=0.8, natural_gas_HHV=39,
                 # 15-40 kWh/tonne feed, https://zhuanlan.zhihu.com/p/30646376322 (accessed 2025-12-23)
                 unit_electricity_SI=30):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.T = T
        self.combustion_eff_SI = combustion_eff_SI
        self.natural_gas_HHV = natural_gas_HHV
        self.unit_electricity_SI = unit_electricity_SI
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
        
        product.copy_flow(vapor, IDs=('Fe3','PO4','Ca2','Mg2','FePO4'), remove=True)
        
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
        natural_gas_heat = (required_heat - provided_heat)/self.combustion_eff_SI
        
        # 50-120 m3/tonne feed, https://zhuanlan.zhihu.com/p/30646376322 (accessed 2025-12-23)
        # the current calculation is around 56 m3/tonne
        natural_gas.ivol['CH4'] = natural_gas_heat/1000/self.natural_gas_HHV
    
    def _design(self):
        self.design_results['Product dry mass flow'] = self.outs[0].F_mass/1000*24
        
        # kW
        self.add_power_utility(self.ins[0].F_mass/1000*self.unit_electricity_SI)
        
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

        # Food waste defined from food sludge ratio and chosen moisture
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
        
        def _total(streams, i):
            return sum(
                np.sum(i * np.array([s.imass[j] for j in IDs], dtype=float))
                for s in streams
            )
        P_in_total = _total((sludge, food_waste), iP)
        COD_in_total = _total((sludge, food_waste), iCOD)

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
        # Remove all mobilized FePO4 from the solid FePO4 pool
        sludge.imol['X_FePO4'] = max(0.0, sludge.imol['X_FePO4'] - dFePO4_mob)
        
        P_used_total_mol = product.imol['X_FePO4']
        f_used = 0.0 if P_released_total <= 0 else min(1.0, P_used_total_mol / P_released_total)

        dFePO4_unused = dFePO4_mob * (1 - f_used)

        # Keep mobilized-but-unused FePO4-P in the sludge liquid as soluble phosphate
        sludge.imass['S_PO4'] += (
            dFePO4_unused * MW_P / sludge.components.S_PO4.i_P
        )

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

        tol = 0.005

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

        # Remaining sludge + generated food waste first go to effluent/supernatant
        effluent.mix_from((sludge, food_waste))

        # Food waste / X_S conversion and partitioning
        XS0 = effluent.imass['X_S']
        effluent.imass['S_A'] += XS0 * self.f_XS_SA
        effluent.imass['S_F'] += XS0 * (self.f_XS_eth + self.f_XS_vfa)
        cake.imass['X_S'] += XS0 * self.f_XS_cake

        gas.imass['S_IC'] += (
        XS0
        * self.f_XS_gas
        * effluent.components.X_S.i_C
        / gas.components.S_IC.i_C
        )
        
        # Tracking COD balance
        COD_lost_to_gas_expected = XS0 * self.f_XS_gas * effluent.components.X_S.i_COD
        
        # Release P associated with X_S converted to non-P-bearing products
        P_from_XS_to_SA = XS0 * self.f_XS_SA * (
            effluent.components.X_S.i_P - effluent.components.S_A.i_P
        )

        P_from_XS_to_gas = XS0 * self.f_XS_gas * effluent.components.X_S.i_P

        P_released_from_XS = P_from_XS_to_SA + P_from_XS_to_gas

        if P_released_from_XS < -1e-12:
            raise RuntimeError(
                f"Negative P release from X_S conversion: {P_released_from_XS:.6e} kg-P/hr."
            )

        effluent.imass['S_PO4'] += P_released_from_XS / effluent.components.S_PO4.i_P

        effluent.imass['X_S'] -= XS0 * (
            self.f_XS_SA + self.f_XS_vfa + self.f_XS_eth + self.f_XS_cake + self.f_XS_gas
        )
        
        # Cake moisture after all cake solids have been assigned
        cake_dry_mass = cake.F_mass - cake.imass['H2O']
        cake.imass['H2O'] = (
            cake_dry_mass / (1 - self.cake_moisture) * self.cake_moisture
        )

        effluent.imass['H2O'] -= cake.imass['H2O']
        if effluent.imass['H2O'] < -1e-9:
            raise ValueError('Negative flow in effluent. Reduce moisture in cake')
            
        # -----------------------------
        # FINAL OUTLET BALANCE CHECKS
        # -----------------------------
        outlets = (product, effluent, cake, gas)
    
        P_out_total = _total(outlets, iP)
        COD_out_total = _total(outlets, iCOD)
    
        P_err_final = (
            0.0 if P_in_total == 0
            else abs((P_out_total - P_in_total) / P_in_total)
        )
    
        COD_err_final = (
            0.0 if COD_in_total == 0
            else abs((COD_out_total + COD_lost_to_gas_expected - COD_in_total) / COD_in_total)
        )
    
        if P_err_final > tol:
            raise RuntimeError(
                f"Final P balance error exceeds tolerance: {P_err_final*100:.2f}%\n"
                f"P_in_total={P_in_total:.6f}, P_out_total={P_out_total:.6f}"
            )
    
        if COD_err_final > tol:
            raise RuntimeError(
                f"Final COD balance error exceeds tolerance: {COD_err_final*100:.2f}%\n"
                f"COD_in_total={COD_in_total:.6f}, "
                f"COD_out_total={COD_out_total:.6f}, "
                f"COD_lost_to_gas_expected={COD_lost_to_gas_expected:.6f}"
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

        conc = {ID: c for ID, c in zip(IDs, C) if c != 0.0 and ID != 'H2O'}

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
        
        if not hasattr(self, '_dyn_tmp_streams'):
            self._dyn_tmp_streams = (
                WasteStream(f'{self.ID}_dyn_sludge_tmp', thermo=thermo),
                WasteStream(f'{self.ID}_dyn_food_tmp', thermo=thermo),
                WasteStream(f'{self.ID}_dyn_product_tmp', thermo=thermo),
                WasteStream(f'{self.ID}_dyn_effluent_tmp', thermo=thermo),
                WasteStream(f'{self.ID}_dyn_cake_tmp', thermo=thermo),
                WasteStream(f'{self.ID}_dyn_gas_tmp', thermo=thermo),
            )

        sludge, food_waste, product, effluent, cake, gas = self._dyn_tmp_streams

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
            y0 = _evaluate_outlet_states(y_ins)
            _state[:] = y0

           
            # y(t+dt) ≈ F(y_ins + dy_ins*dt)
            # yins​(t+dt) ≈ yins​(t) + y˙​ins​(t)*dt
            dt_num = 1e-6  # days # dt_num can be updated in the future.
            y_ins_pert = np.asarray(y_ins, dtype=float) + np.asarray(dy_ins, dtype=float) * dt_num

            # avoid negative concentrations/flows in perturbed inlet states (maybe have errors instead)
            y_ins_pert[:, :-1] = np.maximum(y_ins_pert[:, :-1], 0.0)
            y_ins_pert[:, -1] = np.maximum(y_ins_pert[:, -1], 0.0)

            y1 = _evaluate_outlet_states(y_ins_pert)
            _dstate[:] = (y1 - y0) / dt_num

            _update_state()
            _update_dstate()

        self._AE = y_t
        
class StruviteReactor(SanUnit):
    """Dynamic-capable struvite precipitation and recovery unit.

    The unit has two explicit inlets and two outlets.

    Inputs
    ------
    ins[0]
        Liquid sidestream to be treated (e.g., dewatering centrate).
    ins[1]
        MgCl2 dosing stream. In ``Mg_dosing_mode='external'`` the stream is
        used as supplied. In ``Mg_dosing_mode='auto'`` the unit calculates the
        required Mg dose and populates this stream for reporting.

    Outputs
    -------
    outs[0]
        Recovered wet struvite product. It can also contain captured competing
        precipitates and captured influent impurities.
    outs[1]
        Treated liquid effluent containing uncaptured precipitates.

    Notes
    -----
    * ``eff_PO4_mgL`` remains an empirical soluble-P target.
    * ``precip_yield`` is the fraction of formed struvite captured in the
      product, not the fraction of soluble P precipitated.
    * Concentrations of ``S_PO4``, ``S_NH4``, and ``S_Mg`` are assumed to be on
      elemental P, N, and Mg bases, respectively.
    * Competing precipitation is empirical but stoichiometrically and
      reactant-limited. Set the competition fractions to zero to disable it.
    """

    _N_ins = 2
    _N_outs = 2
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    _neg_tol = -1e-9

    def __init__(
            self, ID='', ins=None, outs=(), thermo=None, *,
            component_ID_NH3='S_NH4',
            component_ID_P='S_PO4',
            component_ID_Mg='S_Mg',
            component_ID_Cl='S_Cl',
            component_ID_Ca='S_Ca',
            component_ID_IC='S_IC',
            component_ID_struvite='X_struv',
            component_ID_ACP='X_ACP',
            component_ID_MgCO3='X_MgCO3',
            precip_yield,
            eff_PO4_mgL,
            Mg_dosing_mode='auto',
            target_Mg_P_molar_ratio=1.1,
            MgCl2_purity=1.0,
            ACP_P_fraction=0.0,
            MgCO3_Mg_fraction=0.0,
            competing_solid_capture_fraction=None,
            impurity_capture_fractions=None,
            product_moisture=0.0,
            product_density_kg_m3=1200.0,
            HRT_hr=1.0,
            pH_ctrl=None,
            strict_component_basis=True,
            component_basis_rtol=0.05,
            check_mass_balance=True,
            mass_balance_rtol=1e-7,
            mass_balance_atol_g_d=1e-3,
            init_with='WasteStream',
            F_BM_default=None,
            isdynamic=True,
            **kwargs,
            ):

        super().__init__(
            ID=ID,
            ins=ins,
            outs=outs,
            thermo=thermo,
            init_with=init_with,
            isdynamic=isdynamic,
            F_BM_default=F_BM_default,
            **kwargs,
        )

        self._component_ID_NH3 = component_ID_NH3
        self._component_ID_P = component_ID_P
        self._component_ID_Mg = component_ID_Mg
        self._component_ID_Cl = component_ID_Cl
        self._component_ID_Ca = component_ID_Ca
        self._component_ID_IC = component_ID_IC
        self._component_ID_struvite = component_ID_struvite
        self._component_ID_ACP = component_ID_ACP
        self._component_ID_MgCO3 = component_ID_MgCO3

        self.precip_yield = float(precip_yield)
        self.eff_PO4_mgL = float(eff_PO4_mgL)
        self.Mg_dosing_mode = str(Mg_dosing_mode).lower()
        self.target_Mg_P_molar_ratio = float(target_Mg_P_molar_ratio)
        self.MgCl2_purity = float(MgCl2_purity)
        self.ACP_P_fraction = float(ACP_P_fraction)
        self.MgCO3_Mg_fraction = float(MgCO3_Mg_fraction)
        self.competing_solid_capture_fraction = (
            self.precip_yield if competing_solid_capture_fraction is None
            else float(competing_solid_capture_fraction)
        )
        self.impurity_capture_fractions = dict(impurity_capture_fractions or {})
        self.product_moisture = float(product_moisture)
        self.product_density_kg_m3 = float(product_density_kg_m3)
        self.HRT_hr = float(HRT_hr)
        self.pH_ctrl = pH_ctrl
        self.strict_component_basis = bool(strict_component_basis)
        self.component_basis_rtol = float(component_basis_rtol)
        self.check_mass_balance = bool(check_mass_balance)
        self.mass_balance_rtol = float(mass_balance_rtol)
        self.mass_balance_atol_g_d = float(mass_balance_atol_g_d)

        self.Mg_dose_kg_hr = 0.0
        self.MgCl2_dose_kg_d = 0.0
        self.dose_MgCl2_kg_d = 0.0  # backward-compatible alias
        self.last_results = {}
        self.last_mass_balance = {}

        self._validate_parameters()

    # --- Component ID properties ---
    @property
    def component_ID_NH3(self): return self._component_ID_NH3

    @property
    def component_ID_P(self): return self._component_ID_P

    @property
    def component_ID_Mg(self): return self._component_ID_Mg

    @property
    def component_ID_Cl(self): return self._component_ID_Cl

    @property
    def component_ID_Ca(self): return self._component_ID_Ca

    @property
    def component_ID_IC(self): return self._component_ID_IC

    @property
    def component_ID_struvite(self): return self._component_ID_struvite

    @property
    def component_ID_ACP(self): return self._component_ID_ACP

    @property
    def component_ID_MgCO3(self): return self._component_ID_MgCO3

    def _validate_parameters(self):
        bounded = {
            'precip_yield': self.precip_yield,
            'ACP_P_fraction': self.ACP_P_fraction,
            'MgCO3_Mg_fraction': self.MgCO3_Mg_fraction,
            'competing_solid_capture_fraction': self.competing_solid_capture_fraction,
        }
        for name, value in bounded.items():
            if not 0.0 <= value <= 1.0:
                raise ValueError(f'{name} must be between 0 and 1; received {value}.')

        for ID, fraction in self.impurity_capture_fractions.items():
            if not 0.0 <= float(fraction) <= 1.0:
                raise ValueError(
                    f'impurity capture fraction for {ID!r} must be between 0 and 1; '
                    f'received {fraction}.'
                )

        if self.eff_PO4_mgL < 0:
            raise ValueError('eff_PO4_mgL cannot be negative.')
        if self.target_Mg_P_molar_ratio < 0:
            raise ValueError('target_Mg_P_molar_ratio cannot be negative.')
        if not 0 < self.MgCl2_purity <= 1:
            raise ValueError('MgCl2_purity must be greater than 0 and no greater than 1.')
        if not 0 <= self.product_moisture < 1:
            raise ValueError('product_moisture must be in [0, 1).')
        if self.product_density_kg_m3 <= 0:
            raise ValueError('product_density_kg_m3 must be positive.')
        if self.Mg_dosing_mode not in ('auto', 'external'):
            raise ValueError("Mg_dosing_mode must be either 'auto' or 'external'.")

    def _validate_components(self, cmps):
        IDs = set(cmps.IDs)
        required = {
            self.component_ID_NH3,
            self.component_ID_P,
            self.component_ID_Mg,
            self.component_ID_struvite,
        }
        if self.Mg_dosing_mode == 'auto':
            required.add(self.component_ID_Cl)
        if self.ACP_P_fraction > 0:
            required.update((self.component_ID_Ca, self.component_ID_ACP))
        if self.MgCO3_Mg_fraction > 0:
            required.update((self.component_ID_IC, self.component_ID_MgCO3))
        required.update(self.impurity_capture_fractions)

        missing = sorted(required.difference(IDs))
        if missing:
            raise ValueError(f'{self.ID}: missing required components: {missing}.')

        if not self.strict_component_basis:
            return

        expected_MW = {
            self.component_ID_P: 30.973761998,   # elemental P basis
            self.component_ID_NH3: 14.0067,      # elemental N basis
            self.component_ID_Mg: 24.305,        # elemental Mg basis
            self.component_ID_struvite: 245.41,  # MgNH4PO4.6H2O
        }
        if self.ACP_P_fraction > 0:
            expected_MW[self.component_ID_Ca] = 40.078
        if self.MgCO3_Mg_fraction > 0:
            expected_MW[self.component_ID_IC] = 12.011

        for ID, expected in expected_MW.items():
            actual = float(cmps[ID].MW)
            if not np.isfinite(actual) or not np.isclose(
                    actual, expected, rtol=self.component_basis_rtol, atol=0.0):
                raise ValueError(
                    f'{self.ID}: component {ID!r} has MW={actual:g}, but this unit '
                    f'expects approximately {expected:g}. This usually means the '
                    'component is not on the assumed elemental/compound mass basis.'
                )

    def _indices_and_MWs(self, cmps):
        IDs = list(cmps.IDs)
        index = {ID: i for i, ID in enumerate(IDs)}
        MW = {ID: float(cmps[ID].MW) for ID in IDs}
        return IDs, index, MW

    def _auto_dose_vector(self, M_feed_g_d, Q_feed_m3_d, index, MW):
        """Return auto-dose mass vector (g/d) and dose-water flow (m3/d)."""
        n = len(M_feed_g_d)
        M_dose = np.zeros(n, dtype=float)

        idx_P = index[self.component_ID_P]
        idx_Mg = index[self.component_ID_Mg]
        idx_Cl = index[self.component_ID_Cl]
        MW_P = MW[self.component_ID_P]
        MW_Mg = MW[self.component_ID_Mg]

        # Use the removable soluble-P load. The competition correction supplies
        # enough Mg to retain the requested Mg:P ratio after empirical MgCO3 loss.
        P_target_g_d = self.eff_PO4_mgL * max(Q_feed_m3_d, 1e-12)
        P_removable_g_d = max(M_feed_g_d[idx_P] - P_target_g_d, 0.0)
        P_for_struvite_g_d = P_removable_g_d * (1.0 - self.ACP_P_fraction)
        mol_P_for_struvite_d = P_for_struvite_g_d / MW_P

        available_fraction = max(1.0 - self.MgCO3_Mg_fraction, 1e-12)
        mol_Mg_target_d = (
            mol_P_for_struvite_d
            * self.target_Mg_P_molar_ratio
            / available_fraction
        )
        Mg_target_g_d = mol_Mg_target_d * MW_Mg
        Mg_add_g_d = max(Mg_target_g_d - M_feed_g_d[idx_Mg], 0.0)

        # Anhydrous MgCl2 molecular weight and chloride stoichiometry.
        MW_Cl = 35.453
        MW_MgCl2 = MW_Mg + 2.0 * MW_Cl
        MgCl2_pure_g_d = Mg_add_g_d * MW_MgCl2 / MW_Mg
        MgCl2_solution_g_d = MgCl2_pure_g_d / self.MgCl2_purity
        carrier_water_g_d = max(MgCl2_solution_g_d - MgCl2_pure_g_d, 0.0)

        M_dose[idx_Mg] = Mg_add_g_d
        M_dose[idx_Cl] = Mg_add_g_d * (2.0 * MW_Cl / MW_Mg)

        if 'H2O' in index:
            M_dose[index['H2O']] = carrier_water_g_d

        self.Mg_dose_kg_hr = Mg_add_g_d / 24.0 / 1000.0
        self.MgCl2_dose_kg_d = MgCl2_solution_g_d / 1000.0
        self.dose_MgCl2_kg_d = self.MgCl2_dose_kg_d

        # Approximate carrier-water volume only; dissolved salt volume is ignored.
        Q_dose_m3_d = carrier_water_g_d / 1e6
        return M_dose, Q_dose_m3_d

    def _set_auto_dose_stream(self, M_dose_g_d):
        dose = self.ins[1]
        dose.mass[:] = 0.0
        dose.mass[:] = np.asarray(M_dose_g_d, dtype=float) / 24.0 / 1000.0
        dose.phase = 'l'
        if self.pH_ctrl is not None:
            try:
                dose.pH = self.pH_ctrl
                dose._pH = self.pH_ctrl
            except AttributeError:
                pass

    def _element_totals(self, M_g_d, index, MW):
        """Return tracked elemental masses (g/d) from the selected species."""
        P_id = self.component_ID_P
        N_id = self.component_ID_NH3
        Mg_id = self.component_ID_Mg
        struv_id = self.component_ID_struvite

        totals = {
            'P': float(M_g_d[index[P_id]]),
            'N': float(M_g_d[index[N_id]]),
            'Mg': float(M_g_d[index[Mg_id]]),
        }

        m_struv = float(M_g_d[index[struv_id]])
        totals['P'] += m_struv * MW[P_id] / MW[struv_id]
        totals['N'] += m_struv * MW[N_id] / MW[struv_id]
        totals['Mg'] += m_struv * MW[Mg_id] / MW[struv_id]

        if self.component_ID_ACP in index:
            acp_id = self.component_ID_ACP
            m_acp = float(M_g_d[index[acp_id]])
            totals['P'] += m_acp * (2.0 * MW[P_id]) / MW[acp_id]
            if self.component_ID_Ca in index:
                ca_id = self.component_ID_Ca
                totals['Ca'] = float(M_g_d[index[ca_id]])
                totals['Ca'] += m_acp * (3.0 * MW[ca_id]) / MW[acp_id]

        if self.component_ID_MgCO3 in index:
            mgco3_id = self.component_ID_MgCO3
            m_mgco3 = float(M_g_d[index[mgco3_id]])
            totals['Mg'] += m_mgco3 * MW[Mg_id] / MW[mgco3_id]
            if self.component_ID_IC in index:
                ic_id = self.component_ID_IC
                totals['C'] = float(M_g_d[index[ic_id]])
                totals['C'] += m_mgco3 * MW[ic_id] / MW[mgco3_id]

        return totals

    def _check_balances(self, M_in_g_d, M_product_g_d, M_eff_g_d, index, MW):
        before = self._element_totals(M_in_g_d, index, MW)
        after = self._element_totals(M_product_g_d + M_eff_g_d, index, MW)
        errors = {}
        failed = []

        for element in sorted(set(before) | set(after)):
            b = float(before.get(element, 0.0))
            a = float(after.get(element, 0.0))
            error = a - b
            errors[element] = {
                'in_g_d': b,
                'out_g_d': a,
                'error_g_d': error,
            }
            tol = self.mass_balance_atol_g_d + self.mass_balance_rtol * max(abs(b), 1.0)
            if abs(error) > tol:
                failed.append((element, error, tol))

        total_in = float(np.sum(M_in_g_d))
        total_out = float(np.sum(M_product_g_d) + np.sum(M_eff_g_d))
        total_error = total_out - total_in
        errors['TotalMass'] = {
            'in_g_d': total_in,
            'out_g_d': total_out,
            'error_g_d': total_error,
        }
        total_tol = (
            self.mass_balance_atol_g_d
            + self.mass_balance_rtol * max(abs(total_in), 1.0)
        )
        if abs(total_error) > total_tol:
            failed.append(('TotalMass', total_error, total_tol))

        self.last_mass_balance = errors
        if self.check_mass_balance and failed:
            detail = ', '.join(
                f'{element}: error={error:.6g} g/d (tol={tol:.6g})'
                for element, error, tol in failed
            )
            raise RuntimeError(f'{self.ID}: elemental mass-balance failure: {detail}.')

    def _calculate(self, M_feed_g_d, Q_feed_m3_d, M_dose_g_d, Q_dose_m3_d,
                   index, MW):
        """Calculate product and effluent mass flows in g/d."""
        M_feed_g_d = np.maximum(np.asarray(M_feed_g_d, dtype=float), 0.0)

        if self.Mg_dosing_mode == 'auto':
            M_dose_g_d, Q_dose_m3_d = self._auto_dose_vector(
                M_feed_g_d, Q_feed_m3_d, index, MW,
            )
            self._set_auto_dose_stream(M_dose_g_d)
        else:
            M_dose_g_d = np.maximum(np.asarray(M_dose_g_d, dtype=float), 0.0)
            self.Mg_dose_kg_hr = M_dose_g_d[index[self.component_ID_Mg]] / 24.0 / 1000.0
            MW_Mg = MW[self.component_ID_Mg]
            MW_Cl = 35.453
            MW_MgCl2 = MW_Mg + 2.0 * MW_Cl
            self.MgCl2_dose_kg_d = (
                M_dose_g_d[index[self.component_ID_Mg]]
                * MW_MgCl2 / MW_Mg
                / self.MgCl2_purity
                / 1000.0
            )
            self.dose_MgCl2_kg_d = self.MgCl2_dose_kg_d

        M_in = M_feed_g_d + M_dose_g_d
        Q_in = max(float(Q_feed_m3_d + Q_dose_m3_d), 1e-12)
        M_eff = M_in.copy()
        M_product = np.zeros_like(M_in)

        P_id = self.component_ID_P
        N_id = self.component_ID_NH3
        Mg_id = self.component_ID_Mg
        struv_id = self.component_ID_struvite
        idx_P = index[P_id]
        idx_N = index[N_id]
        idx_Mg = index[Mg_id]
        idx_struv = index[struv_id]

        if (M_eff[idx_P] < self._neg_tol or
                M_eff[idx_N] < self._neg_tol or
                M_eff[idx_Mg] < self._neg_tol):
            raise ValueError(
                f'{self.ID}: negative reactive inlet load '
                f'(P={M_eff[idx_P]:g}, N={M_eff[idx_N]:g}, '
                f'Mg={M_eff[idx_Mg]:g} g/d).'
            )

        P_target_g_d = self.eff_PO4_mgL * Q_in
        P_removal_request_g_d = max(M_eff[idx_P] - P_target_g_d, 0.0)

        # --- Competing calcium-phosphate precipitation (Ca3(PO4)2 basis) ---
        ACP_formed_g_d = 0.0
        P_to_ACP_g_d = 0.0
        if self.ACP_P_fraction > 0.0:
            Ca_id = self.component_ID_Ca
            ACP_id = self.component_ID_ACP
            idx_Ca = index[Ca_id]
            idx_ACP = index[ACP_id]
            mol_P_target = P_removal_request_g_d * self.ACP_P_fraction / MW[P_id]
            mol_P_Ca_limit = (M_eff[idx_Ca] / MW[Ca_id]) / 1.5
            mol_P_to_ACP = max(min(mol_P_target, mol_P_Ca_limit), 0.0)
            mol_ACP = mol_P_to_ACP / 2.0

            P_to_ACP_g_d = mol_P_to_ACP * MW[P_id]
            Ca_to_ACP_g_d = 3.0 * mol_ACP * MW[Ca_id]
            ACP_formed_g_d = mol_ACP * MW[ACP_id]

            M_eff[idx_P] = max(M_eff[idx_P] - P_to_ACP_g_d, 0.0)
            M_eff[idx_Ca] = max(M_eff[idx_Ca] - Ca_to_ACP_g_d, 0.0)
            M_eff[idx_ACP] += ACP_formed_g_d

        # --- Competing MgCO3 precipitation, limited by S_IC ---
        MgCO3_formed_g_d = 0.0
        if self.MgCO3_Mg_fraction > 0.0:
            IC_id = self.component_ID_IC
            MgCO3_id = self.component_ID_MgCO3
            idx_IC = index[IC_id]
            idx_MgCO3 = index[MgCO3_id]
            mol_Mg_target = (
                M_eff[idx_Mg] / MW[Mg_id]
                * self.MgCO3_Mg_fraction
            )
            mol_IC_limit = M_eff[idx_IC] / MW[IC_id]
            mol_MgCO3 = max(min(mol_Mg_target, mol_IC_limit), 0.0)

            Mg_to_MgCO3_g_d = mol_MgCO3 * MW[Mg_id]
            IC_to_MgCO3_g_d = mol_MgCO3 * MW[IC_id]
            MgCO3_formed_g_d = mol_MgCO3 * MW[MgCO3_id]

            M_eff[idx_Mg] = max(M_eff[idx_Mg] - Mg_to_MgCO3_g_d, 0.0)
            M_eff[idx_IC] = max(M_eff[idx_IC] - IC_to_MgCO3_g_d, 0.0)
            M_eff[idx_MgCO3] += MgCO3_formed_g_d

        # Mineral components are represented on compound-mass bases, while
        # S_PO4/S_NH4/S_Mg/S_Ca/S_IC are elemental-mass bases. The untracked
        # O/H mass incorporated into the minerals is therefore withdrawn from
        # H2O so that the full component mass balance closes.
        mineral_matrix_water_g_d = 0.0
        if ACP_formed_g_d > 0.0:
            mineral_matrix_water_g_d += max(
                ACP_formed_g_d - P_to_ACP_g_d - Ca_to_ACP_g_d, 0.0,
            )
        if MgCO3_formed_g_d > 0.0:
            mineral_matrix_water_g_d += max(
                MgCO3_formed_g_d - Mg_to_MgCO3_g_d - IC_to_MgCO3_g_d, 0.0,
            )

        # --- Struvite formation and capture ---
        P_request_after_ACP_g_d = max(P_removal_request_g_d - P_to_ACP_g_d, 0.0)
        mol_struvite = max(min(
            P_request_after_ACP_g_d / MW[P_id],
            M_eff[idx_N] / MW[N_id],
            M_eff[idx_Mg] / MW[Mg_id],
        ), 0.0)

        P_to_struvite_g_d = mol_struvite * MW[P_id]
        N_to_struvite_g_d = mol_struvite * MW[N_id]
        Mg_to_struvite_g_d = mol_struvite * MW[Mg_id]
        struvite_formed_g_d = mol_struvite * MW[struv_id]
        mineral_matrix_water_g_d += max(
            struvite_formed_g_d
            - P_to_struvite_g_d
            - N_to_struvite_g_d
            - Mg_to_struvite_g_d,
            0.0,
        )

        if mineral_matrix_water_g_d > 0.0:
            if 'H2O' not in index:
                raise ValueError(
                    f'{self.ID}: mineral formation requires an H2O component '
                    'for compound-mass closure.'
                )
            idx_water = index['H2O']
            if M_eff[idx_water] + self._neg_tol < mineral_matrix_water_g_d:
                raise RuntimeError(
                    f'{self.ID}: insufficient H2O for mineral mass closure.'
                )
            M_eff[idx_water] = max(
                M_eff[idx_water] - mineral_matrix_water_g_d, 0.0,
            )

        M_eff[idx_P] = max(M_eff[idx_P] - P_to_struvite_g_d, 0.0)
        M_eff[idx_N] = max(M_eff[idx_N] - N_to_struvite_g_d, 0.0)
        M_eff[idx_Mg] = max(M_eff[idx_Mg] - Mg_to_struvite_g_d, 0.0)
        M_eff[idx_struv] += struvite_formed_g_d

        struvite_available_g_d = M_eff[idx_struv]
        struvite_captured_g_d = struvite_available_g_d * self.precip_yield
        M_product[idx_struv] += struvite_captured_g_d
        M_eff[idx_struv] = struvite_available_g_d - struvite_captured_g_d

        # Capture competing solids as product impurities.
        for solid_ID in (self.component_ID_ACP, self.component_ID_MgCO3):
            if solid_ID not in index:
                continue
            idx = index[solid_ID]
            captured = M_eff[idx] * self.competing_solid_capture_fraction
            M_product[idx] += captured
            M_eff[idx] -= captured

        # Capture user-selected influent/effluent solids as product impurities.
        protected = {
            struv_id,
            self.component_ID_ACP,
            self.component_ID_MgCO3,
            'H2O',
        }
        for impurity_ID, fraction in self.impurity_capture_fractions.items():
            if impurity_ID in protected:
                continue
            idx = index[impurity_ID]
            captured = M_eff[idx] * float(fraction)
            M_product[idx] += captured
            M_eff[idx] -= captured

        # Wet product: remove moisture from the liquid and add it to the product.
        dry_product_g_d = float(M_product.sum())
        product_water_g_d = 0.0
        if self.product_moisture > 0.0 and dry_product_g_d > 0.0:
            if 'H2O' not in index:
                raise ValueError(
                    f'{self.ID}: product_moisture requires an H2O component.'
                )
            idx_water = index['H2O']
            product_water_g_d = (
                dry_product_g_d
                * self.product_moisture
                / (1.0 - self.product_moisture)
            )
            if M_eff[idx_water] + self._neg_tol < product_water_g_d:
                raise RuntimeError(
                    f'{self.ID}: insufficient water to achieve the specified '
                    f'product moisture ({self.product_moisture:.3f}).'
                )
            M_product[idx_water] += product_water_g_d
            M_eff[idx_water] = max(M_eff[idx_water] - product_water_g_d, 0.0)

        wet_product_kg_d = float(M_product.sum()) / 1000.0
        Q_product_m3_d = wet_product_kg_d / self.product_density_kg_m3
        # Only entrained liquid moisture is removed from the liquid hydraulic
        # flow. Dry crystal volume is represented in the product stream but is
        # not subtracted from the aqueous flow used for concentration targets.
        Q_eff_m3_d = max(Q_in - product_water_g_d / 1e6, 1e-12)

        self._check_balances(M_in, M_product, M_eff, index, MW)

        self.last_results = {
            'influent_flow_m3_d': Q_in,
            'effluent_flow_m3_d': Q_eff_m3_d,
            'Mg_dose_kg_hr': self.Mg_dose_kg_hr,
            'MgCl2_dose_kg_d': self.MgCl2_dose_kg_d,
            'P_removal_request_g_d': P_removal_request_g_d,
            'P_to_ACP_g_d': P_to_ACP_g_d,
            'ACP_formed_g_d': ACP_formed_g_d,
            'MgCO3_formed_g_d': MgCO3_formed_g_d,
            'struvite_formed_g_d': struvite_formed_g_d,
            'struvite_captured_g_d': struvite_captured_g_d,
            'dry_product_kg_d': dry_product_g_d / 1000.0,
            'wet_product_kg_d': wet_product_kg_d,
            'product_water_kg_d': product_water_g_d / 1000.0,
            'mineral_matrix_water_kg_d': mineral_matrix_water_g_d / 1000.0,
        }
        return M_product, M_eff, Q_product_m3_d, Q_eff_m3_d

    # ------------------------------------------------------------------
    # Static mass balance
    # ------------------------------------------------------------------
    def _run(self):
        feed, dose = self.ins
        product, effluent = self.outs
        product.phase = 's'
        effluent.phase = 'l'

        cmps = feed.components
        self._validate_components(cmps)
        IDs, index, MW = self._indices_and_MWs(cmps)

        M_feed_g_d = np.asarray(feed.mass, dtype=float) * 24.0 * 1000.0
        M_dose_g_d = np.asarray(dose.mass, dtype=float) * 24.0 * 1000.0
        Q_feed_m3_d = max(float(feed.F_vol) * 24.0, 1e-12)
        Q_dose_m3_d = max(float(dose.F_vol) * 24.0, 0.0)

        M_product, M_eff, _, _ = self._calculate(
            M_feed_g_d, Q_feed_m3_d,
            M_dose_g_d, Q_dose_m3_d,
            index, MW,
        )

        product.mass[:] = M_product / 24.0 / 1000.0
        effluent.mass[:] = M_eff / 24.0 / 1000.0

        if self.pH_ctrl is not None:
            for ws in (product, effluent):
                try:
                    ws.pH = self.pH_ctrl
                    ws._pH = self.pH_ctrl
                except AttributeError:
                    pass

    # ------------------------------------------------------------------
    # Dynamic algebraic evaluation
    # ------------------------------------------------------------------
    # QSDsan natively converts dynamic states back to stream flows for liquid
    # and gas WasteStreams, but not for phase='s'. Therefore, the recovered
    # product uses a solid-specific state convention here:
    #
    #   product_state[:-1] = component mass flows [g/d]
    #   product_state[-1]  = 1.0
    #
    # The product mass vector is also synchronized explicitly in
    # _update_state(), so the recovered stream remains phase='s' and has the
    # correct mass after and during dynamic simulation.
    def _init_state(self):
        product, effluent = self.outs
        n = len(self.components)

        product.phase = 's'
        effluent.phase = 'l'

        self._product_state = np.zeros(n + 1)
        self._product_dstate = np.zeros(n + 1)
        self._effluent_state = np.zeros(n + 1)
        self._effluent_dstate = np.zeros(n + 1)

        # Solid-product convention: direct component mass flows [g/d],
        # followed by a unit scaling factor.
        self._product_state[:-1] = (
            np.asarray(product.mass, dtype=float) * 24.0 * 1000.0
        )
        self._product_state[-1] = 1.0

        # Liquid-effluent convention: concentrations [mg/L] and flow [m3/d].
        self._effluent_state[:-1] = effluent.conc
        self._effluent_state[-1] = max(float(effluent.F_vol) * 24.0, 1e-12)

        # Preserve the conventional public unit state as the liquid effluent.
        self._state = self._effluent_state
        self._dstate = self._effluent_dstate
        self._update_state()
        self._update_dstate()

    def _sync_solid_product_flow(self):
        """Synchronize the phase-'s' product flow from its dynamic state."""
        product = self.outs[0]
        product.phase = 's'

        scale = max(float(self._product_state[-1]), 0.0)
        M_product_g_d = np.maximum(self._product_state[:-1], 0.0) * scale
        product.mass[:] = M_product_g_d / 24.0 / 1000.0

    def _update_state(self):
        product, effluent = self.outs

        if product.state is None:
            product.state = self._product_state.copy()
        else:
            product.state[:] = self._product_state

        # WasteStream._state2flows() does not currently reconstruct phase-'s'
        # streams, so do it explicitly within this SanUnit.
        self._sync_solid_product_flow()

        if effluent.state is None:
            effluent.state = self._effluent_state.copy()
        else:
            effluent.state[:] = self._effluent_state

    def _update_dstate(self):
        product, effluent = self.outs
        if product.dstate is None:
            product.dstate = self._product_dstate.copy()
        else:
            product.dstate[:] = self._product_dstate
        if effluent.dstate is None:
            effluent.dstate = self._effluent_dstate.copy()
        else:
            effluent.dstate[:] = self._effluent_dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        cmps = self.components
        self._validate_components(cmps)
        IDs, index, MW = self._indices_and_MWs(cmps)

        product_state = self._product_state
        product_dstate = self._product_dstate
        effluent_state = self._effluent_state
        effluent_dstate = self._effluent_dstate
        update_state = self._update_state
        update_dstate = self._update_dstate

        def y_t(t, y_ins, dy_ins):
            Q_feed = max(float(y_ins[0, -1]), 1e-12)
            Q_dose = max(float(y_ins[1, -1]), 0.0)
            M_feed = np.maximum(y_ins[0, :-1], 0.0) * Q_feed
            M_dose = np.maximum(y_ins[1, :-1], 0.0) * Q_dose

            M_product, M_eff, _, Q_eff = self._calculate(
                M_feed, Q_feed, M_dose, Q_dose, index, MW,
            )

            # Solid product: store direct component mass flows [g/d].
            product_state[:-1] = M_product
            product_state[-1] = 1.0

            # Liquid effluent: store concentrations [mg/L] and flow [m3/d].
            effluent_state[:-1] = M_eff / max(Q_eff, 1e-12)
            effluent_state[-1] = max(Q_eff, 1e-12)

            # This remains an algebraic unit. The terminal solid product has no
            # independently integrated inventory; its state is recalculated from
            # the current inlet at every AE evaluation.
            product_dstate[:] = 0.0
            effluent_dstate[:-1] = 0.0
            effluent_dstate[-1] = float(dy_ins[:, -1].sum())

            update_state()
            update_dstate()

        self._AE = y_t

    # ------------------------------------------------------------------
    # User-facing diagnostics (no printing)
    # ------------------------------------------------------------------
    def get_recovery_results(self):
        """Return the latest reactor diagnostics as a plain dictionary."""
        results = dict(self.last_results)
        results['mass_balance'] = dict(self.last_mass_balance)
        return results

    def _design(self):
        pass

    def _cost(self):
        pass

class StruviteRedissolution(SanUnit):
    """Dynamic-capable redissolution of residual struvite fines."""

    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    _neg_tol = -1e-9

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 component_ID_P='S_PO4',
                 component_ID_NH3='S_NH4',
                 component_ID_Mg='S_Mg',
                 component_ID_struvite='X_struv',
                 k_max=2.61,
                 d_p=0.3,
                 HRT_min=5.0,
                 isdynamic=True,
                 pH_ctrl=None,
                 **kwargs):

        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with='WasteStream', isdynamic=isdynamic, **kwargs)

        self._component_ID_P = component_ID_P
        self._component_ID_NH3 = component_ID_NH3
        self._component_ID_Mg = component_ID_Mg
        self._component_ID_struvite = component_ID_struvite

        self.k_max = float(k_max)
        self.d_p = float(d_p)
        self.HRT_min = float(HRT_min)
        self.pH_ctrl = pH_ctrl

        R0 = self.d_p / 2.0
        t_diss = R0 / self.k_max
        self._f_diss = float(np.clip(
            1.0 - np.exp(-self.HRT_min / t_diss), 0.0, 1.0,
        ))

    @property
    def component_ID_P(self): return self._component_ID_P

    @property
    def component_ID_NH3(self): return self._component_ID_NH3

    @property
    def component_ID_Mg(self): return self._component_ID_Mg

    @property
    def component_ID_struvite(self): return self._component_ID_struvite

    def _run(self):
        influent, = self.ins
        effluent, = self.outs
        effluent.copy_like(influent)
        if 'H2O' in influent.components.IDs:
            effluent.imass['H2O'] = max(float(effluent.imass['H2O']), 1e-9)

        cmps = influent.components
        MW_P = cmps[self.component_ID_P].MW
        MW_NH4 = cmps[self.component_ID_NH3].MW
        MW_Mg = cmps[self.component_ID_Mg].MW
        MW_STR = cmps[self.component_ID_struvite].MW

        m_struv_in = float(influent.imass[self.component_ID_struvite])
        if m_struv_in < self._neg_tol:
            raise ValueError(
                f'{self.ID}: negative influent struvite ({m_struv_in:g} kg/hr).'
            )

        m_struv_diss = m_struv_in * self._f_diss
        mol_diss = m_struv_diss / MW_STR

        effluent.imass[self.component_ID_P] += mol_diss * MW_P
        effluent.imass[self.component_ID_NH3] += mol_diss * MW_NH4
        effluent.imass[self.component_ID_Mg] += mol_diss * MW_Mg
        if 'H2O' in influent.components.IDs:
            matrix_mass = max(MW_STR - MW_P - MW_NH4 - MW_Mg, 0.0)
            effluent.imass['H2O'] += mol_diss * matrix_mass
        effluent.imass[self.component_ID_struvite] = max(
            m_struv_in - m_struv_diss, 0.0,
        )

        if self.pH_ctrl is not None:
            try:
                effluent.pH = self.pH_ctrl
                effluent._pH = self.pH_ctrl
            except AttributeError:
                pass

    def _init_state(self):
        eff = self.outs[0]
        n = len(self.components)
        self._state = np.empty(n + 1)
        self._state[:-1] = eff.conc
        self._state[-1] = max(eff.F_vol * 24, 1e-12)
        self._dstate = np.zeros(n + 1)
        if eff.state is None:
            eff.state = self._state.copy()
        if eff.dstate is None:
            eff.dstate = self._dstate.copy()
        eff.state[:] = self._state
        eff.dstate[:] = self._dstate

    def _update_state(self):
        eff = self.outs[0]
        if eff.state is None:
            eff.state = self._state.copy()
        else:
            eff.state[:] = self._state

    def _update_dstate(self):
        eff = self.outs[0]
        if eff.dstate is None:
            eff.dstate = self._dstate.copy()
        else:
            eff.dstate[:] = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate

        cmps = self.components
        IDs = list(cmps.IDs)

        idx_P = IDs.index(self.component_ID_P)
        idx_NH4 = IDs.index(self.component_ID_NH3)
        idx_Mg = IDs.index(self.component_ID_Mg)
        idx_struv = IDs.index(self.component_ID_struvite)
        idx_H2O = IDs.index('H2O') if 'H2O' in IDs else None

        MW_P = cmps[self.component_ID_P].MW
        MW_NH4 = cmps[self.component_ID_NH3].MW
        MW_Mg = cmps[self.component_ID_Mg].MW
        MW_STR = cmps[self.component_ID_struvite].MW

        f_diss = self._f_diss
        f_rem = 1.0 - f_diss
        ratio_P = f_diss * MW_P / MW_STR
        ratio_NH4 = f_diss * MW_NH4 / MW_STR
        ratio_Mg = f_diss * MW_Mg / MW_STR
        ratio_H2O = f_diss * max(MW_STR - MW_P - MW_NH4 - MW_Mg, 0.0) / MW_STR

        def y_t(t, y_ins, dy_ins):
            Q_in = max(y_ins[0, -1], 1e-12)
            C_in = y_ins[0, :-1]
            dQ_in = dy_ins[0, -1]
            dC_in = dy_ins[0, :-1]

            c_struv = C_in[idx_struv]
            dc_struv = dC_in[idx_struv]

            C_eff = C_in.copy()
            C_eff[idx_P] += ratio_P * c_struv
            C_eff[idx_NH4] += ratio_NH4 * c_struv
            C_eff[idx_Mg] += ratio_Mg * c_struv
            if idx_H2O is not None:
                C_eff[idx_H2O] += ratio_H2O * c_struv
            C_eff[idx_struv] = f_rem * c_struv

            _state[:-1] = C_eff
            _state[-1] = Q_in

            dC_eff = dC_in.copy()
            dC_eff[idx_P] += ratio_P * dc_struv
            dC_eff[idx_NH4] += ratio_NH4 * dc_struv
            dC_eff[idx_Mg] += ratio_Mg * dc_struv
            if idx_H2O is not None:
                dC_eff[idx_H2O] += ratio_H2O * dc_struv
            dC_eff[idx_struv] = f_rem * dc_struv

            _dstate[:-1] = dC_eff
            _dstate[-1] = dQ_in

            _update_state()
            _update_dstate()

        self._AE = y_t

    def _design(self):
        pass

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