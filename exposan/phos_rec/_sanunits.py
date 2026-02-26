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
from qsdsan import SanUnit
from qsdsan.utils import auom
from qsdsan.sanunits import HXutility, SludgePump
from qsdsan.equipments import VerticalMixer
from biosteam import Stream
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch, pressure_vessel_material_factors as factors
from thermosteam.reaction import ParallelReaction

_C_to_K = 273.15
_316_over_304 = factors['Stainless steel 316'] / factors['Stainless steel 304']
_ton_to_tonne = auom('ton').conversion_factor('tonne')

__all__ = (
    'AcidogenicFermenter',
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
class AcidogenicFermenter(SanUnit):
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
    fe_reduction : float
        Fe3+ reduction ratio during acidogenic fermentation, [-].
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
    HRT = 4*24
    
    # cleaning and unloading, [hr]
    tau_cleaning = 3
    
    # number of fermenter, [-]
    N = 2
    
    # fraction of filled tank to total tank volume, [-]
    V_wf = 0.9
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 food_sludge_ratio=1,
                 food_waste_moisture=0.2,
                 org_to_gas=0.05,
                 org_to_vfa=0.65, 
                 org_to_ethanol=0.02,
                 org_to_residue=0.25,
                 VFA_ratio={'Ac': 0.5, 'Pr': 0.24, 'Bu': 0.23, 'Va': 0.02, 'Lac': 0.01},
                 fe_reduction=0.98,
                 T=37+_C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.food_sludge_ratio = food_sludge_ratio
        self.food_waste_moisture = food_waste_moisture
        self.org_to_gas = org_to_gas
        self.org_to_vfa = org_to_vfa
        self.org_to_ethanol = org_to_ethanol
        self.org_to_residue = org_to_residue
        self.VFA_ratio = VFA_ratio
        self.fe_reduction = fe_reduction
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

        fermentate.imass['Fe2'] = fermentate.imass['Fe3'] * self.fe_reduction
        fermentate.imass['Fe3'] -= fermentate.imass['Fe2']
        
        metal_P_release = {
            0: {'metal': 2.124263, 'P': 11.0693},
            1/3: {'metal': 21.13664, 'P': 44.31886},
            2/3: {'metal': 60.67351, 'P': 71.22673},
            1: {'metal': 83.07559, 'P': 82.30645},
            4/3: {'metal': 85, 'P': 83}
        }
        
        metal_P_to_residue = 0
        
        for metal in ['Fe2','Fe3','Ca2','Mg2']:
            metal_P_to_residue += fermentate.imass[metal]*(1 - metal_P_release[self.food_sludge_ratio]['metal']/100)
            fermentate.imass[metal] *= metal_P_release[self.food_sludge_ratio]['metal']/100
        
        metal_P_to_residue += fermentate.imass['PO4']*(1 - metal_P_release[self.food_sludge_ratio]['P']/100)
        fermentate.imass['PO4'] *= metal_P_release[self.food_sludge_ratio]['P']/100
        
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
class SelectivePrecipitation(SanUnit):
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
    P_recovery : float
        Phosphorus recovery ratio via of Fe3+ and PO43- precipitation, [-].
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
    HRT = 12
    
    # cleaning and unloading, [hr]
    tau_cleaning = 1
    
    # number of precipitation tank, [-]
    N = 2
    
    # fraction of filled tank to total tank volume, [-]
    V_wf = 0.9
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                acid_dose=0.003,
                P_recovery=0.80,
                T=40+_C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_dose = acid_dose
        self.P_recovery = P_recovery
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
        oxidant.imass['H2O2'] = supernatant.imass['Fe2'] / 56 / 2 * 34 * 2
        # assume the volumetic ratio between H2O2 and H2O is 3:7
        oxidant.ivol['H2O'] = oxidant.ivol['H2O2'] / 3 * 7
        
        slurry.mix_from(self.ins)
        
        # Fe2+ oxidation to Fe3+
        slurry.imass['Fe3'] += slurry.imass['Fe2']
        slurry.imass['Fe2'] = 0
        
        # FePO4 precipitation
        P_mol = slurry.imass['PO4'] / 95
        P_precip_mol = P_mol * self.P_recovery
        
        Fe_mol = slurry.imass['Fe3'] / 56
        Fe_reacted_mol = min(Fe_mol, P_precip_mol)
        
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
class HeatDrying(SanUnit):
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
class Sintering(SanUnit):
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
                 product_purity=0.99, f_XS_SA = 0.325, f_XS_vfa = 0.325, f_XS_eth = 0.02, 
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
        
        sludge.phase = 'l'
        food_waste.phase = 'l'
        product.phase = 's'
        effluent.phase = 'l'
        cake.phase = 'l'
        gas.phase = 'g'
        
        self.cmps = self.thermo.chemicals
        
        if self.food_sludge_ratio not in [0, 1/3, 2/3, 1, 4/3]:
            raise RuntimeError('food_sludge_ratio must be one of the follow: 0, 1/3, 2/3, 1, 4/3.')
        
        food_waste.imass['X_S'] = sludge.imass['X_S']*self.food_sludge_ratio
        food_waste.imass['H2O'] = food_waste.imass['X_S']/(1-self.food_waste_moisture) * self.food_waste_moisture
        
        # TODO: this works for at least G1 and H1 for now, but may need further fine-tuning
        # TODO: adjust the Fe dosage in `MetalDosage`
        
        # P released as a function of TP and food_sludge_ratio (probably not the case, therefore commented)
        
        # MW_P = 30.97  # g/mol
        # MW  = np.array([sludge.chemicals[c.ID].MW for c in cmps], dtype=float)  # g/mol
        # nuP = np.array(cmps.i_P, dtype=float) * MW / MW_P                       # mol P / mol comp
        # imol = np.array([sludge.imol[c.ID] for c in cmps], dtype=float)         # kmol/hr
        # P_total = float(imol @ nuP)                                             # kmol-P/hr
        # P_released = P_total * metal_P_release[self.food_sludge_ratio]['P'] / 100
        
        metal_P_release = {
            0: {'metal': 2.124263, 'P': 11.0693},
            1/3: {'metal': 21.13664, 'P': 44.31886},
            2/3: {'metal': 60.67351, 'P': 71.22673},
            1: {'metal': 83.07559, 'P': 82.30645},
            4/3: {'metal': 85, 'P': 83}
        }

        P_pct  = metal_P_release[self.food_sludge_ratio]['P'] / 100
        Fe_pct = metal_P_release[self.food_sludge_ratio]['metal'] / 100

        # FePO4 contains both Fe and P; therefore it would mobilize only once 
        FePO4_pct = min(P_pct, Fe_pct) # This is a conservative assumption 
        dFePO4_mob     = sludge.imol['X_FePO4'] * FePO4_pct   
        dFe_from_FeOH  = sludge.imol['X_FeOH']  * Fe_pct
        dAlPO4_mob     = sludge.imol['X_AlPO4'] * P_pct
        dPO4_from_S    = sludge.imol['S_PO4']   * 1  # completely available since it is soluble 

        Fe_available   = dFe_from_FeOH + dFePO4_mob
        PO4_available  = dPO4_from_S + dFePO4_mob + dAlPO4_mob

        dFePO4_precip = min(Fe_available, PO4_available)
        product.imol['X_FePO4'] = dFePO4_precip

        sludge.imol['X_FeOH']  -= dFe_from_FeOH
        sludge.imol['X_FePO4'] -= dFePO4_mob
        sludge.imol['X_AlPO4'] -= dAlPO4_mob
        sludge.imol['S_PO4']   -= dPO4_from_S

        # Return excess back to sludge
        sludge.imol['S_PO4']  += max(PO4_available - dFePO4_precip, 0.0)
        sludge.imol['X_FeOH'] += max(Fe_available  - dFePO4_precip, 0.0)
        
        # Balance of inert material
        product.imass['X_I'] = product.imass['X_FePO4']/self.product_purity*(1 - self.product_purity)
        
        if (sludge.imass['X_I'] - product.imass['X_I']) < -1e-9: # tol = 1e-9
            raise ValueError(f"Negative X_I mass in cake: {sludge.imass['X_I'] - product.imass['X_I']:.6e} kg/hr. "
                "Check mass balance or upstream removal.")
        cake.imass['X_I'] = max(sludge.imass['X_I'] - product.imass['X_I'], 0.0)
        
        sludge.imass['X_I'] -= (product.imass['X_I'] + cake.imass['X_I'])
        sludge.imass['X_I'] = max(sludge.imass['X_I'], 0.0)
        
        cake.imass['H2O'] = cake.imass['X_I']/(1- self.cake_moisture) * self.cake_moisture
        
        effluent.mix_from(self.ins)
        # Now since effluent = sum of flow of food waste + primary sludge, and cake.imass['H2O'] is already defined.
        effluent.imass['H2O'] -= cake.imass['H2O'] # Because effluent.mix_from(self.ins) 
        if effluent.imass['H2O'] < -1e-9: # tol = 1e-9
            raise ValueError("Negative flow in effluent. Reduce moisture in cake")
        
        effluent.imass['S_A'] += effluent.imass['X_S'] * (self.f_XS_SA) # share of acetate
        effluent.imass['S_F'] += effluent.imass['X_S'] * (self.f_XS_eth + self.f_XS_vfa) # share of ethanol + all non-acetate VFAs (and lactate)
        cake.imass['X_S'] += effluent.imass['X_S']*self.f_XS_cake
        
        #  kg SIC / hr    =    (COD kg/hr)         * f_XS_2_SIC    *          (g C/g COD)        *           (g SIC/g C) (Assuming the unit is NOT anaerobic)
        gas.imass['S_IC'] += effluent.imass['X_S'] * self.f_XS_gas * effluent.components.X_S.i_C * gas.components.S_IC.i_mass # released as CO2
        
        effluent.imass['X_S'] -= effluent.imass['X_S']*(self.f_XS_SA + self.f_XS_vfa + self.f_XS_eth + self.f_XS_cake + self.f_XS_gas)

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
    
    # TODO: check all following functions
    def _init_state(self):
        mixed = qs.WasteStream()
        mixed.mix_from(self.ins)
        self._state = np.empty(len(self.cmps)+1)
        self._state[:-1] = mixed.conc  # first n element be the component concentrations of the mixed stream
        self._state[-1] = mixed.F_vol * 24 # last element be the total volumetric flow
        self._dstate = self._state * 0.
    
    def _update_state(self):
        product, effluent, cake, gas = self.outs
        
        if product.state is None: product.state = self._state.copy()
        
        if effluent.state is None:
            effluent.state = self._state.copy()
        else:
            effluent.state[:-1] = effluent.conc
            effluent.state[-1] = effluent.F_vol * 24
        
        if cake.state is None:
            cake.state = self._state.copy()
        else:
            cake.state[:-1] = cake.conc
            cake.state[-1] = cake.F_vol * 24
            
        if gas.state is None: gas.state = self._state.copy()
        
    def _update_dstate(self):
        product, effluent, cake, gas = self.outs
        
        if product.dstate is None: product.dstate = self._dstate.copy()
        
        if effluent.dstate is None: effluent.dstate = self._dstate.copy()
        
        if cake.dstate is None: cake.dstate = self._dstate.copy()
        
        if gas.dstate is None: gas.dstate = self._dstate.copy()
    
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
        def y_t(t, y_ins, dy_ins):
            Q_ins = y_ins[:,-1]
            C_ins = y_ins[:,:-1]
            _state[-1] = Q = sum(Q_ins)
            _state[:-1] = Q_ins @ C_ins / Q
            _dstate[:] = self._state * 0.
            _update_state()
            _update_dstate()
        self._AE = y_t

# TODO: convert components to modified ASM2d components:
# S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_IC, S_K, S_Mg, S_Ca, S_Na, S_Cl, H2O
# X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, X_CaCO3, X_struv, X_newb, X_ACP, X_MgCO3, X_AlOH, X_AlPO4, X_FeOH, X_FePO4

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