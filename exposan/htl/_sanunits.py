#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd, biosteam as bst
from math import ceil, log, sqrt, pi
from qsdsan import SanUnit
from qsdsan.sanunits import Reactor
from qsdsan.utils import auom
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year, PressureVessel
from biosteam.units.splitting import Splitter

# TODO: add DAP production, anhydrous ammonia production, UAN production

# TODO: adjust the order as needed
__all__ = (
    'AcidExtraction',
    'FuelMixer',
    'HTLmixer',
    'Humidifier',
    'StruvitePrecipitation',
    'UreaSynthesis',
    'WWmixer',
    'WWTP',
    'AmineAbsorption',
    'PreStripper',
    'S2WS'
    )

yearly_operation_hour = 7920 # Jones
_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_Pa_to_psi = auom('Pa').conversion_factor('psi')
_m_to_ft = auom('m').conversion_factor('ft')

# =============================================================================
# Acid Extraction
# =============================================================================

class AcidExtraction(Reactor):
    '''
    H2SO4 is added to hydrochar from HTL to extract phosphorus.
    
    Parameters
    ----------
    ins : Iterable(stream)
        hydrochar, acid.
    outs : Iterable(stream)
        residual, extracted.
    acid_vol: float
        0.5 M H2SO4 to hydrochar ratio: mL/g.
    P_acid_recovery_ratio: float
        The ratio of phosphorus that can be extracted.
        
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
        Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
        2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
        https://doi.org/10.2172/1855835.
    '''
    _N_ins = 2
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', acid_vol=7, P_acid_recovery_ratio=0.8,
                 P=None, tau=2, V_wf=0.8, # tau: [1]
                 length_to_diameter=2, N=1, V=10, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304', # acid condition
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_vol = acid_vol
        self.P_acid_recovery_ratio = P_acid_recovery_ratio
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    def _run(self):
        
        hydrochar, acid = self.ins
        residual, extracted = self.outs
        
        self.HTL = self.ins[0]._source
        
        if hydrochar.F_mass <= 0:
            pass
        else:
            if self.HTL.hydrochar_P <= 0:
                residual.copy_like(hydrochar)
            else: 
                acid.imass['H2SO4'] = hydrochar.F_mass*self.acid_vol*0.5*98.079/1000
                # 0.5 M H2SO4 acid_vol (10 mL/1 g) hydrochar
                # 0.5 M H2SO4 density: 1.03 kg/L 
                # https://www.fishersci.se/shop/products/sulfuric-acid-0-5m-4/11943303 (accessed 2024-08-03)
                acid.imass['H2O'] = hydrochar.F_mass*self.acid_vol*1.03 -\
                                    acid.imass['H2SO4']
                
                residual.imass['Residual'] = hydrochar.F_mass - self.ins[0]._source.\
                                             hydrochar_P*self.P_acid_recovery_ratio
                
                extracted.copy_like(acid)
                extracted.imass['P'] = hydrochar.F_mass - residual.F_mass
                # assume just P can be extracted
                
                residual.phase = 's'
                
                residual.T = extracted.T = hydrochar.T
                residual.P = hydrochar.P
                # H2SO4 reacts with hydrochar to release heat and temperature will increase
            
    @property
    def residual_C(self):
        return self.ins[0]._source.hydrochar_C
    
    @property
    def residual_P(self):
        return self.ins[0]._source.hydrochar_P - self.outs[1].imass['P']
        
    def _design(self):
        self.N = ceil(self.HTL.WWTP.ins[0].F_vol/788.627455/self.V)
        # 1/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 10 m3)
        self.P = self.ins[1].P
        Reactor._design(self)

# =============================================================================
# FuelMixer
# =============================================================================

class FuelMixer(SanUnit):
    '''
    Convert gasoline to diesel or diesel to gasoline based on LHV.
    
    Parameters
    ----------
    ins: Iterable(stream)
        gasoline, diesel
    outs: Iterable(stream)
        fuel
    target: str
        The target can only be 'gasoline' or 'diesel'.
    gasoline_price: float
        Gasoline price, [$/kg].
    diesel_price: float
        Diesel price, [$/kg].
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target='diesel',
                 gasoline_gal_2_kg=2.834894885,
                 diesel_gal_2_kg=3.220628346,
                 gasoline_price=0.9388,
                 diesel_price=0.9722):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target = target
        self.gasoline_gal_2_kg = gasoline_gal_2_kg
        self.diesel_gal_2_kg = diesel_gal_2_kg
        self.gasoline_price = gasoline_price
        self.diesel_price = diesel_price

    def _run(self):
        
        gasoline, diesel = self.ins
        fuel = self.outs[0]
        target = self.target
        
        gasoline_LHV_2_diesel_LHV = (gasoline.LHV/gasoline.F_mass)/(diesel.LHV/diesel.F_mass)
        # KJ/kg gasoline:KJ/kg diesel
        
        if target == 'gasoline':
            fuel.imass['Gasoline'] = gasoline.F_mass + diesel.F_mass/gasoline_LHV_2_diesel_LHV
            fuel.T = gasoline.T
            fuel.P = gasoline.P
        elif target == 'diesel':
            fuel.imass['Diesel'] = diesel.F_mass + gasoline.F_mass*gasoline_LHV_2_diesel_LHV
            fuel.T = diesel.T
            fuel.P = diesel.P
    
    def _cost(self):
        if self.target == 'gasoline':
            self.outs[0].price = self.gasoline_price
        elif self.target == 'diesel':
            self.outs[0].price = self.diesel_price
            
    @property
    def target(self):
        return self._target
    @target.setter
    def target(self, i):
        if i not in ('gasoline', 'diesel'):
            raise ValueError('`target` must be either "gasoline" or "diesel" ',
                             f'the provided "{i}" is not valid.')
        self._target = i

# =============================================================================
# HTL mixer
# =============================================================================

class HTLmixer(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the mixture of HTL
    aqueous and AcidEx effluent.
    
    Parameters
    ----------
    ins : Iterable(stream)
        HTLaqueous, extracted
    outs : Iterable(stream)
        mixture
        
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
        
    def _run(self):
        
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        mixture.empty()
        
        mixture.imass['C'] = self.ins[0]._source.HTLaqueous_C
        mixture.imass['N'] = self.ins[0]._source.HTLaqueous_N
        mixture.imass['P'] = self.ins[0]._source.HTLaqueous_P +\
                             extracted.imass['P']
        mixture.imass['H2O'] = HTLaqueous.F_mass + extracted.F_mass -\
                               mixture.imass['C'] - mixture.imass['N'] -\
                               mixture.imass['P']
        # represented by H2O except C, N, P
        
    @property
    def pH(self):
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]

        base_mol_per_h = HTLaqueous.imass['NaOH']*1000/40
        acid_mol_per_h = extracted.imass['H2SO4']*1000/98*2
        
        volume_L_per_h = mixture.F_vol*1000
        
        base_M = base_mol_per_h/volume_L_per_h
        acid_M = acid_mol_per_h/volume_L_per_h
        
        if base_M > acid_M:
            hydrogen_ion_M = 10**-14/(base_M-acid_M)
        elif base_M == acid_M:
            hydrogen_ion_M = 10**(-7)
        else:
            hydrogen_ion_M = acid_M - base_M

        return -log(hydrogen_ion_M, 10)

# =============================================================================
# Humidifier
# =============================================================================

class Humidifier(SanUnit):
    '''
    A fake unit increases the moisture content of HTL feedstocks to 80%.
    Assume 80% of the water can be recovered from the end. Calcalute makeup water.
        feedstock: X H2O + Y dry matter
        target:    4Y H2O + Y dry matter
        water recovery = 0.8 = (4Y - X - makeup water)/4Y
        makeup water = 0.8Y - X
    
    Parameters
    ----------
    ins : Iterable(stream)
        feedstock, makeup, recycle
    outs : Iterable(stream)
        mixture
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 3
    _N_outs = 1
        
    def _run(self):
        
        feedstock, makeup, recycle = self.ins
        mixture = self.outs[0]
        
        makeup.imass['H2O'] = max(0, 0.8*(feedstock.F_mass - feedstock.imass['H2O']) - feedstock.imass['H2O'])
        
        recycle.imass['H2O'] = (feedstock.F_mass - feedstock.imass['H2O'])/0.2 - feedstock.F_mass - makeup.imass['H2O']

        mixture.mix_from(self.ins)

# =============================================================================
# Struvite Precipitation
# =============================================================================

class StruvitePrecipitation(Reactor):
    '''
    Extracted and HTL aqueous are mixed together before adding MgCl2 for struvite precipitation.
    If mol(N)<mol(P), add NH4Cl to mol(N):mol(P)=1:1.
    
    Parameters
    ----------
    ins : Iterable(stream)
        mixture, supply_MgCl2, supply_NH4Cl, base.
    outs : Iterable(stream)
        struvite, effluent.
    target_pH: float
        Target pH for struvite precipitation.
    Mg_P_ratio: float
        mol(Mg) to mol(P) ratio.   
    P_pre_recovery_ratio: float
        Ratio of phosphorus that can be precipitated out.
    HTLaqueous_NH3_N_2_total_N: float
        Ratio of NH3-N to TN in HTL aqueous phase.
        
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
        Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
        2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
        https://doi.org/10.2172/1855835.
    .. [2] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 4
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', 
                  target_pH = 9,
                  Mg_P_ratio=1,
                  P_pre_recovery_ratio=0.828, # [1]
                  HTLaqueous_NH3_N_2_total_N = 0.853, # [2]
                  P=None, tau=1, V_wf=0.8, # tau: [1]
                  length_to_diameter=2, N=1, V=20, auxiliary=False,
                  mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                  wall_thickness_factor=1,
                  vessel_material='Carbon steel', # basic condition
                  vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_pH = target_pH
        self.Mg_P_ratio = Mg_P_ratio
        self.P_pre_recovery_ratio = P_pre_recovery_ratio
        self.HTLaqueous_NH3_N_2_total_N = HTLaqueous_NH3_N_2_total_N
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
    
        
    def _run(self):
        
        mixture, supply_MgCl2, supply_NH4Cl, base = self.ins
        struvite, effluent = self.outs
        
        self.HTLmixer = self.ins[0]._source
        
        if self.HTLmixer.outs[0].imass['P'] == 0:
            effluent.copy_like(mixture)
        else:
            old_pH = self.HTLmixer.pH
            if old_pH > 9:
                base.imass['MgO'] = 0
            elif old_pH >= 7:
                OH_M = 10**(old_pH-14)
                OH_M_needed = 10**(self.target_pH-14) - OH_M
                base.imass['MgO'] = OH_M_needed/2 * 40.3044/1000*self.ins[0].F_vol*1000
            else:
                neutral_OH_M = 10**(-old_pH)
                to_target_OH_M = 10**(self.target_pH - 14)
                OH_M_needed = neutral_OH_M + to_target_OH_M
                base.imass['MgO'] = OH_M_needed/2 * 40.3044/1000*self.ins[0].F_vol*1000
            
            supply_MgCl2.imass['MgCl2'] = max((mixture.imass['P']/30.973762*self.Mg_P_ratio -\
                                            base.imass['MgO']/40.3044)*95.211, 0)
    
            if mixture.imass['P']/30.973762 > mixture.imass['N']*self.HTLaqueous_NH3_N_2_total_N/14.0067:
            # if P > N, add NH4Cl to make sure N ≥ P
                supply_NH4Cl.imass['NH4Cl'] = (mixture.imass['P']/30.973762 - mixture.imass['N']*\
                                                self.HTLaqueous_NH3_N_2_total_N/14.0067)*53.491

            struvite.imass['Struvite'] = mixture.imass['P']*\
                                          self.P_pre_recovery_ratio/\
                                          30.973762*245.41
            supply_MgCl2.phase = supply_NH4Cl.phase = base.phase = 's'
            
            effluent.copy_like(mixture)
            effluent.imass['P'] -= struvite.imass['Struvite']*30.973762/245.41
            effluent.imass['N'] += supply_NH4Cl.imass['NH4Cl']*14.0067/53.491 -\
                                    struvite.imass['Struvite']*14.0067/245.41
            effluent.imass['H2O'] = self.F_mass_in - struvite.F_mass -\
                                    effluent.imass['C'] - effluent.imass['N'] -\
                                    effluent.imass['P']
            struvite.phase = 's'    
                
            struvite.T = mixture.T
            effluent.T = mixture.T
        
    @property
    def struvite_P(self):
        return self.outs[0].imass['Struvite']*30.973762/245.41

    @property
    def struvite_N(self):
        return self.struvite_P*14.0067/30.973762

    def _design(self):
        self.N = ceil(self.HTLmixer.ins[0]._source.WWTP.ins[0].F_vol*2/788.627455/self.V)
        # 2/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 20 m3)
        self.P = self.ins[0].P
        Reactor._design(self)

# =============================================================================
# UreaSynthesis
# =============================================================================

# TODO: need test

# assume the cost in the reference paper is 2022 dollar
# the installed cost is already included, so BM=1
@cost(basis='Production capacity', ID='Urea synthesizer', units='kg/h',
      cost=4050000, S=1000*1000/365/24,
      CE=CEPCI_by_year[2022], n=0.58, BM=1)
class UreaSynthesis(SanUnit):
    '''
    A black box model of urea synthesis based on [1].
    
    Parameters
    ----------
    ins : Iterable(stream)
        ammonia, carbon_dioxide.
    outs : Iterable(stream)
        urea, water, waste.
    ratio: float
        The overall ratio between NH3 and CO2 as reactants (after considering
        recycling of unconverted reactants).
    efficiency: float
        The overall conversion efficiency (after considering recycling of
        unconverted reactants) of CO2 to urea.
    loss: float
        The loss ratio of unconverted reactants before recycling.
    
    References
    ----------
    [1] Palys, M. J.; Daoutidis, P. Techno-Economic Optimization of Renewable
     Urea Production for Sustainable Agriculture and CO2 Utilization.
     J. Phys. Energy 2023, 6 (1), 015013. https://doi.org/10.1088/2515-7655/ad0ee6.
    '''
    _N_ins = 2
    _N_outs = 3
    _units= {'Production capacity': 'kg/h'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 # TODO: add uncertainty to the following parameters with citations
                 ratio=3.5, # 3-4, Uniform
                 efficiency=0.8, # 0.7-0.9, Uniform
                 loss=0.02): # 0.01-0.03, Uniform
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ratio = ratio
        self.efficiency = efficiency
        self.loss = loss

    def _run(self):
        
        ammonia, carbon_dioxide = self.ins
        urea, water, waste = self.outs
        
        NH3_CO2_molar_ratio = (self.ratio - (self.ratio - 2*self.efficiency)*(1-self.loss))/\
                              (1 - (1 - self.efficiency)*(1-self.loss))
        
        ammonia.imass['NH3'] = carbon_dioxide.imss['CO2']/44.009*NH3_CO2_molar_ratio*17.031
        
        urea.imass['Urea'] = carbon_dioxide.imss['CO2']/44.009*self.efficiency/\
                             (1 - (1 - self.efficiency)*(1-self.loss))*60.06
        water.imass['H2O'] = carbon_dioxide.imss['CO2']/44.009*self.efficiency/\
                             (1 - (1 - self.efficiency)*(1-self.loss))*18.01528
        waste.imass['NH3'] = carbon_dioxide.imss['CO2']/44.009*(self.ratio - 2*self.efficiency)*\
                             self.loss/(1 - (1 - self.efficiency)*(1-self.loss))*17.031
        waste.imass['CO2'] = carbon_dioxide.imss['CO2']/44.009*(1 - self.efficiency)*\
                             self.loss/(1 - (1 - self.efficiency)*(1-self.loss))*44.009
        
        # convert 0.18 MWh/tonne-urea and 0.95 MWh/tonne-urea to kW
        self.power_utility.consumption = (0.18 + 0.95)*1000/1000*urea.imass['urea']
    
    def _design(self):
        
        Design = self.design_results
        Design['Production capacity'] = self.outs[0].F_mass

# =============================================================================
# WWmixer
# =============================================================================

class WWmixer(SanUnit):
    '''
    A fake unit that mixes all wastewater streams and calculates C, N, P, and H2O
    amount.
    Parameters
    ----------
    ins : Iterable(stream)
        supernatant_1, supernatant_2, memdis_ww, ht_ww
    outs : Iterable(stream)
        mixture
    '''
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
    def _run(self):
        
        supernatant, memdis_ww, ht_ww = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        
        HT = self.ins[2]._source.ins[0]._source.ins[0]._source.ins[0]._source.\
             ins[0]._source.ins[0]._source
        
        # only account for C and N from HT if they are not less than 0
        if HT.HTaqueous_C >= 0:
            mixture.imass['C'] += HT.HTaqueous_C
            mixture.imass['H2O'] -= HT.HTaqueous_C
        if HT.HTaqueous_N >=0:
            mixture.imass['N'] += HT.HTaqueous_N
            mixture.imass['H2O'] -= HT.HTaqueous_N

# =============================================================================
# WWTP
# =============================================================================

class WWTP(SanUnit):
    '''
    WWTP is a fake unit that can set up sludge biochemical compositions
    and calculate sludge elemental compositions.
    
    Parameters
    ----------
    ins : Iterable(stream)
        ww.
    outs : Iterable(stream)
        sludge, treated.
    ww_2_dry_sludge: float
        Wastewater-to-dry-sludge conversion factor, [metric ton/day/MGD].
    sludge_moisture: float
        Sludge moisture content.
    sludge_dw_ash: float
        Sludge dry weight ash content.
    sludge_afdw_lipid: float
        Sludge ash free dry weight lipid content.
    sludge_afdw_protein: float
        Sludge ash free dry weight protein content.
    lipid_2_C: float
        Lipid to carbon factor.     
    protein_2_C: float
        Protein to carbon factor.
    carbo_2_C: float
        Carbohydrate to carbon factor.
    C_2_H: float
        Carbon to hydrogen factor.
    protein_2_N: float
        Protein to nitrogen factor.
    N_2_P: float
        Nitrogen to phosphorus factor. 
    operation_hour: float
        Plant yearly operation hour, [hr/yr].
    sludge_distance: float
        Normalized sludge transportation distance, [km].
    biocrude_distance: float
        Distance between WRRFs and oil refineries, [km].
        
    References
    ----------
    .. [1] Metcalf and Eddy, Incorporated. 1991. Wastewater Engineering:
        Treatment Disposal and Reuse. New York: McGraw-Hill.
    .. [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    '''
    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', 
                 ww_2_dry_sludge=0.94, # [1]
                 sludge_moisture=0.99, sludge_dw_ash=0.257, 
                 sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, 
                 lipid_2_C=0.750,
                 protein_2_C=0.545,
                 carbo_2_C=0.400, 
                 lipid_2_H=0.125,
                 protein_2_H=0.068,
                 carbo_2_H=0.067, 
                 protein_2_N=0.159,
                 N_2_P=0.3927,
                 operation_hours=yearly_operation_hour,
                 sludge_distance=100,
                 biocrude_distance=100):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ww_2_dry_sludge = ww_2_dry_sludge
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_ash = sludge_dw_ash
        self.sludge_afdw_lipid = sludge_afdw_lipid
        self.sludge_afdw_protein = sludge_afdw_protein
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.lipid_2_H = lipid_2_H
        self.protein_2_H = protein_2_H
        self.carbo_2_H = carbo_2_H
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.operation_hours = operation_hours
        self.sludge_distance = sludge_distance
        self.biocrude_distance = biocrude_distance
    
    def _run(self):
        
        ww = self.ins[0]
        sludge, treated = self.outs

        self.sludge_afdw_carbo = round(1 - self.sludge_afdw_protein - self.sludge_afdw_lipid, 5)   
        
        if self.sludge_dw_ash >= 1:
            raise Exception ('ash can not be larger than or equal to 1')
        
        if self.sludge_afdw_protein + self.sludge_afdw_lipid > 1:
            raise Exception ('protein and lipid exceed 1')
            
        self.sludge_dw = ww.F_vol*_m3perh_to_MGD*self.ww_2_dry_sludge*1000/24

        sludge.imass['H2O'] = self.sludge_dw/(1-self.sludge_moisture)*self.sludge_moisture
        sludge.imass['Sludge_ash'] = self.sludge_dw*self.sludge_dw_ash

        sludge_afdw = self.sludge_dw*(1 - self.sludge_dw_ash)
        sludge.imass['Sludge_lipid'] = sludge_afdw*self.sludge_afdw_lipid
        sludge.imass['Sludge_protein'] = sludge_afdw*self.sludge_afdw_protein
        sludge.imass['Sludge_carbo'] = sludge_afdw*self.sludge_afdw_carbo

        treated.imass['H2O'] = ww.F_mass - sludge.F_mass
    
    @property
    def sludge_dw_protein(self):
        return self.sludge_afdw_protein*(1-self.sludge_dw_ash)
    
    @property
    def sludge_dw_lipid(self):
        return self.sludge_afdw_lipid*(1-self.sludge_dw_ash)
    
    @property
    def sludge_dw_carbo(self):
        return self.sludge_afdw_carbo*(1-self.sludge_dw_ash)
    
    @property
    def sludge_C_ratio(self):
       return self.sludge_dw_protein*self.protein_2_C + self.sludge_dw_lipid*self.lipid_2_C +\
           self.sludge_dw_carbo*self.carbo_2_C
           
    @property
    def sludge_H_ratio(self):
       return self.sludge_dw_protein*self.protein_2_H + self.sludge_dw_lipid*self.lipid_2_H +\
           self.sludge_dw_carbo*self.carbo_2_H
   
    @property
    def sludge_N_ratio(self):
       return self.sludge_dw_protein*self.protein_2_N
   
    @property
    def sludge_P_ratio(self):
       return self.sludge_N_ratio*self.N_2_P
    
    @property
    def sludge_O_ratio(self):
       return 1 - self.sludge_C_ratio - self.sludge_H_ratio -\
           self.sludge_N_ratio - self.sludge_dw_ash

    @property
    def sludge_C(self):
       return self.sludge_C_ratio*self.sludge_dw
           
    @property
    def sludge_H(self):
       return self.sludge_H_ratio*self.sludge_dw
   
    @property
    def sludge_N(self):
       return self.sludge_N_ratio*self.sludge_dw
   
    @property
    def sludge_P(self):
       return self.sludge_P_ratio*self.sludge_dw
    
    @property
    def sludge_O(self):
       return self.sludge_O_ratio*self.sludge_dw
           
    @property
    def AOSc(self):
       return (3*self.sludge_N_ratio/14.0067 + 2*self.sludge_O_ratio/15.999 -\
               self.sludge_H_ratio/1.00784)/(self.sludge_C_ratio/12.011)

    @property
    def sludge_HHV(self):
       return 100*(0.338*self.sludge_C_ratio + 1.428*(self.sludge_H_ratio -\
              self.sludge_O_ratio/8)) # [2]

    @property
    def H_C_eff(self):
        return (self.sludge_H/1.00784-2*self.sludge_O/15.999)/self.sludge_C*12.011






                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
@cost(basis='Total flow', ID='Absorber', units='kmol/hr',
      cost=4.81e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=4.3)
@cost(basis='Total flow', ID='Stripper', units='kmol/hr',
      cost=4e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=4.3)
@cost(basis='CO2 flow', ID='Pumps', units='kmol/hr',
      # 55595.96 for 613 metric tonne/hr CO2
      kW=55595.96/(613*(1000/44))*(24123*0.1186),
      cost=0.42e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.3)
@cost(basis='Total flow', ID='Condenser', units='kmol/hr',
      cost=0.27e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=4.17)
@cost(basis='Total flow', ID='Reboiler', units='kmol/hr',
      cost=0.53e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.17)
@cost(basis='Total flow', ID='Cross heat exchanger', units='kmol/hr',
      cost=2.28e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.17)
@cost(basis='Total flow', ID='Cooler', units='kmol/hr',
      cost=0.09e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.17)
@cost(basis='Total flow', ID='Makeup tank', units='kmol/hr',
      cost=0.23e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=2.3)
class AmineAbsorption(SanUnit):
    '''
    Create an AmineAbsorption unit for capture of CO2 in the flue gas using
    30 wt% aqueous monoethanolamine (MEA). Capital cost and duty basis are
    from the conventional configuration as detailed in [1]_.
    Cost is extrapolated via the 6/10th rule [2]_. Pump power usage is
    based on [2]_. Bare module factors are based on similar units in BioSTEAM.
    
    Parameters
    ----------
    ins : 
        * [0] Flue gas containing CO2
        * [1] Makeup MEA (neat), updated by the unit
        * [2] Makeup water, updated by the unit
    outs : 
        * [0] CO2-stripped vent
        * [1] Concentrated CO2
    CO2_recovery :
        Percentage of CO2 that can be captured.
    MEA_to_CO2 :
        Net usage of MEA (kg pure MEA/metric tonne CO2 captured).
        The default is 1.5 based on [1]_ and [3]_.
    heat_ratio :
        Unit duty in kJ/kg CO2.
    
    Examples
    --------
    
    >>> import biosteam as bst
    >>> import thermosteam as tmo
    >>> MEA = tmo.Chemical('MEA', search_ID='141-43-5', phase='l')
    >>> chems = tmo.Chemicals(('CO2', 'O2', 'Water', 'N2', MEA))
    >>> tmo.settings.set_thermo(chems)
    >>> flue_gas = tmo.Stream('flue_gas',
    ...                       CO2=2895,
    ...                       N2=17609,
    ...                       O2=3618,
    ...                       units='kmol/hr',
    ...                       T=48+273.15)
    >>> U1 = bst.units.AmineAbsorption('U1',
    ...                                ins=(flue_gas, 'makeup_MEA', 'makeup_water'),
    ...                                outs=('vent', 'CO2'))
    >>> U1.simulate()
    >>> U1.show()
    AmineAbsorption: U1
    ins...
    [0] flue_gas
        phase: 'l', T: 321.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  2.9e+03
                        O2   3.62e+03
                        N2   1.76e+04
    [1] makeup_MEA
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): MEA  2.82
    [2] makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  22.3
    outs...
    [0] vent
        phase: 'l', T: 313.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  290
                        O2   3.62e+03
                        N2   1.76e+04
    [1] CO2
        phase: 'l', T: 313.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  2.61e+03
    >>> U1.results()
    Amine absorption                            Units       U1
    Electricity         Power                      kW 1.23e+03
                        Cost                   USD/hr     96.4
    Low pressure steam  Duty                    kJ/hr 4.36e+08
                        Flow                  kmol/hr 1.12e+04
                        Cost                   USD/hr 2.67e+03
    Design              Total flow            kmol/hr 2.41e+04
                        CO2 flow              kmol/hr 2.61e+03
    Purchase cost       Makeup tank               USD  2.5e+05
                        Cooler                    USD 9.78e+04
                        Cross heat exchanger      USD 2.48e+06
                        Reboiler                  USD 5.76e+05
                        Condenser                 USD 2.94e+05
                        Pumps                     USD  1.2e+05
                        Stripper                  USD 4.35e+06
                        Absorber                  USD 5.23e+06
    Total purchase cost                           USD 1.34e+07
    Utility cost                               USD/hr 2.77e+03
    
    References
    ----------
    .. [1] Karimi et al., Capital Costs and Energy Considerations of Different
        Alternative Stripper Configurations for Post Combustion CO2 Capture.
        Chemical Engineering Research and Design 2011, 89 (8), 1229–1236.
        https://doi.org/10.1016/j.cherd.2011.03.005.
    
    .. [2] Carminati et al., Bioenergy and Full Carbon Dioxide Sinking in
        Sugarcane-Biorefinery with Post-Combustion Capture and Storage:
        Techno-Economic Feasibility. Applied Energy 2019, 254, 113633.
        https://doi.org/10.1016/j.apenergy.2019.113633.
        
    .. [3] Ramezan et al., Carbon Dioxide Capture from Existing Coal-Fired Power
        Plants; DOE/NETL-401/110907; National Energy Technology Laboratory, 2007.
    '''
    
    _N_ins = 3
    _N_outs = 2
    _units = {'Total flow': 'kmol/hr',
              'CO2 flow':'kmol/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', CO2_recovery=0.9, MEA_to_CO2=1.5, heat_ratio=3611):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.CO2_recovery = CO2_recovery
        self.MEA_to_CO2 = MEA_to_CO2
        self.heat_ratio = heat_ratio
    
    def _run(self):
        flue_gas, MEA, water = self.ins
        vent, CO2 = self.outs
        vent.copy_like(flue_gas)
        CO2.imol['CO2'] = self.CO2_recovery * flue_gas.imol['CO2']
        vent.imol['CO2'] = flue_gas.imol['CO2'] - CO2.imol['CO2']
        MEA.imass['MEA'] = self.MEA_to_CO2 * CO2.imass['CO2']/1000
        water.imass['Water'] = MEA.imass['MEA'] / 0.3 * (1-0.3)
        vent.T = CO2.T = 273.15 + 40
        
    def _design(self):
        self.design_results['Total flow'] = self.ins[0].F_mol
        self.design_results['CO2 flow'] = self.outs[1].F_mol
        duty = self.heat_ratio * self.outs[1].F_mass
        self.add_heat_utility(duty, T_in=self.ins[0].T)

# =============================================================================
# PreStripper
# =============================================================================

class PreStripper(SanUnit):
    '''
    Calculate the NH3 concentration in the influent to the stripper.
    
    Parameters
    ----------
    ins : Iterable(stream)
        influent.
    outs : Iterable(stream)
        effluent.
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 influent_pH=8.16, # CHG effluent pH: 8.16 ± 0.25 [1]
                 target_pH=11.25): # 2 unit higher than pKa (9.25)
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.influent_pH = influent_pH
        self.target_pH = target_pH
    
    def _run(self):
        
        influent, base = self.ins
        effluent = self.outs[0]
        
        NaOH_conc = 10**(self.target_pH - 14) - 10**(self.influent_pH - 14)
        NaOH_mol = NaOH_conc*self.ins[0].F_mass
        base.imass['NaOH'] = NaOH_mol*39.997/1000
        
        self.CHG = self.ins[0]._source.ins[0]._source.ins[0]._source
        
        effluent.imass['NH3'] = self.CHG.CHGout_N/14.0067*17.031
        effluent.imass['H2O'] = influent.F_mass + base.F_mass - effluent.imass['NH3'] 
        
# =============================================================================
# S2WS
# =============================================================================
class S2WS(SanUnit):
    '''
    S2WS: Stream to WasteStream.
    A fake unit that enables the calculation of LCA for 'Stream'
    by converting 'Stream' to 'WasteStream'.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream'):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with)

    def _run(self):
        inlet = self.ins[0]
        outlet = self.outs[0]
        # if inlet is a 'MultiStream', for now, just assume the outlet
        # contains water that has the same weight as the inlet
        try:
            outlet.copy_like(inlet)
        except TypeError:
            outlet.imass['H2O'] = inlet.F_mass