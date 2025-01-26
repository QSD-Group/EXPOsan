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

# TODO: add 'ALFTSA' (after chaing the name) and 'CO2ElectrolyzerSystem'
__all__ = (
    'AcidExtraction',
    'FuelMixer',
    'HTLmixer',
    'Humidifier',
    'StruvitePrecipitation',
    'UreaSynthesis',
    'WWmixer',
    'WWTP',
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

# =============================================================================
# ALFTSA
# =============================================================================

# TODO: need test

# TODO: need update before putting the right place

# TODO: change adsorbent and the name of the unit

class ALFTSA(PressureVessel, Splitter):
    '''
    TSA using ALF as adsorbent for CO2 adsorption.
    See biosteam/units/adsorption.py.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    adsorbent : str
        Defaults to 'ALF'.
    adsorbent_cost : dict
        Defaults to 60 $/ft3.
    adsorbent_lifetime : float
        The lifetime of adsorbents. Defaults to 20 years.
    adsorbate_ID : str
        Defaults to 'CO2'.
    split : dict[str, float] or list[float], optional
        Component splits towards the effluent (0th outlet).
        Defaults to dict(O2=0.9964, N2=0.9964, CO2=0.055).
    superficial_velocity : float, optional
        Superficial velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 2160 m/h. Typical velocities are 540 to 2160 m/h for liquids [1]_.
    regeneration_velocity : float, optional
        Mean velocity of the fluid used to regenerate the bed. Defaults to 1332 m/h.
    cycle_time : float, optional
        Time at which the receiving vessel is switched. Defaults to 5 h.
        Typical 2 h [biosteam/units/adsorption.py] - 8 h [2]_.
    rho_adsorbent : float, optional
        The density of ALF. Defaults to 1441 kg/m3, Table S1 [3]_.
    adsorbent_capacity : float, optional
        Amount of CO2 that ALF can hold. Defaults to 0.1188 g/g (2.7 mmol/g), Table S7 [3]_.
    T_regeneration : float, optional
        Temperature during the regeneration phase. Defaults to 418 K, Table S8 [3]_.
    vessel_material : float, optional
        Vessel material. Defaults to 'Stainless steel 316',
    vessel_type : float, optional
        Vessel type. Defaults to 'Vertical'.
    length_unused : float, optional
        Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to 1.219 m.
    treatment_capacity : float
        Flow rate of flue gas to a single TSA unit. Defaults to 10000 m3/h.
        10000 m3/h is an assumed value, at which the diameter of TSA columns is about the half of its length.
        1.74 m3/s for a temperature swing adsorption unit [4]_.
        100000 m3/h at 50 kPa vacuum for a vacuum swing adsorption unit [5]_.
    regen_CO2 : float
        CO2 ratio in the regeneration gas. Defaults to 0.975 (CO2 purity in carbon_dioxide).
    regen_N2 : float
        N2 ratio in the regeneration gas. Defaults to 0.025*0.82/0.87 (flue gas N2:O2 = 82:5).
    regen_O2 : float
        O2 ratio in the regeneration gas. Defaults to 0.025*0.05/0.87 (flue gas N2:O2 = 82:5).
    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] https://www.chemicalprocessing.com/processing-equipment/fluid-handling/\
        article/11302111/select-the-right-valves-for-adsorption-processes-chemical-processing
    [3] Evans, H. A.; Mullangi, D.; Deng, Z.; Wang, Y.; Peh, S. B.; Wei, F.;
        Wang, J.; Brown, C. M.; Zhao, D.; Canepa, P.; Cheetham, A. K.
        Aluminum Formate, Al(HCOO)3: An Earth-Abundant, Scalable, and Highly
        Selective Material for CO2 Capture. Science Advances 2022, 8 (44),
        eade1473. https://doi.org/10.1126/sciadv.ade1473.
    [4] Jareteg, A.; Maggiolo, D.; Larsson, A.; Thunman, H.; Sasic, S.; Ström,
        H. Industrial-Scale Benzene Adsorption: Assessment of a Baseline
        One-Dimensional Temperature Swing Model against Online Industrial Data.
        Ind. Eng. Chem. Res. 2020, 59 (26), 12239–12249.
        https://doi.org/10.1021/acs.iecr.0c01590.
    [5] Webley, P. A. Adsorption Technology for CO2 Separation and Capture:
        A Perspective. Adsorption 2014, 20 (2), 225–231.
        https://doi.org/10.1007/s10450-014-9603-2.
    [6] Ntiamoah, A.; Ling, J.; Xiao, P.; Webley, P. A.; Zhai, Y. CO2 Capture by
        Temperature Swing Adsorption: Use of Hot CO2-Rich Gas for Regeneration.
        Ind. Eng. Chem. Res. 2016, 55 (3), 703–713.
        https://doi.org/10.1021/acs.iecr.5b01384.
    '''
    auxiliary_unit_names = ('heat_exchanger_regeneration',)
    
    _N_ins = 2
    _N_outs = 2
    
    def _init(self, ID='ALF_TSA', ins=(), outs=(),
              adsorbent='ALF',
              adsorbent_cost=60,
              adsorbent_lifetime=20,
              adsorbate_ID='CO2',
              split=dict(O2=0.0036, N2=0.0036, CO2=0.945),
              superficial_velocity=2160,
              regeneration_velocity=1332,
              cycle_time=5,
              rho_adsorbent=1441,
              adsorbent_capacity=0.1188,
              T_regeneration=418,
              vessel_material='Stainless steel 316',
              vessel_type='Vertical',
              length_unused=1.219,
              treatment_capacity=10000,
              regen_CO2=0.975,
              regen_N2=0.025*0.82/0.87,
              regen_O2=0.025*0.05/0.87):
        bst.Splitter._init(self, split=split)
        self.adsorbent = adsorbent
        self.adsorbent_cost = adsorbent_cost
        self.adsorbent_lifetime = adsorbent_lifetime
        self._default_equipment_lifetime = self.equipment_lifetime = {'ALF': adsorbent_lifetime}
        self.adsorbate_ID = adsorbate_ID
        self.superficial_velocity = superficial_velocity
        self.regeneration_velocity = regeneration_velocity
        self.cycle_time = cycle_time
        self.rho_adsorbent = rho_adsorbent
        self.adsorbent_capacity = adsorbent_capacity
        self.T_regeneration = T_regeneration
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_unused = length_unused
        self.treatment_capacity=treatment_capacity
        self.regen_CO2 = regen_CO2
        self.regen_N2 = regen_N2
        self.regen_O2 = regen_O2
        self.heat_exchanger_regeneration = bst.HXutility(None, None, None, thermo=self.thermo)
    
    def _run(self):
        feed, regen = self.ins
        carbon_dioxide, effluent = self.outs 
        regen.empty()
        for i in self.outs: i.empty()
        
        self.N = ceil(feed.F_vol/self.treatment_capacity)

        feed.split_to(carbon_dioxide, effluent, self.split)
        F_vol_feed = feed.F_vol/self.N
        
        superficial_velocity = self.superficial_velocity
        adsorbate_ID = self.adsorbate_ID
        
        F_mass_adsorbate = carbon_dioxide.imass[adsorbate_ID]/self.N
        
        self.diameter = diameter = 2*sqrt(F_vol_feed/(superficial_velocity*pi))
        self.area = area = pi*diameter*diameter/4
        # length of equilibrium section plus unused bed
        total_length = self.cycle_time*F_mass_adsorbate/\
            (self.adsorbent_capacity*self.rho_adsorbent*area) + self.length_unused
        # size of each column
        self.length = length = total_length/2
        self.vessel_volume = length*area
        T_original = regen.T
        # the regen should have the same gas composition as carbon_dioxide, see [6]_
        # carbon dioxide: 97.5% CO2, the rest should be 0.82/0.87 N2 and 0.05/0.87 O2 (see flue gas composition)
        regen.reset_flow(CO2=self.regen_CO2, N2=self.regen_N2, O2=self.regen_O2, phase='g', units='kg/hr')
        carbon_dioxide.T = regen.T = self.T_regeneration
        regen.F_vol = area*self.regeneration_velocity*self.N
        regen.T = T_original
    
    def _design(self):
        feed, regen = self.ins
        design_results = self.design_results
        diameter = self.diameter
        length = self.length
        design_results['Number of sets'] = self.N
        design_results['Number of reactors'] = 3
        design_results.update(self._vessel_design(feed.P*_Pa_to_psi,
                                                  diameter*_m_to_ft,
                                                  length*_m_to_ft))
        hxr = self.heat_exchanger_regeneration
        hxr.ins.empty()
        hxr.outs.empty()
        hxr.ins[0] = regen.copy()
        hxr.outs[0] = regen.copy()
        hxr.T = self.T_regeneration
        hxr.simulate()
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(design_results['Weight'],
                                                                  design_results['Diameter'],
                                                                  design_results['Length']))
        N_sets = design_results['Number of sets']
        N_reactors = design_results['Number of reactors']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_sets*N_reactors
        # 35.3147 ft3/m3
        baseline_purchase_costs[self.adsorbent] = N_sets*N_reactors*35.3147*\
            self.vessel_volume*self.adsorbent_cost

# =============================================================================
# CO2ElectrolyzerSystem
# =============================================================================

# TODO: need test

# TODO: need update before putting the right place

# cost basis year: 2010 (may be later years, but to be conservative)
# see Jouny et al. 2018 SI Excel - tab Data - cell C4
@cost(basis='Electrolyzer area', ID='Electrolyzer', units='m2',
      cost=250.25*1.75*175/1000*10*(1+1/0.65*0.35), S=1, CE=CEPCI_by_year[2010], n=1, BM=1.2)
@cost(basis='Electrolyte flow rate (ethanol)', ID='Distiller (ethanol)', units='L/min',
      cost=4162240, S=1000, CE=CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (formic acid)', ID='Distiller (formic acid)', units='L/min',
      cost=6896190, S=1000, CE=CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (methanol)', ID='Distiller (methanol)', units='L/min',
      cost=4514670, S=1000, CE=CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (propanol)', ID='Distiller (propanol)', units='L/min',
      cost=4687910, S=1000, CE=CEPCI_by_year[2010], n=0.7)
@cost(basis='Total gas flow for PSA', ID='PSA', units='m3/h',
      cost=1989043, S=1000, CE=CEPCI_by_year[2010], n=0.7)
@cost(basis='Gas product and hydrogen flow for PSA', ID='H2_PSA', units='m3/h',
      cost=1989043, S=1000, CE=CEPCI_by_year[2010], n=0.7)
class CO2ElectrolyzerSystem(SanUnit):
    '''
    CO2 electrolyzer system that converts CO2 into reduced 1C, 2C, and 3C products [1]_.
    Compared to [1]_, this system adds a PSA (if gas products) to separate H2.
    Assume the CAPEX and OPEX of the PSA for H2 separation are the same as the PSA for
    CO2 separation.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    target_product : str, optional
        target product of the CO2 reduction system, can only be 'carbon monoxide',
        'ethanol','ethylene','formic acid','methane','methanol', and 'propanol'.
        Defaults to 'formic acid'.
    current_density : float, optional
        Defaults to 0.2 A/cm2.
    cathodic_overpotential : float, optional
        Defaults to 0.454 V.
    cell_voltage : float, optional
        Defaults to 2.3 V.
    product_selectivity : float, optional
        Defaults to 0.9.
    conversion : float, optional
        Defaults to 0.5.
    operating_days_per_year: float/int
        Same to the operating day per year of the system.
    
    References
    ----------
    [1] Jouny, M.; Luc, W.; Jiao, F. General Techno-Economic Analysis of CO2
        Electrolysis Systems. Ind. Eng. Chem. Res. 2018, 57 (6), 2165–2177.
        https://doi.org/10.1021/acs.iecr.7b03514.
    '''
    _N_ins = 2
    _N_outs = 3
    
    _units= {'Electrolyzer area':'m2',
             'Electrolyte flow rate (ethanol)':'L/min',
             'Electrolyte flow rate (formic acid)':'L/min',
             'Electrolyte flow rate (methanol)':'L/min',
             'Electrolyte flow rate (propanol)':'L/min',
             'Total gas flow for PSA':'m3/h',
             'Gas product and hydrogen flow for PSA':'m3/h'}
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None, target_product='formic acid',
                 current_density=0.2, cathodic_overpotential=0.454, cell_voltage=2.3,
                 product_selectivity=0.9, conversion=0.5, operating_days_per_year=350):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.target_product = target_product
        self.current_density = current_density
        self.cathodic_overpotential = cathodic_overpotential
        self.cell_voltage = cell_voltage
        self.product_selectivity = product_selectivity
        self.conversion = conversion
        self.operating_days_per_year = operating_days_per_year
    
    def _run(self):
        carbon_dioxide, water = self.ins
        product, hydrogen, offgas = self.outs
        
        # propanol is n-propanol
        if self.target_product not in ['carbon monoxide','ethanol','ethylene',
                                       'formic acid','methane','methanol','propanol']:
            raise ValueError("target product must be in 'carbon monoxide', 'ethanol', " + 
                             "'ethylene', 'formic acid', 'methane', 'methanol', and 'propanol'")
        
        # density in kg/m3
        product_info = {'chemical': ['carbon monoxide','ethanol','ethylene',
                                     'formic acid','methane','methanol','propanol'],
                        'formula': ['CO','C2H6O','C2H4','HCOOH','CH4','CH4O','C3H8O'],
                        'elec_number': [2, 12, 12, 2, 8, 6, 18],
                        'mole_ratio': [1, 2, 2, 1, 1, 1, 3],
                        'MW': [28, 46, 28, 46, 16, 32, 60],
                        'density': [1.14, 789, 1.18, 1221, 0.656, 792, 803],
                        'state': ['gas','liq','gas','liq','gas','liq','liq'],
                        'potential': [-0.106, 0.084, 0.064, -0.25, 0.169, 0.016, 0.095],
                        'cathode_water_production': [1, 3, 4, 0, 2, 1, 5]}
        
        product_info = pd.DataFrame.from_dict(product_info)
        
        chemical_info = self.chemical_info =\
            product_info[product_info['chemical'] == self.target_product]
        
        # CO2 inlet flow rate [kg/h]
        CO2_inlet_flow_rate = carbon_dioxide.imass['CO2']
        
        # unconverted CO2 will be recycled
        CO2_converted = CO2_inlet_flow_rate*24
        
        # CO2 recycle flow rate [kg/h]
        CO2_recycle_flow_rate_kg_per_h = CO2_inlet_flow_rate*(1-self.conversion)/self.conversion
        
        # CO2 recycle flow rate [m3/h]
        # the density of CO2 is 1.98 kg/m3, Jouny et al. 2018
        self.CO2_recycle_flow_rate_m3_per_h = CO2_recycle_flow_rate_kg_per_h/1.98
        
        # production amount based on inlet CO2 [kg/day]
        product_production = CO2_converted/44/float(chemical_info['mole_ratio'])*\
            float(chemical_info['MW'])
        
        # required current [A]
        current_needed = product_production/24/3600*1000/float(chemical_info['MW'])*\
            float(chemical_info['elec_number'])*96485/self.product_selectivity
        
        # required electrolzer area [m2]
        self.electrolyzer_area = current_needed/self.current_density/10000
        
        # required power [MW]
        self.power_needed = current_needed*self.cell_voltage/1000000
        
        # gas product flow rate [m3/h]
        self.gas_product_flow_rate = 0\
            if chemical_info['state'].to_string(index=False) == 'liq'\
            else product_production/float(chemical_info['density'])/24
        
        # liquid product flow rate [m3/h]
        liquid_product_flow_rate_m3_per_h = 0\
            if chemical_info['state'].to_string(index=False) == 'gas'\
            else product_production/float(chemical_info['density'])/24
        
        # liquid product flow rate [l/min]
        liquid_product_flow_rate_l_per_min = liquid_product_flow_rate_m3_per_h*1000/60
        
        # electrolyte flow rate [l/min]
        # product-rich electrolyte is recycled until a steady-state volume concentration
        # of 10% if reached, Jouny et al. 2018
        self.electrolyte_flow_rate = liquid_product_flow_rate_l_per_min/0.1
        
        # hydrogen flow rate [mol/s]
        # 2 represents 2 e-
        hydrogen_flow_rate_mol_per_s = current_needed*(1-self.product_selectivity)/2/96485
        
        # hydrogen flow rate [m3/h]
        # hydrogen molar mass: 2, hydrogen density: 0.08375 kg/m3, Jouny et al. 2018
        hydrogen_flow_rate_m3_per_h = hydrogen_flow_rate_mol_per_s*2/1000/0.08375*3600
        
        # total gas flow [m3/h]
        self.total_gas_flow = self.CO2_recycle_flow_rate_m3_per_h +\
            self.gas_product_flow_rate + hydrogen_flow_rate_m3_per_h
        
        product.empty()
        product.imass[chemical_info['formula'].to_string(index=False)] = product_production/24
        product.phase = 'g' if chemical_info['state'].to_string(index=False) == 'gas' else 'l'
        
        # note Jouny et al. 2018 calculated process water amount as 'current_needed/4/96485*18/1000*3600'
        # change it to 'current_needed/2/96485*18/1000*3600' as every 18 g H2O can transfer 2 e-
        water.imass['H2O'] = current_needed/2/96485*18/1000*3600
        
        # we added a PSA to separate H2
        # hydrogen density: 0.08375 kg/m3, Jouny et al. 2018
        hydrogen.empty()
        hydrogen.imass['H2'] = hydrogen_flow_rate_m3_per_h*0.08375
        
        offgas.empty()
        offgas.imass['N2'] = carbon_dioxide.imass['N2']
        # H2O after distillation or PSA
        # note formic acid is produced in the cathode with no water produced
        # but still needs distillation since there must be carried-out water
        offgas.imass['H2O'] = product.F_mass/float(chemical_info['MW'])*\
            float(chemical_info['cathode_water_production'])*18
        
        offgas.imass['O2'] = carbon_dioxide.F_mass + water.F_mass -\
            product.F_mass - hydrogen.F_mass - offgas.F_mass
        offgas.phase = 'g'
    
    def _design(self):
        D = self.design_results
        D['Electrolyzer area'] = self.electrolyzer_area
        
        for liquid_product in ['ethanol','formic acid','methanol','propanol']:
            D[f'Electrolyte flow rate ({liquid_product})'] = self.electrolyte_flow_rate\
                if self.chemical_info['chemical'].to_string(index=False) == liquid_product\
                    else 0
        
        D['Total gas flow for PSA'] = 0\
            if self.CO2_recycle_flow_rate_m3_per_h == self.total_gas_flow\
                else self.total_gas_flow
        
        D['Gas product and hydrogen flow for PSA'] = 0\
            if self.gas_product_flow_rate == 0\
            else self.total_gas_flow - self.CO2_recycle_flow_rate_m3_per_h
        
        # PSA operating power: 0.25 kWh/m3, Jouny et al. 2018
        self.add_power_utility(self.power_needed*1000 +\
                               D['Total gas flow for PSA']*0.25 +\
                               D['Gas product and hydrogen flow for PSA']*0.25)
    
    def _cost(self):
        self._decorated_cost()
        
        distillation_OPEX = {'ethanol': 3463310,
                             'formic acid': 11213200,
                             'methanol': 4027820,
                             'propanol': 5610420}
        if self.chemical_info['chemical'].to_string(index=False) in list(distillation_OPEX.keys()):
            self.add_OPEX =\
                {f'{self.chemical_info["chemical"].to_string(index=False)}_distillation_OPEX':
                 self.electrolyte_flow_rate/1000*
                 distillation_OPEX[self.chemical_info['chemical'].to_string(index=False)]/
                 self.operating_days_per_year/24}