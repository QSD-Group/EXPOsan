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

import biosteam as bst
from qsdsan import SanUnit, Construction
import qsdsan.sanunits as qsu
from biosteam.units import IsothermalCompressor, Flash, BinaryDistillation
from exposan.htl._components import create_components
from biosteam.units.design_tools import PressureVessel, flash_vessel_design
from biosteam.units.design_tools.cost_index import CEPCI_by_year as CEPCI
from math import pi, ceil, log, floor
from biosteam.exceptions import DesignError, bounds_warning, DesignWarning
from biosteam import Stream
from biosteam.units.decorators import cost
from thermosteam import indexer, equilibrium
from exposan.htl._process_settings import load_process_settings
from qsdsan.utils import auom, select_pipe
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from warnings import warn

__all__ = ('Reactor',
           'WWTP',
           'HTL',
           'AcidExtraction',
           'HTLmixer',
           'HTLsplitter',
           'StruvitePrecipitation',
           'CHG',
           'MembraneDistillation',
           'HT',
           'HC',
           'WWmixer',
           'PhaseChanger',
           'FuelMixer',
           'HTLpump',
           'HTLHX',
           'HTL_sludge_centrifuge',
           'HTLCHP',
           'HTLCHP_2',
           'HTL_storage_tank',
           'HTLcompressor',
           'HTLflash',
           'HTLdistillation',
           'HTLHXN')

cmps = create_components()

load_process_settings()

yearly_operation_hour = 7920 # Jones

_m3_to_gal = auom('m3').conversion_factor('gallon')

_lb_to_kg = auom('lb').conversion_factor('kg')

_ft3_to_gal = auom('ft3').conversion_factor('gallon')

_ft3_to_m3 = auom('ft3').conversion_factor('m3')

_in_to_ft = auom('in').conversion_factor('ft')

_m_to_ft = auom('m').conversion_factor('ft')

_Pa_to_psi = auom('Pa').conversion_factor('psi')

_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')

_m3perh_to_mmscfd = 1/1177.17 # H2

# =============================================================================
# Reactor
# =============================================================================

class Reactor(SanUnit, PressureVessel, isabstract=True):
    '''
    Create an abstract class for reactor unit, purchase cost of the reactor
    is based on volume calculated by residence time.
    Parameters
    ----------
    ins: stream
        Inlet.
    outs: stream
        Outlet.
    tau: float
        Residence time, [hr].
    V_wf: float
        Fraction of working volume over total volume.
    length_to_diameter : float
        Reactor length to diameter ratio.
    N: int
        Number of reactor.
    V: float
        Volume of reactor, [m3].
    auxiliary: bool
        Whether or not the reactor is an auxiliart unit.      
    mixing_intensity: float
        Mechanical mixing intensity, [/s].
    kW_per_m3: float
        Power usage of agitator
        (converted from 0.5 hp/1000 gal as in [1]).
        If mixing_intensity is provided, this will be calculated based on
        the mixing_intensity and viscosity of the influent mixture as in [2]_
    wall_thickness_factor=1: float
        A safety factor to scale up the calculated minimum wall thickness.
    vessel_material : str, optional
        Vessel material. Default to 'Stainless steel 316'.
    vessel_type : str, optional
        Vessel type. Can only be 'Horizontal' or 'Vertical'.
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Cost Accounting and Capital Cost Estimation. In Product
        and Process Design Principles; Wiley, 2017; pp 470.
    .. [2] Shoener et al. Energy Positive Domestic Wastewater Treatment:
        The Roles of Anaerobic and Phototrophic Technologies.
        Environ. Sci.: Processes Impacts 2014, 16 (6), 1204–1222.
        https://doi.org/10.1039/C3EM00711A.
    '''
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False

    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Single reactor volume': 'm3',
              'Reactor volume': 'm3'}

    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _Vmax = pi/4*(20**2)*40/35.3147
    
    _F_BM_default = PressureVessel._F_BM_default

    def __init__(self, ID='', ins=None, outs=(), *,
                 P=101325, tau=0.5, V_wf=0.8,
                 length_to_diameter=2, N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):

        SanUnit.__init__(self, ID, ins, outs)
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
        
    def _design(self):
        Design = self.design_results
        if not self.auxiliary:
        # if auxiliary, in our system, only K/O drum whose N, and V are provided
        # do not need to deal with self.F_vol_in (auxiliary unit has trouble doing this)
            ins_F_vol = self.F_vol_in
            for i in range(len(self.ins)):
                ins_F_vol -= (self.ins[i].ivol['H2'] +\
                              self.ins[i].ivol['CHG_catalyst'] +\
                              self.ins[i].ivol['HT_catalyst'] +\
                              self.ins[i].ivol['HC_catalyst'])
            # not include gas (e.g. H2)
            V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * _Pa_to_psi # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor

        if self.N:
            if self.V:
                D = (4*self.V/pi/length_to_diameter)**(1/3)
                D *= _m_to_ft # convert from m to ft
                L = D * length_to_diameter

                Design['Total volume'] = self.V*self.N
                Design['Single reactor volume'] = self.V
                Design['Number of reactors'] = self.N
            else:
                V_reactor = V_total/self.N
                D = (4*V_reactor/pi/length_to_diameter)**(1/3)
                D *= _m_to_ft # convert from m to ft
                L = D * length_to_diameter

                Design['Residence time'] = self.tau
                Design['Total volume'] = V_total
                Design['Single reactor volume'] = V_reactor
                Design['Number of reactors'] = self.N
        else:
            N = ceil(V_total/self._Vmax)
            if N == 0:
                V_reactor = 0
                D = 0
                L = 0
            else:
                V_reactor = V_total / N
                D = (4*V_reactor/pi/length_to_diameter)**(1/3)
                D *= _m_to_ft # convert from m to ft
                L = D * length_to_diameter

            Design['Residence time'] = self.tau
            Design['Total volume'] = V_total
            Design['Single reactor volume'] = V_reactor
            Design['Number of reactors'] = N

        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
             Design['Wall thickness'] *= wall_thickness_factor
             # Weight is proportional to wall thickness in PressureVessel design
             Design['Weight'] = round(Design['Weight']*wall_thickness_factor, 2)

        if self.vessel_material == 'Carbon steel':
            self.construction = (
                Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
                )
            self.construction[0].quantity = Design['Weight']*_lb_to_kg
        
        if self.vessel_material == 'Stainless steel 316':
            self.construction = (
                Construction('stainless_steel', linked_unit=self, item='Stainless_steel', quantity_unit='kg'),
                )
            self.construction[0].quantity = Design['Weight']*_lb_to_kg

    def _cost(self):
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs

        if Design['Total volume'] == 0:
            for i, j in purchase_costs.items():
                purchase_costs[i] = 0

        else:
            purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in purchase_costs.items():
                purchase_costs[i] *= Design['Number of reactors']

            self.power_utility(self.kW_per_m3*Design['Total volume'])

    @property
    def BM(self):
        vessel_type = self.vessel_type
        if not vessel_type:
            raise AttributeError('Vessel type not defined')
        elif vessel_type == 'Vertical':
            return self.BM_vertical
        elif vessel_type == 'Horizontal':
            return self.BM_horizontal
        else:
            raise RuntimeError('Invalid vessel type')

    @property
    def kW_per_m3(self):
        G = self.mixing_intensity
        if G is None:
            return self._kW_per_m3
        else:
            mixture = Stream()
            mixture.mix_from(self.ins)
            kW_per_m3 = mixture.mu*(G**2)/1e3
            return kW_per_m3
        
    @kW_per_m3.setter
    def kW_per_m3(self, i):
        if self.mixing_intensity and i is not None:
            raise AttributeError('mixing_intensity is provided, kw_per_m3 will be calculated.')
        else:
            self._kW_per_m3 = i

# =============================================================================
# Sludge Lab
# =============================================================================

class WWTP(SanUnit):
    '''
    WWTP is a fake unit that can set up sludge biochemical compositions
    and calculate sludge elemental compositions.
    Parameters
    ----------
    ins: stream
        ww.
    outs: stream
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
    References
    ----------
    .. [1] Metcalf and Eddy, Incorporated. 1991. Wastewater Engineering:
        Treatment Disposal and Reuse. New York: McGraw-Hill.
    .. [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', 
                 ww_2_dry_sludge=0.94, # [1]
                 sludge_moisture=0.99, sludge_dw_ash=0.257, 
                 sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, 
                 lipid_2_C=0.750, protein_2_C=0.545,
                 carbo_2_C=0.400, C_2_H=0.143,
                 protein_2_N=0.159, N_2_P=0.393,
                 operation_hour=yearly_operation_hour,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ww_2_dry_sludge = ww_2_dry_sludge
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_ash = sludge_dw_ash
        self.sludge_afdw_lipid = sludge_afdw_lipid
        self.sludge_afdw_protein = sludge_afdw_protein
        self.sludge_afdw_carbo = 1 - sludge_afdw_protein - sludge_afdw_lipid
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.C_2_H = C_2_H
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.operation_hour = yearly_operation_hour

    _N_ins = 1
    _N_outs = 2
    
    def _run(self):
        
        ww = self.ins[0]
        sludge, treated = self.outs
        
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
       return self.sludge_C_ratio*self.C_2_H
   
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
        
    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
# KOdrum
# =============================================================================

class KOdrum(Reactor):
    '''
    Konckout drum is a HTL auxiliary unit.
    References
    ----------
    .. [1] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
        Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
        NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
        https://doi.org/10.2172/1111191.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 P=3049.7*6894.76, tau=0, V_wf=0,
                 length_to_diameter=2, N=4, V=None,
                 auxiliary=True,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
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

    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def _run(self):
        pass
        
    def _design(self):
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)

# =============================================================================
# HTL (ignore three phase separator for now, ask Yalin)
# =============================================================================

@cost(basis='Treatment capacity', ID='Solids filter oil/water separator', units='lb/h',
      cost=3945523, S=1219765,
      CE=CEPCI[2011], n=0.68, BM=1.9)
# separator

class HTL(Reactor):
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and biochar
    under elevated temperature (350°C) and pressure. The products percentage
    (wt%) can be evaluated using revised MCA model (Li et al., 2017,
    Leow et al., 2018) with known sludge composition (protein%, lipid%,
    and carbohydrate%, all afdw%).
    Notice that for HTL we just calculate each phases' total mass (except gas)
    and calculate C, N, and P amount in each phase as properties. We don't
    specify components for oil/char since we want to use MCA model to calculate
    C and N amount and it is not necessary to calculate every possible
    components since they will be treated in HT/AcidEx anyway. We also don't
    specify components for aqueous since we want to calculate aqueous C, N, and
    P based on mass balance closure. But later for CHG, HT, and HC, we specify
    each components (except aqueous phase) for the application of flash,
    distillation column, and CHP units.
    Parameters
    ----------
    ins: stream
        dewatered_sludge.
    outs: stream
        biochar, HTLaqueous, biocrude, offgas.
    lipid_2_biocrude: float
        Lipid to biocrude factor.
    protein_2_biocrude: float
        Protein to biocrude factor.
    carbo_2_biocrude: float
        Carbohydrate to biocrude factor.
    protein_2_gas: float
        Protein to gas factor.
    carbo_2_gas: float
        Carbohydrate to gas factor.
    biocrude_C_slope: float
        Biocrude carbon content slope.
    biocrude_C_intercept: float
        Biocrude carbon content intercept.
    biocrude_N_slope: float
        Biocrude nitrogen content slope.
    biocrude_H_slope: float
        Biocrude hydrogen content slope.
    biocrude_H_intercept: float
        Biocrude hydrogen content intercept.
    HTLaqueous_C_slope: float
        HTLaqueous carbon content slope.
    TOC_TC: float   
        HTL TOC/TC.
    biochar_C_slope: float
        Biochar carbon content slope.
    biocrude_moisture_content: float
        Biocrude moisture content.
    biochar_P_recovery_ratio: float
        Biochar phosphorus to total phosphorus ratio.
    gas_composition: dict
        HTL offgas compositions.
    biochar_pre: float
        Biochar pressure, [Pa].
    HTLaqueous_pre: float
        HTL aqueous phase pressure, [Pa].
    biocrude_pre: float
        Biocrude pressure, [Pa].
    offgas_pre: float
        Offgas pressure, [Pa].
    eff_T: float
        HTL effluent temperature, [K].
    CAPEX_factor: float
        Factor used to adjust CAPEX.
    References
    ----------
    .. [1] Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
        Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
        Liquefaction Products from Feedstock Biochemical Composition.
        Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    .. [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    .. [3] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J.
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products.
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727.
        https://doi.org/10.1021/acs.est.8b04035.
    .. [4] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [5] Matayeva, A.; Rasmussen, S. R.; Biller, P. Distribution of Nutrients and
        Phosphorus Recovery in Hydrothermal Liquefaction of Waste Streams.
        BiomassBioenergy 2022, 156, 106323.
        https://doi.org/10.1016/j.biombioe.2021.106323.
    .. [6] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels
        via Liquefaction - Hydrothermal Liquefaction Reactor Design:
        April 5, 2013; NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462,
        1111191. https://doi.org/10.2172/1111191.
    '''
    
    auxiliary_unit_names=('heat_exchanger','kodrum')

    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 lipid_2_biocrude=0.846, # [1]
                 protein_2_biocrude=0.445, # [1]
                 carbo_2_biocrude=0.205, # [1]
                 protein_2_gas=0.074, # [1]
                 carbo_2_gas=0.418, # [1]
                 biocrude_C_slope=-8.37, # [2]
                 biocrude_C_intercept=68.55, # [2]
                 biocrude_N_slope=0.133, # [2]
                 biocrude_H_slope=-2.61, # [2]
                 biocrude_H_intercept=8.20, # [2]
                 HTLaqueous_C_slope=478, # [2]
                 TOC_TC=0.764, # [3]
                 biochar_C_slope = 1.75, # [2]
                 biocrude_moisture_content=0.063, # [4]
                 biochar_P_recovery_ratio=0.86, # [5]
                 gas_composition={'CH4':0.050, 'C2H6':0.032,
                                  'CO2':0.918}, # [4]
                 biochar_pre=3029.7*6894.76, # [4]
                 HTLaqueous_pre=30*6894.76, # [4]
                 biocrude_pre=30*6894.76, # [4]
                 offgas_pre=30*6894.76, # [4]
                 eff_T=60+273.15, # [4]
                 P=None, tau=15/60, V_wf=0.3,
                 length_to_diameter=2, N=4, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lipid_2_biocrude = lipid_2_biocrude
        self.protein_2_biocrude = protein_2_biocrude
        self.carbo_2_biocrude = carbo_2_biocrude
        self.protein_2_gas = protein_2_gas
        self.carbo_2_gas = carbo_2_gas
        self.biocrude_C_slope = biochar_C_slope
        self.biocrude_C_intercept = biocrude_C_intercept
        self.biocrude_N_slope = biocrude_N_slope
        self.biocrude_H_slope = biocrude_H_slope
        self.biocrude_H_intercept = biocrude_H_intercept
        self.HTLaqueous_C_slope = HTLaqueous_C_slope
        self.TOC_TC = TOC_TC
        self.biochar_C_slope = biochar_C_slope
        self.biocrude_moisture_content = biocrude_moisture_content
        self.biochar_P_recovery_ratio = biochar_P_recovery_ratio
        self.gas_composition = gas_composition
        self.biochar_pre = biochar_pre
        self.HTLaqueous_pre = HTLaqueous_pre
        self.biocrude_pre = biocrude_pre
        self.offgas_pre = offgas_pre
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HTLHX(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=eff_T)
        self.kodrum = KOdrum(ID=f'.{ID}_KOdrum')
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
        self.CAPEX_factor = CAPEX_factor

    _N_ins = 1
    _N_outs = 4
    _units= {'Treatment capacity': 'lb/h',
             'Solid filter and separator weight': 'lb'}
    
    def _run(self):
        
        dewatered_sludge = self.ins[0]
        biochar, HTLaqueous, biocrude, offgas = self.outs
        
        self.WWTP = self.ins[0]._source.ins[0]._source.ins[0].\
                         _source.ins[0]._source
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid'] +\
                                dewatered_sludge.imass['Sludge_protein'] +\
                                dewatered_sludge.imass['Sludge_carbo']
        # just use afdw in revised MCA model, other places use dw
        
        afdw_lipid_ratio = self.WWTP.sludge_afdw_lipid
        afdw_protein_ratio = self.WWTP.sludge_afdw_protein
        afdw_carbo_ratio = self.WWTP.sludge_afdw_carbo

        # the following calculations are based on revised MCA model
        biochar.imass['Biochar'] = 0.377*afdw_carbo_ratio*dewatered_sludge_afdw
        
        HTLaqueous.imass['HTLaqueous'] = (0.481*afdw_protein_ratio +\
                                          0.154*afdw_lipid_ratio)*\
                                          dewatered_sludge_afdw
        # HTLaqueous is TDS in aqueous phase
        # 0.377, 0.481, and 0.154 don't have uncertainties becasue they are calculated values
         
        gas_mass = (self.protein_2_gas*afdw_protein_ratio + self.carbo_2_gas*afdw_carbo_ratio)*\
                       dewatered_sludge_afdw
                       
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass*ratio
            
        biocrude.imass['Biocrude'] = (self.protein_2_biocrude*afdw_protein_ratio +\
                                      self.lipid_2_biocrude*afdw_lipid_ratio +\
                                      self.carbo_2_biocrude*afdw_carbo_ratio)*\
                                      dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
                                
        HTLaqueous.imass['H2O'] = dewatered_sludge.F_mass - biochar.F_mass -\
                                  biocrude.F_mass - gas_mass - HTLaqueous.imass['HTLaqueous']
        # assume ash (all soluble based on Jones) goes to water
        
        biochar.phase = 's'
        offgas.phase = 'g'
        
        biochar.P = self.biochar_pre
        HTLaqueous.P = self.HTLaqueous_pre
        biocrude.P = self.biocrude_pre
        offgas.P = self.offgas_pre
        
        for stream in self.outs: stream.T = self.heat_exchanger.T

    @property
    def biocrude_C_ratio(self):
        return (self.WWTP.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100 # [2]
    
    @property
    def biocrude_H_ratio(self):
        return (self.WWTP.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100 # [2]

    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.WWTP.sludge_dw_protein # [2]
    
    @property
    def biocrude_C(self):
        return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.WWTP.sludge_C)


    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.WWTP.sludge_dw_protein*100/1000000/self.TOC_TC,
                   self.WWTP.sludge_C - self.biocrude_C)

    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio

    @property
    def biocrude_N(self):
        return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.WWTP.sludge_N)
    
    @property
    def biocrude_HHV(self):
        return 30.74 - 8.52*self.WWTP.AOSc +\
               0.024*self.WWTP.sludge_dw_protein # [2]
               
    @property
    def energy_recovery(self):
        return self.biocrude_HHV*self.outs[2].imass['Biocrude']/\
               (self.WWTP.outs[0].F_mass -\
               self.WWTP.outs[0].imass['H2O'])/self.WWTP.sludge_HHV # [2]
        
    @property
    def offgas_C(self):
        carbon = 0
        for name in self.gas_composition.keys():
            carbon += self.outs[3].imass[name]*cmps[name].i_C
        return min(carbon, self.WWTP.sludge_C - self.biocrude_C - self.HTLaqueous_C)    
        
    @property
    def biochar_C_ratio(self):
        return min(self.biochar_C_slope*self.WWTP.sludge_dw_carbo, 0.65) # [2]

    @property
    def biochar_C(self):
        return min(self.outs[0].F_mass*self.biochar_C_ratio, self.WWTP.sludge_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)

    @property
    def biochar_P(self):
        return min(self.WWTP.sludge_P*self.biochar_P_recovery_ratio, self.outs[0].F_mass)

    @property
    def HTLaqueous_N(self):
        return self.WWTP.sludge_N - self.biocrude_N
        
    @property
    def HTLaqueous_P(self):
        return self.WWTP.sludge_P*(1 - self.biochar_P_recovery_ratio)

    def _design(self):
        
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.copy_like(hx_ins0)
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx_outs0.T = hx.T
        hx_ins0.P = hx_outs0.P = self.outs[0].P # cooling before depressurized, heating after pressurized
        # in other words, both heating and cooling are performed under relatively high pressure
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

        self.P = self.ins[0].P
        Reactor._design(self)
        Design['Solid filter and separator weight'] = 0.2*Design['Weight'] # assume stainless steel
        # based on [6], case D design table, the purchase price of solid filter and separator to
        # the purchase price of HTL reactor is around 0.2, therefore, assume the weight of solid filter
        # and separator is 0.2*HTL weight
        self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg
        
        self.kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
        # in [6], when knockout drum influent is 1225236 lb/hr, single knockout
        # drum volume is 4230 gal
        
        self.kodrum.simulate()
        
    def _cost(self):
        Reactor._cost(self)
        self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor
        
# =============================================================================
# Acid Extraction
# =============================================================================

class AcidExtraction(Reactor):
    '''
    H2SO4 is added to biochar from HTL to extract phosphorus. 
    Parameters
    ----------
    ins: stream
        biochar, acid.
    outs: stream
        residual, extracted.
    acid_vol: float
        0.5 M H2SO4 to biochar ratio: mL/g.
    P_acid_recovery_ratio: float
        The ratio of phosphorus that can be extracted.
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
        Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
        2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
        https://doi.org/10.2172/1855835.
    '''
    
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', acid_vol=7, P_acid_recovery_ratio=0.8,
                 P=None, tau=2, V_wf=0.8, # tau: [1]
                 length_to_diameter=2, N=1, V=10, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304', # acid condition
                 vessel_type='Vertical',
                 **kwargs):
        
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

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        biochar, acid = self.ins
        residual, extracted = self.outs
        
        self.HTL = self.ins[0]._source
        
        if biochar.F_mass <= 0:
            pass
        else:
            if self.HTL.biochar_P <= 0:
                residual.copy_like(biochar)
            else: 
                acid.imass['H2SO4'] = biochar.F_mass*self.acid_vol*0.5*98.079/1000
                #0.5 M H2SO4 acid_vol (10 mL/1 g) Biochar
                acid.imass['H2O'] = biochar.F_mass*self.acid_vol*1.05 -\
                                    acid.imass['H2SO4']
                # 0.5 M H2SO4 density: 1.05 kg/L 
                # https://www.fishersci.com/shop/products/sulfuric-acid-1n-0-5m-
                # standard-solution-thermo-scientific/AC124240010 (accessed 10-6-2022)
                
                residual.imass['Residual'] = biochar.F_mass - self.ins[0]._source.\
                                             biochar_P*self.P_acid_recovery_ratio
                
                extracted.copy_like(acid)
                extracted.imass['P'] = biochar.F_mass - residual.F_mass
                # assume just P can be extracted
                
                residual.phase = 's'
                
                residual.T = extracted.T = biochar.T
                residual.P = biochar.P
                # H2SO4 reacts with biochar to release heat and temperature will be
                # increased mixture's temperature
            
    @property
    def residual_C(self):
        return self.ins[0]._source.biochar_C
    
    @property
    def residual_P(self):
        return self.ins[0]._source.biochar_P - self.outs[1].imass['P']
        
    def _design(self):
        self.N = ceil(self.HTL.WWTP.ins[0].F_vol/788.627455/self.V)
        # 1/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 10 m3)
        self.P = self.ins[1].P
        Reactor._design(self)
        
    def _cost(self):
        Reactor._cost(self)
    
# =============================================================================
# HTL mixer
# =============================================================================

class HTLmixer(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the mixture of HTL
    aqueous and AcidEx effluent.
    Parameters
    ----------
    ins: stream
        HTLaqueous, extracted
    outs: stream
        mixture
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 2
    _N_outs = 1
        
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
        # assume HTLaqueous pH = 9 [1] (9.08 ± 0.30)
        # extracted pH = 0 (0.5 M H2SO4)
        # since HTLaqueous pH is near to neutral
        # assume pH is dominant by extracted and will be calculate based on dilution
        dilution_factor = self.F_mass_in/self.ins[1].F_mass if self.ins[1].imass['P'] != 0 else 1
        hydrogen_ion_conc = 10**0/dilution_factor
        return -log(hydrogen_ion_conc, 10)
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
# =============================================================================
# HTLsplitter
# =============================================================================

class HTLsplitter(SanUnit):
    '''
    A fake unit that calculates influent based on effluents.
    Parameters
    ----------
    ins: stream
    flow_in
    outs: stream
    flow_out_1, flow_out_2
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo,init_with)

    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        flow_in = self.ins[0]
        flow_out_1, flow_out_2 = self.outs
        
        flow_in.mix_from((flow_out_1, flow_out_2))
    
    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
# Struvite Precipitation
# =============================================================================

class StruvitePrecipitation(Reactor):
    '''
    Extracted and HTL aqueous are mixed together before adding MgCl2 for struvite precipitation.
    If mol(N)<mol(P), add NH4Cl to mol(N):mol(P)=1:1.
    Parameters
    ----------
    ins: stream
        mixture, supply_MgCl2, supply_NH4Cl, base.
    outs: stream
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
                 vessel_type='Vertical',
                 **kwargs):
        
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
        
    _N_ins = 4
    _N_outs = 2
        
    def _run(self):
        
        mixture, supply_MgCl2, supply_NH4Cl, base = self.ins
        struvite, effluent = self.outs
        
        self.HTLmixer = self.ins[0]._source
        
        if self.HTLmixer.outs[0].imass['P'] == 0:
            effluent.copy_like(mixture)
        else:
            old_pH = self.HTLmixer.pH
            neutral_OH_mol = 10**(-old_pH)*self.ins[0].F_mass # ignore solid volume
            to_target_pH = 10**(self.target_pH - 14)*self.ins[0].F_mass # ignore solid volume
            total_OH = neutral_OH_mol + to_target_pH # unit: mol/h
            base.imass['MgO'] = total_OH/2 * 40.3044/1000
            
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
        self.N = ceil(self.HTLmixer.ins[1]._source.HTL.WWTP.ins[0].F_vol*2/788.627455/self.V)
        # 2/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 20 m3)
        self.P = self.ins[0].P
        Reactor._design(self)
        
    def _cost(self):
        Reactor._cost(self)
    
# =============================================================================
# CHG
# =============================================================================

@cost(basis='Treatment capacity', ID='Hydrocyclone', units='lb/h',
      cost=5000000, S=968859,
      CE=CEPCI[2009], n=0.65, BM=2.1)
# hydrocyclone

class CHG(Reactor):
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    Parameters
    ----------
    ins: stream
        chg_in, catalyst_in.
    outs: stream
        chg_out, catalyst_out.
    pump_pressure: float
        CHG influent pressure, [Pa].
    heat_temp: float
        CHG influent temperature, [K].
    cool_temp: float
        CHG effluent temperature, [K].
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    gas_composition: dict
        CHG gas composition.
    gas_C_2_total_C: dict
        CHG gas carbon content to feed carbon content.
    CAPEX_factor: float
        Factor used to adjust CAPEX.
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [2] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.;
        Beckham, G. T.; Humbird, D.; Thompson, D. N.; Roni, M. S. Process
        Design and Economics for the Conversion of Lignocellulosic Biomass
        to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case
        Update; Biochemical Deconstruction and Conversion of Biomass to Fuels
        and Products via Integrated Biorefinery Pathways; NREL/TP--5100-71949,
        1483234; 2018; p NREL/TP--5100-71949, 1483234.
        https://doi.org/10.2172/1483234.
    .. [3] Elliott, D. C.; Neuenschwander, G. G.; Hart, T. R.; Rotness, L. J.;
        Zacher, A. H.; Santosa, D. M.; Valkenburg, C.; Jones, S. B.;
        Rahardjo, S. A. T. Catalytic Hydrothermal Gasification of Lignin-Rich
        Biorefinery Residues and Algae Final Report. 87.
    '''

    auxiliary_unit_names=('pump','heat_ex_heating','heat_ex_cooling')
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Sulfur guard': 2.0}
    _units= {'Treatment capacity': 'lb/h',
             'Hydrocyclone weight': 'lb'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 pump_pressure=3089.7*6894.76,
                 heat_temp=350+273.15,
                 cool_temp=60+273.15,
                 WHSV=3.56,
                 catalyst_lifetime=1*yearly_operation_hour, # 1 year [1]
                 gas_composition={'CH4':0.527,
                                  'CO2':0.432,
                                  'C2H6':0.011,
                                  'C3H8':0.030,
                                  'H2':0.0001}, # [1]
                 gas_C_2_total_C=0.598, # [1]
                 P=None, tau=20/60, void_fraction=0.5, # [2, 3]
                 length_to_diameter=2, N=6, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
        self.pump_pressure = pump_pressure
        self.heat_temp = heat_temp
        self.cool_temp = cool_temp
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.gas_composition = gas_composition
        self.gas_C_2_total_C = gas_C_2_total_C
        pump_in = bst.Stream(f'{ID}_pump_in')
        pump_out = bst.Stream(f'{ID}_pump_out')
        self.pump = HTLpump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=pump_pressure)
        hx_ht_in = bst.Stream(f'{ID}_hx_ht_in')
        hx_ht_out = bst.Stream(f'{ID}_hx_ht_out')
        self.heat_ex_heating = HTLHX(ID=f'.{ID}_hx_ht', ins=hx_ht_in, outs=hx_ht_out, T=heat_temp)
        hx_cl_in = bst.Stream(f'{ID}_hx_cl_in')
        hx_cl_out = bst.Stream(f'{ID}_hx_cl_out')
        self.heat_ex_cooling = HTLHX(ID=f'.{ID}_hx_cl', ins=hx_cl_in, outs=hx_cl_out, T=cool_temp)
        self.P = P
        self.tau = tau
        self.V_wf = void_fraction
        # no headspace, gases produced will be vented, so V_wf = void fraction [2, 3]
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.CAPEX_factor = CAPEX_factor
        
    _N_ins = 2
    _N_outs = 2
    _units= {'Treatment capacity': 'lb/h'} # hydrocyclone
        
    def _run(self):
        
        chg_in, catalyst_in = self.ins
        chg_out, catalyst_out = self.outs
        
        catalyst_in.imass['CHG_catalyst'] = chg_in.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
            
        gas_C_ratio = 0
        for name, ratio in self.gas_composition.items():
            gas_C_ratio += ratio*cmps[name].i_C
            
        gas_mass = chg_in.imass['C']*self.gas_C_2_total_C/gas_C_ratio
        
        for name,ratio in self.gas_composition.items():
            chg_out.imass[name] = gas_mass*ratio
                
        chg_out.imass['H2O'] = chg_in.F_mass - gas_mass
        # all C, N, and P are accounted in H2O here, but will be calculated as properties.
                
        chg_out.T = self.cool_temp
        chg_out.P = self.pump_pressure
        
    @property
    def CHGout_C(self):
        # not include carbon in gas phase
        return self.ins[0].imass['C']*(1 - self.gas_C_2_total_C)
    
    @property
    def CHGout_N(self):
        return self.ins[0].imass['N']
    
    @property
    def CHGout_P(self):
        return self.ins[0].imass['P']
        
    def _design(self):
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx_ht = self.heat_ex_heating
        hx_ht_ins0, hx_ht_outs0 = hx_ht.ins[0], hx_ht.outs[0]
        hx_ht_ins0.copy_like(self.ins[0])
        hx_ht_outs0.copy_like(hx_ht_ins0)
        hx_ht_ins0.T = self.ins[0].T
        hx_ht_outs0.T = hx_ht.T
        hx_ht_ins0.P = hx_ht_outs0.P = pump.P
        hx_ht.simulate_as_auxiliary_exchanger(ins=hx_ht.ins, outs=hx_ht.outs)
            
        hx_cl = self.heat_ex_cooling
        hx_cl_ins0, hx_cl_outs0 = hx_cl.ins[0], hx_cl.outs[0]
        hx_cl_ins0.copy_like(self.outs[0])
        hx_cl_outs0.copy_like(hx_cl_ins0)
        hx_cl_ins0.T = hx_ht.T
        hx_cl_outs0.T = hx_cl.T
        hx_cl_ins0.P = hx_cl_outs0.P = self.outs[0].P
        hx_cl.simulate_as_auxiliary_exchanger(ins=hx_cl.ins, outs=hx_cl.outs)

        self.P = self.pump_pressure
        Reactor._design(self)
        Design['Hydrocyclone weight'] = 0.3*Design['Weight'] # assume stainless steel
        # based on [1], page 54, the purchase price of hydrocyclone to the purchase price of CHG
        # reactor is around 0.3, therefore, assume the weight of hydrocyclone is 0.3*CHG weight
        self.construction[0].quantity += Design['Hydrocyclone weight']*_lb_to_kg
    
    def _cost(self):
        Reactor._cost(self)
        purchase_costs = self.baseline_purchase_costs
        current_cost = 0 # cost w/o sulfur guard
        for item in purchase_costs.keys():
            current_cost += purchase_costs[item]
        purchase_costs['Sulfur guard'] = current_cost*0.05
        self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor
    
# =============================================================================
# Membrane Distillation
# =============================================================================

class MembraneDistillation(SanUnit):
    '''
    Membrane distillation recovers nitrogen as ammonia sulfate based on vapor
    pressure difference across the hydrophobic membrane. Ignore water flux across
    membrane since it will not affect system performance (either TEA or LCA).
    Parameters
    ----------
    ins: stream
        influent, acid, base, mem_in.
    outs: stream
        ammoniumsulfate, ww, mem_out.
    influent_pH: float
        Influent pH.
    target_pH: float
        Target pH for membrane distillation.
    N_S_ratio: float
        mol(N) to mol(S) ratio.
    m2_2_m3: float
        m2 to m3 factor, 1/specific surface area, [m3/m2].
    Dm: float
        NH3 molecular diffusity in air, [m2/s]. 
    porosity: float
        Membrane porosity.
    thickness: float
        Membrane thickness, [m].
    tortuosity: float
        Membrane tortuosity.
    Henry: float
        NH3 Henry constant, [atm*m3/mol].
    Ka: float
        Overall mass transfer coefficient, [m/s].
    capacity: float
        Membrane treatement capacity (permeate flux), [kg/m2/h].
    membrane_price: float
        Membrane price, [$/kg] ([$/m2]).
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    .. [2] Doran, P. M. Chapter 11 - Unit Operations. In Bioprocess Engineering
        Principles (Second Edition); Doran, P. M., Ed.; Academic Press: London,
        2013; pp 445–595. https://doi.org/10.1016/B978-0-12-220851-5.00011-3.
    .. [3] Spiller, L. L. Determination of Ammonia/Air Diffusion Coefficient Using
        Nafion Lined Tube. Analytical Letters 1989, 22 (11–12), 2561–2573.
        https://doi.org/10.1080/00032718908052375.
    .. [4] Scheepers, D. M.; Tahir, A. J.; Brunner, C.; Guillen-Burrieza, E.
        Vacuum Membrane Distillation Multi-Component Numerical Model for Ammonia
        Recovery from Liquid Streams. Journal of Membrane Science
        2020, 614, 118399. https://doi.org/10.1016/j.memsci.2020.118399.
    .. [5] Ding, Z.; Liu, L.; Li, Z.; Ma, R.; Yang, Z. Experimental Study of Ammonia
        Removal from Water by Membrane Distillation (MD): The Comparison of Three
        Configurations. Journal of Membrane Science 2006, 286 (1), 93–103.
        https://doi.org/10.1016/j.memsci.2006.09.015.
    .. [6] Al-Obaidani, S.; Curcio, E.; Macedonio, F.; Di Profio, G.; Al-Hinai, H.;
        Drioli, E. Potential of Membrane Distillation in Seawater Desalination:
        Thermal Efficiency, Sensitivity Study and Cost Estimation.
        Journal of Membrane Science 2008, 323 (1), 85–98.
        https://doi.org/10.1016/j.memsci.2008.06.006.
    .. [7] Kogler, A.; Farmer, M.; Simon, J. A.; Tilmans, S.; Wells, G. F.;
        Tarpeh, W. A. Systematic Evaluation of Emerging Wastewater Nutrient Removal
        and Recovery Technologies to Inform Practice and Advance Resource
        Efficiency. ACS EST Eng. 2021, 1 (4), 662–684.
        https://doi.org/10.1021/acsestengg.0c00253.
    '''
    
    _F_BM_default = {'Membrane': 1}
    
    _units = {'Area': 'm2',
              'Total volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 influent_pH=8.16, # CHG effluent pH: 8.16 ± 0.25 [1]
                 target_pH=10,
                 N_S_ratio=2,
                 # S is excess since not all N can be transfered to form ammonia sulfate
                 # for now, assume N_S_ratio = 2 is ok
                 m2_2_m3=1/1200, # specific surface area, for hollow fiber membrane [2]
                 Dm=2.28*10**(-5), # (2.28 ± 0.12)*10^-5 m^2/s NH3 molecular diffusity in air [3]
                 # (underestimate, this value may be at 15 or 25 C, our feed is 60 C, should be higher)
                 porosity=0.9, # [4]
                 thickness=7*10**(-5), # m [4]
                 tortuosity=1.2, # [4]
                 Henry=1.61*10**(-5), # atm*m3/mol
                 # https://webwiser.nlm.nih.gov/substance?substanceId=315&identifier=\
                 # Ammonia&identifierType=name&menuItemId=81&catId=120#:~:text=The%20\
                 # Henry's%20Law%20constant%20for,m%2Fmole(2). (accessed 11-11-2022)
                 Ka=1.75*10**(-5), # overall mass transfer coefficient 1.2~2.3*10^-5 m/s [5]
                 capacity=6.01, # kg/m2/h [6]
                 # permeate flux, similar values can be found in many other papers ([176], [222] ,[223] in [7])
                 # [177], [223] in [7] show high nitrogen recovery ratio (>85% under optimal conditions)
                 membrane_price=93.29, # $90/m2 2008 [6]
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.influent_pH = influent_pH
        self.target_pH = target_pH
        self.N_S_ratio = N_S_ratio
        self.m2_2_m3 = m2_2_m3
        self.Dm = Dm
        self.porosity = porosity
        self.thickness = thickness
        self.tortuosity = tortuosity
        self.Henry = Henry
        self.Ka = Ka
        self.capacity = capacity
        self.membrane_price = membrane_price

    _N_ins = 4
    _N_outs = 4
    
    def _run(self):
        
        influent, acid, base, mem_in = self.ins
        ammoniumsulfate, ww, mem_out, ammoniumsulfatesolution = self.outs
        
        self.CHG = self.ins[0]._source.ins[0]._source.ins[0]._source
        
        if self.CHG.CHGout_N == 0:
            ww.copy_like(influent)
            self.membrane_area = self.F_vol_in*1000/self.capacity
        else:
            NaOH_conc = 10**(self.target_pH - 14) - 10**(self.influent_pH - 14)
            NaOH_mol = NaOH_conc*self.ins[0].F_mass
            base.imass['NaOH'] = NaOH_mol*39.997/1000
            
            acid.imass['H2SO4'] = self.CHG.CHGout_N/14.0067/self.N_S_ratio*98.079
            acid.imass['H2O'] = acid.imass['H2SO4']*1000/98.079/0.5*1.05 -\
                                acid.imass['H2SO4']
            
            pKa = 9.26 # ammonia pKa
            ammonia_to_ammonium = 10**(-pKa)/10**(-self.target_pH)
            ammonia = self.CHG.CHGout_N*ammonia_to_ammonium/(1 +\
                       ammonia_to_ammonium)*17.031/14.0067
            others = influent.F_mass - ammonia
            
            self.membrane_area = self.F_vol_in*1000/self.capacity
            N2_in_air = self.membrane_area*self.m2_2_m3*self.porosity*0.79*1.204
            O2_in_air = self.membrane_area*self.m2_2_m3*self.porosity*0.21*1.204
            # N2:O2 = 0.79:0.21 in the air, air density is 1.204 kg/m3
            # https://en.wikipedia.org/wiki/Density_of_air#:~:text=It%20also%20\
            # changes%20with%20variation,International%20Standard%20Atmosphere%2\
            # 0(ISA). (accessed 11-14-2022)
            
            imass = indexer.MassFlowIndexer(l=[('H2O', others),
                                               ('NH3', ammonia),
                                               ('N2', 0),
                                               ('O2', 0)],
                                            g=[('H2O', 0),
                                               ('NH3', 0), 
                                               ('N2', N2_in_air),
                                               ('O2', O2_in_air)])
            # N2 amount will be changed based on design, maybe also add O2 (N2:O2 = 4:1)
            
            vle = equilibrium.VLE(imass)
            vle(T=influent.T, P=influent.P)
            X_NH3_f_m = vle.imol['g','NH3']/(vle.imol['g','H2O'] + vle.imol['g','NH3'])
            X_NH3_f = vle.imol['l','NH3']/(vle.imol['l','H2O'] + vle.imol['l','NH3'])
    
            km = self.Dm*self.porosity/self.tortuosity/self.thickness
            
            dimensionless_Henry = self.Henry/8.20575/(10**(-5))/influent.T # H' = H/RT
            # https://www.sciencedirect.com/topics/chemistry/henrys-law#:~:text=\
            # Values%20for%20Henry's%20law%20constants,gas%20constant%20(8.20575%\
            # 20%C3%97%2010 (accessed 11-11-2022)
            
            kf = 1/(1/self.Ka - 1/dimensionless_Henry/km*(1 + 10**(-pKa)*10**\
                 (-self.target_pH)/(10**(-14))))
            
            J = kf*ammonia/influent.F_mass*1000*log(X_NH3_f_m/X_NH3_f)*3600 # in kg/m2/h
            
            NH3_mass_flow = J*self.membrane_area
            
            ammonia_transfer_ratio = min(1, NH3_mass_flow/ammonia)
            
            ammoniumsulfate.imass['NH42SO4'] = ammonia*ammonia_transfer_ratio/34.062*132.14
            ammoniumsulfatesolution.imass['H2O'] = acid.imass['H2O']
            ammoniumsulfatesolution.imass['H2SO4'] = acid.imass['H2SO4'] +\
                                             ammoniumsulfate.imass['NH42SO4']/\
                                             132.14*28.0134 -\
                                             ammoniumsulfate.imass['NH42SO4']
                                            
            ww.copy_like(influent) # ww has the same T and P as influent
            
            ww.imass['N'] = self.CHG.CHGout_N*(1 - ammonia_to_ammonium/(1 +\
                             ammonia_to_ammonium)*ammonia_transfer_ratio)
                                                              
            ww.imass['C'] = self.CHG.CHGout_C
                             
            ww.imass['P'] = self.CHG.CHGout_P
                       
            ww.imass['H2O'] -= (ww.imass['C'] + ww.imass['N'] + ww.imass['P'])
            
            ww.imass['H2O'] += self.ins[2].F_mass
            
            ammoniumsulfate.T = ammoniumsulfatesolution.T = acid.T
            ammoniumsulfate.P = ammoniumsulfatesolution.P = acid.P
            # ammoniumsulfate has the same T and P as acid
            
            mem_in.imass['Membrane'] = 0.15*self.membrane_area/yearly_operation_hour # kg/hr (m2/hr)
            mem_out.copy_like(mem_in)
            # add membrane as streams to include 15% membrane replacement per year [6]
    
    @property
    def N_recovery_ratio(self):
        return 1 - self.outs[1].imass['N']/self.CHG.CHGout_N
        
    def _design(self):
        Design = self.design_results
        Design['Area'] = self.membrane_area
        Design['Total volume'] = Design['Area']*self.m2_2_m3
        
        self.construction = (
            Construction('membrane', linked_unit=self, item='RO', quantity_unit='m2'),
            )
        self.construction[0].quantity = Design['Area']
    
    def _cost(self):
        self.ins[3].price = self.membrane_price
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs
        purchase_costs['Membrane'] = Design['Area']*self.membrane_price

# =============================================================================
# HT
# =============================================================================

class HT(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    Parameters
    ----------
    ins: stream
        biocrude, hydrogen, catalyst_in.
    outs: stream
        ht_out, catalyst_out = self.outs.
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_biocrude: float
        Reacted H2 to biocrude mass ratio.
    hydrogen_excess: float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    hydrocarbon_ratio: float
        Mass ratio of produced hydrocarbon to the sum of biocrude and reacted H2.
    HTin_T: float
        HT influent temperature, [K].
    HTrxn_T: float
        HT effluent (after reaction) temperature, [K].
    HT_composition: dict
        HT effluent composition.
    CAPEX_factor: float
        Factor used to adjust CAPEX.
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [2] Towler, G.; Sinnott, R. Chapter 14 - Design of Pressure Vessels.
        In Chemical Engineering Design (Second Edition); Towler, G., Sinnott, R.,
        Eds.; Butterworth-Heinemann: Boston, 2013; pp 563–629.
        https://doi.org/10.1016/B978-0-08-096659-5.00014-6.
    '''
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    _units = {'Hydrogen': 'mmscfd'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=2*yearly_operation_hour, # 2 years [1]
                 hydrogen_P=1530*6894.76,
                 hydrogen_rxned_to_biocrude=0.046,
                 hydrogen_excess=3,
                 hydrocarbon_ratio=0.875, # 87.5 wt% of biocrude and reacted H2 [1]
                 # spreadsheet HT calculation
                 HTin_T=174+273.15,
                 HTrxn_T=402+273.15, # [1]
                 HT_composition={'CH4':0.02280, 'C2H6':0.02923,
                                 'C3H8':0.01650, 'C4H10':0.00870,
                                 'TWOMBUTAN':0.00408, 'NPENTAN':0.00678,
                                 'TWOMPENTA':0.00408, 'HEXANE':0.00401,
                                 'TWOMHEXAN':0.00408, 'HEPTANE':0.00401,
                                 'CC6METH':0.01020, 'PIPERDIN':0.00408,
                                 'TOLUENE':0.01013, 'THREEMHEPTA':0.01020,
                                 'OCTANE':0.01013, 'ETHCYC6':0.00408,
                                 'ETHYLBEN':0.02040, 'OXYLENE':0.01020,
                                 'C9H20':0.00408, 'PROCYC6':0.00408,
                                 'C3BENZ':0.01020, 'FOURMONAN':0,
                                 'C10H22':0.00240, 'C4BENZ':0.01223,
                                 # C10H22 was originally 0.00203, but it is not
                                 # good for distillation column, the excess amount
                                 # is substracted from HEXANE, HEPTANE, TOLUENE,
                                 # OCTANE, and C9H20, which were originally 0.00408,
                                 # 0.00408, 0.01020, 0.01020, and 0.00408
                                 'C11H24':0.02040, 'C10H12':0.02040,
                                 'C12H26':0.02040, 'OTTFNA':0.01020,
                                 'C6BENZ':0.02040, 'OTTFSN':0.02040,
                                 'C7BENZ':0.02040, 'C8BENZ':0.02040,
                                 'C10H16O4':0.01837, 'C15H32':0.06120,
                                 'C16H34':0.18360, 'C17H36':0.08160, 
                                 'C18H38':0.04080, 'C19H40':0.04080,
                                 'C20H42':0.10200, 'C21H44':0.04080,
                                 'TRICOSANE':0.04080, 'C24H38O4':0.00817,
                                 'C26H42O4':0.01020, 'C30H62':0.00203}, # [1]
                 # spreadsheet HT calculation
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=0.5, void_fraciton=0.4, # [2]
                 length_to_diameter=2, N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_biocrude = hydrogen_rxned_to_biocrude
        self.hydrogen_excess = hydrogen_excess
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HTin_T = HTin_T
        self.HTrxn_T = HTrxn_T
        self.HT_composition = HT_composition
        IC_in = bst.Stream(f'{ID}_IC_in')
        IC_out = bst.Stream(f'{ID}_IC_out')
        self.compressor = HTLcompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HTLHX(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.CAPEX_factor = CAPEX_factor

    _N_ins = 3
    _N_outs = 2
        
    def _run(self):
        
        biocrude, hydrogen, catalyst_in = self.ins
        ht_out, catalyst_out = self.outs
        
        self.HTL = self.ins[0]._source.ins[0]._source
        
        if self.HTL.biocrude_N == 0:
            remove = self.HT_composition['PIPERDIN']
            for chemical in self.HT_composition.keys():  
                self.HT_composition[chemical] /= (1-remove)
            self.HT_composition['PIPERDIN'] = 0
        
        catalyst_in.imass['HT_catalyst'] = biocrude.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen.imass['H2'] = biocrude.imass['Biocrude']*\
                               self.hydrogen_rxned_to_biocrude*self.hydrogen_excess
        hydrogen.phase = 'g'

        hydrocarbon_mass = biocrude.imass['Biocrude']*\
                           (1 + self.hydrogen_rxned_to_biocrude)*\
                           self.hydrocarbon_ratio
        for name, ratio in self.HT_composition.items():
            ht_out.imass[name] = hydrocarbon_mass*ratio
            
        ht_out.imass['H2'] = hydrogen.imass['H2'] -\
                             biocrude.imass['Biocrude']*\
                             self.hydrogen_rxned_to_biocrude
        
        ht_out.imass['H2O'] = biocrude.F_mass + hydrogen.F_mass -\
                              hydrocarbon_mass - ht_out.imass['H2']
        # use water to represent HT aqueous phase,
        # C and N can be calculated base on MB closure.
        
        ht_out.P = biocrude.P
        
        ht_out.T = self.HTrxn_T
        
        if self.HTaqueous_C < -0.1*self.HTL.WWTP.sludge_C:
            raise Exception('carbon mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance
        # should be no C in the aqueous phase, the calculation here is just for MB
        
        if self.HTaqueous_N < -0.1*self.HTL.WWTP.sludge_N:
            raise Exception('nitrogen mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance

        # possibility exist that more carbon is in biooil and gas than in
        # biocrude because we use the biooil/gas compositions to calculate
        # carbon. In this case, the C in HT aqueous phase will be negative.
        # It's OK if the mass balance is within +/- 10% of total carbon in 
        # sludge. Otherwise, an exception will be raised.
        
    @property
    def hydrocarbon_C(self):
        carbon = 0
        for name in self.HT_composition.keys():
            carbon += self.outs[0].imass[name]*cmps[name].i_C
        return carbon

    @property
    def hydrocarbon_N(self):
        nitrogen = 0
        for name in self.HT_composition.keys():
            nitrogen += self.outs[0].imass[name]*cmps[name].i_N
        return nitrogen

    @property
    def HTaqueous_C(self):
        return self.HTL.biocrude_C - self.hydrocarbon_C
    # should be no C in the aqueous phase, the calculation here is just for MB

    @property
    def HTaqueous_N(self):
        return self.HTL.biocrude_N - self.hydrocarbon_N

    def _design(self):
        Design = self.design_results
        Design['Hydrogen'] = self.ins[1].F_vol*_m3perh_to_mmscfd
        
        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(self.ins[1])
        IC_outs0.copy_like(self.ins[1])
        IC_outs0.P = IC.P = self.hydrogen_P
        IC.simulate()
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.ins)
        hx_outs0.copy_like(hx_ins0)
        hx_outs0.T = self.HTin_T
        hx_ins0.P = hx_outs0.P = min(IC_outs0.P, self.ins[0].P)
        # H2 and biocrude have the same pressure
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor

# =============================================================================
# HC
# =============================================================================

class HC(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    Parameters
    ----------
    ins: stream
        heavy_oil, hydrogen, catalyst_in.
    outs: stream
        hc_out, catalyst_out.
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_heavy_oil: float
        Reacted H2 to heavy oil mass ratio.
    hydrogen_excess: float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    hydrocarbon_ratio: float
        Mass ratio of produced hydrocarbon to the sum of heavy oil and reacted H2.
    HCin_T: float
        HC influent temperature, [K].
    HCrxn_T: float
        HC effluent (after reaction) temperature, [K].
    HC_composition: dict
        HC effluent composition.
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    _units = {'Hydrogen': 'mmscfd'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=5*yearly_operation_hour, # 5 years [1]
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_rxned_to_heavy_oil=0.01125,
                 hydrogen_excess=5.556,
                 hydrocarbon_ratio=1, # 100 wt% of heavy oil and reacted H2
                 # nearly all input heavy oils and H2 will be converted to
                 # products [1]
                 # spreadsheet HC calculation
                 HCin_T=394+273.15,
                 HCrxn_T=451+273.15,
                 HC_composition={'CO2':0.03880, 'CH4':0.00630,
                                 'CYCHEX':0.03714, 'HEXANE':0.01111,
                                 'HEPTANE':0.11474, 'OCTANE':0.08125,
                                 'C9H20':0.09086, 'C10H22':0.11756,
                                 'C11H24':0.16846, 'C12H26':0.13198,
                                 'C13H28':0.09302, 'C14H30':0.04643,
                                 'C15H32':0.03250, 'C16H34':0.01923,
                                 'C17H36':0.00431, 'C18H38':0.00099,
                                 'C19H40':0.00497, 'C20H42':0.00033},
                 #combine C20H42 and PHYTANE as C20H42
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=5, void_fraciton=0.4, # Towler
                 length_to_diameter=2, N=None, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_heavy_oil = hydrogen_rxned_to_heavy_oil
        self.hydrogen_excess = hydrogen_excess
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HCin_T = HCin_T
        self.HCrxn_T = HCrxn_T
        self.HC_composition = HC_composition
        IC_in = bst.Stream(f'{ID}_IC_in')
        IC_out = bst.Stream(f'{ID}_IC_out')
        self.compressor = HTLcompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HTLHX(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    _N_ins = 3
    _N_outs = 2
        
    def _run(self):
        
        heavy_oil, hydrogen, catalyst_in = self.ins
        hc_out, catalyst_out = self.outs
        
        catalyst_in.imass['HC_catalyst'] = heavy_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen.imass['H2'] = heavy_oil.F_mass*self.hydrogen_rxned_to_heavy_oil*self.hydrogen_excess
        hydrogen.phase = 'g'

        hydrocarbon_mass = heavy_oil.F_mass*(1 +\
                           self.hydrogen_rxned_to_heavy_oil)*\
                           self.hydrocarbon_ratio

        for name, ratio in self.HC_composition.items():
            hc_out.imass[name] = hydrocarbon_mass*ratio
        
        hc_out.imass['H2'] = hydrogen.imass['H2'] - heavy_oil.F_mass*\
                             self.hydrogen_rxned_to_heavy_oil
        
        hc_out.P = heavy_oil.P
        hc_out.T = self.HCrxn_T
        
        C_in = 0
        total_num = len(list(cmps))
        for num in range(total_num):
            C_in += heavy_oil.imass[str(list(cmps)[num])]*list(cmps)[num].i_C
            
        C_out = self.hydrocarbon_C
        
        if C_out < 0.95*C_in or C_out > 1.05*C_out :
            raise Exception('carbon mass balance is out of +/- 5% for HC')
        # make sure that carbon mass balance is within +/- 5%. Otherwise, an
        # exception will be raised.
        
    @property
    def hydrocarbon_C(self):
        carbon = 0
        for name in self.HC_composition.keys():
            carbon += self.outs[0].imass[name]*cmps[name].i_C
        return carbon

    def _design(self):
        
        Design = self.design_results
        Design['Hydrogen'] = self.ins[1].F_vol*_m3perh_to_mmscfd
        
        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(self.ins[1])
        IC_outs0.copy_like(self.ins[1])
        IC_outs0.P = IC.P = self.hydrogen_P
        IC.simulate()
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.ins)
        hx_outs0.copy_like(hx_ins0)
        hx_outs0.T = self.HCin_T
        hx_ins0.P = hx_outs0.P = min(IC_outs0.P, self.ins[0].P)
        # H2 and biocrude have the same pressure
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)

# =============================================================================
# WWmixer
# =============================================================================

class WWmixer(SanUnit):
    '''
    A fake unit that mix all wastewater streams and calculates C, N, P, and H2O
    amount.
    Parameters
    ----------
    ins: stream
        supernatant_1, supernatant_2, memdis_ww, ht_ww
    outs: stream
        mixture
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 3
    _N_outs = 1
        
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

    def _design(self):
        pass
    
    def _cost(self):
        pass
    
# =============================================================================
# PhaseChanger
# =============================================================================

class PhaseChanger(SanUnit):
    '''
    Correct phase.
    Parameters
    ----------
    ins: stream
        influent
    outs: stream
        effluent
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', phase='l',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.phase = phase

    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
        
    def _run(self):
        
        influent = self.ins[0]
        effluent = self.outs[0]
        
        if self.phase not in ('g','l','s'):
            raise Exception("phase must be 'g','l', or 's'")
        
        effluent.copy_like(influent)
        effluent.phase = self.phase

    def _design(self):
        pass
    
    def _cost(self):
        pass
    
# =============================================================================
# FuelMixer
# =============================================================================

class FuelMixer(SanUnit):
    '''
    Convert gasoline to diesel or diesel to gasoline based on LHV.
    Parameters
    ----------
    ins: stream
        gasoline, diesel
    outs: stream
        fuel
    target: str
        The target can only be 'gasoline' or 'diesel'.
    gasoline_price: float
        Gasoline price, [$/kg].
    diesel_price: float
        Diesel price, [$/kg].
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target='diesel',
                 gasoline_gal_2_kg=2.834894885,
                 diesel_gal_2_kg=3.220628346,
                 gasoline_price=0.9388,
                 diesel_price=0.9722,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target = target
        self.gasoline_gal_2_kg = gasoline_gal_2_kg
        self.diesel_gal_2_kg = diesel_gal_2_kg
        self.gasoline_price = gasoline_price
        self.diesel_price = diesel_price

    _N_ins = 2
    _N_outs = 1

    def _run(self):
        
        gasoline, diesel = self.ins
        fuel = self.outs[0]
        
        if self.target not in ('gasoline','diesel'):
            raise RuntimeError ("target must be either 'gasoline' or 'diesel'")
        
        gasoline_LHV_2_diesel_LHV = (gasoline.LHV/gasoline.F_mass)/(diesel.LHV/diesel.F_mass)
        # KJ/kg gasoline:KJ/kg diesel
        
        if self.target == 'gasoline':
            fuel.imass['Gasoline'] = gasoline.F_mass + diesel.F_mass/gasoline_LHV_2_diesel_LHV
            fuel.T = gasoline.T
            fuel.P = gasoline.P
            fuel.price = self.gasoline_price
        
        if self.target == 'diesel':
            fuel.imass['Diesel'] = diesel.F_mass + gasoline.F_mass*gasoline_LHV_2_diesel_LHV
            fuel.T = diesel.T
            fuel.P = fuel.P
            fuel.price = self.diesel_price
            
    def _design(self):
        pass
    
    def _cost(self):
        if self.target == 'gasoline':
            self.outs[0].price = self.gasoline_price
        if self.target == 'diesel':
            self.outs[0].price = self.diesel_price
            
# =============================================================================
# HTLpump            
# =============================================================================

class HTLpump(qsu.Pump):
    '''
    Pumps used in HTL system
    See qsdsan.sanunits.WWTpump for pipe and pump weight calculation
    See bst.units.Pump for other functions
    All pumps are assumed to be made of stainless steel and specific for sludge.
    Parameters
    ----------
    P : float
        pump pressure
    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
        Valorization of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.
    '''
    
    _N_pump = 1
    _H_ts = 0. # total static head
    _H_p = 0. # total pressure head
    _v = 3 # fluid velocity, [ft/s]
    _C = 110 # Hazen-Williams coefficient for stainless steel (SS)
    _SS_per_pump = 725 * 0.5
    _units = {'Pump pipe stainless steel': 'kg',
              'Pump stainless steel': 'kg'}

    def _design(self):
        super()._design()
        
        pipe, pumps, hdpe = self.design_sludge()

        D = self.design_results
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='Stainless_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = pipe + pumps
    
    def design_sludge(self, Q_mgd=None, N_pump=None, **kwargs):
        '''
        Design pump for handling waste sludge.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or 1

        val_dct = dict(
            L_s=50, # length of suction pipe, [ft]
            L_d=50, # length of discharge pipe, [ft]
            H_ts=0., # H_ds_LIFT (D) - H_ss_LIFT (0)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **val_dct)

        return M_SS_IR_pipe, M_SS_IR_pump, 0
    
    def _design_generic(self, Q_mgd, N_pump=None, L_s=0., L_d=0., H_ts=0., H_p=0.):
        self.Q_mgd = Q_mgd
        self._H_ts = H_ts or self.H_ts
        self._H_p = H_p or self.H_p
        N_pump = N_pump or self.N_pump

        v, C, Q_cfs = self.v, self.C, self.Q_cfs # [ft/s], -, [ft3/s]

        ### Suction side ###
        # Suction pipe (permeate header) dimensions
        OD_s, t_s, ID_s = select_pipe(Q_cfs/N_pump, v) # [in]

        # Suction friction head, [ft]
        self._H_sf = 3.02 * L_s * (v**1.85) * (C**(-1.85)) * ((ID_s/12)**(-1.17))

        ### Discharge side ###
        # Discharge pipe (permeate collector) dimensions
        OD_d, t_d, ID_d = select_pipe(Q_cfs, v)

        # Discharge friction head, [ft]
        self._H_df = 3.02 * L_d * (v**1.85) * (C**(-1.85)) * ((ID_d/12)**(-1.17))

        ### Material usage ###
        # Pipe SS, assume stainless steel, density = 0.29 lbs/in3
        # SS volume for suction, [in3]
        self._N_pump = N_pump
        V_s = N_pump * pi/4*((OD_s)**2-(ID_s)**2) * (L_s*12)
        # SS volume for discharge, [in3]
        V_d = pi/4*((OD_d)**2-(ID_d)**2) * (L_d*12)

        # Total SS mass, [kg]
        M_SS_pipe = 0.29 * (V_s+V_d) * _lb_to_kg
        M_SS_pump = N_pump * self.SS_per_pump
        return M_SS_pipe, M_SS_pump
    
    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
        '''
        return self.F_vol_in*_m3_to_gal*24/1e6
    @Q_mgd.setter
    def Q_mgd(self, i):
        self._Q_mgd = i
        
    @property
    def Q_cfs(self):
        '''[float] Volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal
        
    @property
    def H_ts(self):
        '''[float] Total static head, [ft].'''
        return self._H_ts
    
    @property
    def H_p(self):
        '''[float] Pressure head, [ft].'''
        return self._H_p
    
    @property
    def v(self):
        '''[float] Fluid velocity, [ft/s].'''
        return self._v
    @v.setter
    def v(self, i):
        self._v = i
        
    @property
    def C(self):
        '''[float] Hazen-Williams coefficient to calculate fluid friction.'''
        return self._C
    @C.setter
    def C(self, i):
        self._C = i
        
    @property
    def SS_per_pump(self):
        '''[float] Quantity of stainless steel per pump, [kg/ea].'''
        return self._SS_per_pump
    @SS_per_pump.setter
    def SS_per_pump(self, i):
        self._SS_per_pump = i

# =============================================================================
# HTLHX
# =============================================================================

class HTLHX(qsu.HXutility):
    '''
    Similar to qsdsan.sanunits.HXutility, but can calculate material usage.
    References
    ----------
    .. [1] Seider, W. D., Lewin, D. R., Seader, J. D., Widagdo, S., Gani, R., &
           Ng, M. K. (2017). Product and Process Design Principles. Wiley.
           Chapter 12: Heat Exchanger Design.
    '''
    
    line = qsu.HXutility.line
    _graphics = qsu.HXutility._graphics
    _units = {'Area': 'ft^2',
              'Total tube length': 'ft',
              'Inner pipe weight': 'kg',
              'Outer pipe weight': 'kg',
              'Total steel weight': 'kg',
              'Shell length': 'ft',
              'Shell diameter': 'ft',
              'Shell steel weight': 'kg',
              'Tube weight': 'kg'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    
    def _design(self, duty=None):
        # Set duty and run heat utility
        if duty is None: duty = self.Hnet # Includes heat of formation
        inlet = self.ins[0]
        outlet = self.outs[0] 
        T_in = inlet.T
        T_out = outlet.T
        iscooling = duty < 0.
        if iscooling: # Assume there is a pressure drop before the heat exchanger
            if T_out > T_in: T_in = T_out
        else:
            if T_out < T_in: T_out = T_in
        self.add_heat_utility(duty, T_in, T_out, 
                              heat_transfer_efficiency=self.heat_transfer_efficiency,
                              hxn_ok=True)
        bst.units.HX._design(self)

        D = self.design_results
        if D['Area'] < 150: # double pipe
            # Assume use 1 1/4 nominal size of inner tube, based on [1] page 365
            # Table 12.3, when use Schedule 40, surface area per foot is 0.435 ft2
            # and weight is 2.28 lb steel per foot
            D['Total tube length'] = D['Area']/0.435
            D['Inner pipe weight'] = D['Total tube length']*2.28*_lb_to_kg
            # Assume use 2 nominal size of outer tube, same length as inner tube
            # the weight is 3.66 lb steel per foot
            D['Outer pipe weight'] = D['Total tube length']*3.66*_lb_to_kg
            D['Total steel weight'] = D['Inner pipe weight'] + D['Outer pipe weight']
        else: # shell and tube
            # assume all tubes are 16 ft long and 3/4 in O.D., 16 BWG
            D['Shell length'] = 16
            # D['Total tube length'] = 'N/A'
            single_shell_area = D['Area']/self.N_shells
            if single_shell_area <= 100:
                D['Shell diameter'] = 1
            elif single_shell_area <= 400:
                D['Shell diameter'] = 2
            elif single_shell_area < 1100:
                D['Shell diameter'] = 3
            else:
                D['Shell diameter'] = 3
                self.N_shells = ceil(D['Area']/1100) # ceil shell number
                # if increase the number of N_shells, ft (correction factor) will increase (max = 1),
                # then the required area will decrease, so the calculation here is conservative.
                single_shell_area = D['Area']/self.N_shells

            Shell_design = self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, D['Shell diameter'], D['Shell length'])
            D['Shell steel weight'] = Shell_design['Weight']*_lb_to_kg*self.N_shells
            
            single_tube_area = pi*(3/4)*_in_to_ft*D['Shell length']
            
            D['Total tube length'] = D['Shell length']*single_shell_area/single_tube_area*self.N_shells
            
            # according to [1] page 367, Table 12.4, 3/4 in O.D., 16 BWG tube: 0.520 lb/ft
            D['Tube weight'] = D['Total tube length']*0.520*_lb_to_kg
            
            D['Total steel weight'] = D['Shell steel weight'] + D['Tube weight']
            
        self.construction = (
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = D['Total steel weight']
            
    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3['Carbon steel']
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design

# =============================================================================
# HTL_sludge_centrifuge
# =============================================================================
    
class HTL_sludge_centrifuge(qsu.SludgeThickening, bst.units.SolidsCentrifuge):
    '''
    Similar to qsdsan.sanunits.SludgeCentrifuge, but can calculate material usage.
    References
    ----------
    .. [1] https://dolphincentrifuge.com/wastewater-centrifuge/ (accessed 12-4-2022).
    '''
    
    _units = {'Total pump stainless steel': 'kg',
              'Total pipe stainless steel': 'kg',
              'Centrifige stainless steel': 'kg',
              'Total stainless steel': 'kg'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 sludge_moisture=0.8, solids=(),
                 centrifuge_type='scroll_solid_bowl'):
        qsu.SludgeThickening.__init__(self, ID, ins, outs, thermo, init_with,
                                sludge_moisture=sludge_moisture,
                                solids=solids)
        self.centrifuge_type = centrifuge_type
        ID = self.ID
        eff = self.outs[0].proxy(f'{ID}_eff')
        sludge = self.outs[1].proxy(f'{ID}_sludge')
        self.effluent_pump = HTLpump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
        self.sludge_pump = HTLpump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)

    _run = qsu.SludgeThickening._run

    def _design(self):
        bst.units.SolidsCentrifuge._design(self)
        D = self.design_results
        self.effluent_pump.simulate()
        self.sludge_pump.simulate()
        D['Total pump stainless steel'] = self.effluent_pump.design_results['Pump stainless steel'] +\
                                          self.sludge_pump.design_results['Pump stainless steel']
        D['Total pipe stainless steel'] = self.effluent_pump.design_results['Pump pipe stainless steel'] +\
                                          self.sludge_pump.design_results['Pump pipe stainless steel']
        # based on [1]:
        # when rated capacity <= 80 GPM: weight = 2500 lb
        # when rated capacity (80, 170]: weight = 4000 lb
        # when rated capacity > 170 GPM, use a combination of large and small centrifuges
        
        D['Number of large centrifuge'] = floor(self.F_vol_in*_m3_to_gal/60/170)
        D['Number of small centrifuge'] = 0
        if self.F_vol_in*_m3_to_gal/60 - D['Number of large centrifuge']*170 <= 80:
            D['Number of small centrifuge'] = 1
        else:
            D['Number of large centrifuge'] += 1
        
        D['Centrifige stainless steel'] = (4000*D['Number of large centrifuge'] + 2500*D['Number of small centrifuge'])*_lb_to_kg
        D['Total stainless steel'] = D['Total pump stainless steel'] + D['Total pipe stainless steel'] + D['Centrifige stainless steel']
        
        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='Stainless_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = D['Total stainless steel']
        
    def _cost(self):
        qsu.SludgeThickening._cost(self)

# =============================================================================
# HTLCHP
# =============================================================================
    
class HTLCHP(qsu.CHP):
    '''
    Similar to qsdsan.sanunits.CHP, but can calculate material usage.
    References
    ----------
    .. [1] Havukainen, J.; Nguyen, M. T.; Väisänen, S.; Horttanainen, M.
           Life Cycle Assessment of Small-Scale Combined Heat and Power Plant:
           Environmental Impacts of Different Forest Biofuels and Replacing
           District Heat Produced from Natural Gas. Journal of Cleaner
           Production 2018, 172, 837–846.
           https://doi.org/10.1016/j.jclepro.2017.10.241.
    '''
    
    _units = {'Steel': 'kg',
              'Furnace': 'kg',
              'Concrete': 'kg',
              'Reinforcing steel': 'kg'}
    
    def _design(self):
        super()._design()
        D = self.design_results
        
        # material calculation based on [1], linearly scaled on power (kW)
        # in [1], a 580 kW CHP:
        # steel: 20098 kg
        # furnace: 12490 kg
        # reinforced concrete: 15000 kg (concrete + reinforcing steel)
        # 1 m3 reinforced concrete: 98 v/v% concrete with a density of 2500 kg/m3 (2450 kg)
        #                            2 v/v% reinforcing steel with a density of 7850 kg/m3 (157kg)
        factor = self.H_net_feeds/3600/580
        D['Steel'] = factor*20098
        D['Furnace'] = factor*12490
        D['Concrete'] = factor*15000*2450/(2450 + 157)
        D['Reinforcing steel'] = factor*15000*157/(2450 + 157)
        
        self.construction = (
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            Construction('furnace', linked_unit=self, item='Furnace', quantity_unit='kg'),
            Construction('concrete', linked_unit=self, item='concrete', quantity_unit='kg'),
            Construction('reinforcing_steel', linked_unit=self, item='Reinforcing_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = D['Steel']
        self.construction[1].quantity = D['Furnace']
        self.construction[2].quantity = D['Concrete']
        self.construction[3].quantity = D['Reinforcing steel']
        
        
        
        
        
        
        
# =============================================================================
# HTLCHP_2
# =============================================================================

from thermosteam.reaction import ParallelReaction

class HTLCHP_2(qsu.CHP):
    '''
    Similar to qsdsan.sanunits.CHP, but can calculate material usage.
    References
    ----------
    .. [1] Havukainen, J.; Nguyen, M. T.; Väisänen, S.; Horttanainen, M.
           Life Cycle Assessment of Small-Scale Combined Heat and Power Plant:
           Environmental Impacts of Different Forest Biofuels and Replacing
           District Heat Produced from Natural Gas. Journal of Cleaner
           Production 2018, 172, 837–846.
           https://doi.org/10.1016/j.jclepro.2017.10.241.
    '''
    _N_ins = 2
    _N_outs = 2
    _units = {'Steel': 'kg',
              'Furnace': 'kg',
              'Concrete': 'kg',
              'Reinforcing steel': 'kg'}
    
    def _design(self):
        feed, air = self.ins
        emission, ash = self.outs
        for i in (air, ash):
            i.empty()
        feed.phase = air.phase = emission.phase = 'g'
        ash.phase = 's'
        emission.P = ash.P = 101325
        emission.T = ash.T = 298.15

        cmps = self.components
        rxns = []
        for cmp in cmps:
            if cmp.locked_state in ('l', 's') and (not cmp.organic or cmp.degradability=='Undegradable'):
                continue
            rxn = cmp.get_combustion_reaction()
            if rxn:
                rxns.append(rxn)
        combustion_rxns = self.combustion_reactions = ParallelReaction(rxns)

        def react(natural_gas_flow=0):
            emission.copy_flow(feed)
            combustion_rxns.force_reaction(emission.mol)
            air.imol['O2'] = -emission.imol['O2']
            emission.imol['N2'] = air.imol['N2'] = air.imol['O2']/0.21*0.79
            emission.imol['O2'] = 0
            H_net_feed = feed.H + feed.HHV - emission.H # subtracting the energy in emission
            return H_net_feed

        self.H_net_feeds = react(0)
        
        pu = self.power_utility
        pu.production = self.H_net_feeds/3600*self.combined_eff

        ash_IDs = [i.ID for i in cmps if not i.formula]
        ash.copy_flow(emission, IDs=tuple(ash_IDs), remove=True)
        
        D = self.design_results
        
        # material calculation based on [1], linearly scaled on power (kW)
        # in [1], a 580 kW CHP:
        # steel: 20098 kg
        # furnace: 12490 kg
        # reinforced concrete: 15000 kg (concrete + reinforcing steel)
        # 1 m3 reinforced concrete: 98 v/v% concrete with a density of 2500 kg/m3 (2450 kg)
        #                            2 v/v% reinforcing steel with a density of 7850 kg/m3 (157kg)
        factor = self.H_net_feeds/3600/580
        D['Steel'] = factor*20098
        D['Furnace'] = factor*12490
        D['Concrete'] = factor*15000*2450/(2450 + 157)
        D['Reinforcing steel'] = factor*15000*157/(2450 + 157)
        
        self.construction = (
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            Construction('furnace', linked_unit=self, item='Furnace', quantity_unit='kg'),
            Construction('concrete', linked_unit=self, item='concrete', quantity_unit='kg'),
            Construction('reinforcing_steel', linked_unit=self, item='Reinforcing_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = D['Steel']
        self.construction[1].quantity = D['Furnace']
        self.construction[2].quantity = D['Concrete']
        self.construction[3].quantity = D['Reinforcing steel']        













# =============================================================================
# HTL_storage_tank
# =============================================================================
        
class HTL_storage_tank(qsu.StorageTank):
    '''
    Similar to qsdsan.sanunits.StorageTank, but can calculate material usage.
    '''
    
    _units = {'Diameter': 'ft',
              'Length': 'ft',
              'Wall thickness': 'in',
              'Weight': 'kg'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  vessel_type=None, tau=None, V_wf=None,
                  vessel_material=None, kW_per_m3=0.,
                  init_with='WasteStream', F_BM_default=None, length_to_diameter=2):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                      vessel_type=vessel_type, tau=tau, V_wf=V_wf,
                      vessel_material=vessel_material, kW_per_m3=kW_per_m3,
                      init_with=init_with, F_BM_default=F_BM_default)
        self.length_to_diameter = length_to_diameter
    
    def _design(self):
        super()._design()
        D = self.design_results
        
        Diameter = (4*D['Total volume']/pi/self.length_to_diameter)**(1/3)
        Diameter *= _m_to_ft # convert from m to ft
        L = Diameter * self.length_to_diameter # ft

        Tank_design = self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, Diameter, L)
        
        D['Diameter'] = Diameter
        D['Length'] = L
        D['Wall thickness'] = Tank_design['Wall thickness']
        D['Material'] = self.vessel_material
        D['Weight'] = Tank_design['Weight']*_lb_to_kg
        
        self.construction = (
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = D['Weight']
        
    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3[self.vessel_material]
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design

# =============================================================================
# HTLcompressor
# =============================================================================

class HTLcompressor(IsothermalCompressor):
    '''
    Similar to biosteam.units.IsothermalCompressor, but can calculate number of units.
    '''
    
    def _design(self):
        super()._design()
        D = self.design_results
        power = D['Ideal power']/D['Driver efficiency']
        D['Number of 300 kW unit'] = floor(power/300)
        D['Number of 4 kW unit'] = 0
        if (power - D['Number of 300 kW unit']*300) <= 60:
            # according to Ecoinvent 3: the impact of at most 15 4 kW unit is smaller than 1 300 kW unit
            # therefore, if the rest of power is smaller than 60 kW, use multiple small units
            # else, add one large unit
            D['Number of 4 kW unit'] = ceil((power - D['Number of 300 kW unit']*300)/4)
        else:
            D['Number of 300 kW unit'] += 1
        
        self.construction = (
            Construction('compressor_4kW', linked_unit=self, item='Compressor_4kW', quantity_unit='ea'),
            Construction('compressor_300kW', linked_unit=self, item='Compressor_300kW', quantity_unit='ea')
            )
        self.construction[0].quantity = D['Number of 4 kW unit']
        self.construction[1].quantity = D['Number of 300 kW unit']

# =============================================================================
# HTLflash
# =============================================================================

class HTLflash(Flash):
    '''
    Similar to biosteam.units.Flash, but includes construction.
    '''
    
    def _design(self):
        super()._design()
        D = self.design_results
        
        self.construction = (
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            )
        if 'Weight' in D.keys():
            self.construction[0].quantity = D['Weight']*_lb_to_kg

# =============================================================================
# HTLdistillation
# =============================================================================

class HTLdistillation(BinaryDistillation):
    '''
    Similar to biosteam.units.BinaryDistillation, but includes construction.
    '''
    
    def _design(self):
        super()._design()
        D = self.design_results
        
        self.construction = (
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            )
        self.construction[0].quantity = (D['Rectifier weight'] + D['Stripper weight'])*_lb_to_kg

# =============================================================================
# HTLHXN
# =============================================================================

class HTLHXN(qsu.HeatExchangerNetwork):
    '''
    Similar to qsdsan.sanunits.HeatExchangerNetwork, but enable LCA.
    '''
    
    def __init__(self, ID=''):
        super().__init__(ID=ID)
        self.construction = ()
        self.transportation = ()