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

References:

(1) Snowden-Swan, L. J.; Li, S.; Jiang, Y.; Thorson, M. R.; Schmidt, A. J.;
    Seiple, T. E.; Billing, J. M.; Santosa, D. M.; Hart, T. R.; Fox, S. P.;
    Cronin, D.; Kallupalayam Ramasamy, K.; Anderson, D. B.; Hallen, R. T.;
    Fonoll-Almansa, X.; Norton, J. Wet Waste Hydrothermal Liquefaction and
    Biocrude Upgrading to Hydrocarbon Fuels: 2021 State of Technology;
    PNNL-32731; Pacific Northwest National Lab. (PNNL), Richland, WA
    (United States), 2022. https://doi.org/10.2172/1863608.

(2) Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
    Quantitative Evaluation of an Integrated System for Valorization of
    Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
    Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
    https://doi.org/10.1021/acs.est.8b04035.
    
(3) Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
    Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
    Liquefaction Products from Feedstock Biochemical Composition.
    Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    
(4) Snowden-Swan, L. J.; Zhu, Y.; Jones, S. B.; Elliott, D. C.; Schmidt, A. J.; 
    Hallen, R. T.; Billing, J. M.; Hart, T. R.; Fox, S. P.; Maupin, G. D. 
    Hydrothermal Liquefaction and Upgrading of Municipal Wastewater Treatment 
    Plant Sludge: A Preliminary Techno-Economic Analysis; 
    PNNL--25464, 1258731; 2016; https://doi.org/10.2172/1258731.

(5) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to
    Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
    PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.
    
(6) Hao, S.; Choi, Y.-J.; Wu, B.; Higgins, C. P.; Deeb, R.; Strathmann, T. J.
    Hydrothermal Alkaline Treatment for Destruction of Per- and Polyfluoroalkyl
    Substances in Aqueous Film-Forming Foam. Environ. Sci. Technol.
    2021, 55(5), 3283–3295. https://doi.org/10.1021/acs.est.0c06906.
    
(7)	Matayeva, A.; Rasmussen, S. R.; Biller, P. Distribution of Nutrients and
    Phosphorus Recovery in Hydrothermal Liquefaction of Waste Streams.
    BiomassBioenergy 2022, 156, 106323.
    https://doi.org/10.1016/j.biombioe.2021.106323.
    
(8) Doran, P. M. Chapter 11 - Unit Operations. In Bioprocess Engineering
    Principles (Second Edition); Doran, P. M., Ed.; Academic Press: London,
    2013; pp 445–595. https://doi.org/10.1016/B978-0-12-220851-5.00011-3.
    
(9) Elliott, D. C.; Neuenschwander, G. G.; Hart, T. R.; Rotness, L. J.;
    Zacher, A. H.; Santosa, D. M.; Valkenburg, C.; Jones, S. B.;
    Rahardjo, S. A. T. Catalytic Hydrothermal Gasification of Lignin-Rich
    Biorefinery Residues and Algae Final Report. 87.

(10) Spiller, L. L. Determination of Ammonia/Air Diffusion Coefficient Using
    Nafion Lined Tube. Analytical Letters 1989, 22 (11–12), 2561–2573.
    https://doi.org/10.1080/00032718908052375.

(11) Scheepers, D. M.; Tahir, A. J.; Brunner, C.; Guillen-Burrieza, E.
    Vacuum Membrane Distillation Multi-Component Numerical Model for Ammonia
    Recovery from Liquid Streams. Journal of Membrane Science
    2020, 614, 118399. https://doi.org/10.1016/j.memsci.2020.118399.

(12) Ding, Z.; Liu, L.; Li, Z.; Ma, R.; Yang, Z. Experimental Study of Ammonia
    Removal from Water by Membrane Distillation (MD): The Comparison of Three
    Configurations. Journal of Membrane Science 2006, 286 (1), 93–103.
    https://doi.org/10.1016/j.memsci.2006.09.015.
    
(13) Towler, G.; Sinnott, R. Chapter 14 - Design of Pressure Vessels.
    In Chemical Engineering Design (Second Edition); Towler, G., Sinnott, R.,
    Eds.; Butterworth-Heinemann: Boston, 2013; pp 563–629.
    https://doi.org/10.1016/B978-0-08-096659-5.00014-6.
    
(14) Al-Obaidani, S.; Curcio, E.; Macedonio, F.; Di Profio, G.; Al-Hinai, H.;
    Drioli, E. Potential of Membrane Distillation in Seawater Desalination:
    Thermal Efficiency, Sensitivity Study and Cost Estimation.
    Journal of Membrane Science 2008, 323 (1), 85–98.
    https://doi.org/10.1016/j.memsci.2008.06.006.
    
(15) Kogler, A.; Farmer, M.; Simon, J. A.; Tilmans, S.; Wells, G. F.;
    Tarpeh, W. A. Systematic Evaluation of Emerging Wastewater Nutrient Removal
    and Recovery Technologies to Inform Practice and Advance Resource
    Efficiency. ACS EST Eng. 2021, 1 (4), 662–684.
    https://doi.org/10.1021/acsestengg.0c00253.

(16) Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
    Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
    2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
    https://doi.org/10.2172/1855835.
'''

import biosteam as bst
from qsdsan import SanUnit
from qsdsan.sanunits import HXutility, Pump
from biosteam.units import IsothermalCompressor
from exposan.htl._components import create_components
from biosteam.units.design_tools import PressureVessel
from biosteam.units.design_tools.cost_index import CEPCI_by_year as CEPCI
from math import pi, ceil, log
from biosteam.exceptions import DesignError
from biosteam import Stream
from biosteam.units.decorators import cost
from thermosteam import indexer, equilibrium
from exposan.htl._process_settings import load_process_settings

__all__ = ('Reactor',
           'SludgeLab',
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
           'FuelMixer')

cmps = create_components()

load_process_settings()

yearly_operation_hour = 7920 # Jones

# =============================================================================
# Reactor
# =============================================================================

class Reactor(SanUnit, PressureVessel, isabstract=True):
    '''
    Create an abstract class for reactor unit, purchase cost of the reactor
    is based on volume calculated by residence time.
    Parameters
    ----------
    ins : stream
        Inlet.
    outs : stream
        Outlet.
    tau : float
        Residence time [hr].
    V_wf : float
        Fraction of working volume over total volume.
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
        `<https://doi.org/10.1039/C3EM00711A>`_.
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
            ins_F_vol = self.F_vol_in
            for i in range(len(self.ins)):
                ins_F_vol -= self.ins[i].ivol['H2']
            # not include gas (e.g. H2)
            V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor

        if self.N:
            if self.V:
                D = (4*self.V/pi/length_to_diameter)**(1/3)
                D *= 3.28084 # convert from m to ft
                L = D * length_to_diameter

                Design['Total volume'] = self.V*self.N
                Design['Single reactor volume'] = self.V
                Design['Number of reactors'] = self.N
            else:
                V_reactor = V_total/self.N
                D = (4*V_reactor/pi/length_to_diameter)**(1/3)
                D *= 3.28084 # convert from m to ft
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
                D *= 3.28084 # convert from m to ft
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

class SludgeLab(SanUnit):
    
    '''
    SludgeLab is a fake unit that can set up sludge biochemical compositions
    and calculate sludge elemental compositions.
    
    Model method: just _run, no _design or _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        ww
    outs: Iterable (stream)
        sludge, treated
    '''

    _m3perh_2_MGD = 1/157.725491

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', 
                 ww_2_dry_sludge=0.94, # ww (MGD) to dw sludge (metric ton/day) # SS PNNL 2022
                 # how much metric ton/day sludge can be produced by 1 MGD of ww
                 sludge_moisture=0.99, sludge_dw_ash=0.257, 
                 sludge_afdw_protein=0.463, sludge_afdw_lipid=0.204,
                 lipid_2_C=0.750, protein_2_C=0.545,
                 carbo_2_C=0.400, C_2_H=0.143,
                 protein_2_N=0.159, N_2_P=0.393,
                 operation_hour=yearly_operation_hour,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ww_2_dry_sludge = ww_2_dry_sludge
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_ash = sludge_dw_ash
        self.sludge_afdw_protein = sludge_afdw_protein
        self.sludge_afdw_lipid = sludge_afdw_lipid
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
            
        self.sludge_dw = ww.F_vol*self._m3perh_2_MGD*self.ww_2_dry_sludge*1000/24

        sludge.imass['H2O'] = self.sludge_dw/(1-self.sludge_moisture)*self.sludge_moisture
        sludge.imass['Sludge_ash'] = self.sludge_dw*self.sludge_dw_ash

        sludge_afdw = self.sludge_dw*(1 - self.sludge_dw_ash)
        sludge.imass['Sludge_protein'] = sludge_afdw*self.sludge_afdw_protein
        sludge.imass['Sludge_lipid'] = sludge_afdw*self.sludge_afdw_lipid
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
              self.sludge_O_ratio/8)) # Li 2017, in MJ/kg
        
    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
# KOdrum
# =============================================================================

class KOdrum(Reactor):
    
    '''
    konckout drum: a HTL auxiliary unit
    
    Model method: just design and cost.
    
    Design: for each HTL reactor (4 in total), add a KOdrum to separate liquid
    from gas. V = 4230 gal based on NREL 2013.
    
    Parameters
    ----------
    ins: Iterable (stream)
        dewatered_sludge
    outs: Iterable (stream)
        biochar, HTLaqueous, biocrude, offgas
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 P=3049.7*6894.76, tau=0, V_wf=0,
                 length_to_diameter=2, N=4,
                 V=4230*0.00378541, # 4230 gal to m3, NREL 2013
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
    
    Model method: empirical model (based on MCA model and experimental data).
    
    Design: tau = 15/60 based on Jones 2014. V_wf = 0.3 calculate based on data
    in NREL 2013. N = 4 based on NREL 2013. Instead of a long PFR, use length_to_
    diameter = 2.
    
    Parameters
    ----------
    ins: Iterable (stream)
        dewatered_sludge
    outs: Iterable (stream)
        biochar, HTLaqueous, biocrude, offgas
    '''
    
    auxiliary_unit_names=('heat_exchanger','kodrum')
    
    _kg_2_lb = 2.20462
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 lipid_2_biocrude = 0.846, # revised MCA
                 protein_2_biocrude = 0.445, # revised MCA
                 carbo_2_biocrude = 0.205, # revised MCA
                 protein_2_gas = 0.074, # revised MCA
                 carbo_2_gas = 0.418, # revised MCA
                 biocrude_C_slope = -8.37, # MCA
                 biocrude_C_intercept = 68.55, # MCA
                 biocrude_N_slope = 0.133, # MCA
                 biocrude_H_slope = -2.61, # MCA
                 biocrude_H_intercept = 8.20, # MCA
                 biochar_C_slope = 1.75, # MCA
                 biocrude_moisture_content=0.063, # Jones
                 biochar_P_recovery_ratio=0.86, # Matayeva: 0.84-0.88
                 gas_composition={'CH4':0.050, 'C2H6':0.032,
                                  'CO2':0.918}, # Jones
                 biochar_pre=3029.7*6894.76, # Jones
                 HTLaqueous_pre=30*6894.76, # Jones
                 biocrude_pre=30*6894.76, # Jones
                 offgas_pre=30*6894.76, # Jones
                 eff_T=60+273.15, # Jones
                 
                 P=None, tau=15/60, V_wf=0.3,
                 length_to_diameter=2, N=4, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',         
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
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=eff_T)
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

    _N_ins = 1
    _N_outs = 4
    _units= {'Treatment capacity': 'lb/h'} # separator
    
    def _run(self):
        
        dewatered_sludge = self.ins[0]
        biochar, HTLaqueous, biocrude, offgas = self.outs
        
        self.sludgelab = self.ins[0]._source.ins[0]._source.ins[0].\
                         _source.ins[0]._source.ins[0]._source
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid'] +\
                                dewatered_sludge.imass['Sludge_protein'] +\
                                dewatered_sludge.imass['Sludge_carbo']
        # just use afdw in revised MCA model, other places use dw
        
        afdw_protein_ratio = self.sludgelab.sludge_afdw_protein
        afdw_lipid_ratio = self.sludgelab.sludge_afdw_lipid
        afdw_carbo_ratio = self.sludgelab.sludge_afdw_carbo

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
                                
        HTLaqueous.imass['H2O'] = dewatered_sludge.imass['H2O'] -\
                                  biocrude.imass['H2O'] +\
                                  dewatered_sludge.imass['Sludge_ash']
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
        return (self.sludgelab.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100
        # revised MCA model
    
    @property
    def biocrude_H_ratio(self):
        return (self.sludgelab.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100
        # revised MCA model

    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.sludgelab.sludge_dw_protein # revised MCA model
    
    @property
    def biocrude_C(self):
        return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.sludgelab.sludge_C)

    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio

    @property
    def biocrude_N(self):
        return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.sludgelab.sludge_N)
    
    @property
    def biocrude_HHV(self):
        return 30.74 - 8.52*self.sludgelab.AOSc +\
               0.024*self.sludgelab.sludge_dw_protein # Li 2017, in MJ/kg
               
    @property
    def energy_recovery(self):
        return self.biocrude_HHV*self.outs[2].imass['Biocrude']/\
               (self.sludgelab.outs[0].F_mass -\
               self.sludgelab.outs[0].imass['H2O'])/self.sludgelab.sludge_HHV
        # Li 2017
        
    @property
    def biochar_C_ratio(self):
        return min(self.biochar_C_slope*self.sludgelab.sludge_dw_carbo, 0.65)
    # revised MCA model

    @property
    def biochar_C(self):
        return min(self.outs[0].F_mass*self.biochar_C_ratio, self.sludgelab.sludge_C - self.biocrude_C)

    @property
    def biochar_P(self):
        return self.sludgelab.sludge_P*self.biochar_P_recovery_ratio

    @property
    def offgas_C(self):
        carbon = 0
        for name in self.gas_composition.keys():
            carbon += self.outs[3].imass[name]*cmps[name].i_C
        return min(carbon, self.sludgelab.sludge_C - self.biocrude_C - self.biochar_C)
               
    # C and N in aqueous phase are calculated based on mass balance closure
    # in the case that HTLaqueous_C is less than 0, which
    # is likely to happen when the carbo% is high,
    # we assume the HTLaqueous 0 and excess biochar_C is ok
    
    @property
    def HTLaqueous_C(self):
        return self.sludgelab.sludge_C - self.biocrude_C - self.biochar_C - self.offgas_C
    # !!! keep in mind, in reality, HTLaqueous_C cannot be 0.
    # when inputs for HTL become extreme, the production will deviate from MCA model

    @property
    def HTLaqueous_N(self):
        return self.sludgelab.sludge_N - self.biocrude_N
        
    @property
    def HTLaqueous_P(self):
        return self.sludgelab.sludge_P*(1 - self.biochar_P_recovery_ratio)

    def _design(self):
        
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass*self._kg_2_lb
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.copy_like(hx_ins0)
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx_outs0.T = hx.T
        hx_ins0.P = hx_outs0.P = self.outs[1].P
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        self.kodrum.simulate()

        self.P = self.ins[0].P
        Reactor._design(self)
        
    def _cost(self):
        Reactor._cost(self)
        self._decorated_cost()
        
# =============================================================================
# Acid Extraction
# =============================================================================

class AcidExtraction(Reactor):
    
    '''
    H2SO4 is added to biochar from HTL to extract P. 
    
    Model method: assume P recovery ratio, add filters in _design and _cost.
    
    Design: tau, V_wf, and kW_per_m3 use MixTank default values. Material uses 304.
    
    Parameters
    ----------
    ins: Iterable (stream)
        biochar, acid
    outs: Iterable (stream)
        residual, extracted
    '''
    
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', acid_vol=7, P_acid_recovery_ratio=0.8,
                 P=None, tau=2, V_wf=0.8, # tau: Zhu
                 length_to_diameter=2, N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0.0985, # use MixTank default value
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
            if self.HTL.biochar_P == 0:
                residual.copy_like(biochar)
                residual.price = 0.1 # $/kg made-up value based on online info, will change later
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
                residual.price = 0.1 # made-up value
                
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
        if self.ins[0].F_mass <= 0 or self.HTL.biochar_P == 0:
            pass
        else:
            self.P = self.ins[1].P
            Reactor._design(self)
        
    def _cost(self):
        if self.ins[0].F_mass <= 0 or self.HTL.biochar_P == 0:
            pass
        else:
            Reactor._cost(self)
    
# =============================================================================
# HTL mixer
# =============================================================================

class HTLmixer(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the mixture of HTL
    aqueous and AcidEx effluent.
    
    Model method: elements separation, don't need _design and _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        HTLaqueous, extracted
    outs: Iterable (stream)
        mixture
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
        # assume HTLaqueous pH = 9 Li 2018 SI (9.08 ± 0.30)
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
    
    Model method: use experimental data.
    
    Parameters
    ----------
    ins: Iterable (stream)
    flow_in
    outs: Iterable (stream)
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
    extracted_P and HTL aqueous are mixed together (Mixer) before adding
    MgCl2 and struvite precipitation.
    
    Model method: P recovery rate with uncertainty from literature data.
    If mol(N)<mol(P), add NH4Cl to mol(N):mol(P)=1:1
    
    Design: tau, V_wf, and kW_per_m3 use MixTank default values. Material uses carbon steel.
    
    Parameters
    ----------
    ins: Iterable (stream)
        mixture, supply_MgCl2, supply_NH4Cl
    outs: Iterable (stream)
        struvite, effluent
    '''
    
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', 
                 target_pH = 9,
                 Mg_P_ratio=2.5,
                 P_pre_recovery_ratio=0.99, # Zhu
                 P_in_struvite=0.126,
                 HTLaqueous_NH3_N_2_total_N = 0.853, # Jones
                 
                 P=None, tau=1, V_wf=0.8, # tau: Zhu
                 length_to_diameter=2, N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0.0985, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Carbon steel', # basic condition
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_pH = target_pH
        self.Mg_P_ratio = Mg_P_ratio
        self.P_pre_recovery_ratio = P_pre_recovery_ratio
        self.P_in_struvite = P_in_struvite
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
            base.price = 0.2
    
            if mixture.imass['P']/30.973762 > self.Mg_P_ratio*mixture.imass['N']*\
                                              self.HTLaqueous_NH3_N_2_total_N/14.0067:
            # N:P >= Mg:P
                supply_NH4Cl.imass['NH4Cl'] = (mixture.imass['P']/30.973762 -\
                                               self.Mg_P_ratio*mixture.imass['N']*\
                                               self.HTLaqueous_NH3_N_2_total_N/14.0067)*53.491
            supply_NH4Cl.price = 0.13
            # make sure N:P >= 1:1
            
            supply_MgCl2.imass['MgCl2'] = (mixture.imass['P']/30.973762*self.Mg_P_ratio -\
                                           base.imass['MgO']/40.3044)*95.211
                                           # Mg:P = 2:1 (1.5:1 - 4:1)
            supply_MgCl2.price = 0.5452
            struvite.imass['Struvite'] = mixture.imass['P']*\
                                         self.P_pre_recovery_ratio/\
                                         self.P_in_struvite
            struvite.price = 0.6784
            supply_MgCl2.phase = supply_NH4Cl.phase = base.phase = 's'
            
            effluent.copy_like(mixture)
            effluent.imass['P'] -= struvite.imass['Struvite']*self.P_in_struvite
            effluent.imass['N'] += supply_NH4Cl.imass['NH4Cl']*14.0067/53.491 -\
                                   struvite.imass['Struvite']*\
                                   self.P_in_struvite/30.973762*14.0067
            effluent.imass['H2O'] += (supply_MgCl2.imass['MgCl2'] +\
                                      supply_NH4Cl.imass['NH4Cl'] +\
                                      base.imass['NaOH'] -\
                                      struvite.imass['Struvite']*\
                                      (1 - self.P_in_struvite*\
                                      (1+14.0067/30.973762)))
            struvite.phase = 's'    
                
            struvite.T = mixture.T
            effluent.T = mixture.T
        
    @property
    def struvite_P(self):
        return self.outs[0].imass['Struvite']*self.P_in_struvite

    @property
    def struvite_N(self):
        return self.struvite_P*14.0067/30.973762

    def _design(self):
        if self.HTLmixer.outs[0].imass['P'] == 0:
            pass
        else:
            self.P = self.ins[0].P
            Reactor._design(self)
        
    def _cost(self):
        if self.HTLmixer.outs[0].imass['P'] == 0:
            pass
        else:
            Reactor._cost(self)
    
# =============================================================================
# CHG
# =============================================================================

@cost(basis='Treatment capacity', ID='Hydrocyclone', units='lb/h',
      cost=5000000, S=968859,
      CE=CEPCI[2009], n=0.65, BM=2.1)
# hydrocyclone

class CHG(Reactor, SludgeLab):
   
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    
    Model method: use experimental data, assume no NH3 loss for now.
    
    Design: tau = 1/3 h and V_wf = 0.5 based on Davis NREL 2018. Assuming no
    headspace. N=6 based on Jones 2014.
    
    Parameters
    ----------
    ins: Iterable (stream)
        chg_in, catalyst_in
    outs: Iterable (stream)
        chg_out, catalyst_out
    '''
    
    auxiliary_unit_names=('pump','heat_ex_heating','heat_ex_cooling')
    
    _kg_2_lb = 2.20462
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Sulfur guard': 2.0}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 
                 pump_pressure=3089.7*6894.76,
                 heat_temp=350+273.15,
                 cool_temp=60+273.15,
                 
                 WHSV=3.56, # wt./hr per wt. catalyst
                 catalyst_lifetime=1*yearly_operation_hour, # 1 year Jones
                 catalyst_price=134.53,
                 gas_composition={'CH4':0.527,
                                  'CO2':0.432,
                                  'C2H6':0.011,
                                  'C3H8':0.030,
                                  'H2':0.0001},
                 # Jones
                 # will not be a variable in uncertainty/sensitivity analysis
                 gas_C_2_total_C=0.598, # Jones
                 
                 P=None, tau=20/60, void_fraction=0.5,
                 length_to_diameter=2, N=6, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
        self.pump_pressure = pump_pressure
        self.heat_temp = heat_temp
        self.cool_temp = cool_temp
        
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_price = catalyst_price
        self.gas_composition = gas_composition
        self.gas_C_2_total_C = gas_C_2_total_C
        
        pump_in = bst.Stream(f'{ID}_pump_in')
        pump_out = bst.Stream(f'{ID}_pump_out')
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=pump_pressure)
        
        hx_ht_in = bst.Stream(f'{ID}_hx_ht_in')
        hx_ht_out = bst.Stream(f'{ID}_hx_ht_out')
        self.heat_ex_heating = HXutility(ID=f'.{ID}_hx_ht', ins=hx_ht_in, outs=hx_ht_out, T=heat_temp)

        hx_cl_in = bst.Stream(f'{ID}_hx_cl_in')
        hx_cl_out = bst.Stream(f'{ID}_hx_cl_out')
        self.heat_ex_cooling = HXutility(ID=f'.{ID}_hx_cl', ins=hx_cl_in, outs=hx_cl_out, T=cool_temp)
        
        self.P = P
        self.tau = tau
        self.V_wf = void_fraction # no headspace, gases produced will be vented, so V_wf = void fraction Elliott
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
    _units= {'Treatment capacity': 'lb/h'} # hydrocyclone
        
    def _run(self):
        
        chg_in, catalyst_in = self.ins
        chg_out, catalyst_out = self.outs
        
        if self.ins[0].imass['C'] == 0:
            chg_out.copy_like(chg_in)
        else:
            catalyst_in.imass['CHG_catalyst'] = chg_in.F_mass/self.WHSV/self.catalyst_lifetime
            catalyst_in.phase = 's'
            catalyst_out.copy_like(catalyst_in)
            # catalysts amount is quite low compared to the main stream, therefore do not consider
            # heating/cooling of catalysts
            catalyst_in.price = self.catalyst_price
            
            gas_C_ratio = 0
            for name, ratio in self.gas_composition.items():
                gas_C_ratio += ratio*cmps[name].i_C
            
            gas_mass = chg_in.imass['C']*self.gas_C_2_total_C/gas_C_ratio
            
            for name,ratio in self.gas_composition.items():
                chg_out.imass[name] = gas_mass*ratio
                
            chg_out.imass['H2O'] = chg_in.F_mass - gas_mass
            # all C, N, and P are accounted in H2O here, but will be calculated as
            # properties.
                
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
        if self.ins[0].imass['C'] == 0:
            pass
        else:
            Design = self.design_results
            Design['Treatment capacity'] = self.ins[0].F_mass*self._kg_2_lb
        
            
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

            self.P = self.ins[0].P
            Reactor._design(self)
    
    def _cost(self):
        if self.ins[0].imass['C'] == 0:
            pass
        else:
            Reactor._cost(self)
            purchase_costs = self.baseline_purchase_costs
            current_cost = 0 # cost w/o sulfur guard
            for item in purchase_costs.keys():
                current_cost += purchase_costs[item]
            purchase_costs['Sulfur guard'] = current_cost*0.05
            self._decorated_cost()
    
# =============================================================================
# Membrane Distillation
# =============================================================================

class MembraneDistillation(SanUnit):
    
    '''
    Membrane distillation recovers nitrogen as ammonia sulfate based on vapor
    pressure difference across the hydrophobic membrane.
    
    Model method: 
        1. calculate NH3 in feed based on pH
        2. calculate NH3/H2O VLE
        3. calculate NH3 flux and removal rate
        
    Design: based on Al-Obaidani et al.

    Parameters
    ----------
    ins: Iterable (stream)
        influent, acid, base
    outs: Iterable (stream)
        ammoniasulfate, ww
    '''
    
    _F_BM_default = {'Membrane': 1}
    
    _units = {'Area': 'm2',
              'Total volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 inffluent_pH=8.16, # CHG effluent pH: 8.16 ± 0.25 Li 2018
                 target_pH=10,
                 N_S_ratio=2,
                 # S is excess since not all N can be transfered to form ammonia sulfate
                 # for now, assume N_S_ratio = 2 is ok
                 m2_2_m3=1/1200, # specific surface area, for hollow fiber membrane, Doran
                 Dm=2.28*10**(-5), # (2.28 ± 0.12)*10^-5 m^2/s NH3 molecular diffusity in air Spiller
                 # (underestimate, this value may be at 15 or 25 C, our feed is 60 C, should be higher)
                 porosity=0.9, # Scheepers
                 thickness=7*10**(-5), # m, Scheepers
                 tortuosity=1.2, # Scheepers
                 Henry=1.61*10**(-5), # atm*m3/mole
                 # https://webwiser.nlm.nih.gov/substance?substanceId=315&identifier=\
                 # Ammonia&identifierType=name&menuItemId=81&catId=120#:~:text=The%20\
                 # Henry's%20Law%20constant%20for,m%2Fmole(2). (accessed 11-11-2022)
                 Ka=1.75*10**(-5), # overall mass transfer coefficient 1.2~2.3*10^-5 m/s Ding
                 capacity=6.01, # kg/m2/h, Al-Obaidani
                 # permeate flux, similar values can be found in many other papers ([176], [222] ,[223] in Kogler Stanford)
                 # [177], [223] in Kogler Stanford show high nitrogen recovery ratio (>85% under optimal conditions)
                 membrane_price=93.29,
                 ammonium_sulfate_price=0.3236, # convert to solution price
                 # Al-Obaidani: $90/m2, $2008
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.inffluent_pH = inffluent_pH
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
        self.ammonium_sulfate_price = ammonium_sulfate_price

    _N_ins = 4
    _N_outs = 3
    
    def _run(self):
        
        influent, acid, base, mem_in = self.ins
        ammoniumsulfate, ww, mem_out = self.outs
        
        
        self.CHG = self.ins[0]._source.ins[0]._source.ins[0]._source
        
        if self.CHG.CHGout_N == 0:
            ww.copy_like(influent)
        else:
            NaOH_conc = 10**(self.target_pH - 14) - 10**(self.inffluent_pH - 14)
            NaOH_mol = NaOH_conc*self.ins[0].F_mass
            base.imass['NaOH'] = NaOH_mol*39.997/1000
            base.price = 0.526
            
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
            ammoniumsulfate.imass['H2O'] = acid.imass['H2O']
            ammoniumsulfate.imass['H2SO4'] = acid.imass['H2SO4'] +\
                                             ammoniumsulfate.imass['NH42SO4']/\
                                             132.14*28.0134 -\
                                             ammoniumsulfate.imass['NH42SO4']
            ammoniumsulfate.price = self.ammonium_sulfate_price*ammoniumsulfate.imass['NH42SO4']/ammoniumsulfate.F_mass
            # made-up value based on online price, will change later
                                            
            ww.copy_like(influent) # ww has the same T and P as influent
            
            ww.imass['N'] = self.CHG.CHGout_N*(1 - ammonia_to_ammonium/(1 +\
                             ammonia_to_ammonium)*ammonia_transfer_ratio)
                                                              
            ww.imass['C'] = self.CHG.CHGout_C
                             
            ww.imass['P'] = self.CHG.CHGout_P
                       
            ww.imass['H2O'] -= (ww.imass['C'] + ww.imass['N'] + ww.imass['P'])
            
            ww.imass['H2O'] += self.ins[2].F_mass
            
            ammoniumsulfate.T = acid.T
            ammoniumsulfate.P = acid.P
            # ammoniumsulfate has the same T and P as acid
            
            mem_in.imass['Membrane'] = 0.15*self.membrane_area/yearly_operation_hour # m2/hr
            mem_out.copy_like(mem_in)
            mem_in.price = self.membrane_price
            # add membrane as streams to include 15% membrane replacement per year (Al-Obaidani)
    
    @property
    def N_recovery_ratio(self):
        return 1 - self.outs[1].imass['N']/self.CHG.CHGout_N
        
    def _design(self):
        if self.CHG.CHGout_N == 0:
            pass
        else:
            Design = self.design_results
            Design['Area'] = self.membrane_area
            Design['Total volume'] = Design['Area']*self.m2_2_m3
    
    def _cost(self):
        if self.CHG.CHGout_N == 0:
            pass
        else:
            Design = self.design_results
            purchase_costs = self.baseline_purchase_costs
            purchase_costs['Membrane'] = Design['Area']*self.membrane_price

# =============================================================================
# HT
# =============================================================================

class HT(Reactor):
    
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-products include fuel gas and
    char. The amount of biooil and fuel gas can be estimated using values from
    Li et al., 2018.
    The amount of char can be calculated based on mass closure.
    
    Model method: use experimental data.
    
    Design: tau assumes 2 times HTL tau. V_wf calculated based on
    LHSV = 0.5 L/L/hr and H2 and biocrude volume ratio.
    Void fraction of catalyst = 0.4 (the same as HC, Chemical Engineering Design).
    
    Parameters
    ----------
    ins: Iterable (stream)
    biocrude, hydrogen, catalyst_in
    outs: Iterable (stream)
    ht_out, catalyst_out
    '''
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _m3perhr_2_mmscfd = 1/1177.17
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst Jones
                 catalyst_lifetime=2*yearly_operation_hour, # 2 years Jones
                 catalyst_price=38.79,
                 hydrogen_P=1530*6894.76,
                 hydrogen_rxned_to_biocrude=0.046,
                 hydrocarbon_ratio=0.875, # 87.5 wt% of biocrude and reacted H2
                 # Jones et al., 2014
                 # spreadsheet HT calculation
                 HTin_T=174+273.15,
                 HTrxn_T=402+273.15, # Jones 2014
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
                                 'C26H42O4':0.01020, 'C30H62':0.00203},
                 # Jones et al., 2014
                 # spreadsheet HT calculation
                 # will not be a variable in uncertainty/sensitivity analysis
                 
                 P=None, tau=0.5, void_fraciton=0.4, # Towler
                 length_to_diameter=2, N=None, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_price = catalyst_price
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_biocrude = hydrogen_rxned_to_biocrude
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HTin_T = HTin_T
        self.HTrxn_T = HTrxn_T
        self.HT_composition = HT_composition
        
        IC_in = bst.Stream(f'{ID}_IC_in')
        IC_out = bst.Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        
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
        catalyst_in.price = self.catalyst_price
        
        hydrogen.imass['H2'] = biocrude.imass['Biocrude']*\
                               self.hydrogen_rxned_to_biocrude
        hydrogen.phase = 'g'
        hydrogen.price = 1.61

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
        
        if self.HTaqueous_C < -0.1*self.HTL.sludgelab.sludge_C:
            raise Exception('carbon mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance
        # should be no C in the aqueous phase, the calculation here is just for MB
        
        if self.HTaqueous_N < -0.1*self.HTL.sludgelab.sludge_N:
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
        Design['Hydrogen'] = self.ins[1].F_vol*self._m3perhr_2_mmscfd
        
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
        
        V_H2 = self.ins[1].F_vol*101325/self.hydrogen_P
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)

# =============================================================================
# HC
# =============================================================================

class HC(Reactor):
    
    '''
    Hydrocracking further cracks down heavy part in HT biooil to diesel and
    gasoline.
    
    Model method: use experimental data.
    
    Design: tau assumes 10 times HT tau. V_wf calculated based on
    LHSV = 0.5 L/L/hr and H2 and biocrude volume ratio.
    Void fraction of catalyst = 0.4 (Chemical Engineering Design).
    wall_thickness_factor = 1.5 to account for harsher conditions in HC than HT.
    (https://www.differencebetween.com/difference-between-hydrocracking-and-\
    hydrotreating/#:~:text=The%20key%20difference%20between%20hydrocracking,\
    processes%20in%20petroleum%20oil%20refining. (accessed 11-19-2022))
    
    Parameters
    ----------
    ins: Iterable (stream)
    heavy_oil, hydrogen, catalyst_in
    outs: Iterable (stream)
    ht_out, catalyst_out
    '''
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _m3perhr_2_mmscfd = 1/1177.17
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst Jones
                 catalyst_lifetime=5*yearly_operation_hour, # 5 years Jones
                 catalyst_price=38.79,
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_rxned_to_heavy_oil=0.01125,
                 hydrocarbon_ratio=1, # 100 wt% of heavy oil and reacted H2
                 # nearly all input heavy oils and H2 will be converted to
                 # products
                 # Jones et al., 2014
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
        
        SanUnit.__init__(self, ID, ins, outs, thermo,init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_price = catalyst_price
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_heavy_oil = hydrogen_rxned_to_heavy_oil
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HCin_T = HCin_T
        self.HCrxn_T = HCrxn_T
        self.HC_composition = HC_composition
        
        IC_in = bst.Stream(f'{ID}_IC_in')
        IC_out = bst.Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        
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
        catalyst_in.price = self.catalyst_price
        
        hydrogen.imass['H2'] = heavy_oil.F_mass*self.hydrogen_rxned_to_heavy_oil
        hydrogen.phase = 'g'
        hydrogen.price = 1.61

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
        Design['Hydrogen'] = self.ins[1].F_vol*self._m3perhr_2_mmscfd
        
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
        
        V_H2 = self.ins[1].F_vol*101325/self.hydrogen_P
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
    
    Model method: elements calculation, don't need _design and _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        supernatant_1, supernatant_2, memdis_ww, ht_ww
    outs: Iterable (stream)
        mixture
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 4
    _N_outs = 1
        
    def _run(self):
        
        supernatant_1, supernatant_2, memdis_ww, ht_ww = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        
        HT = self.ins[3]._source.ins[0]._source.ins[0]._source.ins[0]._source.\
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
    
    Model method: just correct phase, don't need _design and _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        influent
    outs: Iterable (stream)
        effluent
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', phase='l',
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
    
    Model method: Fake unit.
    
    Parameters
    ----------
    ins: Iterable (stream)
        gasoline, diesel
    outs: Iterable (stream)
        fuel
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', target='diesel',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target = target

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
        
        if self.target == 'diesel':
            fuel.imass['Diesel'] = diesel.F_mass + gasoline.F_mass*gasoline_LHV_2_diesel_LHV
            fuel.T = diesel.T
            fuel.P = fuel.P
            
    def _design(self):
        pass
    
    def _cost(self):
        pass