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
    
(1) Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
    Quantitative Evaluation of an Integrated System for Valorization of
    Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
    Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
    https://doi.org/10.1021/acs.est.8b04035.
    
(2) Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
    Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
    Liquefaction Products from Feedstock Biochemical Composition.
    Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    
(3) Snowden-Swan, L. J.; Zhu, Y.; Jones, S. B.; Elliott, D. C.; Schmidt, A. J.; 
    Hallen, R. T.; Billing, J. M.; Hart, T. R.; Fox, S. P.; Maupin, G. D. 
    Hydrothermal Liquefaction and Upgrading of Municipal Wastewater Treatment 
    Plant Sludge: A Preliminary Techno-Economic Analysis; 
    PNNL--25464, 1258731; 2016; https://doi.org/10.2172/1258731.

(4) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to
    Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
    PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.
    
(5) Hao, S.; Choi, Y.-J.; Wu, B.; Higgins, C. P.; Deeb, R.; Strathmann, T. J.
    Hydrothermal Alkaline Treatment for Destruction of Per- and Polyfluoroalkyl
    Substances in Aqueous Film-Forming Foam. Environ. Sci. Technol.
    2021, 55(5), 3283–3295. https://doi.org/10.1021/acs.est.0c06906.
    
(6)	Matayeva, A.; Rasmussen, S. R.; Biller, P. Distribution of Nutrients and
    Phosphorus Recovery in Hydrothermal Liquefaction of Waste Streams.
    BiomassBioenergy 2022, 156, 106323.
    https://doi.org/10.1016/j.biombioe.2021.106323.
'''

import biosteam as bst
from qsdsan import SanUnit
from qsdsan.sanunits import HXutility
from biosteam.units import IsothermalCompressor
from exposan.htl._components_design import create_components
from biosteam.units.design_tools import PressureVessel
from biosteam.units.design_tools.cost_index import CEPCI_by_year as CEPCI
from math import pi, ceil
from biosteam.exceptions import DesignError
from biosteam import Stream
from biosteam.units.decorators import cost

__all__ = ('SludgeLab',
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
          'PhaseChanger')

cmps = create_components()







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
              'Reactor volume': 'm3'}

    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _Vmax = pi/4*(20**2)*40/35.3147
    
    _F_BM_default = PressureVessel._F_BM_default

    def __init__(self, ID='', ins=None, outs=(), *,
                 P=101325, tau=0.5, V_wf=0.8,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):

        SanUnit.__init__(self, ID, ins, outs)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _design(self):
        Design = self.design_results
        
        ins_F_vol = self.F_vol_in
        for i in range(len(self.ins)):
            ins_F_vol -= self.ins[i].ivol['H2']
        # not include gas (e.g. H2)
        V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor

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
        fake_sludge
    outs: Iterable (stream)
        real_sludge
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', 
                 sludge_moisture=0.99, sludge_dw_protein=0.341,
                 sludge_dw_lipid=0.226, sludge_dw_carbo=0.167,
                 sludge_P_ratio = 0.019, # data are from SS PNNL 2021
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_protein = sludge_dw_protein
        self.sludge_dw_carbo = sludge_dw_carbo
        self.sludge_dw_lipid = sludge_dw_lipid
        self.sludge_dw_ash = 1 - sludge_dw_protein - sludge_dw_carbo -\
                             sludge_dw_lipid
        self.sludge_P_ratio = sludge_P_ratio
        # set P as an independent variable, assume S (sulfur) is 0
        
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        
        fake_sludge = self.ins[0]
        real_sludge = self.outs[0]
        
        real_sludge.imass['H2O'] = fake_sludge.F_mass * self.sludge_moisture
        sludge_dw = fake_sludge.F_mass*(1 - self.sludge_moisture)
        real_sludge.imass['Sludge_protein'] = sludge_dw*self.sludge_dw_protein
        real_sludge.imass['Sludge_carbo'] = sludge_dw*self.sludge_dw_carbo
        real_sludge.imass['Sludge_lipid'] = sludge_dw*self.sludge_dw_lipid
        real_sludge.imass['Sludge_ash'] = sludge_dw*self.sludge_dw_ash
          
    # all sludge elemental analysis are based on empirical equation
    @property
    def sludge_C_ratio(self):
       return self.sludge_dw_carbo*0.44 + self.sludge_dw_lipid*0.75 +\
           self.sludge_dw_protein*0.53
    # https://pubmed.ncbi.nlm.nih.gov/2061559/ (accessed 2022-10-27)
    # https://encyclopedi/a2.thefreedictionary.com/Proteins
    # (accessed 2022-10-27)
    # add uncertainty after serious calibration!
    
    @property
    def sludge_H_ratio(self):
       return self.sludge_C_ratio*0.143
    # based on SS PNNL 2021 data, H ~ C/7
    # add uncertainty after serious calibration!
   
    @property
    def sludge_N_ratio(self):
       return self.sludge_dw_protein*0.16 
    # or change according to Leow 2018: afdw pr/4.78
    # but here, we use dw instead of afdw
    # https://www.fao.org/3/y5022e/y5022e03.htm#:~:text=On%20the%20basis%20of%
    # 20early,is%20confounded%20by%20two%20considerations (accessed 2022-10-27)
    # add uncertainty after serious calibration!
   
    @property
    def sludge_P_ratio(self):
       return self._sludge_P_ratio
    # set P as an indepedent variable since hard to find any association with
    # sludge biochemical compositions
    # add uncertainty after serious calibration!
    
    @sludge_P_ratio.setter
    def sludge_P_ratio(self, i):
        if not 0 <= i <= 1:
            raise AttributeError('`sludge_P` must be within [0, 1], '
                                f'the provided value {i} is outside the\
                                range.')
        self._sludge_P_ratio = i
    
    @property
    def sludge_O_ratio(self):
       return 1 - self.sludge_C_ratio - self.sludge_H_ratio -\
           self.sludge_N_ratio - self.sludge_P_ratio - self.sludge_dw_ash*0.75
    # sludge_O is calculated based on mass balance closure and * 0.75 since
    # double count some elements. 0.75 is based on SS PNNL 2021
    # add uncertainty after serious calibration!
    
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
    
    Parameters
    ----------
    ins: Iterable (stream)
        dewatered_sludge
    outs: Iterable (stream)
        biochar, HTLaqueous, biocrude, offgas
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
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
                 NaOH_2_water = 0.2, # Hao: 5 M NaOH
                 biocrude_moisture_content=0.056, # Jones
                 biochar_P_ratio=0.86, # Matayeva: 0.84-0.88
                 gas_composition={'CH4':0.050, 'C2H6':0.032,
                                  'CO2':0.918}, # Jones
                 biochar_pre=3029.7*6894.76, # Jones
                 HTLaqueous_pre=30*6894.76, # Jones
                 biocrude_pre=30*6894.76, # Jones
                 offgas_pre=30*6894.76, # Jones
                 eff_T=60+273.15, # Jones
                 
                 P=None, tau=1, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
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
        self.NaOH_2_water = NaOH_2_water
        self.biocrude_moisture_content = biocrude_moisture_content
        self.biochar_P_ratio = biochar_P_ratio
        self.gas_composition = gas_composition
        self.biochar_pre = biochar_pre
        self.HTLaqueous_pre = HTLaqueous_pre
        self.biocrude_pre = biocrude_pre
        self.offgas_pre = offgas_pre
        self.eff_T = eff_T
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    _N_ins = 2
    _N_outs = 4
    _units= {'Treatment capacity': 'lb/h'} # separator
    
    def _run(self):
        
        dewatered_sludge, base = self.ins
        biochar, HTLaqueous, biocrude, offgas = self.outs
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid'] +\
                                dewatered_sludge.imass['Sludge_protein'] +\
                                dewatered_sludge.imass['Sludge_carbo']
        # just use afdw in revised MCA model, other places use dw
        
        # NaOH added here. target is 5 M (Shilai Hao EST)
        # for water solution: 5 M NaOH: 20 g NaOH / 100 mL H2O
        # (0.2 kg NaOH / 1 kg H2O)
        # here, the calculation is based on the water amount in the dewatered
        # sludge (assume the initial pH = 7, and solids don't affect pH)
        base.imass['NaOH'] = dewatered_sludge.imass['H2O']*self.NaOH_2_water
        
        lipid_ratio = dewatered_sludge.imass['Sludge_lipid']/\
                      dewatered_sludge_afdw
        protein_ratio = dewatered_sludge.imass['Sludge_protein']/\
                        dewatered_sludge_afdw
        carbo_ratio = dewatered_sludge.imass['Sludge_carbo']/\
                      dewatered_sludge_afdw

        # the following calculations are based on revised MCA model
        biochar.imass['Biochar'] = 0.377*carbo_ratio*dewatered_sludge_afdw  
        
        HTLaqueous.imass['HTLaqueous'] = (0.154*lipid_ratio +\
                                          0.481*protein_ratio)*\
                                          dewatered_sludge_afdw
        # HTLaqueous is TDS in aqueous phase
         
        gas_mass = (self.protein_2_gas*protein_ratio + self.carbo_2_gas*carbo_ratio)*\
                       dewatered_sludge_afdw
                       
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass*ratio
            
        biocrude.imass['Biocrude'] = (self.lipid_2_biocrude*lipid_ratio + self.protein_2_biocrude*protein_ratio\
                                      + self.carbo_2_biocrude*carbo_ratio)*\
                                      dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
                                
        HTLaqueous.imass['H2O'] = dewatered_sludge.imass['H2O'] -\
                                  biocrude.imass['H2O'] +\
                                  dewatered_sludge.imass['Sludge_ash'] +\
                                  base.imass['NaOH']
        # assume ash (all soluble based on Jones) goes to water
        # all NaOH also goes to water to maintain pH for membrane distillation
        
        biochar.phase = 's'
        offgas.phase = 'g'
        
        biochar.P = self.biochar_pre
        HTLaqueous.P = self.HTLaqueous_pre
        biocrude.P = self.biocrude_pre
        offgas.P = self.offgas_pre
        
        for stream in self.outs: stream.T = self.eff_T
        
        self.sludgelab = self.ins[0]._source.ins[0]._source.ins[0].\
                         _source.ins[0]._source.ins[0]._source
        
        # in the case that HTLaqueous_C or HTLaqueous_N is less than 0, which
        # is not likely to happen, we add the following two exceptions.
        if self.HTLaqueous_C < 0:
            raise Exception('double check sludge composition, HTLaqueous_C '\
                            'is not likely to be less than 0, otherwise, we '\
                            'do not need CHG')
        if self.HTLaqueous_N < 0:
            raise Exception('double check sludge composition, HTLaqueous_N '\
                            'is not likely to be less than 0, otherwise, we '\
                            'do not need membrane distillation')
            
        # self.HTLaqueous_P is always larger than 0, since an constraint is
        # added when calculated biochar_P

    @property
    def biochar_C_ratio(self):
        return min(self.biochar_C_slope*self.sludgelab.sludge_dw_carbo, 0.65)
    # revised MCA model

    @property
    def biochar_C(self):
        return self.outs[0].F_mass*self.biochar_C_ratio

    @property
    def biochar_P(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*\
                self.sludgelab.sludge_P_ratio*self.biochar_P_ratio
        
    @property
    def biorude_P(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*\
                self.sludgelab.sludge_P_ratio*(1 - self.biochar_P_ratio)

    @property
    def biocrude_C_ratio(self):
        return (self.sludgelab.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100 # revised MCA model
    
    @property
    def biocrude_H_ratio(self):
        return (self.sludgelab.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100 # revised MCA model

    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.sludgelab.sludge_dw_protein # revised MCA model
    
    @property
    def biocrude_C(self):
        return self.outs[2].F_mass*self.biocrude_C_ratio

    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio

    @property
    def biocrude_N(self):
        return self.outs[2].F_mass*self.biocrude_N_ratio
    
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
    def offgas_C(self):
        carbon = 0
        for name in self.gas_composition.keys():
            carbon += self.outs[3].imass[name]*cmps[name].i_C
        return carbon   
               
    # C and N in aqueous phase are calculated based on mass balance closure
    @property
    def HTLaqueous_C(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*\
                self.sludgelab.sludge_C_ratio - self.biochar_C -\
                self.biocrude_C - self.offgas_C

    @property
    def HTLaqueous_N(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*\
                self.sludgelab.sludge_N_ratio - self.biocrude_N

    def _design(self):
        
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass*self._kg_2_lb
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx.T = hx_outs0.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

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
    
    Parameters
    ----------
    ins: Iterable (stream)
        biochar, acid
    outs: Iterable (stream)
        residual, extracted
    '''
    
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', acid_vol=10, P_acid_recovery_ratio=0.95,
                 
                 P=None, tau=1, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_vol = acid_vol
        self.P_acid_recovery_ratio = P_acid_recovery_ratio
        
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
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
        
        acid.imass['H2SO4'] = biochar.F_mass*self.acid_vol*0.5*98.079/1000
        #0.5 M H2SO4 acid_vol (10 mL/1 g) Biochar
        acid.imass['H2O'] = biochar.F_mass*self.acid_vol*1.05 -\
                            acid.imass['H2SO4']
        # 0.5 M H2SO4 density: 1.05 kg/L 
        # https://www.fishersci.com/shop/products/sulfuric-acid-1n-0-5m-
        # standard-solution-thermo-scientific/AC124240010 (accessed 10-6-2022)
        
        residual.imass['Residual'] = biochar.F_mass*(1 - self.ins[0].
                                     _source.biochar_P_ratio*self.
                                     P_acid_recovery_ratio)
        
        extracted.copy_like(acid)
        extracted.imass['P'] = biochar.F_mass - residual.F_mass
        # assume just P can be extracted
        
        residual.phase = 's'
        
        residual.T = extracted.T = biochar.T
        # H2SO4 reacts with biochar to release heat and temperature will be
        # increased mixture's temperature
        
    @property
    def residual_C(self):
        return self.ins[0]._source.biochar_C
    
    @property
    def residual_N(self):
        return self.ins[0]._source.biochar_N

    @property
    def residual_P(self):
        return self.ins[0]._source.biochar_P - self.outs[1].imass['P']
        
    def _design(self):
        
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
        
        mixture.imass['C'] = self.ins[0]._source.HTLaqueous_C
        mixture.imass['N'] = self.ins[0]._source.HTLaqueous_N
        mixture.imass['P'] = self.ins[0]._source.HTLaqueous_P +\
                             extracted.imass['P']
        mixture.imass['H2O'] = HTLaqueous.F_mass + extracted.F_mass -\
                               mixture.imass['C'] - mixture.imass['N'] -\
                               mixture.imass['P']
        # represented by H2O except C, N, P
        
        mixture.T = extracted.T
        mixture.P = extracted.P
        
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
    
    Parameters
    ----------
    ins: Iterable (stream)
        mixture, supply_MgCl2, supply_NH4Cl
    outs: Iterable (stream)
        struvite, effluent
    '''
    
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', Mg_P_ratio=1,
                 P_pre_recovery_ratio=0.95, P_in_struvite=0.127,
                 
                 P=None, tau=1, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.Mg_P_ratio = Mg_P_ratio
        self.P_pre_recovery_ratio = P_pre_recovery_ratio
        self.P_in_struvite = P_in_struvite
        
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    _N_ins = 3
    _N_outs = 2
        
    def _run(self):
        
        mixture, supply_MgCl2, supply_NH4Cl = self.ins
        struvite, effluent = self.outs
        
        if mixture.imass['P']/30.973762 > mixture.imass['N']/14.0067:
            supply_NH4Cl.imass['NH4Cl'] = (mixture.imass['P']/30.973762 -\
                                           mixture.imass['N']/14.0067)*53.491
        # make sure N:P >= 1:1
        
        supply_MgCl2.imass['MgCl2'] = mixture.imass['P']/30.973762*95.211*\
                                      self.Mg_P_ratio # Mg:P = 1:1
        struvite.imass['Struvite'] = mixture.imass['P']*\
                                     self.P_pre_recovery_ratio/\
                                     self.P_in_struvite
        supply_MgCl2.phase = 's'
        
        effluent.copy_like(mixture)
        effluent.imass['P'] -= struvite.imass['Struvite']*self.P_in_struvite
        effluent.imass['N'] += supply_NH4Cl.imass['NH4Cl']*14.0067/53.491 -\
                               struvite.imass['Struvite']*\
                               self.P_in_struvite/30.973762*14.0067
        effluent.imass['H2O'] += (supply_MgCl2.imass['MgCl2'] +\
                                  supply_NH4Cl.imass['NH4Cl'] -\
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
        
        self.P = self.ins[0].P
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)

# =============================================================================
# CHG
# =============================================================================

class CHG(Reactor):
   
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    
    Model method: use experimental data, assume no NH3 loss for now.
    
    Parameters
    ----------
    ins: Iterable (stream)
        CHGfeed
    outs: Iterable (stream)
        CHGfuelgas, effluent
    '''
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 gas_composition={'CH4':0.527,
                                  'CO2':0.432,
                                  'C2H6':0.011,
                                  'C3H8':0.030,
                                  'H2':0.0001},
                 # Jones
                 # will not be a variable in uncertainty/sensitivity analysis
                 gas_c_to_total_c=0.764*0.262, # Li EST
                 # Jones 2014: pressure before flash
                 
                 P=None, tau=1, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.gas_composition = gas_composition
        self.gas_c_to_total_c = gas_c_to_total_c
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    _N_ins = 1
    _N_outs = 1
        
    def _run(self):
        
        CHGfeed = self.ins[0]
        CHGout = self.outs[0]
        
        gas_C_ratio = 0
        for name, ratio in self.gas_composition.items():
            gas_C_ratio += ratio*cmps[name].i_C
        
        gas_mass = CHGfeed.imass['C']*self.gas_c_to_total_c/gas_C_ratio
        
        for name,ratio in self.gas_composition.items():
            CHGout.imass[name] = gas_mass*ratio
            
        CHGout.imass['H2O'] = CHGfeed.F_mass - gas_mass
        # all C, N, and P are accounted in H2O here, but will be calculated as
        # properties.
            
        CHGout.T = CHGfeed.T
        
        CHGout.P = CHGfeed.P
        
    @property
    def CHGout_C(self):
        # not include carbon in gas phase
        return self.ins[0].imass['C']*(1 - self.gas_c_to_total_c)
    
    @property
    def CHGout_N(self):
        return self.ins[0].imass['N']
    
    @property
    def CHGout_P(self):
        return self.ins[0].imass['P']
        
    def _design(self):
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.outs)
        hx_outs0.mix_from(self.outs)
        hx_ins0.T = self.ins[0].T # temperature before/after CHG are similar
        hx.T = hx_outs0.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        Reactor._Vmax /= 2 # so that there are 6 CHG tanks
        self.P = self.ins[0].P
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)
        purchase_costs = self.baseline_purchase_costs
        current_cost = 0 # cost w/o sulfur guard
        for item in purchase_costs.keys():
            current_cost += purchase_costs[item]
        purchase_costs['sulfur guard'] = current_cost*0.05
    
# =============================================================================
# Membrane Distillation
# =============================================================================

class MembraneDistillation(Reactor):
    
    '''
    Membrane distillation recovers nitrogen as ammonia sulfate based on vapor
    pressure difference across the hydrophobic membrane.
    
    Model method: 
        1. Feed pH = 10, permeate pH = 1.5 (0.5 M H2SO4)
        2. All N in the feed are NH4+/NH3 (Jones PNNL 2014)
        3. 95% NH3 in feed can be transfered to permeate (assume 95% for now,
           use literature data to find a conservative assumpation later)
        4. All NH3 in permeate can form (NH4)2SO4 (which makes sense since
           just water evaporates)
        5. _design and _cost refer to
           A.A. et al., Membrane distillation: A comprehensive review

    Parameters
    ----------
    ins: Iterable (stream)
        influent, acid
    outs: Iterable (stream)
        ammoniasulfate, ww
    '''
    
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 N_S_ratio=2, pH=10, ammonia_transfer_ratio=0.95,
                 
                 P=None, tau=1, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.N_S_ratio = N_S_ratio
        self.pH = pH
        self.ammonia_transfer_ratio = ammonia_transfer_ratio
        
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        influent, acid = self.ins
        ammoniumsulfate, ww = self.outs
        
        influent.imass['C'] = self.ins[0]._source.ins[0]._source.ins[0].\
                              _source.CHGout_C
        influent.imass['N'] = self.ins[0]._source.ins[0]._source.ins[0].\
                              _source.CHGout_N
        influent.imass['P'] = self.ins[0]._source.ins[0]._source.ins[0].\
                              _source.CHGout_P
        influent.imass['H2O'] -= (influent.imass['C'] + influent.imass['N'] +\
                                  influent.imass['P'])
        
        acid.imass['H2SO4'] = influent.imass['N']/14.0067/self.N_S_ratio*98.079
        acid.imass['H2O'] = acid.imass['H2SO4']*1000/98.079/0.5*1.05 -\
                            acid.imass['H2SO4']
        
        pKa = 9.26 # ammonia pKa
        ammonia_to_ammonium = 10**(-pKa)/10**(-self.pH)
        ammonia_in_feed = influent.imass['N']/14.0067*ammonia_to_ammonium/(1 +\
                          ammonia_to_ammonium)*17.031

        ammoniumsulfate.imass['NH42SO4'] = ammonia_in_feed*\
                                          self.ammonia_transfer_ratio/34.062*\
                                          132.14
        ammoniumsulfate.imass['H2O'] = acid.imass['H2O']
        ammoniumsulfate.imass['H2SO4'] = acid.imass['H2SO4'] +\
                                         ammoniumsulfate.imass['NH42SO4']/\
                                         132.14*28.0134 -\
                                         ammoniumsulfate.imass['NH42SO4']
                                        
        ww.copy_like(influent) # ww has the same T and P as influent
        ww.imass['N'] -= influent.imass['N']*self.ammonia_transfer_ratio*\
                         ammonia_to_ammonium/(1 + ammonia_to_ammonium)
                         
        ammoniumsulfate.T = acid.T
        ammoniumsulfate.P = acid.P
        # ammoniumsulfate has the same T and P as acid
        
    def _design(self):
        
        self.P = self.ins[1].P
        Reactor._design(self)
    
    def _cost(self):
        Reactor._cost(self)

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
    
    Parameters
    ----------
    ins: Iterable (stream)
    biocrude, hydrogen
    outs: Iterable (stream)
    HTaqueous, fuel_gas, gasoline, diesel, heavy_oil
    '''
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _m3perhr_2_mmscfd = 1/1177.17
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 hydrogen_P=1530*6894.76,
                 hydrogen_to_biocrude=0.138,
                 hydrogen_rxned_to_biocrude=0.046,
                 hydrocarbon_ratio=0.875,
                 # Jones et al., 2014
                 # spreadsheet HT calculation
                 HTin_T=174+273.15,
                 HTrxn_T=402+273.15, # Jones 2014
                 HT_composition={'CH4':0.0228, 'C2H6':0.0292,
                                 'C3H8':0.0165, 'C4H10':0.0087,
                                 'TWOMBUTAN':0.0041, 'NPENTAN':0.0068,
                                 'TWOMPENTA':0.0041, 'HEXANE':0.0041,
                                 'TWOMHEXAN':0.0041, 'HEPTANE':0.0041,
                                 'CC6METH':0.0102, 'PIPERDIN':0.0041,
                                 'TOLUENE':0.0102, 'THREEMHEPTA':0.0102,
                                 'OCTANE':0.0102, 'ETHCYC6':0.0041,
                                 'ETHYLBEN':0.0204, 'OXYLENE':0.0102,
                                 'C9H20':0.0041, 'PROCYC6':0.0041,
                                 'C3BENZ':0.0102, 'FOURMONAN':0,
                                 'C10H22':0.0204, 'C4BENZ':0.0122,
                                 'C11H24':0.0204, 'C10H12':0.0204,
                                 'C12H26':0.0204, 'OTTFNA':0.0102,
                                 'C6BENZ':0.0204, 'OTTFSN':0.0204,
                                 'C7BENZ':0.0204, 'C8BENZ':0.0204,
                                 'C10H16O4':0.0184, 'C15H32':0.0612,
                                 'C16H34':0.1836, 'C17H36':0.0816, 
                                 'C18H38':0.0408, 'C19H40':0.0408,
                                 'C20H42':0.1020, 'C21H44':0.0408,
                                 'TRICOSANE':0.0408, 'C24H38O4':0.0082,
                                 'C26H42O4':0.0102, 'C30H62':0.0020},
                 # Jones et al., 2014
                 # spreadsheet HT calculation
                 # will not be a variable in uncertainty/sensitivity analysis
                 
                 P=None, tau=2, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.hydrogen_P = hydrogen_P
        self.hydrogen_to_biocrude = hydrogen_to_biocrude
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
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        
        biocrude, hydrogen = self.ins
        ht_out = self.outs[0]
        
        hydrogen.imass['H2'] = biocrude.imass['Biocrude']*\
                               self.hydrogen_to_biocrude
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
        
        self.HTL = self.ins[0]._source.ins[0]._source
        
        if self.HTaqueous_C < -0.1*self.HTL.biocrude_C:
            raise Exception('carbon mass balance is out of +/- 10%')
        # allow +/- 10% out of mass balance
        
        if self.HTaqueous_N < -0.1*self.HTL.biocrude_N:
            raise Exception('nitrogen mass balance is out of +/- 10%')
        # allow +/- 10% out of mass balance
        
        # possibility exist that more carbon is in biooil and gas than in
        # biocrude because we use the biooil/gas compositions to calculate
        # carbon. In this case, the C in HT aqueous phase will be negative.
        # It's OK if the mass balance is within +/- 10%. Otherwise, an
        # exception will be raised.
        
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
        hx_outs0.mix_from(self.ins)
        hx_outs0.T = self.HTin_T
        hx_ins0.P = hx_outs0.P = min(IC_outs0.P, self.ins[0].P)
        # H2 and biocrude have the same pressure
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
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
    
    Parameters
    ----------
    ins: Iterable (stream)
    heavy_oil, hydrogen
    outs: Iterable (stream)
    gasoline, diesel, off_gas
    '''
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _m3perhr_2_mmscfd = 1/1177.17
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_to_heavy_oil=0.0625,
                 hydrogen_rxned_to_heavy_oil=0.01125,
                 hydrocarbon_ratio=1,
                 # nearly all input heavy oils and H2 will be converted to
                 # products
                 # Jones et al., 2014
                 # spreadsheet HC calculation
                 HCin_T=394+273.15,
                 HC_rxn_T=451+273.15,
                 HC_composition={'CO2':0.0388, 'CH4':0.0063,
                                 'CYCHEX':0.0389,'HEXANE':0.0116,
                                 'HEPTANE':0.1202, 'OCTANE':0.0851,
                                 'C9H20':0.0952, 'C10H22':0.1231,
                                 'C11H24':0.1764, 'C12H26':0.1382,
                                 'C13H28':0.0974, 'C14H30':0.0486,
                                 'C15H32':0.0340, 'C16H34':0.0201,
                                 'C17H36':0.0045, 'C18H38':0.0010,
                                 'C19H40':0.0052, 'C20H42':0.0003},
                 #combine C20H42 and PHYTANE as C20H42
                 # will not be a variable in uncertainty/sensitivity analysis
                 
                 P=None, tau=2, V_wf=0.5,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo,init_with)
        self.hydrogen_P = hydrogen_P
        self.hydrogen_to_heavy_oil = hydrogen_to_heavy_oil
        self.hydrogen_rxned_to_heavy_oil = hydrogen_rxned_to_heavy_oil
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HCin_T = HCin_T

        self.HC_rxn_T = HC_rxn_T

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
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        
        heavy_oil, hydrogen = self.ins
        hc_out = self.outs[0]
        
        hydrogen.imass['H2'] = heavy_oil.F_mass*self.hydrogen_to_heavy_oil
        hydrogen.phase = 'g'

    
        hydrocarbon_mass = heavy_oil.F_mass*(1 +\
                           self.hydrogen_rxned_to_heavy_oil)*\
                           self.hydrocarbon_ratio

        for name, ratio in self.HC_composition.items():
            hc_out.imass[name] = hydrocarbon_mass*ratio
        
        hc_out.imass['H2'] = hydrogen.imass['H2'] - heavy_oil.F_mass*\
                             self.hydrogen_rxned_to_heavy_oil
        
        hc_out.P = heavy_oil.P
        hc_out.T = self.HC_rxn_T
        
        C_in = 0
        total_num = len(list(cmps))
        for num in range(total_num):
            C_in += heavy_oil.imass[str(list(cmps)[num])]*list(cmps)[num].i_C
            
        C_out = self.hydrocarbon_C
        
        if C_out < 0.95*C_in or C_out > 1.05*C_out :
            raise Exception('carbon mass balance is out of +/- 5%')
        # make sure that carbon mass balance is within +/- 10%. Otherwise, an
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
        hx_outs0.mix_from(self.ins)
        hx_outs0.T = self.HCin_T
        hx_ins0.P = hx_outs0.P = min(IC_outs0.P, self.ins[0].P)
        # H2 and biocrude have the same pressure
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
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
             ins[0]._source
        
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