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
    Quantitative Evaluation of an Integrated System for Valorization of Wastewater 
    Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
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
    Process Design and Economics for the Conversion of Algal Biomass to Hydrocarbons: 
    Whole Algae Hydrothermal Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.
'''

import biosteam as bst, qsdsan as qs
from qsdsan import SanUnit
from biosteam.units import StorageTank
from qsdsan.sanunits import HXutility

__all__ = (
    'SludgeLab',
    'HTL',
    'HT',
    'HC',
    'AcidExtraction',
    'HTLmixer',
    'StruvitePrecipitation',
    'CHG',
    'MembraneDistillation',
    'GasMixer')


class SludgeLab(SanUnit):
    '''
    SludgeLab is a fake unit that makes it easier to incorperate sludge compositions
    and moisture content into uncertainty analysis.
    
    Parameters
    ----------
    ins: Iterable (stream)
        fake_sludge
    outs: Iterable (stream)
        real_sludge
    '''

    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream', 
                 sludge_moisture=0.99,sludge_dw_protein=0.341,
                 sludge_dw_lipid=0.226,sludge_dw_carbo=0.167,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.sludge_moisture=sludge_moisture
        self.sludge_dw_protein=sludge_dw_protein
        self.sludge_dw_carbo=sludge_dw_carbo
        self.sludge_dw_lipid=sludge_dw_lipid
        
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        fake_sludge=self.ins[0]
        real_sludge=self.outs[0]
        real_sludge.imass['H2O']=fake_sludge.F_mass*self.sludge_moisture
        real_sludge.imass['Sludge_protein']=fake_sludge.F_mass*(1-self.sludge_moisture)*\
            self.sludge_dw_protein
        real_sludge.imass['Sludge_carbo']=fake_sludge.F_mass*(1-self.sludge_moisture)*\
            self.sludge_dw_carbo
        real_sludge.imass['Sludge_lipid']=fake_sludge.F_mass*(1-self.sludge_moisture)*\
            self.sludge_dw_lipid
        real_sludge.imass['Sludge_ash']=fake_sludge.F_mass*(1-self.sludge_moisture)*\
            (1-self.sludge_dw_protein-self.sludge_dw_carbo-self.sludge_dw_lipid)

        
    def _design(self):
        pass
    
    def _cost(self):
        pass

class HTL(SanUnit):
    
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and biochar under
    elevated temperature (350°C) and pressure. The products percentage (wt%) can be evaluated
    using revised MCA model (Li et al., 2017, Leow et al., 2018) with known sludge
    composition (protein%, lipid%, and carbohydrate%)
    
    Parameters
    ----------
    ins: Iterable (stream)
        dewatered_sludge
    outs: Iterable (stream)
        biochar, HTLaqueous, biocrude, offgas
    '''
    auxiliary_unit_names=('heat_exchanger',)

    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream', 
                 biocrude_moisture_content=0.044,sludge_C_ratio=0.411,
                 sludge_P_ratio=0.019,sludge_H_ratio=0.058,sludge_S_ratio=0.01,
                 sludge_N_ratio=0.05,sludge_O_ratio=0.261,
                 biochar_C_N_ratio=15.5,biochar_C_P_ratio=2.163,
                 eff_T=60+273.15, #Jones 2014
                 **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.biocrude_moisture_content=biocrude_moisture_content
        self.sludge_C_ratio=sludge_C_ratio
        self.sludge_P_ratio=sludge_P_ratio
        self.sludge_H_ratio=sludge_H_ratio
        self.sludge_S_ratio=sludge_S_ratio
        self.sludge_N_ratio=sludge_N_ratio
        self.sludge_O_ratio=sludge_O_ratio
        self.biochar_C_N_ratio=biochar_C_N_ratio
        self.biochar_C_P_ratio=biochar_C_P_ratio
        self.eff_T = eff_T
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out, T=eff_T)


    _N_ins = 1
    _N_outs = 4
    

    def _run(self):
        dewatered_sludge=self.ins[0]
        biochar,HTLaqueous,biocrude,offgas= outs = self.outs
        for s in outs: s.T = dewatered_sludge.T
        
        
        self.dewatered_sludge_afdw = dewatered_sludge_afdw =\
            dewatered_sludge.imass['Sludge_lipid']\
            + dewatered_sludge.imass['Sludge_protein']\
            + dewatered_sludge.imass['Sludge_carbo']
            
        self.dewatered_sludge_dw = dewatered_sludge_afdw + dewatered_sludge.imass['Sludge_ash']
                               
        lipid_ratio = dewatered_sludge.imass['Sludge_lipid']/dewatered_sludge_afdw
        protein_ratio = dewatered_sludge.imass['Sludge_protein']/dewatered_sludge_afdw
        carbo_ratio = dewatered_sludge.imass['Sludge_carbo']/dewatered_sludge_afdw
        
        biocrude.imass['Biocrude'] = (0.846*lipid_ratio + 0.445*protein_ratio\
            + 0.205*carbo_ratio) * dewatered_sludge_afdw
        biochar.imass['Biochar'] = 0.377*carbo_ratio * dewatered_sludge_afdw     
        HTLaqueous.imass['HTLaqueous'] = (0.154*lipid_ratio + 0.481*protein_ratio)\
            * dewatered_sludge_afdw #HTLaqueous is TDS in aqueous phase
        offgas.imass['CO2'] = (0.074*protein_ratio + 0.418*carbo_ratio)\
            * dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1-self.biocrude_moisture_content)-\
            biocrude.imass['Biocrude']
        HTLaqueous.imass['H2O'] = dewatered_sludge.imass['H2O']-biocrude.imass['H2O']+\
            dewatered_sludge.imass['Sludge_ash'] #assume ash goes to water
        
        biochar.phase='s'
        offgas.phase='g'
        for stream in outs: stream.T = self.eff_T
        
        biochar.P = 3029.7*6894.76 #Jones 2014: 3029.7 psia
        HTLaqueous.P = 30*6894.76 #Jones 2014: 30 psia
        biocrude.P = 30*6894.76
        offgas.P = 30*6894.76

    
    def _design(self):
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.outs)
        hx_ins0.T = self.ins[0].T
        hx_outs0.copy_flow(hx_ins0)
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
    
    @property
    def AOSc(self):
       return (3*self.sludge_N_ratio/14+2*self.sludge_O_ratio/16-\
               self.sludge_H_ratio/1)/(self.sludge_C_ratio/12)
    
    @property
    def sludge_carbo_ratio_dw(self):
        return self.ins[0].imass['Sludge_carbo']/self.dewatered_sludge_dw
    
    @property
    def sludge_protein_ratio_dw(self):
        return self.ins[0].imass['Sludge_protein']/self.dewatered_sludge_dw
    
    @property
    def biochar_C_ratio(self):
        return min(1.75*self.sludge_carbo_ratio_dw,0.65)
    
    @property
    def biochar_N_ratio(self):
        return self.biochar_C_ratio/self.biochar_C_N_ratio
    
    @property
    def biochar_P_ratio(self):
        return self.biochar_C_ratio/self.biochar_C_P_ratio

    @property
    def biochar_C(self):
        return self.outs[0].F_mass*self.biochar_C_ratio

    @property
    def biochar_N(self):
        return self.outs[0].F_mass*self.biochar_N_ratio

    @property
    def biochar_P(self):
        return min((self.ins[0].F_mass-self.ins[0].imass['H2O'])*self.sludge_P_ratio,
                   self.outs[0].F_mass*self.biochar_P_ratio)

    @property
    def biocrude_C_ratio(self):
        return (self.AOSc*(-8.37)+68.55)/100

    @property
    def biocrude_N_ratio(self):
        return 0.133*self.sludge_protein_ratio_dw

    @property
    def biocrude_C(self):
        return self.outs[2].F_mass*self.biocrude_C_ratio

    @property
    def biocrude_N(self):
        return self.outs[2].F_mass*self.biocrude_N_ratio

    @property
    def offgas_C(self):
        return self.outs[3].F_mass*12/44

    @property
    def HTLaqueous_C(self):
        return self.dewatered_sludge_dw*self.sludge_C_ratio -\
            self.biochar_C - self.biocrude_C - self.offgas_C

    @property
    def HTLaqueous_N(self):
        return self.dewatered_sludge_dw*self.sludge_N_ratio -\
            self.biochar_N - self.biocrude_N

    @property
    def HTLaqueous_P(self):
        return self.dewatered_sludge_dw*self.sludge_P_ratio -\
            self.biochar_P

    def _design(self):
        pass
    
    def _cost(self):
        pass


class AcidExtraction(SanUnit):
    
    '''
    H2SO4 is added to biochar from HTL to extract P. 
    
    Parameters
    ----------
    ins: Iterable (stream)
        biochar, acid
    outs: Iterable (stream)
        residual, extracted
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',acid_vol=10,
                 P_acid_recovery_ratio=0.95,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.acid_vol=acid_vol
        self.P_acid_recovery_ratio=P_acid_recovery_ratio

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        biochar, acid = self.ins
        residual, extracted = self.outs
        
        acid.imass['H2SO4']=biochar.F_mass*self.acid_vol*0.5*98/1000 #0.5 M H2SO4 acid_vol (10 mL/1 g) Biochar
        acid.imass['H2O']=biochar.F_mass*self.acid_vol*1.05-acid.imass['H2SO4'] #0.5 M H2SO4 density: 1.05 kg/L 
        #https://www.fishersci.com/shop/products/sulfuric-acid-1n-0-5m-standard-solution-thermo-
        #scientific/AC124240010 (accessed 10-6-2022)
        acid.T = biochar.T
        acid.P = self.ins[0]._source.outs[1].P
        
        residual.imass['Residual']=biochar.F_mass*\
            (1-self.ins[0]._source.biochar_P_ratio*self.P_acid_recovery_ratio)
        
        extracted.copy_like(acid)
        extracted.imass['P']=biochar.F_mass-residual.F_mass
        
        residual.phase='s'
        
        residual.T = biochar.T
        extracted.T = acid.T
        
    @property
    def residual_C(self):
        return self.ins[0]._source.biochar_C
    
    @property
    def residual_N(self):
        return self.ins[0]._source.biochar_N

    @property
    def residual_P(self):
        return self.ins[0]._source.biochar_P-self.outs[1].imass['P']
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
    
class HTLmixer(SanUnit):
    '''
    
    Parameters
    ----------
    ins: Iterable (stream)
        HTLaqueous, extracted
    outs: Iterable (stream)
        mixture
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',
                 **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]
        
        mixture.imass['C']=self.ins[0]._source.HTLaqueous_C
        mixture.imass['N']=self.ins[0]._source.HTLaqueous_N
        mixture.imass['P']=self.ins[0]._source.HTLaqueous_P+extracted.imass['P']
        mixture.imass['H2O']=HTLaqueous.F_mass+extracted.F_mass-\
            mixture.imass['C']-mixture.imass['N']-mixture.imass['P'] #Represented by H2O except C, N, P
        
        mixture.T = extracted.T
        mixture.P = HTLaqueous.P
        
    def _design(self):
        pass
    
    def _cost(self):
        pass


class StruvitePrecipitation(SanUnit):
    '''
    extracted_P and HTL aqueous are mixed together (Mixer) before adding MgCl2 and struvite precipitation.
    
    Parameters
    ----------
    ins: Iterable (stream)
        mixture, supply_MgCl2
    outs: Iterable (stream)
        struvite, effluent
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',Mg_P_ratio=1,
                 P_pre_recovery_ratio=0.95,P_in_struvite=0.127,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.Mg_P_ratio=Mg_P_ratio
        self.P_pre_recovery_ratio=P_pre_recovery_ratio
        self.P_in_struvite=P_in_struvite

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        mixture, supply_MgCl2 = self.ins
        struvite, effluent = self.outs
        
        supply_MgCl2.imass['MgCl2']=mixture.imass['P']/31*95.211*self.Mg_P_ratio
        struvite.imass['Struvite']=mixture.imass['P']*self.P_pre_recovery_ratio/self.P_in_struvite
        
        supply_MgCl2.phase='s'
        
        effluent.copy_like(mixture)
        effluent.imass['P'] -= struvite.imass['Struvite']*self.P_in_struvite
        effluent.imass['N'] -= struvite.imass['Struvite']*self.P_in_struvite/31*14
        effluent.imass['H2O'] += (supply_MgCl2.imass['MgCl2']-\
                                       struvite.imass['Struvite']*(1-self.P_in_struvite*(1+14/31)))
        
        struvite.phase='s'    
            
        struvite.T = mixture.T
        effluent.T = mixture.T
        
    @property
    def struvite_P(self):
        return self.outs[0].imass['Struvite']*self.P_in_struvite

    @property
    def struvite_N(self):
        return self.struvite_P*14/31

    def _design(self):
        pass
    
    def _cost(self):
        pass
    

class CHG(SanUnit):
   
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel gas 
    under elevated temperature (350°C) and pressure.
    
    Parameters
    ----------
    ins: Iterable (stream)
        CHGfeed
    outs: Iterable (stream)
        CHGfuelgas, effluent
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',ch4_ratio=0.244,
                 co_ratio=0.029,co2_ratio=0.150,c2h6_ratio=0.043,toc_tc_ratio=0.764,
                 toc_to_gas_c_ratio=0.262,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.ch4_ratio=ch4_ratio
        self.co_ratio=co_ratio
        self.co2_ratio=co2_ratio
        self.c2h6_ratio=c2h6_ratio
        self.toc_tc_ratio=toc_tc_ratio
        self.toc_to_gas_c_ratio=toc_to_gas_c_ratio
        self.heat_exchanger = HXutility(ID=None, ins=(None,), outs=None, thermo=self.thermo)
        self._mixed_higherT = qs.Stream(f'{self.ID}_mixed_higherT')
        self._mixed_lowerT = qs.Stream(f'{self.ID}_mixed_lowerT')
        
    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        CHGfeed = self.ins[0]
        CHGfuelgas, effluent = outs = self.outs
        
        for s in outs: s.T = CHGfeed.T
        
        for i in (self.ch4_ratio,self.co_ratio,self.co2_ratio,self.c2h6_ratio,
                  1-self.ch4_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio):
            if i < 0: i = 0
        
        CHGfuel_gas_mass = CHGfeed.imass['C']*self.toc_tc_ratio*self.toc_to_gas_c_ratio/(self.\
            ch4_ratio*12/16+self.co_ratio*12/28+self.co2_ratio*12/44+self.c2h6_ratio*24/30)

        CHGfuelgas_composition = {
            'CH4':self.ch4_ratio,
            'CO':self.co_ratio,
            'CO2':self.co2_ratio,
            'C2H6':self.c2h6_ratio,
            'H2':1-self.ch4_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio
            }
        
        for name,ratio in CHGfuelgas_composition.items():
            CHGfuelgas.imass[name] = CHGfuel_gas_mass * ratio
        
        effluent.copy_like(CHGfeed)
        effluent.imass['C']*=(1-self.toc_tc_ratio*self.toc_to_gas_c_ratio)
        effluent.imass['H2O']=CHGfeed.F_mass-CHGfuel_gas_mass-effluent.imass['C']-\
            effluent.imass['N']-effluent.imass['P']
        
        CHGfuelgas.phase='g'
        
        mixed_higherT, mixed_lowerT = self._mixed_higherT, self._mixed_lowerT
        mixed_higherT.mix_from(outs)

        CHGfuelgas.T = 60+273.15 #Jones 2014
        effluent.T = 60+273.15
        
        mixed_lowerT.mix_from(outs)

        hx = self.heat_exchanger
        hx.simulate_as_auxiliary_exchanger(inlet=mixed_higherT, outlet=mixed_lowerT)
        
        CHGfuelgas.P = 50*6894.76 #Jones 2014: 50 psia
        effluent.P = 50*6894.76
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    

class MembraneDistillation(SanUnit):
    
    '''
    Membrane distillation will be modeled using vapor-liquid-liquid equilibrium later.
    Here is simplified.
    
    Parameters
    ----------
    ins: Iterable (stream)
        influent, acid
    outs: Iterable (stream)
        ammoniasulfate, ww
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',N_S_ratio=2,
                 N_recovery_rate=0.825,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.N_S_ratio=N_S_ratio
        self.N_recovery_rate=N_recovery_rate

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        influent, acid = self.ins
        ammoniasulfate, ww = self.outs
        
        acid.imass['H2SO4'] = influent.imass['N']/14/self.N_S_ratio*98
        acid.imass['H2O'] = acid.imass['H2SO4']*1000/98/0.5*1.05-acid.imass['H2SO4']
        ammoniasulfate.imass['NH42SO4']=influent.imass['N']*self.N_recovery_rate/28*132
        ammoniasulfate.imass['H2O']=acid.imass['H2O']
        ww.copy_like(influent)
        ww.imass['N']*=(1-self.N_recovery_rate)
        ww.imass['H2O']+=ww.imass['N']*self.N_recovery_rate+acid.imass['H2SO4']-\
            ammoniasulfate.imass['NH42SO4']
        ammoniasulfate.T = acid.T
        ww.T = acid.T
        ww.P = acid.P
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    

class HT(SanUnit):
    
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C) and pressure
    to produce upgraded biooil. Co-products include fuel gas and char. The amount of
    biooil and fuel gas can be estimated using values from Li et al., 2018.
    The amount of char can be calculated based on mass closure.
    
    Parameters
    ----------
    ins: Iterable (stream)
    biocrude,hydrogen_gas
    outs: Iterable (stream)
    HTaqueous, fuel_gas, gasoline, diesel, heavy_oil
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',
                 biooil_ratio=0.77,gas_ratio=0.07,#Jones et al., 2014
                 biocrude_h2_ratio = 18.14,
                 gasoline_ratio = 0.103, heavy_oil_ratio = 0.113,
                 co_ratio=0.128,co2_ratio=0.007,
                 c2h6_ratio=0.188,c3h8_ratio=0.107,c4h10_ratio=0.09,
                 biooil_C_ratio=0.855,biooil_N_ratio=0.01,
                 **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.biooil_ratio=biooil_ratio
        self.gas_ratio=gas_ratio
        self.biocrude_h2_ratio=biocrude_h2_ratio
        self.gasoline_ratio=gasoline_ratio
        self.heavy_oil_ratio=heavy_oil_ratio  
        self.c4h10_ratio=c4h10_ratio
        self.co_ratio=co_ratio
        self.co2_ratio=co2_ratio
        self.c2h6_ratio=c2h6_ratio
        self.c3h8_ratio=c3h8_ratio
        self.biooil_C_ratio = biooil_C_ratio
        self.biooil_N_ratio = biooil_N_ratio
        self.heat_exchanger = HXutility(ID=None, ins=(None,), outs=None, thermo=self.thermo)
        self._mixed_higherT = qs.Stream(f'{self.ID}_mixed_higherT')
        self._mixed_lowerT = qs.Stream(f'{self.ID}_mixed_lowerT')

    _N_ins = 2
    _N_outs = 5
        
    def _run(self):
        
        biocrude,hydrogen_gas = self.ins
        HTaqueous,fuel_gas,gasoline,diesel,heavy_oil = outs = self.outs
        
        for s in outs: s.T = biocrude.T
        
        hydrogen_gas.imass['H2'] = biocrude.F_mass/self.biocrude_h2_ratio
        hydrogen_gas.phase='g'
        hydrogen_gas.T = 117.5+273.15 #Jones 2014
        hydrogen_gas.P = 1530.0*6894.76
        
        biooil_total_mass = biocrude.F_mass*self.biooil_ratio
        
        for i in (self.gasoline_ratio,self.heavy_oil_ratio):
            if i<0: i=0
            
        gasoline.imass['Gasoline'] = biooil_total_mass*self.gasoline_ratio
        diesel.imass['Diesel'] = biooil_total_mass*(1-self.gasoline_ratio-self.heavy_oil_ratio)
        heavy_oil.imass['Heavy_oil'] = biooil_total_mass*self.heavy_oil_ratio
        
        gas_mass = biocrude.F_mass*self.gas_ratio
        
        for i in (self.c4h10_ratio,self.co_ratio,self.co2_ratio,self.c2h6_ratio,self.c3h8_ratio,
                    1-self.c4h10_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio-\
                    self.c3h8_ratio):
            if i<0: i=0
    
        fuelgas_HT_composition = {
            'CO':self.co_ratio,
            'CO2':self.co2_ratio,
            'C2H6':self.c2h6_ratio,
            'C3H8':self.c3h8_ratio,
            'C4H10':self.c4h10_ratio,
            'CH4':1-self.c4h10_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio-\
                self.c3h8_ratio
            }
            
        for name,ratio in fuelgas_HT_composition.items():
            fuel_gas.imass[name] = gas_mass * ratio
        
        HTaqueous.imass['HTaqueous'] = biocrude.F_mass+hydrogen_gas.F_mass-biooil_total_mass-gas_mass
        #HTaqueous is liquid waste from HT
        
        fuel_gas.phase='g'
        
        mixed_higherT, mixed_lowerT = self._mixed_higherT, self._mixed_lowerT
        mixed_higherT.mix_from(outs)

        HTaqueous.T = 47+273.15 #Jones 2014
        fuel_gas.T = 44+273.15
        gasoline.T = 109+273.15
        diesel.T = 265+273.15
        heavy_oil.T = 381+273.15
        
        mixed_lowerT.mix_from(outs)

        hx = self.heat_exchanger
        hx.simulate_as_auxiliary_exchanger(inlet=mixed_higherT, outlet=mixed_lowerT)
        
        HTaqueous.P = 55*6894.76 #Jones 2014: 55 psia
        fuel_gas.P = 20*6894.76
        gasoline.P = 25*6894.76
        diesel.P = 18.7*6894.76
        heavy_oil.P = 18.7*6894.76
        
    @property
    def biooil_C(self):
        return min(self.outs[2].F_mass*self.biooil_C_ratio,self.ins[0]._source.biocrude_C)
        
    @property
    def HTfuel_gas_C(self):
        HTfuel_gas_C = 0
        fuelgas_carbo_ratio = {
            'C4H10':60/72,
            'CO':12/28,
            'CO2':12/44,
            'C2H6':24/30,
            'C3H8':36/44,
            'CH4':12/16
              }
        for name, ratio in fuelgas_carbo_ratio.items():
            HTfuel_gas_C+=self.outs[1].imass[name]*ratio
        return min(HTfuel_gas_C,self.ins[0]._source.biocrude_C-self.biooil_C)
    
    @property
    def biooil_N(self):
        return self.outs[2].F_mass*self.biooil_N_ratio

    biocrude_C_ratio = HTL.biocrude_C_ratio

    @property
    def HTaqueous_C(self):
        return max(0,self.ins[0]._source.biocrude_C - self.HTfuel_gas_C - self.biooil_C)

    @property
    def HTaqueous_N(self):
        return self.ins[0]._source.biocrude_N - self.biooil_N

    def _design(self):
        pass
    
    def _cost(self):
        pass
   
   
class HC(SanUnit):
    
    '''
    Hydrocracking further cracks down heavy part in HT biooil to diesel and gasoline.
    
    Parameters
    ----------
    ins: Iterable (stream)
    heavy_oil, hydrogen_gas
    outs: Iterable (stream)
    gasoline, diesel, off_gas
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',
                 heavy_oil_h2_ratio=45.17,gasoline_hc_ratio=0.303,off_gas_ratio=0.0394,
                 **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.heavy_oil_h2_ratio=heavy_oil_h2_ratio
        self.gasoline_hc_ratio=gasoline_hc_ratio
        self.off_gas_ratio=off_gas_ratio
        self.heat_exchanger = HXutility(ID=None, ins=(None,), outs=None, thermo=self.thermo)
        self._mixed_higherT = qs.Stream(f'{self.ID}_mixed_higherT')
        self._mixed_lowerT = qs.Stream(f'{self.ID}_mixed_lowerT')
        
    _N_ins = 2
    _N_outs = 3
        
    def _run(self):
        
        heavy_oil, hydrogen_gas = self.ins
        gasoline, diesel, off_gas = outs = self.outs
        
        for s in outs: s.T = heavy_oil.T
        
        hydrogen_gas.imass['H2'] = heavy_oil.F_mass/self.heavy_oil_h2_ratio
        hydrogen_gas.phase='g'
        hydrogen_gas.T = 128+273.15 #Jones 2014
        hydrogen_gas.P = 1039.7*6894.76
        
        for i in [self.gasoline_hc_ratio,self.off_gas_ratio]:
            if i < 0: i=0
        
        gasoline.imass['Gasoline'] = (heavy_oil.F_mass+hydrogen_gas.F_mass)*self.gasoline_hc_ratio

        diesel.imass['Diesel'] = (heavy_oil.F_mass+hydrogen_gas.F_mass)*\
            (1-self.gasoline_hc_ratio-self.off_gas_ratio)
        
        off_gas.imass['CO2'] = (heavy_oil.F_mass+hydrogen_gas.F_mass)*self.off_gas_ratio
        
        off_gas.phase='g'
        
        mixed_higherT, mixed_lowerT = self._mixed_higherT, self._mixed_lowerT
        mixed_higherT.mix_from(outs)

        gasoline.T = 60+273.15 #Jones 2014
        diesel.T = 220+273.15
        off_gas.T = 45+273.15
        
        mixed_lowerT.mix_from(outs)

        hx = self.heat_exchanger
        hx.simulate_as_auxiliary_exchanger(inlet=mixed_higherT, outlet=mixed_lowerT)
        
        gasoline.P = 20*6894.76 #Jones 2014: 20 psia
        diesel.P = 24*6894.76
        off_gas.P = 20*6894.76
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
   

class FuelMixer(SanUnit):
    
    '''
    Mix fuels.
    
    Parameters
    ----------
    ins: Iterable (stream)
        ht_fuel, hc_fuel
    outs: Iterable (stream)
        mixed_fuel
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        
        ht_fuel, hc_fuel = ins = self.ins
        mixed_fuel = self.outs[0]
        
        mixed_fuel.mix_from(ins)
        # mixed_fuel.copy_like(ht_fuel)
        # mixed_fuel.imass['Gasoline'] += hc_fuel.imass['Gasoline']
        # mixed_fuel.imass['Diesel'] += hc_fuel.imass['Diesel']
        
    def _design(self):
        pass
    
    def _cost(self):
        pass 

    
class GasMixer(SanUnit):
    
    '''
    Mix all gases.
    
    Parameters
    ----------
    ins: Iterable (stream)
        htl_off_gas, ht_fuel_gas, hc_off_gas, chg_fuel_gas
    outs: Iterable (stream)
        fuel_gas
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='Stream',**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

    _N_ins = 4
    _N_outs = 1
        
    def _run(self):
        
        htl_off_gas, ht_fuel_gas, hc_off_gas, chg_fuel_gas = ins = self.ins
        fuel_gas = self.outs[0]
        
        fuel_gas.mix_from(ins)
        fuel_gas.phase='g'
        # fuel_gas.copy_like(ht_fuel_gas)
        # for gas in ['CH4','CO','CO2','C2H6','H2','C3H8','C4H10']:
        #     fuel_gas.imass[gas] += (htl_off_gas.imass[gas]+hc_off_gas.imass[gas]+
        #                             chg_fuel_gas.imass[gas])
            
    def _design(self):
        pass
    
    def _cost(self):
        pass