#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import biosteam as bst
from qsdsan import SanUnit,Construction
from qsdsan.sanunits import SludgeThickening
from qsdsan.utils import ospath, data_path, load_data #!!!not sure if ... is used correctly
vr_su_data_path = ospath.join(data_path, 'sanunit_data/vr')

#%%

vr_filter_press_path = ospath.join(vr_su_data_path, '_vr_filter_press.tsv')

class VolumeReductionFilterPress(SludgeThickening): 
    #!!!SludgeThickening or SludgeSeparator?
    '''
    A Filter Press unit for the dewatering of mixed excreta in volume reduction
    generation II reinveted toilet [1]
    
    The 0th outs is the water-rich supernatant (effluent) and
    the 1st outs is the solid-rich sludge.
    
    The following componenet should be included in system thermo object for simulation:
    Water.
    
    The following impact items should be pre-constructed for life cycle assessment:
    Stainless steel.
    
    Parameters
    ----------
    ins: Iterable(stream)
      Pasteurized solids waste for dewatering treatment
    outs: Iterable(stream)
      Liquids and solids produced from filter press.
    sludge_moisture: float
      Moisture content of the solids cake after filter press [wt% water].  
    References
    -----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
    https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1
    
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', 
                 ins=None, outs=(), thermo=None, 
                 init_with='WasteStream',
                 sludge_moisture = 0.45, 
                 solids = (),
                 **kwargs):
        SludgeThickening.__init__(self, ID, ins, outs, thermo, init_with,
                                sludge_moisture=sludge_moisture,
                                solids=solids,
                                **kwargs)
                
        data = load_data(path=vr_filter_press_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _init_lca(self): #!!!use "_init_lca" or "_refres_lca"?
        self.construction = [
            Construction('stainless steel',linked_unit=self, 
                         item='Stainless steel', 
                         quantity_unit='kg'),
            ]
    
    def _run(self): #!!! not sure if this is needed, super class already has all
        waste = self.ins
        liq, cake_sol = self.outs
        SludgeThickening._run(self)
        
    def _design(self):
        design = self.design_results #!!!not sure where is "design_result" from.
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.filterpress_ss_weight
        self.add_construction(add_cost = False) #!!!not sure if add cost or not, if so, how to check ID inventory?
        
    def _cost(self):
        C = self.baseline_purchase_costs  #!!! where is baseline_purchase_costs from?
        C["Filter Press"] = self.filterpress_purchase_cost * self.price_ratio #USD
        
        self.power_utility(self.filterpress_energy_persolids*self.outs[1].F_mass) # kW
        
        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost() #USD/hr
        
    def _calc_replacement_cost(self): #USD/hr
        filter_media_replacement_cost = (
               self.filterpress_filter_media_price*self.filterpress_filter_media_area/self.filterpress_filter_media_lifetime         
            )
        return filter_media_replacement_cost / (365*24) * self.price_ratio
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        filter_press_maintenance_labor_cost= (self.filterpress_maintenance * self.wages)
        return filter_press_maintenance_labor_cost / (365*24) * self.price_ratio
        
