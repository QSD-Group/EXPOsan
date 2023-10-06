
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Bright Elijah <be05055@georgiasouthern.edu>
    

Department of Civil Engineering and Construction
Georgia Southern University


This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.


'''
# %%


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
from qsdsan.sanunits._decay import Decay
import os

from math import ceil, pi
#from . import Decay
#from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path

__all__ = ('AgNP_CWF',)

#path to csv with all the inputs
#data_path += '\sanunit_data/_pou_chlorination.tsv'
agnp_path = ospath.join(data_path, 'sanunit_data/_AgNP_CWF_2.csv')


class AgNP_CWF(SanUnit):
    '''mh 
    Point of use water treatment technology: Desinfection through Silver Nanoparticles
    
    
    Reference documents
    ------------------- 
    N/A
    
    Parameters
    ----------
    ins : Raw water 
          AgNP
    outs : Treated water
    
    
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', number_of_households=1, **kwargs,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.number_of_households = number_of_households


      
        data = load_data(path=agnp_path)
     
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):

        raw_water,  = self.ins
        treated_water = self.outs[0]
        
        # give treated water all the properties and cmps of raw water
        # these will be changed below
        treated_water.copy_like(self.ins[0])
        
        
                
        #set the equation for log reduction following the Chick Watson Kinetic Model for log removal. Log(N/No) = -K*co*t
        #This was adopted from <https://pubs.acs.org/doi/full/10.1021/es4026084>
        #Here the log reduction value was set in the san sunit data (this wil be improved to accomodate changes in the raw water quality)
        
        No = raw_water.imass['Ecoli']/raw_water.F_vol * 10**-4  # E coli CFU/ mL ICC/ml (intact cells counts) used by (Cheswick et al., 2020)
        
        log_reduction = self.log_reduction_cwf
        
        log_N = log_reduction + np.log10(No) #CFU/mL
        N = 10**(log_N)
        
        #set conditional statement for simulation to capture the impact of raw water quality parameters

        self.hardness = raw_water.imass['Ca'] + raw_water.imass['Mg'] 
        
        if raw_water._turbidity <= 10 and self.hardness <= 60:
            self.AgNP_lifetime = 2
        elif raw_water._turbidity > 10 or self.hardness > 60:
            self.AgNP_lifetime = 1.5
        elif raw_water._turbidity > 10 and self.hardness > 60:
            self.AgNP_lifetime = 0.5
            



        
        
   
      
                            
        
     
    
  

             
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        design['CWFClay'] = CWF_clay = self.number_of_households * self.Clay_mass
        design['SilverNP'] =  AgNP = self.number_of_households * self.AgNP
        design['Sawdust'] =  Sawdust = self.number_of_households *self.Sawdust_mass
        design['PE'] = Container = self.number_of_households * self.PE_in_container/2
        

        self.construction = (
            Construction(item='CWFClay', quantity = CWF_clay, quantity_unit = 'kg'),
            Construction(item='SilverNP', quantity = AgNP, quantity_unit = 'kg', lifetime = self.AgNP_lifetime, lifetime_unit='yr'),
            Construction(item='Sawdust', quantity = Sawdust, quantity_unit = 'kg'),
            Construction(item='PE', quantity = Container, quantity_unit = 'kg')
            )
        self.add_construction(add_cost=False)
        
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        #breakpoint()
        #self.baseline_purchase_costs['cwf'] = (self.CWF_cost)*self.number_of_households
        self.baseline_purchase_costs['brush'] = (self.Brush_cost)*self.number_of_households
        self.baseline_purchase_costs['clay'] = (self.CWF_clay_cost)*self.number_of_households
        self.baseline_purchase_costs['grog'] = (self.CWF_grog_cost)*self.number_of_households
        self.baseline_purchase_costs['sawdust'] = (self.CWF_sawdust_cost)*self.number_of_households
        self.baseline_purchase_costs['water'] = (self.CWF_water_cost)*self.number_of_households
        self.baseline_purchase_costs['wood'] = (self.CWF_wood_cost)*self.number_of_households
        self.baseline_purchase_costs['labor'] = (self.CWF_labor_cost)*self.number_of_households
        self.baseline_purchase_costs['AgNP'] = (self.CWF_AgNP_cost*self.AgNP/self.Argenol_AgNP_content)*self.number_of_households
        self.baseline_purchase_costs['Bucket'] = (self.CWF_bucket)*self.number_of_households
        self.baseline_purchase_costs['Lid'] = (self.CWF_lid)*self.number_of_households
        self.baseline_purchase_costs['Spout'] = (self.CWF_spout)*self.number_of_households

        
        
        
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        # AgNP_lifetime = 0.1
        
        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with  the cost of the labor to replace them
        # USD/yr
        replacement_cost = (self.CWF_AgNP_cost*self.AgNP/self.Argenol_AgNP_content) * self.number_of_households / self.AgNP_lifetime
        
        #self.add_OPEX = self.chlorine_rate / self.NaClO_density / self.container_vol * self.operator_refill_cost  # USD/hr (all items are per hour)
        self.add_OPEX = replacement_cost / (365*24)
        
        
        


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Bright Elijah <be05055@georgiasouthern.edu & brightcarlelijah@gmail.com>
    
Department of Civil Engineering and Construction
Georgia Southern University

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.


'''
# %%


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
from qsdsan.sanunits._decay import Decay
import os

from math import ceil, pi
from . import Decay
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path

__all__ = ('POUChlorination',)

#path to csv with all the inputs
#data_path += '\sanunit_data/_pou_chlorination.tsv'
poucl_path = ospath.join(data_path, 'sanunit_data/_pou_chlorination.csv')


class POUChlorination(SanUnit):
    '''mh 
    Point of use water treatment technology: Desinfection through chlorination
    
    Reference documents
    ------------------- 
    Hussein, M.; Brown, J.; Njee, R. M.; Clasen, T.; Malebo, H. M.; Mbuligwe, S. Point-of-Use Chlorination of Turbid Water: Results from a Field Study in Tanzania. Journal of Water and Health 2015, 13 (2), 544–552. http://dx.doi.org/10.2166/wh.2014.001.
    
    Tamene, A. A Qualitative Analysis of Factors Influencing Household Water Treatment Practices Among Consumers of Self-Supplied Water in Rural Ethiopia. Risk Management and Healthcare Policy 2021, 14, 1129–1139. https://doi.org/10.2147/RMHP.S299671.
    
    Parameters
    ----------
    ins : Raw water 
        Chlorine 
    outs : Treated water
    
    
    References
    ----------
    .. N/A
    
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', number_of_households=1, **kwargs,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.number_of_households = number_of_households

# load data from tsv each name will be self.name  

        data = load_data(path=poucl_path)
     
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
# define the number of influent and effluent streams    
    _N_ins = 3
    _N_outs = 1

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):

        raw_water, chlorine, Cl_bottle = self.ins
        treated_water = self.outs[0]
        
        # give treated water all the properties and cmps of raw water
        # these will be changed below
        treated_water.copy_like(self.ins[0])
        chlorine.phase = 'l'
        Cl_bottle.phase = 's'
        
        # add chlorine to the treated_water
        #breakpoint()

        if raw_water.turbidity <= 10:
            pass
        elif raw_water.turbidity > 10:
            self.NaClO_dose *= 2
      
        
        chlorine.imass['NaClO'] = self.NaClO_dose * treated_water.F_vol / 1000 # kg NaClO/hr
        self.chlorine_rate = chlorine.imass['NaClO']
        Cl_bottle.imass['Polyethylene'] = chlorine.imass['NaClO'] * self.PE_to_NaClO
        
        # disinfect bacteria from treated_water
        
        Cl_Concentration = self.NaClO_dose #mg/L
      
        No = raw_water.imass['Ecoli']/raw_water.F_vol * 10**-4  # E coli CFU/ mL ICC/ml (intact cells counts) used by (Cheswick et al., 2020)
             
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items are in the _impacts_items.xlsx data
        design['PE'] = Container = self.number_of_households * self.PE_in_container


        self.construction = (
            Construction(item='PE', quantity = Container, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)       
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['WaterContainer'] = (self.container_cost)*self.number_of_households
       
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        

        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with  the cost of the labor to replace them

        #self.add_OPEX = self.chlorine_rate / self.NaClO_density / self.container_vol * self.operator_refill_cost  # USD/hr (all items are per hour)



#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Bright Elijah <be05055@georgiasouthern.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.


'''
# %%


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
from qsdsan.sanunits._decay import Decay
import os

from math import ceil, pi
#from . import Decay
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path

__all__ = ('POU_UV',)

#path to csv with all the inputs
#data_path += '\sanunit_data/_pou_chlorination.tsv'
pou_uv_path = ospath.join(data_path, 'sanunit_data/_pou_uv.csv')


class POU_UV(SanUnit):
    '''mh 
    Point of use water treatment technology: Desinfection through POU UV
    

    
    Reference documents
    ------------------- 
    N/A
    
    Parameters
    ----------
    ins : Raw water 
        
    outs : Treated water
    
        

        
    References
    ----------
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', number_of_households=1, **kwargs,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.number_of_households = number_of_households

# load data from tsv each name will be self.name  
 

        data = load_data(path=pou_uv_path)
     
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)


        
# define the number of influent and effluent streams  
# There is one influent streams which is  the raw water (from the raw water san unit)
# The effluent stream is the treated water
  
    _N_ins = 1
    _N_outs = 1

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
     
        raw_water,  = self.ins
        treated_water = self.outs[0]
        
        # give treated water all the properties and cmps of raw water
        # these will be changed below
        treated_water.copy_like(self.ins[0])
        
        if raw_water.F_vol > (self.uv_flow*60):
            breakpoint()
            
        self.run_time = raw_water.F_vol/(self.uv_flow*60)
        

            
    
            
    
        
        

    
      
        No = raw_water.imass['Ecoli']/raw_water.F_vol * 10**-4  # E coli CFU/ mL ICC/ml (intact cells counts) used by (Cheswick et al., 2020)


        # if self.raw_water.uv_turbidty <= 18 :
        #     log_removal =  5.62
        # else:
        #     print("high turbidity")
        
        #set the equation for log reduction following the Chick Watson Kinetic Model for log removal. Log(N/No) = -K*co*t
        #This was adopted from <https://escholarship.org/uc/item/3p76b9gb>
        # Chatterley, C.; Linden, K. Demonstration and Evaluation of Germicidal UV-LEDs for Point-of-Use Water Disinfection. Journal of Water and Health 2010, 8 (3), 479–486. https://doi.org/10.2166/wh.2010.124.
        #Here the log reduction value was set in the san sunit data
        
        log_removal =  self.uv_slope*self.uv_dose + self.uv_intercept
        log_N = log_removal + np.log10(No) #CFU/mL
        N = 10**(log_N)
        # N = np.exp(log_red  uction + np.log10(No))
         
        ################ Log (N0/N) = Kd × UV dose where kd = inactivation rate constant
        if raw_water._turbidity <= 10:
            self.lamp_lifespan_factor = 1
        elif raw_water._turbidity > 10:
            self.lamp_lifespan_factor = 0.5
        
        self.lamp_life_span *= self.lamp_lifespan_factor
        
        ## determine the flowrate
        
        # if raw_water.UVT <= 80:
        #     self.flowrate = 11.1 #ml/min
        # elif raw_water.UVT > 80:
        #     self.flowrate = 20   #ml/min
        # else:
        #     pass
        
        # ## determine the time
        
        # self.residence_time = (self.unit_volume*1000) / self.flowrate /60    #time in hr
        
        ## factor this time into power demand cost and into lamp replacement lifetime.
        
        
        
        
                
    

             
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        design['PE'] = uv_storage = self.number_of_households * self.storage_PE*2
        design['PVC'] = uv_pvc = self.number_of_households * self.uv_PVC 
        design['Uvlamp'] = uv_lamp_mecury = self.number_of_households * self.number_of_uv_lamps
        design['Aluminum'] =  uv_aluminum_foil = self.number_of_households * self.uv_aluminum_foil
        

        self.construction = (
            Construction(item='PE', quantity = uv_storage, quantity_unit = 'kg'),
            Construction(item='PVC', quantity = uv_pvc, quantity_unit = 'kg'),
            Construction(item='Uvlamp', quantity = uv_lamp_mecury, quantity_unit = 'kg', lifetime = (self.lamp_life_span ), lifetime_unit='hr'),
            Construction(item='Aluminum', quantity = uv_aluminum_foil, quantity_unit = 'kg')
            )
        self.add_construction(add_cost=False)
        
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        #breakpoint()
        self.baseline_purchase_costs['uv_unit'] = (self.uv_unit_cost)*self.number_of_households
        self.baseline_purchase_costs['uv_storage'] = (self.uv_storage_cost*2)*self.number_of_households
        self.add_OPEX['uvlamp'] = (self.uv_lamp_cost/(self.lamp_life_span)) 
         
        
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        power_demand = (self.uv_electric_demand / 1000) * self.run_time / self.lamp_lifespan_factor
        self.power_utility(power_demand)
        


        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with  the cost of the labor to replace them

        #self.add_OPEX = self.chlorine_rate / self.NaClO_density / self.container_vol * self.operator_refill_cost  # USD/hr (all items are per hour)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Bright Elijah <be05055@georgiasouthern.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.


'''
# %%


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
from qsdsan.sanunits._decay import Decay
import os

from math import ceil, pi
#from . import Decay
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path

__all__ = ('UV_LED',)

#path to csv with all the inputs
#data_path += '\sanunit_data/_pou_chlorination.tsv'
uv_led_path = ospath.join(data_path, 'sanunit_data/_uv_led.csv')


class UV_LED(SanUnit):
    '''mh 
    Point of use water treatment technology: Desinfection through UV LED
    

    
    Reference documents
    ------------------- 
    N/A
    
    Parameters
    ----------
    ins : Raw water 
        
    outs : Treated water
    
        

        
    References
    ----------
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', number_of_households=1, **kwargs,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.number_of_households = number_of_households

# load data from tsv each name will be self.name  
 
        data = load_data(path = uv_led_path)
     
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)


        
# define the number of influent and effluent streams  
# There is one influent streams which is  the raw water (from the raw water san unit)
# The effluent stream is the treated water
  
    _N_ins = 1
    _N_outs = 1

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
     
        raw_water,  = self.ins
        treated_water = self.outs[0]
        
        # give treated water all the properties and cmps of raw water
        # these will be changed below
        treated_water.copy_like(self.ins[0])
        
        
    
      
        No = raw_water.imass['Ecoli']/raw_water.F_vol * 10**-4  # E coli CFU/ mL ICC/ml (intact cells counts) used by (Cheswick et al., 2020)


        # if self.raw_water.uv_turbidty <= 18 :
        #     log_removal =  5.62
        # else:
        #     print("high turbidity")
        
        #set the equation for log reduction following the Chick Watson Kinetic Model for log removal. Log(N/No) = -K*co*t
        #This was adopted from <https://pubmed.ncbi.nlm.nih.gov/26179637/>
        # Chatterley, C.; Linden, K. Demonstration and Evaluation of Germicidal UV-LEDs for Point-of-Use Water Disinfection. Journal of Water and Health 2010, 8 (3), 479–486. https://doi.org/10.2166/wh.2010.124.
        #Here the log reduction value was set in the san sunit data

        log_removal =  self.uv_led_slope*self.uv_led_dose + self.uv_led_intercept
        log_N = log_removal + np.log10(No) #CFU/mL
        N = 10**(log_N)
        # N = np.exp(log_reduction + np.log10(No))
        
        if raw_water._turbidity <= 10:
            self.led_lifespan_factor = 1
        elif raw_water._turbidity > 10:
            self.led_lifespan_factor = 0.5
            
        self.uv_led_lifespan *= self.led_lifespan_factor
        
        
        
        if raw_water.F_vol > (self.uv_led_flow*60):
            breakpoint()
            
        self.run_time = raw_water.F_vol/(self.uv_led_flow*60)
        
        
        ## determine the flowrate
        
        # if raw_water.UVT <= 80:
        #     self.flowrate = 11.1 #ml/min
        # elif raw_water.UVT > 80:
        #     self.flowrate = 20   #ml/min
        # else:
        #     pass
        
        # ## determine the time
        
        # self.residence_time = (self.unit_volume*1000) / self.flowrate /60    #time in hr
        
        ## factor this time into power demand cost and into lamp replacement lifetime.
         

             
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        design['Quartz'] = uv_led_quartz = self.number_of_households * self.uv_quartz
        design['PE'] = uv_led_storage = self.number_of_households * self.uv_led_storage_PE*2
        design['StainlessSteel'] =  uv_led_steel = self.number_of_households * self.StainlessSteel
        design['LED'] = uv_led = self.number_of_households * self.uv_led_weight

        

        self.construction = (
            Construction(item='Quartz', quantity = uv_led_quartz, quantity_unit = 'kg'),
            Construction(item='PE', quantity = uv_led_storage, quantity_unit = 'kg'),
            Construction(item='StainlessSteel', quantity = uv_led_steel, quantity_unit = 'kg'),
            Construction(item='LED', quantity = uv_led, quantity_unit = 'kg', lifetime = (self.uv_led_lifespan), lifetime_unit='hr'),
            )
        self.add_construction(add_cost=False)
        
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        #breakpoint()
        self.baseline_purchase_costs['uv_led_unit'] = (self.uv_led_unit_cost)*self.number_of_households
        self.baseline_purchase_costs['uv_led_storage'] = (self.uv_led_storage_cost*2)*self.number_of_households
        #self.baseline_purchase_costs['uv_pump'] = (self.uv_led_pump_cost)*self.number_of_households
        self.add_OPEX['LED'] = (self.uv_led_cost/(self.uv_led_lifespan)) 
         
        
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        
        power_demand = (self.led_electricity_demand / 1000) * self.run_time /self.led_lifespan_factor
        self.power_utility(power_demand)        

        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with  the cost of the labor to replace them

        #self.add_OPEX = self.chlorine_rate / self.NaClO_density / self.container_vol * self.operator_refill_cost  # USD/hr (all items are per hour)