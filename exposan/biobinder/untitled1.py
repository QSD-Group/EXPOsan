# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 07:19:56 2024

@author: aliah
"""
import qsdsan as qs
from qsdsan import WasteStream, Unit
# Define a basic WasteStream object for feedstock input
feedstock = qs.WasteStream('feedstock')


# Set flow by concentrations
concentrations = {
    'lipids': 15,        # in kg/hr
    'proteins': 5,       # in kg/hr
    'carbohydrates': 10, # in kg/hr
    'ash': 5,            # in kg/hr
    'water': 65          # in kg/hr
}


# Set the flow using set_flow_by_concentration method
feedstock.set_flow_by_concentration(flow_tot=100, concentrations=concentrations, units='kg/hr')