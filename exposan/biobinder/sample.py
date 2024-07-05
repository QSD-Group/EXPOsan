# -*- coding: utf-8 -*-

import biosteam as bst
import qsdsan as qs

from qsdsan import Component
from qsdsan import sanunits as qsu
from qsdsan import WasteStream, Unit

# Define components
from qsdsan import Component

# Define biocrude oil component
biocrude = qs.Component('biocrude', formula='C10H20O2', phase='l', particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)


# Define NaCl (salt) component
nacl = qs.Component('nacl', formula='NaCl', phase='s', particle_size = 'Particulate', degradability = 'Biological',
                     organic = False)

# Define water component
water = qs.Component('water', formula='H2O', phase='l', particle_size = 'Particulate', degradability = 'Biological',
                     organic = False)

# Define lipids component
lipids = qs.Component('lipids', formula='C10H20O2', phase='l', particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)

# Define proteins component
proteins = qs.Component('proteins', formula='C6H12O6N2', phase='s',particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)

# Define carbohydrates component
carbs = qs.Component('carbs', formula='C6H12O6', phase='s', particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)

# Define ash component
ash = qs.Component('ash', phase='s', particle_size = 'Particulate', degradability = 'Biological',
                     organic = False)

# Define WasteStreams
feedstock = WasteStream('feedstock')
feedstock.set_flow_by_concentration(flow_tot=100, concentrations={'lipids': 15, 'proteins': 5, 'carbohydrates': 10, 'ash': 5, 'water': 65}, units='kg/hr')

# Define Units
from qsdsan import Mixer, Splitter, StorageTank

mixer = Mixer('mixer')
splitter = Splitter('splitter', split={'lipids': 0.5, 'proteins': 0.3, 'carbohydrates': 0.2})
storage_tank = StorageTank('storage_tank', volume=100, phase='l')

# Connect Units
feedstock.source = mixer
mixer-0-0 == splitter
splitter-0 == storage_tank

# Simulate the system
from qsdsan import System

system = System('example_system', path='.')
system.add(feedstock, mixer, splitter, storage_tank)
system.simulate()

# Display results
print(storage_tank.ins[0].show(flow='kg/hr'))

# Define biocrude oil component
biocrude = qs.Component('biocrude', formula='C10H20O2', phase='l', particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)


# Define NaCl (salt) component
nacl = qs.Component('nacl', formula='NaCl', phase='s', particle_size = 'Particulate', degradability = 'Biological',
                     organic = False)

# Define water component
water = qs.Component('water', formula='H2O', phase='l', particle_size = 'Particulate', degradability = 'Biological',
                     organic = False)

# Define lipids component
lipids = qs.Component('lipids', formula='C10H20O2', phase='l', particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)

# Define proteins component
proteins = qs.Component('proteins', formula='C6H12O6N2', phase='s',particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)

# Define carbohydrates component
carbs = qs.Component('carbs', formula='C6H12O6', phase='s', particle_size = 'Particulate', degradability = 'Biological',
                     organic = True)

# Define ash component
ash = qs.Component('ash', phase='s', particle_size = 'Particulate', degradability = 'Biological',
                     organic = False)


# Create a new BioSTEAM system
bst.main_flowsheet.set_flowsheet('biocrude_processing')

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


                      

# Deasher

# Assuming an ash component in biocrude
ash_content = qs.Chemical('Ash', formula='Ash', phase='s', CAS='123-45-6')
biocrude.imass['Ash'] = 5  # Example: 5 kg/hr of ash in biocrude

# Create a de-ashing unit using Biosteam


class DeashingUnit(bst.Unit):
    def __init__(self, ins, outs):
        super().__init__(ins, outs)
        
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        self.outs[0].remove_chemicals(['Ash'])

# Initialize the de-ashing unit
deashing_unit = DeashingUnit(ins='biocrude', outs='deashed_biocrude')

# Add the de-ashing unit to the flowsheet
bst.main_flowsheet.add_unit(deashing_unit)

# Desalter

# Assume an initial salt content in biocrude
biocrude.imass['NaCl'] = 2  # 2 kg/hr of NaCl in biocrude

# Create a desalting unit using Biosteam

class DesaltingUnit(bst.Unit):
    def __init__(self, ins, outs):
        super().__init__(ins, outs)
        
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        self.outs[0].remove_chemicals(['NaCl'])

# Initialize the desalting unit
desalting_unit = DesaltingUnit(ins='deashed_biocrude', outs='desalted_biocrude')

# Add the desalting unit to the flowsheet
bst.main_flowsheet.add_unit(desalting_unit)

# Dewater
# Assume an initial water content in biocrude
biocrude.imass['Water'] = 10  # Example: 10 kg/hr of water in biocrude

# Create a dewatering unit using Biosteam

class DewateringUnit(bst.Unit):
    def __init__(self, ins, outs):
        super().__init__(ins, outs)
        
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        self.outs[0].remove_chemicals(['Water'])

# Initialize the dewatering unit
dewatering_unit = DewateringUnit(ins='desalted_biocrude', outs='dewatered_biocrude')

# Add the dewatering unit to the flowsheet
bst.main_flowsheet.add_unit(dewatering_unit)

# Create a simple simulation
bst.main_flowsheet.simulate()

# Print results
print(deashing_unit.results())
print(desalting_unit.results())
print(dewatering_unit.results())




