# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from exposan.metab_mock import data_path
from qsdsan.utils import ospath
from bw2qsd import DataDownloader as ddld, CFgetter, remove_setups_pickle

# dld = ddld()
# dld.available_databases

def export_ei_CFs(save_to=''):
    ei = CFgetter('ei')
    ei.load_database('ecoinvent_cutoff38')
    ei.load_indicators(method='TRACI', method_exclude='no LT', add=True)
    
    # *********** construction *************
    # reactor (maybe gas holding tank also)
    # consider a double-wall design for better insulation of H2E
    ei.load_activities('steel production, 18/8', limit=1, add=True)
    ei.load_activities('polystyrene foam slab, insulation', limit=1, add=True)  # for insulation
    ei.load_activities('"stone wool production, packed"', limit=1, add=True)
    ei.load_activities('cellulose fibre production', add=True)
    ei.load_activities('concrete', limit=2, add=True)
    ei.load_activities('concrete production North America', limit=3, add=True)
    ei.load_activities('aluminium production, ingot', limit=2, add=True)
    ei.load_activities('sheet rolling, aluminium', add=True, limit=1)
    
    # bead materials
    ei.load_activities('"ethylene glycol production"', add=True)
    ei.load_activities('"methacrylic acid"', filter=dict(location='RoW'), limit=1, add=True)
    ei.load_activities('acrylamide', add=True)
    ei.load_activities('oxidation of methanol', limit=1, add=True) # produces formaldehyde
    ei.load_activities('"copper sulfate production"', add=True)
    ei.load_activities('sulfuric acid production', limit=1, add=True)
    # ei.load_activities('ethylenediamine', limit=1, add=True)
    ei.load_activities('"dimethylamine production"', add=True)
    ei.load_activities('"ethylene dichloride production"', limit=1, add=True)
    ei.load_activities('"potassium hydroxide production"', add=True)
    
    ei.load_activities('"activated carbon production"', limit=1, add=True)
    ei.load_activities('sodium persulfate', limit=1, add=True)
    
    # degassing membrane
    ei.load_activities('Polydimethylsiloxane production', add=True)
    ei.load_activities('polycarbonate production', add=True)
    ei.load_activities('polypropylene production', limit=1, add=True)
    ei.load_activities('polyurethane production', limit=1, add=True)
    
    ei.load_activities('"polyethylene production"', limit=3, add=True) # HDPE, LDPE
    # ei.load_activities('propylene production', limit=1, add=True)
    ei.load_activities('pvc', limit=3, add=True)
    ei.load_activities('epoxy resin', limit=1, add=True)
    ei.load_activities('polysulfone production', add=True)
    
    ei.load_activities('fibre reinforced plastic production, polyester', limit=1, add=True)
    ei.load_activities('"synthetic rubber production"', add=True)
    
    ei.load_activities('ptfe', add=True)
    ei.load_activities('polyvinylfluoride production film', add=True)
    
    ei.load_activities('injection moulding CA-QC', limit=2, add=True)
    ei.load_activities('extrusion plastic pipes CA-QC', add=True)
    
    # biogas holding tank
    ei.load_activities('"storage production, 650 l"', add=True)
    ei.load_activities('textile production polyester', limit=1, add=True)
    ei.load_activities('varnish production', add=True)
    # biogas pipe
    ei.load_activities('steel pipe production', limit=1, add=True)   
    ei.load_activities('compressor production, 4kW', filter=dict(location='RER'), add=True)
    # biogas pretreatment -- iron sponge (wood chips + iron oxide)
    ei.load_activities('ferric oxide', limit=1, add=True)
    ei.load_activities('"wood chips production"', add=True)
    ei.load_activities('soda ash, dense, neutralising agent', add=True)
    ei.load_activities('calcium carbonate', limit=1, add=True)
    ei.load_activities('"heat pump production, for mini CHP plant"', add=True)
    # CHP
    ei.load_activities('"mini CHP plant construction"', add=True)
    ei.load_activities('"mini CHP plant production"', add=True)
    ei.load_activities('gas boiler production', limit=1, add=True) # optional, only when not using CHP
    
    # pumping & piping
    ei.load_activities('"pump production"', limit=2, add=True)
    ei.load_activities('polyethylene pipe production', add=True)
    
    # *********** operation *************
    ei.load_activities('electricity production mix US-MRO', limit=1, add=True)
    ei.load_activities('heat production, natural gas, boiler <100kW', limit=1, add=True)
    
    # ei.get_CFs(show=False, path=save_to)
    return ei

#%%
if __name__ == '__main__':
    path = ospath.join(data_path, 'CFs.xlsx')
    ei = export_ei_CFs(path)
    remove_setups_pickle()