# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from exposan.metab import data_path
from qsdsan.utils import ospath
from bw2qsd import DataDownloader as ddld, CFgetter, remove_setups_pickle
from bw2qsd.utils import format_name

# dld = ddld()
# dld.available_databases

def export_ei_CFs(save_to=''):
    ei = CFgetter('ei')
    # ei.load_database('ecoinvent_cutoff38')
    ei.load_database('ecoinvent_apos38')
    ei.load_indicators(method='TRACI', method_exclude='no LT', add=True)
    
    # *********** construction *************
    # reactor (maybe gas holding tank also)
    # consider a double-wall design for better insulation of H2E
    # ei.load_activities('steel production, 18/8', limit=1, add=True)t
    ei.load_activities('market for steel, chromium steel 18/8', add=True)
    # ei.load_activities('polystyrene foam slab, insulation', limit=1, add=True)  # for insulation
    ei.load_activities('market for stone wool', limit=2, add=True)
    # ei.load_activities('cellulose fibre production', add=True)
    ei.load_activities('market for concrete slab', limit=2, add=True)
    ei.load_activities('market for concrete, normal, RNA', add=True, limit=1)
    ei.load_activities('market for concrete, medium strength, RNA', add=True, limit=1)
    ei.load_activities('aluminium production, ingot, RoW', limit=2, add=True)
    ei.load_activities('sheet rolling, aluminium', add=True, limit=1)
    
    # bead materials
    ei.load_activities('chemical factory, organics', limit=2, add=True)
    ei.load_activities('"ethylene glycol production", RoW', add=True)
    ei.load_activities('"methacrylic acid production", RoW', add=True)
    ei.load_activities('polyacrylamide', add=True)
    ei.load_activities('"oxidation of methanol", RoW', add=True) # produces formaldehyde
    ei.load_activities('copper sulfate', limit=2, add=True)
    ei.load_activities('sulfuric acid, RoW', limit=2, add=True)
    ei.load_activities('dimethylamine, RoW', limit=2, add=True)
    ei.load_activities('ethylene dichloride, RoW', limit=2, add=True)
    ei.load_activities('potassium hydroxide', limit=3, add=True)
    
    ei.load_activities('activated carbon, granular', limit=2, add=True)
    ei.load_activities('sodium persulfate', add=True)
    
    # degassing membrane
    # ei.load_activities('Polydimethylsiloxane production', add=True)
    # ei.load_activities('polycarbonate production', add=True)
    # ei.load_activities('polypropylene production', limit=1, add=True)
    # ei.load_activities('polyurethane production', limit=1, add=True)
    
    ei.load_activities('polypropylene, granulate', limit=2, add=True)
    ei.load_activities('pvc, RoW', limit=2, add=True)
    ei.load_activities('epoxy resin, liquid, RoW', limit=2, add=True)
    ei.load_activities('polysulfone', add=True)
    
    # ei.load_activities('fibre reinforced plastic production, polyester', limit=1, add=True)
    # ei.load_activities('"synthetic rubber production"', add=True)
    
    # ei.load_activities('ptfe', add=True)
    # ei.load_activities('polyvinylfluoride production film', add=True)
    
    ei.load_activities('injection moulding, RoW', limit=1, add=True)
    ei.load_activities('extrusion, plastic, RoW', limit=2, add=True)
    
    # biogas holder
    ei.load_activities('textile production polyester', limit=1, add=True)
    ei.load_activities('acrylic varnish, RoW', add=True)
    ei.load_activities('air compressor, 4kW, GLO', add=True)
    ei.load_activities('air compressor, 4kW, RoW', add=True)
    
    # biogas pretreatment -- iron sponge (wood chips + iron oxide)
    ei.load_activities('portafer production, RoW', limit=1, add=True)   # iron oxide
    ei.load_activities('"wood chips production"', add=True)
    ei.load_activities('soda ash, dense', limit=2, add=True)
    ei.load_activities('calcium carbonate, RoW', limit=2, add=True)
    ei.load_activities('reinforcing steel production, RoW', limit=1, add=True)
    
    # CHP
    # ei.load_activities('"mini CHP plant construction"', add=True)
    # ei.load_activities('"mini CHP plant production"', add=True)
    # ei.load_activities('gas boiler production', limit=1, add=True) # optional, only when not using CHP
    
    # pumping & piping
    ei.load_activities('market for pump, 40W', add=True)
    ei.load_activities('market for pump, 22kW', add=True)
    ei.load_activities('polyethylene, high density, granulate', limit=1, add=True)
    ei.load_activities('"polyethylene production, high density, granulate", RoW', limit=1, add=True)
    ei.load_activities('heat and power co-generation, natural gas', add=True)
    ei.load_activities('heat and power co-generation, oil, natural gas, US-WECC', add=True)
    ei.load_activities('plastic processing factory, GLO', limit=1, add=True)
    
    # *********** operation *************
    ei.load_activities('market for electricity, production mix, US-WECC', limit=4, add=True)
    ei.load_activities('heat, central or small-scale, natural gas, RoW, <100kW', limit=7, add=True)
    ei.load_activities('"transport, freight, lorry 3.5-7.5 metric ton", RoW', add=True)
    ei.load_activities('natural gas, low pressure, RoW', limit=1, add=True)  # offset
    
    # membrane cleaning
    ei.load_activities('sodium hypochlorite, RoW', limit=1, add=True)
    ei.load_activities('citric acid GLO', limit=1, add=True)
    
    if save_to: ei.get_CFs(show=False, path=save_to)
    return ei

def format_name_TRACI(name):
    name = name.split(' (')[0]
    name = name.replace(': ', ' ')
    return format_name(name)

def format_alias_TRACI(ind):
    name = ind[2].split(' (')
    if len(name) > 1:
        return name[1].rstrip(')')
    else:
        name = name[0]
        for i in (', ', '-', ': '): name = name.replace(i, ' ')
        return ''.join([i[0].capitalize() for i in name.split(' ')])



#%%
if __name__ == '__main__':
    # path = ospath.join(data_path, 'ecoinvent38apos_CFs.xlsx')
    # ei = export_ei_CFs(path)
    # df = ei.export_indicators(name_formatter=format_name_TRACI, 
    #                           alias_formatter=format_alias_TRACI,
    #                           path=ospath.join(data_path, 'TRACI_indicators.xlsx'))
    ei = export_ei_CFs()
    remove_setups_pickle()