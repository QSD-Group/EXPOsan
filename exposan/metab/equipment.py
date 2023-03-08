# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from biosteam import Stream
from qsdsan import Equipment, Construction
from qsdsan.utils import auom
from exposan.metab_mock.utils import encap_lci, er_lci
from math import pi

__all__ = ('Beads',
           'IronSpongeTreatment',
           'DoubleMembraneGasHolder')

#%% Beads
class Beads(Equipment):
    
    _units = {
        'EG': 'kg',
        'MAA': 'kg',
        'PAM': 'kg',
        'FMD': 'kg',
        'CuSO4': 'kg',
        'H2SO4': 'kg',
        'DMA': 'kg',
        'EDCl': 'kg',
        'KOH': 'kg',
        'NaPS': 'kg',
        'GAC': 'kg'
        }
    
    def __init__(self, F_BM=1.1, lifetime=1, **kwargs):
        super().__init__(F_BM=F_BM, lifetime=lifetime, 
                         units=Beads._units, **kwargs)
        const = []
        # breakpoint()
        for k, v in self._units.items():
            const.append(
                Construction(ID=f'{self.ID}_{k}',
                             item=k, lifetime=lifetime)
                )
        for k in self._manufacturing_unit_input.keys():
            const.append(
                Construction(ID=f'{self.ID}_{k}', 
                             item=k, lifetime=lifetime)
                )
        self.construction = const
        
    # encapsulation recipe
    _recipe = dict(
        n_bead=45,
        d_bead=0.4,         # inch
        PEGDMA_1000=4.5e-3, # kg
        BIS=2.25e-4,        # kg
        TEMED=6e-5,         # L
        APS=3e-5,           # kg
        PAC=3e-4            # kg        
        )

    _price = {
        'PEGDMA_1000': 21,            # 20-22 USD/kg; https://www.echemi.com/produce/pr2210112485-polyethyleneglycoldimethacrylate-99-colourless-liquid-c8h4na2o4-pharmacy-grade-senwayer.html
        'BIS': 10,                    # 5-15 USD/kg; https://www.alibaba.com/product-detail/Best-price-N-N-Methylenebisacrylamide-CAS_1600724924581.html?spm=a2700.galleryofferlist.normal_offer.d_title.58d341d0GVNAh1
        'TEMED': 16*0.775,            # 3-30 USD/kg, 0.775 kg/L; https://www.alibaba.com/product-detail/TMEDA-N-N-N-N-tetramethylethylenediamine_1600669191140.html?spm=a2700.galleryofferlist.normal_offer.d_title.78635f0fXo3DvL
        'APS': 2.8,                   # 1.6-4.0 USD/kg; https://www.alibaba.com/product-detail/Persulfate-Ammonium-Molecular-formula-NH4-2S2O8_1600618995452.html?spm=a2700.galleryofferlist.normal_offer.d_title.17b42c417WytQh    
        'PAC': 1.5,                   # 1.38-1.60 USD/kg; https://www.alibaba.com/product-detail/Best-Sale-325mesh-Wood-Based-Powder_1600694829290.html?spm=a2700.galleryofferlist.normal_offer.d_title.62c1efdbTj3aHI&s=p
        }
    
    _bead_density = 1860    # kg/m3
    _manufacturing_unit_input = {
        'chemical_factory': (4e-10, ''),    # unit/kg
        'electricity': (0.02, 'kWh'),       # kWh/kg
        'heat': (1.6e-6, 'MJ'),             # MJ/kg
        'trucking': (15e-3, 'tonne*km')     # km, assume transport distance is always 15 km
        }
    
    def _design(self):
        linked_unit = self.linked_unit
        V_beads = linked_unit.design_results['Bead volume']
        m_beads = V_beads * self._bead_density # kg
        D = self._design_results
        D.update(encap_lci.encap_material_input(V_beads, **self._recipe))
        creg = Construction.registry
        get = getattr
        for k, v in D.items():
            const = get(creg, f'{self.ID}_{k}')
            const.quantity = v
        for k, v in self._manufacturing_unit_input.items():
            qt, qu = v
            const = get(creg, f'{self.ID}_{k}')
            const.quantity = m_beads*qt
        return D
        
    def _cost(self):
        V_beads = self.linked_unit.design_results['Bead volume']
        C = encap_lci.encap_material_cost(V_beads, **self._recipe, 
                                          unit_prices=self._price.values())
        self._baseline_purchase_costs = C
        return C
    
    def update_lifetime(self, lt):
        self.lifetime = int(lt)
        for const in self.construction:
            const.lifetime = lt


#%% IronSpongeTreatment

class IronSpongeTreatment(Equipment):

    # https://projects.sare.org/wp-content/uploads/3b-Iron-Sponge-design-Considerations.pdf
    
    _default_lifetime = {
        'Vessel': int(10),
        'Compressor': int(10),
        'Control system': int(10),
        'Iron sponge': 1,
        }
    
    def __init__(self, influent_H2S_ppmv=5000, reaction_efficiency=0.7, 
                 empty_contact_time=90.0, blower_efficiency=0.6, 
                 design_psia=15, lifetime=None, **kwargs):
        lt = lifetime or self._default_lifetime
        super().__init__(lifetime=lt, **kwargs)
        self._mixed = Stream(phase='g')
        self.influent_H2S_ppmv = influent_H2S_ppmv
        self.reaction_efficiency = reaction_efficiency
        self.empty_contact_time = empty_contact_time
        self.blower_efficiency = blower_efficiency
        self.design_psia = design_psia
        linked_unit = self.linked_unit
        self.construction = [
            Construction(ID=f'{self.ID}_carbon_steel', linked_unit=linked_unit,
                         item='carbon_steel', lifetime=lt['Vessel']),
            Construction(ID=f'{self.ID}_iron_sponge', linked_unit=linked_unit,
                         item='iron_sponge', lifetime=lt['Iron sponge']),
            Construction(ID=f'{self.ID}_compressor', linked_unit=linked_unit,
                         item='air_compressor', lifetime=lt['Compressor']),
            ]
    
    _vessel_thickness = 5            # mm
    _carbon_steel_density = 7850     # kg/m3
    _iron_sponge_bulk_density = 800  # kg/m3
    _Fe2O3_content = 15              # lb Fe2O3/bushel
    _min_contact_time = 60           # s
    _compressibility = 1.
    
    _prices = {
        'Control system': 104,        # $ per
        'Compressor': 5180,           # $/kW
        'Vessel': 3110,               # $/m2, 15 psig
        'Iron sponge': 2.07           # $/kg
        }
    
    _units = {
        'Vessel volume': 'm3',
        'Diameter': 'm',
        'Bed height': 'm',
        'Vessel surface area': 'm2',
        'Carbon steel': 'kg',
        'Iron sponge': 'kg',
        'Media lifetime': 'd',
        'Compressor': 'kW',
        }
    
    def get_F_mol(self):
        'Total molar flow in kmol/hr.'
        mixed = self._mixed
        mf = self.influent_H2S_ppmv/1e6
        cmps = mixed.chemicals
        return sum(mixed.mass * cmps.i_mass / cmps.chem_MW)/(1-mf)
       
    def _design(self):
        D = self._design_results
        linked_unit = self.linked_unit
        lifetime = self.lifetime
        mixed = self._mixed
        product_bgs = [ws for ws in linked_unit._system.products\
                       if ws.phase == 'g']
        mixed.mix_from(product_bgs)
        T, P = mixed.T, auom('Pa').convert(mixed.P, 'psi')
        psia = self.design_psia
        Z = self._compressibility
        kmolph = self.get_F_mol()
        Qg = auom('m3/hr').convert(kmolph * 22.4, 'ft3/d') * 1e-6  # million cubic feet per day
        degree_R = T*1.8
        _V = Qg*degree_R*Z/psia
        ID_min = (360*_V)**(1/2)
        ID_max = 5**(1/2) * ID_min
        ID_min = max(ID_min, (5.34*Qg*self.influent_H2S_ppmv)**(1/2))
        dia = (ID_min + ID_max)/2
        d = D['Diameter'] = auom('inch').convert(dia, 'm')              
        tau = max(self.empty_contact_time, self._min_contact_time)
        H = tau * 60 * _V/dia**2
        h = D['Bed height'] = auom('ft').convert(H, 'm')
        Bu = 4.4e-3 * dia**2 * H   # US bushel
        V = D['Vessel volume'] = auom('bushel').convert(Bu, 'm3')
        D['Iron sponge'] = V * self._iron_sponge_bulk_density
        Fe = self._Fe2O3_content
        t_replace = D['Media lifetime'] = 3.14e-8 * Fe * dia**2 * H * self.reaction_efficiency \
            /(Qg*self.influent_H2S_ppmv*1e-6)
        yr_replace = int(t_replace/365)
        if yr_replace == 0: yr_replace = 1
        D['Iron sponge'] *= yr_replace*365/t_replace
        lifetime['Iron sponge'] = yr_replace
        hp = 144*P*(Qg*1e6/24/60)*1.41/(33000*0.41)*((psia/P)**(0.41/1.41)-1)   # https://www.engineeringtoolbox.com/horsepower-compressed-air-d_1363.html
        kW = D['Compressor'] = auom('hp').convert(hp/self.blower_efficiency, 'kW')
        linked_unit.power_utility.consumption += kW
        S = D['Vessel surface area'] = pi*d*(d/2+h) * 1.05
        D['Carbon steel'] = S * self._vessel_thickness/1e3 * self._carbon_steel_density
        
        const = self.construction
        const[0].quantity = D['Carbon steel']
        const[1].quantity = D['Iron sponge']
        const[1].lifetime = yr_replace
        const[2].quantity = (kW/4)**0.6
        return D
    
    def _cost(self):
        D, C = self._design_results, self._baseline_purchase_costs
        _p = self._prices
        C['Vessel'] = D['Vessel surface area'] * _p['Vessel']
        C['Control system'] = _p['Control system']
        C['Compressor'] = D['Compressor'] * _p['Compressor']
        C['Iron sponge'] = D['Iron sponge'] * _p['Iron sponge']
        return C

#%% DoubleMembraneGasHolder

class DoubleMembraneGasHolder(Equipment):
    
    def __init__(self, max_holding_time=12, T=293.15, P=101325, 
                 slab_concrete_unit_cost=582.48,
                 membrane_unit_cost=1.88, 
                 F_BM=1.2, **kwargs):
        lt = kwargs.pop('lifetime', {})
        super().__init__(F_BM=F_BM, **kwargs)
        self.lifetime = lt
        self._mixed = Stream(phase='g')
        self.max_holding_time = max_holding_time
        self.T = T
        self.P = P
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.membrane_unit_cost = membrane_unit_cost
        const = []
        linked_unit = self.linked_unit
        for i in ('PE', 'PVC', 'Varnish', 'Slab concrete'):
            name = i.lower().replace(' ', '_') if len(i) > 3 else i
            const.append(
                Construction(ID=f'{self.ID}_{name}', linked_unit=linked_unit, item=name)
                )
        self.construction = const

    _density = {
        'PE': 1300,  # 1230-1380 kg/m3
        'PVC': 1380, # kg/m3
        'varnish': 900, # https://www.industrialcoatingsltd.com/app/uploads/2020/10/Selett-Clear-Varnish-Data-Sheet.pdf
        }
    
    d_pe = 1e-3     # assume membrane thickness 1 mm
    d_pvc = 75e-6   # assume surface treatment thickness = 75 um
    d_slab = 0.2    # assume 20 cm
    # https://www.industrialcoatingsltd.com/app/uploads/2020/10/Selett-Clear-Varnish-Data-Sheet.pdf
    varnish_spreading_rate = 10.65e3 # 11.6 - 9.7 m2/L
    
    _units =  {
        'PE': 'kg',
        'PVC': 'kg',
        'Varnish': 'kg',
        'Slab concrete': 'm3',
        'Capacity': 'm3',
        'Diameter': 'm',
        'Membrane surface area': 'm2'
        }
    
    def get_F_mol(self):
        'Total molar flow in kmol/hr.'
        mixed = self._mixed
        cmps = mixed.chemicals
        return sum(mixed.mass * cmps.i_mass / cmps.chem_MW)
    
    def _design(self):
        D = self._design_results
        linked_unit = self.linked_unit
        mixed = self._mixed
        product_bgs = [ws for ws in linked_unit._system.products\
                       if ws.phase == 'g']
        mixed.mix_from(product_bgs)
        Q = self.get_F_mol()*1e3 * 8.314 * self.T/self.P  # m3/hr
        V_max = D['Capcity'] = Q * self.max_holding_time
        D.update(er_lci.gas_holder_input(V_max, self.d_pe, self.d_pvc, 
                                         self.d_slab, self.varnish_spreading_rate))
        r = er_lci.cap_radius_from_V(V_max)
        D['Diameter'] = r*2
        D['Membrane surface area'] = er_lci.A_cap(r)
        const = self.construction
        for i, key in enumerate(('PE', 'PVC', 'Varnish', 'Slab concrete')):
            const[i].quantity = D[key]
        return D
    
    def _cost(self):
        D, C = self._design_results, self._baseline_purchase_costs
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        A = D['Membrane surface area']
        C['Membrane'] = A*2*self.membrane_unit_cost
        return C
        
