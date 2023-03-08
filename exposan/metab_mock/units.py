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
from qsdsan import Construction
from qsdsan.sanunits import AnaerobicCSTR, Pump, HXutility
from qsdsan.utils import auom
from exposan.metab.equipment import Beads
from exposan.metab.utils import (
    pipe_design, pipe_friction_head, hdpe_price, 
    heat_transfer_U, 
    stainless_steel_wall_thickness as wt_ssteel, 
    add_prefix,
    _construct_water_pump,
    )
        
from math import pi
from collections import defaultdict

__all__ = ('METAB_AnCSTR',)


#%% METAB_AnCSTR
_fts2mhr = auom('ft/s').conversion_factor('m/hr')
_cmph_2_gpm = auom('m3/hr').conversion_factor('gpm')

class METAB_AnCSTR(AnaerobicCSTR):
    
    auxiliary_unit_names = ('heat_exchanger', )
    
    def __init__(self, ID='', lifetime=30, bead_lifetime=10,
                 encapsulate_concentration=25,
                 reactor_height_to_diameter=1.5,
                 mixing_power=10,                   # hp/1000 gallons
                 wall_concrete_unit_cost=1081.73,   # $850/m3 in 2014 USD, converted to 2021 USD with concrete PPI
                 slab_concrete_unit_cost=582.48,    # $458/m3 in 2014 USD 
                 stainless_steel_unit_cost=1.8,     # https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
                 rockwool_unit_cost=0.59,           # https://www.alibaba.com/product-detail/mineral-wool-insulation-price-mineral-wool_60101640303.html?spm=a2700.7724857.0.0.262334d1rZXb48
                 carbon_steel_unit_cost=0.5,        # https://www.alibaba.com/product-detail/ASTM-A106-Ss400-Q235-Standard-Ms_1600406694387.html?s=p
                 **kwargs):
        equip = kwargs.pop('equipment', [])
        equip.append(Beads(ID=f'{ID}_beads', lifetime=bead_lifetime))
        super().__init__(ID, lifetime=lifetime, equipment=equip, **kwargs)
        self.mixing_power=mixing_power
        self.encapsulate_concentration = encapsulate_concentration
        self.reactor_height_to_diameter = reactor_height_to_diameter
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.stainless_steel_unit_cost = stainless_steel_unit_cost
        self.rockwool_unit_cost = rockwool_unit_cost
        self.carbon_steel_unit_cost = carbon_steel_unit_cost
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
        for i in ('Wall concrete', 'Slab concrete', 'Stainless steel', 
                  'Rockwool', 'Carbon steel', 'HDPE pipes'):
            name = i.lower().replace(' ', '_')
            self.construction.append(
                Construction(ID=name, linked_unit=self, item=name)
                )
        for aux in self.auxiliary_units:
            self.construction += aux.construction
        for equip in self.equipment:
            self.construction += equip.construction
    
    @property
    def bead_lifetime(self):
        for equip in self.equipment:
            if isinstance(equip, Beads): return equip.lifetime
            
    @bead_lifetime.setter
    def bead_lifetime(self, lt):
        for equip in self.equipment:
            if isinstance(equip, Beads): 
                equip.update_lifetime(lt)
    
    def _setup(self):
        hasfield = hasattr
        setfield = setattr
        for i, ws in enumerate(self.ins):
            field = f'Pump_ins{i}'
            if not hasfield(self, field):
                pump = Pump(ws.ID+'_Pump', ins=Stream(f'{ws.ID}_proxy'))
                setfield(self, field, pump)
                self.construction += [
                    Construction(ID='22kW', linked_unit=pump, item='pump_22kW'),
                    Construction(ID='40W', linked_unit=pump, item='pump_40W')
                    ]
            self.auxiliary_unit_names = tuple({*self.auxiliary_unit_names, field})
        super()._setup()
        
    _steel_insulate_thickness = 50     # mm
    _concrete_cover_thickness = 100    # mm
    _concrete_wall_thickness = 150     # mm
    _concrete_base_thickness = 160     # mm
    _cncr_insulate_thickness = 25      # mm
    _facing_thickness = 3              # mm
    
    _density = {
        'Aluminum': 2710,        # 2,640 - 2,810 kg/m3
        'Stainless steel': 7930, # kg/m3, 18/8 Chromium
        'Rockwool': 100,         # kg/m3
        'Carbon steel': 7840
        }
    
    T_air = 273.15 + 20
    T_earth = 273.15 + 20
    
    _l_min_velocity = 3*_fts2mhr
    _g_min_velocity = 10*_fts2mhr
    
    _units = {
        'Volume': 'm3',
        'Height': 'm',
        'Outer diameter': 'm',
        'Area': 'm2',
        'Wall concrete': 'm3',
        'Slab concrete': 'm3',
        'Stainless steel': 'kg',
        'Rockwool': 'kg',
        'Carbon steel': 'kg',
        'HDPE pipes': 'kg',
        'Bead volume': 'm3',
        'Agitator': 'hp'
        }
    
    def _design(self):
        D = self.design_results
        den = self._density
        V = D['Volume'] = self.V_liq + self.V_gas
        dia = (V*4/self.reactor_height_to_diameter/pi) ** (1/3)
        h = D['Height'] = dia * self.reactor_height_to_diameter
        tface = self._facing_thickness/1e3
        if V >= 5:
            twall = self._concrete_wall_thickness/1e3
            tcover = self._concrete_cover_thickness/1e3
            tbase = self._concrete_base_thickness/1e3
            tinsl = self._cncr_insulate_thickness/1e3
            D['Outer diameter'] = OD = dia + twall * 2
            S_wall = pi*OD*h
            S_base = D['Area'] = pi*(OD/2)**2
            D['Wall concrete'] = S_wall * twall
            D['Slab concrete'] = S_base*(tcover + tbase)
            D['Stainless steel'] = 0
            D['Rockwool'] = S_wall * tinsl * den['Rockwool']
            D['Carbon steel'] = S_wall * tface * den['Carbon steel']
            Uwall, Ucover, Ubase = heat_transfer_U(twall, tinsl, tface, tbase, concrete=True)
        else:
            P_liq = self.external_P * 1e5 + self._mixed.rho * 9.80665 * h * self.V_liq/V
            twall = tbase = wt_ssteel(P_liq, dia, h)
            tcover = twall/2
            tinsl = self._steel_insulate_thickness/1e3
            D['Outer diameter'] = OD = dia + twall*2
            S_wall = pi*OD*h
            S_base = D['Area'] = pi*(OD/2)**2
            D['Wall concrete'] = 0
            D['Slab concrete'] = 0
            D['Stainless steel'] = (S_wall * twall + S_base*(tcover + tbase)) * den['Stainless steel']
            D['Rockwool'] = (S_base+S_wall) * tinsl * den['Rockwool']
            D['Carbon steel'] = (S_base+S_wall) * tface * den['Carbon steel']
            Uwall, Ucover, Ubase = heat_transfer_U(twall, tinsl, tface, tbase, concrete=False)
        
        # Calculate needed heating
        T = self.T
        hx = self.heat_exchanger
        mixed = self._mixed
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        mixed.mix_from(self.ins)
        hx_ins0.copy_flow(mixed)
        hx_outs0.copy_flow(mixed)
        hx_ins0.T = mixed.T
        hx_outs0.T = T
        hx_ins0.P = hx_outs0.P = mixed.T
        
        # Heat loss
        wall_loss = Uwall * S_wall * (T-self.T_air) # [W]
        base_loss = Ubase * S_base * (T-self.T_earth) # [W]
        cover_loss = Ucover * S_base * (T-self.T_air) # [W]
        duty = (wall_loss+base_loss+cover_loss)*60*60/1e3 # kJ/hr
        hx.H = hx_ins0.H + duty # stream heating and heat loss
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        if T <= mixed.T: hx.heat_utilities[0].cost = 0
        
        # Piping
        L_inlets = OD * 1.25
        L_outlets = h + OD*0.25
        L_gas = h + OD
        pipe_IDs = []
        HDPE_pipes = []
        for ws in self.ins:
            _inch, _kg_per_m = pipe_design(ws.F_vol, self._l_min_velocity)
            pipe_IDs.append(_inch)
            HDPE_pipes.append(_kg_per_m*L_inlets)
        for ws in self.outs:
            if ws.phase == 'g':
                D['Stainless steel'] += pipe_design(ws.F_vol, self._g_min_velocity, True)[1] * L_gas                
            else:
                _inch, _kg_per_m = pipe_design(ws.F_vol, self._l_min_velocity)
                pipe_IDs.append(_inch)
                HDPE_pipes.append(_kg_per_m*L_outlets)
        D['HDPE pipes'] = sum(HDPE_pipes)
        self._hdpe_ids, self._hdpe_kgs = pipe_IDs, HDPE_pipes
        
        # Vertical turbine agitator
        D['Agitator'] = hp = self.mixing_power * auom('m3').convert(self.V_liq, 'gallon')/1000
        self.add_power_utility(auom('hp').convert(hp, 'kW'))
        
        # Pumps
        flowsheet_ID = self.system.flowsheet.ID
        getfield = getattr
        creg = Construction.registry
        for i, ws, in enumerate(self.ins):
            ID = pipe_IDs[i]
            hf = pipe_friction_head(ws.F_vol*_cmph_2_gpm, L_inlets, ID)  # friction head loss
            TDH = hf + h # in m, assume suction head = 0, discharge head = reactor height
            field = f'Pump_ins{i}'
            pump = getfield(self, field)
            pump.ins[0].copy_flow(ws)
            pump.dP_design = TDH * 9804.14  # in Pa
            pump.simulate()
            if self.include_construction:
                q22, q40 = _construct_water_pump(pump)
                p22 = getfield(creg, f'{flowsheet_ID}_{pump.ID}_22kW')
                p40 = getfield(creg, f'{flowsheet_ID}_{pump.ID}_40W')
                p22.quantity = q22
                p40.quantity = q40
        
        if self.include_construction:
            for i in ('Wall concrete', 'Slab concrete', 'Stainless steel', 
                      'Rockwool', 'Carbon steel', 'HDPE pipes'):
                name = i.lower().replace(' ', '_')
                const = getfield(creg, f'{flowsheet_ID}_{self.ID}_{name}')
                const.quantity = D[i]
                        
        # Beads
        cmps = mixed.components
        conc = self._state[:len(cmps)]        
        retained_concentration = sum(conc * cmps.i_mass * (self._f_retain > 0))    # kg mass/m3
        V_beads = D['Bead volume'] = self.V_liq * retained_concentration / self.encapsulate_concentration
        if V_beads >= self.V_liq:
            raise RuntimeError('retained biomass concentration > design encapsualation concentration')
        self.add_equipment_design()

    _NG_price = 0.85*auom('kJ').conversion_factor('therm') # [USD/kJ] 5.47 2021USD/Mcf vs. 4.19 2016USD/Mcf    

    def _cost(self):
        bg = self.outs[0]
        cmps = bg.components
        KJ_per_kg = cmps.i_mass/cmps.chem_MW*cmps.LHV
        self.add_OPEX['NG_offset'] = -sum(bg.mass*KJ_per_kg)*self._NG_price # kJ/hr * USD/kJ = USD/hr
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Wall concrete'] = D['Wall concrete']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        C['Stainless steel'] = D['Stainless steel']*self.stainless_steel_unit_cost
        C['Rockwool'] = D['Rockwool']*self.rockwool_unit_cost
        C['Carbon steel'] = D['Carbon steel']*self.carbon_steel_unit_cost
        C['HDPE pipes'] = sum(hdpe_price(inch)*kg for inch, kg in zip(self._hdpe_ids, self._hdpe_kgs))
        hp = D['Agitator']
        D['Stainless steel'] += 29.792*hp**0.4785   # Empirical function based on data on https://www.agitadoresfluidmix.com/wp-content/uploads/2023/01/4.Data%20sheet%20VTA%20ENG.pdf
        if hp < 2: C['Agitator'] = 2610*hp**0.34 * 708/567 # Seider et al., 2017, pp481; adjusted to 2021$ using CEPCI
        else: C['Agitator'] = 3730*hp**0.54 * 708/567
        self.add_equipment_cost()
    
    def add_equipment_design(self):
        unit_design = self.design_results
        unit_units = self._units
        isa = isinstance
        get = getattr
        if isa(self.equipment_lifetime, int):
            lt = self.equipment_lifetime
            self.equipment_lifetime = defaultdict(lambda: lt)
        F_BM, F_D, F_P, F_M, lifetime = \
            self.F_BM, self.F_D, self.F_P, self.F_M, self.equipment_lifetime

        for equip in self.equipment:
            equip_ID = equip.ID
            prefix = f'{equip.__class__.__name__} {equip_ID}'
            equip_design = equip._design_results = equip._design()
            equip_design = {} if not equip_design else equip_design
            unit_design.update(add_prefix(equip_design, prefix))

            equip_units = {} if not equip.units else equip.units
            unit_units.update(add_prefix(equip_units, prefix))
            for unit_attr, equip_attr in zip(
                    (F_BM, F_D, F_P, F_M, lifetime),
                    ('F_BM', 'F_D', 'F_P', 'F_M', 'lifetime'),
                    ):
                equip_attr = get(equip, equip_attr)
                if isa(equip_attr, dict):
                    unit_attr.update(add_prefix(equip_attr, prefix))
                else:
                    unit_attr[equip_ID] = equip_attr

