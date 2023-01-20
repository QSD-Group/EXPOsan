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
from qsdsan import SanUnit, Equipment
from qsdsan.sanunits import AnaerobicCSTR, CSTR, Pump, HXutility
from qsdsan.utils import auom, ospath, load_data
from exposan.metab_mock import data_path
from exposan.metab_mock.utils import encap_lci, dm_lci, er_lci
import numpy as np
from warnings import warn
from math import pi, ceil

__all__ = ('rhos_adm1_ph_ctrl',
           'Beads',
           'DegassingMembrane',
           'METAB_AnCSTR',
           'IronSpongeTreatment',
           'DoubleMembraneGasHolder')

#%% rhos_adm1_ph_ctrl
rhos = np.zeros(22) # 22 kinetic processes
Cs = np.empty(19)

from qsdsan.processes._adm1 import (
    R,
    mass2mol_conversion,
    T_correction_factor,
    acid_base_rxn,
    substr_inhibit,
    Hill_inhibit,
    non_compet_inhibit
    )

def rhos_adm1_ph_ctrl(state_arr, params):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    KHb = params['K_H_base']
    Kab = params['Ka_base']
    KH_dH = params['K_H_dH']
    Ka_dH = params['Ka_dH']
    kLa = params['kLa']
    T_base = params['T_base']
    root = params['root']
    

    Cs[:8] = state_arr[12:20]
    Cs[8:12] = state_arr[19:23]
    Cs[12:] = state_arr[16:23]
    substrates = state_arr[:8]
    S_va, S_bu, S_h2, S_IN = state_arr[[3,4,7,10]]
    unit_conversion = mass2mol_conversion(cmps)
    cmps_in_M = state_arr[:27] * unit_conversion
    weak_acids = cmps_in_M[[24, 25, 10, 9, 6, 5, 4, 3]]

    T_op, pH = state_arr[-2:]   #!!! change in state_arr
    biogas_S = state_arr[7:10].copy()
    biogas_p = R * T_op * state_arr[27:30]
    Kas = Kab * T_correction_factor(T_base, T_op, Ka_dH)
    KH = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[7:10]

    rhos[:-3] = ks * Cs
    Monod = substr_inhibit(substrates, Ks)
    rhos[4:12] *= Monod
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)

    h = 10**(-pH)
    delta = acid_base_rxn(h, weak_acids, Kas)
    S_cat = weak_acids[0] - delta

    nh3 = Kas[1] * weak_acids[2] / (Kas[1] + h)
    co2 = weak_acids[3] - Kas[2] * weak_acids[3] / (Kas[2] + h)
    biogas_S[-1] = co2 / unit_conversion[9]
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Iin = substr_inhibit(S_IN, KS_IN)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    rhos[4:12] *= Iph * Iin
    rhos[6:10] *= Ih2
    rhos[10] *= Inh3
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    root.data = {
        'pH':pH, 
        'Iph':Iph, 
        'Ih2':Ih2, 
        'Iin':Iin, 
        'Inh3':Inh3,
        'Monod':Monod,
        'rhos':rhos[4:12].copy(),
        'gas_transfer':rhos[-3:].copy(),
        'S_cat':S_cat
        }
    return rhos

#%% Beads
class Beads(Equipment):
    
    def __init__(self, F_BM=1.1, lifetime=1, **kwargs):
        super().__init__(F_BM=F_BM, lifetime=lifetime, **kwargs)
        
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
        'PEGDMA_1000': 1017.00,            # USD/kg; https://www.polysciences.com/default/polyethylene-glycol-dimethacrylate-pegdma-1000
        'BIS': 137/0.5,                    # USD/kg; https://www.sigmaaldrich.com/US/en/product/sial/146072
        'TEMED': 169.00,                   # USD/L;  https://www.sigmaaldrich.com/US/en/product/mm/808742
        'APS': 237/2.5,                    # USD/kg; https://www.sigmaaldrich.com/US/en/product/sigald/215589    
        'PAC': 393/5                       # USD/kg; https://www.sigmaaldrich.com/US/en/product/sigald/161551
        }
    
    _units = {k:v[1] for k,v in encap_lci.encap_items.items()}
    
    def _design(self):
        V_beads = self.linked_unit.design_results['Bead volume']
        D = self._design_results
        D.update(encap_lci.encap_material_input(V_beads, **self._recipe))
        return D
        
    def _cost(self):
        V_beads = self.linked_unit.design_results['Bead volume']
        C = encap_lci.encap_material_cost(V_beads, **self._recipe, 
                                          unit_prices=self._price.values())
        self._baseline_purchase_costs = C
        return C

#%% DegassingMembrane
class DegassingMembrane(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    auxiliary_unit_names = ('vacuum_pump', 'water_pump',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=True,
                 tau=0.01, vacuum_pressure=6e4, water_pressure=6e5,
                 H2_degas_efficiency=0.85, CH4_degas_efficiency=0.85, 
                 CO2_degas_efficiency=0.05, gas_IDs=('S_h2', 'S_ch4', 'S_IC'),
                 design_liquid_flow=(1,11), # m3/hr, DuPont Ligasep LDM-040
                 unit_price=4126):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.tau = tau
        self.vacuum_pressure = vacuum_pressure
        self.water_pressure = water_pressure
        self.H2_degas_efficiency = H2_degas_efficiency
        self.CH4_degas_efficiency = CH4_degas_efficiency
        self.CO2_degas_efficiency = CO2_degas_efficiency
        self.gas_IDs = gas_IDs
        self._split = np.zeros(len(self.components))
        self._gas_idx = self.components.indices(gas_IDs)
        self.design_liquid_flow = design_liquid_flow
        self.unit_price = unit_price
    
    @property
    def H2_degas_efficiency(self):
        return self._h2_ermv
    
    @H2_degas_efficiency.setter
    def H2_degas_efficiency(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'degassing efficiency must be within [0, 1], not {e}')
        self._h2_ermv = e
        
    @property
    def CH4_degas_efficiency(self):
        return self._ch4_ermv
    
    @CH4_degas_efficiency.setter
    def CH4_degas_efficiency(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'degassing efficiency must be within [0, 1], not {e}')
        self._ch4_ermv = e
    
    @property
    def CO2_degas_efficiency(self):
        return self._co2_ermv
    
    @CO2_degas_efficiency.setter
    def CO2_degas_efficiency(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'degassing efficiency must be within [0, 1], not {e}')
        self._co2_ermv = e
    
    @property
    def split(self):
        s = self._split * 0
        s[self._gas_idx] = [self._h2_ermv, self._ch4_ermv, self._co2_ermv]
        return s   
    
    def _setup(self):
        hasfield = hasattr
        inf, = self.ins
        gas = self.outs[0]
        aux = self.auxiliary_unit_names
        if not hasfield(self, 'vacuum_pump'):
            self.vacuum_pump = Pump('VacPump', ins=Stream(f'{gas.ID}_proxy'),
                                    dP_design=self.vacuum_pressure)
        if not hasfield(self, 'water_pump'):
            self.water_pump = Pump(f'{inf.ID}_Pump', ins=Stream(f'{inf.ID}_proxy'),
                                   P=self.water_pressure)
        self.auxiliary_unit_names = tuple({*aux, 'vacuum_pump', 'water_pump'})
        super()._setup()
    
    def _run(self):
        inf, = self.ins
        gas, liquid = self.outs
        s = self.split
        inf.split_to(gas, liquid, s)
        gas.phase = 'g'

    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.
        
    def _update_state(self):
        arr = self._state
        gas, liquid = self.outs
        s = self.split
        Q_liq = arr[-1]
        if gas.state is None: gas.state = arr*0.0
        gas.state[:-1] = s * arr[:-1] * Q_liq
        gas.state[-1] = 1
        if liquid.state is None: liquid.state = arr*0.0    
        liquid.state[:-1] = (1-s) * arr[:-1]
        liquid.state[-1] = Q_liq

    def _update_dstate(self):
        arr = self._dstate
        gas, liquid = self.outs
        s = self.split
        Q_liq = self._state[-1]
        C_liq = self._state[:-1]
        if gas.dstate is None: gas.dstate = arr*0.0
        gas.dstate[:-1] = s * (arr[:-1] * Q_liq + C_liq * arr[-1])
        gas.dstate[-1] = 0
        if liquid.dstate is None: liquid.dstate = arr*0.0    
        liquid.dstate[:-1] = (1-s) * arr[:-1]
        liquid.dstate[-1] = arr[-1]

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE
    
    def _compile_ODE(self):
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        tau = self.tau
        def dy_dt(t, QC_ins, QC, dQC_ins):
            _dstate[:] = (QC_ins[0] - QC)/tau
            _dstate[-1] = dQC_ins[0,-1]
            _update_dstate()
        self._ODE = dy_dt
    
    _units = {k:v[1] for k,v in dm_lci.DuPont_items.items()}
    
    # https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/MDG-Ligasep-LDM-040-PDS-45-D00501-en.pdf
    _DuPont_specs = dict(
        od_fiber = 210e-6,           # 180-240 um
        dw_fiber = 35e-6,            # 30-40 um
        l_fiber = 536e-3,            # 536 mm
        od_shell = 165e-3,           # 165 mm
        dw_shell = 2.5e-3,           # assume 2.5 mm thick housing plastic
        l_shell = 536e-3,            # 536 mm
        od_potting = 180e-3,         # 180 mm
        l_potting = 36e-3 * 2,       # assume 36 mm of potting length on each end
        od_pipe = 108e-3,            # 108 mm
        h_pipe = (159-165/2)*1e-3,   # mm
        V_liq = 6.5e-3,              # 6.5 L
        total_mass = 10,             # kg
        )

    def _design(self):
        D = self.design_results
        inf, = self.ins
        D['Number'] = self.parallel['self'] = ceil(inf.F_vol/self.design_liquid_flow[1])
        dm_specs = self._DuPont_specs
        D.update(dm_lci.DuPont_input(**dm_specs))
        vac, wat = self.vacuum_pump, self.water_pump        
        vac.ins[0].copy_like(self.outs[0])
        vac.dP_design = self.vacuum_pressure
        vac.simulate()
        wat.ins[0].copy_like(self.ins[0])
        wat.P = self.water_pressure
        wat.simulate()
    
    def _cost(self):
        self.baseline_purchase_costs['Each'] = self.unit_price

#%% METAB_AnCSTR

_fts2mhr = auom('ft/s').conversion_factor('m/hr')
_m2in = auom('m').conversion_factor('inch')
_cmph_2_gpm = auom('m3/hr').conversion_factor('gpm')
_ft2m = auom('ft').conversion_factor('m')

def pipe_design(F_vol, vmin, data):
    ID = (F_vol/vmin/pi)**(1/2) * 2 * _m2in
    df = data[(data.ID <= ID).to_numpy()]
    f_unit = auom(data.Weight.columns[0]).conversion_factor('kg/m')
    if len(df) == 0: 
        ID = data.ID.iloc[0,0] # inch
        kg_per_m = data.Weight.iloc[0,0] * f_unit
    else: 
        ID = df.ID.iloc[-1,0]
        kg_per_m = df.Weight.iloc[-1,0] * f_unit
    return ID, kg_per_m

def pipe_friction_head(q, L, c, ID):
    '''Hazen-Williams equation. 
    https://www.engineeringtoolbox.com/hazen-williams-water-d_797.html'''
    return 2.083e-3 * (100*q / c)**1.852 / ID**4.8655 * L * _ft2m

def hdpe_price(ID):
    '''Price in [USD/kg] as a function of inner diameter [inch],
    projection based on prices in https://hdpesupply.com/hdpe-straight-length-pipe/'''
    return 9.625*ID**(-0.368)

class METAB_AnCSTR(AnaerobicCSTR):
    
    auxiliary_unit_names = ('heat_exchanger', )
    
    def __init__(self, ID='', lifetime=20,
                 encapsulate_concentration=25, 
                 wall_concrete_unit_cost=1081.73,
                 slab_concrete_unit_cost=582.48,
                 stainless_steel_unit_cost=4.19,
                 rockwool_unit_cost=0.59,
                 aluminum_unit_cost=15.56,
                 **kwargs):
        equip = kwargs.pop('equipment', [])
        equip.append(Beads(ID=f'{ID}_beads'))
        super().__init__(ID, lifetime=lifetime, equipment=equip, **kwargs)
        self.encapsulate_concentration = encapsulate_concentration
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.stainless_steel_unit_cost = stainless_steel_unit_cost
        self.rockwool_unit_cost = rockwool_unit_cost
        self.aluminum_unit_cost = aluminum_unit_cost
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
    
    def _setup(self):
        hasfield = hasattr
        setfield = setattr
        for i, ws in enumerate(self.ins):
            field = f'Pump_ins{i}'
            if not hasfield(self, field):
                setfield(self, field, Pump(ws.ID+'_Pump', ins=Stream(f'{ws.ID}_proxy')))
            self.auxiliary_unit_names = tuple({*self.auxiliary_unit_names, field})
        super()._setup()
        
    
    _reactor_height_to_diameter = 1.5
    _steel_separator_thickness = 15    # mm
    _steel_wall_thickness = 30         # mm
    _steel_base_thickness = 40         # mm
    _steel_insulate_thickness = 50     # mm
    _concrete_cover_thickness = 100    # mm
    _concrete_wall_thickness = 300     # mm
    _concrete_base_thickness = 300     # mm
    _cncr_insulate_thickness = 25      # mm
    _alum_facing_thickness = 3         # mm
    _gas_separator_r_frac = 0.75       # cone radius to reactor radius
    _gas_separator_h2r = 1/3**(1/2)    # cone height to cone radius
    _baffle_slope = 2/3**(1/2)         # slant/base of the baffles on the side wall
    
    _density = {
        'Aluminum': 2710,        # 2,640 - 2,810 kg/m3
        'Stainless steel': 7930, # kg/m3, 18/8 Chromium
        'Rockwool': 100,         # kg/m3
        }
    
    _h_air = 37.5       # W/m2/K, convective heat transfer coefficient, assume at free air relative velocity = 10 m/s, https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
    _h_water = 3000     # W/m2/K, assume moderate forced flow, https://www.engineersedge.com/heat_transfer/convective_heat_transfer_coefficients__13378.htm
    _k_rockwool = 0.038 # 0.035 – 0.040 W/m/K, thermal conductivity
    _k_ssteel = 20      # 16-24 W/m/K
    _k_alum = 230       # 205-250
    
    T_air = 273.15 + 20
    T_earth = 273.15 + 20
    
    _l_min_velocity = 3*_fts2mhr
    _g_min_velocity = 10*_fts2mhr
    
    _HDPE_pipe_path = ospath.join(data_path, 'HDPE_pipe_chart.xlsx')               # https://www.piping-designer.com/index.php/datasheets/piping-datasheets/1651-pipe-hdpe-ansi-dr-11-0-ips-in
    _ssteel_pipe_path = ospath.join(data_path, 'stainless_steel_pipe_chart.xlsx')  # https://amerpipe.com/wp-content/uploads/2015/10/APP-chart-v7-web.pdf
    
    _heat_transfer_coefficients = {
        'Concrete wall, insulated': 0.7, # 0.6-0.8 W/m2/K
        'Concrete floor': 1.7,
        'Concrete cover, insulated': 1.4, # 1.2-1.6 W/m2/K
        }
    
    _Hazen_Williams_coefficients = {
        'Stainless steel': 110,
        'HDPE': 140
        }
    
    _units = {
        'Volume': 'm3',
        'Height': 'm',
        'Outer diameter': 'm',
        'Wall concrete': 'm3',
        'Slab concrete': 'm3',
        'Stainless steel': 'kg',
        'Rockwool': 'kg',
        'Aluminum sheet': 'kg',
        'HDPE pipes': 'kg',
        'Bead volume': 'm3'
        }
    
    def _design(self):
        D = self.design_results
        den = self._density
        V = D['Volume'] = self.V_liq + self.V_gas
        dia = (V*4/self._reactor_height_to_diameter/pi) ** (1/3)
        h = D['Height'] = dia * self._reactor_height_to_diameter
        r_cone = dia/2*self._gas_separator_r_frac
        Vg = 1/3*pi*r_cone**3*self._gas_separator_h2r
        if Vg < 1.5*self.V_gas:
            Vg = 1.5*self.V_gas
            h_cone = Vg/(1/3*pi*r_cone**2)
        S_cone = pi*r_cone*(r_cone + (r_cone**2 + h_cone**2)**(1/2))
        S_baffle = (pi*(dia/2)**2 - pi*(dia/2*(self._gas_separator_r_frac-1))**2)\
            *self._baffle_slope*2
        tface = self._alum_facing_thickness/1e3
        tsep = self._steel_separator_thickness/1e3
        if V >= 25:
            twall = self._concrete_wall_thickness/1e3
            tcover = self._concrete_cover_thickness/1e3
            tbase = self._concrete_base_thickness/1e3
            tinsl = self._cncr_insulate_thickness/1e3
            D['Outer diameter'] = OD = dia + twall * 2
            S_wall = pi*OD*h
            S_base = D['Area'] = pi*(OD/2)**2
            D['Wall concrete'] = S_wall * twall
            D['Slab concrete'] = S_base*(tcover + tbase)
            D['Stainless steel'] = (S_cone+S_baffle) * tsep * den['Stainless steel']
            D['Rockwool'] = (S_base+S_wall) * tinsl * den['Rockwool']
            D['Aluminum sheet'] = (S_base+S_wall) * tface * den['Aluminum']
            U = self._heat_transfer_coefficients
            Uwall = U['Concrete wall, insulated']
            Ucover = U['Concrete cover, insulated']
            Ubase = U['Concrete floor']
        else:
            twall = tcover = self._steel_wall_thickness/1e3
            tbase = self._steel_base_thickness/1e3
            tinsl = self._steel_insulate_thickness/1e3
            D['Outer diameter'] = OD = dia + twall*2
            S_wall = pi*OD*h
            S_base = D['Area'] = pi*(OD/2)**2
            V_stainless = S_wall * twall + S_base*(tcover + tbase)
            D['Wall concrete'] = 0
            D['Slab concrete'] = 0
            D['Stainless steel'] = ((S_cone+S_baffle) * tsep + V_stainless) * den['Stainless steel']
            D['Rockwool'] = (S_base*2+S_wall) * tinsl * den['Rockwool']
            D['Aluminum sheet'] = (S_base*2+S_wall) * tface * den['Aluminum']
            Uwall = Ucover = 1/(1/self._h_water + 1/self._h_air \
                                + twall/self._k_ssteel \
                                + tinsl/self._k_rockwool \
                                + tface/self._k_alum)
            Ubase = 1/(1/self._h_water \
                       + tbase/self._k_ssteel \
                       + tinsl/self._k_rockwool \
                       + tface/self._k_alum)
        
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
        
        # Piping
        df_l = load_data(self._HDPE_pipe_path, header=[0,1], index_col=None)
        df_g = load_data(self._ssteel_pipe_path, header=[0,1], index_col=None)
        L_inlets = OD * 1.25
        L_outlets = h + OD*0.25
        L_gas = h + OD
        pipe_IDs = []
        HDPE_pipes = []
        for ws in self.ins:
            _inch, _kg_per_m = pipe_design(ws.F_vol, self._l_min_velocity, df_l)
            pipe_IDs.append(_inch)
            HDPE_pipes.append(_kg_per_m*L_inlets)
        for ws in self.outs:
            if ws.phase == 'g':
                D['Stainless steel'] += pipe_design(ws.F_vol, self._g_min_velocity, df_g)[1] * L_gas                
            else:
                _inch, _kg_per_m = pipe_design(ws.F_vol, self._l_min_velocity, df_l)
                pipe_IDs.append(_inch)
                HDPE_pipes.append(_kg_per_m*L_outlets)
        D['HDPE pipes'] = sum(HDPE_pipes)
        self._hdpe_ids, self._hdpe_kgs = pipe_IDs, HDPE_pipes
        
        # Pumps
        HWc = self._Hazen_Williams_coefficients
        getfield = getattr
        for i, ws, in enumerate(self.ins):
            ID = pipe_IDs[i]
            hf = pipe_friction_head(ws.F_vol*_cmph_2_gpm, L_inlets, HWc['HDPE'], ID)  # friction head loss
            #!!! consider adding velocity head to promote mixing?
            TDH = hf + h # in m, assume suction head = 0, discharge head = reactor height
            pump = getfield(self, f'Pump_ins{i}')
            pump.ins[0].copy_flow(ws)
            pump.dP_design = TDH * 9804.14  # in Pa
            pump.simulate()
        
        # Beads
        cmps = mixed.components
        conc = self._state[:len(cmps)]        
        retained_concentration = sum(conc * cmps.i_mass * (self._f_retain > 0))    # kg mass/m3
        V_beads = D['Bead volume'] = self.V_liq * retained_concentration / self.encapsulate_concentration
        if V_beads >= self.V_liq:
            raise RuntimeError('retained biomass concentration > design encapsualation concentration')
        self.add_equipment_design()
    
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Wall concrete'] = D['Wall concrete']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        C['Stainless steel'] = D['Stainless steel']*self.stainless_steel_unit_cost
        C['Rockwool'] = D['Rockwool']*self.rockwool_unit_cost
        C['Aluminum sheet'] = D['Aluminum sheet']*self.aluminum_unit_cost
        C['HDPE pipes'] = sum(hdpe_price(inch)*kg for inch, kg in zip(self._hdpe_ids, self._hdpe_kgs))
        self.add_equipment_cost()

#%% IronSpongeTreatment

class IronSpongeTreatment(Equipment):

    # https://projects.sare.org/wp-content/uploads/3b-Iron-Sponge-design-Considerations.pdf
    
    _default_lifetime = {
        'Vessel': 10,
        'Compressor': 10,
        'Control system': 10,
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
    
    _vessel_thickness = 5            # mm
    _carbon_steel_density = 7850     # kg/m3
    _iron_sponge_bulk_density = 800  # kg/m3
    _Fe2O3_content = 15              # lb Fe2O3/bushel
    _min_contact_time = 60           # s
    _compressibility = 1.
    
    _prices = {
        'Control system': 10400,      # $ per
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
        mixed = self._mixed
        product_bgs = [ws for ws in self.linked_unit._system.products\
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
        self.lifetime['Iron sponge'] = t_replace / 365
        hp = 144*P*(Qg*1e6/24/60)*1.41/(33000*0.41)*((psia/P)**(0.41/1.41)-1)   # https://www.engineeringtoolbox.com/horsepower-compressed-air-d_1363.html
        kW = D['Compressor'] = auom('hp').convert(hp/self.blower_efficiency, 'kW')
        self.linked_unit.power_utility.consumption += kW
        S = D['Vessel surface area'] = pi*d*(d/2+h) * 1.05
        D['Carbon steel'] = S * self._vessel_thickness/1e3 * self._carbon_steel_density
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
        super().__init__(F_BM=F_BM, **kwargs)
        self._mixed = Stream(phase='g')
        self.max_holding_time = max_holding_time
        self.T = T
        self.P = P
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.membrane_unit_cost = membrane_unit_cost

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
        mixed = self._mixed
        product_bgs = [ws for ws in self.linked_unit._system.products\
                       if ws.phase == 'g']
        mixed.mix_from(product_bgs)
        Q = self.get_F_mol()*1e3 * 8.314 * self.T/self.P  # m3/hr
        V_max = D['Capcity'] = Q * self.max_holding_time
        D.update(er_lci.gas_holder_input(V_max, self.d_pe, self.d_pvc, 
                                         self.d_slab, self.varnish_spreading_rate))
        r = er_lci.cap_radius_from_V(V_max)
        D['Diameter'] = r*2
        D['Membrane surface area'] = er_lci.A_cap(r)
        return D
    
    def _cost(self):
        D, C = self._design_results, self._baseline_purchase_costs
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        A = D['Membrane surface area']
        C['Membrane'] = A*2*self.membrane_unit_cost
        return C
        