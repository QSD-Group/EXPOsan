# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from biosteam import Stream, VacuumSystem
from qsdsan import SanStream, WasteStream, CompiledProcesses, SanUnit, Construction
from qsdsan.sanunits import AnaerobicCSTR, Pump, HXutility
from qsdsan.utils import auom
from exposan.metab.equipment import Beads
from exposan.metab.utils import (
    dm_lci, 
    pipe_design, pipe_friction_head, hdpe_price, 
    heat_transfer_U, 
    stainless_steel_wall_thickness as wt_ssteel, 
    UASB_sizing,
    add_prefix,
    _construct_water_pump,
    _construct_vacuum_pump,
    _F_mass, _F_vol
    )
import numpy as np, flexsolve as flx
from math import pi, ceil
from warnings import warn
from collections import defaultdict

__all__ = ('DegassingMembrane',
           'UASB',
           'METAB_FluidizedBed',
           'METAB_PackedBed',
           'METAB_BatchExp')

_fts2mhr = auom('ft/s').conversion_factor('m/hr')
_cmph_2_gpm = auom('m3/hr').conversion_factor('gpm')

#%% DegassingMembrane
class DegassingMembrane(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    auxiliary_unit_names = ('vacuum_pump', 'water_pump',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=True,
                 vacuum_pressure=7e4, H2_degas_efficiency=0.625, CH4_degas_efficiency=0.455, 
                 CO2_degas_efficiency=0.13, gas_IDs=('S_h2', 'S_ch4', 'S_IC'),
                 design_liquid_flow=(1,11), # m3/hr, DuPont Ligasep LDM-040
                 unit_price=4126):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.vacuum_pressure = vacuum_pressure
        self.H2_degas_efficiency = H2_degas_efficiency
        self.CH4_degas_efficiency = CH4_degas_efficiency
        self.CO2_degas_efficiency = CO2_degas_efficiency
        self.gas_IDs = gas_IDs
        self._split = np.zeros(len(self.components))
        self._gas_idx = self.components.indices(gas_IDs)
        self.design_liquid_flow = design_liquid_flow
        self.unit_price = unit_price
        self.construction += [
            Construction(ID=i, linked_unit=self, item=i)\
                for i, u in self._material_units.items()
            ]
        # just to account for env impacts & additional OPEX
        self.NaOCl = SanStream(f'{ID}_NaOCl', H2O=1)
        self.citric_acid = SanStream(f'{ID}_citric_acid', H2O=1)
    
    @property
    def tau(self):
        '''HRT in d.'''
        return self._DuPont_specs['V_liq'] / self.design_liquid_flow[1] / 24
    
    @property
    def pressure_drop(self):
        '''Pressure drop in Pa.'''
        specs = self._DuPont_specs
        # Q = self.design_liquid_flow[1]
        # d = specs['od_fiber'] - specs['dw_fiber'] * 2
        # tau = specs['V_liq'] / Q * 3600
        L = specs['l_fiber']
        # v = L/tau
        # mu = self.ins[0].mu
        # return 32*mu*L*v/d**2
        dP = 2*6895 # 2 psi pressure drop, shell side liquid, PermSelect, max liquid flowrate, fiber length 14.2 cm
        return dP * L/14.2e-2

    
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
        # gas = self.outs[0]
        aux = self.auxiliary_unit_names
        if not hasfield(self, 'vacuum_pump'):
        #     pump = self.vacuum_pump = Pump(f'{self.ID}_VacPump', ins=Stream(f'{gas.ID}_proxy'),
        #                                     dP_design=self.vacuum_pressure)
        #     self.construction.append(
        #         Construction(ID='surrogate', linked_unit=pump, item='air_compressor')
        #         )
            self.construction.append(
                Construction(ID='VacPump_surrogate', linked_unit=self, item='air_compressor')
                )
        if not hasfield(self, 'water_pump'):
            pump = self.water_pump = Pump(f'{inf.ID}_Pump', ins=Stream(f'{inf.ID}_proxy'))
            self.construction += [
                Construction(ID='22kW', linked_unit=pump, item='pump_22kW'),
                Construction(ID='40W', linked_unit=pump, item='pump_40W')
                ]
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
    
    _material_units = {
        'PP': 'kg',
        'PVC': 'kg',
        'PS': 'kg',
        'epoxy': 'kg',
        'electricity': 'kWh',
        'molding': 'kg',
        'extrusion': 'kg',
        }
    
    _units = {
        **_material_units,
        'Number': '',
        'Cleaning frequncy': 'month^(-1)'
        }
    
    _NG_price = 0.85*auom('kJ').conversion_factor('therm') # [USD/kJ] 5.47 2021USD/Mcf vs. 4.19 2016USD/Mcf
    
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
    
    _min_cleaning_frequency = 0.5    # times per month, assuming operated at 11 m3/h, TSS = 1 ppm, TOC = 1 ppm
       
    def _calc_NaOCl(self, freq):
        V_sol = self._DuPont_specs['V_liq'] * 2 * freq * 12  # L/yr
        m_NaOCl = V_sol * 1.22 * 0.125   # kg/yr, 1.22 kg/L density of 12.5% solution
        return m_NaOCl/365/24   # kg/hr pure NaOCl
    
    def _calc_citric_acid(self, freq):
        V_sol = self._DuPont_specs['V_liq'] * 2 * freq * 12  # L/yr
        m_ca = V_sol * 0.03 * 193/1000   # kg/yr, 30mM solution, MW ~ 193
        return m_ca/365/24      # kg/hr pure citric acid

    def _design(self):
        D = self.design_results
        inf, = self.ins
        D['Number'] = ceil(inf.F_vol/self.design_liquid_flow[1])
        dm_specs = self._DuPont_specs
        D.update(dm_lci.DuPont_input(**dm_specs))
        # vac, wat = self.vacuum_pump, self.water_pump        
        # vac.ins[0].copy_like(self.outs[0])
        # vac.dP_design = self.vacuum_pressure
        # vac.simulate()
        gas = self.outs[0]
        V_gas = dm_lci.V_lumen(**dm_specs)
        self.vacuum_pump = vac = VacuumSystem(
            self, F_mass=_F_mass(gas), F_vol=_F_vol(gas), 
            P_suction=101325-self.vacuum_pressure, vessel_volume=V_gas
            )
        wat = self.water_pump
        wat.ins[0].copy_like(self.ins[0])
        wat.dP_design = self.pressure_drop
        wat.simulate()
        freq = D['Cleaning frequncy'] = \
            (max(inf.get_TSS(), inf.composite('C', organic=True))\
            * inf.F_vol/11)**0.6 * self._min_cleaning_frequency # mg/L * m3/hr = g/h
        creg = Construction.registry
        get = getattr
        if self.include_construction:
            flowsheet_ID = self.system.flowsheet.ID
            for i in self._material_units.keys():
                const = get(creg, f'{flowsheet_ID}_{self.ID}_{i}')
                const.quantity = D[i]
            q22, q40 = _construct_water_pump(wat)
            p22 = get(creg, f'{flowsheet_ID}_{wat.ID}_22kW')
            p40 = get(creg, f'{flowsheet_ID}_{wat.ID}_40W')
            p22.quantity = q22
            p40.quantity = q40
            qvac = _construct_vacuum_pump(vac)
            # pvac = get(creg, f'{flowsheet_ID}_{vac.ID}_surrogate')
            pvac = get(creg, f'{flowsheet_ID}_{self.ID}_VacPump_surrogate')
            pvac.quantity = qvac
        self.NaOCl.F_mass = self._calc_NaOCl(freq)
        self.citric_acid.F_mass = self._calc_citric_acid(freq)
        
    def _cost(self):
        bg = self.outs[0]
        cmps = bg.components
        KJ_per_kg = cmps.i_mass/cmps.chem_MW*cmps.LHV
        bg.price = sum(bg.mass*KJ_per_kg)/bg.F_mass*self._NG_price # kJ/kg * USD/kJ = USD/kg
        D, C = self.design_results, self.baseline_purchase_costs
        C['Module'] = self.unit_price * D['Number']
        self.add_OPEX['NaOCl'] = self.NaOCl.F_mass/0.125 * 0.78 # USD/hr, $0.78/kg 12.5% solution, https://www.alibaba.com/product-detail/wholesale-sodium-hypochlorite-NaClO-15-Industrial_1600307294563.html?spm=a2700.galleryofferlist.normal_offer.d_title.21145d84U7uilV
        self.add_OPEX['citric_acid'] = self.citric_acid.F_mass * 0.75 # USD/hr, $0.75/kg https://www.alibaba.com/product-detail/Best-selling-powder-lemon-acid-price_1600657054188.html?spm=a2700.galleryofferlist.0.0.7141505dSfKzKN

#%% UASB
class UASB(AnaerobicCSTR):
    
    auxiliary_unit_names = ('heat_exchanger', )
    def __init__(self, ID='', lifetime=30, T=295.15,
                 fraction_retain=0.963, pH_ctrl=False,
                 max_depth_to_diameter=4,
                 design_upflow_velocity=0.5,        # m/h
                 wall_concrete_unit_cost=1081.73,   # $850/m3 in 2014 USD, converted to 2021 USD with concrete PPI
                 slab_concrete_unit_cost=582.48,    # $458/m3 in 2014 USD 
                 stainless_steel_unit_cost=1.8,     # https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
                 rockwool_unit_cost=0.59,           # https://www.alibaba.com/product-detail/mineral-wool-insulation-price-mineral-wool_60101640303.html?spm=a2700.7724857.0.0.262334d1rZXb48
                 carbon_steel_unit_cost=0.5,        # https://www.alibaba.com/product-detail/ASTM-A106-Ss400-Q235-Standard-Ms_1600406694387.html?s=p
                 **kwargs):
        
        super().__init__(ID, lifetime=lifetime, T=T, **kwargs)
        self._f_retain = self.thermo.chemicals.x * fraction_retain
        self.pH_ctrl = pH_ctrl
        self.max_depth_to_diameter = max_depth_to_diameter
        self.design_upflow_velocity = design_upflow_velocity
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
        self._cached_state = None
    
    def _setup(self):
        hasfield = hasattr
        setfield = setattr
        if self.fixed_headspace_P and not hasfield(self, 'vacuum_pump'):
            # gas = self.outs[0]
            # dP = max(0, (self._P_atm-self._P_gas)*1e5)
            # pump = self.vacuum_pump = Pump(gas.ID+'_VacPump', ins=Stream(f'{gas.ID}_proxy'),
            #                                dP_design=dP)
            # self.construction.append(
            #     Construction(ID='surrogate', linked_unit=pump, item='air_compressor')
            #     )
            self.construction.append(
                Construction(ID='VacPump_surrogate', linked_unit=self, item='air_compressor')
                )
            self.auxiliary_unit_names = tuple({*self.auxiliary_unit_names, 'vacuum_pump'})
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
        
    _steel_separator_thickness = 5     # mm
    _steel_insulate_thickness = 50     # mm
    _concrete_cover_thickness = 100    # mm
    _concrete_wall_thickness = 150     # mm
    _concrete_base_thickness = 160     # mm
    _cncr_insulate_thickness = 25      # mm
    _facing_thickness = 3              # mm
    _gas_separator_r_frac = 0.75       # cone radius to reactor radius
    _gas_separator_h2r = 1/3**(1/2)    # cone height to cone radius
    _baffle_slope = 2/3**(1/2)         # slant/base of the baffles on the side wall
    
    _density = {
        'Aluminum': 2710,        # 2,640 - 2,810 kg/m3
        'Stainless steel': 7930, # kg/m3, 18/8 Chromium
        'Rockwool': 100,         # kg/m3
        'Carbon steel': 7840
        }
    
    T_air = 273.15 + 22
    T_earth = 273.15 + 22
    
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
        }
    
    def _run(self):
        super()._run()
        if self._mixed.T > self.T:
            if self.T < self.T_air: self.T_air = self.T
            self._correct_T()
    
    def _correct_T(self):
        U = 1.2e-3      # kW/m2
        c = 4.186       # kJ/kg/C
        m = self._mixed.F_mass/3600 # kg/s
        V, h, dia = UASB_sizing(self._mixed.F_vol*24, self.V_liq, self.V_gas,
                                self.max_depth_to_diameter, 
                                self.design_upflow_velocity)
        S = pi*dia*h + pi*dia**2/2  # m2
        T_in = self._mixed.T
        T_ext = self.T_air
        self.T = (m*c*T_in + U*S*T_ext)/(m*c+U*S)
    
    def _cache_state(self):
        self._cached_state = self._state[:-1].copy()
    
    def _init_state(self):
        if self._cached_state is not None:
            Q = self._mixed.F_vol * 24
            self._state = np.append(self._cached_state, Q)
            self._dstate = self._state * 0.
        else:
            super()._init_state()
        
    def _compile_ODE(self):
        cmps = self.components
        f_rtn = self._f_retain
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        T = self.T
        _params = self.model.rate_function._params
        _f_rhos = lambda state_arr: self.model.flex_rate_function(
            state_arr, _params, T_op=T, pH=self.pH_ctrl, gas_transfer=True
            )
        M_stoichio = self.model.stoichio_eval()
        n_cmps = len(cmps)
        n_gas = self._n_gas
        V_liq = self.V_liq
        V_gas = self.V_gas
        gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]
        if self._fixed_P_gas:
            f_qgas = self.f_q_gas_fixed_P_headspace
        else:
            f_qgas = self.f_q_gas_var_P_headspace
        def dy_dt(t, QC_ins, QC, dQC_ins):
            #!!! to avoid accumulation of floating error due to limited precision
            QC[np.abs(QC) < 2.22044604925e-16] = 0
            S_liq = QC[:n_cmps]
            S_gas = QC[n_cmps: (n_cmps+n_gas)]
            Q_ins = QC_ins[:, -1]
            S_ins = QC_ins[:, :-1] * 1e-3  # mg/L to kg/m3
            Q = sum(Q_ins)
            rhos = _f_rhos(QC)
            #!!! to avoid accumulation of floating error due to limited precision
            rhos[np.abs(rhos) < 2.22044604925e-16] = 0
            _dstate[:n_cmps] = (Q_ins @ S_ins - Q*S_liq*(1-f_rtn))/V_liq \
                + np.dot(M_stoichio.T, rhos)
            q_gas = f_qgas(rhos[-3:], S_gas, T)
            _dstate[n_cmps: (n_cmps+n_gas)] = - q_gas*S_gas/V_gas \
                + rhos[-3:] * V_liq/V_gas * gas_mass2mol_conversion
            _dstate[-1] = dQC_ins[0,-1]
            _update_dstate()
        self._ODE = dy_dt
    
    def biomass_tss(self, biomass_IDs):
        y = self._state
        cmps = self.components
        bm_idx = cmps.indices(biomass_IDs)
        return sum(y[bm_idx] * cmps.i_mass[bm_idx])
    
    def get_retained_mass(self, biomass_IDs):
        return self.biomass_tss(biomass_IDs) * self.V_liq
    
    def _design(self):
        D = self.design_results
        den = self._density
        Q = self._mixed.F_vol * 24
        V, h, dia = UASB_sizing(Q, self.V_liq, self.V_gas,
                                self.max_depth_to_diameter, 
                                self.design_upflow_velocity)
        D['Volume'] = V
        D['Height'] = h
        r_cone = dia/2*self._gas_separator_r_frac
        Vg = 1/3*pi*r_cone**3*self._gas_separator_h2r
        if Vg < 1.5*self.V_gas:
            Vg = 1.5*self.V_gas
        h_cone = Vg/(1/3*pi*r_cone**2)
        S_cone = pi*r_cone*(r_cone + (r_cone**2 + h_cone**2)**(1/2))
        S_baffle = (pi*(dia/2)**2 - pi*(dia/2*(self._gas_separator_r_frac-1))**2)\
            *self._baffle_slope*2
        tface = self._facing_thickness/1e3
        tsep = self._steel_separator_thickness/1e3
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
            D['Stainless steel'] = (S_cone+S_baffle) * tsep * den['Stainless steel']
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
            V_stainless = S_wall * twall + S_base*(tcover + tbase)
            D['Wall concrete'] = 0
            D['Slab concrete'] = 0
            D['Stainless steel'] = ((S_cone+S_baffle) * tsep + V_stainless) * den['Stainless steel']
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
        
        # Pumps
        flowsheet_ID = self.system.flowsheet.ID
        getfield = getattr
        creg = Construction.registry
        if self.fixed_headspace_P:
            # vac = self.vacuum_pump
            # vac.ins[0].copy_like(self.outs[0])
            # vac.dP_design = (self.external_P - self.headspace_P) * 1e5
            # vac.simulate()
            gas = self.outs[0]
            self.vacuum_pump = vac = VacuumSystem(
                self, F_mass=_F_mass(gas), F_vol=_F_vol(gas), 
                P_suction=self.headspace_P * 1e5,
                vessel_volume=self.V_gas
                )
            qvac = _construct_vacuum_pump(vac)
            # pvac = getfield(creg, f'{flowsheet_ID}_{vac.ID}_surrogate')
            pvac = getfield(creg, f'{flowsheet_ID}_{self.ID}_VacPump_surrogate')
            pvac.quantity = qvac
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
        
        self.add_equipment_design()  
        
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
    
    _NG_price = 0.85*auom('kJ').conversion_factor('therm') # [USD/kJ] 5.47 2021USD/Mcf vs. 4.19 2016USD/Mcf

    def _cost(self):
        bg = self.outs[0]
        if bg.F_mass > 0:
            cmps = bg.components
            KJ_per_kg = cmps.i_mass/cmps.chem_MW*cmps.LHV
            # self.add_OPEX['NG_offset'] = -sum(bg.mass*KJ_per_kg)*self._NG_price # kJ/hr * USD/kJ = USD/hr
            bg.price = sum(bg.mass*KJ_per_kg)/bg.F_mass*self._NG_price # kJ/kg * USD/kJ = USD/kg
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Wall concrete'] = D['Wall concrete']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        C['Stainless steel'] = D['Stainless steel']*self.stainless_steel_unit_cost
        C['Rockwool'] = D['Rockwool']*self.rockwool_unit_cost
        C['Carbon steel'] = D['Carbon steel']*self.carbon_steel_unit_cost
        C['HDPE pipes'] = sum(hdpe_price(inch)*kg for inch, kg in zip(self._hdpe_ids, self._hdpe_kgs))
        self.add_equipment_cost()

#%% METAB_FluidizedBed
def Ergun_equation(u, rho_p, rho, d_p, mu, e=0.4):
    return 1.75*rho*u**2/d_p/e**3 + 150*mu*(1-e)*u/d_p**2/e**3 - 9.81*(rho_p-rho)

class METAB_FluidizedBed(AnaerobicCSTR):
    
    auxiliary_unit_names = ('heat_exchanger', )
    
    def __init__(self, ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_liq=3400, V_gas=300, 
                 voidage=0.6, bead_diameter=2, n_layer=5,
                 boundary_layer_thickness=0.01, diffusivity=None,
                 f_diff=0.55, max_encapsulation_tss=16, model=None,
                 pH_ctrl=False, T=295.15, headspace_P=0.1, external_P=1.013, 
                 pipe_resistance=5.0e4, fixed_headspace_P=False,
                 isdynamic=True, exogenous_vars=(), lifetime=30, bead_lifetime=10,
                 reactor_height_to_diameter=1.5, recirculation_ratio=None,
                 wall_concrete_unit_cost=1081.73,   # $850/m3 in 2014 USD, converted to 2021 USD with concrete PPI
                 slab_concrete_unit_cost=582.48,    # $458/m3 in 2014 USD 
                 stainless_steel_unit_cost=1.8,     # https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
                 rockwool_unit_cost=0.59,           # https://www.alibaba.com/product-detail/mineral-wool-insulation-price-mineral-wool_60101640303.html?spm=a2700.7724857.0.0.262334d1rZXb48
                 carbon_steel_unit_cost=0.5,        # https://www.alibaba.com/product-detail/ASTM-A106-Ss400-Q235-Standard-Ms_1600406694387.html?s=p
                 **kwargs):   
        equip = kwargs.pop('equipment', [])
        equip.append(Beads(ID=f'{ID}_beads', lifetime=bead_lifetime))
        SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, isdynamic=isdynamic, 
                         exogenous_vars=exogenous_vars, lifetime=lifetime, 
                         equipment=equip, **kwargs)
        self.split = split
        self.V_gas = V_gas
        self.V_liq = V_liq
        self.voidage = voidage
        self.bead_diameter = bead_diameter
        self.n_layer = n_layer
        self.boundary_layer_thickness = boundary_layer_thickness
        self.diffusivity = diffusivity
        self.f_diff = f_diff
        self.max_encapsulation_tss = max_encapsulation_tss
        
        self.pH_ctrl = pH_ctrl
        self.T = T
        self._q_gas = 0
        self._n_gas = None
        self._gas_cmp_idx = None
        self._state_keys = None
        self._S_vapor = None
        self._biogas = WasteStream(phase='g')
        self.headspace_P = headspace_P
        self.external_P = external_P
        self.pipe_resistance = pipe_resistance
        self.fixed_headspace_P = fixed_headspace_P
        self._mixed = WasteStream()
        self._tempstate = []
        
        self.reactor_height_to_diameter = reactor_height_to_diameter
        self.recirculation_ratio = recirculation_ratio
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.stainless_steel_unit_cost = stainless_steel_unit_cost
        self.rockwool_unit_cost = rockwool_unit_cost
        self.carbon_steel_unit_cost = carbon_steel_unit_cost
        self.model = model
        self._cached_state = None
        
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
    
    def min_fluidizing_velocity(self, bead_density):
        '''
        Estimate minimum fluidizing velocity based on bead size, bead density, 
        voidage, and fluid properties.

        Parameters
        ----------
        bead_density : float
            Overall density of encapsulation beads, in kg/m3.

        Returns
        -------
        [float] Minimal upflow velocity to fluidize beads, in m/h.
        '''
        mixed = self._mixed
        mixed.mix_from(self.ins)
        rho = mixed.rho
        if bead_density <= rho: return 0.
        mu = mixed.mu
        d = self.bead_diameter * 1e-3
        u_min = flx.IQ_interpolation(f=Ergun_equation, x0=0, x1=1e3, 
                                     args=(bead_density, rho, d, mu))
        return u_min*3600
    
    @property
    def recirculation_ratio(self):
        '''
        [float] Recirculation ratio, as in recirculation flowrate over influent flowrate. 
        
        .. note::
            Only used for design and costing. If recirculation streams have been 
            considered in process simulation, should set recirculation_ratio to 0 
            or None to avoid double counting.
        
        '''
        if not self._rQ: 
            u_min = self.min_fluidizing_velocity(Beads._bead_density)
            if u_min <= 0: return 0
            A_bed = (pi*self.V_bed**2/4/self.reactor_height_to_diameter**2)**(1/3)
            A_liq = A_bed * 0.4
            return max(0, u_min * A_liq/self._mixed.F_vol - 1)
        return self._rQ
    @recirculation_ratio.setter
    def recirculation_ratio(self, r):
        if r:
            u_min = self.min_fluidizing_velocity(Beads._bead_density)
            A_bed = (pi*self.V_bed**2/4/self.reactor_height_to_diameter**2)**(1/3)
            A_liq = A_bed * 0.4
            u = self._mixed.F_vol * (1+r)/A_liq
            if u_min and u < u_min:
                warn(f'Recirculation ratio {r} is too low to fluidize beads, '
                     f'current estimated flow velocity through bed is {u:.2f} m/h, '
                     f'minimal fluidizing velocity is {u_min:.2f} m/h. ')
        self._rQ = r
        
    @property
    def V_liq(self):
        '''[float] The liquid-phase volume, in m^3.'''
        return self._V_liq
    @V_liq.setter
    def V_liq(self, V):
        self._V_liq = V

    @property
    def V_gas(self):
        '''[float] The total gas phase volume, in m^3.'''
        return self._V_gas
    @V_gas.setter
    def V_gas(self, V):
        self._V_gas = V        

    @property
    def voidage(self):
        '''[float] Void fraction of the fluidized bed, unitless.'''
        return self.f_void
    @voidage.setter
    def voidage(self, f):
        if f <= 0 or f >= 1:
            raise ValueError(f'voidage must be in (0,1), not {f}')
        self.f_void = f
    
    @property
    def V_beads(self):
        '''[float] Total encapsulant bead volume, in m^3'''
        f_void = self.voidage
        return self.V_liq / f_void * (1-f_void)
    @property
    def V_bed(self):
        '''[float] Total reactor volume, in m^3'''
        return self.V_liq / self.voidage + self.V_gas
    
    @property
    def bead_diameter(self):
        '''[float] Encapsulation bead diameter in mm.'''
        return self.r_beads * 2e3
    @bead_diameter.setter
    def bead_diameter(self, db):
        self.r_beads = db/2 * 1e-3
    
    @property
    def n_layer(self):
        '''[int] Number of layers for discretization of encapsulation beads along radius.'''
        return self.n_dz
    @n_layer.setter
    def n_layer(self, n):
        self.n_dz = int(n)
    
    @property
    def boundary_layer_thickness(self):
        '''[float] Thickness of liquid boundary layer around beads, in mm.'''
        return self.l_bl * 1e3
    @boundary_layer_thickness.setter
    def boundary_layer_thickness(self, l):
        self.l_bl = min(l, self.r_beads/10)
    
    _diffusivities = {
        'S_su': 4.56e-06,
        'S_aa': 8.62e-06,
        'S_fa': 5.33e-06,
        'S_va': 5e-06,
        'S_bu': 5.04e-06,
        'S_pro': 6e-06,
        'S_ac': 6.48e-06,
        'S_h2': 0.000402,
        'S_ch4': 0.000136,
        'S_IC': 0.000171,
        'S_IN': 0.000152,
        'S_I': 6e-06,
        'S_cat': 0.000117,
        'S_an': 0.000117
        }
    
    
    @property
    def diffusivity(self):
        '''[numpy.ndarray] Diffusivities of components in pure water, in m^2/d.'''
        return self._diff
    @diffusivity.setter
    def diffusivity(self, arr):
        isa = isinstance
        if arr is None:
            self._set_diffusivities(**self._diffusivities)
        elif isa(arr, (list, tuple, np.ndarray)):
            arr = np.asarray(arr)
            if arr.shape == (len(self.thermo.chemicals)):
                self._diff = arr
            else: raise ValueError(f'diffusivity should be an array of the same length'
                                   f'as the components, not of shape {arr.shape}')
        elif isa(arr, dict):
            self._set_diffusivities(**arr)
        else:
            raise TypeError(f'diffusivity must be array-like or a dict, if not None, not {type(arr)}')
        
    def _set_diffusivities(self, **kwargs):
        cmps = self.thermo.chemicals            
        idx = cmps.indices(kwargs.keys())
        if not hasattr(self, '_diff'):
            self._diff = np.zeros(len(cmps))
        self._diff[idx] = list(kwargs.values())
        
    @property
    def f_diff(self):
        '''[float] Encapsulant/water diffusivity ratio, unitless.'''
        return self._f_diff
    @f_diff.setter
    def f_diff(self, f):
        self._f_diff = f
    
    @property
    def k_bl(self):
        '''[float] Mass transfer coefficients through the liquid boundary layer, in m/d.'''
        return self._diff / self.l_bl
    
    @property
    def D(self):
        '''[float] Diffusivities through encapsulant, in m^2/d'''
        return self._diff * self.f_diff
        
    @property
    def max_encapsulation_tss(self):
        '''[float] Maximum biomass encapsulation density, in gTSS/L'''
        return self.K_tss * 2
    @max_encapsulation_tss.setter
    def max_encapsulation_tss(self, tss):
        self.K_tss = tss / 2
    
    @property
    def detachment_half_saturation_tss(self):
        '''[float] Half saturation coefficient for particulate detachment from 
        the encapsulation matrix, in gTSS/L'''
        return self.K_tss
    
    @property
    def model(self):
        '''[:class:`CompiledProcesses`] Anaerobic digestion model.'''
        return self._model
    @model.setter
    def model(self, model):
        if not isinstance(model, CompiledProcesses): 
            raise TypeError(f'model must be a CompiledProesses, not {type(model)}')
        self._model = model
        self._prep_model()

    def _prep_model(self):
        model = self.model
        self._S_vapor = self.ideal_gas_law(p=self.p_vapor())
        self._n_gas = len(model._biogas_IDs)
        cmps = self.thermo.chemicals
        layers = [*range(self.n_dz), 'bulk']
        self._state_keys = [f'{cmp}-{i}' for i in layers for cmp in cmps.IDs] \
            + [ID+'_gas' for ID in self.model._biogas_IDs] \
            + ['Q']
        self._gas_cmp_idx = cmps.indices(self.model._biogas_IDs)
        self._state_header = self._state_keys
        self._gas_state_idx = dict(zip(self.model._biogas_IDs, range(self._n_gas)))  
    
    def _setup(self):
        hasfield = hasattr
        setfield = setattr
        if self.fixed_headspace_P and not hasfield(self, 'vacuum_pump'):
            self.construction.append(
                Construction(ID='VacPump_surrogate', linked_unit=self, item='air_compressor')
                )
            self.auxiliary_unit_names = tuple({*self.auxiliary_unit_names, 'vacuum_pump'})
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
            
    def _run(self):
        super()._run()
        if self._mixed.T > self.T:
            if self.T < self.T_air: self.T_air = self.T
            self._correct_T()
    
    def _correct_T(self):
        U = 1.2e-3      # kW/m2
        c = 4.186       # kJ/kg/C
        m = self._mixed.F_mass/3600 # kg/s
        h2d = self.reactor_height_to_diameter
        dia = (4*self.V_bed/pi/h2d)**(1/3)
        h = dia * h2d
        S = pi*dia*h + pi*dia**2/2  # m2
        T_in = self._mixed.T
        T_ext = self.T_air
        self.T = (m*c*T_in + U*S*T_ext)/(m*c+U*S)
    
    def set_init_conc(self, arr=None, **kwargs):
        '''set the initial concentrations [kg/m3] of components in the fluidized bed, 
        applies uniformly to all layers in beads and bulk liquid unless an array is input.'''
        cmps = self.thermo.chemicals
        cmpx = cmps.index
        n_dz = self.n_dz
        if arr is None:
            Cs = np.zeros(len(cmps))
            for k, v in kwargs.items(): Cs[cmpx(k)] = v
            self._concs = np.tile(Cs, n_dz+1)
        else:
            arr = np.asarray(arr)
            if arr.shape != (len(cmps)*(n_dz+1),):
                raise ValueError(f'arr must be None or a 1d array of length {len(cmps)*(n_dz+1)}, '
                                 f'not {arr.shape}')
            self._concs = arr
    
    def _cache_state(self):
        self._cached_state = self._state[:-(self._n_gas+1)].copy()
    
    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        if self._cached_state is not None: Cs = self._cached_state
        elif self._concs is not None: Cs = self._concs
        else: Cs = np.tile(mixed.conc, (self.n_dz+1))
        self._state = np.append(Cs, [0]*self._n_gas + [Q]).astype('float64')
        self._dstate = self._state * 0.

    def _update_state(self):
        cmps = self.components
        n_cmps = len(cmps)
        n_gas = self._n_gas
        y = self._state[-(n_cmps+n_gas+1):]
        i_mass = cmps.i_mass
        chem_MW = self.components.chem_MW
        Cs = y[:n_cmps]*1e3 # kg/m3 to mg/L
        if self.split is None:
            gas, liquid = self._outs
            if liquid.state is None:
                liquid.state = np.append(Cs, y[-1])
            else:
                liquid.state[:n_cmps] = Cs
                liquid.state[-1] = y[-1]
        else:
            gas = self._outs[0]
            liquids = self._outs[1:]
            for liquid, spl in zip(liquids, self.split):
                if liquid.state is None:
                    liquid.state = np.append(Cs, y[-1]*spl)
                else:
                    liquid.state[:n_cmps] = Cs
                    liquid.state[-1] = y[-1]*spl        
        if gas.state is None:
            gas.state = np.zeros(n_cmps+1)
        gas.state[self._gas_cmp_idx] = y[n_cmps:(n_cmps + n_gas)]
        gas.state[cmps.index('H2O')] = self._S_vapor
        gas.state[-1] = self._q_gas
        gas.state[:n_cmps] = gas.state[:n_cmps] * chem_MW / i_mass * 1e3 # i.e., M biogas to mg (measured_unit) / L

    def _update_dstate(self):
        self._tempstate = self.model.rate_function._params['root'].data.copy()
        n_cmps = len(self.components)
        n_gas = self._n_gas
        dy = self._dstate[-(n_cmps+n_gas+1):]
        dCs = dy[:n_cmps]*1e3
        if self.split is None:
            gas, liquid = self._outs
            if liquid.dstate is None:
                liquid.dstate = np.append(dCs, dy[-1])
            else:
                liquid.dstate[:n_cmps] = dCs
                liquid.dstate[-1] = dy[-1]
        else:
            gas = self._outs[0]
            liquids = self._outs[1:]
            for liquid, spl in zip(liquids, self.split):
                if liquid.dstate is None:
                    liquid.dstate = np.append(dCs, dy[-1]*spl)
                else:
                    liquid.dstate[:n_cmps] = dCs
                    liquid.dstate[-1] = dy[-1]*spl
        if gas.dstate is None:
            # contains no info on dstate
            gas.dstate = np.zeros(n_cmps+1)
    
    @property
    def ODE(self):
        if self._ODE is None:
            self._prep_model()
            self._compile_ODE()
        return self._ODE
    
    def _compile_ODE(self):
        cmps = self.components
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        n_cmps = len(cmps)
        n_dz = self.n_layer
        n_gas = self._n_gas
        T = self.T
        
        _f_rhos = self.model.flex_rate_function
        _params = self.model.rate_function.params
        stoi_bk = self.model.stoichio_eval()
        stoi_en = stoi_bk[:-n_gas]  # no liquid-gas transfer
        Rho_bk = lambda state_arr: _f_rhos(state_arr, _params, 
                                           T_op=T, pH=self.pH_ctrl, 
                                           gas_transfer=True)
        Rho_en = lambda state_arr: _f_rhos(state_arr, _params, 
                                           T_op=T, pH=self.pH_ctrl, 
                                           gas_transfer=False)
        gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]

        V_liq = self.V_liq
        V_gas = self.V_gas
        V_beads = self.V_beads
        r_beads = self.r_beads
        A_beads = 3 * V_beads / r_beads # m2, total bead surface area

        dz = r_beads / n_dz
        zs = np.linspace(dz, r_beads, n_dz)
        dV = 4/3*np.pi*(zs)**3
        dV[1:] -= dV[:-1]
        V_bead = (4/3*np.pi*r_beads**3)
        
        D = self.D          # Diffusivity in beads
        k = self.k_bl       # mass transfer coeffient through liquid boundary layer
        K_tss = self.K_tss
        S_idx = list(*(D.nonzero()))
        D_ov_dz2 = D[S_idx]/(dz**2)     # (n_soluble,)
        D_ov_dz = D[S_idx]/dz
        _1_ov_z = 1/zs                  # (n_dz,)
        
        if self._fixed_P_gas:
            f_qgas = self.f_q_gas_fixed_P_headspace
        else:
            f_qgas = self.f_q_gas_var_P_headspace
        
        def dy_dt(t, y_ins, y, dy_ins):
            S_gas = y[-(n_gas+1):-1]
            Cs_bk = y[-(n_cmps+n_gas+1):-(n_gas+1)]          # bulk liquid concentrations
            Cs_en = y[:n_dz*n_cmps].reshape((n_dz, n_cmps))  # each row is one control volume
            Q_ins = y_ins[:, -1]
            Cs_ins = y_ins[:, :-1] * 1e-3  # mg/L to kg/m3
            Q = sum(Q_ins)
            
            # Transformation
            rhos_en = np.apply_along_axis(Rho_en, 1, Cs_en)
            Rs_en = rhos_en @ stoi_en       # n_dz * n_cmps
            rhos_bk = Rho_bk(y[-(n_cmps+n_gas+1):-1])
            Rs_bk = np.dot(stoi_bk.T, rhos_bk) # n_cmps+5
            q_gas = f_qgas(rhos_bk[-n_gas:], S_gas, T)
            gas_transfer = - q_gas*S_gas/V_gas + rhos_bk[-n_gas:] * V_liq/V_gas * gas_mass2mol_conversion
            
            # Detachment -- particulates
            tss = np.sum(Cs_en * (cmps.x*cmps.i_mass), axis=1)
            x_net_growth = np.sum(Rs_en * cmps.x, axis=1)/np.sum(Cs_en * cmps.x, axis=1) # d^(-1), equivalent to k_de
            u_de = 1/(1+np.exp(K_tss-tss)) * np.maximum(x_net_growth, 0) * (tss > 0)
            de_en = np.diag(u_de) @ (Cs_en * cmps.x)
            tot_de = np.sum(np.diag(dV) @ de_en, axis=0) / V_bead  # detachment per unit volume of beads

            #!!! Mass transfer (centered differences) -- MOL; solubles only
            C_lf = Cs_en[-1]
            J_lf = k*(Cs_bk - C_lf)
            S_en = Cs_en[:, S_idx]
            M_transfer = np.zeros_like(Cs_en)
            M_transfer[1:-1, S_idx] = D_ov_dz2 * (S_en[2:] - 2*S_en[1:-1] + S_en[:-2])\
                + D_ov_dz * (np.diag(_1_ov_z[1:-1]) @ (S_en[2:] - S_en[:-2]))
            M_transfer[0, S_idx] = 2 * D_ov_dz2 * (S_en[1] - S_en[0])
            M_transfer[-1, S_idx] = 2 * D_ov_dz2 * (S_en[-2] - S_en[-1])\
                + 2 * (1/dz + _1_ov_z[-1]) * J_lf[S_idx]
 
            # Mass balance
            dCdt_en = M_transfer + Rs_en - de_en
            dCdt_bk = (Q_ins @ Cs_ins - Q*Cs_bk)/V_liq \
                - A_beads/V_liq*J_lf + Rs_bk + V_beads*tot_de/V_liq
    
            _dstate[:n_dz*n_cmps] = dCdt_en.flatten()
            _dstate[-(n_cmps+n_gas+1):-(n_gas+1)] = dCdt_bk
            _dstate[-(n_gas+1):-1] = gas_transfer
            _dstate[-1] = sum(dy_ins[:,-1])
            _update_dstate()
        
        self._ODE = dy_dt

    def biomass_tss(self, biomass_IDs):
        '''Returns a 2-tuple of biomass TSS [kg/m3] in bulk and in encapsulation matrix (on average)'''
        y = self._state
        cmps = self.components
        n_cmps = len(cmps)
        n_dz = self.n_dz
        bm_idx = cmps.indices(biomass_IDs)
        en_bm = np.sum(y[:n_dz*n_cmps].reshape((n_dz, n_cmps))[:,bm_idx] * cmps.i_mass[bm_idx], axis=1)
        bk_bm = np.sum((y[n_dz*n_cmps: ((n_dz+1)*n_cmps)] * cmps.i_mass)[bm_idx])
        
        dz = self.r_beads / n_dz
        zs = np.linspace(dz, self.r_beads, n_dz)
        dV = 4/3*np.pi*(zs)**3
        V_bead = dV[-1]
        dV[1:] -= dV[:-1]
        C_en_avg = np.dot(en_bm, dV)/V_bead
        return bk_bm, C_en_avg
    
    def get_retained_mass(self, biomass_IDs):
        bk_bm, C_en_avg = self.biomass_tss(biomass_IDs)
        return self.V_liq * bk_bm + self.V_beads * C_en_avg

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
    
    T_air = 273.15 + 22
    T_earth = 273.15 + 22
    
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
        }
    
    def _static_lift_equivalent(self):
        dia = (self.V_bed*4/self.reactor_height_to_diameter/pi) ** (1/3)
        h = dia * self.reactor_height_to_diameter
        return h
    
    def _design(self):
        D = self.design_results
        den = self._density
        V = D['Volume'] = self.V_bed
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
        rQ = self.recirculation_ratio or 0
        for ws in self.ins:
            _inch, _kg_per_m = pipe_design(ws.F_vol, self._l_min_velocity)
            pipe_IDs.append(_inch)
            HDPE_pipes.append(_kg_per_m*L_inlets)
        if rQ > 0:
            _inch, _kg_per_m = pipe_design(mixed.F_vol*rQ, self._l_min_velocity)
            pipe_IDs.append(_inch)
            HDPE_pipes.append(_kg_per_m*(L_inlets+L_outlets))
        for ws in self.outs:
            if ws.phase == 'g':
                D['Stainless steel'] += pipe_design(ws.F_vol, self._g_min_velocity, True)[1] * L_gas                
            else:
                _inch, _kg_per_m = pipe_design(ws.F_vol, self._l_min_velocity)
                pipe_IDs.append(_inch)
                HDPE_pipes.append(_kg_per_m*L_outlets)
        D['HDPE pipes'] = sum(HDPE_pipes)
        self._hdpe_ids, self._hdpe_kgs = pipe_IDs, HDPE_pipes
        
        # Pumps
        flowsheet_ID = self.system.flowsheet.ID
        getfield = getattr
        creg = Construction.registry
        if self.fixed_headspace_P:
            gas = self.outs[0]
            self.vacuum_pump = vac = VacuumSystem(
                self, F_mass=_F_mass(gas), F_vol=_F_vol(gas), 
                P_suction=self.headspace_P * 1e5,
                vessel_volume=self.V_gas
                )
            qvac = _construct_vacuum_pump(vac)
            pvac = getfield(creg, f'{flowsheet_ID}_{self.ID}_VacPump_surrogate')
            pvac.quantity = qvac
        h_lift = self._static_lift_equivalent()
        for i, ws, in enumerate(self.ins):
            ID = pipe_IDs[i]
            hf = pipe_friction_head(ws.F_vol*_cmph_2_gpm, L_inlets, ID)  # friction head loss
            TDH = hf + h_lift # in m, assume suction head = 0, discharge head = reactor height
            field = f'Pump_ins{i}'
            pump = getfield(self, field)
            pump.ins[0].copy_flow(ws)
            # pump.ins[0].set_total_flow(ws.F_vol*(1+rQ), 'm3/h')
            pump.dP_design = TDH * 9804.14  # in Pa
            pump.simulate()
            if self.include_construction:
                q22, q40 = _construct_water_pump(pump)
                p22 = getfield(creg, f'{flowsheet_ID}_{pump.ID}_22kW')
                p40 = getfield(creg, f'{flowsheet_ID}_{pump.ID}_40W')
                p22.quantity = q22
                p40.quantity = q40
        if rQ > 0:
            field = 'Pump_recirculation'
            if not hasattr(self, field):
                pump = Pump(self.ID+'recir_Pump', ins=Stream(f'{self.ID}_recir'))
                setattr(self, field, pump)
                self.construction += [
                    Construction(ID='22kW', linked_unit=pump, item='pump_22kW'),
                    Construction(ID='40W', linked_unit=pump, item='pump_40W')
                    ]
                self.auxiliary_unit_names = tuple({*self.auxiliary_unit_names, field})
            ID = pipe_IDs[len(self.ins)]
            hf = pipe_friction_head(mixed.F_vol*rQ*_cmph_2_gpm, L_inlets+L_outlets, ID)  # friction head loss
            TDH = hf # in m, assume suction head = 0, static lift = 0
            pump = self.Pump_recirculation
            pump.ins[0].copy_flow(self.outs[1])
            pump.ins[0].set_total_flow(mixed.F_vol*rQ, 'm3/h')
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
        D['Bead volume'] = self.V_beads
        self.add_equipment_design()

    _NG_price = 0.85*auom('kJ').conversion_factor('therm') # [USD/kJ] 5.47 2021USD/Mcf vs. 4.19 2016USD/Mcf    

    _cost = UASB._cost
    add_equipment_design = UASB.add_equipment_design

#%% METAB_PackedBed

class METAB_PackedBed(METAB_FluidizedBed):
    
    def __init__(self, ID='', voidage=0.39, n_cstr=None, **kwargs):
        self._n_cstr = n_cstr
        self._checked_ncstr = False
        V_liq = kwargs.pop('V_liq', 3400)
        V_gas = kwargs.pop('V_gas', 300)
        V_tot = V_liq/0.38 + V_gas
        V_gas = V_tot - V_liq/voidage
        if V_gas <= 0: V_gas = 0.1
        super().__init__(ID=ID, V_liq=V_liq, V_gas=V_gas, voidage=voidage, **kwargs)

    voidage = property(METAB_FluidizedBed.voidage.fget)
    @voidage.setter
    def voidage(self, f):
        if f < 0.35 or f > 0.45:
            raise ValueError(f'voidage must be in [0.35, 0.45], not {f}')
        if hasattr(self, 'f_void'):
            f_old = self.f_void
            Vg = self.V_gas
            Vg_subtract = self.V_liq/(1/f - 1/f_old)
            self.V_gas = max(0.1, Vg-Vg_subtract)
        self.f_void = f
    
    @property
    def recirculation_ratio(self):
        return self._rQ or 0.
    @recirculation_ratio.setter
    def recirculation_ratio(self, r):
        self._rQ = r
    
    @property
    def n_cstr(self):
        if not self._checked_ncstr:
            n = self._n_cstr
            self.n_cstr = n
        return self._n_cstr
    @n_cstr.setter
    def n_cstr(self, n):
        if self.recirculation_ratio:
            if self.recirculation_ratio >= 1:
                self._n_cstr = 1
                warn('n_cstr is forced to be 1 due to recirculation ratio >= 1')
        else:
            n = n or ceil(self.reactor_height_to_diameter)
            self._n_cstr = int(max(n, 2))
        self._checked_ncstr = True
        
    def _prep_model(self):
        model = self.model
        self._S_vapor = self.ideal_gas_law(p=self.p_vapor())
        self._n_gas = len(model._biogas_IDs)
        cmps = self.thermo.chemicals
        layers = [*range(self.n_dz), 'bulk']
        self._state_keys = [f'{cmp}-R{i}-{j}' for i in range(self.n_cstr) 
                            for j in layers for cmp in cmps.IDs] \
            + [ID+'_gas' for ID in self.model._biogas_IDs] \
            + ['Q']
        self._gas_cmp_idx = cmps.indices(self.model._biogas_IDs)
        self._state_header = self._state_keys
        self._gas_state_idx = dict(zip(self.model._biogas_IDs, range(self._n_gas)))
    
    def set_init_conc(self, arr=None, **kwargs):
        '''set the initial concentrations [kg/m3] of components in the packed bed, 
        applies uniformly to all layers in beads and bulk liquid and all well-mixed
        compartments unless an array is input.'''
        cmps = self.thermo.chemicals
        cmpx = cmps.index
        n_dz = self.n_dz
        n_cstr = self.n_cstr
        if arr is None:
            Cs = np.zeros(len(cmps))
            for k, v in kwargs.items(): Cs[cmpx(k)] = v
            self._concs = np.tile(Cs, (n_dz+1)*n_cstr)
        else:
            arr = np.asarray(arr)
            if arr.shape != (len(cmps)*(n_dz+1)*n_cstr,):
                raise ValueError(f'arr must be None or a 1d array of length {len(cmps)*(n_dz+1)}, '
                                 f'not {arr.shape}')
            self._concs = arr
    
    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        if self._cached_state is not None: Cs = self._cached_state
        elif self._concs is not None: Cs = self._concs
        else: Cs = np.tile(mixed.conc, (self.n_dz+1)*self.n_cstr)
        self._state = np.append(Cs, [0]*self._n_gas + [Q]).astype('float64')
        self._dstate = self._state * 0.
    
    def f_q_gas_fixed_P_headspace(self, rhoTs, S_gas, T):
        cmps = self.components
        gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]
        q_gas = self._R*T/(self._P_gas-self.p_vapor(convert_to_bar=True))\
                                *self.V_liq/self.n_cstr*sum(rhoTs*gas_mass2mol_conversion)
        return q_gas
    
    def f_q_gas_var_P_headspace(self, rhoTs, S_gas, T):
        p_gas = S_gas * self._R * T
        self._P_gas = P = sum(p_gas) + self.p_vapor(convert_to_bar=True) 
        q_gas = max(0, self._k_p * (P - self._P_atm))
        return q_gas
    
    def _compile_ODE(self):
        cmps = self.components
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        n_cmps = len(cmps)
        n_cstr = self.n_cstr
        n_dz = self.n_layer
        n_gas = self._n_gas
        T = self.T
        
        _f_rhos = self.model.flex_rate_function
        _params = self.model.rate_function.params
        stoi_bk = self.model.stoichio_eval()
        stoi_en = stoi_bk[:-n_gas]  # no liquid-gas transfer
        Rho_bk = lambda state_arr: _f_rhos(state_arr, _params, 
                                           T_op=T, pH=self.pH_ctrl, 
                                           gas_transfer=True)
        Rho_en = lambda state_arr: _f_rhos(state_arr, _params, 
                                           T_op=T, pH=self.pH_ctrl, 
                                           gas_transfer=False)
        gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]

        V_liq = self.V_liq/n_cstr
        V_gas = self.V_gas
        V_beads = self.V_beads/n_cstr
        r_beads = self.r_beads
        A_beads = 3 * V_beads / r_beads # m2, total bead surface area per cstr

        dz = r_beads / n_dz
        zs = np.linspace(dz, r_beads, n_dz)
        dV = 4/3*np.pi*(zs)**3
        dV[1:] -= dV[:-1]
        V_bead = (4/3*np.pi*r_beads**3)
        
        D = self.D          # Diffusivity in beads
        k = self.k_bl       # mass transfer coeffient through liquid boundary layer
        K_tss = self.K_tss
        S_idx = list(*(D.nonzero()))
        D_ov_dz2 = D[S_idx]/(dz**2)     # (n_soluble,)
        D_ov_dz = D[S_idx]/dz
        _1_ov_z = 1/zs                  # (n_dz,)
        rhos_gas = np.zeros((n_cstr, n_gas))

        if self._fixed_P_gas:
            f_qgas = self.f_q_gas_fixed_P_headspace
        else:
            f_qgas = self.f_q_gas_var_P_headspace
        
        def dy_dt(t, y_ins, y, dy_ins):
            S_gas = y[-(n_gas+1):-1]
            Cs_all = y[:-(n_gas+1)].reshape((n_cstr, n_cmps*(n_dz+1)))
            Q_ins = y_ins[:, -1]
            Cs_ins = y_ins[:, :-1] * 1e-3  # mg/L to kg/m3
            flow_in = Q_ins @ Cs_ins
            Q = sum(Q_ins)
            for i, yi in enumerate(Cs_all):
                Cs_bk = yi[-n_cmps:]                                # bulk liquid concentrations
                Cs_en = yi[:n_dz*n_cmps].reshape((n_dz, n_cmps))     # each row is one control volume

                # Transformation
                rhos_en = np.apply_along_axis(Rho_en, 1, Cs_en)
                Rs_en = rhos_en @ stoi_en                           # n_dz * n_cmps
                rhos_bk = Rho_bk(np.append(Cs_bk, S_gas))
                Rs_bk = np.dot(stoi_bk.T, rhos_bk)
                rhos_gas[i,:] = rhos_bk[-n_gas:].copy()
            
                # Detachment -- particulates
                tss = np.sum(Cs_en * (cmps.x*cmps.i_mass), axis=1)
                x_net_growth = np.sum(Rs_en * cmps.x, axis=1)/np.sum(Cs_en * cmps.x, axis=1) # d^(-1), equivalent to k_de
                u_de = 1/(1+np.exp(K_tss-tss)) * np.maximum(x_net_growth, 0) * (tss > 0)
                de_en = np.diag(u_de) @ (Cs_en * cmps.x)
                tot_de = np.sum(np.diag(dV) @ de_en, axis=0) / V_bead  # detachment per unit volume of beads
    
                #!!! Mass transfer (centered differences) -- MOL; solubles only
                C_lf = Cs_en[-1]
                J_lf = k*(Cs_bk - C_lf)
                S_en = Cs_en[:, S_idx]
                M_transfer = np.zeros_like(Cs_en)
                M_transfer[1:-1, S_idx] = D_ov_dz2 * (S_en[2:] - 2*S_en[1:-1] + S_en[:-2])\
                    + D_ov_dz * (np.diag(_1_ov_z[1:-1]) @ (S_en[2:] - S_en[:-2]))
                M_transfer[0, S_idx] = 2 * D_ov_dz2 * (S_en[1] - S_en[0])
                M_transfer[-1, S_idx] = 2 * D_ov_dz2 * (S_en[-2] - S_en[-1])\
                    + 2 * (1/dz + _1_ov_z[-1]) * J_lf[S_idx]
     
                # Mass balance
                dCdt_en = M_transfer + Rs_en - de_en
                dCdt_bk = (flow_in - Q*Cs_bk)/V_liq \
                    - A_beads/V_liq*J_lf + Rs_bk + V_beads*tot_de/V_liq
                
                start = i*((n_dz+1)*n_cmps)
                stop_en = start + n_dz*n_cmps
                stop_bk = stop_en + n_cmps
                _dstate[start:stop_en] = dCdt_en.flatten()
                _dstate[stop_en:stop_bk] = dCdt_bk
                
                flow_in = Q * Cs_bk             # for next CSTR
            
            q_gas = np.apply_along_axis(f_qgas, 1, rhos_gas, S_gas, T)
            self._q_gas = sum(q_gas)
            d_gas = - sum(q_gas)*S_gas/V_gas + np.sum(rhos_gas*V_liq, axis=0)/V_gas * gas_mass2mol_conversion
            _dstate[-(n_gas+1):-1] = d_gas
            _dstate[-1] = sum(dy_ins[:,-1])
            _update_dstate()
        
        self._ODE = dy_dt

    def biomass_tss(self, biomass_IDs):
        '''Returns a 2-tuple of biomass TSS [kg/m3] in bulk and in encapsulation matrix (on average)'''
        y = self._state
        cmps = self.components
        n_cmps = len(cmps)
        n_dz = self.n_dz
        n_cstr = self.n_cstr
        bm_idx = cmps.indices(biomass_IDs)
        outs = []
        dz = self.r_beads / n_dz
        zs = np.linspace(dz, self.r_beads, n_dz)
        dV = 4/3*np.pi*(zs)**3
        V_bead = dV[-1]
        dV[1:] -= dV[:-1]
        for i in range(n_cstr):
            start = i*((n_dz+1)*n_cmps)
            stop_en = start + n_dz*n_cmps
            stop_bk = stop_en + n_cmps
            en_bm = np.sum(y[start: stop_en].reshape((n_dz, n_cmps))[:,bm_idx] * cmps.i_mass[bm_idx], axis=1)
            bk_bm = np.sum((y[stop_en: stop_bk] * cmps.i_mass)[bm_idx])
            C_en_avg = np.dot(en_bm, dV)/V_bead
            outs.append([bk_bm, C_en_avg])
        outs = np.asarray(outs)
        return np.mean(outs, axis=0)

    def _static_lift_equivalent(self):
        dia = (self.V_bed*4/self.reactor_height_to_diameter/pi) ** (1/3)
        L = dia * self.reactor_height_to_diameter
        void = self.voidage
        mixed = self._mixed
        rho = mixed.rho
        mu = mixed.mu
        d = self.bead_diameter * 1e-3
        A_bed = (pi*self.V_bed**2/4/self.reactor_height_to_diameter**2)**(1/3)
        A_liq = A_bed * void
        u = mixed.F_vol/A_liq/3600     # m/s
        dP = L * (150 * mu * (1-void)**2 * u / (void**3 * d**2) \
            + 1.75 * (1-void) * rho * u**2 / (void**3 * d))
        return dP / (9.81 * rho) + L # Pa to m
    
#%% Batch experiment

class METAB_BatchExp(METAB_FluidizedBed):
    _N_ins = 0
    _N_outs = 1
    
    def __init__(self, ID='', outs=(), thermo=None,
                 init_with='WasteStream', V_liq=25, V_gas=30, 
                 V_beads=10, bead_diameter=10, n_layer=5,
                 boundary_layer_thickness=0.01, diffusivity=None,
                 f_diff=0.55, max_encapsulation_tss=16, model=None,
                 pH_ctrl=False, T=295.15, headspace_P=1.013, external_P=1.013, 
                 pipe_resistance=5.0e-1, fixed_headspace_P=True,
                 isdynamic=True, exogenous_vars=(), **kwargs):   
        
        SanUnit.__init__(self, ID=ID, ins=None, outs=outs, thermo=thermo,
                         init_with=init_with, isdynamic=isdynamic, 
                         exogenous_vars=exogenous_vars, **kwargs)
        self.V_gas = V_gas
        self.V_liq = V_liq
        self.V_beads = V_beads
        self.bead_diameter = bead_diameter
        self.n_layer = n_layer
        self.boundary_layer_thickness = boundary_layer_thickness
        self.diffusivity = diffusivity
        self.f_diff = f_diff
        self.max_encapsulation_tss = max_encapsulation_tss
        
        self.pH_ctrl = pH_ctrl
        self.T = T
        self._q_gas = 0
        self._n_gas = None
        self._gas_cmp_idx = None
        self._state_keys = None
        self._S_vapor = None
        self._biogas = WasteStream(phase='g')
        self.headspace_P = headspace_P
        self.external_P = external_P
        self.pipe_resistance = pipe_resistance
        self.fixed_headspace_P = fixed_headspace_P
        self._tempstate = []
        
        self.model = model
        self._cached_state = None

    @property
    def V_beads(self):
        return self._V_beads
    @V_beads.setter
    def V_beads(self, V):
        self._V_beads = V

    def _setup(self):
        AnaerobicCSTR._setup(self)
            
    def _run(self):
        pass

    def _set_init_Cs(self, arr=None, **kwargs):
        cmps = self.thermo.chemicals
        cmpx = cmps.index
        if arr is None:
            Cs = np.zeros(len(cmps))
            for k, v in kwargs.items(): Cs[cmpx(k)] = v
            return Cs
        else:
            arr = np.asarray(arr)
            if arr.shape != (len(cmps),):
                raise ValueError(f'arr must be None or a 1d array of length {len(cmps)}, '
                                 f'not {arr.shape}')
            return arr

    def set_bulk_init(self, arr=None, **kwargs):
        self._bulk_concs = self._set_init_Cs(arr, **kwargs)
    
    def set_encap_init(self, arr=None, **kwargs):
        n_dz = self.n_dz
        Cs = self._set_init_Cs(arr, **kwargs)
        self._encap_concs = np.tile(Cs, n_dz)

    def _init_state(self):
        if self._cached_state is not None: 
            Cs = self._cached_state
            Cs[-len(self._bulk_concs):] = self._bulk_concs
        else: Cs = np.append(self._encap_concs, self._bulk_concs)
        self._state = np.append(Cs, [0]*self._n_gas + [0]).astype('float64')
        self._dstate = self._state * 0.
        
    def _update_state(self):
        cmps = self.components
        n_cmps = len(cmps)
        n_gas = self._n_gas
        y_gas = self._state[-(n_gas+1):-1]
        i_mass = cmps.i_mass
        chem_MW = cmps.chem_MW        
        gas, = self._outs
        if gas.state is None:
            gas.state = np.zeros(n_cmps+1)
        gas.state[self._gas_cmp_idx] = y_gas
        gas.state[cmps.index('H2O')] = self._S_vapor
        gas.state[-1] = self._q_gas
        gas.state[:n_cmps] = gas.state[:n_cmps] * chem_MW / i_mass * 1e3 # i.e., M biogas to mg (measured_unit) / L

    def _update_dstate(self):
        self._tempstate = self.model.rate_function._params['root'].data.copy()        
        gas, = self._outs
        if gas.dstate is None:
            # contains no info on dstate
            n_cmps = len(self.components)
            gas.dstate = np.zeros(n_cmps+1)
            
    def _compile_ODE(self):
        cmps = self.components
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        n_cmps = len(cmps)
        n_dz = self.n_layer
        n_gas = self._n_gas
        T = self.T
        
        _f_rhos = self.model.flex_rate_function
        _params = self.model.rate_function.params
        stoi_bk = self.model.stoichio_eval()
        stoi_en = stoi_bk[:-n_gas]  # no liquid-gas transfer
        Rho_bk = lambda state_arr: _f_rhos(state_arr, _params, 
                                           T_op=T, pH=self.pH_ctrl, 
                                           gas_transfer=True)
        Rho_en = lambda state_arr: _f_rhos(state_arr, _params, 
                                           T_op=T, pH=self.pH_ctrl, 
                                           gas_transfer=False)
        gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]

        V_liq = self.V_liq
        V_gas = self.V_gas
        V_beads = self.V_beads
        r_beads = self.r_beads
        A_beads = 3 * V_beads / r_beads # m2, total bead surface area

        dz = r_beads / n_dz
        zs = np.linspace(dz, r_beads, n_dz)
        dV = 4/3*np.pi*(zs)**3
        dV[1:] -= dV[:-1]
        V_bead = (4/3*np.pi*r_beads**3)
        
        D = self.D          # Diffusivity in beads
        k = self.k_bl       # mass transfer coeffient through liquid boundary layer
        K_tss = self.K_tss
        S_idx = list(*(D.nonzero()))
        D_ov_dz2 = D[S_idx]/(dz**2)     # (n_soluble,)
        D_ov_dz = D[S_idx]/dz
        _1_ov_z = 1/zs                  # (n_dz,)
        
        if self._fixed_P_gas:
            f_qgas = self.f_q_gas_fixed_P_headspace
        else:
            f_qgas = self.f_q_gas_var_P_headspace
        
        def dy_dt(t, y_ins, y, dy_ins):
            S_gas = y[-(n_gas+1):-1]
            Cs_bk = y[-(n_cmps+n_gas+1):-(n_gas+1)]          # bulk liquid concentrations
            Cs_en = y[:n_dz*n_cmps].reshape((n_dz, n_cmps))  # each row is one control volume
            
            # Transformation
            rhos_en = np.apply_along_axis(Rho_en, 1, Cs_en)
            Rs_en = rhos_en @ stoi_en       # n_dz * n_cmps
            rhos_bk = Rho_bk(y[-(n_cmps+n_gas+1):-1])
            Rs_bk = np.dot(stoi_bk.T, rhos_bk) # n_cmps+5
            q_gas = f_qgas(rhos_bk[-n_gas:], S_gas, T)
            gas_transfer = - q_gas*S_gas/V_gas + rhos_bk[-n_gas:] * V_liq/V_gas * gas_mass2mol_conversion
            
            # Detachment -- particulates
            tss = np.sum(Cs_en * (cmps.x*cmps.i_mass), axis=1)
            x_net_growth = np.sum(Rs_en * cmps.x, axis=1)/np.sum(Cs_en * cmps.x, axis=1) # d^(-1), equivalent to k_de
            u_de = 1/(1+np.exp(K_tss-tss)) * np.maximum(x_net_growth, 0) * (tss > 0)
            de_en = np.diag(u_de) @ (Cs_en * cmps.x)
            tot_de = np.sum(np.diag(dV) @ de_en, axis=0) / V_bead  # detachment per unit volume of beads

            #!!! Mass transfer (centered differences) -- MOL; solubles only
            C_lf = Cs_en[-1]
            J_lf = k*(Cs_bk - C_lf)
            S_en = Cs_en[:, S_idx]
            M_transfer = np.zeros_like(Cs_en)
            M_transfer[1:-1, S_idx] = D_ov_dz2 * (S_en[2:] - 2*S_en[1:-1] + S_en[:-2])\
                + D_ov_dz * (np.diag(_1_ov_z[1:-1]) @ (S_en[2:] - S_en[:-2]))
            M_transfer[0, S_idx] = 2 * D_ov_dz2 * (S_en[1] - S_en[0])
            M_transfer[-1, S_idx] = 2 * D_ov_dz2 * (S_en[-2] - S_en[-1])\
                + 2 * (1/dz + _1_ov_z[-1]) * J_lf[S_idx]
 
            # Mass balance
            dCdt_en = M_transfer + Rs_en - de_en
            dCdt_bk = - A_beads/V_liq*J_lf + Rs_bk + V_beads*tot_de/V_liq
    
            _dstate[:n_dz*n_cmps] = dCdt_en.flatten()
            _dstate[-(n_cmps+n_gas+1):-(n_gas+1)] = dCdt_bk
            _dstate[-(n_gas+1):-1] = gas_transfer
            _update_dstate()
        
        self._ODE = dy_dt

    def _design(self):
        pass
    
    def _cost(self):
        pass