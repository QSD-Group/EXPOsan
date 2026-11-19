# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Focused native BioSTEAM extensions for conversion and product separation."""
from math import ceil, pi
import numpy as np
import biosteam as bst
from biosteam.units.design_tools import PressureVessel, CEPCI_by_year
from ._thermal_basis import CostBasis
from ._thermal_chemistry import ELEMENTS, atom_vector, gasify, complete_gas_balance, FEED_FORMULAS

PNNL = 'https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-23227.pdf'
NREL = 'https://digital.library.unt.edu/ark:/67531/metadc1015731/m2/1/high_res_d/993332.pdf'
SEIDER = 'https://biosteam.readthedocs.io/en/latest/_modules/biosteam/units/design_tools/pressure_vessel.html'

# https://www.osti.gov/servlets/purl/1081663 : Appendix E attributes Wright's
# equipment basis to 2007 USD; NREL Appendix A Table A-1: 500 dry metric t/day
# per fluid bed, purchase 836000 USD, installed/purchase = 3.02.
PYROLYZER = CostBasis(836000., 500000/24, 'kg dry feed/h', 2007,
    525.4, NREL, 'Appendix A, Table A-1/A-2, pyrolysis fluid beds',
    .6, installation_factor=3.02,
    scaling_basis='engineering six-tenths scaling assumption; source gives reference equipment only',
    scope='reactor only; excludes cyclones, condensation, feed preparation and utilities')

# https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-23227.pdf
# Appendix B, PSA-310: 1.75 million USD purchase, 10 mmscfd hydrogen, 2004 USD.
HYDROGEN_PSA = CostBasis(1750000., 10., 'mmscfd H2 at 60 F, 1 atm', 2004,
    444.2, PNNL, 'Appendix B, PSA-310, hydrogen recycle PSA',
    .8, installation_factor=2.47,
    scope='hydrogen separation package; product compression is separate')

# https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-23227.pdf
# Appendix B SP240A, Kreb's Feb 2009 quote: 5 million USD at 968859 lb/h.
HYDROCYCLONE = CostBasis(5e6, 968859*.45359237, 'kg slurry/h', 2009,
    521.9, PNNL, 'Appendix B SP240A hydrocyclone, Kreb Engineering quote',
    .65, installation_factor=2.1,
    scope='high-pressure solids-removal proxy, not a validated supercritical salt separator')


def equipment_basis(config, key, default):
    basis = config.cost_overrides.get(key,default)
    if isinstance(basis,dict): basis = CostBasis(**basis)
    if not isinstance(basis,CostBasis) or basis.size_unit != default.size_unit:
        raise ValueError(f'{key} cost size_unit must be {default.size_unit!r}; convert the source size explicitly')
    return basis


def water_h(T, P):
    """Pure-water enthalpy, kJ/kg, with real-fluid pressure dependence."""
    from CoolProp.CoolProp import PropsSI
    return PropsSI('Hmass', 'T', float(T), 'P', float(P), 'Water')/1000


def sensible_h(stream, T=None):
    """ThermoSTEAM mixture H with a CoolProp correction to its water term."""
    T = stream.T if T is None else T
    clone = stream.copy()
    clone.T = T
    if not stream.imass['Water']:
        return float(clone.H)
    pure = bst.Stream(thermo=stream.thermo, Water=1., units='kg/hr', T=T, P=stream.P,
                      phase=stream.phase if len(stream.phase) == 1 else 'l')
    # https://coolprop.org/fluid_properties/fluids/Water.html : IAPWS water EOS.
    # Reference at 298.15 K/1 atm aligns the sensible-enthalpy zero.
    if stream.phase == 'g' and T < 647.096:
        from CoolProp.CoolProp import PropsSI
        # https://coolprop.org/fluid_properties/fluids/Water.html : wet-gas
        # water is vapor at its partial pressure, not compressed liquid at total
        # pressure. Saturated-vapor h(T) is the low-temperature gas approximation.
        if stream.P >= float(stream.chemicals.Water.Psat(T)):
            real_water = PropsSI('Hmass','T',T,'Q',1,'Water')/1000
        else:
            real_water = water_h(T,stream.P)
    else:
        real_water = water_h(T, stream.P)
    real_water -= water_h(298.15,101325.)
    return float(clone.H)+float(stream.imass['Water'])*(real_water-float(pure.H))


class WaterHX(bst.HXutility):
    """Native exchanger/furnace sizing with real-fluid water heat duty."""
    def _design(self, duty=None):
        if duty is None:
            duty = sensible_h(self.outs[0])-sensible_h(self.ins[0])
        super()._design(duty=duty)


class OilCondenser(WaterHX):
    """Native exchanger budget including a declared oil latent-heat surrogate."""
    def _init(self, latent_kj_kg, **kwargs):
        self.latent_kj_kg = latent_kj_kg
        super()._init(**kwargs)

    def _design(self, duty=None):
        latent = float(self.ins[0].imass['PyrolysisOil'])*self.latent_kj_kg if 'PyrolysisOil' in self.chemicals else 0.
        duty = sensible_h(self.outs[0])-sensible_h(self.ins[0])-latent
        super()._design(duty=duty)
        self.design_results['Organic condensation duty kJ/h'] = latent


class BoundedCentrifuge(bst.SolidsCentrifuge):
    """Native centrifuge with cake moisture bounded by available feed water."""
    def _run(self):
        requested = self.moisture_content
        water = sum(float(s.imass['Water']) for s in self.ins)
        dry = sum(sum(float(s.imass[c.ID])*float(self.isplit[c.ID])
                  for c in self.chemicals if c.ID != 'Water') for s in self.ins)
        attainable = water/(water+dry) if water+dry else 0.
        # https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/separations.html
        # Bound the specified cake water before native separation; a dry incoming
        # feed cannot gain water from an unspecified source. Strict mode avoids
        # the native non-strict moisture correction's negative-permeate branch.
        self.moisture_content = min(requested,attainable*(1-1e-12))
        try:
            super()._run()
        finally:
            self.moisture_content = requested


class WaterRecuperator(bst.HXprocess):
    """Countercurrent heat recovery; native LMTD/area/cost implementation."""
    # https://www.psi.ch/sites/default/files/import/cpe/PosterCPEEN/Poster_Erlangen.pdf
    # 75% exchanger efficiency is a process-study assumption, not measured UA.
    recovery_fraction = .75

    def _run(self):
        from scipy.optimize import brentq
        cold, hot = self.ins
        co, ho = self.outs
        co.copy_like(cold)
        ho.copy_like(hot)
        if not cold.F_mass or not hot.F_mass or hot.T <= cold.T+self.dT:
            self.total_heat_transfer = 0.
            return
        hc, hh = sensible_h(cold), sensible_h(hot)
        q = self.recovery_fraction*min(sensible_h(cold,hot.T-self.dT)-hc,
                                      hh-sensible_h(hot,cold.T+self.dT))
        q = max(float(q),0.)
        if q:
            co.T = brentq(lambda T:sensible_h(cold,T)-hc-q, cold.T,hot.T-self.dT)
            ho.T = brentq(lambda T:sensible_h(hot,T)-hh+q, cold.T+self.dT,hot.T)
        self.total_heat_transfer = q


class Conversion(bst.Unit, PressureVessel):
    """Element-conserving conversion scenario with native pressure-vessel budget."""
    _N_ins = 1
    _N_outs = 1
    _units = {**PressureVessel._units, 'Total volume':'m3', 'Single volume':'m3',
              'Residence time':'hr', 'Dry throughput':'kg/hr'}
    _F_BM_default = dict(PressureVessel._F_BM_default)
    _vessel_material = 'Stainless steel 316'
    _vessel_type = 'Vertical'

    def _init(self, technology, config, composition, products, **kwargs):
        self.technology = technology
        self.config = config
        self.reference_composition = composition
        self.products = products
        self.vessel_material = 'Stainless steel 316'
        self.cost_audit = []

    def _run(self):
        feed, = self.ins
        out, = self.outs
        out.copy_like(feed)
        out.T, out.P = self.config.temperature_k, self.config.pressure_pa
        actual = np.array([float(feed.imass[k]) for k in FEED_FORMULAS])
        reference = np.array([self.reference_composition.get(k,0.) for k in FEED_FORMULAS])
        if actual.sum() and not np.allclose(actual/actual.sum(),reference/reference.sum(),atol=1e-10,rtol=1e-8):
            raise ValueError('Changed organic composition: rebuild the thermochemical system/property package')
        atoms = sum((float(feed.imol[k])*atom_vector(self.chemicals[k]) for k in FEED_FORMULAS),np.zeros(5))
        water = float(feed.imol['Water'])
        if self.technology == 'htg':
            converted = atoms*self.config.gasification_conversion
            for k in FEED_FORMULAS:
                out.imol[k] *= 1-self.config.gasification_conversion
            gases = gasify(converted, self.config.methanation_completion,self.config.co_carbon_fraction)
            remaining_water = water-gases.pop('reaction_water')
            if remaining_water < -1e-9:
                raise ValueError('Insufficient reaction water for selected HTG conversion')
            out.imol['Water'] = max(remaining_water,0.)
            for name, value in gases.items():
                out.imol[name] = value
        else:
            for k in FEED_FORMULAS:
                out.imol[k] = 0.
            residual = atoms.copy()
            for ID, retention in (
                ('PyrolysisOil',self.config.oil_retention),
                ('BiocharMatrix' if self.technology == 'fast_pyrolysis' else 'HydrocharMatrix', self.config.char_retention),
                ('AqueousOrganics',self.config.aqueous_retention)):
                if ID not in self.chemicals:
                    continue
                retained = atoms*np.array([retention.get(e,0.) for e in ELEMENTS])
                formula = atom_vector(self.chemicals[ID])
                index = int(np.argmax(formula))
                out.imol[ID] = retained[index]/formula[index]
                residual -= retained
            for name, value in complete_gas_balance(residual,water).items():
                out.imol[name] = value
        if out.mol.min() < -1e-9:
            raise ValueError('Negative component flow after conversion')

    def _design(self):
        c = self.config
        feed = self.ins[0]
        d = self.design_results
        d['Dry throughput'] = float(feed.F_mass-feed.imass['Water'])
        if self.technology != 'fast_pyrolysis':
            from CoolProp.CoolProp import PropsSI
            rho = PropsSI('Dmass','T',c.temperature_k,'P',c.pressure_pa,'Water')
            hot_volume = float(feed.imass['Water'])/rho
            hot_volume += sum(float(feed.imass[k])/1400 for k in FEED_FORMULAS)
            total = hot_volume*c.residence_time_h/c.working_volume_fraction
            number = max(ceil(total/c.max_vessel_m3),1)
            volume = total/number
            # https://biosteam.readthedocs.io/en/latest/_modules/biosteam/units/design_tools/pressure_vessel.html
            # L/D=3 is an engineering geometry scenario; native vessel thickness
            # uses operating pressure plus its own design margin and corrosion.
            diameter = (4*volume/(pi*3))**(1/3)
            length = 3*diameter
            d.update(self._vessel_design(c.pressure_pa/6894.757293, diameter/0.3048,length/0.3048))
            d.update({'Total volume':total,'Single volume':volume,'Number of vessels':number,
                      'Residence time':c.residence_time_h,'Hot water density kg/m3':rho})
        else:
            d['Residence time'] = c.residence_time_h
        # Heat of reaction and real-fluid water correction are included once.
        duty = sensible_h(self.outs[0])-sensible_h(feed)+self.outs[0].Hf-feed.Hf
        if self.technology == 'fast_pyrolysis':
            # https://chemicals.readthedocs.io/chemicals.phase_change.html
            # Supply organic vaporization once here and remove the same latent
            # energy once in OilCondenser; the lump remains an empirical species.
            latent = (float(self.outs[0].imass['PyrolysisOil'])*self.products['PyrolysisOil']['latent_kj_kg']
                      if 'PyrolysisOil' in self.products else 0.)
            d['Organic vaporization duty kJ/h'] = latent
            duty += latent
        self.add_heat_utility(duty,c.temperature_k,c.temperature_k)

    def _cost(self):
        self.cost_audit = []
        if self.technology == 'fast_pyrolysis':
            basis = equipment_basis(self.config,'pyrolyzer',PYROLYZER)
            audit = basis.estimate(self.design_results['Dry throughput'],bst.CE,self.config.cost_year)
            self.cost_audit.append(audit)
            if audit['purchase_usd'] is None:
                self.installed_costs['Pyrolyzer'] = audit['installed_usd']
            else:
                self.baseline_purchase_costs['Pyrolyzer'] = audit['purchase_usd']
                self.F_BM['Pyrolyzer'] = basis.installation_factor
        else:
            d = self.design_results
            if 'reactor' in self.config.cost_overrides:
                basis = self.config.cost_overrides['reactor']
                if isinstance(basis,dict): basis = CostBasis(**basis)
                if not isinstance(basis,CostBasis) or basis.size_unit != 'm3 per vessel':
                    raise ValueError('Reactor quote size_unit must be m3 per vessel')
                audit = basis.estimate(d['Single volume'],bst.CE,self.config.cost_year)
                number = d['Number of vessels']
                audit.update(number_of_vessels=number,installed_total_usd=audit['installed_usd']*number)
                self.cost_audit = [audit]
                if audit['purchase_usd'] is None:
                    self.installed_costs['Quoted reactor package'] = audit['installed_usd']*number
                else:
                    self.baseline_purchase_costs['Quoted reactor package'] = audit['purchase_usd']*number
                    self.F_BM['Quoted reactor package'] = basis.installation_factor
                return
            costs = self._vessel_purchase_cost(d['Weight'],d['Diameter'],d['Length'])
            self.baseline_purchase_costs.update({k:v*d['Number of vessels'] for k,v in costs.items()})


class SaltSeparator(bst.Unit):
    """Specified ash/brine recovery with a published high-pressure solids budget."""
    _N_ins = 1
    _N_outs = 2
    def _init(self, config):
        self.config = config
        self.cost_audit = []

    def _run(self):
        feed, = self.ins
        clean, brine = self.outs
        clean.copy_like(feed)
        brine.empty()
        brine.T, brine.P = feed.T, feed.P
        brine.imass['Ash'] = float(feed.imass['Ash'])*self.config.salt_capture
        brine.imass['Water'] = min(float(feed.imass['Water']),
            float(brine.imass['Ash'])*self.config.brine_water_per_kg_ash)
        clean.mol[:] -= brine.mol

    def _design(self):
        self.design_results['Slurry kg/h'] = float(self.ins[0].F_mass)

    def _cost(self):
        basis = equipment_basis(self.config,'salt_separator',HYDROCYCLONE)
        audit = basis.estimate(self.design_results['Slurry kg/h'],bst.CE,self.config.cost_year)
        self.cost_audit = [audit]
        if audit['purchase_usd'] is None:
            self.installed_costs['High-pressure solids separator'] = audit['installed_usd']
        else:
            self.baseline_purchase_costs['High-pressure solids separator'] = audit['purchase_usd']
            self.F_BM['High-pressure solids separator'] = basis.installation_factor


class GasLiquidSeparator(bst.Flash):
    """Explicit empirical cooled-gas split with native flash-vessel sizing/cost.

    Water vapour uses pure-water Psat, not zero-water or a universal dry gas.
    NH3 and H2S partition are scenario inputs, not electrolyte equilibrium.
    """
    def _run(self):
        feed = self.ins[0]
        gas, liquid = self.outs
        liquid.copy_like(feed)
        gas.empty()
        gas.phase = 'g'
        gas.T = liquid.T = self.T
        gas.P = liquid.P = self.P
        gas.copy_flow(liquid,('CO2','CO','CH4','H2','N2','O2'),remove=True)
        # https://www.psi.ch/en/cpe/hydropilot : nitrogenous water/gas cleaning.
        # 1% NH3 and 90% H2S to gas are uncalibrated partition scenarios, not
        # electrolyte-equilibrium predictions. They must be replaced for design.
        for name, split in (('NH3',.01),('H2S',.90)):
            amount = float(liquid.imol[name])*split
            gas.imol[name] = amount
            liquid.imol[name] -= amount
        ywater = min(float(self.chemicals.Water.Psat(self.T))/self.P,.95)
        vapour = min(float(liquid.imol['Water']),float(gas.F_mol)*ywater/(1-ywater))
        gas.imol['Water'] = vapour
        liquid.imol['Water'] -= vapour

    def _design(self):
        self._size_flash_vessel()
        # Empirical isothermal separation may vaporize a small amount of water.
        # Native auxiliary exchanger supplies this phase-change duty once.
        duty = sum(sensible_h(s)+s.Hf for s in self.outs)-sensible_h(self.ins[0])-self.ins[0].Hf
        if abs(duty) > 1e-6:
            self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=self.ins,outs=self.outs,
                duty=duty,vle=False)


class HotCyclone(bst.Splitter):
    """Specified hot-char capture with the NREL vapor-cyclone cost anchor."""
    def _init(self, config, dry_feed):
        self.config = config
        self.dry_feed = dry_feed
        self.cost_audit = []
        super()._init(split={'BiocharMatrix':config.solids_capture,'Ash':config.solids_capture})

    def _run(self):
        super()._run()
        self.outs[0].phase = 's'
        self.outs[1].phase = 'g'

    def _design(self):
        self.design_results['Dry feed kg/h'] = float(self.dry_feed.F_mass-self.dry_feed.imass['Water'])

    def _cost(self):
        # https://digital.library.unt.edu/ark:/67531/metadc1015731/m2/1/high_res_d/993332.pdf
        # Appendix A Table A-1: one cyclone package at 2000 dry metric t/day.
        default = CostBasis(1262000.,2000000/24,
            'kg dry feed/h',2007,525.4,NREL,'Appendix A Table A-1, pyrolysis vapor cyclones',
            .6,installation_factor=3.02,scaling_basis='engineering six-tenths assumption',
            scope='hot vapor cyclone package; cost year qualified by ORNL/TM-2010/120 Appendix E')
        basis = equipment_basis(self.config,'cyclone',default)
        audit = basis.estimate(self.design_results['Dry feed kg/h'],bst.CE,self.config.cost_year)
        self.cost_audit = [audit]
        if audit['purchase_usd'] is None:
            self.installed_costs['Hot cyclone'] = audit['installed_usd']
        else:
            self.baseline_purchase_costs['Hot cyclone'] = audit['purchase_usd']
            self.F_BM['Hot cyclone'] = basis.installation_factor


class ProductSeparator(bst.Splitter):
    """Specified component recoveries; costs are separate from purity assumptions."""
    def _init(self, target, config, route):
        self.target = target
        self.config = config
        self.route = route
        split = {k:config.impurity_recovery for k in self.chemicals.IDs}
        split[target] = config.gas_recovery
        split['Water'] = 0.
        super()._init(split=split)
        self.cost_audit = []

    def _design(self):
        gas = self.ins[0]
        # 60 F/1 atm standard volume (not 0 C normal volume).
        self.design_results['H2 mmscfd'] = float(gas.imol['H2'])*23.690*24/0.028316846592/1e6
        self.design_results['Raw gas scfm'] = float(gas.F_mol)*23.690/0.028316846592/60

    def _cost(self):
        self.cost_audit = []
        if self.route == 'hydrogen':
            basis = equipment_basis(self.config,'gas_purification',HYDROGEN_PSA)
            size = self.design_results['H2 mmscfd']
        else:
            # https://www.epa.gov/system/files/documents/2024-01/pdh_full.pdf
            # Table 4-10: 2020 USD, treatment+compression 8300 USD/scfm at1000
            # and6200 at6000 scfm. Installed LFG-to-RNG package is a service proxy.
            from math import log
            exponent = 1+log(6200/8300)/log(6)
            default = CostBasis(
                8300000.,1000.,'scfm raw gas at 60 F, 1 atm',2020,596.2,
                'https://www.epa.gov/system/files/documents/2024-01/pdh_full.pdf',
                'Table 4-10, LFG-to-RNG treatment and compression',exponent,kind='installed',
                scope='bundled treatment and compression; no second compressor capital')
            basis = equipment_basis(self.config,'gas_purification',default)
            size = self.design_results['Raw gas scfm']
        if size <= 0:
            return
        audit = basis.estimate(size,bst.CE,self.config.cost_year)
        self.cost_audit = [audit]
        if audit['purchase_usd'] is None:
            self.installed_costs['Gas purification package'] = audit['installed_usd']
        else:
            self.baseline_purchase_costs['Gas purification package'] = audit['purchase_usd']
            self.F_BM['Gas purification package'] = basis.installation_factor


class BundledCompressor(bst.IsothermalCompressor):
    """Native compressor work; capital belongs to the RNG package once."""
    def _cost(self):
        pass


class MinimumReciprocatingCompressor(bst.IsothermalCompressor):
    """Native work, with an installed minimum-size reciprocating package budget."""
    def _cost(self):
        if not self.ins[0].F_mass:
            return
        # https://biosteam.readthedocs.io/en/latest/_modules/biosteam/units/compressor.html
        # Seider reciprocating correlation:100-20000 hp,5-7000 acfm.
        # Budget a minimum100 hp/5 acfm package below that range, rather than
        # extrapolating its purchase equation below the stated lower limit.
        # Actual operating power remains the native thermodynamic requirement;
        # part-load efficiency/turndown still requires a vendor performance curve.
        d = self.design_results
        alg = self.baseline_cost_algorithms['Reciprocating']
        number = d['Compressors in parallel']
        power = float(self.power_utility.get_property('consumption','hp'))/number
        rating = max(power,alg.hp_bounds[0])
        d['Cost-rated power per compressor hp'] = rating
        d['Minimum rated inlet acfm'] = max(float(self.ins[0].get_total_flow('cfm'))/number,alg.acfm_bounds[0])
        self.baseline_purchase_costs['Compressor(s)'] = number*bst.CE/alg.CE*alg.cost(rating)
        self.F_D['Compressor(s)'] = self.design_factors[d['Driver']]


class OffgasOxidizer(bst.ThermalOxidizer):
    """Native combustion; Seider fired-heater cost as an explicit service proxy."""
    _F_BM_default = {'Offgas furnace':2.19}

    def _run(self):
        if not self.ins[0].F_mass:
            self.ins[1].empty()
            self.ins[2].empty()
            self.outs[0].copy_like(self.ins[0])
        else:
            super()._run()

    def _design(self):
        if self.outs[0].F_mass:
            super()._design()

    def _cost(self):
        from biosteam.units.heat_exchange import compute_furnace_purchase_cost
        release = sum(float(s.LHV) for s in self.ins)
        if release > 0:
            self.baseline_purchase_costs['Offgas furnace'] = compute_furnace_purchase_cost(release*.947817,bst.CE)
