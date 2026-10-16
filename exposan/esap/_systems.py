#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    BlueS928

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from qsdsan import unit_operations as qsu
from biosteam.units import (
    DrumDryer,
    IsothermalCompressor,
    MolecularSieve,
    MultiEffectEvaporator,
    Stripper,
    )
from exposan.esap._components import (
    NH3_solution_molarity,
    ammonium_sulfate_solution_molarity,
    create_components,
    )
from exposan.esap._sanunits import ESAPBatchCrystallizer
from exposan.esap._process_settings import (
    acid_prices as default_acid_prices,
    electricity_price,
    natural_gas_price,
    natural_gas_price_per_kmol,
    )
from exposan.htl import _sanunits as htl_su

__all__ = ('create_system',)

_products = {
    'ammonium_sulfate': {
        'solute': 'NH42SO4',
        'crystal': 'NH42SO4Crystal',
        'feed': 'ammonium_sulfate_solution',
        'product': 'ammonium_sulfate_solids',
        },
    'anhydrous_ammonia': {
        'solute': 'NH3',
        'feed': 'ammonia_solution',
        'product': 'anhydrous_ammonia',
        },
    'ammonium_citrate': {
        'solute': 'NH4Citrate',
        'crystal': 'NH4CitrateCrystal',
        'feed': 'ammonium_citrate_solution',
        'product': 'ammonium_citrate_dibasic_solids',
        },
    'ammonium_lactate': {
        'solute': 'NH4Lactate',
        'feed': 'ammonium_lactate_solution',
        'product': 'ammonium_lactate_solution_60wt',
        },
    'ammonium_maleate': {
        'solute': 'NH4Maleate',
        'crystal': 'NH4MaleateCrystal',
        'feed': 'ammonium_maleate_solution',
        'product': 'ammonium_maleate_powder',
        },
    }

_product_aliases = {
    'as': 'ammonium_sulfate',
    'ammonium sulfate': 'ammonium_sulfate',
    'ammonium_sulfate_solids': 'ammonium_sulfate',
    'nh3': 'anhydrous_ammonia',
    'ammonia': 'anhydrous_ammonia',
    'anhydrous ammonia': 'anhydrous_ammonia',
    'ammonium citrate': 'ammonium_citrate',
    'ammonium_citrate_dibasic': 'ammonium_citrate',
    'ammonium_citrate_dibasic_solids': 'ammonium_citrate',
    'ammonium lactate': 'ammonium_lactate',
    'ammonium_lactate_60wt': 'ammonium_lactate',
    'ammonium_lactate_solution_60wt': 'ammonium_lactate',
    'ammonium maleate': 'ammonium_maleate',
    'ammonium_maleate_powder': 'ammonium_maleate',
    }

_default_feed_molarities = {
    'NH42SO4': ammonium_sulfate_solution_molarity,
    'NH3': NH3_solution_molarity,
    'NH4Citrate': ammonium_sulfate_solution_molarity,
    'NH4Lactate': 2 * ammonium_sulfate_solution_molarity,
    'NH4Maleate': ammonium_sulfate_solution_molarity,
    }

_acid_solution_molarity = 1
_acid_solution_density_kg_per_L = 1

_acid_components = {
    'ammonium_sulfate': 'H2SO4',
    'ammonium_citrate': 'CitricAcid',
    'ammonium_lactate': 'LacticAcid',
    'ammonium_maleate': 'MaleicAcid',
    }

_default_acid_molarities = {
    'ammonium_sulfate': _acid_solution_molarity,
    'ammonium_citrate': _acid_solution_molarity,
    'ammonium_lactate': _acid_solution_molarity,
    'ammonium_maleate': _acid_solution_molarity,
    }

_ammonium_per_product = {
    'ammonium_sulfate': 2,
    'ammonium_citrate': 2,
    'ammonium_lactate': 1,
    'ammonium_maleate': 2,
    }

_acid_per_ammonium = {
    'ammonium_sulfate': 0.5,
    'ammonium_citrate': 0.5,
    'ammonium_lactate': 1,
    'ammonium_maleate': 0.5,
    }

_evaporator_pressures = (101325, 73581, 50892, 32777, 20000)

_default_crystallizer_parameters = {
    'ammonium_sulfate': {
        'crystallization_recovery': 0.85,
        'crystal_capture_efficiency': 0.98,
        'wet_cake_moisture': 0.20,
        'heat_of_crystallization': 0,
        'separator_kWhr_per_m3': 1.40,
        },
    'ammonium_citrate': {
        'crystallization_recovery': 0.80,
        'crystal_capture_efficiency': 0.98,
        'wet_cake_moisture': 0.25,
        'heat_of_crystallization': 0,
        'separator_kWhr_per_m3': 1.40,
        },
    'ammonium_maleate': {
        'crystallization_recovery': 0.80,
        'crystal_capture_efficiency': 0.98,
        'wet_cake_moisture': 0.20,
        'heat_of_crystallization': 0,
        'separator_kWhr_per_m3': 1.40,
        },
    }


def _normalize_product(product):
    product = product.lower().strip().replace('-', '_')
    product = _product_aliases.get(product, product)
    if product not in _products:
        raise ValueError(
            '`product` must be one of '
            f'{tuple(_products)}, not {product!r}.'
            )
    return product


def _set_solution_flow(stream, solute_ID, molarity, flow_L_hr,
                       solution_density_kg_per_L, components):
    solute = getattr(components, solute_ID)
    solute_mass_kg_per_L = molarity * solute.MW / 1000
    water_mass_kg_per_L = solution_density_kg_per_L - solute_mass_kg_per_L
    if water_mass_kg_per_L < 0:
        raise ValueError(
            f'Solution density of {solution_density_kg_per_L} kg/L is too low '
            f'for {molarity} mol/L {solute_ID}.'
            )

    stream.empty()
    stream.phase = 'l'
    stream.imol[solute_ID] = molarity * flow_L_hr / 1000
    stream.imass['H2O'] = water_mass_kg_per_L * flow_L_hr


def _get_evaporation_fraction(solute_ID, molarity, flow_L_hr,
                              target_solute_wt_frac,
                              solution_density_kg_per_L, components):
    solute = getattr(components, solute_ID)
    solute_mass = molarity * flow_L_hr * solute.MW / 1000
    water_mass = solution_density_kg_per_L * flow_L_hr - solute_mass
    target_water_mass = solute_mass * (1 - target_solute_wt_frac) / target_solute_wt_frac
    water_removed = max(water_mass - target_water_mass, 0)
    water_removed_kmol = water_removed / components.H2O.MW
    total_kmol = molarity * flow_L_hr / 1000 + water_mass / components.H2O.MW
    return min(water_removed_kmol / total_kmol, 0.95) if total_kmol else 0


def _create_influent_streams(components, flow_L_hr, feed_molarities,
                             solution_density_kg_per_L):
    influents = {}
    for product_ID, data in _products.items():
        feed = qs.WasteStream(data['feed'], phase='l')
        influents[product_ID] = feed
        _set_solution_flow(
            feed,
            data['solute'],
            feed_molarities[data['solute']],
            flow_L_hr,
            solution_density_kg_per_L,
            components,
            )
    return influents


def _create_acid_stream(product, feed, solute_ID, components,
                        acid_molarities=None, acid_prices=None,
                        acid_solution_density_kg_per_L=(
                            _acid_solution_density_kg_per_L)):
    acid_molarities = dict(_default_acid_molarities, **(acid_molarities or {}))
    acid_prices = acid_prices or {}
    acid_ID = _acid_components[product]
    acid = qs.WasteStream(f'{acid_ID}_solution', phase='l')
    acid_molarity = acid_molarities[product]
    if acid_molarity <= 0:
        raise ValueError('Acid molarity must be greater than zero.')

    ammonium_flow = feed.imol[solute_ID] * _ammonium_per_product[product]
    acid_flow = ammonium_flow * _acid_per_ammonium[product]
    flow_L_hr = acid_flow * 1000 / acid_molarity
    _set_solution_flow(
        acid,
        acid_ID,
        acid_molarity,
        flow_L_hr,
        acid_solution_density_kg_per_L,
        components,
        )
    dry_acid_price = acid_prices.get(product, acid_prices.get(acid_ID, 0))
    if dry_acid_price < 0:
        raise ValueError('Dry-acid price cannot be negative.')
    if acid.F_mass:
        acid.price = dry_acid_price * acid.imass[acid_ID] / acid.F_mass
    return acid


def _add_acid_accounting_unit(product, acid, units):
    """Add a pass-through path that includes stoichiometric acid in TEA feeds."""
    AcidAccounting = htl_su.StreamTypeConverter(
        'AcidAccounting',
        ins=acid,
        outs=f'{_acid_components[product]}_consumed',
        init_with='WasteStream',
        )
    units.append(AcidAccounting)


def _create_solid_product_system(product, feed, solute_ID, crystal_ID, product_ID,
                                 evaporator_V, dryer_moisture_content,
                                 crystallizer_parameters, acid=None):
    units = []
    if acid:
        _add_acid_accounting_unit(product, acid, units)

    MEE = MultiEffectEvaporator(
        'MEE',
        ins=feed,
        outs=('concentrated_solution', 'evaporated_water'),
        P=_evaporator_pressures,
        V=evaporator_V,
        V_definition='Overall',
        )

    Crystallizer = ESAPBatchCrystallizer(
        'Crystallizer',
        ins=MEE-0,
        outs=('wet_crystal_cake', 'mother_liquor'),
        solute_ID=solute_ID,
        crystal_ID=crystal_ID,
        tau=4,
        N=2,
        **crystallizer_parameters,
        )

    Dryer = DrumDryer(
        'Dryer',
        ins=(Crystallizer-0, 'dry_gas', 'natural_gas'),
        outs=('dried_solids', 'hot_gas', 'dryer_missions'),
        split={solute_ID: 0, crystal_ID: 0, 'H2O': 1},
        moisture_content=dryer_moisture_content,
        moisture_ID='H2O',
        utility_agent='Natural gas',
        gas_composition=(('N2', 0.78), ('O2', 0.22)),
        )

    Product = htl_su.StreamTypeConverter(
        'ProductConverter',
        ins=Dryer-0,
        outs=product_ID,
        init_with='WasteStream',
        )

    return qs.System.from_units(
        f'esap_{product}',
        units=units + [MEE, Crystallizer, Dryer, Product],
        )


def _create_lactate_system(feed, evaporator_V, acid=None):
    units = []
    if acid:
        _add_acid_accounting_unit('ammonium_lactate', acid, units)

    MEE = MultiEffectEvaporator(
        'MEE',
        ins=feed,
        outs=('ammonium_lactate_solution_60wt_raw', 'evaporated_water'),
        P=_evaporator_pressures,
        V=evaporator_V,
        V_definition='Overall',
        )

    Product = htl_su.StreamTypeConverter(
        'ProductConverter',
        ins=MEE-0,
        outs='ammonium_lactate_solution_60wt',
        init_with='WasteStream',
        )

    return qs.System.from_units(
        'esap_ammonium_lactate',
        units=units + [MEE, Product],
        )


def _create_anhydrous_ammonia_system(feed, steam_kg_per_hr,
                                     target_pH, ammonia_pKa):
    PreStripper = qsu.Mixer(
        'PreStripper',
        ins=(feed, 'NaOH'),
        outs='NH3_solution',
        init_with='WasteStream',
        )

    @PreStripper.add_specification(run=True)
    def adjust_NaOH():
        influent, base = PreStripper.ins
        OH_excess_kmol_per_m3 = 10**(target_pH - 14)
        base.empty()
        base.phase = 's'
        base.imol['NaOH'] = OH_excess_kmol_per_m3 * influent.F_vol

    water_steam = qs.Stream('water_steam', H2O=steam_kg_per_hr, units='kg/hr',
                            phase='l', T=298.15)

    boiler = qsu.HXutility(
        'boiler',
        ins=water_steam,
        outs='steam',
        T=390,
        init_with='Stream',
        rigorous=True,
        )

    NH3Stripper = Stripper(
        'NH3Stripper',
        N_stages=2,
        ins=(PreStripper-0, boiler-0),
        outs=('vapor', 'liquid'),
        solute='NH3',
        )

    NH3MS1 = MolecularSieve(
        'NH3MS1',
        ins=NH3Stripper-0,
        outs=('NH3_rich_1', 'water_rich_1'),
        split=dict(Water=0.16, NH3=0.98),
        )

    NH3MS2 = MolecularSieve(
        'NH3MS2',
        ins=NH3MS1-0,
        outs=('NH3_rich_2', 'water_rich_2'),
        split=dict(Water=0.08, NH3=0.98),
        )

    NH3MS3 = MolecularSieve(
        'NH3MS3',
        ins=NH3MS2-0,
        outs=('NH3_rich_3', 'water_rich_3'),
        split=dict(Water=0, NH3=0.98),
        )

    NH3Compressor = IsothermalCompressor(
        'NH3Compressor',
        ins=NH3MS3-0,
        outs='anhydrous_ammonia_gas',
        P=2e6,
        eta=1,
        vle=True,
        )

    NH3Cooler = qsu.HXutility(
        'NH3Cooler',
        ins=NH3Compressor-0,
        outs='anhydrous_ammonia_cooled',
        T=298.15,
        init_with='Stream',
        rigorous=True,
        )

    Product = htl_su.StreamTypeConverter(
        'ProductConverter',
        ins=NH3Cooler-0,
        outs='anhydrous_ammonia',
        init_with='WasteStream',
        )

    return qs.System.from_units(
        'esap_anhydrous_ammonia',
        units=[PreStripper, boiler, NH3Stripper, NH3MS1, NH3MS2, NH3MS3,
               NH3Compressor, NH3Cooler, Product],
        )


def create_system(product='ammonium_sulfate', flowsheet=None,
                  flow_L_hr=1000, feed_molarities=None,
                  solution_density_kg_per_L=1,
                  acid_molarities=None, acid_prices=None,
                  solid_evaporator_solute_wt_frac=0.5,
                  lactate_product_solute_wt_frac=0.6,
                  dryer_moisture_content=0.01,
                  crystallizer_parameters=None,
                  target_pH=11.25,
                  ammonia_pKa=9.25,
                  steam_kg_per_hr=None,
                  simulate=True,
                  acid_solution_density_kg_per_L=(
                      _acid_solution_density_kg_per_L)):
    '''
    Create a downstream ESAP fertilizer production system.

    Sets shared BioSTEAM/QSDsan electricity and natural-gas prices from
    ``_process_settings``. Electricity is 0.1429 USD/kWh; natural gas is
    converted from 4.58 USD/thousand standard cubic feet using ideal methane
    at 60 deg F and 14.73 psia. Both the dryer fuel price (USD/kg) and the
    natural-gas heating-agent price (USD/kmol) are updated. These global
    settings also affect other systems in the same Python process.

    Parameters
    ----------
    product : str, optional
        Nitrogen fertilizer product to create. Options are:
        ``'ammonium_sulfate'``, ``'anhydrous_ammonia'``,
        ``'ammonium_citrate'``, ``'ammonium_lactate'``, and
        ``'ammonium_maleate'``. Common aliases such as ``'NH3'`` and
        ``'anhydrous ammonia'`` are also accepted.
    flowsheet : qsdsan.Flowsheet, optional
        Flowsheet for registering created streams and units. If not provided,
        a new flowsheet named ``'esap_{product}'`` is created.
    flow_L_hr : float, optional
        Volumetric flow rate of each influent solution, [L/hr]. This is used
        for all five product systems.
    feed_molarities : dict[str, float], optional
        Custom solute molarities for influent streams, [mol/L]. Keys should be
        component IDs, e.g., ``'NH42SO4'``, ``'NH3'``, ``'NH4Citrate'``,
        ``'NH4Lactate'``, and ``'NH4Maleate'``. Provided values override the
        defaults in ``_default_feed_molarities``.
    solution_density_kg_per_L : float, optional
        Assumed density of each influent solution, [kg/L]. This is used to
        convert molarity and volumetric flow into solute and water mass flows.
        Currently one density is applied to all influent solutions.
    acid_molarities : dict[str, float], optional
        Custom acid molarities for acid makeup streams, [mol/L]. Keys should be
        product IDs: ``'ammonium_sulfate'``, ``'ammonium_citrate'``,
        ``'ammonium_lactate'``, and ``'ammonium_maleate'``. All acids default
        to 1 mol/L. Acid consumption is calculated stoichiometrically from the
        ammonium fertilizer production rate. Used only for ammonium salt
        systems, not ``'anhydrous_ammonia'``.
    acid_prices : dict[str, float], optional
        Dry-acid prices, [USD/kg pure acid]. Keys may be product IDs or acid
        component IDs (e.g., ``'ammonium_sulfate'`` or ``'H2SO4'``). The
        effective price assigned to each aqueous acid stream charges only its
        pure-acid mass, not dilution water. If None, use the purity-adjusted
        wholesale quote midpoints in ``_process_settings.acid_prices``.
        An explicit dictionary replaces these defaults; missing entries then
        have zero price. Pass an empty dictionary for zero acid costs.
    acid_solution_density_kg_per_L : float, optional
        Assumed density of the acid makeup solution, [kg/L]. This is used to
        convert its volumetric flow into acid and water mass flows. The acid
        stream is an upstream consumption-accounting stream and is not mixed
        into the already-formed fertilizer solution.
    solid_evaporator_solute_wt_frac : float, optional
        Target solute mass fraction after multi-effect evaporation, [kg solute
        / kg solution]. Used only for solid-product systems:
        ``'ammonium_sulfate'``, ``'ammonium_citrate'``, and
        ``'ammonium_maleate'``.
    lactate_product_solute_wt_frac : float, optional
        Target ammonium lactate mass fraction after multi-effect evaporation,
        [kg ammonium lactate / kg solution]. Used only for
        ``'ammonium_lactate'``. The default is 0.6 for a 60 wt% solution.
    dryer_moisture_content : float, optional
        Target water mass fraction in dried product, [kg water / kg dried
        product]. Used only for solid-product systems:
        ``'ammonium_sulfate'``, ``'ammonium_citrate'``, and
        ``'ammonium_maleate'``.
    crystallizer_parameters : dict[str, float], optional
        Overrides for the selected solid product's black-box crystallizer and
        separator. Supported keys are ``'crystallization_recovery'``,
        ``'crystal_capture_efficiency'``, ``'wet_cake_moisture'``, and
        ``'heat_of_crystallization'``, and ``'separator_kWhr_per_m3'``.
        Baseline
        values are assumptions for preliminary TEA: sulfate uses 0.85, 0.98,
        and 0.20; citrate uses 0.80, 0.98, and 0.25; maleate uses 0.80, 0.98,
        and 0.20, respectively. Heat of crystallization defaults to zero kJ/kg
        for all salts. These values are not equilibrium predictions and should
        be replaced or varied when product-specific data are found. Separator
        electricity defaults to 1.40 kWh/m3 of crystallizer effluent, matching
        BioSTEAM's ``SolidsCentrifuge`` default.
    target_pH : float, optional
        Target pH before ammonia stripping. Used only for
        ``'anhydrous_ammonia'`` to estimate NaOH addition. 
        2 units higher than ammonia pKa.
    ammonia_pKa : float, optional
        Acid dissociation constant of ammonium/ammonia, [-]. Used only for
        ``'anhydrous_ammonia'`` to estimate the free-ammonia fraction at
        ``target_pH``.
    steam_kg_per_hr : float, optional
        Steam flow rate to the ammonia stripper, [kg/hr]. Used only for
        ``'anhydrous_ammonia'``. If not provided, this is set to 20% of the
        ammonia solution influent mass flow.
    simulate : bool, optional
        Whether to simulate the system immediately after creation.

    Returns
    -------
    sys : qsdsan.System
        Created downstream production system for the selected product.

    '''
    product = _normalize_product(product)
    qs.PowerUtility.price = electricity_price
    qs.settings.stream_prices['Natural gas'] = natural_gas_price
    qs.HeatUtility.get_agent('natural_gas').regeneration_price = natural_gas_price_per_kmol
    if acid_prices is None:
        acid_prices = default_acid_prices
    flowsheet_ID = f'esap_{product}'

    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
    flowsheet = flowsheet or qs.Flowsheet(flowsheet_ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)

    components = create_components(set_thermo=True)
    feed_molarities = dict(_default_feed_molarities, **(feed_molarities or {}))
    influents = _create_influent_streams(
        components,
        flow_L_hr,
        feed_molarities,
        solution_density_kg_per_L,
        )

    data = _products[product]
    solute_ID = data['solute']
    feed = influents[product]
    acid = None
    if product in _acid_components:
        acid = _create_acid_stream(
            product,
            feed,
            solute_ID,
            components,
            acid_molarities,
            acid_prices,
            acid_solution_density_kg_per_L,
            )

    if product == 'anhydrous_ammonia':
        if steam_kg_per_hr is None:
            steam_kg_per_hr = 0.2 * feed.F_mass
        sys = _create_anhydrous_ammonia_system(
            feed,
            steam_kg_per_hr,
            target_pH,
            ammonia_pKa,
            )
    elif product == 'ammonium_lactate':
        evaporator_V = _get_evaporation_fraction(
            solute_ID,
            feed_molarities[solute_ID],
            flow_L_hr,
            lactate_product_solute_wt_frac,
            solution_density_kg_per_L,
            components,
            )
        sys = _create_lactate_system(feed, evaporator_V, acid)
    else:
        crystallizer_parameters = dict(
            _default_crystallizer_parameters[product],
            **(crystallizer_parameters or {}),
            )
        evaporator_V = _get_evaporation_fraction(
            solute_ID,
            feed_molarities[solute_ID],
            flow_L_hr,
            solid_evaporator_solute_wt_frac,
            solution_density_kg_per_L,
            components,
            )
        sys = _create_solid_product_system(
            product,
            feed,
            solute_ID,
            data['crystal'],
            data['product'],
            evaporator_V,
            dryer_moisture_content,
            crystallizer_parameters,
            acid,
            )

    if simulate: sys.simulate()
    return sys
