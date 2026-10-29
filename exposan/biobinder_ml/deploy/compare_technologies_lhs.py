# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Paired HTL/AD/composting LHS simulations, checkpointed under deploy/results.

Run directly with --help. Composition draws are shared across technologies;
each recipe component keeps its wet blending weight and inherent moisture.
All native simulation objects are rebuilt per sample. No global feedstock data
are mutated, no HTL RF retraining is performed, and missing LCA factors are
reported as unavailable rather than silently interpreted as zero impacts.
"""
import argparse
import contextlib
import faulthandler
from datetime import datetime
import hashlib
import importlib.util
import inspect
from itertools import combinations
import json
from math import isfinite
import multiprocessing
import os
from pathlib import Path
import sys
import time
import traceback
from zoneinfo import ZoneInfo

sys.dont_write_bytecode = True
DEPLOY = Path(__file__).resolve().parent
RESULTS = DEPLOY/'results'
DRY_COMPONENTS = ('Carbohydrates','Proteins','Ash','Lipids')
TECHNOLOGIES = ('htl','ad','ad_rng','composting')
_MODULES = {}


class LiveLog:
    def __init__(self,stream):
        self.stream = stream

    def write(self,text):
        value = self.stream.write(text)
        if '\n' in text:
            self.stream.flush()
        return value

    def flush(self):
        self.stream.flush()


def bootstrap(extra_paths=()):
    repo = DEPLOY.parents[2]
    paths = [repo,DEPLOY,DEPLOY/'systems',*map(Path,extra_paths)]
    for name in ('biosteam','thermosteam','QSDsan'):
        candidate = repo.parent/name
        if candidate.is_dir():
            paths.append(candidate)
    for path in reversed(paths):
        if str(path) not in sys.path:
            sys.path.insert(0,str(path))


def load_file(path, name):
    path = Path(path).resolve()
    if not path.is_relative_to(DEPLOY):
        raise ValueError('Factory implementation must be inside deploy')
    if name not in _MODULES:
        spec = importlib.util.spec_from_file_location(name,path,
            submodule_search_locations=[str(path.parent)] if path.name=='__init__.py' else None)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        _MODULES[name] = module
    return _MODULES[name]


def recovery_tools():
    return load_file(DEPLOY/'technologies/__init__.py','_deployment_recovery_technologies')


def clean(value):
    """Strict JSON: nonfinite quantities are missing, never JSON NaN."""
    if isinstance(value,dict):
        return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [clean(v) for v in value]
    if isinstance(value,Path):
        return str(value)
    if hasattr(value,'item'):
        return clean(value.item())
    if isinstance(value,float) and not isfinite(value):
        return None
    return value


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(clean(value),indent=2,allow_nan=False),encoding='utf-8')
    temporary.replace(path)


def fingerprint(config):
    basis = {k:v for k,v in config.items() if k not in ('workers','python_paths')}
    return hashlib.sha256(json.dumps(basis,sort_keys=True).encode()).hexdigest()


def capacity_basis(config,composition):
    """Translate legacy dry short tons/day to operating kg/h and wet metric t/day."""
    water = composition['Water']
    if not 0<=water<1:
        raise ValueError('Capacity conversion requires water fraction in [0,1)')
    if 'capacity_dry_tpd' in config:
        if 'capacity_tpd' in config:
            raise ValueError('Specify capacity_dry_tpd or legacy wet capacity_tpd, not both')
        value = config['capacity_dry_tpd']
        uptime = config.get('capacity_uptime_ratio',.9)
        if not isfinite(value) or value<=0 or not isfinite(uptime) or not 0<uptime<=1:
            raise ValueError('Dry capacity must be positive; capacity_uptime_ratio must be in (0,1]')
        # Match exposan.saf._process_settings: short tons * 907.185 / (24 * uptime).
        dry_kg_h = value*907.185/(24*uptime)
        basis = 'legacy dry US short tons/day; hourly feed divided by capacity_uptime_ratio'
    else:
        value = config.get('capacity_tpd')
        if value is None or not isfinite(value) or value<=0:
            raise ValueError('Specify positive capacity_dry_tpd')
        dry_kg_h = value*1000/24*(1-water)
        uptime = None
        basis = 'backward-compatible wet metric tonnes per operating day'
    return dict(input_basis=basis,input_value=value,capacity_uptime_ratio=uptime,
        dry_kg_h=dry_kg_h,dry_metric_tpd_operating=dry_kg_h*24/1000,
        wet_metric_tpd_operating=dry_kg_h*24/1000/(1-water))


def validate(config):
    recipe = config['recipe_wet_mass']
    if (not recipe or any(not isfinite(v) or v < 0 for v in recipe.values()) or
        abs(sum(recipe.values())-1)>1e-9):
        raise ValueError('recipe_wet_mass must be nonnegative and sum to one')
    if not isinstance(config['n_samples'],int) or config['n_samples']<1:
        raise ValueError('n_samples must be a positive integer')
    if not 0 <= config['variation'] < 1:
        raise ValueError('variation must be in [0,1)')
    if config['variation_mode'] not in ('relative','absolute'):
        raise ValueError('variation_mode must be relative or absolute')
    capacity_basis(config,{'Water':0.})
    if not isinstance(config['workers'],int) or config['workers']<1:
        raise ValueError('workers must be a positive integer')
    timeout = config.get('max_technology_seconds',300)
    if not isfinite(timeout) or timeout<=0:
        raise ValueError('max_technology_seconds must be positive and finite')
    for key in ('worker_threads','worker_tasks_before_restart'):
        if not isinstance(config.get(key,1 if key=='worker_threads' else 20),int) or config.get(key,1 if key=='worker_threads' else 20)<1:
            raise ValueError(f'{key} must be a positive integer')
    interval = config.get('diagnostic_stack_interval_seconds',0)
    if not isfinite(interval) or interval<0:
        raise ValueError('diagnostic_stack_interval_seconds must be nonnegative')
    tools = recovery_tools()
    for feed in recipe:
        tools.Feedstock.from_name(feed)
    forbidden = {'feedstock','capacity_tpd','economics','characterization','config'}
    finance_keys = {'target_cepci','electricity_buy_usd_kwh','electricity_sell_usd_kwh',
                   'lifetime','start_year','discount_rate','operating_days','income_tax',
                   'depreciation','construction_schedule','working_capital_fraction',
                   'finance_interest','finance_years','finance_fraction','startup_months',
                   'startup_FOCfrac','startup_VOCfrac','startup_salesfrac'}
    for technology in TECHNOLOGIES[1:]:
        if forbidden.intersection(config[technology].get('system_kwargs',{})):
            raise ValueError(f'{technology}.system_kwargs overrides the paired comparison basis')
        if finance_keys.intersection(config[technology].get('economics',{})):
            raise ValueError('Use tea_overrides for shared finance instead of technology-specific finance')


def make_sample(config, sample_id, unit_draw):
    import numpy as np
    tools = recovery_tools()
    predictor = sys.modules[tools.__name__+'._cornell_ad'] if tools.__name__+'._cornell_ad' in sys.modules else None
    if predictor is None:
        predictor = load_file(DEPLOY/'technologies/_cornell_ad.py','_comparison_cornell_stoichiometry')
    perturbation = (2*np.asarray(unit_draw)-1)*config['variation']
    wet = dict.fromkeys((*DRY_COMPONENTS,'Water'),0.)
    component_inputs = {}
    potential = corrected = 0.
    provenance = []
    for name,weight in config['recipe_wet_mass'].items():
        if not weight:
            continue
        feed = tools.Feedstock.from_name(name)
        water = feed.composition['Water']
        baseline = np.array([feed.composition[k]/(1-water) for k in DRY_COMPONENTS])
        raw = (baseline*(1+perturbation) if config['variation_mode']=='relative' else
               np.maximum(baseline+perturbation,0.))
        if raw.sum()<=0:
            raise ValueError('Composition perturbation leaves no dry mass')
        dry = raw/raw.sum()
        composition = {k:float(v*(1-water)) for k,v in zip(DRY_COMPONENTS,dry)}
        composition['Water'] = water
        component_inputs[name] = composition
        for key in wet:
            wet[key] += weight*composition[key]
        prediction = predictor.predict(composition,1.,1.)
        # Cornell fD is methane-COD/feed-COD. The theoretical methane potential
        # is proportional to total stoichiometric COD (0.350 Nm3 CH4/kg COD).
        # Thus these are COD weights, not wet-mass or VS weights.
        cod_potential = weight*prediction['theoretical_methane_nm3_kg_wet']/.350
        potential += cod_potential
        corrected += cod_potential*feed.ad_biodegradability
        provenance.append(dict(feedstock=name,wet_weight=weight,
            theoretical_cod_kg_kg_blend=cod_potential,
            biodegradability=feed.ad_biodegradability,
            provenance=feed.ad_biodegradability_metadata))
    if potential<=0 or abs(sum(wet.values())-1)>1e-9:
        raise ValueError('Blend composition or COD potential does not close')
    sample = dict(sample_id=sample_id,composition=wet,
                  component_compositions=component_inputs,
                  ad_biodegradability=corrected/potential,
                  biodegradability_provenance=provenance,
                  lhs_quantiles=dict(zip(DRY_COMPONENTS,map(float,unit_draw))),
                  perturbations=dict(zip(DRY_COMPONENTS,map(float,perturbation))))
    sample['dry_composition'] = {k:wet[k]/(1-wet['Water']) for k in DRY_COMPONENTS}
    sample['capacity'] = capacity_basis(config,wet)
    return sample


def generate_samples(config):
    from scipy.stats import qmc
    draws = qmc.LatinHypercube(d=4,seed=config['seed']).random(config['n_samples'])
    return [make_sample(config,i+1,draw) for i,draw in enumerate(draws)]


def native_irr(tea):
    diagnostics = load_file(DEPLOY/'systems/ua_column_diagnostics.py','_comparison_ua_diagnostics')
    try:
        value,reason = diagnostics.solve_IRR_aitken_limit(tea)
        return clean(value*100),reason
    except Exception as error:
        return None,f'{type(error).__name__}: {error}'


def shared_basis(config):
    """Resolve the configured HTL finance/price basis even if its solver fails."""
    print('Common basis: loading recovery helpers',flush=True)
    tools = recovery_tools()
    relative,_ = config['htl']['factory'].rsplit(':',1)
    print('Common basis: importing HTL factory and simulation libraries',flush=True)
    module = load_file(DEPLOY/relative,'_comparison_htl_'+Path(relative).stem)
    print('Common basis: applying HTL process settings',flush=True)
    module._load_process_settings()
    import biosteam as bst
    e = tools.Economics.from_htl()
    fields = ('operating_days','lifetime','discount_rate','income_tax','start_year','depreciation',
              'construction_schedule','startup_months','startup_FOCfrac','startup_VOCfrac',
              'startup_salesfrac','working_capital_fraction','finance_interest','finance_years',
              'finance_fraction','insurance_fraction','labor_burden','maintenance_fraction')
    finance = {k:getattr(e,k) for k in fields}
    electricity_price = config['htl'].get('system_kwargs',{}).get('electricity_price')
    if electricity_price is None:
        electricity_price = bst.PowerUtility.price
    finance.update(target_cepci=float(bst.CE),electricity_buy_usd_kwh=electricity_price,
                   electricity_sell_usd_kwh=electricity_price)
    overrides = dict(config['htl'].get('system_kwargs',{}).get('tea_overrides',{}))
    overrides.update(config.get('tea_overrides',{}))
    aliases = dict(IRR='discount_rate',WC_over_FCI='working_capital_fraction',
                   property_insurance='insurance_fraction',maintenance='maintenance_fraction')
    for key,value in overrides.items():
        if key=='duration':
            finance.update(start_year=value[0],lifetime=value[1]-value[0])
        elif aliases.get(key,key) in finance:
            finance[aliases.get(key,key)] = value
    price = config.get('feedstock_price_usd_t')
    if price is None:
        price = sum(weight*module.feedprice_dct[name]*1000 for name,weight in config['recipe_wet_mass'].items())
    tools.Economics(**finance,feedstock_price_usd_t=price)
    return finance,price


def htl_equipment_audit(system):
    """Snapshot unit breakdowns; unavailable reporting fields cannot fail a case."""
    import biosteam as bst
    osbl = system.TEA.OSBL_units
    hours = system.operating_hours
    records = []
    for unit in system.units:
        errors = []
        def read(attribute, convert):
            try:
                value = getattr(unit,attribute)
                return convert(value) if value is not None else None
            except Exception as error:
                errors.append(f'{attribute}: {type(error).__name__}: {error}')
                return None
        installed = read('installed_cost',float)
        power = read('power_utility',lambda p:float(p.rate))
        utility = read('utility_cost',lambda v:float(v)*hours)
        records.append(dict(
            name=unit.ID,unit_id=unit.ID,unit_class=type(unit).__name__,
            cost_area='OSBL' if unit in osbl else 'ISBL',
            included_in_system_capital=unit in system.cost_units,
            installed_usd=installed,native_installed_usd=installed,
            purchase_cost_usd=read('purchase_cost',float),
            baseline_purchase_costs_usd=read('baseline_purchase_costs',dict),
            purchase_costs_usd=read('purchase_costs',dict),
            installed_costs_usd=read('installed_costs',dict),
            design_results=read('design_results',dict),design_units=read('_units',dict),
            annual_utility_cost_usd=utility,
            utility_cost_status='available' if utility is not None else 'unavailable',
            net_electricity_kw=power,
            annual_electricity_cost_usd=power*bst.PowerUtility.price*hours if power is not None else None,
            target_cepci=float(bst.CE),reporting_errors=errors,
            source='Native HTL unit costs; ISBL/OSBL split matches comparison_reporting.py',
            cost_basis='Native unit capital; excludes TEA indirect costs and land'))
    total = sum(r['native_installed_usd'] or 0. for r in records if r['included_in_system_capital'])
    expected = float(system.installed_equipment_cost)
    reconciled = abs(total-expected)<=max(1e-6,abs(expected)*1e-8)
    for record in records:
        record.update(system_installed_usd=expected,
            audit_installed_total_usd=total,capital_reconciliation_status='OK' if reconciled else 'MISMATCH')
    return records


def htl_equipment_snapshot(system):
    """Return reporting status independently of simulation/TEA acceptance."""
    try:
        rows = htl_equipment_audit(system)
        status = 'PARTIAL' if any(r['reporting_errors'] or r['capital_reconciliation_status']!='OK' for r in rows) else 'OK'
        return dict(equipment_audit=rows,equipment_audit_status=status)
    except Exception as error:
        return dict(equipment_audit=[],equipment_audit_status='FAILED',
                    equipment_audit_error=f'{type(error).__name__}: {error}')


def htl_case(config,sample):
    timing_started = time.monotonic()
    import qsdsan as qs
    relative, function = config['htl']['factory'].rsplit(':',1)
    module = load_file(DEPLOY/relative,'_comparison_htl_'+Path(relative).stem)
    factory = getattr(module,function)
    if 'feedstock_composition_override' not in inspect.signature(factory).parameters:
        raise TypeError('HTL factory must explicitly accept feedstock_composition_override')
    kwargs = dict(config['htl'].get('system_kwargs',{}))
    for reserved in ('feedstock_composition_override','central_dry_flowrate','HTL_yields_override'):
        if reserved in kwargs:
            raise ValueError(f'{reserved} would override the sampled feed or RF; configure it in the paired basis instead')
    wet = sample['composition']
    capacity = capacity_basis(config,wet)
    kwargs.update(feedstock_composition_override=wet,
        feedstock_id=next(iter(config['recipe_wet_mass'])),
        central_dry_flowrate=capacity['dry_kg_h'])
    overrides = dict(kwargs.pop('tea_overrides',{}))
    overrides.update(config.get('tea_overrides',{}))
    kwargs['tea_overrides'] = overrides
    # Execution scheduling only: preserve the trained RF and its predictions.
    # A single model thread also avoids nested multiprocessing inside LHS jobs.
    if hasattr(module,'_get_rf_model'):
        model = module._get_rf_model()
        if hasattr(model,'n_jobs'):
            model.n_jobs = config['htl'].get('rf_n_jobs',1)
    print('HTL: building sampled native system',flush=True)
    setup_seconds = time.monotonic()-timing_started
    build_started = time.monotonic()
    system = factory(**kwargs)
    build_seconds = time.monotonic()-build_started
    tea = system.TEA
    system.operating_hours = tea.operating_days*24
    prices = module.feedprice_dct
    feed_price = config.get('feedstock_price_usd_t')
    if feed_price is None:
        feed_price = sum(weight*prices[name]*1000 for name,weight in config['recipe_wet_mass'].items())
    system.flowsheet.stream.scaled_feedstock.price = feed_price/1000
    for name,price in config['htl'].get('stream_prices_usd_kg',{}).items():
        getattr(system.flowsheet.stream,name).price = price
    lca = system.LCA
    lca.lifetime = tea.duration[1]-tea.duration[0]
    lca.uptime_ratio = system.operating_hours/(365*24)
    # Replace the inherited avoided-landfill credit with the shared scenario.
    item = qs.ImpactItem.get_item('feedstock_item')
    if item is None:
        raise ValueError('HTL has no feedstock_item; cannot align the avoided-disposal boundary')
    for indicator in lca.indicators:
        if indicator.ID in item.CFs:
            item.CFs[indicator.ID] = -config['lca']['avoided_disposal_kg_co2e_t']/1000
    print('HTL: simulating sampled native system',flush=True)
    simulate_started = time.monotonic()
    system.simulate(**config['htl'].get('simulate_kwargs',{}))
    simulate_seconds = time.monotonic()-simulate_started
    print('HTL: collecting native TEA/LCA',flush=True)
    feed = system.flowsheet.stream.scaled_feedstock
    realized_tpd = feed.F_mass*24/1000
    requested_tpd = capacity['wet_metric_tpd_operating']
    if abs(realized_tpd-requested_tpd)>max(1e-6,requested_tpd*1e-6):
        raise ValueError(f'HTL realized {realized_tpd} wet t/day; requested {requested_tpd}')
    for key,value in wet.items():
        if abs(float(feed.imass[key])/feed.F_mass-value)>1e-7:
            raise ValueError(f'HTL realized a different sampled {key} fraction')
    annual_tonnes = feed.F_mass*system.operating_hours/1000
    irr,reason = native_irr(tea)
    annual_impacts = lca.get_total_impacts(annual=True)
    gwp = annual_impacts.get('GWP')
    reporting = load_file(DEPLOY/'systems/comparison_reporting.py','_lhs_comparison_reporting')
    snapshot = reporting.collect(system,'HTL',{
        'Native_cut_fractions':kwargs.get('cutoff_fracs'),
        'Intermediate_route':kwargs.get('heavy_route')})
    metrics = dict(NPV_usd=float(tea.NPV),IRR_percent=irr,IRR_status=reason,
        FCI_usd=float(tea.FCI),TCI_usd=float(tea.TCI),
        annual_feed_tonnes=annual_tonnes,realized_capacity_tpd=realized_tpd,capacity_basis=capacity,
        GWP_kgCO2e_year=gwp,GWP_kgCO2e_t_wet=gwp/annual_tonnes if gwp is not None else None,
        LCA_status='native HTL boundary; review against alternatives',
        annual_impacts=annual_impacts,
        operating_hours=system.operating_hours,
        net_electricity_kw=float(system.power_utility.rate),
        product_kg_h={s.ID:float(s.F_mass) for s in system.products},
        product_prices_usd_kg={s.ID:float(s.price) for s in system.products},
        HTL_yields=dict(system.flowsheet.unit.HTL.dw_yields))
    metrics.update(htl_equipment_snapshot(system))
    metrics.update(snapshot['metrics'])
    # Preserve reporter titles, row order, duplicate items, and native table columns.
    metrics['reporting_blocks'] = {
        title: json.loads(table.to_json(orient='split',date_format='iso'))
        for title,table in snapshot['blocks']}
    metrics['distillation_columns'] = snapshot['columns']
    metrics['distillation_trials'] = snapshot['trials']
    metrics.update(cashflow_usd=tea.cashflow_array.tolist(),
                   discount_periods=tea._get_duration_array().tolist(),
                   execution_timing_seconds=dict(RF_and_setup=setup_seconds,
                       system_build=build_seconds,system_simulate=simulate_seconds))
    tools = recovery_tools()
    economics = tools.Economics.from_htl(tea,feedstock_price_usd_t=feed_price)
    finance = {k:getattr(economics,k) for k in (
        'operating_days','lifetime','discount_rate','income_tax','start_year','depreciation',
        'construction_schedule','startup_months','startup_FOCfrac','startup_VOCfrac',
        'startup_salesfrac','working_capital_fraction','finance_interest','finance_years',
        'finance_fraction','insurance_fraction','labor_burden','maintenance_fraction')}
    metrics['financial_basis'] = finance
    import biosteam as bst
    finance.update(target_cepci=float(bst.CE),electricity_buy_usd_kwh=bst.PowerUtility.price,
                   electricity_sell_usd_kwh=bst.PowerUtility.price)
    metrics['feedstock_price_usd_t'] = feed_price
    if hasattr(module,'_INTEGRATION_ASSUMPTIONS'):
        module._INTEGRATION_ASSUMPTIONS.pop(system,None)
    return metrics,finance,feed_price


def alternative_lca_factors(technology,config,inventory):
    """Reuse matching HTL utility factors; retain explicit first-run exclusions."""
    settings = config['lca']
    factors = {}
    sources = {}
    if settings.get('reuse_htl_utilities',False):
        from exposan.saf import gwp_dct
        overrides = config['htl'].get('system_kwargs',{})
        electricity = overrides.get('electricity_GHG')
        if electricity is None:
            electricity = gwp_dct['electricity']
        for key in ('electricity_import_kwh','electricity_export_kwh'):
            factors[key] = electricity
            sources[key] = 'HTL electricity factor; same grid mix for import/displacement'
        for key in ('heat_import_mj','heat_export_mj'):
            factors[key] = gwp_dct['steam']
            sources[key] = 'HTL steam factor; exported heat assumes displacement of the same steam'
    configured = dict(settings.get('factors',{}).get('ad',{})) if technology == 'ad_rng' else {}
    configured.update(settings.get('factors',{}).get(technology,{}))
    for key,value in configured.items():
        if value is not None:
            factors[key] = value
            sources[key] = 'Explicit configured factor'
    missing = [k for k,v in inventory.items() if v and k not in factors]
    exclusions = dict(settings.get('excluded_activities',{}).get('ad',{})) if technology == 'ad_rng' else {}
    exclusions.update(settings.get('excluded_activities',{}).get(technology,{}))
    excluded = {k:exclusions[k] for k in missing if k in exclusions}
    unresolved = [k for k in missing if k not in excluded]
    audit = [dict(activity=k,annual_quantity=v,factor_kgCO2e_per_unit=factors.get(k),
        status='CHARACTERIZED' if k in factors else 'EXCLUDED' if k in excluded else 'MISSING' if v else 'INACTIVE',
        source_or_reason=sources.get(k,excluded.get(k,'')),
        annual_kgCO2e=v*factors[k] if k in factors else None) for k,v in inventory.items()]
    return factors,missing,excluded,unresolved,audit


def alternative_case(technology,config,sample,finance,feed_price):
    tools = recovery_tools()
    settings = config[technology]
    wet = sample['composition']
    capacity = capacity_basis(config,wet)
    wet_tpd = capacity['wet_metric_tpd_operating']
    feed = tools.Feedstock(name='+'.join(f'{k}:{v:g}' for k,v in config['recipe_wet_mass'].items()),
        total_solids=1-wet['Water'],vs_fraction=1-wet['Ash']/(1-wet['Water']),
        composition=wet,ad_biodegradability=sample['ad_biodegradability'],
        ad_biodegradability_metadata=dict(status='COD-weighted literature mixture',
            components=sample['biodegradability_provenance']),
        product_price_usd_t=settings.get('product_price_usd_t',0.),
        marketable_fraction=settings.get('marketable_fraction',1.))
    economic_kwargs = dict(finance)
    economic_kwargs.update(settings.get('economics',{}))
    economic_kwargs['feedstock_price_usd_t'] = feed_price
    economics = tools.Economics(**economic_kwargs)
    cfg = dict(settings.get('config',{}))
    kwargs = dict(settings.get('system_kwargs',{}))
    if technology in ('ad','ad_rng'):
        cfg['route'] = 'rng' if technology == 'ad_rng' else 'chp'
        result = tools.create_ad_system(wet_tpd,feedstock=feed,
                    config=tools.ADConfig(**cfg),economics=economics,**kwargs)
    else:
        if isinstance(cfg.get('pactea'),dict):
            cfg['pactea'] = tools.PaCTEAConfig(**cfg['pactea'])
        bulking = settings.get('bulking_composition')
        if bulking is None and settings.get('bulking_feedstock'):
            bulking = tools.Feedstock.from_name(settings['bulking_feedstock']).composition
        result = tools.create_composting_system(wet_tpd,feedstock=feed,
                    config=tools.CompostConfig(**cfg),economics=economics,
                    bulking_composition=bulking,**kwargs)
    factors,missing,excluded,unresolved,audit = alternative_lca_factors(technology,config,result.inventory_annual)
    # Missing/excluded activities use zero only inside the native adapter.
    # Export a partial total only when every active omission is explicitly excluded.
    diagnostic_factors = {k:factors.get(k,0.) for k in result.inventory_annual}
    native,lca,streams = result.to_qsdsan(diagnostic_factors,
        avoided_disposal_kg_co2e_t=config['lca']['avoided_disposal_kg_co2e_t'])
    native.simulate()
    analysis,tea = result.tea(native_system=native,return_tea=True)
    equipment_audit = [item for unit in native.units for item in getattr(unit,'cost_basis_audit',[])]
    irr,reason = native_irr(tea)
    metrics = dict(NPV_usd=analysis['NPV'],IRR_percent=irr,IRR_status=reason,
        FCI_usd=analysis['fixed_capital'],TCI_usd=analysis['total_investment'],
        annual_opex_usd=analysis['annual_opex'],annual_revenue_usd=analysis['annual_revenue'],
        EBITDA_usd=analysis['EBITDA'],annual_feed_tonnes=result.annual_tonnes,
        realized_capacity_tpd=result.capacity_tpd,operating_hours=result.hours,capacity_basis=capacity,
        land_area_acres=result.land_acres,land_area_ha=result.land_acres*0.40468564224,
        net_electricity_kw=float(native.power_utility.rate),
        product_kg_h={s.ID:float(s.F_mass) for s in native.products},
        annual_inventory=result.inventory_annual,design=result.design,
        equipment_audit=equipment_audit,
        native_conventions=dict(unit_classes=[type(u).__mro__[1].__name__ for u in native.units],
            installed_equipment_cost_usd=float(native.installed_equipment_cost),
            target_cepci=economics.target_cepci,purchase_cost_split_available=False,
            mass_balance_error_kg_h=sum(float(u.mass_balance_error()) for u in native.units),
            net_heat_utility_kj_h=sum(float(hu.duty) for u in native.units for hu in u.heat_utilities),
            process_basis='composition-dependent aggregate correlation; rebuild for changed feed'),
        financial_basis=analysis['financial_assumptions'],feedstock_price_usd_t=feed_price,
        annual_costs=analysis['annual_costs'],annual_revenues=analysis['annual_revenues'],
        LCA_missing_factors=missing,LCA_excluded_activities=excluded,LCA_factor_audit=audit,notes=result.notes,
        cashflow_usd=analysis['cashflows'],discount_periods=analysis['discount_periods'],
        GWP_kgCO2e_year=None,GWP_kgCO2e_t_wet=None)
    if unresolved or (excluded and config['lca'].get('require_complete_factors')):
        metrics['LCA_status'] = 'MISSING_FACTORS'
        if config['lca'].get('require_complete_factors'):
            raise ValueError(f'{technology} LCA missing factors: {missing}')
    else:
        impacts = lca.get_total_impacts(annual=True)
        annual = float(next(iter(impacts.values())))
        expected = result.lca(diagnostic_factors,avoided_disposal_kg_co2e_t=config['lca']['avoided_disposal_kg_co2e_t'])['annual_kg_co2e']
        if abs(annual-expected)>max(1e-6,abs(expected)*1e-8):
            raise ValueError('Native QSDsan annual inventory does not reconcile')
        metrics.update(LCA_status='PARTIAL_EXCLUSIONS; utility-only screening' if excluded else 'OK; explicit alternative inventory factors',
                       GWP_kgCO2e_year=annual,GWP_kgCO2e_t_wet=annual/result.annual_tonnes)
    if abs(float(tea.NPV)-analysis['NPV'])>1e-6:
        raise ValueError('Native TEA NPV does not reconcile')
    return metrics


def _run_technology(config,sample,run_dir,technology):
    """One isolated worker job; exceptions remain explicit per technology."""
    bootstrap(config.get('python_paths',()))
    rows = []
    finance = feed_price = None
    for technology in (technology,):
        started = time.monotonic()
        log = Path(run_dir)/'logs'/f'{sample["sample_id"]:05d}_{technology}.log'
        row = dict(sample_id=sample['sample_id'],technology=technology,
                   status='FAILED',log=str(log),metrics={},error='')
        with log.open('w',encoding='utf-8') as stream, (Path(run_dir)/'logs'/f'{sample["sample_id"]:05d}_{technology}.stack.log').open('w') as stacks,contextlib.redirect_stdout(LiveLog(stream)),contextlib.redirect_stderr(LiveLog(stream)):
            faulthandler.enable(file=stacks,all_threads=True)
            interval = config.get('diagnostic_stack_interval_seconds',0)
            if interval:
                faulthandler.dump_traceback_later(interval,repeat=True,file=stacks)
            try:
                print(f'Worker PID={os.getpid()} Python={sys.executable} NUMBA_DISABLE_JIT={os.environ.get("NUMBA_DISABLE_JIT","unset")}',flush=True)
                basis_started = time.monotonic()
                finance,feed_price = shared_basis(config)
                basis_seconds = time.monotonic()-basis_started
                if technology=='htl':
                    metrics,finance,feed_price = htl_case(config,sample)
                else:
                    if finance is None:
                        raise RuntimeError('Common financial/feed-price basis unavailable because HTL failed; alternative not run with unmatched assumptions')
                    metrics = alternative_case(technology,config,sample,finance,feed_price)
                metrics.setdefault('execution_timing_seconds',{})['common_basis_setup'] = basis_seconds
                row.update(status='OK',metrics=clean(metrics),worker_pid=os.getpid())
            except Exception:
                row['error'] = traceback.format_exc()
                print(row['error'])
            finally:
                faulthandler.cancel_dump_traceback_later()
                faulthandler.disable()
                # Workers own these registries. Snapshot scalar results before
                # releasing objects; no live UI/session system is modified.
                if technology=='htl':
                    for module in _MODULES.values():
                        if hasattr(module,'_INTEGRATION_ASSUMPTIONS'):
                            module._INTEGRATION_ASSUMPTIONS.clear()
        row['elapsed_seconds'] = time.monotonic()-started
        rows.append(row)
    return rows


def _persistent_worker(job_path,config,run_dir):
    """Keep imports/RF warm while rebuilding each native system from its sample."""
    threads = str(config.get('worker_threads',1))
    for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name] = threads
    bootstrap(config.get('python_paths',()))
    job_path = Path(job_path)
    try:
        while True:
            if not job_path.exists():
                time.sleep(.05)
                continue
            task = json.loads(job_path.read_text())
            job_path.unlink()
            if task is None:
                return
            sample,technology = task
            row = _run_technology(config,sample,run_dir,technology)[0]
            result_path = Path(run_dir)/'logs'/f'{sample["sample_id"]:05d}_{technology}.result.json'
            write_json(result_path,row)
    finally:
        job_path.unlink(missing_ok=True)


def _technology_process(config,sample,run_dir,technology,result_path):
    """Compatibility entry point for already-running isolated-worker batches."""
    write_json(result_path,_run_technology(config,sample,run_dir,technology)[0])


def run_sample(config,sample,run_dir):
    """Run a preflight technology set through the same supervised worker scheduler."""
    iterator = completed_samples(config,[sample],run_dir)
    try:
        return next(iterator)[1]
    finally:
        iterator.close()


def comparison_gaps(row):
    """Open evidence gaps; native API compliance alone does not validate a model."""
    tech,m = row['technology'],row.get('metrics',{})
    gaps = []
    def add(area,detail):
        gaps.append(dict(sample_id=row['sample_id'],technology=tech,area=area,status='OPEN',detail=detail))
    if row['status']!='OK':
        add('Execution','No accepted native simulation/TEA/LCA result: '+row['status'])
        return gaps
    if tech in ('ad','ad_rng','composting'):
        add('Process resolution','Aggregate composition model; detailed reactor, separator and utility-network validation remains outstanding.')
        add('Capital scope','Bundled installed-plant allocations lack equipment FOB quotations and installation/BM splits; reconcile OSBL, civil, engineering, startup and contingency with HTL.')
        add('Utility design','Electricity/heat/diesel intensities are screening assumptions; confirm pump, blower, mixer, separation and CHP performance and heat quality/offtake.')
        add('Mass/energy validation','Total mass closure checked with pseudo-components; full elemental and enthalpy balances are not established.')
        add('Operating costs','Confirm labor, maintenance, disposal, utility contracts and replacement schedules on a common cost-year and scope basis.')
        add('Product market','Confirm marketable fraction, price, permitted outlet and product specification; zero default sale price does not establish quality.')
        if tech in ('ad','ad_rng'):
            add('AD conversion','Literature BMP corrections require full-scale reactor/HRT, inhibition and co-digestion validation.')
            add('AD capital basis',m.get('design',{}).get('capital_basis',{}).get('status','Unverified AD cost basis'))
            add('AD effluent','Confirm digestate/centrate quality, treatment equipment and disposal/application boundary and costs.')
        else:
            add('Compost sizing','Capital scales with total throughput; active/curing volume and area do not yet independently size/cost the pad, aeration, biofilter and screen.')
            add('Compost quality','Resolve carbohydrate pools and confirm moisture, maturity, pathogens, contaminants and nutrient availability for this sludge-rich blend.')
        for item in m.get('equipment_audit',[]):
            ratio = item['size']/item.get('reference_size',item['size'])
            if not .1<=ratio<=10:
                add('Cost extrapolation',f'{item["name"]}: {ratio:g} times reference size; obtain an applicable cost correlation.')
    else:
        add('HTL scenario','Review RF applicability, cut fractions, surrogate composition, product outlets and convergence for this blend.')
    for activity in m.get('LCA_missing_factors',[]):
        add('LCA factor',activity+': '+m.get('LCA_excluded_activities',{}).get(activity,'Factor unavailable'))
    add('Lifecycle boundary','Reconcile transport, construction, storage, product use, nutrient credits and soil carbon across all technologies.')
    return gaps


def export_system_workbook(path,summary,rows):
    """Keep native reporting tables on separate sheets, stacked by sample."""
    import re
    import pandas as pd
    from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE
    tables = {}
    for row in rows:
        identity = {key:row[key] for key in ('sample_id','technology','status')}
        metrics = row.get('metrics',{})
        for title,table in metrics.get('reporting_blocks',{}).items():
            if title == 'TEA TOTALS':
                continue  # These scalar totals already appear on System results.
            frame = pd.DataFrame(table['data'],columns=table['columns'])
            for key,value in reversed(list(identity.items())):
                frame.insert(0,key,value)
            tables.setdefault(title,[]).append(frame)
        for key,title in (('distillation_columns','Distillation'),
                          ('distillation_trials','Distillation trials')):
            records = metrics.get(key,[])
            if records:
                tables.setdefault(title,[]).append(pd.DataFrame([
                    dict(record,**identity) for record in records]))
    def excel_values(frame):
        def excel_value(value):
            if isinstance(value,(dict,list,tuple)):
                value = json.dumps(clean(value),allow_nan=False,default=str)
            return ILLEGAL_CHARACTERS_RE.sub('',value) if isinstance(value,str) else value
        return frame.apply(lambda col: col.map(excel_value))
    with pd.ExcelWriter(path,engine='openpyxl') as writer:
        excel_values(summary).to_excel(writer,sheet_name='System results',index=False)
        used = {'system results'}
        for title,frames in tables.items():
            name = re.sub(r'[\\/*?:\[\]]','_',title)[:31]
            base = name
            counter = 2
            while name.lower() in used:
                suffix = f' {counter}'
                name = base[:31-len(suffix)]+suffix
                counter += 1
            used.add(name.lower())
            excel_values(pd.concat(frames,ignore_index=True)).to_excel(
                writer,sheet_name=name,index=False)
        for sheet in writer.book:
            sheet.freeze_panes = 'D2'
            sheet.auto_filter.ref = sheet.dimensions
            for column in sheet.columns:
                sheet.column_dimensions[column[0].column_letter].width = min(48,max(16,
                    max(len(str(cell.value or '')) for cell in column[:50])+2))


def export_results(run_dir,rows,config):
    import numpy as np
    import pandas as pd
    # Summary rows contain scalar metrics only; detailed records stay in JSONL
    # and dedicated workbook sheets rather than expanding the CSV horizontally.
    summary_rows = [dict(row,metrics={key:value for key,value in row.get('metrics',{}).items()
        if not isinstance(value,(dict,list,tuple))}) for row in rows]
    flat = pd.json_normalize(summary_rows,sep='.')
    flat.to_csv(run_dir/'system_results.csv',index=False)
    export_system_workbook(run_dir/'system_results.xlsx',flat,rows)
    metrics = ('NPV_usd','IRR_percent','FCI_usd','annual_opex_usd','annual_revenue_usd','GWP_kgCO2e_t_wet')
    paired = pd.DataFrame(index=sorted(set(row['sample_id'] for row in rows)))
    paired.index.name = 'sample_id'
    gap_rows = [gap for row in rows for gap in comparison_gaps(row)]
    pd.DataFrame(gap_rows).to_csv(run_dir/'comparison_gaps.csv',index=False)
    equipment_rows = [dict(sample_id=row['sample_id'],technology=row['technology'],**item)
        for row in rows for item in row.get('metrics',{}).get('equipment_audit',[])]
    pd.DataFrame(equipment_rows).to_csv(run_dir/'equipment_cost_audit.csv',index=False)
    pd.DataFrame([dict(sample_id=row['sample_id'],technology=row['technology'],
        status=row['status'],**row.get('metrics',{}).get('native_conventions',{})) for row in rows]
        ).to_csv(run_dir/'native_conventions.csv',index=False)
    paired['comparison_status'] = 'NOT_READY: inspect comparison_gaps.csv'
    summary = []
    for technology in TECHNOLOGIES:
        selected = [row for row in rows if row['technology']==technology]
        summary.append(dict(technology=technology,attempted=len(selected),
                            succeeded=sum(row['status']=='OK' for row in selected),
                            failed=sum(row['status']!='OK' for row in selected),metric='status'))
        for metric in metrics:
            values = {row['sample_id']:row['metrics'].get(metric) if row['status']=='OK' else None for row in selected}
            paired[f'{technology}.{metric}'] = pd.Series(values,dtype=float)
            finite = np.array([v for v in values.values() if v is not None and np.isfinite(v)])
            if finite.size:
                summary.append(dict(technology=technology,metric=metric,count=len(finite),
                    mean=float(finite.mean()),p05=float(np.quantile(finite,.05)),
                    median=float(np.median(finite)),p95=float(np.quantile(finite,.95))))
    for first,second in combinations(TECHNOLOGIES,2):
        for metric in ('NPV_usd','GWP_kgCO2e_t_wet'):
            if gap_rows and not config.get('comparison',{}).get('allow_screening_differences',False):
                continue
            if metric.startswith('GWP') and not config['lca'].get('boundaries_confirmed_comparable'):
                continue
            paired[f'{first}_minus_{second}.{metric}'] = paired[f'{first}.{metric}']-paired[f'{second}.{metric}']
    paired.to_csv(run_dir/'paired_comparison.csv')
    pd.DataFrame(summary).to_csv(run_dir/'summary.csv',index=False)
    inventories = []
    factor_audit = []
    templates = {key:{} for key in TECHNOLOGIES[1:]}
    for row in rows:
        metrics_dict = row.get('metrics',{})
        factor_audit.extend(dict(sample_id=row['sample_id'],technology=row['technology'],**entry)
            for entry in metrics_dict.get('LCA_factor_audit',[]))
        for key,value in metrics_dict.get('annual_inventory',{}).items():
            inventories.append(dict(sample_id=row['sample_id'],technology=row['technology'],
                                    activity=key,annual_quantity=value))
            if row['technology'] in templates:
                audit_entry = next((a for a in metrics_dict.get('LCA_factor_audit',[]) if a['activity']==key),{})
                templates[row['technology']][key] = audit_entry.get('factor_kgCO2e_per_unit')
    pd.DataFrame(inventories,columns=['sample_id','technology','activity','annual_quantity']).to_csv(run_dir/'alternative_inventories.csv',index=False)
    write_json(run_dir/'lca_factors_required.json',templates)
    pd.DataFrame(factor_audit).to_csv(run_dir/'lca_factor_audit.csv',index=False)


def completed_samples(config,pending,run_dir):
    """Supervise persistent processes; replace only a crashed/timed-out worker."""
    if not pending:
        return
    context = multiprocessing.get_context('spawn')
    tasks = iter((sample,technology) for sample in pending for technology in TECHNOLOGIES)
    slots = [None]*min(config['workers'],len(pending)*len(TECHNOLOGIES))
    finished = {}
    exhausted = False
    timeout = config.get('max_technology_seconds',300)

    def start():
        job = Path(run_dir)/'logs'/f'.worker_{os.getpid()}_{time.time_ns()}.job.json'
        process = context.Process(target=_persistent_worker,args=(str(job),config,str(run_dir)))
        process.start()
        return dict(process=process,job=job,task=None,completed=0)

    def stop(slot):
        if slot is None:
            return
        process = slot['process']
        if process.is_alive() and slot['task'] is None:
            write_json(slot['job'],None)
            process.join(2)
        if process.is_alive():
            process.terminate()
            process.join(5)
            if process.is_alive():
                process.kill()
                process.join(5)
        slot['job'].unlink(missing_ok=True)
        process.close()

    try:
        while True:
            active = False
            for index,slot in enumerate(slots):
                if slot is not None and slot['task'] is not None:
                    active = True
                    row = None
                    if slot['result_path'].exists():
                        try:
                            row = json.loads(slot['result_path'].read_text(encoding='utf-8'))
                        except (PermissionError, FileNotFoundError):
                        # Temporary Windows access conflict; retry next scheduler pass.
                            row = None
                            row = json.loads(slot['result_path'].read_text())
                    expired = time.monotonic()-slot['started']>timeout
                    dead = not slot['process'].is_alive()
                    if row is None and (expired or dead):
                        sample,technology = slot['task']
                        code = slot['process'].exitcode
                        code_text = f'{code} (0x{code & 0xffffffff:08X})' if code is not None else 'unavailable'
                        row = dict(sample_id=sample['sample_id'],technology=technology,
                            status='TIMEOUT' if expired and not dead else 'FAILED',metrics={},
                            log=str(Path(run_dir)/'logs'/f'{sample["sample_id"]:05d}_{technology}.log'),
                            worker_pid=slot['process'].pid,
                            error=f'Exceeded {timeout} seconds; no result accepted.' if expired and not dead else
                                  f'Worker exited with code {code_text}; inspect the plain log and fatal stack log.')
                    if row is not None:
                        sample,technology = slot['task']
                        row['wall_elapsed_seconds'] = time.monotonic()-slot['started']
                        row.setdefault('elapsed_seconds',row['wall_elapsed_seconds'])
                        finished.setdefault(sample['sample_id'],{})[technology] = row
                        slot['completed'] += 1
                        slot['task'] = None if not (expired or dead) else slot['task']
                        if expired or dead or slot['completed']>=config.get('worker_tasks_before_restart',20):
                            stop(slot)
                            slots[index] = slot = None
                        write_json(Path(run_dir)/'logs'/f'{sample["sample_id"]:05d}_{technology}.result.json',row)
                        if len(finished[sample['sample_id']])==len(TECHNOLOGIES):
                            records = finished.pop(sample['sample_id'])
                            yield sample,[records[t] for t in TECHNOLOGIES]
                if slot is None or slot['task'] is None:
                    if not exhausted:
                        task = next(tasks,None)
                        if task is None:
                            exhausted = True
                        else:
                            if slot is None:
                                slots[index] = slot = start()
                            sample,technology = task
                            result = Path(run_dir)/'logs'/f'{sample["sample_id"]:05d}_{technology}.result.json'
                            result.unlink(missing_ok=True)
                            slot['task'] = task
                            slot['started'] = time.monotonic()
                            slot['result_path'] = result
                            write_json(slot['job'],task)
                            active = True
            if exhausted and not any(s is not None and s['task'] is not None for s in slots):
                break
            if active:
                time.sleep(.05)
    finally:
        for slot in slots:
            stop(slot)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=DEPLOY/'comparison_lhs_config.json')
    parser.add_argument('--samples',type=int)
    parser.add_argument('--variation',type=float)
    parser.add_argument('--variation-mode',choices=('relative','absolute'))
    parser.add_argument('--workers',type=int)
    parser.add_argument('--case-timeout',type=float,help='Per-technology time limit in seconds; default 300')
    parser.add_argument('--limit',type=int,help='Run only the first N sampled cases; retain the complete design')
    parser.add_argument('--prepare-only',action='store_true',help='Save the LHS design without loading HTL')
    parser.add_argument('--preflight-only',action='store_true',help='Run the nominal blend once, retaining the complete LHS design')
    parser.add_argument('--resume',type=Path,help='Resume a matching checkpoint folder under deploy/results')
    parser.add_argument('--continue-on-baseline-failure',action='store_true')
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit<1:
        parser.error('--limit must be positive')
    config = json.loads(args.config.read_text(encoding='utf-8-sig'))
    if args.case_timeout is not None:
        if not isfinite(args.case_timeout) or args.case_timeout<=0:
            parser.error('--case-timeout must be positive and finite')
        config['max_technology_seconds'] = args.case_timeout
    for key in ('samples','variation','variation_mode','workers'):
        value = getattr(args,key)
        if value is not None:
            config['n_samples' if key=='samples' else key] = value
    bootstrap(config.get('python_paths',()))
    validate(config)
    stamp = datetime.now(ZoneInfo('America/New_York')).strftime('%Y%m%d_%H%M%S_%f')
    run_dir = (args.resume if args.resume else RESULTS/f'lhs_comparison_{stamp}').resolve()
    if not run_dir.is_relative_to(RESULTS.resolve()):
        raise ValueError('All results must stay under deploy/results')
    if args.resume and not (run_dir/'manifest.json').is_file():
        raise ValueError('Resume folder has no manifest')
    run_dir.mkdir(parents=True,exist_ok=True)
    (run_dir/'logs').mkdir(exist_ok=True)
    manifest_path = run_dir/'manifest.json'
    if args.resume:
        saved_manifest = json.loads(manifest_path.read_text())
        if saved_manifest['config_sha256']!=fingerprint(config):
            raise ValueError('Resume configuration differs from the saved run')
    else:
        write_json(manifest_path,dict(config=config,config_sha256=fingerprint(config),
            created_new_york=datetime.now(ZoneInfo('America/New_York')).isoformat(),
            python=sys.executable,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            execution='Supervised persistent spawn workers; fresh native systems per technology/sample',
            runtime_environment={key:os.environ.get(key) for key in
                ('NUMBA_DISABLE_JIT','PYTHONPATH','PYTHONNOUSERSITE')},
            feedstock_data_sha256=hashlib.sha256((DEPLOY/'feedstocks.py').read_bytes()).hexdigest(),
            sampling='Four-dimensional LHS before per-component dry normalization; fixed wet blend weights and moisture',
            pairwise_lca_enabled=config['lca'].get('boundaries_confirmed_comparable',False)))
    samples = generate_samples(config)
    sample_text = ''.join(json.dumps(sample,allow_nan=False)+'\n' for sample in samples)
    design_hash = hashlib.sha256(sample_text.encode()).hexdigest()
    saved_manifest = json.loads(manifest_path.read_text())
    if args.resume:
        if saved_manifest.get('sample_design_sha256',design_hash)!=design_hash:
            raise ValueError('Regenerated sample design differs; start a new run')
        current_feed_hash = hashlib.sha256((DEPLOY/'feedstocks.py').read_bytes()).hexdigest()
        if saved_manifest.get('feedstock_data_sha256',current_feed_hash)!=current_feed_hash:
            raise ValueError('Deployment feedstock data changed; start a new run')
    saved_manifest['sample_design_sha256'] = design_hash
    saved_manifest['nominal_capacity'] = capacity_basis(config,make_sample(config,-1,[.5]*4)['composition'])
    write_json(manifest_path,saved_manifest)
    with (run_dir/'samples.jsonl').open('w',encoding='utf-8') as stream:
        for sample in samples:
            stream.write(json.dumps(sample,allow_nan=False)+'\n')
    import pandas as pd
    pd.DataFrame([dict(sample_id=s['sample_id'],Water=s['composition']['Water'],
        **{'capacity_'+key:value for key,value in s['capacity'].items()},
        ad_biodegradability=s['ad_biodegradability'],
        **{'wet_'+k:v for k,v in s['composition'].items() if k!='Water'},
        **{'dry_'+k:v for k,v in s['dry_composition'].items()},
        **{'lhs_'+k:v for k,v in s['lhs_quantiles'].items()},
        **{'perturbation_'+k:v for k,v in s['perturbations'].items()}) for s in samples]).to_csv(run_dir/'samples.csv',index=False)
    print(f'Saved {len(samples)} paired LHS inputs: {run_dir}',flush=True)
    write_json(run_dir/'run_status.json',dict(status='PREPARED',planned_samples=len(samples)))
    if args.prepare_only:
        return run_dir
    baseline = make_sample(config,-1,[.5]*4)
    baseline_path = run_dir/'baseline_results.json'
    previous = json.loads(baseline_path.read_text()) if baseline_path.exists() else []
    baseline_rows = (previous if args.resume and not args.preflight_only and previous and
                     all(r['status']=='OK' for r in previous) else run_sample(config,baseline,run_dir))
    write_json(run_dir/'baseline_results.json',baseline_rows)
    export_results(run_dir,baseline_rows,config)
    failures = [r for r in baseline_rows if r['status']!='OK']
    print('Baseline: '+', '.join(r['technology']+'='+r['status'] for r in baseline_rows),flush=True)
    write_json(run_dir/'run_status.json',dict(status='PREFLIGHT_FAILED' if failures else 'PREFLIGHT_OK',
        planned_samples=len(samples),baseline_status={r['technology']:r['status'] for r in baseline_rows}))
    if args.preflight_only:
        return run_dir
    if failures and not args.continue_on_baseline_failure:
        raise RuntimeError('Baseline failed; full simulation stopped to avoid repeating an invalid setup. Inspect baseline_results.json and logs; the full LHS design is saved.')
    checkpoint = run_dir/'results.jsonl'
    history = [json.loads(line) for line in checkpoint.read_text().splitlines() if line.strip()] if checkpoint.exists() else []
    latest = {(r['sample_id'],r['technology']):r for r in history}
    rows = list(latest.values())
    complete = {sample['sample_id'] for sample in samples if
                all(any(r['sample_id']==sample['sample_id'] and r['technology']==t for r in rows) for t in TECHNOLOGIES)}
    selected = samples[:args.limit] if args.limit is not None else samples
    pending = [s for s in selected if s['sample_id'] not in complete]
    write_json(run_dir/'run_status.json',dict(status='RUNNING',planned_samples=len(samples),selected_samples=len(selected)))
    for number,(sample,records) in enumerate(completed_samples(config,pending,run_dir),1):
        with checkpoint.open('a',encoding='utf-8') as stream:
            for record in records:
                stream.write(json.dumps(clean(record),allow_nan=False)+'\n')
            stream.flush()
        latest.update({(r['sample_id'],r['technology']):r for r in records})
        rows = list(latest.values())
        print(f'{number}/{len(pending)} samples: '+', '.join(r['technology']+'='+r['status'] for r in records),flush=True)
        if number%10==0:
            export_results(run_dir,rows,config)
    export_results(run_dir,rows,config)
    write_json(run_dir/'run_status.json',dict(status='COMPLETE' if len({r['sample_id'] for r in rows})==len(samples) else 'PARTIAL_COMPLETE',planned_samples=len(samples),
        selected_samples=len(selected),
        completed_samples=len({r['sample_id'] for r in rows}),completed_technology_cases=len(rows),
        failed_technology_cases=sum(r['status']!='OK' for r in rows)))
    print(f'Results saved: {run_dir}',flush=True)
    return run_dir


if __name__=='__main__':
    multiprocessing.freeze_support()
    main()
